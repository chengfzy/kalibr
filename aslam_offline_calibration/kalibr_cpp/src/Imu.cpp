#include "cc/Imu.h"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Imu.h>
#include <iostream>
#include "aslam/backend/BSplineMotionError.hpp"
#include "aslam/backend/MEstimatorPolicies.hpp"
#include "aslam/splines/EuclideanBSplineDesignVariable.hpp"
#include "cc/Heading.hpp"
#include "cc/Util.h"

using namespace cc;
using namespace std;
using namespace fmt;
using namespace aslam;
using namespace Eigen;
using namespace rosbag;
using namespace aslam::backend;
using namespace kalibr_errorterms;

Imu::Imu(const std::string& bagFile, const ImuParameters& imuParams, const ros::Time& startTime,
         const ros::Time& endTime)
    : params_(imuParams), gyroBias(3), accBias(3) {
    cout << Paragraph("Initialize IMU");
    cout << format("\tDataset: {}", bagFile) << endl;
    cout << format("\tTopic: {}", imuParams.topic) << endl;
    cout << format("\tTime: [{}, {}]", startTime, endTime) << endl;

    // read IMU data and add to IMU measurements
    Bag bag(bagFile);
    View view(bag, TopicQuery(imuParams.topic), ros::Time(startTime), ros::Time(endTime));
    cout << format("bag message size = {}, t = [{}, {}]", view.size(), view.getBeginTime(), view.getEndTime()) << endl;
    Matrix3d accR = imuParams.getAccNoiseDiscrete() * imuParams.getAccNoiseDiscrete() * Matrix3d::Identity();
    Matrix3d gyroR = imuParams.getGyroNoiseDiscrete() * imuParams.getGyroNoiseDiscrete() * Matrix3d::Identity();
    for (auto it = view.begin(); it != view.end(); ++it) {
        sensor_msgs::Imu::ConstPtr v = it->instantiate<sensor_msgs::Imu>();
        if (v) {
            data.emplace_back(ImuMeasurement(
                Time(v->header.stamp.toSec()),
                Vector3d(v->linear_acceleration.x, v->linear_acceleration.y, v->linear_acceleration.z),
                Vector3d(v->angular_velocity.x, v->angular_velocity.y, v->angular_velocity.z), accR, gyroR));
        }
    }
    cout << format("obtain {}/{} IMU data", data.size(), view.size()) << endl;
}

void Imu::initBiasSpline(const bsplines::BSplinePose& poseSpline, int splineOrder, int biasKnotsPerSecond) {
    double start = poseSpline.tMin();
    double end = poseSpline.tMax();
    int knots = round((end - start) * biasKnotsPerSecond);
    cout << Section(format("Initializing the bias Splines with {} knots", knots)) << endl;

    gyroBias = bsplines::BSpline(splineOrder);
    gyroBias.initConstantSpline(start, end, knots, gyroBiasPrior);
    accBias = bsplines::BSpline(splineOrder);
    accBias.initConstantSpline(start, end, knots, Vector3d::Zero());
}

void Imu::addDesignVariable(calibration::OptimizationProblem& problem) {
    gyroBiasDesignVar = boost::make_shared<splines::EuclideanBSplineDesignVariable>(gyroBias);
    for (size_t i = 0; i < gyroBiasDesignVar->numDesignVariables(); ++i) {
        auto v = calibration::OptimizationProblem::DesignVariableSP(gyroBiasDesignVar->designVariable(i),
                                                                    sm::null_deleter());
        v->setActive(true);
        problem.addDesignVariable(v, kHelperGroupId);
    }

    accBiasDesignVar = boost::make_shared<splines::EuclideanBSplineDesignVariable>(accBias);
    for (size_t i = 0; i < accBiasDesignVar->numDesignVariables(); ++i) {
        auto v =
            calibration::OptimizationProblem::DesignVariableSP(accBiasDesignVar->designVariable(i), sm::null_deleter());
        v->setActive(true);
        problem.addDesignVariable(v, kHelperGroupId);
    }

    // q_IB, the rotation for multiple IMUs
    qIBDesignVar = boost::make_shared<backend::RotationQuaternion>(qIBPrior);
    qIBDesignVar->setActive(false);
    problem.addDesignVariable(qIBDesignVar, kHelperGroupId);

    // rB, the translation for multiple IMUs
    rBDesignVariable = boost::make_shared<backend::EuclideanPoint>(Vector3d(0, 0, 0));
    rBDesignVariable->setActive(false);
    problem.addDesignVariable(rBDesignVariable, kHelperGroupId);
}

void Imu::addAccelerometerErrorTerms(calibration::OptimizationProblem& problem,
                                     const boost::shared_ptr<splines::BSplinePoseDesignVariable>& poseDesignVariable,
                                     const backend::EuclideanExpression& gravity, const double& sigma,
                                     const double& accNoiseScale) {
    cout << Section("Add Accelerometer Error Terms");
    double weight = 1.0 / accNoiseScale;

    // M estimator
    boost::shared_ptr<backend::MEstimator> estimator;
    if (sigma > 0.0) {
        estimator = boost::make_shared<backend::HuberMEstimator>(sigma);
    } else {
        estimator = boost::make_shared<backend::NoMEstimator>();
    }

    // add error terms
    size_t skipedNum{0};
    for (auto& v : data) {
        double t = v.stamp.toSec();
        if (poseDesignVariable->spline().tMin() < t && t < poseDesignVariable->spline().tMax()) {
            RotationExpression Rbw = poseDesignVariable->orientation(t).inverse();
            EuclideanExpression aw = poseDesignVariable->linearAcceleration(t);
            EuclideanExpression bi = accBiasDesignVar->toEuclideanExpression(t, 0);
            EuclideanExpression wb = poseDesignVariable->angularVelocityBodyFrame(t);
            EuclideanExpression wDotB = poseDesignVariable->angularAccelerationBodyFrame(t);
            RotationExpression Rib = qIBDesignVar->toExpression();
            EuclideanExpression rB = rBDesignVariable->toExpression();
            // Rib, rB is the rotation and position for multiple IMUs
            EuclideanExpression a = Rib * (Rbw * (aw - gravity) + wDotB.cross(rB) + wb.cross(wb.cross(rB)));
            auto error = boost::make_shared<EuclideanError>(v.acc, v.accInvR * weight, a + bi);
            error->setMEstimatorPolicy(estimator);
            accErrors.emplace_back(error);
            problem.addErrorTerm(error);
        } else {
            ++skipedNum;
        }
    }

    cout << format("\tAdd {} of {} accelerometer error terms (skipped {} out of bounds measurements)",
                   data.size() - skipedNum, data.size(), skipedNum);
}

void Imu::addGyroscopeErrorTerms(calibration::OptimizationProblem& problem,
                                 const boost::shared_ptr<splines::BSplinePoseDesignVariable>& poseDesignVariable,
                                 const double& sigma, const double& gyroNoiseScale) {
    cout << Section("Add Gyroscope Error Terms");
    double weight = 1.0 / gyroNoiseScale;

    // M estimator
    boost::shared_ptr<backend::MEstimator> estimator;
    if (sigma > 0.0) {
        estimator = boost::make_shared<backend::HuberMEstimator>(sigma);
    } else {
        estimator = boost::make_shared<backend::NoMEstimator>();
    }

    // add error terms
    size_t skipedNum{0};
    for (auto& v : data) {
        double t = v.stamp.toSec();
        if (poseDesignVariable->spline().tMin() < t && t < poseDesignVariable->spline().tMax()) {
            EuclideanExpression wb = poseDesignVariable->angularVelocityBodyFrame(t);
            EuclideanExpression bi = gyroBiasDesignVar->toEuclideanExpression(t, 0);
            RotationExpression Rib = qIBDesignVar->toExpression();
            EuclideanExpression w = Rib * wb;
            auto error = boost::make_shared<EuclideanError>(v.gyro, v.gyroInvR * weight, w + bi);
            error->setMEstimatorPolicy(estimator);
            gyroErrors.emplace_back(error);
            problem.addErrorTerm(error);
        } else {
            ++skipedNum;
        }
    }

    cout << format("\tAdd {} of {} gyroscope error terms (skipped {} out of bounds measurements)",
                   data.size() - skipedNum, data.size(), skipedNum);
}

void Imu::addBiasMotionTerms(aslam::calibration::OptimizationProblem& problem) {
    cout << Section("Add IMU Bias Motion Terms");
    Matrix3d gyroW = Matrix3d::Identity() / (params_.gyroRandomWalk * params_.gyroRandomWalk);
    auto gyroBiasMotionError = boost::make_shared<BSplineMotionError<splines::EuclideanBSplineDesignVariable>>(
        gyroBiasDesignVar.get(), gyroW, 1);
    problem.addErrorTerm(gyroBiasMotionError);

    Matrix3d accW = Matrix3d::Identity() / (params_.accRandomWalk * params_.accRandomWalk);
    auto accBiasMotionError = boost::make_shared<BSplineMotionError<splines::EuclideanBSplineDesignVariable>>(
        accBiasDesignVar.get(), accW, 1);
    problem.addErrorTerm(accBiasMotionError);
}