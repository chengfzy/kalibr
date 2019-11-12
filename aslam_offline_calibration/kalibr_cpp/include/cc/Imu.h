#pragma once
#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include "aslam/Time.hpp"
#include "aslam/backend/EuclideanPoint.hpp"
#include "aslam/backend/RotationQuaternion.hpp"
#include "aslam/calibration/core/OptimizationProblem.h"
#include "aslam/splines/BSplinePoseDesignVariable.hpp"
#include "aslam/splines/EuclideanBSplineDesignVariable.hpp"
#include "bsplines/BSplinePose.hpp"
#include "cc/ImuParameters.h"
#include "kalibr_errorterms/AccelerometerError.hpp"
#include "kalibr_errorterms/EuclideanError.hpp"
#include "kalibr_errorterms/GyroscopeError.hpp"

namespace cc {

class ImuMeasurement {
  public:
    ImuMeasurement(const aslam::Time& stamp, const Eigen::Vector3d& acc, const Eigen::Vector3d& gyro,
                   const Eigen::Matrix3d& accR, const Eigen::Matrix3d& gyroR)
        : stamp(stamp), acc(acc), gyro(gyro), accR(accR), gyroR(gyroR) {
        accInvR = accR.inverse();
        gyroInvR = gyroR.inverse();
    }

  public:
    aslam::Time stamp;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
    Eigen::Matrix3d accR;
    Eigen::Matrix3d accInvR;
    Eigen::Matrix3d gyroR;
    Eigen::Matrix3d gyroInvR;
};

class Imu {
  public:
    Imu(const std::string& bagFile, const ImuParameters& imuParams, const ros::Time& startTime = ros::TIME_MIN,
        const ros::Time& endTime = ros::TIME_MAX);

    void initBiasSpline(const bsplines::BSplinePose& poseSpline, int splineOrder, int biasKnotsPerSecond);

    void addDesignVariable(aslam::calibration::OptimizationProblem& problem);

    void addAccelerometerErrorTerms(
        aslam::calibration::OptimizationProblem& problem,
        const boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable>& poseDesignVariable,
        const aslam::backend::EuclideanExpression& gravity, const double& sigma, const double& accNoiseScale);

    void addGyroscopeErrorTerms(aslam::calibration::OptimizationProblem& problem,
                                const boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable>& poseDesignVariable,
                                const double& sigma, const double& gyroNoiseScale);

    void addBiasMotionTerms(aslam::calibration::OptimizationProblem& problem);

  public:
    ImuParameters params_;             // IMU parameters
    std::vector<ImuMeasurement> data;  // IMU measurements

    bsplines::BSpline gyroBias;                                                           // gyro bias spline
    bsplines::BSpline accBias;                                                            // acc bias spline
    boost::shared_ptr<aslam::splines::EuclideanBSplineDesignVariable> gyroBiasDesignVar;  // gyro bias design variable
    boost::shared_ptr<aslam::splines::EuclideanBSplineDesignVariable> accBiasDesignVar;   // acc bias design variable

    Eigen::Vector3d gyroBiasPrior = Eigen::Vector3d::Zero();

    // qIB and rB is the relative rotation and position w.r.t anchor IMU, used for multiple IMU system
    Eigen::Vector4d qIBPrior = Eigen::Vector4d(0, 0, 0, 1);
    boost::shared_ptr<aslam::backend::RotationQuaternion> qIBDesignVar;
    boost::shared_ptr<aslam::backend::EuclideanPoint> rBDesignVariable;

    bool isReference = true;
    bool estimateTimeDelay = false;

    std::vector<boost::shared_ptr<kalibr_errorterms::EuclideanError>> accErrors;
    std::vector<boost::shared_ptr<kalibr_errorterms::EuclideanError>> gyroErrors;
};

}  // namespace cc
