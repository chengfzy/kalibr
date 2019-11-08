#include "cc/Calibrator.h"
#include <fmt/format.h>
#include "aslam/backend/BlockCholeskyLinearSystemSolver.hpp"
#include "aslam/backend/LevenbergMarquardtTrustRegionPolicy.hpp"
#include "aslam/calibration/core/IncrementalEstimator.h"
#include "cc/Heading.hpp"
#include "cc/Util.h"

using namespace std;
using namespace cc;
using namespace fmt;
using namespace Eigen;
using namespace aslam;
using namespace aslam::backend;
using namespace aslam::calibration;

Calibrator::Calibrator(const cc::Camera& camera, const cc::Imu& imu) : camera(camera), imu(imu) {}

void Calibrator::buildProblem(int splineOrder, int poseKnotsPerSecond, int biasKnotsPerSecond, bool doPoseMotionError,
                              const double& mrTranslationVariance, const double& mrRotationVariance,
                              bool doBiasMotionError, int blakeZisserCam, const double& huberAccel,
                              const double& huberGyro, bool noTimeCalibration, bool, int maxIterations,
                              const double& gyroNoiseScale, const double& accelNoiseScale,
                              const double& timeOffsetPadding, bool) {
    cout << Section("Build Problem");
    cout << format("\tSpline order: {}", splineOrder) << endl;
    cout << format("\tPose knots per second: {}", poseKnotsPerSecond) << endl;
    cout << format("\tDo pose motion regularization: {}", doPoseMotionError) << endl;
    cout << format("\t\txddot translation variance: {}", mrTranslationVariance) << endl;
    cout << format("\t\txddot rotation variance: {}", mrRotationVariance) << endl;
    cout << format("\tBias knots per second: {}", biasKnotsPerSecond) << endl;
    cout << format("\tDo bias motion regularization: {}", doBiasMotionError) << endl;
    cout << format("\tBlake-Zisserman on reprojection errors: {}", blakeZisserCam) << endl;
    cout << format("\tAcceleration Huber width(sigma): {}", huberAccel) << endl;
    cout << format("\tGyroscope Huber width(sigma): {}", huberGyro) << endl;
    cout << format("\tDo time calibration: {}", !noTimeCalibration) << endl;
    cout << format("\tMax iteration: {}", maxIterations) << endl;
    cout << format("\tTime offset padding: {}", timeOffsetPadding) << endl;

    camera.findTimeShiftCameraImuPrior(imu);
    camera.findOrientationPriorCameraToImu(imu);
    const Vector3d& gravity = camera.gravity;

    // init optimization problem
    // initialize a pose spline using the camera poses in the cam chain
    bsplines::BSplinePose poseSpline =
        camera.initPoseSplineFromCamera(splineOrder, poseKnotsPerSecond, timeOffsetPadding);

    // initialize bias spines for IMU
    imu.initBiasSpline(poseSpline, splineOrder, biasKnotsPerSecond);

    // init design variables
    poseDesignVar = boost::make_shared<splines::BSplinePoseDesignVariable>(poseSpline);
    for (size_t i = 0; i < poseDesignVar->numDesignVariables(); ++i) {
        auto v =
            calibration::OptimizationProblem::DesignVariableSP(poseDesignVar->designVariable(i), sm::null_deleter());
        v->setActive(true);
        problem.addDesignVariable(v, kHelperGroupId);
    }
    // add the calibration target orientation design variable
    gravityDesignVar = boost::make_shared<backend::EuclideanDirection>(gravity);
    gravityExpression = gravityDesignVar->toExpression();
    gravityDesignVar->setActive(true);
    problem.addDesignVariable(gravityDesignVar, kHelperGroupId);
    // add design variable for IMU
    imu.addDesignVariable(problem);
    camera.addDesignVariable(problem, noTimeCalibration);

    // add calibration target reprojection error terms
    camera.addErrorTerms(problem, poseDesignVar, camera.TcbDesignVar->toExpression(), blakeZisserCam,
                         timeOffsetPadding);

    // add IMU errors
    imu.addAccelerometerErrorTerms(problem, poseDesignVar, gravityExpression, huberAccel, accelNoiseScale);
    imu.addGyroscopeErrorTerms(problem, poseDesignVar, huberGyro, gyroNoiseScale);
    if (doBiasMotionError) {
        imu.addBiasMotionTerms(problem);
    }

    // add pose motion terms
    if (doPoseMotionError) {
        cout << Section("Add Pose Motion Error");
        // pass
    }
}

void Calibrator::optimize(int maxIterations, bool recoverCov) {
    // options
    optimizerOptions.convergenceDeltaX = 1e-5;
    optimizerOptions.convergenceDeltaJ = 1e-2;
    optimizerOptions.maxIterations = maxIterations;
    optimizerOptions.trustRegionPolicy = boost::make_shared<LevenbergMarquardtTrustRegionPolicy>(10);
    optimizerOptions.linearSystemSolver = boost::make_shared<BlockCholeskyLinearSystemSolver>();

    // run the optimization
    optimizer.options() = optimizerOptions;
    optimizer.setProblem(boost::shared_ptr<OptimizationProblem>(&problem, sm::null_deleter()));
    optimizer.optimize();

    // recover covariance
    if (recoverCov) {
        cout << Section("Recovering Covariance");
        IncrementalEstimator estimator(kCalibrationGroupId);
        estimator.addBatch(boost::shared_ptr<OptimizationProblem>(&problem, sm::null_deleter()), true);
        VectorXd estimatedStd = estimator.getSigma2Theta().diagonal();

        // split and store the variance
        transVar = estimatedStd.head<6>();
        timesVar = estimatedStd.tail(6);
    }
}

void Calibrator::printErrorStatistics() {
    // print out normalized residuals
    cout << Paragraph("Normalized Residuals");
    // evaluate reprojection error
    VectorXd reprojectNormalErrors(camera.reprojectionErrors.size());
    VectorXd reprojectErrors(camera.reprojectionErrors.size());
    for (size_t i = 0; i < camera.reprojectionErrors.size(); ++i) {
        reprojectNormalErrors[i] = sqrt(camera.reprojectionErrors[i]->evaluateError());
        reprojectErrors[i] = camera.reprojectionErrors[i]->error().norm();
    }
    cout << format("Reprojection error:\tmean: {}, std: {}", reprojectNormalErrors.mean(),
                   sqrt((reprojectNormalErrors.rowwise() - reprojectNormalErrors.colwise().mean())
                            .rowwise()
                            .squaredNorm()
                            .mean()))
         << endl;

    // evaluate gyroscope error
    VectorXd gyroNormalErrors(imu.gyroErrors.size());
    VectorXd gyroErrors(imu.gyroErrors.size());
    for (size_t i = 0; i < gyroNormalErrors.size(); ++i) {
        gyroNormalErrors[i] = sqrt(imu.gyroErrors[i]->evaluateError());
        gyroErrors[i] = imu.gyroErrors[i]->error().norm();
    }
    cout << format(
                "Gyroscope error:\tmean: {}, std: {}", gyroNormalErrors.mean(),
                sqrt((gyroNormalErrors.rowwise() - gyroNormalErrors.colwise().mean()).rowwise().squaredNorm().mean()))
         << endl;

    // evaluate acceleration error
    VectorXd accNormalErrors(imu.accErrors.size());
    VectorXd accErrors(imu.accErrors.size());
    for (size_t i = 0; i < accNormalErrors.size(); ++i) {
        accNormalErrors[i] = sqrt(imu.accErrors[i]->evaluateError());
        accErrors[i] = imu.accErrors[i]->error().norm();
    }
    cout << format("Accelerometer error:\tmean: {}, std: {}", accNormalErrors.mean(),
                   sqrt((accNormalErrors.rowwise() - accNormalErrors.colwise().mean()).rowwise().squaredNorm().mean()))
         << endl;

    // print out residuals
    cout << Paragraph("Residuals");
    cout << format("Reprojection error[px]:\tmean: {}, std: {}", reprojectErrors.mean(),
                   sqrt((reprojectErrors.rowwise() - reprojectErrors.colwise().mean()).rowwise().squaredNorm().mean()))
         << endl;
    cout << format("Gyroscope error[rad/s]:\tmean: {}, std: {}", gyroErrors.mean(),
                   sqrt((gyroErrors.rowwise() - gyroErrors.colwise().mean()).rowwise().squaredNorm().mean()))
         << endl;
    cout << format("Accelerometer error[m/s^2]:\tmean: {}, std: {}", accErrors.mean(),
                   sqrt((accErrors.rowwise() - accErrors.colwise().mean()).rowwise().squaredNorm().mean()))
         << endl;
}

void Calibrator::printResult(bool) {
    cout << Section("Calibration Result");
    cout << "Transforation T_cam0_imu0 (imu0 to cam0, Tci):" << endl;
    cout << camera.getResultTransformationImuToCam().T() << endl;
    cout << endl;
    cout << format("cam0 to imu0 time(t_imu = t_cam + shift): {}[s]", camera.getResultTimeShift()) << endl;
}