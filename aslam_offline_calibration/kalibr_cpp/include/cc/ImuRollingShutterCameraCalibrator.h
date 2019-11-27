#pragma once
#include "Imu.h"
#include "RollingShutterCamera.h"
#include "aslam/backend/EuclideanDirection.hpp"
#include "aslam/backend/Optimizer2.hpp"
#include "aslam/backend/Optimizer2Options.hpp"
#include "aslam/calibration/core/OptimizationProblem.h"
#include "aslam/splines/BSplinePoseDesignVariable.hpp"

namespace cc {

class ImuRollingShutterCameraCalibrator {
  public:
    ImuRollingShutterCameraCalibrator(const RollingShutterCamera& camera, const Imu& imu);

  public:
    void buildProblem(int splineOrder = 6, int poseKnotsPerSecond = 100, int biasKnotsPerSecond = 50,
                      bool doPoseMotionError = false, const double& mrTranslationVariance = 1e6,
                      const double& mrRotationVariance = 1e5, bool doBiasMotionError = true, int blakeZisserCam = -1,
                      const double& huberAccel = -1, const double& huberGyro = -1, bool noTimeCalibration = false,
                      bool noChainExtrinsics = true, int maxIterations = 30, const double& gyroNoiseScale = 1.0,
                      const double& accelNoiseScale = 1.0, const double& timeOffsetPadding = 0.03,
                      bool verbose = false);

    void optimize(int maxIterations = 30, bool recoverCov = false);

    void printErrorStatistics();

    void printResult(bool withCov = false);

  private:
    RollingShutterCamera camera;  // camera
    Imu imu;                      // IMU

    aslam::calibration::OptimizationProblem problem;                             // problem
    boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable> poseDesignVar;  // pose design variable
    boost::shared_ptr<aslam::backend::EuclideanDirection> gravityDesignVar;      // gravity design variable
    aslam::backend::EuclideanExpression gravityExpression;

    aslam::backend::Optimizer2 optimizer;
    aslam::backend::Optimizer2Options optimizerOptions;
    Eigen::Matrix<double, 6, 1> transVar;
    Eigen::MatrixXd timesVar;
};

}  // namespace cc
