#pragma once
#include <ros/ros.h>
#include <Eigen/Core>
#include "aslam/PinholeUndistorter.hpp"
#include "aslam/Time.hpp"
#include "aslam/backend/CameraDesignVariable.hpp"
#include "aslam/backend/CovarianceReprojectionError.hpp"
#include "aslam/calibration/core/OptimizationProblem.h"
#include "aslam/cameras.hpp"
#include "aslam/cameras/GridCalibrationTargetAprilgrid.hpp"
#include "aslam/cameras/GridDetector.hpp"
#include "aslam/splines/BSplinePoseDesignVariable.hpp"
#include "cc/AprilTargetParameters.h"
#include "cc/CameraParameters.h"

namespace cc {

/**
 * @brief Calibration for rolling shutter camera with pin hole projection
 */
class RollingShutterCameraCalibrator {
  public:
    struct Options {
        double deltaX = 1e-8;
        double deltaJ = 1e-4;
        int maxIterationNumber = 30;
        int maxKnotPlacementIterations = 10;
        bool adaptiveKnotPlacement = true;
        double timeOffsetConstantSparsityPattern = 0.08;
        double inverseFeatureCovariance = 1. / 0.26;
        int splineOrder = 4;
        double timeOffsetPadding = 0.05;
        int numberOfKnots = -1;
        int frameRate = 30;
    };

  public:
    RollingShutterCameraCalibrator(const std::string& bagFile, const CameraParameters& cameraParams,
                                   const AprilTargetParameters& targetParams, const Options& options,
                                   const ros::Time& startTime = ros::TIME_MIN,
                                   const ros::Time& endTime = ros::TIME_MAX);

  public:
    // calibration
    void calibrate();

  private:
    // extract observations from ros bag image
    void extractObservations();

    // init intrinsics
    void initIntrinsics();

    // init extrinsics
    void initExtrinsics();

    // init pose B-Spline
    bsplines::BSplinePose initPoseBSpline(int splineOrder, const double& timeOffsetPadding, int numberOfKnots,
                                          const double& frameRate);

    // build the optimization problem
    void buildProblem(const bsplines::BSplinePose& poseSpline, const Eigen::Matrix<double, 6, 6>& W);

    // solve problem
    void solve();

    // print results
    void printResult();

  public:
    // type define
    using CameraGeometry = aslam::cameras::DistortedPinholeRsCameraGeometry;
    using Frame = aslam::Frame<CameraGeometry>;
    using Distortion = aslam::cameras::RadialTangentialDistortion;
    using Projection = aslam::cameras::PinholeProjection<Distortion>;
    using Shutter = aslam::cameras::RollingShutter;
    using CameraDesignVariable = aslam::backend::CameraDesignVariable<CameraGeometry>;
    using AdaptiveCovarianceReprojectionError =
        aslam::backend::CovarianceReprojectionError<aslam::Frame<CameraGeometry>>;

    Options options;  // options

    std::string bagFile;                 // ros bag file
    ros::Time startTime;                 // start time for ros bag
    ros::Time endTime;                   // end time for ros bag
    CameraParameters cameraParams;       // camera parameters
    AprilTargetParameters targetParams;  // target parameters

    CameraGeometry cameraGeometry;                                                    // camera geometry
    boost::shared_ptr<CameraDesignVariable> cameraDeignVariable;                      // camera design variable
    aslam::cameras::GridDetector detector;                                            // grid detector
    std::vector<aslam::cameras::GridCalibrationTargetObservation> observations;       // observations
    std::vector<boost::shared_ptr<Frame>> frames;                                     // frames
    boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable> poseDesignVariable;  // pose spline design variable
    std::vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>> reprojectionErrors;  // reprojection errors

    boost::shared_ptr<aslam::calibration::OptimizationProblem> problem;  // optimization problem
};

}  // namespace cc