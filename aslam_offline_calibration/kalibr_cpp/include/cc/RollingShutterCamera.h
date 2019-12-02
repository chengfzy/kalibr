#pragma once
#include <ros/ros.h>
#include <rosbag/view.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/filesystem.hpp>
#include <vector>
#include "aslam/PinholeUndistorter.hpp"
#include "aslam/Time.hpp"
#include "aslam/backend/CameraDesignVariable.hpp"
#include "aslam/backend/CovarianceReprojectionError.hpp"
#include "aslam/backend/Scalar.hpp"
#include "aslam/backend/SimpleReprojectionError.hpp"
#include "aslam/cameras.hpp"
#include "aslam/cameras/GridCalibrationTargetAprilgrid.hpp"
#include "aslam/cameras/GridDetector.hpp"
#include "aslam/splines/BSplinePoseDesignVariable.hpp"
#include "bsplines/BSplinePose.hpp"
#include "cc/AprilTargetParameters.h"
#include "cc/CameraParameters.h"
#include "cc/Imu.h"
#include "cc/TransformationDesignVariable.h"
#include "sm/kinematics/Transformation.hpp"

namespace cc {

class RollingShutterCamera {
  public:
    RollingShutterCamera(const std::string& bagFile, const CameraParameters& cameraParams,
                         const AprilTargetParameters& targetParams, const ros::Time& startTime = ros::TIME_MIN,
                         const ros::Time& endTime = ros::TIME_MAX);

  public:
    // initialize a pose spine using camera poses(pose spine = T_TB)
    bsplines::BSplinePose initPoseSplineFromCamera(int splineOrder = 6, int poseKnotsPerSecond = 100,
                                                   const double& timeOffsetPadding = 0.02);

    // estimates the timeshift between camera and imu using cross correlation approach
    void findTimeShiftCameraImuPrior(const Imu& imu);

    // estimate IMU-Camera rotation prior
    void findOrientationPriorCameraToImu(Imu& imu);

    void addDesignVariable(aslam::calibration::OptimizationProblem& problem, bool noTimeCalibration = true);

    void addErrorTerms(aslam::calibration::OptimizationProblem& problem,
                       const boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable>& poseDesignVariable,
                       aslam::backend::TransformationExpression& Tcb, int blakeZisserCam,
                       const double& timeOffsetPadding);

    sm::kinematics::Transformation getResultTransformationImuToCam();

    double getResultTimeShift();

  private:
    // detect observation(corners) from ros bag
    void detectObservations(rosbag::View& view, const AprilTargetParameters& targetParams,
                            const boost::filesystem::path& obsPath);

  public:
    // type define
    using CameraGeometry = aslam::cameras::DistortedPinholeRsCameraGeometry;
    using CameraDesignVariable = aslam::backend::CameraDesignVariable<CameraGeometry>;
    using Frame = aslam::Frame<CameraGeometry>;
    using KeyPoint = aslam::Keypoint<2>;
    using ReprojectionError = aslam::backend::SimpleReprojectionError<Frame>;
    using AdaptiveCovarianceReprojectionError =
        aslam::backend::CovarianceReprojectionError<aslam::Frame<CameraGeometry>>;

    CameraParameters cameraParams;   // camera parameters
    double cornerUncertainty = 1.0;  // corner uncertainty
    // NOTE by CC: I'm not sure the extrinsic is T_BC or T_CB, please note that kalibr using JPL conversion, just regard
    // it as what you think it is.
    sm::kinematics::Transformation extrinsic = sm::kinematics::Transformation();  // extrinsic, T_BC set to default
    double timeshiftCameraToImuPrior = 0;                                         // timeshift between camera and IMU
    Eigen::Vector3d gravity = Eigen::Vector3d(9.80655, 0., 0.);                   // gravity

    CameraGeometry geometry;                                                     // camera geometry
    std::vector<aslam::cameras::GridCalibrationTargetObservation> observations;  // observations

    boost::shared_ptr<cc::TransformationDesignVariable> TcbDesignVar;    // extrinsic TCB design variable
    boost::shared_ptr<aslam::backend::Scalar> cameraTimeToImuDesignVar;  // time delay design variable
    std::vector<boost::shared_ptr<Frame>> frames;                        // frames
    std::vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>> reprojectionErrors;  // reprojection error
};

}  // namespace cc