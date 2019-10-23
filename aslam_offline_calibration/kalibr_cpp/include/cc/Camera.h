#pragma once
#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <aslam/cameras/GridDetector.hpp>
#include <vector>
#include "aslam/PinholeUndistorter.hpp"
#include "aslam/Time.hpp"
#include "aslam/backend/Scalar.hpp"
#include "aslam/backend/SimpleReprojectionError.hpp"
#include "aslam/cameras.hpp"
#include "aslam/cameras/GridCalibrationTargetAprilgrid.hpp"
#include "aslam/splines/BSplinePoseDesignVariable.hpp"
#include "bsplines/BSplinePose.hpp"
#include "cc/AprilTargetParameters.h"
#include "cc/CameraParameters.h"
#include "cc/Imu.h"
#include "cc/TransformationDesignVariable.h"
#include "sm/kinematics/Transformation.hpp"

namespace cc {

class Camera {
  public:
    Camera(const std::string& bagFile, const CameraParameters& cameraParams, const AprilTargetParameters& targetParams,
           const ros::Time& startTime = ros::TIME_MIN, const ros::Time& endTime = ros::TIME_MAX);

  public:
    // initialize a pose spine using camera poses(pose  spine = T_wb)
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

  public:
    // type define
    using Frame = aslam::Frame<aslam::cameras::DistortedPinholeCameraGeometry>;
    using KeyPoint = aslam::Keypoint<2>;
    using ReprojectionError = aslam::backend::SimpleReprojectionError<Frame>;
    using Undistorter = aslam::PinholeUndistorter<aslam::cameras::RadialTangentialDistortion, aslam::cameras::NoMask>;

    CameraParameters cameraParams;                                                // camera parameters
    double cornerUncertainty = 1.0;                                               // corner uncertainty
    sm::kinematics::Transformation extrinsic = sm::kinematics::Transformation();  // extrinsic, T_CB set to default
    double timeshiftCameraToImuPrior = 0;                                         // timeshift between camera and IMU
    Eigen::Vector3d gravity = Eigen::Vector3d(9.80655, 0., 0.);                   // gravity

    aslam::cameras::DistortedPinholeCameraGeometry geometry;                     // camera geometry
    aslam::cameras::GridDetector detector;                                       // grid detector
    std::vector<aslam::cameras::GridCalibrationTargetObservation> observations;  // observations

    boost::shared_ptr<cc::TransformationDesignVariable> TcbDesignVar;      // extrinsic TCB design variable
    boost::shared_ptr<aslam::backend::Scalar> cameraTimeToImuDesignVar;    // time delay design variable
    std::vector<boost::shared_ptr<ReprojectionError>> reprojectionErrors;  // reprojection error
};

}  // namespace cc