#include "cc/RollingShutterCamera.h"
#include <cv_bridge/cv_bridge.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <glog/logging.h>
#include <rosbag/bag.h>
#include <sensor_msgs/Image.h>
#include <aslam/backend/BlockCholeskyLinearSystemSolver.hpp>
#include <aslam/backend/RotationQuaternion.hpp>
#include <iostream>
#include <numeric>
#include "aslam/backend/EuclideanExpression.hpp"
#include "aslam/backend/EuclideanPoint.hpp"
#include "aslam/backend/MEstimatorPolicies.hpp"
#include "aslam/backend/OptimizationProblem.hpp"
#include "aslam/backend/Optimizer2.hpp"
#include "cc/Heading.hpp"
#include "cc/Util.h"
#include "kalibr_errorterms/GyroscopeError.hpp"
#include "sm/kinematics/RotationVector.hpp"
#include "sm/kinematics/Transformation.hpp"
#include "sm/kinematics/transformations.hpp"

using namespace std;
using namespace cc;
using namespace rosbag;
using namespace fmt;
using namespace Eigen;
using namespace aslam;
using namespace aslam::backend;

RollingShutterCamera::RollingShutterCamera(const string& bagFile, const CameraParameters& cameraParams,
                                           const AprilTargetParameters& targetParams, const ros::Time& startTime,
                                           const ros::Time& endTime)
    : cameraParams(cameraParams) {
    cout << Paragraph("Initialize Camera");
    cout << format("\tDataset: {}", bagFile) << endl;
    cout << format("\tTopic: {}", cameraParams.topic) << endl;
    cout << format("\tTime: [{}, {}]", startTime, endTime) << endl;

    // init ros bag reader
    Bag bag(bagFile);
    View view(bag, TopicQuery(cameraParams.topic), ros::Time(startTime), ros::Time(endTime));
    cout << format("bag message size = {}, t = [{}, {}]", view.size(), view.getBeginTime(), view.getEndTime()) << endl;

    // init camera model
    cameras::RadialTangentialDistortion dist(cameraParams.d[0], cameraParams.d[1], cameraParams.d[2], cameraParams.d[3],
                                             cameraParams.d[4]);
    cameras::PinholeProjection<cameras::RadialTangentialDistortion> proj(
        cameraParams.f[0], cameraParams.f[1], cameraParams.c[0], cameraParams.c[1], cameraParams.resolution[0],
        cameraParams.resolution[1], dist);
    geometry = CameraGeometry(proj, cameras::RollingShutter(cameraParams.lineDelay));

    // check observations data is in folder, if is, then load observations from folder, otherwise detect
    // observations(corners) from rosbag image
    boost::filesystem::path dataPath(bagFile);
    boost::filesystem::path obsPath = dataPath.parent_path() / "imuCameraObservations.xml";
    if (boost::filesystem::exists(obsPath)) {
        // load observations(corners) from folder
        cout << format("load observations from file \"{}\"", obsPath.string()) << endl;
        observations.clear();
        sm::boost_serialization::load_xml(observations, obsPath);
    } else {
        // detect observations(corners) from rosbag image
        detectObservations(view, targetParams, obsPath);
    }
}

// initialize a pose spine using camera poses(pose spine = T_TB)
bsplines::BSplinePose RollingShutterCamera::initPoseSplineFromCamera(int splineOrder, int poseKnotsPerSecond,
                                                                     const double& timeOffsetPadding) {
    // FIXME, the extrinsic is T_BC, not T_CB, even though they are the same during initialization
    Matrix4d Tcb = extrinsic.T();
    bsplines::BSplinePose pose(splineOrder, boost::make_shared<sm::kinematics::RotationVector>());

    // get the checkerboard times
    const size_t kN = observations.size();
    VectorXd times(kN + 2);
    MatrixXd curve(6, kN + 2);
    for (size_t i = 0; i < kN; ++i) {
        times[i + 1] = observations[i].time().toSec() + timeshiftCameraToImuPrior;
#if defined(DebugTest) && false
        cout << format("time = {:.10f}", observations[i].time().toSec()) << endl;
        cout << format("T_t_c.T() = \n{}", observations[i].T_t_c().T()) << endl;
        cout << format("Tcb = \n{}", Tcb) << endl;
        cout << format("T_t_c.T() * Tcb = \n{}", observations[i].T_t_c().T() * Tcb) << endl;
#endif
        curve.col(i + 1) = pose.transformationToCurveValue(observations[i].T_t_c().T() * Tcb);
    }
    CHECK(!curve.array().isNaN().any()) << "NaNs in curve values";
    // add 2 seconds on either end to allow the spline to slide during optimization
    times[0] = times[1] - 2.0 * timeOffsetPadding;
    times[kN + 1] = times[kN] + 2.0 * timeOffsetPadding;
    curve.col(0) = curve.col(1);
    curve.col(kN + 1) = curve.col(kN);
#if defined(DebugTest) && false
    // print first and last 10 rows in times and curve
    for (size_t i = 0; i < 10; ++i) {
        cout << format("times[{}] = {:.10f}", i, times[i]) << endl;
    }
    for (size_t i = kN + 2 - 10; i < kN + 2; ++i) {
        cout << format("times[{}] = {:.10f}", i, times[i]) << endl;
    }
    for (size_t i = 0; i < 10; ++i) {
        cout << format("curve[:, {}] = [{:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}]", i, curve(0, i),
                       curve(1, i), curve(2, i), curve(3, i), curve(4, i), curve(5, i))
             << endl;
    }
    for (size_t i = kN + 2 - 10; i < kN + 2; ++i) {
        cout << format("curve[{}] = [{:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}]", i, curve(0, i),
                       curve(1, i), curve(2, i), curve(3, i), curve(4, i), curve(5, i))
             << endl;
    }
#endif
    // make sure the rotation vector doesn't flip
    // NOTE by CC: adjust the angle, and select the min distance to previous rotation vector
    for (size_t i = 1; i < kN + 2; ++i) {
        Vector3d preRotationVector = curve.block<3, 1>(3, i - 1);
        Vector3d r = curve.block<3, 1>(3, i);
        double angle = r.norm();
        Vector3d axis = r / angle;
        Vector3d bestR = r;
        double bestDist = (bestR - preRotationVector).norm();

        for (int s = -3; s < 4; ++s) {
            Vector3d aa = axis * (angle + 2.0 * s * M_PI);
            double dist = (aa - preRotationVector).norm();
            if (dist < bestDist) {
                bestR = aa;
                bestDist = dist;
            }
        }

        curve.block<3, 1>(3, i) = bestR;
    }

    double seconds = times[kN + 1] - times[0];
    int knots = round(seconds * poseKnotsPerSecond);
    cout << format("Initializing a pose spline with {} knots ({} knots per second over {} seconds)", knots,
                   poseKnotsPerSecond, seconds)
         << endl;
    // int knots = round(seconds * cameraParams.frameRate / 3);
    // cout << format("Initializing a pose spline with {} knots ({} knots per second over {} seconds)", knots,
    //                knots / seconds, seconds)
    //      << endl;

    // note by CC: seems like the curve fitting. The first parameters is the timestamp, the second is the pose(position
    // + rotation, size = 6), third parameters is the knots number, and the last one don't understand yet
    pose.initPoseSplineSparse(times, curve, knots, 1e-4);
    return pose;
}

void RollingShutterCamera::findTimeShiftCameraImuPrior(const Imu& imu) {
    cout << SubSection("Estimating time shift camera to IMU");
    bsplines::BSplinePose poseSpline = initPoseSplineFromCamera(initSplineOrder, poseKnotsPerSecond, 0);

    // predict time shift prior
    vector<double> t;
    vector<double> measureGyroNorm;
    vector<double> predictGyroNorm;
    for (auto& m : imu.data) {
        double tk = m.stamp.toSec();
        if (poseSpline.tMin() < tk && tk < poseSpline.tMax()) {
            // get IMU measurements and spline from camera
            const Vector3d& measure = m.gyro;
            backend::EuclideanExpression predict(poseSpline.angularVelocityBodyFrame(tk));

            // calculate norm
            t.emplace_back(tk);
            measureGyroNorm.emplace_back(measure.norm());
            predictGyroNorm.emplace_back(predict.toEuclidean().norm());
        }
    }
    CHECK(!measureGyroNorm.empty() && !predictGyroNorm.empty())
        << "The time ranges of the camera and IMU do not overlap. "
           "Please make sure that your sensors are synchronized correctly.";

    // get the time shift, the correlation between angular velocity of IMU and camera
    vector<double> corr = correlate(predictGyroNorm, measureGyroNorm);
    auto itMax = max_element(corr.begin(), corr.end());
    int discreteShift = distance(corr.begin(), itMax) - measureGyroNorm.size() + 1;
    // calculate the average time difference between IMU measurement
    vector<double> times;
    transform(imu.data.begin(), imu.data.end(), back_inserter(times),
              [](const ImuMeasurement& m) { return m.stamp.toSec(); });
    double sum{0};
    for (size_t i = 1; i < times.size(); ++i) {
        sum += times[i] - times[i - 1];
    }
    double dT = sum / (times.size() - 1);
    timeshiftCameraToImuPrior = -discreteShift * dT;
    cout << format("Time shift camera to IMU (t_imu = t_cam + shift) = {:.10f}", timeshiftCameraToImuPrior) << endl;
}

// estimate IMU-Camera rotation prior
void RollingShutterCamera::findOrientationPriorCameraToImu(Imu& imu) {
    cout << SubSection("Estimate IMU-Camera Rotation Prior");

    // build the problem
    boost::shared_ptr<backend::OptimizationProblem> problem = boost::make_shared<backend::OptimizationProblem>();

    // add rotation R_BC as design variable
    boost::shared_ptr<backend::RotationQuaternion> qCB = boost::make_shared<backend::RotationQuaternion>(extrinsic.q());
    qCB->setActive(true);
    problem->addDesignVariable(qCB);

    // add the gyro bias b_g as design variable
    boost::shared_ptr<backend::EuclideanPoint> gyroBias = boost::make_shared<backend::EuclideanPoint>(Vector3d::Zero());
    gyroBias->setActive(true);
    problem->addDesignVariable(gyroBias);

    // initialize a pose spline using the camera poses
    bsplines::BSplinePose poseSpline = initPoseSplineFromCamera(initSplineOrder, poseKnotsPerSecond, 0.0);
    cout << format("pose spline time = [{:.10f}, {:.10f}]", poseSpline.tMin(), poseSpline.tMax()) << endl;
    for (auto& m : imu.data) {
        double tk = m.stamp.toSec();
        if (poseSpline.tMin() < tk && tk < poseSpline.tMax()) {
            backend::EuclideanExpression gyroPredict =
                qCB->toExpression() * backend::EuclideanExpression(poseSpline.angularVelocityBodyFrame(tk));
            problem->addErrorTerm(boost::make_shared<kalibr_errorterms::GyroscopeError>(m.gyro, m.gyroInvR, gyroPredict,
                                                                                        gyroBias->toExpression()));
#if defined(DebugTest) && false
            // cout << format("qIC = \n{}", qIC->toExpression().toRotationMatrix()) << endl;
            // cout << format("angularVelocity = {}", poseSpline.angularVelocityBodyFrame(tk).transpose()) << endl;
            cout << format("t= {:.10f}, measured = {}, predict = {}", tk, m.gyro.transpose(),
                           gyroPredict.toValue().transpose())
                 << endl;
            // cout << format("gyroInvR = \n{}", m.gyroInvR) << endl;
            // cout << format("gyroBias = {}", gyroBias->toExpression().toValue().transpose()) << endl;
#endif
        }
    }

    // check error terms number
    if (problem->numErrorTerms() == 0) {
        LOG(FATAL)
            << "failed to obtain orientation prior. Please make sure that your sensors are synchronized correctly.";
    }

    // optimization
    backend::Optimizer2Options options;
    options.verbose = true;
    options.linearSystemSolver = boost::make_shared<backend::BlockCholeskyLinearSystemSolver>();
    options.nThreads = 2;
    options.convergenceDeltaX = 1e-4;
    options.convergenceDeltaJ = 1;
    options.maxIterations = 50;
    backend::Optimizer2 optimizer(options);
    optimizer.setProblem(problem);
    optimizer.optimize();

    // overwrite the external rotation prior
    Matrix3d Rcb = qCB->toRotationMatrix().transpose();
    extrinsic = sm::kinematics::Transformation(sm::kinematics::rt2Transform(Rcb, extrinsic.t()));
    cout << format("Orientation prior camera-IMU found as (T_CB) = \n{}", Rcb) << endl;
    cout << format("estimated extrinsic = \n{}", extrinsic.T()) << endl;

    // estimate gravity in the world coordinate frame as the mean specific force
    cout << SubSection("Estimate Gravity");
    vector<Vector3d> aw;
    for (auto& m : imu.data) {
        double tk = m.stamp.toSec();
        if (poseSpline.tMin() < tk && tk < poseSpline.tMax()) {
            //- R_WB * a(t) = R_WC * R_CB * a(t) = -R_WB * R_CB * a(t), note R_WB = R_WC * R_CB before initialization
            aw.emplace_back(-poseSpline.orientation(tk) * Rcb * m.acc);
        }
    }
    Vector3d awMean = accumulate(aw.begin(), aw.end(), Vector3d(0, 0, 0),
                                 [](const Vector3d& sum, const Vector3d& v) { return sum + v; }) /
                      aw.size();
    gravity = awMean / awMean.norm() * 9.80655;
    cout << format("gravity was initialized to {} m/s^2", gravity.transpose()) << endl;

    // set the gyro bias prior
    imu.gyroBiasPrior = gyroBias->toExpression().toEuclidean();
    cout << format("gyro bias prior found as: {}", imu.gyroBiasPrior.transpose()) << endl;
}

void RollingShutterCamera::addDesignVariable(aslam::calibration::OptimizationProblem& problem, bool noTimeCalibration) {
    // add IMU-Camera extrinsics design variable
    TcbDesignVar = boost::make_shared<TransformationDesignVariable>(extrinsic, true, true);
    for (int i = 0; i < TcbDesignVar->numDesignVariables(); ++i) {
        problem.addDesignVariable(TcbDesignVar->designVariable(i), kCalibrationGroupId);
    }

    // add time delay design variable
    cameraTimeToImuDesignVar = boost::make_shared<backend::Scalar>(0.0);
    cameraTimeToImuDesignVar->setActive(!noTimeCalibration);
    problem.addDesignVariable(cameraTimeToImuDesignVar, kCalibrationGroupId);
}

void RollingShutterCamera::addErrorTerms(
    aslam::calibration::OptimizationProblem& problem,
    const boost::shared_ptr<aslam::splines::BSplinePoseDesignVariable>& poseDesignVariable,
    aslam::backend::TransformationExpression& Tcb, int blakeZisserCam, const double& timeOffsetPadding) {
    cout << Paragraph("Adding Camera Error Terms");

    // add camera design variable
    CameraDesignVariable cameraDeignVariable(boost::shared_ptr<CameraGeometry>(&geometry, sm::null_deleter()));
    problem.addDesignVariable(cameraDeignVariable.shutterDesignVariable(), kCalibrationGroupId);
    problem.addDesignVariable(cameraDeignVariable.projectionDesignVariable(), kCalibrationGroupId);
    problem.addDesignVariable(cameraDeignVariable.distortionDesignVariable(), kCalibrationGroupId);

    for (auto& v : observations) {
        // as we are applying an initial time shift outside the optimization, so we need to make sure that we don't add
        // data outside the spline definition
        backend::ScalarExpression frameTime =
            cameraTimeToImuDesignVar->toExpression() + v.time().toSec() + timeshiftCameraToImuPrior;
        double frameTimeScalar = frameTime.toScalar();
        if (frameTimeScalar <= poseDesignVariable->spline().tMin() ||
            frameTimeScalar >= poseDesignVariable->spline().tMax()) {
            continue;
        }

        // get the image and target points corresponding to the frame
        vector<cv::Point2f> imageCorners;
        vector<cv::Point3f> targetCorners;
        v.getCornersImageFrame(imageCorners);
        v.getCornersTargetFrame(targetCorners);

        // setup a frame to handle the distortion
        auto frame = boost::make_shared<Frame>();
        frame->setGeometry(boost::shared_ptr<CameraGeometry>(&geometry, sm::null_deleter()));
        frame->setTime(v.time());
        frames.emplace_back(frame);

        // add error terms fro each observed corner
        for (size_t n = 0; n < imageCorners.size(); ++n) {
            // keypoint time offset by line delay
            Vector2d keypoint(imageCorners[n].x, imageCorners[n].y);
            backend::ScalarExpression keypointTime = cameraDeignVariable.keypointTime(Time(frameTimeScalar), keypoint);

            // transformation from IMU to world T_WB
            TransformationExpression Twb =
                poseDesignVariable->transformationAtTime(keypointTime, timeOffsetPadding, timeOffsetPadding);
            TransformationExpression Tbw = Twb.inverse();
            // calibration target from world frame to camera frame
            TransformationExpression Tcw = Tcb * Tbw;

            // transform target point to camera frame
            Vector4d targetPoint(targetCorners[n].x, targetCorners[n].y, targetCorners[n].z, 1);
            // HomogeneousExpression TransformationExpression::operator*(const HomogeneousExpression& rhs) const
            HomogeneousExpression p = Tcw * HomogeneousExpression(targetPoint);

            // add key point frame
            size_t keypointId = frame->numKeypoints();
            KeyPoint k;
            k.setMeasurement(keypoint);
            k.setInverseMeasurementCovariance(Matrix2d::Identity() * cornerUncertainty * cornerUncertainty);
            frame->addKeypoint(k);

            // create reprojection error
            auto error = boost::make_shared<AdaptiveCovarianceReprojectionError>(
                frame.get(), keypointId, p, cameraDeignVariable, poseDesignVariable.get());
            // add blake-zisserman M-estimator
            if (blakeZisserCam > 0) {
                auto estimator = boost::make_shared<BlakeZissermanMEstimator>(blakeZisserCam);
                error->setMEstimatorPolicy(estimator);
            }

            reprojectionErrors.emplace_back(error);
            problem.addErrorTerm(error);
        }
    }

    cout << format("\tAdd {} camera errors terms", observations.size()) << endl;
}

sm::kinematics::Transformation RollingShutterCamera::getResultTransformationImuToCam() {
    return sm::kinematics::Transformation(TcbDesignVar->T());
}

double RollingShutterCamera::getResultTimeShift() {
    return cameraTimeToImuDesignVar->toScalar() + timeshiftCameraToImuPrior;
}

// detect observation(corners) from ros bag
void RollingShutterCamera::detectObservations(rosbag::View& view, const AprilTargetParameters& targetParams,
                                              const boost::filesystem::path& obsPath) {
    // setup target detector
    cameras::GridCalibrationTargetAprilgrid::AprilgridOptions gridOptions;
    gridOptions.showExtractionVideo = false;
    gridOptions.minTagsForValidObs = max(targetParams.rows, targetParams.cols) + 1;
    auto grid = boost::make_shared<aslam::cameras::GridCalibrationTargetAprilgrid>(
        targetParams.rows, targetParams.cols, targetParams.size, targetParams.spacing, gridOptions);
    cameras::GridDetector::GridDetectorOptions detectorOptions;
    detectorOptions.imageStepping = false;
    detectorOptions.plotCornerReprojection = false;
    detectorOptions.filterCornerOutliers = true;
    aslam::cameras::GridDetector detector =
        cameras::GridDetector(cameras::CameraGeometryBase::Ptr(&geometry, sm::null_deleter()), grid, detectorOptions);

    // extract corners from data set
    size_t index{0};
    cout << format("\tExtract corners for image {}/{}", index, view.size()) << flush;
    for (auto it = view.begin(); it != view.end(); ++it) {
        sensor_msgs::Image::ConstPtr v = it->instantiate<sensor_msgs::Image>();
        if (v) {
            cout << format("\r\tExtract corners for image {}/{}", index, view.size()) << flush;
            ++index;
            Time time(v->header.stamp.sec, v->header.stamp.nsec);
            cv_bridge::CvImageConstPtr cvPtr = cv_bridge::toCvShare(v, v->encoding);
            cameras::GridCalibrationTargetObservation obs;
            if (detector.findTarget(cvPtr->image, time, obs)) {
#if defined(DebugTest) && false
                cout << endl;
                cout << format("time = {}, {}", obs.time().sec, obs.time().nsec) << endl;
                // cout << format("obs.points = \n{}", obs.points()) << endl;
                cout << format("obs.T_t_c().T() = \n{}", obs.T_t_c().T()) << endl;
#endif

                // clear images
                obs.clearImage();

                observations.emplace_back(obs);
            }
        }
    }
    cout << endl;
    cout << format("Extracted corners for {}/{} images", observations.size(), view.size()) << endl;

    // save observations to file
    cout << format("serialize observations to file \"{}\"", obsPath.string()) << endl;
    sm::boost_serialization::save_xml(observations, obsPath);
}
