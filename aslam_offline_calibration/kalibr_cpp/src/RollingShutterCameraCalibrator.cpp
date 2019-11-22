#include "cc/RollingShutterCameraCalibrator.h"
#include <cv_bridge/cv_bridge.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <glog/logging.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <sensor_msgs/Image.h>
#include <aslam/backend/BlockCholeskyLinearSystemSolver.hpp>
#include <aslam/backend/DogLegTrustRegionPolicy.hpp>
#include "aslam/Keypoint.hpp"
#include "aslam/backend/BSplineMotionError.hpp"
#include "aslam/backend/HomogeneousPoint.hpp"
#include "aslam/backend/Optimizer2.hpp"
#include "aslam/backend/Optimizer2Options.hpp"
#include "cc/Heading.hpp"
#include "cc/KnotSequenceUpdateStrategy.h"
#include "cc/Util.h"
#include "sm/kinematics/RotationVector.hpp"

using namespace std;
using namespace cc;
using namespace Eigen;
using namespace rosbag;
using namespace fmt;
using namespace aslam;

RollingShutterCameraCalibrator::RollingShutterCameraCalibrator(const string& bagFile,
                                                               const CameraParameters& cameraParams,
                                                               const AprilTargetParameters& targetParams,
                                                               const Options& options, const ros::Time& startTime,
                                                               const ros::Time& endTime)
    : options(options),
      bagFile(bagFile),
      startTime(startTime),
      endTime(endTime),
      cameraParams(cameraParams),
      targetParams(targetParams) {}

// calibration
void RollingShutterCameraCalibrator::calibrate() {
    // init
    extractObservations();
    initIntrinsics();
    initExtrinsics();
    // set the motion model prior
    Matrix<double, 6, 6> W = Matrix<double, 6, 6>::Identity();
    W.block<3, 3>(0, 0) *= 1e-5;
    W.block<3, 3>(3, 3) *= 1e-2;

    cameraDeignVariable = boost::make_shared<CameraDesignVariable>(
        boost::shared_ptr<CameraGeometry>(&cameraGeometry, sm::null_deleter()));

    // build problem and solve
    auto poseSpline =
        initPoseBSpline(options.splineOrder, options.timeOffsetPadding, options.numberOfKnots, options.frameRate);
    buildProblem(poseSpline, W);
    solve();

    // continue with knot replacement
    if (options.adaptiveKnotPlacement) {
        KnotSequenceUpdateStrategy knotSequenceUpdateStrategy(options.frameRate);

        for (int i = 0; i < options.maxKnotPlacementIterations; ++i) {
            cout << Section(format("[{}] Adaptive Knot Placement", i));
            // generate the ne knots list
            vector<double> knots;
            bool requireUpdate =
                knotSequenceUpdateStrategy.generateKnots(reprojectionErrors, poseDesignVariable->spline(), knots);
            // if no new knots was generated, break
            if (!requireUpdate) {
                break;
            }

            // otherwise update the spline dv and rebuild the problem
            poseSpline =
                knotSequenceUpdateStrategy.getUpdatedSpline(poseDesignVariable->spline(), knots, options.splineOrder);
            buildProblem(poseSpline, W);
            solve();
        }
    }

    printResult();
}

// extract observations from ros bag image
void RollingShutterCameraCalibrator::extractObservations() {
    cout << SubSection("Extract Observations from ROS Bag");
    cout << format("\tDataset: {}", bagFile) << endl;
    cout << format("\tTopic: {}", cameraParams.topic) << endl;
    cout << format("\tTime: [{}, {}]", startTime, endTime) << endl;

    // init ros bag reader
    Bag bag(bagFile);
    View view(bag, TopicQuery(cameraParams.topic), ros::Time(startTime), ros::Time(endTime));
    cout << format("bag message size = {}, t = [{}, {}]", view.size(), view.getBeginTime(), view.getEndTime()) << endl;

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
    detector = cameras::GridDetector(cameras::CameraGeometryBase::Ptr(&cameraGeometry, sm::null_deleter()), grid,
                                     detectorOptions);

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
            if (detector.findTargetNoTransformation(cvPtr->image, time, obs)) {
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
}

// init intrinsics
void RollingShutterCameraCalibrator::initIntrinsics() {
    cout << SubSection("Initialize Intrinsics");
    int sensorRows = observations.begin()->imRows();
    cameraGeometry.shutter().setParameters(Matrix<double, 1, 1>(1.0 / options.frameRate / sensorRows));
    cameraGeometry.initializeIntrinsics(observations);

    Matrix<double, 8, 1> params;
    params << cameraParams.f, cameraParams.c, cameraParams.d;
    cameraGeometry.setParameters(params, true, true, false);

    // print parameters after init
    MatrixXd initParams;
    cameraGeometry.getParameters(initParams, true, true, true);
    cout << format("init camera parameters: {}", initParams.transpose()) << endl;
}

// init extrinsics
void RollingShutterCameraCalibrator::initExtrinsics() {
    cout << SubSection("Initialize Extrinsics");
    for (size_t i = 0; i < observations.size(); ++i) {
        sm::kinematics::Transformation T_TC;
        if (cameraGeometry.estimateTransformation(observations[i], T_TC)) {
            observations[i].set_T_t_c(T_TC);
        } else {
            LOG(ERROR) << format("Could not estimate T_TC for observation at [{}]", i);
        }
    }
}

// init pose B-Spline
bsplines::BSplinePose RollingShutterCameraCalibrator::initPoseBSpline(int splineOrder, const double& timeOffsetPadding,
                                                                      int numberOfKnots, const double& frameRate) {
    cout << SubSection("Initialize Pose B-Spline");
    bsplines::BSplinePose poseSpline(splineOrder, boost::make_shared<sm::kinematics::RotationVector>());

    // get the checkerboard times
    const size_t kN = observations.size();
    VectorXd times(kN + 2);
    MatrixXd curve(6, kN + 2);
    for (size_t i = 0; i < kN; ++i) {
        times[i + 1] = observations[i].time().toSec();
        curve.col(i + 1) = poseSpline.transformationToCurveValue(observations[i].T_t_c().T());
    }
    CHECK(!curve.array().isNaN().any()) << "NaNs in curve values";
    // add 2 seconds on either end to allow the spline to slide during optimization
    times[0] = times[1] - 2.0 * timeOffsetPadding;
    times[kN + 1] = times[kN] + 2.0 * timeOffsetPadding;
    curve.col(0) = curve.col(1);
    curve.col(kN + 1) = curve.col(kN);

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
    int knots = numberOfKnots;
    if (knots == -1) {
        knots = round(seconds * frameRate / 3);
    }
    cout << format("Initializing a pose spline with {} knots ({} knots per second over {} seconds)", knots,
                   knots / seconds, seconds)
         << endl;

    // note by CC: seems like the curve fitting. The first parameters is the timestamp, the second is the pose(position
    // + rotation, size = 6), third parameters is the knots number, and the last one don't understand yet
    poseSpline.initPoseSplineSparse(times, curve, knots, 1e-4);

    return poseSpline;
}

// build the optimization problem
void RollingShutterCameraCalibrator::buildProblem(const bsplines::BSplinePose& poseSpline,
                                                  const Eigen::Matrix<double, 6, 6>& W) {
    cout << SubSection("Build Optimization Problem");

    problem = boost::make_shared<aslam::calibration::OptimizationProblem>();

    // get the design variable representation of the pose spline and add them to the problem
    poseDesignVariable = boost::make_shared<splines::BSplinePoseDesignVariable>(poseSpline);
    for (size_t i = 0; i < poseDesignVariable->numDesignVariables(); ++i) {
        auto v = calibration::OptimizationProblem::DesignVariableSP(poseDesignVariable->designVariable(i),
                                                                    sm::null_deleter());
        v->setActive(true);
        problem->addDesignVariable(v, kCalibrationGroupId);
    }

    // get the keypoint index of target corner
    vector<unsigned int> keypointIndex;
    vector<cv::Point3f> targets;  // targets, all 3D points for april target
    for (auto& v : observations) {
        if (static_cast<int>(v.getCornersIdx(keypointIndex)) == targetParams.keyPointNumber()) {
            v.getCornersTargetFrame(targets);
            break;
        }
    }
    // check the index and targets number is the maximum
    if (static_cast<int>(keypointIndex.size()) != targetParams.keyPointNumber()) {
        LOG(FATAL) << format("the maximum number of keypoint {} don't match the target keypoint size {}",
                             keypointIndex.size(), targetParams.keyPointNumber());
    }

    // add target point to problem
    std::unordered_map<unsigned int, backend::HomogeneousExpression> targetExpression;
    for (size_t i = 0; i < targets.size(); ++i) {
        Vector4d targetPoint(targets[i].x, targets[i].y, targets[i].z, 1);
        auto targetDesignVariable = boost::make_shared<backend::HomogeneousPoint>(targetPoint);
        targetDesignVariable->setActive(false);
        targetExpression[keypointIndex[i]] = targetDesignVariable->toExpression();
        problem->addDesignVariable(targetDesignVariable, kCalibrationGroupId);
    }

    // add camera design variable and activate
    problem->addDesignVariable(cameraDeignVariable->shutterDesignVariable(), kCalibrationGroupId);
    problem->addDesignVariable(cameraDeignVariable->projectionDesignVariable(), kCalibrationGroupId);
    problem->addDesignVariable(cameraDeignVariable->distortionDesignVariable(), kCalibrationGroupId);
    cameraDeignVariable->setActive(false, false, true);

    // regularization term(motion prior
    auto motionError = boost::make_shared<backend::BSplineMotionError<splines::BSplinePoseDesignVariable>>(
        poseDesignVariable.get(), W);
    problem->addErrorTerm(motionError);

    // add observation to problem
    frames.clear();
    reprojectionErrors.clear();
    for (auto& v : observations) {
        if (v.hasSuccessfulObservation()) {
            // add a frame
            boost::shared_ptr<Frame> frame = boost::make_shared<Frame>();
            frame->setGeometry(boost::shared_ptr<CameraGeometry>(&cameraGeometry, sm::null_deleter()));
            frame->setTime(v.time());
            frames.emplace_back(frame);

            // get the keypoint index of current image(observation)
            vector<unsigned int> keypointIds;
            v.getCornersIdx(keypointIds);

            // add an error term for every observed corner
            vector<cv::Point2f> corners;
            v.getCornersImageFrame(corners);
            for (size_t i = 0; i < corners.size(); ++i) {
                Vector2d point(corners[i].x, corners[i].y);
                // keypoint time offset by line delay as expression type
                backend::ScalarExpression keypointTime = cameraDeignVariable->keypointTime(frame->time(), point);

                // from target to world transformation
                backend::TransformationExpression T_WT = poseDesignVariable->transformationAtTime(
                    keypointTime, options.timeOffsetConstantSparsityPattern, options.timeOffsetConstantSparsityPattern);
                backend::TransformationExpression T_TW = T_WT.inverse();

                // transform target point to camera frame
                backend::HomogeneousExpression point_T = T_TW * targetExpression[keypointIds[i]];

                // create the keypoint
                size_t keypointId = frame->numKeypoints();
                Keypoint<2> keypoint;
                keypoint.setMeasurement(point);
                keypoint.setInverseMeasurementCovariance(Matrix2d::Identity() * options.inverseFeatureCovariance);
                keypoint.setLandmarkId(LandmarkId(keypointId));
                frame->addKeypoint(keypoint);

                // create reprojection error
                auto reprojectionError = boost::make_shared<AdaptiveCovarianceReprojectionError>(
                    frame.get(), keypointId, point_T, *cameraDeignVariable, poseDesignVariable.get());
                reprojectionErrors.emplace_back(reprojectionError);
                problem->addErrorTerm(reprojectionError);
            }
        }
    }
}

// solve problem
void RollingShutterCameraCalibrator::solve() {
    cout << Section("Solve Problem");
    printResult();

    // optimizer
    backend::Optimizer2Options optimizerOptions;
    optimizerOptions.verbose = true;
    optimizerOptions.linearSystemSolver = boost::make_shared<backend::BlockCholeskyLinearSystemSolver>();
    optimizerOptions.doSchurComplement = true;
    optimizerOptions.maxIterations = options.maxIterationNumber;
    optimizerOptions.convergenceDeltaX = options.deltaX;
    optimizerOptions.convergenceDeltaJ = options.deltaJ;
    optimizerOptions.trustRegionPolicy = boost::make_shared<backend::DogLegTrustRegionPolicy>();

    // create optimizer and solve
    backend::Optimizer2 optimizer(optimizerOptions);
    optimizer.setProblem(problem);
    optimizer.optimize();
}

// print results
void RollingShutterCameraCalibrator::printResult() {
    cout << Paragraph("Results");
    cout << format("LineDelay: {}", cameraDeignVariable->shutterDesignVariable()->value().lineDelay()) << endl;
    MatrixXd projection, distortion;
    cameraDeignVariable->projectionDesignVariable()->value().getParameters(projection);
    cameraDeignVariable->distortionDesignVariable()->value().getParameters(distortion);
    cout << format("Projection: {}", projection.transpose()) << endl;
    cout << format("Distortion: {}", distortion.transpose()) << endl;
}
