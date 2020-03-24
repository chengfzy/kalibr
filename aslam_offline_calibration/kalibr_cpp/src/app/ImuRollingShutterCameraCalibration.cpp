#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include "cc/Heading.hpp"
#include "cc/ImuRollingShutterCameraCalibrator.h"

using namespace std;
using namespace Eigen;
using namespace fmt;
using namespace cc;

DEFINE_string(bagFile, "../data/data.bag", "bag file");
DEFINE_int32(poseKnotsPerSecond, 100, "B-Spline pose knots number per second");

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // print heading
    cout << Section("", false) << Section("IMU Rolling Shutter Camera Calibration", false) << Section("", false);
    cout << format("bag file: \"{}\"", FLAGS_bagFile) << endl;
    cout << format("poseKnotsPerSecond : {}", FLAGS_poseKnotsPerSecond) << endl;

    // load and set parameters
    cout << Section("Load and Set Parameters");
    // ROS
    // double startTime{1568798143.76}, endTime{1568798153.76};
    // IMU parameters
    cout << Paragraph("IMU Parameters");
    ImuParameters imuParameters;
    imuParameters.topic = "/imuSbg";
    imuParameters.updateRate = 100;
    imuParameters.accNoiseDensity = 0.00568;
    imuParameters.accRandomWalk = 8.97e-5;
    imuParameters.gyroNoiseDensity = 0.000488;
    imuParameters.gyroRandomWalk = 3.19e-5;
    cout << imuParameters << endl;
    // target parameters
    cout << Paragraph("Target Parameters");
    AprilTargetParameters targetParameters;
    targetParameters.type = "AprilGrid";
    targetParameters.cols = 7;
    targetParameters.rows = 7;
    targetParameters.size = 0.09;
    targetParameters.spacing = 0.3;
    cout << targetParameters << endl;
    // camera parameters: SensingJAX52202
    cout << Paragraph("Camera Parameters");
    CameraParameters cameraParameters;
    cameraParameters.topic = "/camSensing/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(826.999140224398, 827.942754081805);
    cameraParameters.c = Vector2d(682.225164651297, 311.932742151825);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d = Vector4d(-0.314814115131056, 0.0945857861790417, 0.000176515325718172, 0.000678933511638204);
    cameraParameters.lineDelay = 3.68e-5;
    cameraParameters.resolution = Vector2i(1936, 1216);
    cameraParameters.frameRate = 30;
    cout << cameraParameters << endl;

    // init
    cout << Section("Initialization");
    Imu imu(FLAGS_bagFile, imuParameters);
    RollingShutterCamera camera(FLAGS_bagFile, cameraParameters, targetParameters);
    // Note, set extrinsics of TCB,  and rotation using JPL expression
    // Vector4d qCB(0.00000, 0.71541, 0.69864, 0.00951);  // [qx, qy, qz, qw]
    // Vector3d pCB(0, 0.03, -0.03);
    // camera.extrinsic = sm::kinematics::Transformation(qCB, pCB);
    // cout << "Initial Tcb, Transformation T_cam0_imu0 (imu0 to cam0):" << endl;
    // cout << camera.extrinsic.T() << endl;
    // calibrator
    ImuRollingShutterCameraCalibrator::Options options;
    options.poseKnotsPerSecond = FLAGS_poseKnotsPerSecond;
    options.splineOrder = 6;
    options.timeOffsetConstantSparsityPattern = 0.08;
    options.timeOffsetPadding = 0.5;
    ImuRollingShutterCameraCalibrator calibrator(camera, imu);
    calibrator.options = options;
    calibrator.camera.initSplineOrder = options.splineOrder;            // set spline order
    calibrator.camera.poseKnotsPerSecond = options.poseKnotsPerSecond;  // set spline order
    calibrator.buildProblem();

    cout << Section("Before Optimization");
    calibrator.printErrorStatistics();

    cout << Section("Optimizing");
    calibrator.optimize(30, false);

    cout << Section("After Optimization");
    calibrator.printErrorStatistics();
    calibrator.printResult();

    cout << Section("Finish IMU-Rolling Shutter Camera Calibration");

    google::ShutdownGoogleLogging();
    google::ShutDownCommandLineFlags();
    return 0;
}