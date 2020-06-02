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
    // imuParameters.topic = "/imuSbg";
    // imuParameters.updateRate = 100;
    // imuParameters.accNoiseDensity = 0.00568;
    // imuParameters.accRandomWalk = 8.97e-5;
    // imuParameters.gyroNoiseDensity = 0.000488;
    // imuParameters.gyroRandomWalk = 3.19e-5;
    // the parameters of SparkFun and SanChi consider the same
    imuParameters.topic = "/imuSparkFun0";
    imuParameters.updateRate = 100;
    imuParameters.accNoiseDensity = 0.00447;
    imuParameters.accRandomWalk = 7.071e-5;
    imuParameters.gyroNoiseDensity = 0.063246;
    imuParameters.gyroRandomWalk = 0.001;
    // below parameter don't seem reasonable, ref: https://confluence.ygomi.com:8443/display/RF/SBG+IMU
    // imuParameters.topic = "/imuSparkFun0";
    // imuParameters.updateRate = 100;
    // imuParameters.accNoiseDensity = 0.0562;
    // imuParameters.accRandomWalk = 0.0049;
    // imuParameters.gyroNoiseDensity = 0.0050;
    // imuParameters.gyroRandomWalk = 1.1569e-5;
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
    // camera parameters
    cout << Paragraph("Camera Parameters");
    CameraParameters cameraParameters;
    // Sensing JAX52202
    // cameraParameters.topic = "/camSensing/image_raw";
    // cameraParameters.model = "pinhole";
    // cameraParameters.f = Vector2d(826.999140224398, 827.942754081805);
    // cameraParameters.c = Vector2d(682.225164651297, 311.932742151825);
    // cameraParameters.distortModel = "radtan";
    // cameraParameters.d[0] = -0.314814115131056;
    // cameraParameters.d[1] = 0.0945857861790417;
    // cameraParameters.d[2] = 0.000176515325718172;
    // cameraParameters.d[3] = 0.000678933511638204;
    // cameraParameters.d[4] = 0;
    // cameraParameters.lineDelay = 3.68e-5;
    // Sensing JAX80802
    cameraParameters.topic = "/camNormal1/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(829.708551330165, 828.779399078419);
    cameraParameters.c = Vector2d(652.639155240132, 348.883447318945);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d[0] = -0.345600280288694;
    cameraParameters.d[1] = 0.171489212385776;
    cameraParameters.d[2] = -3.44338669628556e-05;
    cameraParameters.d[3] = -1.49180538008743e-05;
    cameraParameters.d[4] = -0.0525549437324129;
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