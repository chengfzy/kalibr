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

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // print heading
    cout << Section("", false) << Section("IMU Rolling Shutter Camera Calibration", false) << Section("", false);

    // load and set parameters
    cout << Section("Load and Set Parameters");
    // ROS
    // double startTime{1568798143.76}, endTime{1568798153.76};
    // IMU parameters
    cout << Paragraph("IMU Parameters");
    ImuParameters imuParameters;
    imuParameters.topic = "/imuSanChi";
    imuParameters.updateRate = 100;
    imuParameters.accNoiseDensity = 0.00447;
    imuParameters.accRandomWalk = 7.071e-5;
    imuParameters.gyroNoiseDensity = 0.063246;
    imuParameters.gyroRandomWalk = 0.001;
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
    cameraParameters.topic = "/camNormal/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(876.338148402960, 876.461971047948);
    cameraParameters.c = Vector2d(633.483570573810, 332.245456171106);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d = Vector4d(-0.324894923827232, 0.0946332635762246, -0.000641748136769180, -0.00204087263576566);
    cameraParameters.resolution = Vector2i(1280, 720);
    cameraParameters.lineDelay = 3.65e-5;
    cout << cameraParameters << endl;

    // init
    cout << Section("Initialization");
    Imu imu(FLAGS_bagFile, imuParameters);
    RollingShutterCamera camera(FLAGS_bagFile, cameraParameters, targetParameters);
    ImuRollingShutterCameraCalibrator calibrator(camera, imu);
    calibrator.buildProblem();

    cout << Section("Before Optimization");
    calibrator.printErrorStatistics();

    cout << Section("Optimizing...");
    calibrator.optimize(30, false);

    cout << Section("After Optimization");
    calibrator.printErrorStatistics();
    calibrator.printResult();

    cout << Section("Finish IMU-Rolling Shutter Camera Calibration");

    google::ShutdownGoogleLogging();
    google::ShutDownCommandLineFlags();
    return 0;
}