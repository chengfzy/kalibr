#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include "cc/Calibrator.h"
#include "cc/Heading.hpp"

using namespace std;
using namespace Eigen;
using namespace fmt;
using namespace cc;

DEFINE_string(bagFile, "../data/data.bag", "bag file");

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // load and set parameters
    cout << Section("Load and Set Parameters");
    // ROS
    // double startTime{1568798143.76}, endTime{1568798153.76};
    double startTime{1568797000.69}, endTime{1568798317.18};
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
    cameraParameters.topic = "/camIds/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(1395.55335893991, 1396.18082713661);
    cameraParameters.c = Vector2d(996.891838559857, 575.901405887637);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d = Vector4d(-0.160213305777107, 0.0873377403779236, 0.000591627579404801, -1.01770731705533e-05);
    cameraParameters.resolution = Vector2i(1936, 1216);
    cout << cameraParameters << endl;

    // init
    cout << Section("Initialization");
    Imu imu(FLAGS_bagFile, imuParameters, ros::Time(startTime), ros::Time(endTime));
    Camera camera(FLAGS_bagFile, cameraParameters, targetParameters, ros::Time(startTime), ros::Time(endTime));
    Calibrator calibrator(camera, imu);
    calibrator.buildProblem();

    cout << Section("Before Optimization");
    calibrator.printErrorStatistics();

    cout << Section("Optimizing...");
    calibrator.optimize(30, false);

    cout << Section("After Optimization");
    calibrator.printErrorStatistics();
    calibrator.printResult();

    cout << Section("Finish IMU-Camera Calibration");

    google::ShutdownGoogleLogging();
    google::ShutDownCommandLineFlags();
    return 0;
}