#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include "cc/Heading.hpp"
#include "cc/ImuCameraCalibrator.h"

using namespace std;
using namespace Eigen;
using namespace fmt;
using namespace cc;

DEFINE_string(bagFile, "../data/data.bag", "bag file");

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    cout << Section("", false) << Section("IMU Camera Calibration", false) << Section("", false);

    // load and set parameters
    cout << Section("Load and Set Parameters");
    // ROS
    // double startTime{1568798143.76}, endTime{1568798153.76};
    // double startTime{1568797000.69}, endTime{1568798317.18};
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
    // camera parameters: SensingJAX52202
    cout << Paragraph("Camera Parameters");
    CameraParameters cameraParameters;
    cameraParameters.topic = "/camNormal/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(823.906785813718, 824.704794976756);
    cameraParameters.c = Vector2d(648.608565867332, 314.810494971340);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d[0] = -0.312098601430490;
    cameraParameters.d[1] = 0.0928470270344407;
    cameraParameters.d[2] = -1.93958495467811e-05;
    cameraParameters.d[3] = -0.000132104569851275;
    cameraParameters.d[4] = 0.0;
    cameraParameters.resolution = Vector2i(1280, 720);
    cameraParameters.lineDelay = 3.66e-5;
    cout << cameraParameters << endl;

    // init
    cout << Section("Initialization");
    Imu imu(FLAGS_bagFile, imuParameters);
    Camera camera(FLAGS_bagFile, cameraParameters, targetParameters);
    ImuCameraCalibrator calibrator(camera, imu);
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