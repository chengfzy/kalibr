#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include "cc/Heading.hpp"
#include "cc/RollingShutterCameraCalibrator.h"

using namespace std;
using namespace Eigen;
using namespace cc;

DEFINE_string(bagFile, "../data/data.bag", "bag file");
DEFINE_int32(frameRate, 30, "Approximate framerate of the camera");
DEFINE_double(featureVariance, 1.0, "Estimated variance of the feature detector");
DEFINE_int32(maxIter, 30, "Max iterations");

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    cout << Section("", false) << Section("Rolling Shutter Camera Calibration", false) << Section("", false);

    // load and set parameters
    cout << Section("Load and Set Parameters");
    // target parameters
    cout << Paragraph("Target Parameters");
    AprilTargetParameters targetParameters;
    targetParameters.type = "AprilGrid";
    targetParameters.cols = 7;
    targetParameters.rows = 7;
    targetParameters.size = 0.09;
    targetParameters.spacing = 0.3;
    cout << targetParameters << endl;
    // camera parameters: Sensing Camera JAX52102
    cout << Paragraph("Camera Parameters");
    CameraParameters cameraParameters;
    cameraParameters.topic = "/camSensing/image_raw";
    cameraParameters.model = "pinhole";
    cameraParameters.f = Vector2d(826.999140224398, 827.942754081805);
    cameraParameters.c = Vector2d(682.225164651297, 311.932742151825);
    cameraParameters.distortModel = "radtan";
    cameraParameters.d[0] = -0.314814115131056;
    cameraParameters.d[1] = 0.0945857861790417;
    cameraParameters.d[2] = 0.000176515325718172;
    cameraParameters.d[3] = 0.000678933511638204;
    cameraParameters.d[4] = 0.0;
    cameraParameters.resolution = Vector2i(1936, 1216);
    cout << cameraParameters << endl;

    // calibration
    RollingShutterCameraCalibrator::Options options;
    options.inverseFeatureCovariance = 1. / FLAGS_featureVariance;
    options.maxIterationNumber = FLAGS_maxIter;
    options.frameRate = FLAGS_frameRate;
    options.timeOffsetConstantSparsityPattern = 0.08;
    options.timeOffsetPadding = 0.5;
    RollingShutterCameraCalibrator calibrator(FLAGS_bagFile, cameraParameters, targetParameters, options);
    calibrator.calibrate();

    google::ShutdownGoogleLogging();
    google::ShutDownCommandLineFlags();
    return 0;
}