// Bring in gtest
#include <gtest/gtest.h>
#include <aslam/cameras/DepthCameraGeometry.hpp>
#include <aslam/cameras/test/CameraGeometryTestHarness.hpp>
#include <sm/eigen/gtest.hpp>

TEST(AslamCamerasTestSuite, testDepthCameraGeometry) {
    using namespace aslam::cameras;
    CameraGeometryTestHarness<DepthCameraGeometry> harness(1e-2);
    SCOPED_TRACE("");
    harness.testAll();
}
