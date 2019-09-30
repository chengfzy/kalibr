// Bring in gtest
#include <gtest/gtest.h>
#include <aslam/cameras.hpp>
#include <aslam/cameras/test/CameraGeometryTestHarness.hpp>
#include <sm/eigen/gtest.hpp>

TEST(AslamCamerasTestSuite, testDoubleSphereCameraGeometry) {
    using namespace aslam::cameras;
    CameraGeometryTestHarness<DoubleSphereCameraGeometry> harness(0.1);

    SCOPED_TRACE("double sphere camera");
    harness.testAll();
}
