#include "cc/CameraParameters.h"
#include <fmt/format.h>

using namespace std;
using namespace fmt;

namespace cc {

std::ostream& operator<<(std::ostream& os, const cc::CameraParameters& cameraParams) {
    os << format("\tCamera model: {}", cameraParams.model) << endl;
    os << format("\tFocal length: [{}, {}]", cameraParams.f[0], cameraParams.f[1]) << endl;
    os << format("\tPrincipal point: [{}, {}]", cameraParams.c[0], cameraParams.c[1]) << endl;
    os << format("\tDistortion model: {}", cameraParams.distortModel) << endl;
    os << format("\tDistortion coefficients: [{}, {}, {}, {}]", cameraParams.d[0], cameraParams.d[1], cameraParams.d[2],
                 cameraParams.d[3])
       << endl;
    os << format("\tROS topic: {}", cameraParams.topic) << endl;
    os << format("\tResolution: [{}, {}]", cameraParams.resolution[0], cameraParams.resolution[1]) << endl;
    return os;
}

}  // namespace cc