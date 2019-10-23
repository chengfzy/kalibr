#include "cc/ImuParameters.h"
#include <fmt/format.h>

using namespace std;
using namespace fmt;

namespace cc {

std::ostream& operator<<(std::ostream& os, const cc::ImuParameters& imuParams) {
    os << format("\tUpdate rate: {:.1f}", imuParams.updateRate) << endl;
    os << "\tAccelerometer" << endl;
    os << format("\t\tNoise density: {}", imuParams.accNoiseDensity) << endl;
    os << format("\t\tNoise density (discrete): {}", imuParams.getAccNoiseDiscrete()) << endl;
    os << format("\t\tRandom walk: {}", imuParams.accRandomWalk) << endl;
    os << "\tGyroscope" << endl;
    os << format("\t\tNoise density: {}", imuParams.gyroNoiseDensity) << endl;
    os << format("\t\tNoise density (discrete): {}", imuParams.getGyroNoiseDiscrete()) << endl;
    os << format("\t\tRandom walk: {}", imuParams.gyroRandomWalk) << endl;
    os << format("\tROS topic: {}", imuParams.topic) << endl;
    return os;
}

}  // namespace cc