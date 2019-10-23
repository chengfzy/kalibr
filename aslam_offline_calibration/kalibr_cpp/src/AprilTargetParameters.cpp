#include "cc/AprilTargetParameters.h"
#include <fmt/format.h>

using namespace std;
using namespace cc;
using namespace fmt;

namespace cc {

std::ostream& operator<<(std::ostream& os, const AprilTargetParameters& targetParameters) {
    os << format("\tType: {}", targetParameters.type) << endl;
    os << "\tTags" << endl;
    os << format("\t\tRows: {}", targetParameters.rows) << endl;
    os << format("\t\tCols: {}", targetParameters.cols) << endl;
    os << format("\t\tSize: {} m", targetParameters.size) << endl;
    os << format("\t\tSpacing: {} m", targetParameters.spacing) << endl;

    return os;
}

}  // namespace cc