#pragma once
#include <Eigen/Core>
#include <string>

namespace cc {

class AprilTargetParameters {
  public:
    friend std::ostream& operator<<(std::ostream& os, const AprilTargetParameters& targetParameters);

  public:
    std::string type;  // target type
    std::size_t cols;  // target cols
    std::size_t rows;  // target rows
    double size;       // size, m
    double spacing;    // ratio of space between tags to tag size
};

}  // namespace cc