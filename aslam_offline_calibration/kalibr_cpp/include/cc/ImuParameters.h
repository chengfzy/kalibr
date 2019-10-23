#pragma once
#include <Eigen/Core>
#include <string>

namespace cc {

class ImuParameters {
  public:
    friend std::ostream& operator<<(std::ostream& os, const ImuParameters& imuParams);

    double getAccNoiseDiscrete() const { return accNoiseDensity / std::sqrt(1.0 / updateRate); }

    double getGyroNoiseDiscrete() const { return gyroNoiseDensity / std::sqrt(1.0 / updateRate); }

  public:
    std::string topic;        // ros topic
    double updateRate;        // update rate(Hz)
    double accNoiseDensity;   // accelerometer noise density in continuous time
    double accRandomWalk;     // accelerometer random walk
    double gyroNoiseDensity;  // gyroscope noise density in continuous time
    double gyroRandomWalk;    // gyroscope random walk
};

}  // namespace cc