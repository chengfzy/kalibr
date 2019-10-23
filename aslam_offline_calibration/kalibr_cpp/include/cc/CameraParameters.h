#pragma once
#include <Eigen/Core>
#include <string>

namespace cc {

class CameraParameters {
  public:
    friend std::ostream& operator<<(std::ostream& os, const CameraParameters& cameraParams);

  public:
    std::string topic;           // ros topic
    std::string model;           // camera model
    std::string distortModel;    // distortion model
    Eigen::Vector2d f;           // focal length
    Eigen::Vector2d c;           // principal point
    Eigen::Vector4d d;           // distortion coefficients
    Eigen::Vector2i resolution;  // camera resolution
};

}  // namespace cc