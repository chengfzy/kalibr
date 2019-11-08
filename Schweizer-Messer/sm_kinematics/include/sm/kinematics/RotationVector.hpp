#ifndef _ROTATIONVECTOR_H_
#define _ROTATIONVECTOR_H_

#include "RotationalKinematics.hpp"

namespace sm {
namespace kinematics {

class RotationVector : public RotationalKinematics {
  public:
    ~RotationVector() override = default;

    Eigen::Matrix3d parametersToRotationMatrix(const Eigen::Vector3d& parameters,
                                               Eigen::Matrix3d* S = nullptr) const override;
    Eigen::Vector3d rotationMatrixToParameters(const Eigen::Matrix3d& rotationMatrix) const override;

    // to left Jacobian, note kalibr using JPL conversion instead of Hamilton ones, so there exist much different
    Eigen::Matrix3d parametersToSMatrix(const Eigen::Vector3d& parameters) const override;

    Eigen::Vector3d angularVelocityAndJacobian(const Eigen::Vector3d& p, const Eigen::Vector3d& pdot,
                                               Eigen::Matrix<double, 3, 6>* Jacobian) const override;
    // to Inverse of left Jacobian
    Eigen::Matrix3d parametersToInverseSMatrix(const Eigen::Vector3d& parameters) const;
};

}  // namespace kinematics
}  // namespace sm

#endif /* _ROTATIONVECTOR_H_ */
