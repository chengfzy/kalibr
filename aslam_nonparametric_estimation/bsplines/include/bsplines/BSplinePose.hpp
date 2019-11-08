/**
 * @file   BSpline.hpp
 * @author Paul Furgale <paul.furgale@utoronto.ca>
 * @date   Fri Feb 11 13:51:57 2011
 *
 * @brief  A class to facilitate state estimation for vehicles in 3D
 *         space using B-splines.
 *
 *
 */

#ifndef _BSPLINE_POSE_HPP
#define _BSPLINE_POSE_HPP

#include <sm/kinematics/RotationalKinematics.hpp>
#include "BSpline.hpp"

namespace bsplines {

/**
 * @brief A class to facilitate state estimation for vehicles in 3D space using B-Splines.
 * The spline represents a pose with respect to some navigation frame \f$ \mathbf F_n \f$.
 *
 * Variable x = [p, r] = [p, phi]
 */
class BSplinePose : public BSpline {
  public:
    /**
     * @brief Create a spline of the specified order. The resulting B-spline will be a series of piecewise polynomials
     * of degree splineOrder - 1.
     *
     * @param splineOrder The order of the spline.
     */
    BSplinePose(int splineOrder, const sm::kinematics::RotationalKinematics::Ptr& rotationalKinematics);

    /**
     * A destructor
     */
    ~BSplinePose();

    Eigen::Matrix4d transformation(double tk) const;
    Eigen::Matrix4d transformationAndJacobian(double tk, Eigen::MatrixXd* J = nullptr,
                                              Eigen::VectorXi* coefficientIndices = nullptr) const;

    Eigen::Matrix4d inverseTransformationAndJacobian(double tk, Eigen::MatrixXd* J = nullptr,
                                                     Eigen::VectorXi* coefficientIndices = nullptr) const;
    Eigen::Matrix4d inverseTransformation(double tk) const;

    Eigen::Vector4d transformVectorAndJacobian(double tk, const Eigen::Vector4d& v, Eigen::MatrixXd* J = nullptr,
                                               Eigen::VectorXi* coefficientIndices = nullptr) const;
    /**
     * @brief Position at certain time. p = x[:3]
     * @param tk Time t
     * @return Position at certain time
     */
    Eigen::Vector3d position(double tk) const;

    /**
     * @brief Orientation at certain time, phi = x[3:], R = Exp(phi)
     * @param tk Time t
     * @return Orientation at certain time
     */
    Eigen::Matrix3d orientation(double tk) const;
    Eigen::Matrix3d orientationAndJacobian(double tk, Eigen::MatrixXd* J, Eigen::VectorXi* coefficientIndices) const;

    /**
     * @brief Inverse orientation at certain time, phi = x[3:], R^-1 = R^T = Exp(-phi)
     * @param tk Time t
     * @return Inverse orientation at certain time
     */
    Eigen::Matrix3d inverseOrientation(double tk) const;
    Eigen::Matrix3d inverseOrientationAndJacobian(double tk, Eigen::MatrixXd* J,
                                                  Eigen::VectorXi* coefficientIndices) const;

    /**
     * @brief Velocity v_W = dp / dt
     * @param tk
     * @return
     */
    Eigen::Vector3d linearVelocity(double tk) const;
    // v_B = R_BW *

    /**
     * @brief v_B = R_BW * v_W
     * @param tk
     * @return
     */
    Eigen::Vector3d linearVelocityBodyFrame(double tk) const;

    Eigen::Vector3d linearAcceleration(double tk) const;
    Eigen::Vector3d linearAccelerationBodyFrame(double tk) const;
    Eigen::Vector3d linearAccelerationAndJacobian(double tk, Eigen::MatrixXd* J,
                                                  Eigen::VectorXi* coefficientIndices) const;

    Eigen::Vector3d angularVelocity(double tk) const;
    Eigen::Vector3d angularVelocityBodyFrame(double tk) const;
    Eigen::Vector3d angularVelocityBodyFrameAndJacobian(double tk, Eigen::MatrixXd* J,
                                                        Eigen::VectorXi* coefficientIndices) const;

    Eigen::Vector3d angularVelocityAndJacobian(double tk, Eigen::MatrixXd* J,
                                               Eigen::VectorXi* coefficientIndices) const;

    // Fabio (not used and not tested)
    Eigen::Vector3d angularAccelerationBodyFrame(double tk) const;
    Eigen::Vector3d angularAccelerationAndJacobian(double tk, Eigen::MatrixXd* J,
                                                   Eigen::VectorXi* coefficientIndices) const;
    Eigen::Vector3d angularAccelerationBodyFrameAndJacobian(double tk, Eigen::MatrixXd* J,
                                                            Eigen::VectorXi* coefficientIndices) const;

    void initPoseSpline(double t0, double t1, const Eigen::Matrix4d& T_n_t0, const Eigen::Matrix4d& T_n_t);
    void initPoseSpline2(const Eigen::VectorXd& times, const Eigen::Matrix<double, 6, Eigen::Dynamic>& poses,
                         int numSegments, double lambda);
    void initPoseSpline3(const Eigen::VectorXd& times, const Eigen::Matrix<double, 6, Eigen::Dynamic>& poses,
                         int numSegments, double lambda);
    void initPoseSplineSparse(const Eigen::VectorXd& times, const Eigen::Matrix<double, 6, Eigen::Dynamic>& poses,
                              int numSegments, double lambda);
    void initPoseSplineSparseKnots(const Eigen::VectorXd& times, const Eigen::MatrixXd& interpolationPoints,
                                   const Eigen::VectorXd& knots, double lambda);

    void addPoseSegment(double tk, const Eigen::Matrix4d& T_n_tk);
    void addPoseSegment2(double tk, const Eigen::Matrix4d& T_n_tk, double lambda);

    Eigen::Matrix4d curveValueToTransformation(const Eigen::VectorXd& c) const;
    Eigen::VectorXd transformationToCurveValue(const Eigen::Matrix4d& T) const;

    sm::kinematics::RotationalKinematics::Ptr rotation() const;

    Eigen::Matrix4d curveValueToTransformationAndJacobian(const Eigen::VectorXd& c, Eigen::MatrixXd* J) const;

  private:
    sm::kinematics::RotationalKinematics::Ptr rotation_;
};

}  // namespace bsplines

#endif /* _BSPLINE_POSE_HPP */
