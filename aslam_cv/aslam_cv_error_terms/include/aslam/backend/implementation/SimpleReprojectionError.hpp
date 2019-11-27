namespace aslam {
namespace backend {

template <typename F>
SimpleReprojectionError<F>::SimpleReprojectionError() {}

template <typename F>
SimpleReprojectionError<F>::SimpleReprojectionError(const frame_t* frame, int keypointIndex,
                                                    const HomogeneousExpression& point)
    : _point(point) {
    SM_ASSERT_TRUE(aslam::Exception, frame, "The frame must not be null");
    _y = frame->keypoint(keypointIndex).y();
    _geometry = &frame->geometry();
    parent_t::setInvR(frame->keypoint(keypointIndex).invR());
    JacobianContainer::set_t dvs;
    point.getDesignVariables(dvs);
    parent_t::setDesignVariablesIterator(dvs.begin(), dvs.end());
}

template <typename F>
SimpleReprojectionError<F>::SimpleReprojectionError(const measurement_t& y, const inverse_covariance_t& invR,
                                                    const HomogeneousExpression& point,
                                                    const camera_geometry_t& geometry)
    : _y(y), _point(point) {
    _geometry = &geometry;
    parent_t::setInvR(invR);
    JacobianContainer::set_t dvs;
    point.getDesignVariables(dvs);
    parent_t::setDesignVariablesIterator(dvs.begin(), dvs.end());
}

template <typename F>
double SimpleReprojectionError<F>::evaluateErrorImplementation() {
    // the reprojection process
    Eigen::Vector4d p = _point.toHomogeneous();
    measurement_t hat_y;
    _geometry->homogeneousToKeypoint(p, hat_y);
    parent_t::setError(_y - hat_y);
    return parent_t::error().dot(parent_t::invR() * parent_t::error());
}

template <typename F>
void SimpleReprojectionError<F>::evaluateJacobiansImplementation(aslam::backend::JacobianContainer& _jacobians) const {
    Eigen::Vector4d p = _point.toHomogeneous();
    typename camera_geometry_t::jacobian_homogeneous_t J;
    measurement_t hat_y;

    // reproject point from world(target) frame to camera frame, and calculate the Jacobians
    _geometry->homogeneousToKeypoint(p, hat_y, J);

    // call HomogeneousExpressionNodeMultiply::evaluateJacobiansImplementation(JacobianContainer& outJacobians, const
    // Eigen::MatrixXd& applyChainRule) const, but why -J?
    _point.evaluateJacobians(_jacobians, -J);
}

}  // namespace backend
}  // namespace aslam
