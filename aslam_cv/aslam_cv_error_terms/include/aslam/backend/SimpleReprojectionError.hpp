#ifndef ASLAM_BACKEND_SIMPLE_REPROJECTION_ERROR_HPP
#define ASLAM_BACKEND_SIMPLE_REPROJECTION_ERROR_HPP

#include <aslam/backend/ErrorTerm.hpp>
#include <aslam/backend/HomogeneousExpression.hpp>
#include <boost/shared_ptr.hpp>

namespace aslam {
namespace backend {

template <typename FRAME_T>
class SimpleReprojectionError : public ErrorTermFs<FRAME_T::KeypointDimension> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef FRAME_T frame_t;
    typedef typename frame_t::keypoint_t keypoint_t;
    typedef typename frame_t::camera_geometry_t camera_geometry_t;
    enum {
        // The dimension of the keypoint associated with this geometry policy
        KeypointDimension = frame_t::KeypointDimension
    };

    typedef Eigen::Matrix<double, KeypointDimension, 1> measurement_t;
    typedef Eigen::Matrix<double, KeypointDimension, KeypointDimension> inverse_covariance_t;
    typedef ErrorTermFs<KeypointDimension> parent_t;

    SimpleReprojectionError();

    /**
     * @brief Construct with frame and target point, the measurement and covariance matrix will be extracted from frame
     * by index
     * @param frame         Frame, include the measurement
     * @param keypointIndex Measurement point index in frame
     * @param point         Target point in camera frame, expressed in homogeneous coordinates
     */
    SimpleReprojectionError(const frame_t* frame, int keypointIndex, const HomogeneousExpression& point);

    /**
     * @brief Construct with measurement, target point and covariance matrix
     * @param y         Measurement, the point(corner) extract from image
     * @param invR      Point covariance matrix
     * @param point     Target point in camera frame, expressed in homogeneous coordinates
     * @param geometry  Camera geometry
     */
    SimpleReprojectionError(const measurement_t& y, const inverse_covariance_t& invR,
                            const HomogeneousExpression& point, const camera_geometry_t& geometry);

    virtual ~SimpleReprojectionError() = default;

  protected:
    /**
     * @brief Evaluate the error term
     * @return Error
     */
    virtual double evaluateErrorImplementation();

    /**
     * @brief Evaluate the jacobian
     * @param _jacobians  Jacobians
     */
    virtual void evaluateJacobiansImplementation(aslam::backend::JacobianContainer& _jacobians) const;

  protected:
    measurement_t _y;                    // measurement, the point in the frame(image)
    const camera_geometry_t* _geometry;  // camera geometry
    HomogeneousExpression _point;        // the homogeneous target point expressed in the camera frame
};

}  // namespace backend
}  // namespace aslam

#include "implementation/SimpleReprojectionError.hpp"

#endif /* ASLAM_BACKEND_SIMPLE_REPROJECTION_ERROR_HPP */
