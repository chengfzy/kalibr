#ifndef ASLAM_CAMERAS_PINHOLE_PROJECTION_HPP
#define ASLAM_CAMERAS_PINHOLE_PROJECTION_HPP

#include <aslam/cameras/GridCalibrationTargetObservation.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/eigen.hpp>
#include <sm/PropertyTree.hpp>
#include <sm/boost/serialization.hpp>
#include <sm/kinematics/Transformation.hpp>
#include <sm/logging.hpp>
#include "StaticAssert.hpp"

namespace aslam {
namespace cameras {

/**
 * @brief
 *
 * @note NOTE by CC
 *  The frame defined in `kalibr`
 *    (1) Keypoint frame: the image frame F_uv
 *    (2) Euclidean frame: the camera frame F_C, sometime it means the normalized camera frame
 */
template <typename DISTORTION_T>
class PinholeProjection {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    enum { KeypointDimension = 2 };

    enum { IntrinsicsDimension = 4 };
    enum { DesignVariableDimension = IntrinsicsDimension };

    typedef DISTORTION_T distortion_t;
    typedef Eigen::Matrix<double, KeypointDimension, 1> keypoint_t;
    typedef Eigen::Matrix<double, KeypointDimension, IntrinsicsDimension> jacobian_intrinsics_t;

    /// \brief Default constructor
    PinholeProjection();

    PinholeProjection(double focalLengthU, double focalLengthV, double imageCenterU, double imageCenterV,
                      int resolutionU, int resolutionV, distortion_t distortion);

    PinholeProjection(double focalLengthU, double focalLengthV, double imageCenterU, double imageCenterV,
                      int resolutionU, int resolutionV);

    PinholeProjection(const sm::PropertyTree& config);

    /// \brief destructor.
    virtual ~PinholeProjection();

    template <typename DERIVED_P, typename DERIVED_K>
    bool euclideanToKeypoint(const Eigen::MatrixBase<DERIVED_P>& p,
                             const Eigen::MatrixBase<DERIVED_K>& outKeypoint) const;

    /**
     * @brief This is this reprojection procedure, reproject point from world frame(p, 3x1 vector) to image
     * frame(outKeypoint, 2x1 vector), include distortion and projection. And then calculate the Jacobian of outKeypoint
     * w.r.t to p. The Jacobian calculation using the chain rule, first calculate the Jacobian(J1) of distorted point
     * w.r.t p, then calculate the Jacobians(J2) of outKeypoint w.r.t distorted point, so the final Jacobian J = J2 * J1
     */
    template <typename DERIVED_P, typename DERIVED_K, typename DERIVED_JP>
    bool euclideanToKeypoint(const Eigen::MatrixBase<DERIVED_P>& p, const Eigen::MatrixBase<DERIVED_K>& outKeypoint,
                             const Eigen::MatrixBase<DERIVED_JP>& outJp) const;

    template <typename DERIVED_P, typename DERIVED_K>
    bool homogeneousToKeypoint(const Eigen::MatrixBase<DERIVED_P>& p,
                               const Eigen::MatrixBase<DERIVED_K>& outKeypoint) const;

    template <typename DERIVED_P, typename DERIVED_K, typename DERIVED_JP>
    bool homogeneousToKeypoint(const Eigen::MatrixBase<DERIVED_P>& p, const Eigen::MatrixBase<DERIVED_K>& outKeypoint,
                               const Eigen::MatrixBase<DERIVED_JP>& outJp) const;

    // NOTE by CC: convert point in image frame p_uv to point in normalized camera frame. This procedure contains two
    // steps: (1) Back projection (2) Undistortion
    template <typename DERIVED_K, typename DERIVED_P>
    bool keypointToEuclidean(const Eigen::MatrixBase<DERIVED_K>& keypoint,
                             const Eigen::MatrixBase<DERIVED_P>& outPoint) const;

    template <typename DERIVED_K, typename DERIVED_P, typename DERIVED_JK>
    bool keypointToEuclidean(const Eigen::MatrixBase<DERIVED_K>& keypoint, const Eigen::MatrixBase<DERIVED_P>& outPoint,
                             const Eigen::MatrixBase<DERIVED_JK>& outJk) const;

    template <typename DERIVED_K, typename DERIVED_P>
    bool keypointToHomogeneous(const Eigen::MatrixBase<DERIVED_K>& keypoint,
                               const Eigen::MatrixBase<DERIVED_P>& outPoint) const;

    template <typename DERIVED_K, typename DERIVED_P, typename DERIVED_JK>
    bool keypointToHomogeneous(const Eigen::MatrixBase<DERIVED_K>& keypoint,
                               const Eigen::MatrixBase<DERIVED_P>& outPoint,
                               const Eigen::MatrixBase<DERIVED_JK>& outJk) const;

    template <typename DERIVED_P, typename DERIVED_JI>
    void euclideanToKeypointIntrinsicsJacobian(const Eigen::MatrixBase<DERIVED_P>& p,
                                               const Eigen::MatrixBase<DERIVED_JI>& outJi) const;

    template <typename DERIVED_P, typename DERIVED_JD>
    void euclideanToKeypointDistortionJacobian(const Eigen::MatrixBase<DERIVED_P>& p,
                                               const Eigen::MatrixBase<DERIVED_JD>& outJd) const;

    template <typename DERIVED_P, typename DERIVED_JI>
    void homogeneousToKeypointIntrinsicsJacobian(const Eigen::MatrixBase<DERIVED_P>& p,
                                                 const Eigen::MatrixBase<DERIVED_JI>& outJi) const;

    template <typename DERIVED_P, typename DERIVED_JD>
    void homogeneousToKeypointDistortionJacobian(const Eigen::MatrixBase<DERIVED_P>& p,
                                                 const Eigen::MatrixBase<DERIVED_JD>& outJd) const;

    /**
     * @brief Whether the point is in image
     *
     * @tparam DERIVED_K
     * @param keypoint
     * @return
     */
    template <typename DERIVED_K>
    bool isValid(const Eigen::MatrixBase<DERIVED_K>& keypoint) const;

    template <typename DERIVED_P>
    bool isEuclideanVisible(const Eigen::MatrixBase<DERIVED_P>& p) const;

    template <typename DERIVED_P>
    bool isHomogeneousVisible(const Eigen::MatrixBase<DERIVED_P>& ph) const;

    // aslam::backend compatibility
    void update(const double* v);
    int minimalDimensions() const;
    void getParameters(Eigen::MatrixXd& P) const;
    void setParameters(const Eigen::MatrixXd& P);
    Eigen::Vector2i parameterSize() const;

    enum { CLASS_SERIALIZATION_VERSION = 0 };
    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template <class Archive>
    void load(Archive& ar, const unsigned int version);
    template <class Archive>
    void save(Archive& ar, const unsigned int version) const;

    // \brief creates a random valid keypoint.
    virtual Eigen::VectorXd createRandomKeypoint() const;

    // \brief creates a random visible point. Negative depth means random between 0 and 100 meters.
    virtual Eigen::Vector3d createRandomVisiblePoint(double depth = -1.0) const;

    bool isProjectionInvertible() const { return false; }

    void setDistortion(const distortion_t& distortion) { _distortion = distortion; }
    distortion_t& distortion() { return _distortion; };
    const distortion_t& distortion() const { return _distortion; };

    Eigen::Matrix3d getCameraMatrix() {
        Eigen::Matrix3d K;
        K << _fu, 0.0, _cu, 0.0, _fv, _cv, 0.0, 0.0, 1.0;
        return K;
    };

    double focalLengthCol() const { return _fu; }
    double focalLengthRow() const { return _fv; }
    double opticalCenterCol() const { return _cu; }
    double opticalCenterRow() const { return _cv; }

    /// \brief The horizontal focal length in pixels.
    double fu() const { return _fu; }
    /// \brief The vertical focal length in pixels.
    double fv() const { return _fv; }
    /// \brief The horizontal image center in pixels.
    double cu() const { return _cu; }
    /// \brief The vertical image center in pixels.
    double cv() const { return _cv; }
    /// \brief The horizontal resolution in pixels.
    int ru() const { return _ru; }
    /// \brief The vertical resolution in pixels.
    int rv() const { return _rv; }
    /// \brief The horizontal resolution in pixels.
    int width() const { return _ru; }
    /// \brief The vertical resolution in pixels.
    int height() const { return _rv; }

    int keypointDimension() const { return KeypointDimension; }

    bool isBinaryEqual(const PinholeProjection<distortion_t>& rhs) const;

    static PinholeProjection<distortion_t> getTestProjection();

    /// \brief resize the intrinsics based on a scaling of the image.
    void resizeIntrinsics(double scale);

    /// \brief Get a set of border rays
    void getBorderRays(Eigen::MatrixXd& rays);

    /// \brief initialize the intrinsics based on a list of views of a gridded calibration target
    /// \return true on success
    bool initializeIntrinsics(const std::vector<GridCalibrationTargetObservation>& observations);

    /// \brief compute the reprojection error based on a checkerboard observation.
    /// \return the number of corners successfully observed and projected
    size_t computeReprojectionError(const GridCalibrationTargetObservation& obs,
                                    const sm::kinematics::Transformation& T_target_camera, double& outErr) const;

    /// \brief estimate the transformation of the camera with respect to the calibration target
    ///        On success out_T_t_c is filled in with the transformation that takes points from
    ///        the camera frame to the target frame
    /// \return true on success
    bool estimateTransformation(const GridCalibrationTargetObservation& obs,
                                sm::kinematics::Transformation& out_T_t_c) const;

  private:
    void updateTemporaries();

    /// \brief The horizontal focal length in pixels.
    double _fu;
    /// \brief The vertical focal length in pixels.
    double _fv;
    /// \brief The horizontal image center in pixels.
    double _cu;
    /// \brief The vertical image center in pixels.
    double _cv;
    /// \brief The horizontal resolution in pixels.
    int _ru;
    /// \brief The vertical resolution in pixels.
    int _rv;

    /// \brief A computed value for speeding up computation.
    double _recip_fu;
    double _recip_fv;
    double _fu_over_fv;

    distortion_t _distortion;
};

}  // namespace cameras
}  // namespace aslam

#include "implementation/PinholeProjection.hpp"

SM_BOOST_CLASS_VERSION_T1(aslam::cameras::PinholeProjection);

#endif /* ASLAM_CAMERAS_PINHOLE_PROJECTION_HPP */
