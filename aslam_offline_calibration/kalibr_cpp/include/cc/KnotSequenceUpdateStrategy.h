#pragma once
#include <aslam/backend/CovarianceReprojectionError.hpp>
#include "aslam/cameras.hpp"

namespace cc {

/**
 * @brief Knot Sequence Update Strategies need two method:
 * 1) generateKnotList: which generates a new knot list given the current spline and reprojection errors. Returns a
 * boolean flag if the updated knot sequence needs another step of optimization.
 * 2) getUpdatedSpline: given a knot list,
 * spline order and spline generates a new spline initialized with the values of the given spline and the given knot
 * sequence.
 */
class KnotSequenceUpdateStrategy {
  public:
    // type define
    using CameraGeometry = aslam::cameras::DistortedPinholeRsCameraGeometry;
    using Frame = aslam::Frame<CameraGeometry>;
    using AdaptiveCovarianceReprojectionError =
        aslam::backend::CovarianceReprojectionError<aslam::Frame<CameraGeometry>>;

  public:
    explicit KnotSequenceUpdateStrategy(const double& frameRate);

    /**
     * @brief Generate a new knots list given the current spline and reprojection errors
     * @param reprojectionErrors    Current reprojection errors
     * @param poseSpline            Pose spline
     * @return  True if the updated knot sequence need another step of optimization
     */
    bool generateKnots(const std::vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>>& reprojectionErrors,
                       const bsplines::BSplinePose& poseSpline, std::vector<double>& updatedKnots);

    // Get a spline with the ne knots build on the poses of the old spline
    bsplines::BSplinePose getUpdatedSpline(const bsplines::BSplinePose& poseSpline, std::vector<double>& knots,
                                           int splineOrder);

  private:
    // extract the timestamps and reprojection error values
    void getErrorAndTimestamp(
        const std::vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>>& reprojectionErrors,
        std::vector<double>& times, std::vector<double>& errors);

    // get the total error per segment and number of error terms per segment
    void getErrorPerSegment(const std::vector<double>& times, const std::vector<double>& errors,
                            const std::vector<double>& knots, std::vector<double>& errorPerSegment,
                            std::vector<double>& errorTermsPerSegment);

    // get the knot section [n0, n1) the time locate in
    std::pair<int, int> timeToKnotSection(const double& t, const std::vector<double>& knots,
                                          const std::pair<int, int>& segment);

    // remove segment without any improvement, return filtered knots
    void removeSegmentWithoutImprovement(const std::vector<double>& times, const std::vector<double>& errors);

    // remove segment without observations
    std::vector<double> removeSegmentWithoutObservations(const std::vector<double>& knots,
                                                         const std::vector<double>& errorPerSegment);

    // check whether the segment is disabled
    bool isSegmentDisabled(const double& t);

  private:
    double frameRate_;          // frame rate
    double maxKnotsPerSecond_;  // max knot per second

    std::vector<double> previousKnots_;                            // previous knots sequence
    std::vector<double> previousErrors_;                           // previous error terms
    std::vector<std::pair<double, double>> disabledTimeSegments_;  // disabled time segments
};

}  // namespace cc