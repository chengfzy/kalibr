#include "cc/KnotSequenceUpdateStrategy.h"
#include <fmt/format.h>
#include <glog/logging.h>
#include "cc/Heading.hpp"

using namespace std;
using namespace cc;
using namespace fmt;
using namespace aslam;
using namespace Eigen;

// const double kEpsilon = 1e-8;
const double kEpsilon = 0;

KnotSequenceUpdateStrategy::KnotSequenceUpdateStrategy(const double& frameRate)
    : frameRate_(frameRate), maxKnotsPerSecond_(1. / 2. / frameRate) {}

// Generate a new knots list given the current spline and reprojection errors
bool KnotSequenceUpdateStrategy::generateKnots(
    const vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>>& reprojectionErrors,
    const bsplines::BSplinePose& poseSpline, vector<double>& updatedKnots) {
    // evaluate the timestamp and reprojection error values from the input reprojection error terms
    vector<double> times, errors;
    getErrorAndTimestamp(reprojectionErrors, times, errors);

    // take a copy of the old knots
    vector<double> knots = poseSpline.knots();

    // get the errors and errors term numbers in each pose spline segment
    vector<double> errorPerSegment;
    vector<double> errorTermsPerSegment;
    getErrorPerSegment(times, errors, knots, errorPerSegment, errorTermsPerSegment);

    // remove segments
    removeSegmentWithoutImprovement(times, errors);
    vector<double> filteredKnots = removeSegmentWithoutObservations(knots, errorPerSegment);
    vector<double> errorPerSegmentFiltered;
    vector<double> errorTermsPerSegmentFiltered;
    getErrorPerSegment(times, errors, filteredKnots, errorPerSegmentFiltered, errorTermsPerSegmentFiltered);

    // generate new knots
    updatedKnots.clear();
    for (size_t i = 0; i < errorPerSegmentFiltered.size(); ++i) {
        // FIXME by CC: expected error should be 2N because the dimension of reprojection error is 2. see paper
        double expectedNormalError = errorTermsPerSegmentFiltered[i];
        if (errorPerSegmentFiltered[i] > expectedNormalError && i < knots.size() - 1) {
            double newKnot = (knots[i] + knots[i + 1]) / 2.0;
            double deltaT = newKnot - knots[i];
            if (deltaT + kEpsilon <= maxKnotsPerSecond_) {
                // reach the max number of knots per second: do not split
                updatedKnots.emplace_back(knots[i]);
            } else if (!disabledTimeSegments_.empty() && isSegmentDisabled(newKnot)) {
                // segment is disabled: don't split
                updatedKnots.emplace_back(knots[i]);
            } else {
                // split
                updatedKnots.emplace_back(knots[i]);
                updatedKnots.emplace_back(newKnot);
            }
        } else {
            updatedKnots.emplace_back(knots[i]);
        }
    }

    bool requiresUpdate{false};
    if (previousKnots_.empty()) {
        requiresUpdate = true;
    } else {
        // require at least a 1% increase in knots for a next iteration being worth the effort
        requiresUpdate = updatedKnots.size() > 1.01 * previousKnots_.size();
    }

    //  keep a copy of the knot sequence
    previousKnots_ = updatedKnots;
    previousErrors_ = errorPerSegmentFiltered;

    return requiresUpdate;
}

// Get a spline with the ne knots build upon the poses of the old spline
bsplines::BSplinePose KnotSequenceUpdateStrategy::getUpdatedSpline(const bsplines::BSplinePose& poseSpline,
                                                                   vector<double>& knots, int splineOrder) {
    // linear sample the old spline
    size_t kN = knots.size();
    vector<double> times(kN);
    double dt = (poseSpline.tMax() - poseSpline.tMin()) / (static_cast<int>(kN) - 1);
    // generate(times.begin(), times.end(), [&, t = poseSpline.tMin()]() mutable { return t + dt; });
    for (size_t i = 0; i < kN; ++i) {
        times[i] = poseSpline.tMin() + i * dt;
    }

    // spline poses
    MatrixXd splinePoses(6, kN);
    for (size_t i = 0; i < kN; ++i) {
        splinePoses.col(i) = poseSpline.eval(times[i]);
    }

    // guarantee that begin and end times of the spline remain unchanged
    vector<double> oldKnots = poseSpline.knots();
    int i{0};
    while (oldKnots[i] + kEpsilon < knots[0]) {
        ++i;
    }
    knots.insert(knots.begin(), oldKnots.begin(), oldKnots.begin() + i);
    i = static_cast<int>(oldKnots.size()) - 1;
    while (oldKnots[i] - kEpsilon > *knots.rbegin()) {
        --i;
    }
    knots.insert(knots.end(), oldKnots.begin() + i, oldKnots.end());

    // generate new pose spline
    bsplines::BSplinePose newPoseSpline(splineOrder, poseSpline.rotation());
    newPoseSpline.initPoseSplineSparseKnots(Map<VectorXd>(times.data(), times.size()), splinePoses,
                                            Map<VectorXd>(knots.data(), knots.size()), 1e-6);

    return newPoseSpline;
}

// evaluate the timestamp and reprojection error values from the input reprojection error terms, and sort it by
// timestamp
void KnotSequenceUpdateStrategy::getErrorAndTimestamp(
    const vector<boost::shared_ptr<AdaptiveCovarianceReprojectionError>>& reprojectionErrors, vector<double>& times,
    vector<double>& errors) {
    vector<double> rawTimes;
    vector<double> rawErrors;
    for (auto& v : reprojectionErrors) {
        rawTimes.emplace_back(v->observationTime());
        rawErrors.emplace_back(v->evaluateError());
    }

    // sort by time
    times.resize(rawTimes.size());
    errors.resize(rawTimes.size());
    vector<size_t> index(rawTimes.size());
    iota(index.begin(), index.end(), 0);
    sort(index.begin(), index.end(), [&](const size_t& i0, const size_t& i1) { return rawTimes[i0] < rawTimes[i1]; });
    transform(index.begin(), index.end(), times.begin(), [&](const size_t& i) { return rawTimes[i]; });
    transform(index.begin(), index.end(), errors.begin(), [&](const size_t& i) { return rawErrors[i]; });
}

// get the errors and errors term numbers in each pose spline segment
void KnotSequenceUpdateStrategy::getErrorPerSegment(const std::vector<double>& times, const std::vector<double>& errors,
                                                    const std::vector<double>& knots, vector<double>& errorPerSegment,
                                                    vector<double>& errorTermsPerSegment) {
    errorPerSegment.resize(knots.size(), 0);       // error in each spline segment
    errorTermsPerSegment.resize(knots.size(), 0);  // error term numbers in each spline segment
    pair<int, int> segment(-1, -1);
    for (size_t n = 0; n < times.size(); ++n) {
        segment = timeToKnotSection(times[n], knots, segment);
        CHECK(segment.first != -1) << format("the found segment = ({}, {}) is error", segment.first, segment.second);
        errorPerSegment[segment.first] += errors[n];
        errorTermsPerSegment[segment.first] += 1;
    }
}

// get the knot section [n0, n1) the time located in
std::pair<int, int> KnotSequenceUpdateStrategy::timeToKnotSection(const double& t, const std::vector<double>& knots,
                                                                  const std::pair<int, int>& segment) {
    int n = segment.first;
    if (-1 == n) {
        n = 0;
    }

    while (n < static_cast<int>(knots.size() - 1)) {
        if (knots[n] + kEpsilon <= t && t + kEpsilon < knots[n + 1]) {
            return {n, n + 1};
        }
        ++n;
    }

    return {-1, -1};
}

// remove segment without any improvement
void KnotSequenceUpdateStrategy::removeSegmentWithoutImprovement(const std::vector<double>& times,
                                                                 const std::vector<double>& errors) {
    // first compare the reprojection error in the previous segments, if we don't observe a significant drop, stop
    // adding errors
    if (!previousKnots_.empty() && !previousErrors_.empty()) {
        vector<pair<double, double>> timeSegments;
        vector<double> errorPerSegment;  // error in each segment in current solution
        errorPerSegment.resize(previousKnots_.size(), 0);
        pair<int, int> segment(-1, -1);
        // analyze each section of the knot sequence
        for (size_t i = 0; i < times.size(); ++i) {
            segment = timeToKnotSection(times[i], previousKnots_, segment);
            CHECK(segment.first != -1) << format("the found segment = ({}, {}) is error", segment.first,
                                                 segment.second);
            errorPerSegment[segment.first] += errors[i];
            timeSegments.emplace_back(make_pair(previousKnots_[segment.first], previousKnots_[segment.second]));
        }

        // compare to previous errors, if we don't observe a significant drop, stop adding errors
        for (size_t i = 0; i < previousErrors_.size(); ++i) {
            if (previousErrors_[i] * 0.8 < errorPerSegment[i]) {
                disabledTimeSegments_.emplace_back(timeSegments[i]);
            }
        }
    }
}

// remove segment without observations, return filtered knots
vector<double> KnotSequenceUpdateStrategy::removeSegmentWithoutObservations(
    const std::vector<double>& knots, const std::vector<double>& errorPerSegment) {
    vector<double> filteredKnots;  // after removed errors

    // remove segment with consecutive 0-valued errors
    double pError{0};  // previous error
    for (size_t i = 0; i < errorPerSegment.size(); ++i) {
        // this should depend on the spline order, FIXME by CC: seems the value should be 4 instead of 6?
        if (pError == 0 && errorPerSegment[i] == 0 && 6 < i && i < errorPerSegment.size() - 6) {
            // add the segment between the previous and the next knot to the blacklist
            disabledTimeSegments_.emplace_back(knots[i - 1], knots[i + 1]);
        } else {
            filteredKnots.emplace_back(knots[i]);
        }
        pError = errorPerSegment[i];
    }
    return filteredKnots;
}

// check whether the segment is disabled
bool KnotSequenceUpdateStrategy::isSegmentDisabled(const double& t) {
    for (auto& s : disabledTimeSegments_) {
        if (s.first <= t && t < s.second) {
            return true;
        }
    }
    return false;
}