#pragma once
#include <vector>

namespace cc {

const std::size_t kCalibrationGroupId = 0;
const std::size_t kHelperGroupId = 1;

/**
 * @brief Cross correlation of two 1D sequences.
 *      c[k] = sum(a[n+k] * v[n])
 *  Ref: https://docs.scipy.org/doc/numpy/reference/generated/numpy.correlate.html
 * @param a First input sequence
 * @param v Second input sequence
 * @return Cross correlation of a and v
 */
std::vector<double> correlate(const std::vector<double>& a, const std::vector<double>& v);

}  // namespace cc
