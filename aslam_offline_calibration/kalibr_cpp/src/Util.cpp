#include "cc/Util.h"
#include <glog/logging.h>

using namespace std;

namespace cc {

// Cross correlation of two 1D sequences
vector<double> correlate(const vector<double>& a, const vector<double>& v) {
    CHECK(a.size() == v.size()) << "input sequence size should be the same";
    const size_t kN = a.size();
    vector<double> c;
    for (size_t n = 0; n < 2 * kN - 1; ++n) {
        // start and end index for a and v
        size_t i0{0}, i1{0}, j0{0};
        if (n < kN) {
            i0 = 0;
            i1 = n + 1;
            j0 = kN - n - 1;
            // j1 = kN;
        } else {
            i0 = n - kN + 1;
            i1 = kN;
            j0 = 0;
            // j1 = 2 * kN - n - 1;
        }

        // sum
        double sum{0};
        for (size_t i = i0, j = j0; i < i1; ++i, ++j) {
            sum += a[i] * v[j];
        }
        c.emplace_back(sum);
    }

    return c;
}

}  // namespace cc