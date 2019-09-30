#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <numpy_eigen/boost_python_headers.hpp>

void import_bspline_python();
void import_bspline_pose_python();

BOOST_PYTHON_MODULE(libbsplines_python) {
    import_bspline_python();
    import_bspline_pose_python();
}
