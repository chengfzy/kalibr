#include <Eigen/Core>

#include <numpy_eigen/boost_python_headers.hpp>
Eigen::Matrix<boost::int64_t, 1, 3> test_long_1_3(const Eigen::Matrix<boost::int64_t, 1, 3>& M) { return M; }
void export_long_1_3() { boost::python::def("test_long_1_3", test_long_1_3); }
