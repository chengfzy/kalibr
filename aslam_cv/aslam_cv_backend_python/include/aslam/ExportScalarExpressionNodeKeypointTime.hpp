// It is extremely important to use this header
// if you are using the numpy_eigen interface
#include <aslam/backend/ScalarExpressionNodeKeypointTime.hpp>
#include <numpy_eigen/boost_python_headers.hpp>
#include <sm/python/stl_converters.hpp>

namespace aslam {
namespace python {

template <typename CAMERA_T>
void exportScalarExpressionNodeKeypointTime(std::string name) {
    using namespace boost::python;
    using namespace aslam;

    class_<backend::ScalarExpressionNodeKeypointTime<CAMERA_T>,
           boost::shared_ptr<backend::ScalarExpressionNodeKeypointTime<CAMERA_T> > >(
        (name + "ScalarExpressionNodeKeypointTime").c_str(),
        init<const aslam::Time &, const Eigen::VectorXd &,
             boost::shared_ptr<aslam::backend::CameraDesignVariable<CAMERA_T> > >(
            "ScalarExpressionNodeKeypointTime(frame timestamp, keypoint, camera design variable)"))
        .def("toScalar", &backend::ScalarExpressionNodeKeypointTime<CAMERA_T>::toScalar);
}
}  // namespace python
}  // namespace aslam
