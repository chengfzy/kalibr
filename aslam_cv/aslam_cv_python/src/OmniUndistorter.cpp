#include <aslam/OmniUndistorter.hpp>
#include <aslam/python/exportOmniUndistorter.hpp>
#include <numpy_eigen/boost_python_headers.hpp>

using namespace boost::python;
using namespace aslam;

void exportOmniUndistorter() {
    aslam::cameras::exportOmniUndistorter<aslam::cameras::NoMask>("OmniUndistorterNoMask");
    aslam::cameras::exportOmniUndistorter<aslam::cameras::ImageMask>("OmniUndistorterMasked");
}
