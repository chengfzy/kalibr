// This file automatically generated by create_export_module.py
#define NO_IMPORT_ARRAY

#include <NumpyEigenConverter.hpp>

#include <boost/cstdint.hpp>

void import_1_D_uchar() {
    NumpyEigenConverter<Eigen::Matrix<boost::uint8_t, 1, Eigen::Dynamic> >::register_converter();
}
