// This file automatically generated by create_export_module.py
#define NO_IMPORT_ARRAY

#include <NumpyEigenConverter.hpp>

#include <boost/cstdint.hpp>

void import_D_6_uchar() {
    NumpyEigenConverter<Eigen::Matrix<boost::uint8_t, Eigen::Dynamic, 6> >::register_converter();
}
