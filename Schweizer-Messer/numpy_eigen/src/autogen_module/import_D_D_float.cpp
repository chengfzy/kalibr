// This file automatically generated by create_export_module.py
#define NO_IMPORT_ARRAY

#include <NumpyEigenConverter.hpp>

#include <boost/cstdint.hpp>

void import_D_D_float() {
    NumpyEigenConverter<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> >::register_converter();
}
