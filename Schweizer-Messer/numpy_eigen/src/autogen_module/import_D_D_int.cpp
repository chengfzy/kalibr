// This file automatically generated by create_export_module.py
#define NO_IMPORT_ARRAY

#include <NumpyEigenConverter.hpp>

#include <boost/cstdint.hpp>

void import_D_D_int() {
    NumpyEigenConverter<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> >::register_converter();
}
