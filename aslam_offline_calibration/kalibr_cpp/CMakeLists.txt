cmake_minimum_required(VERSION 3.2)
project(kalibr_cpp)

set(CMAKE_CXX_STANDARD 11)

find_package(catkin REQUIRED COMPONENTS roscpp rosbag cv_bridge cmake_modules)
catkin_simple()

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})

add_definitions(-DDebugTest)

##################################
# error terms (+python export)
##################################
include_directories(./include)

# IMU-Camera Calibration
file(GLOB_RECURSE FILE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.*)
file(GLOB_RECURSE FILE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.*)
list(REMOVE_ITEM FILE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/RollingShutterCameraCalibration.cpp)
add_executable(ImuCameraCalibration ${FILE_INCLUDES} ${FILE_SRCS})
target_link_libraries(ImuCameraCalibration ${catkin_LIBRARIES} fmt::fmt-header-only glog gflags)

# Rolling Camera Calibration
file(GLOB_RECURSE FILE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.*)
file(GLOB_RECURSE FILE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.*)
list(REMOVE_ITEM FILE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/ImuCameraCalibration.cpp)
add_executable(RollingShutterCameraCalibration ${FILE_INCLUDES} ${FILE_SRCS})
target_link_libraries(RollingShutterCameraCalibration ${catkin_LIBRARIES} fmt::fmt-header-only glog gflags)

##################################
# tests
##################################


##################################
# EXPORT
##################################
#catkin_python_setup()
cs_install()
cs_export()
