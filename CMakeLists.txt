cmake_minimum_required(VERSION 3.16)
project(bspline_interportation)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(/usr/local/include/eigen3)

find_package(Eigen3)
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)

aux_source_directory(. DIRS_SRCS)

add_executable(${PROJECT_NAME} ${DIRS_SRCS} BSpline.cpp BSpline.h)

target_link_libraries(${PROJECT_NAME} Qt5::Widgets)
target_link_libraries(${PROJECT_NAME} Qt5::Core)