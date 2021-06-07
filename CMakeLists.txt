cmake_minimum_required(VERSION 3.17)
project(sensor_spectra)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/include/eigen3/)
include_directories(/usr/include/glog/)
include_directories(.)
include_directories(./libs)

link_directories(/usr/local/lib/)

find_package(Ceres)

add_executable(sensor_spectra sensor_spectra.cpp rgb2spec.cpp spectrum.cpp xy_conversion.cpp xy_conversion.h sampler.cpp sampler.h camera_xy_lookup.cpp camera_xy_lookup.h ppm_image.cpp ppm_image.h spectral_image.cpp spectral_image.h k_means.cpp k_means.h rbf.cpp rbf.h color_correction.h lookup_2d.cpp lookup_2d.h)
target_link_libraries(sensor_spectra PRIVATE glog ceres)

add_executable(sensor_images sensor_images.cpp
        k_means.h k_means.cpp
        xy_conversion.h xy_conversion.cpp
        rbf.h rbf.cpp
        spectrum.h spectrum.cpp
        spectral_image.h spectral_image.cpp
        ppm_image.h ppm_image.cpp
        camera_xy_lookup.h camera_xy_lookup.cpp
        color_correction.h lookup_2d.cpp lookup_2d.h)

add_executable(sensor_correction sensor_correction.cpp
        k_means.h k_means.cpp
        xy_conversion.h xy_conversion.cpp
        rbf.h rbf.cpp
        spectrum.h spectrum.cpp
        spectral_image.h spectral_image.cpp
        ppm_image.h ppm_image.cpp
        camera_xy_lookup.h camera_xy_lookup.cpp
        color_correction.h lookup_2d.cpp lookup_2d.h)

add_executable(image_chromatic image_chromatic.cpp
        ppm_image.h ppm_image.cpp)

add_executable(binary_spectra binary_spectra.cpp
        spectrum.h spectrum.cpp
        rgb2spec.h)

add_custom_command(TARGET sensor_spectra PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/tables/ $<TARGET_FILE_DIR:sensor_spectra>/tables/)
add_custom_command(TARGET sensor_spectra PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/resources/ $<TARGET_FILE_DIR:sensor_spectra>/resources/)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -march=native -O3")
