cmake_minimum_required (VERSION 3.5)

project(point_data_lib)

add_library(point_data INTERFACE)
target_compile_definitions(point_data
INTERFACE LIBRARY_HEADER_ONLY)

install(
    FILES
    scr/PointData.hpp
    scr/PointData_0D.hpp
    scr/PointData_1D.hpp
    scr/PointData_2D.hpp
    scr/PointData_3D.hpp
    DESTINATION include
)

install(
    TARGETS point_data
    EXPORT library_targets
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)

install(
    EXPORT library_targets
    FILE point_data-config.cmake
    NAMESPACE liton::
    DESTINATION lib/cmake/point_data
)
