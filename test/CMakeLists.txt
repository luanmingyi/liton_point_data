cmake_minimum_required (VERSION 3.5)

project(liton_point_data)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(tinyxml STATIC "../scr/tinyxml2/tinyxml2.cpp")

add_library(ordered_tec STATIC "../scr/liton_ordered_tec/ordered_tec.cpp")
target_compile_definitions(ordered_tec
PUBLIC "$<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>"
PUBLIC OT_TINYXML)
target_link_libraries(ordered_tec tinyxml)

add_library(point_data INTERFACE)
target_include_directories(point_data INTERFACE "../scr/liton_point_data")

enable_testing()

add_subdirectory(Test_D0_PointData_C)

add_subdirectory(Test_D1_RangeT)
add_subdirectory(Test_D1_SizeT)
add_subdirectory(Test_D1_PointData_C)
add_subdirectory(Test_D1_PointData_H)
add_subdirectory(Test_D1_tecfile)

add_subdirectory(Test_D2_RangeT)
add_subdirectory(Test_D2_SizeT)
add_subdirectory(Test_D2_PointData_C)
add_subdirectory(Test_D2_PointData_H)
add_subdirectory(Test_D2_tecfile)

add_subdirectory(Test_D3_RangeT)
add_subdirectory(Test_D3_SizeT)
add_subdirectory(Test_D3_PointData_C)
add_subdirectory(Test_D3_PointData_H)
add_subdirectory(Test_D3_tecfile)

add_subdirectory(Example_D1_heat)
add_subdirectory(Example_D1_integration)
add_subdirectory(Example_D2_heat)
add_subdirectory(Test_HPP)
