cmake_minimum_required (VERSION 3.5)

project(liton_point_data)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(tinyxml2 REQUIRED)
find_package(snippets REQUIRED)
find_package(ordered_tec REQUIRED)

add_library(point_data INTERFACE)
target_include_directories(point_data INTERFACE "../scr")
target_compile_definitions(point_data INTERFACE 
"$<$<CONFIG:DEBUG>:_CHECK_POINTDATA_RANGE>"
"$<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>"
OT_TINYXML)
target_link_libraries(point_data INTERFACE tinyxml2::tinyxml2 liton::snippets liton::ordered_tec)


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
