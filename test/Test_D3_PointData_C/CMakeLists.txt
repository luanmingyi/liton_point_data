project (Test_D3_PointData_C)

add_executable("_${PROJECT_NAME}" "${PROJECT_NAME}.cpp")
target_link_libraries("_${PROJECT_NAME}" point_data)

add_test(NAME "${PROJECT_NAME}_" COMMAND "_${PROJECT_NAME}")
