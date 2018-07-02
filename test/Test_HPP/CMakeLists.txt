project (Test_HPP)

add_library(sub STATIC "Test_HPP_sub.cpp")
target_link_libraries(sub point_data)

add_executable("_${PROJECT_NAME}" "${PROJECT_NAME}.cpp")
target_link_libraries("_${PROJECT_NAME}" point_data sub)

add_test(NAME "${PROJECT_NAME}_" COMMAND "_${PROJECT_NAME}")
