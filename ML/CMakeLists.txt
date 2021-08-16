cmake_minimum_required(VERSION 3.0)
project(ml_demo)

set(CMAKE_MODULE_PATH ../modules/)
find_package(GMP REQUIRED)

find_library(PRICOMP NAMES pricomp HINTS ${CMAKE_SOURCE_DIR}/../lib/)
find_library(SMC_TWO NAMES smc_two HINTS ${CMAKE_SOURCE_DIR}/../lib/)

add_subdirectory(ml_protocol)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/ ./ ${GMP_INCLUDE_DIR})

add_executable(train_secure_nb train_secure_nb.cpp)
add_executable(test_secure_nb test_secure_nb.cpp)
add_executable(train_secure_kmeans train_secure_kmeans.cpp)
add_executable(test_secure_kmeans test_secure_kmeans.cpp)
add_executable(demo_knn demo_knn.cpp)

target_link_libraries(test_secure_nb ${SMC_TWO} ${PRICOMP} ml_protocol ${GMP_LIBRARIES})
target_link_libraries(train_secure_nb ${SMC_TWO} ${PRICOMP} ml_protocol ${GMP_LIBRARIES})
target_link_libraries(demo_knn ${SMC_TWO} ${PRICOMP} ml_protocol ${GMP_LIBRARIES})
target_link_libraries(train_secure_kmeans ${SMC_TWO} ${PRICOMP} ml_protocol ${GMP_LIBRARIES})
target_link_libraries(test_secure_kmeans ${SMC_TWO} ${PRICOMP} ml_protocol ${GMP_LIBRARIES})
