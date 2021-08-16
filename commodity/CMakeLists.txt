cmake_minimum_required(VERSION 2.8)
project(commodity)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../modules/)
find_package(GMP REQUIRED)

set(CMAKE_C_COMPILER c++)

find_library(PRICOMP NAMES pricomp HINTS ${CMAKE_SOURCE_DIR}/../lib/)

find_library(SMC_TWO NAMES smc_two HINTS ${CMAKE_SOURCE_DIR}/../lib/)

include_directories(../src/ ./)

add_executable(commodity commodity.c commodity.h)
add_executable(genOfflineRandomBits genOfflineRandomBits.c commodity.h)
add_executable(genCommodityRandomBits genCommodityRandomBits.c)

target_link_libraries(commodity ${PRICOMP} ${SMC_TWO} ${GMP_LIBRARIES})
