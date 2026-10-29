cmake_minimum_required(VERSION 3.24)
project(cfar_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)

# pocketfft is a single header with no CMake project of its own, so point
# SOURCE_SUBDIR at a directory that does not exist and wrap it ourselves.
FetchContent_Declare(pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG 4d865cba24359b05b37da4e75aee24498c3b78e7 # branch cpp
  SOURCE_SUBDIR no-cmake)
FetchContent_Declare(tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
  GIT_TAG v3.4.0)
FetchContent_MakeAvailable(pocketfft tomlplusplus)

add_library(pocketfft INTERFACE)
target_include_directories(pocketfft SYSTEM INTERFACE ${pocketfft_SOURCE_DIR})

find_package(Threads REQUIRED) # pocketfft uses std::thread for its plan cache

set(CFAR_WARNINGS -Wall -Wextra -Wpedantic)

# Detector library (deliverable D-09). Must not depend on cfar_sim.
add_library(cfar src/spectrum.cpp)
target_include_directories(cfar PUBLIC include)
target_link_libraries(cfar PRIVATE pocketfft Threads::Threads)
target_compile_options(cfar PRIVATE ${CFAR_WARNINGS})

# Simulation code: signal and noise models, TOML config.
add_library(cfar_sim src/sim/signal.cpp src/sim/noise.cpp src/sim/config.cpp)
target_link_libraries(cfar_sim PUBLIC cfar PRIVATE tomlplusplus::tomlplusplus)
target_compile_options(cfar_sim PRIVATE ${CFAR_WARNINGS})

add_executable(sim_main apps/sim_main.cpp)
target_link_libraries(sim_main PRIVATE cfar_sim)
target_compile_options(sim_main PRIVATE ${CFAR_WARNINGS})

enable_testing()
find_package(GTest REQUIRED)
add_executable(test_spectrum tests/test_spectrum.cpp)
target_link_libraries(test_spectrum PRIVATE cfar_sim GTest::gtest_main)
target_compile_options(test_spectrum PRIVATE ${CFAR_WARNINGS})
target_compile_definitions(test_spectrum PRIVATE CFAR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
include(GoogleTest)
gtest_discover_tests(test_spectrum)
