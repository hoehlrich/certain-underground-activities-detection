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

# core: pure DSP and detection. Links nothing else in this repo; no file I/O.
add_library(core core/spectrum.cpp core/cfar.cpp core/config.cpp core/detector.cpp)
target_include_directories(core PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(core PRIVATE pocketfft Threads::Threads)
target_compile_options(core PRIVATE ${CFAR_WARNINGS})

# io: config loading; later sample sources and loggers. Depends on core only.
add_library(io io/config_loader.cpp)
target_link_libraries(io PUBLIC core PRIVATE tomlplusplus::tomlplusplus)
target_compile_options(io PRIVATE ${CFAR_WARNINGS})

# sim: signal and noise models, sim config, Monte Carlo runners, detection theory.
add_library(sim sim/signal.cpp sim/noise.cpp sim/config.cpp sim/pfa.cpp sim/theory.cpp
  sim/pd.cpp)
target_link_libraries(sim PUBLIC core io PRIVATE tomlplusplus::tomlplusplus)
target_compile_options(sim PRIVATE ${CFAR_WARNINGS})

foreach(app sim_main pfa_check pd_sweep roc)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} PRIVATE sim)
  target_compile_options(${app} PRIVATE ${CFAR_WARNINGS})
endforeach()

enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)
foreach(test test_spectrum test_cfar test_pd)
  add_executable(${test} tests/${test}.cpp)
  target_link_libraries(${test} PRIVATE sim GTest::gtest_main)
  target_compile_options(${test} PRIVATE ${CFAR_WARNINGS})
  target_compile_definitions(${test} PRIVATE CFAR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  gtest_discover_tests(${test} DISCOVERY_TIMEOUT 30)
endforeach()

# Dependency rules, checked on the sources: core includes nothing from io, sim or
# apps; io includes nothing from sim or apps.
add_test(NAME boundary.core COMMAND ${CMAKE_COMMAND} -DDIR=${PROJECT_SOURCE_DIR}/core
  "-DFORBIDDEN=sim;io;apps" -P ${PROJECT_SOURCE_DIR}/cmake/check_includes.cmake)
add_test(NAME boundary.io COMMAND ${CMAKE_COMMAND} -DDIR=${PROJECT_SOURCE_DIR}/io
  "-DFORBIDDEN=sim;apps" -P ${PROJECT_SOURCE_DIR}/cmake/check_includes.cmake)
