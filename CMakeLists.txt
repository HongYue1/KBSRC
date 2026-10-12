cmake_minimum_required(VERSION 3.20)
project(kbsrc VERSION 1.2.0 DESCRIPTION "Kaiser-Bessel fold-back sample-rate converter" LANGUAGES CXX)

if(MSVC)
  message(FATAL_ERROR "kbsrc needs GCC or Clang (uses __int128 and x87 long double)")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(KBSRC_NATIVE "Optimize for the build machine (-march=native). OFF targets x86-64-v3 (AVX2 + FMA)." ON)
option(KBSRC_BUILD_TESTS "Register the Python test suite with CTest" ON)

find_package(Threads REQUIRED)

add_library(kbsrc_options INTERFACE)
if(KBSRC_NATIVE)
  target_compile_options(kbsrc_options INTERFACE -march=native)
else()
  target_compile_options(kbsrc_options INTERFACE -march=x86-64-v3 -mtune=generic)
endif()
# Required: the double-double arithmetic relies on exact IEEE rounding of every operation.
target_compile_options(kbsrc_options INTERFACE -ffp-contract=off -Wall -Wextra -Wno-ignored-attributes)
if(MINGW)
  # GCC on Windows cannot guarantee 32-byte stack alignment for AVX spills.
  target_compile_options(kbsrc_options INTERFACE -Wa,-muse-unaligned-vector-move)
endif()

add_library(kbsrc_lib STATIC
  src/resampler.cpp
  src/design.cpp
  src/lpc.cpp
)
add_library(kbsrc::kbsrc ALIAS kbsrc_lib)
set_target_properties(kbsrc_lib PROPERTIES OUTPUT_NAME kbsrc)
target_include_directories(kbsrc_lib
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
  PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_compile_definitions(kbsrc_lib PRIVATE KBSRC_VERSION="${PROJECT_VERSION}")
target_link_libraries(kbsrc_lib PRIVATE kbsrc_options PUBLIC Threads::Threads)

add_library(kbsrc_io STATIC src/io/wav.cpp)
target_include_directories(kbsrc_io PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(kbsrc_io PRIVATE kbsrc_options)

add_executable(kbsrc app/main.cpp)
target_link_libraries(kbsrc PRIVATE kbsrc_lib kbsrc_io kbsrc_options)
if(WIN32)
  target_link_libraries(kbsrc PRIVATE shell32)
endif()
if(MINGW)
  target_link_options(kbsrc PRIVATE -static)
endif()

include(GNUInstallDirs)
install(TARGETS kbsrc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(KBSRC_BUILD_TESTS)
  enable_testing()
  find_package(Python3 COMPONENTS Interpreter)
  set(KBSRC_TEST_DATA "${PROJECT_SOURCE_DIR}/tests/data" CACHE PATH "Directory with the 64-bit float test WAVs")
  if(Python3_FOUND)
    add_test(NAME rates COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_rates.py --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME formats COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_formats.py --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME reference COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_reference.py
             --bin $<TARGET_FILE:kbsrc> --data ${KBSRC_TEST_DATA})
  endif()
endif()
