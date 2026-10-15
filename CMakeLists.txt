cmake_minimum_required(VERSION 3.20)
project(kbsrc VERSION 1.3.0 DESCRIPTION "Kaiser-Bessel sample-rate converter" LANGUAGES CXX)

if(MSVC)
  message(FATAL_ERROR "kbsrc needs GCC or Clang (uses __int128 and x87 long double)")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(KBSRC_ARCH "native" CACHE STRING
    "native: tune for the build CPU; portable: any CPU of this architecture (x86-64: AVX2 engine picked at run time)")
set_property(CACHE KBSRC_ARCH PROPERTY STRINGS native portable)
option(KBSRC_BUILD_TESTS "Register the Python test suite with CTest" ON)

find_package(Threads REQUIRED)

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" KBSRC_CPU)
set(KBSRC_X86 OFF)
if(KBSRC_CPU MATCHES "^(x86_64|amd64|x64)$")
  set(KBSRC_X86 ON)
endif()

# Flags shared by every target. Required: the double-double arithmetic relies on exact IEEE rounding of every
# operation, so no contraction into FMA and no fast-math.
add_library(kbsrc_options INTERFACE)
target_compile_options(kbsrc_options INTERFACE -ffp-contract=off -Wall -Wextra -Wno-ignored-attributes)
if(MINGW)
  # GCC on Windows cannot guarantee 32-byte stack alignment for AVX spills.
  target_compile_options(kbsrc_options INTERFACE -Wa,-muse-unaligned-vector-move)
endif()
if(KBSRC_ARCH STREQUAL "native")
  if(KBSRC_CPU MATCHES "^(arm|aarch64|arm64)")
    target_compile_options(kbsrc_options INTERFACE -mcpu=native)
  else()
    target_compile_options(kbsrc_options INTERFACE -march=native)
  endif()
elseif(NOT KBSRC_ARCH STREQUAL "portable")
  message(FATAL_ERROR "KBSRC_ARCH must be native or portable")
endif()

# The engine is built once per instruction set, each copy in its own namespace (src/isa.hpp).
set(KBSRC_ENGINE_SOURCES src/resampler.cpp src/design.cpp src/lpc.cpp)
function(kbsrc_engine name)
  add_library(kbsrc_engine_${name} OBJECT ${KBSRC_ENGINE_SOURCES})
  target_include_directories(kbsrc_engine_${name} PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)
  target_compile_definitions(kbsrc_engine_${name} PRIVATE KBSRC_ISA=${name})
  target_compile_options(kbsrc_engine_${name} PRIVATE ${ARGN})
  target_link_libraries(kbsrc_engine_${name} PRIVATE kbsrc_options)
endfunction()
kbsrc_engine(base)
set(KBSRC_ENGINE_OBJECTS $<TARGET_OBJECTS:kbsrc_engine_base>)
if(KBSRC_ARCH STREQUAL "portable" AND KBSRC_X86)
  kbsrc_engine(avx2 -march=x86-64-v3 -mtune=generic)
  list(APPEND KBSRC_ENGINE_OBJECTS $<TARGET_OBJECTS:kbsrc_engine_avx2>)
endif()

add_library(kbsrc_lib STATIC src/dispatch.cpp src/parallel.cpp ${KBSRC_ENGINE_OBJECTS})
add_library(kbsrc::kbsrc ALIAS kbsrc_lib)
set_target_properties(kbsrc_lib PROPERTIES OUTPUT_NAME kbsrc)
target_include_directories(kbsrc_lib
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
  PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_compile_definitions(kbsrc_lib PRIVATE KBSRC_VERSION="${PROJECT_VERSION}"
  $<$<TARGET_EXISTS:kbsrc_engine_avx2>:KBSRC_ENGINE_AVX2>)
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
    add_test(NAME rates COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_rates.py
             --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME formats COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_formats.py
             --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME robust COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_robust.py
             --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME alias COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_alias.py
             --bin $<TARGET_FILE:kbsrc>)
    add_test(NAME reference COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_reference.py
             --bin $<TARGET_FILE:kbsrc> --data ${KBSRC_TEST_DATA})
    set_tests_properties(reference PROPERTIES SKIP_RETURN_CODE 77)
    # Engine identity and rounding-mode independence of Exact; the driver calls the library with the caller in a
    # chosen rounding mode, which the CLI cannot do.
    add_executable(kbsrc_fenv_driver tests/fenv_driver.cpp)
    target_link_libraries(kbsrc_fenv_driver PRIVATE kbsrc_lib kbsrc_options)
    if(MINGW)
      target_link_options(kbsrc_fenv_driver PRIVATE -static)
    endif()
    add_test(NAME engines COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_engines.py
             --bin $<TARGET_FILE:kbsrc> --driver $<TARGET_FILE:kbsrc_fenv_driver>)
    # A portable x86-64 build must run on CPUs without AVX: no AVX instruction outside the AVX2 engine.
    if(TARGET kbsrc_engine_avx2)
      add_test(NAME isa COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/test_isa.py
               --bin $<TARGET_FILE:kbsrc> --objdump "${CMAKE_OBJDUMP}")
      set_tests_properties(isa PROPERTIES SKIP_RETURN_CODE 77)
    endif()
  endif()
endif()
