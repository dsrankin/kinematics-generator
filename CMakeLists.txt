cmake_minimum_required(VERSION 3.16)
project(kinegen VERSION 0.2.0 LANGUAGES CXX)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(KG_BUILD_TESTS "Build unit tests" ON)
option(KG_BUILD_PYTHON "Build pybind11 Python bindings" OFF)
option(KG_BUILD_EXAMPLES "Build C++ examples" ON)
set(KG_WITH_ROOT "AUTO" CACHE STRING "ROOT output support: AUTO, ON or OFF")

include(FetchContent)

# --- yaml-cpp (system package if available, otherwise fetched) -------------
find_package(yaml-cpp 0.7 QUIET)
if(NOT yaml-cpp_FOUND)
  message(STATUS "yaml-cpp not found: fetching 0.8.0")
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
  set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  # yaml-cpp 0.8.0 declares an old cmake_minimum_required that CMake >= 4 rejects.
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG 0.8.0
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(yaml-cpp)
endif()
if(TARGET yaml-cpp::yaml-cpp)
  set(KG_YAML_TARGET yaml-cpp::yaml-cpp)
else()
  set(KG_YAML_TARGET yaml-cpp)
endif()

# --- core library (no ROOT dependency) ---------------------------------------
add_library(kinegen
  src/Expression.cc
  src/Distribution.cc
  src/Param.cc
  src/Kinematics.cc
  src/Detector.cc
  src/Features.cc
  src/Prior.cc
  src/Task.cc
  src/Truth.cc
  src/Observation.cc
  src/Serialization.cc
  src/Split.cc
  src/Legacy.cc)
target_include_directories(kinegen PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(kinegen PUBLIC ${KG_YAML_TARGET})
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(kinegen PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

# --- optional ROOT output ------------------------------------------------------
set(KG_HAVE_ROOT OFF)
if(NOT KG_WITH_ROOT STREQUAL "OFF")
  if(KG_WITH_ROOT STREQUAL "ON")
    find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Physics MathCore)
  else()
    find_package(ROOT QUIET COMPONENTS Core RIO Tree Physics MathCore)
  endif()
  if(ROOT_FOUND)
    set(KG_HAVE_ROOT ON)
    add_library(kinegen_root src/RootWriter.cc)
    target_link_libraries(kinegen_root PUBLIC kinegen ROOT::Core ROOT::RIO ROOT::Tree)
    message(STATUS "ROOT ${ROOT_VERSION} found: building kinegen_root")
  endif()
endif()

# --- command-line tools --------------------------------------------------------
add_executable(kg-convert-card tools/kg_convert_card.cc)
target_link_libraries(kg-convert-card PRIVATE kinegen)
add_executable(kg-generate tools/kg_generate.cc)
target_link_libraries(kg-generate PRIVATE kinegen)
if(KG_HAVE_ROOT)
  target_link_libraries(kg-generate PRIVATE kinegen_root)
  target_compile_definitions(kg-generate PRIVATE KG_HAVE_ROOT=1)
endif()

# --- examples --------------------------------------------------------------------
if(KG_BUILD_EXAMPLES)
  add_executable(example_stream examples/stream_example.cc)
  target_link_libraries(example_stream PRIVATE kinegen)
  target_compile_definitions(example_stream PRIVATE KG_EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples")
endif()

# --- tests --------------------------------------------------------------------------
if(KG_BUILD_TESTS)
  enable_testing()
  add_executable(kg_tests
    tests/test_main.cc
    tests/test_kinematics.cc
    tests/test_params.cc
    tests/test_validation.cc
    tests/test_truth.cc
    tests/test_detector.cc
    tests/test_features.cc
    tests/test_reproducibility.cc
    tests/test_legacy.cc
    tests/test_isr_correlations.cc
    tests/test_composite.cc)
  target_link_libraries(kg_tests PRIVATE kinegen)
  target_compile_definitions(kg_tests PRIVATE KG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME kg_tests COMMAND kg_tests)
  if(KG_HAVE_ROOT)
    add_executable(kg_root_tests tests/test_main.cc tests/test_root.cc)
    target_link_libraries(kg_root_tests PRIVATE kinegen_root ROOT::Physics ROOT::MathCore)
    target_compile_definitions(kg_root_tests PRIVATE KG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME kg_root_tests COMMAND kg_root_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()

# --- Python bindings -----------------------------------------------------------------
if(KG_BUILD_PYTHON)
  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  set(PYBIND11_FINDPYTHON ON) # reuse the interpreter found above
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found: fetching v2.13.6")
    FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.13.6
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(pybind11)
  endif()
  pybind11_add_module(kinegen_py python/bindings.cc)
  set_target_properties(kinegen_py PROPERTIES OUTPUT_NAME kinegen)
  target_link_libraries(kinegen_py PRIVATE kinegen)
  if(KG_BUILD_TESTS)
    add_test(NAME python_tests
      COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(python_tests PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:kinegen_py>")
  endif()
endif()
