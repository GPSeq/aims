cmake_minimum_required(VERSION 3.24)

project(AIMS
  VERSION 0.1.0
  DESCRIPTION "Adaptive Information-per-Memory-access Seed Index research system"
  LANGUAGES CXX)

option(AIMS_ENABLE_ASAN "Enable AddressSanitizer and UndefinedBehaviorSanitizer in Debug builds" ON)
option(AIMS_ENABLE_LTO "Enable interprocedural optimization for Release builds when supported" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIPOSupported)
include(CTest)

find_package(Doxygen QUIET)
find_program(LATEXMK_EXECUTABLE latexmk)

add_library(aims
  src/build/index_builder.cpp
  src/codecs/delta_varint.cpp
  src/codecs/posting_block_codec.cpp
  src/index/exact_dictionary.cpp
  src/index/postings.cpp
  src/instrumentation/metrics.cpp
  src/io/fasta.cpp
  src/io/mapped_file.cpp
  src/query/candidate_accumulator.cpp
  src/query/kmer_search.cpp
  src/query/query_planner.cpp
  src/seeds/kmer_generator.cpp
  src/serialization/index_format.cpp)

target_include_directories(aims PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_compile_features(aims PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(aims PRIVATE /W4 /permissive-)
else()
  target_compile_options(aims PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

if(AIMS_ENABLE_ASAN AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT MSVC)
  target_compile_options(aims PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(aims PUBLIC -fsanitize=address,undefined)
endif()

if(AIMS_ENABLE_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
  check_ipo_supported(RESULT aims_ipo_supported OUTPUT aims_ipo_error)
  if(aims_ipo_supported)
    set_property(TARGET aims PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO/LTO not supported: ${aims_ipo_error}")
  endif()
endif()

if(DOXYGEN_FOUND)
  add_custom_target(aims_doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating AIMS Doxygen HTML documentation"
    VERBATIM)
endif()

if(LATEXMK_EXECUTABLE)
  add_custom_target(aims_algorithm_pdf
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/texmf-var
    COMMAND ${CMAKE_COMMAND} -E env
      TEXMFVAR=${CMAKE_CURRENT_BINARY_DIR}/texmf-var
      ${LATEXMK_EXECUTABLE}
      -lualatex
      -interaction=nonstopmode
      -halt-on-error
      -outdir=${CMAKE_CURRENT_SOURCE_DIR}/docs
      ${CMAKE_CURRENT_SOURCE_DIR}/docs/algorithm.tex
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating AIMS algorithm PDF"
    VERBATIM)
endif()

foreach(app aims_build aims_query aims_dump aims_bench aims_validate)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} PRIVATE aims)
  target_compile_features(${app} PRIVATE cxx_std_20)
endforeach()

if(BUILD_TESTING)
  add_executable(aims_core_tests
    tests/test_core.cpp)
  target_link_libraries(aims_core_tests PRIVATE aims)
  add_test(NAME aims_core_tests COMMAND aims_core_tests)

  add_test(NAME aims_cli_smoke
    COMMAND ${CMAKE_COMMAND}
      -DAIMS_BUILD_EXE=$<TARGET_FILE:aims_build>
      -DAIMS_QUERY_EXE=$<TARGET_FILE:aims_query>
      -DAIMS_DUMP_EXE=$<TARGET_FILE:aims_dump>
      -DAIMS_VALIDATE_EXE=$<TARGET_FILE:aims_validate>
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_smoke.cmake)

  add_test(NAME aims_bench_smoke
    COMMAND $<TARGET_FILE:aims_bench>
      --ref ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/kmer_refs.fa
      --query ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/kmer_queries.fa
      --truth ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/kmer_truth.tsv
      --k 5,7
      --topk 3
      --dataset synthetic_kmer_fixture)
  set_tests_properties(aims_bench_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "\"top1_recall\"")

  add_test(NAME aims_generated_fixture_smoke
    COMMAND ${CMAKE_COMMAND}
      -DPYTHON_EXECUTABLE=python3
      -DAIMS_BENCH_EXE=$<TARGET_FILE:aims_bench>
      -DAIMS_GENERATOR=${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_synthetic_kmer_dataset.py
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/generated_fixture_smoke.cmake)
endif()
