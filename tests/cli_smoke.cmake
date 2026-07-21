set(workdir "${CMAKE_CURRENT_BINARY_DIR}/aims_cli_smoke")
file(MAKE_DIRECTORY "${workdir}")

file(WRITE "${workdir}/refs.fa" ">ref0\nAACCGGTTAACC\n>ref1\nTTTTAACCGGTT\n")
file(WRITE "${workdir}/queries.fq" "@q0\nAACCGGTT\n+\nFFFFFFFF\n")

execute_process(
  COMMAND "${AIMS_BUILD_EXE}" --ref "${workdir}/refs.fa" --out "${workdir}/idx.aims" --k 4,6
  RESULT_VARIABLE build_result
  OUTPUT_VARIABLE build_stdout
  ERROR_VARIABLE build_stderr)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR "aims_build failed: ${build_stderr}")
endif()

execute_process(
  COMMAND "${AIMS_QUERY_EXE}" --index "${workdir}/idx.aims" --query "${workdir}/queries.fq" --topk 2 --emit jsonl --out "${workdir}/query_results.jsonl" --mmap --posting-cache-blocks 2 --threads 2
  RESULT_VARIABLE query_result
  OUTPUT_VARIABLE query_stdout
  ERROR_VARIABLE query_stderr)
if(NOT query_result EQUAL 0)
  message(FATAL_ERROR "aims_query failed: ${query_stderr}")
endif()

file(READ "${workdir}/query_results.jsonl" query_output_file)

string(FIND "${query_output_file}" "\"comparison_stage\":\"exact_retrieval\"" stage_pos)
if(stage_pos EQUAL -1)
  message(FATAL_ERROR "query output missing comparison_stage: ${query_output_file}")
endif()

string(FIND "${query_output_file}" "\"postings_decoded\":" postings_pos)
if(postings_pos EQUAL -1)
  message(FATAL_ERROR "query output missing postings_decoded: ${query_output_file}")
endif()

string(FIND "${query_output_file}" "\"topk_results\":" topk_pos)
if(topk_pos EQUAL -1)
  message(FATAL_ERROR "query output missing topk_results: ${query_output_file}")
endif()

execute_process(
  COMMAND "${AIMS_VALIDATE_EXE}" --index "${workdir}/idx.aims"
  RESULT_VARIABLE validate_result
  OUTPUT_VARIABLE validate_stdout
  ERROR_VARIABLE validate_stderr)
if(NOT validate_result EQUAL 0)
  message(FATAL_ERROR "aims_validate failed: ${validate_stderr}")
endif()

string(FIND "${validate_stdout}" "format=KmerExactIndex" validate_pos)
if(validate_pos EQUAL -1)
  message(FATAL_ERROR "validate output missing format: ${validate_stdout}")
endif()

execute_process(
  COMMAND "${AIMS_DUMP_EXE}" --index "${workdir}/idx.aims"
  RESULT_VARIABLE dump_result
  OUTPUT_VARIABLE dump_stdout
  ERROR_VARIABLE dump_stderr)
if(NOT dump_result EQUAL 0)
  message(FATAL_ERROR "aims_dump failed: ${dump_stderr}")
endif()

string(FIND "${dump_stdout}" "KmerExactIndex" dump_pos)
if(dump_pos EQUAL -1)
  message(FATAL_ERROR "dump output missing KmerExactIndex: ${dump_stdout}")
endif()
