set(workdir "${CMAKE_CURRENT_BINARY_DIR}/aims_generated_fixture")
file(REMOVE_RECURSE "${workdir}")
file(MAKE_DIRECTORY "${workdir}")

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "${AIMS_GENERATOR}" --out-dir "${workdir}" --documents 8 --length 240 --queries 12 --query-length 60
  RESULT_VARIABLE generator_result
  OUTPUT_VARIABLE generator_stdout
  ERROR_VARIABLE generator_stderr)
if(NOT generator_result EQUAL 0)
  message(FATAL_ERROR "synthetic generator failed: ${generator_stderr}")
endif()

execute_process(
  COMMAND "${AIMS_BENCH_EXE}" --ref "${workdir}/refs.fa" --query "${workdir}/queries.fa" --truth "${workdir}/truth.tsv" --k 9,13 --topk 5 --mmap --posting-cache-blocks 4 --repeats 2 --query-metrics-out "${workdir}/query_metrics.jsonl" --dataset generated_fixture
  RESULT_VARIABLE bench_result
  OUTPUT_VARIABLE bench_stdout
  ERROR_VARIABLE bench_stderr)
if(NOT bench_result EQUAL 0)
  message(FATAL_ERROR "generated fixture bench failed: ${bench_stderr}")
endif()

string(FIND "${bench_stdout}" "\"query_runs\":24" query_runs_pos)
if(query_runs_pos EQUAL -1)
  message(FATAL_ERROR "bench output missing expected query_runs: ${bench_stdout}")
endif()

if(NOT EXISTS "${workdir}/query_metrics.jsonl")
  message(FATAL_ERROR "query metrics output was not written")
endif()
