execute_process(
  COMMAND "${AIMS_BENCH_EXE}"
    --ref "${AIMS_REF}"
    --query "${AIMS_QUERY}"
    --truth "${AIMS_TRUTH}"
    --k 5,7
    --topk 3
    --dataset synthetic_kmer_fixture
  RESULT_VARIABLE bench_result
  OUTPUT_VARIABLE bench_stdout
  ERROR_VARIABLE bench_stderr)

if(NOT bench_result EQUAL 0)
  message(FATAL_ERROR "aims_bench failed: ${bench_stderr}")
endif()

string(JSON stage ERROR_VARIABLE json_error GET "${bench_stdout}" comparison_stage)
if(json_error OR NOT stage STREQUAL "exact_retrieval")
  message(FATAL_ERROR "invalid benchmark JSON/stage: ${json_error}: ${bench_stdout}")
endif()

string(JSON cache_mode ERROR_VARIABLE cache_error GET "${bench_stdout}" cache_mode)
if(cache_error OR NOT cache_mode STREQUAL "disabled")
  message(FATAL_ERROR "unexpected benchmark cache mode: ${cache_error}: ${bench_stdout}")
endif()

string(JSON top1_recall ERROR_VARIABLE recall_error GET "${bench_stdout}" metrics top1_recall)
if(recall_error OR top1_recall LESS 0 OR top1_recall GREATER 1)
  message(FATAL_ERROR "invalid benchmark recall: ${recall_error}: ${bench_stdout}")
endif()
