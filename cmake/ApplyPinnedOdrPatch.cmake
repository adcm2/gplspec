# Apply an ODR-only patch to an owned FetchContent source tree. A patch is
# accepted only when the complete patch applies or the complete reverse patch
# applies (meaning it was already applied). Partial/unexpected states fail
# without writing source files.
function(gplspec_apply_odr_patch source patch build_root)
  get_filename_component(_source_real "${source}" REALPATH)
  get_filename_component(_build_real "${build_root}" REALPATH)
  set(_owned_prefix "${_build_real}/_deps/")
  string(FIND "${_source_real}/" "${_owned_prefix}" _owned_index)
  if(NOT "${_owned_index}" STREQUAL "0")
    message(FATAL_ERROR
      "Refusing to patch dependency outside build-owned FetchContent tree: ${_source_real}; expected below ${_owned_prefix}")
  endif()

  execute_process(
    COMMAND git apply --check "${patch}"
    WORKING_DIRECTORY "${_source_real}"
    RESULT_VARIABLE _forward_check
    OUTPUT_QUIET ERROR_QUIET)
  if("${_forward_check}" STREQUAL "0")
    execute_process(
      COMMAND git apply "${patch}"
      WORKING_DIRECTORY "${_source_real}"
      RESULT_VARIABLE _apply_result
      OUTPUT_QUIET ERROR_QUIET)
    if(NOT "${_apply_result}" STREQUAL "0")
      message(FATAL_ERROR
        "Failed to apply complete ODR patch ${patch} to ${_source_real} (git result: ${_apply_result})")
    endif()
    return()
  endif()

  execute_process(
    COMMAND git apply --reverse --check "${patch}"
    WORKING_DIRECTORY "${_source_real}"
    RESULT_VARIABLE _reverse_check
    OUTPUT_QUIET ERROR_QUIET)
  if("${_reverse_check}" STREQUAL "0")
    return()
  endif()

  message(FATAL_ERROR
    "Refusing partial or unexpected dependency source for ODR patch ${patch}: complete forward and reverse checks failed (git results: ${_forward_check}; ${_reverse_check}); source left unchanged")
endfunction()
