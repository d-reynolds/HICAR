# git_clone_retry.cmake — resilient git clone for flaky upstream hosts.
#
# Used as a DOWNLOAD_COMMAND from FindSNOWPACK.cmake. CMake's built-in git
# downloader retries the clone 3 times back-to-back; when the upstream host
# (git.wsl.ch) is unreachable for a few minutes those 3 attempts all time out
# together and the whole build fails. This wrapper instead retries with a
# growing backoff so attempts are spread over minutes and a short outage can
# recover between tries.
#
# Required -D args: GIT_EXECUTABLE, REPO, TAG, DEST
# Optional -D args: RETRIES (default 5), DEPTH (default 1)

if(NOT DEFINED RETRIES)
  set(RETRIES 5)
endif()
if(NOT DEFINED DEPTH)
  set(DEPTH 1)
endif()

# Already populated (e.g. an incremental reconfigure) -> nothing to do.
if(EXISTS "${DEST}/.git")
  message(STATUS "git_clone_retry: ${DEST} already cloned, skipping")
  return()
endif()

set(_ok FALSE)
set(_attempt 1)
while(NOT _ok AND _attempt LESS_EQUAL ${RETRIES})
  # A previous failed attempt may have left a partial directory behind; git
  # clone requires the target to be absent or empty.
  file(REMOVE_RECURSE "${DEST}")
  message(STATUS "git_clone_retry: cloning ${REPO} (${TAG}) -> ${DEST} [attempt ${_attempt}/${RETRIES}]")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" clone --branch "${TAG}" --single-branch
            --depth ${DEPTH} --recurse-submodules "${REPO}" "${DEST}"
    RESULT_VARIABLE _rc
  )
  if(_rc EQUAL 0)
    set(_ok TRUE)
  else()
    if(_attempt LESS ${RETRIES})
      # Backoff: 15s, 30s, 45s, ... — spreads attempts over minutes so a brief
      # host outage recovers between tries (unlike CMake's back-to-back retries).
      math(EXPR _sleep "15 * ${_attempt}")
      message(STATUS "git_clone_retry: attempt ${_attempt} failed (rc=${_rc}); retrying in ${_sleep}s")
      execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep ${_sleep})
    endif()
    math(EXPR _attempt "${_attempt} + 1")
  endif()
endwhile()

if(NOT _ok)
  message(FATAL_ERROR "git_clone_retry: failed to clone ${REPO} (${TAG}) after ${RETRIES} attempts")
endif()
