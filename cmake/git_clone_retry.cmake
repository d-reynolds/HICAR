# git_clone_retry.cmake — resilient git clone/update for the flaky upstream host
# (git.wsl.ch). Used as the DOWNLOAD_COMMAND for SNOWPACK and MeteoIO.
#
# Two design goals:
#   1. Keep tracking the TIP of a moving branch TAG (the project wants HEAD on
#      `master` / `fortran-bindings`, not a pinned commit).
#   2. A git.wsl.ch stall or outage must NOT fail the build once a previously
#      fetched copy is available to fall back on.
#
# How the no-fail guarantee works — a persistent FALLBACK copy:
#   FALLBACK (optional) is a STABLE directory outside the build tree (CI keeps it
#   warm with actions/cache; see .github/actions/snowpack-src-cache). Each run:
#     * if DEST is absent but FALLBACK exists, DEST is seeded from FALLBACK;
#     * DEST is then UPDATED to the tip of TAG (so we track HEAD) — or, if the
#       host is unreachable, we BUILD AGAINST the seeded copy instead of failing;
#     * after a good DEST is obtained, FALLBACK is refreshed to it (latest-good).
#   Keeping FALLBACK in a stable home dir (not the build dir) means it survives
#   `rm -rf build/*` reconfigures and ephemeral CI runners, and needs no change
#   to the build scripts.
#
# Only a COLD build (no DEST and no FALLBACK) can still hard-fail — there is
# nothing to fall back to yet; the first successful build seeds FALLBACK.
#
# Network hardening: a stalled transfer (< ~1 KB/s for LOW_SPEED_TIME s) is
# aborted so the backoff/retry kicks in promptly instead of the job hanging to
# its global timeout; each git invocation also has a hard wall-clock TIMEOUT.
#
# Required -D args: GIT_EXECUTABLE, REPO, TAG, DEST
# Optional -D args: FALLBACK, RETRIES (default 8), DEPTH (default 1),
#                   ATTEMPT_TIMEOUT (s, default 600), LOW_SPEED_TIME (s, default 30)

if(NOT DEFINED RETRIES)
  set(RETRIES 8)
endif()
if(NOT DEFINED DEPTH)
  set(DEPTH 1)
endif()
if(NOT DEFINED ATTEMPT_TIMEOUT)
  set(ATTEMPT_TIMEOUT 600)
endif()
if(NOT DEFINED LOW_SPEED_TIME)
  set(LOW_SPEED_TIME 30)
endif()

# git -c flags: abort a hung transfer (< 1000 B/s sustained for LOW_SPEED_TIME)
# instead of blocking forever, so a stalled connection fails fast and retries.
set(_net -c "http.lowSpeedLimit=1000" -c "http.lowSpeedTime=${LOW_SPEED_TIME}")

# --- Seed DEST from the persistent fallback, if we have one and DEST is cold. --
if(DEFINED FALLBACK AND EXISTS "${FALLBACK}/.git" AND NOT EXISTS "${DEST}/.git")
  message(STATUS "git_clone_retry: seeding ${DEST} from fallback ${FALLBACK}")
  file(REMOVE_RECURSE "${DEST}")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_directory "${FALLBACK}" "${DEST}")
endif()

set(_have_dest FALSE)

if(EXISTS "${DEST}/.git")
  # ----- Update path: bring the existing copy to the tip of TAG; on a network
  #       failure, keep the existing copy rather than failing the build. -------
  message(STATUS "git_clone_retry: ${DEST} present — updating to tip of ${TAG}")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" ${_net} -C "${DEST}"
            fetch --force --depth ${DEPTH} origin "${TAG}"
    RESULT_VARIABLE _rc
    TIMEOUT ${ATTEMPT_TIMEOUT})
  if(_rc EQUAL 0)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${DEST}" reset --hard FETCH_HEAD)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" ${_net} -C "${DEST}"
              submodule update --init --recursive --depth ${DEPTH}
      TIMEOUT ${ATTEMPT_TIMEOUT})
    message(STATUS "git_clone_retry: ${DEST} updated to tip of ${TAG}")
  else()
    message(WARNING
      "git_clone_retry: could not reach ${REPO} (rc=${_rc}). "
      "Building against the EXISTING copy in ${DEST} (possibly slightly stale) "
      "instead of failing the build.")
  endif()
  set(_have_dest TRUE)
else()
  # ----- Clone path: no copy yet. Retry with a growing backoff so a multi-minute
  #       outage can recover between attempts. Only this path can hard-fail. ----
  set(_attempt 1)
  while(NOT _have_dest AND _attempt LESS_EQUAL ${RETRIES})
    # A previous failed attempt may have left a partial directory; git clone
    # requires the target to be absent or empty.
    file(REMOVE_RECURSE "${DEST}")
    message(STATUS "git_clone_retry: cloning ${REPO} (${TAG}) -> ${DEST} [attempt ${_attempt}/${RETRIES}]")
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" ${_net} clone --branch "${TAG}" --single-branch
              --depth ${DEPTH} --recurse-submodules "${REPO}" "${DEST}"
      RESULT_VARIABLE _rc
      TIMEOUT ${ATTEMPT_TIMEOUT})
    if(_rc EQUAL 0)
      set(_have_dest TRUE)
    else()
      if(_attempt LESS ${RETRIES})
        # Backoff 15s, 30s, 45s, ... spreads attempts over minutes so a brief
        # host outage can recover between tries (unlike git's back-to-back retries).
        math(EXPR _sleep "15 * ${_attempt}")
        message(STATUS "git_clone_retry: attempt ${_attempt} failed (rc=${_rc}); retrying in ${_sleep}s")
        execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep ${_sleep})
      endif()
      math(EXPR _attempt "${_attempt} + 1")
    endif()
  endwhile()
endif()

if(NOT _have_dest)
  message(FATAL_ERROR
    "git_clone_retry: failed to clone ${REPO} (${TAG}) after ${RETRIES} attempts, "
    "and no fallback copy was available. (Once one build succeeds, CI caches a "
    "fallback so later outages no longer fail the build.)")
endif()

# --- Refresh the persistent fallback to the copy we ended up with (HEAD when
#     the host was reachable, otherwise the last-good copy). ------------------
if(DEFINED FALLBACK AND EXISTS "${DEST}/.git")
  file(REMOVE_RECURSE "${FALLBACK}")
  get_filename_component(_fb_parent "${FALLBACK}" DIRECTORY)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E make_directory "${_fb_parent}")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_directory "${DEST}" "${FALLBACK}")
endif()
