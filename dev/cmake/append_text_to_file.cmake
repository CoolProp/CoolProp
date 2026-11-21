# Append the contents of DATA_FILE to TARGET_FILE. Run as a standalone
# `cmake -P` process (see CMakeLists.txt's Octave ~/.octaverc registration)
# so that a write failure here raises only this process's own exit code
# instead of aborting the parent `cmake --install` run: `cmake -E touch`
# cannot be used as a write-access probe for that purpose, since on POSIX
# a file's owner can update its mtime without write permission on its
# contents, so the probe can report success while this append still fails.
if(NOT DEFINED TARGET_FILE OR NOT DEFINED DATA_FILE)
  message(FATAL_ERROR "append_text_to_file.cmake requires -DTARGET_FILE=... -DDATA_FILE=...")
endif()
file(READ "${DATA_FILE}" contents)
file(APPEND "${TARGET_FILE}" "${contents}")
