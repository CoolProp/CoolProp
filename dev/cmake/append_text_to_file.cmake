# Append LINE (preceded by SEP, usually empty or a single newline) to
# TARGET_FILE. Run as a standalone `cmake -P` process (see CMakeLists.txt's
# Octave ~/.octaverc registration) so that a write failure here raises only
# this process's own exit code instead of aborting the parent
# `cmake --install` run: `cmake -E touch` cannot be used as a write-access
# probe for that purpose, since on POSIX a file's owner can update its
# mtime without write permission on its contents, so the probe can report
# success while this append still fails. Everything that can fail,
# including the write itself, happens in this child process; the parent
# passes LINE as a plain `-D` argument (no shell is involved, so no
# quoting/escaping concerns) instead of writing it to an intermediate file,
# which would reintroduce the same abort-the-install risk this exists to
# avoid.
if(NOT DEFINED TARGET_FILE OR NOT DEFINED LINE)
  message(FATAL_ERROR "append_text_to_file.cmake requires -DTARGET_FILE=... -DLINE=...")
endif()
if(NOT DEFINED SEP)
  set(SEP "")
endif()
file(APPEND "${TARGET_FILE}" "${SEP}${LINE}\n")
