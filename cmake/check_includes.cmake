# Fails if any .hpp or .cpp in DIR includes a header from a FORBIDDEN module.
# Usage: cmake -DDIR=<dir> "-DFORBIDDEN=sim;io" -P check_includes.cmake

if(NOT DIR OR NOT FORBIDDEN)
  message(FATAL_ERROR "check_includes: set DIR and FORBIDDEN")
endif()

list(JOIN FORBIDDEN "|" modules)
set(pattern "^[ \t]*#[ \t]*include[ \t]*[<\"](\\.\\.?/)*(${modules})/")

file(GLOB_RECURSE sources "${DIR}/*.hpp" "${DIR}/*.cpp")
set(violations "")
foreach(src IN LISTS sources)
  # One list element per line; semicolons in code would otherwise split lines.
  file(READ "${src}" text)
  string(REPLACE ";" "," text "${text}")
  string(REPLACE "\n" ";" lines "${text}")
  set(n 0)
  foreach(line IN LISTS lines)
    math(EXPR n "${n} + 1")
    if(line MATCHES "${pattern}")
      string(APPEND violations "\n  ${src}:${n}: ${line}")
    endif()
  endforeach()
endforeach()

if(violations)
  message(FATAL_ERROR "includes across a module boundary (${FORBIDDEN}):${violations}")
endif()
list(LENGTH sources count)
message(STATUS "${DIR}: ${count} files, no includes from ${FORBIDDEN}")
