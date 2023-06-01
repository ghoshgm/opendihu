
find_path(BASE64_INCLUDE_DIRS NAMES base64.h)

if(EXISTS ${BASE64_INCLUDE_DIRS})
  message(STATUS "Found Base64: ${BASE64_INCLUDE_DIRS}/base64.h")
endif()

add_library(
  base64 SHARED 
  ${BASE64_INCLUDE_DIRS}/base64.h
)

set_target_properties(base64 PROPERTIES LINKER_LANGUAGE CXX)

set(BASE64_LIBRARIES "base64")