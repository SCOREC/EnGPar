
set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")


include(GetGitRevisionDescription)
get_git_head_revision(ENGPAR_GIT_REFSPEC ENGPAR_GIT_SHA1)
configure_file(engpar_version.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/engpar_version.cpp)

set(SOURCES
  engpar_support.cpp
  engpar_log.cpp
  engpar_message.cpp
  engpar_version.cpp
  engpar_memory.cpp
  )

set(SUPPORT_HEADERS
  engpar_support.h
  engpar_version.h
  )

add_library(engpar_support ${SOURCES})

target_include_directories(engpar_support
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

engpar_export_lib(engpar_support "${SUPPORT_HEADERS}")

bob_end_subdir()
