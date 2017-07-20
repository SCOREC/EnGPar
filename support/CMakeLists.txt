
#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
#project(EnGPar)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")


include(GetGitRevisionDescription)
get_git_head_revision(ENGPAR_GIT_REFSPEC ENGPAR_GIT_SHA1)
configure_file(engpar_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/engpar_version.h)

set(SOURCES
  engpar_support.cpp
  )

set(SUPPORT_HEADERS
  engpar_support.h
  )

add_library(engpar_support ${SOURCES})

engpar_export_lib(engpar_support "${SUPPORT_HEADERS}")
