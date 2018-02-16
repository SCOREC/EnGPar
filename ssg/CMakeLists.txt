set(PACKAGE_VERSION "1.0")
project(ssg)

set(SSG_SOURCES
  ssgConvert.cpp
  )

set(SSG_HEADERS
  sellCSigma.h
  )

add_library(ssg ${SSG_SOURCES})

target_include_directories(ssg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(ssg PUBLIC engpar_support agi)

engpar_export_lib(ssg "${SSG_HEADERS}")
