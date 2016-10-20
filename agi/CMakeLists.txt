set(PACKAGE_VERSION "1.0")
project(agi)

set(AGI_SOURCES
  agi.cpp
  )

set(AGI_HEADERS
  agi.h
  )

add_library(agi ${AGI_SOURCES})
