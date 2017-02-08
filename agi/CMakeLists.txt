set(PACKAGE_VERSION "1.0")
project(agi)

set(AGI_SOURCES
  ngraph.cpp
  agiBalancer.cpp
  )

set(AGI_HEADERS
  agi.h
  ngraph.h
  agiBalancer.h
  )

add_library(agi ${AGI_SOURCES})
