set(PACKAGE_VERSION "1.0")
project(agi)

set(AGI_SOURCES
  ngraph.cpp
  agiBalancer.cpp
  migrate.cpp
  )

set(AGI_HEADERS
  ngraph.h
  )

add_library(agi ${AGI_SOURCES})

engpar_export_lib(agi "${AGI_HEADERS}")
