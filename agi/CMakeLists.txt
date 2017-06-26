set(PACKAGE_VERSION "1.0")
project(agi)

include_directories(
  ../support
  )
set(AGI_SOURCES
  ngraph.cpp
  agiIterate.cpp
  agiBalancer.cpp
  agiMigrate.cpp
  agiIO.cpp
  )

set(AGI_HEADERS
  ngraph.h
  )

add_library(agi ${AGI_SOURCES})

target_link_libraries(agi PUBLIC engpar_support)
engpar_export_lib(agi "${AGI_HEADERS}")
