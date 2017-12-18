option(ENABLE_PARMETIS "Enable ParMETIS interface [ON|OFF]" OFF)
message(STATUS "ENABLE_PARMETIS: " ${ENABLE_PARMETIS})

add_subdirectory(Diffusive)

add_subdirectory(Multilevel)

if (ENABLE_PARMETIS)
  find_package(Parmetis MODULE REQUIRED)
endif()

set(ENGPAR_SOURCES
  engpar_split.cpp
  engpar_input.cpp
  )

set(ENGPAR_HEADERS
  engpar_split.h
  )

add_library(engpar ${ENGPAR_SOURCES})

target_include_directories(engpar PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(engpar PUBLIC agi diffusive multilevel)


if (ENABLE_PARMETIS)
  message(${PARMETIS_INCLUDE_DIRS})
  target_include_directories(engpar PRIVATE
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(engpar PUBLIC
    ${PARMETIS_LIBRARIES}
    )
endif()

engpar_export_lib(engpar "${ENGPAR_HEADERS}")
