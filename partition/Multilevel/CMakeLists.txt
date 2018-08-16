option(ENABLE_PARMETIS "Enable ParMETIS iterface [ON|OFF]" OFF)
message(STATUS "ENABLE_PARMETIS: " ${ENABLE_PARMETIS})
if (ENABLE_PARMETIS)
  add_definitions(-DHAS_PARMETIS)
endif()

option(ENABLE_ZOLTAN "Enable Zoltan iterface [ON|OFF]" OFF)
message(STATUS "ENABLE_ZOLTAN: " ${ENABLE_ZOLTAN})
if (ENABLE_ZOLTAN)
  add_definitions(-DHAS_ZOLTAN)
endif()


set(PARMETIS_SOURCE
  engpar_split.cpp
  engpar_parmetis.cpp
  engpar_zoltan.cpp
  )
set(HEADERS
  engpar_split_input.h
  )

add_library(multilevel ${PARMETIS_SOURCE})

target_link_libraries(multilevel PUBLIC engpar_metrics)

target_include_directories(multilevel PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_PARMETIS)
  find_package(Parmetis MODULE REQUIRED)
  target_compile_definitions(multilevel INTERFACE -DHAS_PARMETIS)
  target_include_directories(multilevel PRIVATE
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(multilevel PUBLIC
    ${PARMETIS_LIBRARIES}
    )
endif()

if (ENABLE_ZOLTAN)
  find_package(Zoltan MODULE REQUIRED)
  target_compile_definitions(multilevel INTERFACE -DHAS_ZOLTAN)
  target_include_directories(multilevel PRIVATE
    ${ZOLTAN_INCLUDE_DIRS}
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(multilevel PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    )
endif()

engpar_export_lib(multilevel ${HEADERS})

bob_end_subdir()
