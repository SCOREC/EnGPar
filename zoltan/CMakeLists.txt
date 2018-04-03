set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

if (ENABLE_ZOLTAN)
  set(ZOLTAN_SOURCE
    ZoltanCallbacks.cpp
    ZoltanCutVertex.cpp
    )
else()
  set(ZOLTAN_SOURCE
    ZoltanShell.cpp
    )
endif()


add_library(agi_zoltan ${ZOLTAN_SOURCE})

target_link_libraries(agi_zoltan PUBLIC agi)

if (ENABLE_ZOLTAN)
  target_include_directories(agi_zoltan PRIVATE
    ${ZOLTAN_INCLUDE_DIRS}
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(agi_zoltan PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    )
endif()

engpar_export_lib(agi_zoltan "")
