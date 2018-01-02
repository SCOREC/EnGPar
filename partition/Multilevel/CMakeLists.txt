if (ENABLE_PARMETIS)
  add_definitions(-DHAS_PARMETIS)
endif()


set(PARMETIS_SOURCE
  engpar_parmetis.cpp
  )


add_library(multilevel ${PARMETIS_SOURCE})

target_link_libraries(multilevel PUBLIC agi)

target_include_directories(multilevel PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_PARMETIS)
  target_include_directories(multilevel PRIVATE
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(multilevel PUBLIC
    ${PARMETIS_LIBRARIES}
    )
endif()
