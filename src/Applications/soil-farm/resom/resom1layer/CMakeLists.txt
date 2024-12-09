set(resom1LAYER_SOURCES
  BgcresomIndexType.F90
  BgcresomCompetType.F90
  BgcresomDecompType.F90
  BgcresomNitDenType.F90
  BgcresomSOMType.F90
  BgcresomType.F90
  BgcresomMath.F90
  BgcresomDebType.F90
  DebGrowMod.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/resom/resomPara)                
  include(add_betr_library)
  add_betr_library(resom1layer ${resom1LAYER_SOURCES})

  set(BETR_LIBRARIES resom1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES resom1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS resom1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/resom/resom1layer)
  endif()

