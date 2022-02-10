target_include_directories(arlequin PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MATERIAL_HEADERS
    TPZMatPoissonCS.h
  )

set(MATERIAL_SRCS
    TPZMatPoissonCS.cpp
  )


set(MESH_HEADERS

  )

set(MESH_SRCS

  )

# let us create a config file with the directory containing the meshes
# create folder with meshes
set(MESHESDIR ${PROJECT_SOURCE_DIR}/mesh)
set (CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/arlequin_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/arlequin_config.h.in" ${CONFIG_H})

set(UTIL_HEADERS
  
  )

set(UTIL_SRCS
  
  )

target_sources(arlequin PRIVATE
  ${MATERIAL_HEADERS} ${MATERIAL_SRCS}
  ${MESH_HEADERS} ${MESH_SRCS}
  ${UTIL_HEADERS} ${UTIL_SRCS})