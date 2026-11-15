#
# PROJECT_NAME='daqana'
#
# objects written to art files; the ROOT dictionary is generated from classes_def.xml
#
message("--- daqana/art_obj: PROJECT_NAME: ${PROJECT_NAME}")

set(subpkg "art_obj"                       )
set(dict   "${PROJECT_NAME}_${subpkg}_dict")
set(lib    "${PROJECT_NAME}_${subpkg}"     )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
  ${PROJECT_SOURCE_DIR}
)

add_library(${lib} SHARED
  TrkSegment.cc
)
add_library(daqana::art_obj ALIAS ${lib})

target_link_libraries(${lib} PUBLIC
  canvas::canvas
  Offline::RecoDataProducts
  Offline::TrackerGeom
  Offline::Mu2eUtilities
  ROOT::Core
)

art_dictionary( NO_CHECK_CLASS_VERSION NO_EXPORT
  CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
  CLASSES_H       ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
  LIB_TARGET      ${dict}
  DICTIONARY_LIBRARIES ${lib}
)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib} ${dict})
install(DIRECTORY ${subpkg} DESTINATION include/${PROJECT_NAME})
