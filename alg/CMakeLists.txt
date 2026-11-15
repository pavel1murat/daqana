#
# PROJECT_NAME='daqana'
#
# algorithms working on the objects of daqana/art_obj
#
message("--- daqana/alg: PROJECT_NAME: ${PROJECT_NAME}")

set(subpkg "alg"                          )
set(dict   "${PROJECT_NAME}_${subpkg}_dict")
set(lib    "${PROJECT_NAME}_${subpkg}"    )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
  ${PROJECT_SOURCE_DIR}
)

ROOT_GENERATE_DICTIONARY(${dict}
  SegmentFit.hh
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  SegmentFit.cc
  ${dict}.cxx
)
add_library(daqana::alg ALIAS ${lib})

target_link_libraries(${lib} PUBLIC daqana_art_obj Offline::Mu2eUtilities ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
install(DIRECTORY ${subpkg} DESTINATION include/${PROJECT_NAME})
