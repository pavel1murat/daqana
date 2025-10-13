#
# PROJECT_NAME='daqana'
#
message("--- daqana/obj: PROJECT_NAME: ${PROJECT_NAME}")

set(subpkg "obj"                  )
set(dict   "${PROJECT_NAME}_${subpkg}_dict")
set(lib    "${PROJECT_NAME}_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
  $ENV{SPACK_ENV}/${PROJECT_NAME} 
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  DaqEvent.hh DaqTrack.hh DaqStrawDigi.hh DaqStrawHit.hh DaqComboHit.hh DaqTimeCluster.hh
  TrkSegment.hh
  LINKDEF obj/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  DaqEvent.cc DaqTrack.cc DaqStrawDigi.cc DaqStrawHit.cc DaqComboHit.cc DaqTimeCluster.cc
  TrkSegment.cc SegmentFit.cc
  ${dict}.cxx
)

target_link_libraries(${lib}   Offline::Mu2eUtilities ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${PROJECT_NAME})
