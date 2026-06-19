#
# PROJECT_NAME='daqana'
#
message("--- daqana/ana: PROJECT_NAME: ${PROJECT_NAME}")

set(subpkg "ana"                  )
set(dict   "${PROJECT_NAME}_${subpkg}_dict")
set(lib    "${PROJECT_NAME}_${subpkg}"     )

include_directories(
  $ENV{SPACK_VIEW}/include
  $ENV{ROOTSYS}/include/root
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
  $ENV{SPACK_ENV}/${PROJECT_NAME} 
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  plot_n001_sd.hh booking.hh 
  LINKDEF ana/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  plot_n001_sd.cc
  ${dict}.cxx
)

target_link_libraries(${lib}   daqana::obj Offline::Mu2eUtilities ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${PROJECT_NAME})
