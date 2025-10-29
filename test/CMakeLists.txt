#
# PROJECT_NAME='daqana'
#
message("--- daqana/test: PROJECT_NAME: ${PROJECT_NAME}")

set(subpkg "test"                          )
set(dict   "${PROJECT_NAME}_${subpkg}_dict")
set(lib    "${PROJECT_NAME}_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
  $ENV{SPACK_ENV}/${PROJECT_NAME} 
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict}
  test3.hh
  LINKDEF ${subpkg}/Linkdef.h)

# set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
  test3.cc read_combohits.cc
  ${dict}.cxx
)

target_link_libraries(${lib}  daqana_obj Offline::Mu2eUtilities ROOT::Core)

cet_add_compiler_flags(-g -O0)

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/${PROJECT_NAME})
