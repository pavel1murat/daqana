#
set(subpkg "obj"                  )
set(dict   "daqana_${subpkg}_dict")
set(lib    "daqana_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
#  $ENV{SPACK_ENV}/daqana
)

# link_directories($ENV{ROOTSYS}/lib)

ROOT_GENERATE_DICTIONARY(${dict} DaqEvent.hh DaqStrawDigi.hh LINKDEF Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

add_library(${lib} SHARED
            DaqEve nt.cc DaqStrawDigi.cc ${dict}.cxx
)

target_link_libraries(${lib} ROOT::Core)

cet_add_compiler_flags(-g -O0)

message("--- daqana/obj: PROJECT_NAME: ${PROJECT_NAME}"

install(TARGETS ${lib})
# install_headers()
# install_source()
install(DIRECTORY ${subpkg} DESTINATION include/daqana)
