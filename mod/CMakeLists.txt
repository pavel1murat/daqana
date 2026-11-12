# PM: use cet
set (subpkg "mod"                  )
set (dict   "daqana_${subpkg}_dict")
set (lib    "daqana_${subpkg}"     )

include_directories(
  $ENV{ROOTSYS}/include/root
  $ENV{SPACK_VIEW}/include
  ${CMAKE_SOURCE_DIR}/
  ${CMAKE_CURRENT_SOURCE_DIR}/${subpkg}
)

# link_directories($ENV{ROOTSYS}/lib)

set (include_files
  PSegFinder_module.hh
  PSegFinderHolder.hh
#  TrkFragmentAna_module.hh
)

ROOT_GENERATE_DICTIONARY(${dict}
  ${include_files}
  MODULE ${dict}
  LINKDEF ${subpkg}/Linkdef.h)

set_property(SOURCE ${dict}.cxx PROPERTY generated)

# dictionary of the modules listed in Linkdef.h
add_library(${lib} SHARED
  ${dict}.cxx
)
set_target_properties(${lib} PROPERTIES OUTPUT_NAME ${dict})
install(TARGETS ${lib})

# the PSegFinder plugin is a MODULE library: link it by name, so the dependency is recorded
# as libmod_PSegFinder_module.so and found through the RPATH
add_dependencies(${lib} mod_PSegFinder_module)
set_target_properties(${lib} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
target_link_libraries(${lib}
  -L$<TARGET_FILE_DIR:mod_PSegFinder_module>
  -lmod_PSegFinder_module
  daqana_obj
  Stntuple_mod
  Stntuple_print
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  Offline::TrackerGeom
  Offline::ProditionsService
  art::Framework_Core
  ROOT::Core
)

# cet_make_library(LIBRARY_NAME daqana_mod
#   SOURCE
#   TrkSegment.cc
#   LIBRARIES PUBLIC
#   ROOT::Geom
# )

cet_build_plugin(TrkFragmentAna art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(FragmentAna art::module LIBRARIES REG
  Stntuple_mod
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(FileCatalog art::module LIBRARIES REG
  Stntuple_mod
  #   Stntuple_gui
  Stntuple_print
  Stntuple_geom
  Stntuple_obj
  Stntuple_base
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(StationAna art::module LIBRARIES REG
  Stntuple_mod Stntuple_print Stntuple_geom Stntuple_obj Stntuple_base
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist ROOT::Tree ROOT::Core ROOT::RIO ROOT::Gui
)

target_include_directories(mod_StationAna_module     PRIVATE ${Stntuple_INCLUDE_DIR})
target_include_directories(mod_FileCatalog_module    PRIVATE ${Stntuple_INCLUDE_DIR})
target_include_directories(mod_TrkFragmentAna_module PRIVATE ${Stntuple_INCLUDE_DIR})
target_include_directories(mod_FragmentAna_module    PRIVATE ${Stntuple_INCLUDE_DIR})

target_link_directories(mod_StationAna_module     PRIVATE ${Stntuple_LIBRARY_DIR})
target_link_directories(mod_TrkFragmentAna_module PRIVATE ${Stntuple_LIBRARY_DIR})
target_link_directories(mod_FileCatalog_module    PRIVATE ${Stntuple_LIBRARY_DIR})
target_link_directories(mod_FragmentAna_module    PRIVATE ${Stntuple_LIBRARY_DIR})

cet_build_plugin(MakeDigiNtuple art::module LIBRARIES REG
  # daqana::mod
  daqana_obj
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist ROOT::Tree ROOT::RIO ROOT::Gui ROOT::Geom ROOT::Matrix ROOT::GenVector ROOT::Core
)

cet_build_plugin(PSegFinder art::module LIBRARIES REG
  daqana_obj
  Stntuple_mod
  Stntuple_print
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  Offline::TrackerGeom
  Offline::ProditionsService
)

cet_build_plugin(TCLineFilter art::module LIBRARIES REG
  daqana_obj
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(DaqanaFilter art::module LIBRARIES REG
  daqana_obj
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::Core
  ROOT::RIO
  ROOT::Gui
)

cet_build_plugin(MixFragments art::module LIBRARIES REG
  # daqana::mod
  daqana_obj
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art::Framework_IO_ProductMix
  art_root_io::art_root_io
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::RIO
  ROOT::Gui
  ROOT::Geom
  ROOT::Matrix
  ROOT::GenVector
  ROOT::Core
)

cet_build_plugin(FragmentNormalizeEB art::module LIBRARIES REG
  # daqana::mod
  daqana_obj
  Offline::Mu2eUtilities
  Offline::RecoDataProducts
  art::Framework_IO_ProductMix
  art_root_io::art_root_io
  art_root_io::TFileService_service
  artdaq-core-mu2e::Data
  ROOT::Hist
  ROOT::Tree
  ROOT::RIO
  ROOT::Gui
  ROOT::Geom
  ROOT::Matrix
  ROOT::GenVector
  ROOT::Core
)

cet_add_compiler_flags(-g -O0)

# install(TARGETS ${lib})
# install_headers()
# install_source()
