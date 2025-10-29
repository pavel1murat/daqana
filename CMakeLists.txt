# ======================================================================
#  daqana
# ======================================================================
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)
find_package(cetmodules)

project(daqana VERSION 25.04.21)

include(CetCMakeEnv)
cet_cmake_env()
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

cet_set_compiler_flags(
        DIAGS VIGILANT
        # WERROR
        NO_UNDEFINED
        EXTRA_FLAGS -fpermissive -Wno-unused-parameter -pthread -g  -O0 # -O0 -v Equivalent to debug qualifier, for testing
)
#------------------------------------------------------------------------------
# it looks that these packages are looked for only after cet_cmake_env
#------------------------------------------------------------------------------
find_package(Offline           10.00.00 REQUIRED)
find_package(Stntuple          25.04.21 REQUIRED)

message("---------- daqana: PRINT CMAKE VARIABLES " )
# 
# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#   message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
  
# if (DEFINED ENV{MIDASSYS})
#    message(-- P.Murat: otsdaq-mu2e-tracker: MIDASSYS=.$ENV{MIDASSYS}. )
#    find_package(Midas REQUIRED)
# else()
#    message("-- P.Murat: no MIDASSYS")
# endif()

include_directories($ENV{BTRK_INC} $ENV{KINKAL_INC})
link_directories($ENV{BTRK_LIB} $ENV{KINKAL_LIB} $ENV{SPACK_VIEW}/lib)

#string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
#if( ${BTYPE_UC} MATCHES "DEBUG" )
#  cet_add_compiler_flags(-fsanitize=address)
#endif()

#cet_report_compiler_flags()

#string(REGEX MATCHALL "[0-9][0-9]?[0-9]?" ART_VERSION_LIST "$ENV{ART_VERSION}")
#list(GET ART_VERSION_LIST 2 ART_PATCH_VERSION)

#message("-DART_HEX_VERSION=0x${ART_MAJOR_VERSION}${ART_MINOR_VERSION}${ART_PATCH_VERSION}")
#set(ART_HEX_VERSION 0x${ART_MAJOR_VERSION}${ART_MINOR_VERSION}${ART_PATCH_VERSION})
#add_definitions(-DART_HEX_VERSION=0x${ART_MAJOR_VERSION}${ART_MINOR_VERSION}${ART_PATCH_VERSION})

include(ArtDictionary)
include(BuildPlugins)

# Add a custom target to run a script
add_custom_target(create_links
    COMMAND scripts/create_links.sh
#    DEPENDS my_app
)


# Code
add_subdirectory(mod)
add_subdirectory(obj)
add_subdirectory(test)
#add_subdirectory(ana)

cet_cmake_config()
