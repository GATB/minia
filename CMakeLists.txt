project(minia)

cmake_minimum_required (VERSION 2.6)

################################################################################
# The version number.
################################################################################
# The default version number is the latest official build
SET (gatb-tool_VERSION_MAJOR 3)
SET (gatb-tool_VERSION_MINOR 2)
SET (gatb-tool_VERSION_PATCH 1)

# But, it is possible to define another release number during a local build
IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

#############################
# getting git version
#from http://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake
exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "rev-parse --short HEAD"
    OUTPUT_VARIABLE VERSION_SHA1 )

if ( VERSION_SHA1 MATCHES "^[0-9a-f]+$") # might return "git commit fatal: not a git repository"
    message("Git version: ${VERSION_SHA1}")
    add_definitions( -DGIT_SHA1="${VERSION_SHA1}" )
endif()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")
# by default, it's compiling up to 128.

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*)
add_executable(${PROJECT_NAME} ${ProjectFiles})
target_link_libraries(${PROJECT_NAME} ${gatb-core-libraries})

# merci
set (MERCI_SOURCE_DIR ${PROJECT_SOURCE_DIR}/merci)
include_directories (${MERCI_SOURCE_DIR})
file (GLOB_RECURSE  MerciFiles  ${MERCI_SOURCE_DIR}/*.cpp)
add_executable("merci"  ${MerciFiles})
target_link_libraries("merci" ${gatb-core-libraries})


# that file will hold the version information as a #define
configure_file (
${PROJECT_SOURCE_DIR}/src/build_info.hpp.in
${PROJECT_SOURCE_DIR}/src/build_info.hpp
)

################################################################################
#  DOC
################################################################################
set (SKIP_DOC 1)
if (NOT DEFINED SKIP_DOC)
    add_subdirectory (doc)
endif()

################################################################################
#  PACKAGING
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwanted files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES
    "^${PROJECT_SOURCE_DIR}/\\.git/"     ;
    "^${PROJECT_SOURCE_DIR}/\\.gitmodules" ;
    "^${PROJECT_SOURCE_DIR}/\\.gitignore" ;
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/\\.cproject" ;
    "^${GATB_CORE_HOME}/\\.git/"     ;
    "^${GATB_CORE_HOME}/\\.project"  ;
    "^${GATB_CORE_HOME}/\\.gitignore"
)

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
INSTALL (DIRECTORY "${PROJECT_SOURCE_DIR}/test/" DESTINATION test)
INSTALL (FILES LICENSE README.md DESTINATION bin/..)

# We include the "bin" tag into binary archive file name
set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME})

include (CPack)
