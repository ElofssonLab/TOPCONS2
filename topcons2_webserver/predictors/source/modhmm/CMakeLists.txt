cmake_minimum_required(VERSION 2.6.0)
project(modhmm)


ENABLE_TESTING()

SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})

OPTION(BUILD_DOCBOOK "build Docbook html documentation" OFF )
OPTION(WITH_LIBXML "build xml functionality. This adds a dependeny for libxml" OFF )
OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle. This is  will download gengetopt, lixml2, zlib and build them. The libxml2-library ( and therefore indirectly zlib ) will be statically linked in" OFF )

SET(PERMS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE  WORLD_READ WORLD_EXECUTE )



macro(my_find_program progname)
  find_program(${progname}_EXECUTABLE ${progname})
  IF(NOT ${progname}_EXECUTABLE)
    message(FATAL_ERROR  "${progname} not found" )
  ENDIF(NOT ${progname}_EXECUTABLE)
endmacro(my_find_program)

my_find_program(xsltproc)
my_find_program(gengetopt)


add_subdirectory(src/c)
add_subdirectory(src/docbook)
add_subdirectory(src/data)
add_subdirectory(src/shellscripts)
add_subdirectory(src/xslt)






# after running cmake, we can generate source tarballs ( in tar.gz and zip ) by running
# make package_source
#

#SET(CPACK_GENERATOR RPM;DEB)

SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libxml2 (>= 2.6.16 )")


SET(CPACK_SOURCE_GENERATOR TGZ;ZIP)
# SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
SET(CPACK_PACKAGE_CONTACT "erik.sjolund@sbc.su.se") 
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "modhmm is a software package for building, training and scoring hidden Markov models (HMMs)")
SET(CPACK_PACKAGE_VENDOR "Hakan Viklund, Erik Sjölund")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
SET(CPACK_RPM_PACKAGE_LICENSE "GPL")
SET(CPACK_RPM_PACKAGE_RELEASE "1")


SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
INCLUDE(CPack)
