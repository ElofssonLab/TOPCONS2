IF(UNIX)	
	IF(CMAKE_COMPILER_IS_GNUCC)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
	ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF(UNIX)





IF(STATIC)
 SET(LIBXML2_INCLUDE_DIR  "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2" )
 SET(LIBXML2_LIBRARIES  "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a" )

# The next line seems to be need for macosx
 SET(ZLIB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" )
ENDIF(STATIC)

SET(LIBXML2_VERSION 2.7.2 )
SET(ZLIB_VERSION 1.2.3 )

SET(DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/downloads" INTERNAL )
SET(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/build INTERNAL )
SET(INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/install INTERNAL )

file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/download )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/build )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/install )

ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_CURRENT_BINARY_DIR}/script.nsi" DEPENDS fastdist fnj )
ADD_CUSTOM_TARGET( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" )
ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" COMMAND  wget --directory-prefix=external/download  "ftp://xmlsoft.org/libxml2/libxml2-${LIBXML2_VERSION}.tar.gz" )
ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz" COMMAND  wget --directory-prefix=external/download  "http://www.gzip.org/zlib/zlib-${ZLIB_VERSION}.tar.gz" )

ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a" "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_libxml2.sh" ${LIBXML2_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install"  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a"  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build" )
ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_zlib.sh" ${ZLIB_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install"  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz"  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build" )


set(commentout_xml_on_stdin_ggo "#")
IF(WITH_LIBXML)
  include(FindLibXml2)
  IF(NOT LIBXML2_FOUND)
  message(FATAL_ERROR  "libxml2 devel not found" )
  ENDIF(NOT LIBXML2_FOUND)

  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  set(commentout_xml_on_stdin_ggo "")
ENDIF(WITH_LIBXML)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/script.nsi.cmake ${CMAKE_CURRENT_BINARY_DIR}/script.nsi)

INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)

IF(STATIC AND WITH_LIBXML)
 set( LIBXML_EXTRA_COMPILE_FLAGS " -DLIBXML_STATIC " )
ENDIF(STATIC AND WITH_LIBXML)


IF(WITH_LIBXML)
 IF(STATIC)
  SET(FNJ_XML_SRCS programs/fnj/XmlInputStream.cpp  "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" )
 ELSE(STATIC)
  SET(FNJ_XML_SRCS programs/fnj/XmlInputStream.cpp )
 ENDIF(STATIC)
ENDIF(WITH_LIBXML)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/${LIBXML2_INCLUDE_DIR}/libxml/xmlreader.h PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )

SET (COMMON_H_FILES funcs.h
                     structs.h)
SET (COMMON_C_FILES readhmm.c
                     readhmm_multialpha.c
                     std_funcs.c
                     savehmm.c
                     savehmm_multialpha.c
                     debug_funcs.c
                     readseqs_multialpha.c)
SET (SCORING_C_FILES core_algorithms_multialpha.c
                     std_calculation_funcs.c)

ADD_LIBRARY(commonc STATIC ${COMMON_C_FILES})


IF(NOT TARGETS)
 SET(TARGETS general octopus spoctopus scampi)
ENDIF(NOT TARGETS)

SET(PROGRAMS modhmms modhmmt modseqalign add_alphabet) 

FOREACH(j ${TARGETS} )
 if(${j} STREQUAL "general")
  set(suffix "")
  set(commentout_userscore_ggo "#")
 else(${j} STREQUAL "general")
  set(suffix "_${j}")
  set(commentout_userscore_ggo "")
 endif(${j} STREQUAL "general")
 if(${j} STREQUAL "scampi")
  set(commentout_globmemscan_ggo "")
 else(${j} STREQUAL "scampi")
  set(commentout_globmemscan_ggo "#")
 endif(${j} STREQUAL "scampi")

 ADD_LIBRARY(scoring${suffix} STATIC ${SCORING_C_FILES})
 SET_TARGET_PROPERTIES(scoring${suffix} PROPERTIES COMPILE_FLAGS "-DTARGET${j}")

 FOREACH(i ${PROGRAMS} )
  set(currenttarget ${i})
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/gengetopt/cmdline_${i}.ggo.cmake ${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.ggo @ONLY)
  add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.c   ${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.h 
                     COMMAND ${gengetopt_EXECUTABLE} -i ${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.ggo --file-name=cmdline_${i} --output-dir=${CMAKE_CURRENT_BINARY_DIR}/${j}
                     DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.ggo  
                     VERBATIM)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.c PROPERTIES GENERATED true)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/${j}/cmdline_${i}.h PROPERTIES GENERATED true)
 ENDFOREACH(i)

 add_executable (modhmms${suffix} modhmms.c ${j}/cmdline_modhmms.c  training_algorithms_multialpha.c)
 add_executable (modhmmt${suffix}  modhmmt.c ${j}/cmdline_modhmmt.c training_algorithms_multialpha.c)
 add_executable (modseqalign${suffix} modseqalign.c ${j}/cmdline_modseqalign.c )
 add_executable (add_alphabet${suffix} add_alphabet.c ${j}/cmdline_add_alphabet.c )

 FOREACH(i ${PROGRAMS} )
  SET_TARGET_PROPERTIES(${i}${suffix} PROPERTIES COMPILE_FLAGS "-DTARGET${j} -I${CMAKE_CURRENT_BINARY_DIR}/${j}" )
 ENDFOREACH(i)

 FOREACH(i ${PROGRAMS} )
  IF(STATIC)
   set_property(TARGET ${i}${suffix} PROPERTY LINK_SEARCH_END_STATIC ON )
  ENDIF(STATIC)
  target_link_libraries (${i}${suffix} commonc m ${LIBXML2_LIBRARIES})
 ENDFOREACH(i)
 target_link_libraries (modhmms${suffix} scoring${suffix})
 target_link_libraries (modhmmt${suffix} scoring${suffix})
 target_link_libraries (modseqalign${suffix} scoring${suffix})
 FOREACH(i ${PROGRAMS} )
  install(TARGETS ${i}${suffix} DESTINATION bin)
  get_target_property(${i}${suffix}_exe ${i}${suffix} LOCATION)
 ENDFOREACH(i)
 file(GLOB testdir test/${j}/[^.]*)
 FOREACH(k ${testdir} )
   add_subdirectory(${k})
 ENDFOREACH(k ${testdir} )
ENDFOREACH(j)



