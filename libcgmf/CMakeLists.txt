# Shared or Static Library Build
option(cgmf.shared_library
  "Build cgmf as a shared library" ${BUILD_SHARED_LIBS})

if(cgmf.shared_library)
  add_library(cgmf SHARED "")
else()
  add_library(cgmf STATIC "")
endif()

target_include_directories(cgmf PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<INSTALL_INTERFACE:include/cgmf>)

option(cgmf.devutil
  "Preprocessor DEVUTIL option" OFF)

if(cgmf.devutil)
  target_compile_definitions(cgmf PUBLIC DEVUTIL)
endif()

if(MSVC_VERSION GREATER_EQUAL 1910)
  target_compile_options(cgmf PUBLIC
    "$<$<COMPILE_LANGUAGE:CXX>:/permissive->"
  )
endif()

add_subdirectory(src)

set(BUILD_DATADIR ${CMAKE_CURRENT_SOURCE_DIR}/../data)
set(INSTALL_DATADIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/data/cgmf-${CGMF_VER})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cgmf_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cgmf_config.h)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_LIST_DIR} include)

if(cgmf.installation)
  install(DIRECTORY include/
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cgmf-${CGMF_VER}"
          FILES_MATCHING PATTERN "*.h"
  )
  install(TARGETS cgmf EXPORT cgmf 
          DESTINATION "${CMAKE_INSTALL_LIBDIR}/cgmf-${CGMF_VER}"
  )
endif()
