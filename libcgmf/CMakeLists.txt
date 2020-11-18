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

add_subdirectory(src)

get_filename_component(ABS_INSTALL_PATH ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
set(BUILD_DATADIR ${CMAKE_CURRENT_SOURCE_DIR}/../data)
set(INSTALL_DATADIR ${ABS_INSTALL_PATH}/cgmf/data)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cgmf_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cgmf_config.h)

include_directories(${CMAKE_CURRENT_LIST_DIR} include)

install(DIRECTORY include DESTINATION cgmf)
install(TARGETS cgmf EXPORT cgmf DESTINATION lib)
