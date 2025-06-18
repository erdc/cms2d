cmake_minimum_required(VERSION 3.15)
project(XMDF C)

find_package(ZLIB REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(src/include)

SET(SRC_LIST
  src/Xmdf.c
  src/Xmdf_f.c
  src/xmdf_private.c
  src/xmdf_timestep.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
target_link_libraries(${PROJECT_NAME} HDF5::HDF5)

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  target_compile_options(${PROJECT_NAME} PRIVATE "-g")
#  target_compile_options(${PROJECT_NAME} PRIVATE "-fdefault-integer-8")
endif()