cmake_minimum_required(VERSION 3.12)

PROJECT("SBL-Thickness")

find_package(Armadillo REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx)

foreach(SRC_FILE ${SRC_FILES})
	get_filename_component(EXE_FILE ${SRC_FILE} NAME_WE)
	add_executable(${EXE_FILE} ${SRC_FILE})
	target_link_libraries(${EXE_FILE} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${ARMADILLO_LIBRARIES})
endforeach(SRC_FILE ${SRC_FILES})
