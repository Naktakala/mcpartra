cmake_minimum_required(VERSION 3.2)

SET(TARGET mcpartra)
project(${TARGET} C CXX)

include(~/Desktop/ChiTech/chi-tech/ChiResources/Macros/Downstream.cmake)
include(~/Desktop/ChiTech/chi-tech/ChiResources/Macros/Filter.cmake)

set(SOURCES "mcpartra.cc")
file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

FOREACH(item ${MORE_SOURCES})
    IF(${item} MATCHES "lua/mcpartra_maindoc_page.cc")
        LIST(REMOVE_ITEM MORE_SOURCES ${item})
    ENDIF(${item} MATCHES "lua/mcpartra_maindoc_page.cc")
ENDFOREACH(item)

set(SOURCES ${SOURCES} ${MORE_SOURCES})

set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O3")

add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})