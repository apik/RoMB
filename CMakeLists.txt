CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

## Here comes the name of your project:

SET(PROJECT_NAME "RoMB")

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

## Change 'hg-tip' to the current version number of your project if you wish.
## Optionally, you can leave it as is and set PROJECT_VERSION from
## the cmake-gui when you make a release.
## The last parameter is a help string displayed by CMAKE.

SET(PROJECT_VERSION "hg-tip"
    CACHE STRING "${PROJECT_NAME} version string")

## Do not edit this.
PROJECT(${PROJECT_NAME})

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

## The next part looks for LEMON. Typically, you don't want to modify it.
##
## First, it tries to use LEMON as a CMAKE subproject by looking for
## it in the 'lemon' or 'deps/lemon' subdirectories or in directory
## given by the LEMON_SOURCE_ROOT_DIR variable.
## If LEMON isn't there, then CMAKE will try to find an installed
## version of LEMON. If it is installed at some non-standard place,
## then you must tell its location in the LEMON_ROOT_DIR CMAKE config
## variable. (Do not hard code it into your config! Others may keep
## LEMON at different places.)

#FIND_PATH(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
#  PATHS ${CMAKE_SOURCE_DIR}/lemon ${CMAKE_SOURCE_DIR}/deps/lemon
#  NO_DEFAULT_PATH
#  DOC "Location of LEMON source as a CMAKE subproject")

#IF(EXISTS ${LEMON_SOURCE_ROOT_DIR})
#  ADD_SUBDIRECTORY(${LEMON_SOURCE_ROOT_DIR} deps/lemon)
#  SET(LEMON_INCLUDE_DIRS
#    ${LEMON_SOURCE_ROOT_DIR}
#    ${CMAKE_BINARY_DIR}/deps/lemon
#  )
#  SET(LEMON_LIBRARIES lemon)
#  UNSET(LEMON_ROOT_DIR CACHE)
#  UNSET(LEMON_DIR CACHE)
#  UNSET(LEMON_INCLUDE_DIR CACHE)
#  UNSET(LEMON_LIBRARY CACHE)
#ELSE()
#  FIND_PACKAGE(LEMON QUIET NO_MODULE)
#  FIND_PACKAGE(LEMON REQUIRED)
#ENDIF()

## This line finds doxygen (for document creation)

FIND_PACKAGE(Doxygen)

## These are the include directories used by the compiler.

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
#  ${LEMON_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/cuba
)

## Here we define an executable target. Its name is 'lemon-project' and
## is compiled from 'main.cc'. You can add more source files separated
## with whitespaces (including newlines). If you want to build more
## executables, simple repeat (and edit) the following ADD_EXECUTABLE and
## TARGET_LINK_LIBRARIES statements.

LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${PROJECT_SOURCE_DIR}/cuba/ )
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( demos )


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src )

#ADD_EXECUTABLE(romb_test main.cpp )
#TARGET_LINK_LIBRARIES(romb_test romb cuba ginac glpk boost_filesystem ${LEMON_LIBRARIES})

