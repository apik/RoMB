CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

## Here comes the name of your project:
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

## Change 'hg-tip' to the current version number of your project if you wish.
## Optionally, you can leave it as is and set PROJECT_VERSION from
## the cmake-gui when you make a release.
## The last parameter is a help string displayed by CMAKE.

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


## These are the include directories used by the compiler.

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${LEMON_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/cuba
)

## Here we define an executable target. Its name is 'lemon-project' and
## is compiled from 'main.cc'. You can add more source files separated
## with whitespaces (including newlines). If you want to build more
## executables, simple repeat (and edit) the following ADD_EXECUTABLE and
## TARGET_LINK_LIBRARIES statements.

#LINK_DIRECTORIES( ${LINK_DIRECTORIES} ./cuba/ )
#ADD_SUBDIRECTORY( src )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src )
SET(ROMB_LIBS romb cuba ginac glpk boost_filesystem romb ${LEMON_LIBRARIES})

## One-loop B0 function example
ADD_EXECUTABLE(B0_test B0.cpp )
TARGET_LINK_LIBRARIES(B0_test ${ROMB_LIBS})

## One-loop massless box 
ADD_EXECUTABLE(boxm0_test boxm0.cpp )
TARGET_LINK_LIBRARIES(boxm0_test ${ROMB_LIBS})

## One-loop massive box
ADD_EXECUTABLE(boxM_test boxM.cpp )
TARGET_LINK_LIBRARIES(boxM_test ${ROMB_LIBS})

## Two-loop sunset example
ADD_EXECUTABLE(sunset2loop_test sunset2loop.cpp )
TARGET_LINK_LIBRARIES(sunset2loop_test ${ROMB_LIBS})

## One-loop B0 function example
#ADD_EXECUTABLE(B0_test B0.cpp )
#TARGET_LINK_LIBRARIES(B0_test ${ROMB_LIBS})

## One-loop B0 function example
#ADD_EXECUTABLE(B0_test B0.cpp )
#TARGET_LINK_LIBRARIES(B0_test ${ROMB_LIBS})

## One-loop B0 function example
#ADD_EXECUTABLE(B0_test B0.cpp )
#TARGET_LINK_LIBRARIES(B0_test ${ROMB_LIBS})