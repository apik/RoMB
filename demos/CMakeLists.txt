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

LINK_DIRECTORIES( ${LINK_DIRECTORIES} ${PROJECT_SOURCE_DIR}/cuba/ )
#ADD_SUBDIRECTORY( src )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src )
SET(ROMB_LIBS romb cuba ginac glpk boost_filesystem gsl gslcblas romb ${LEMON_LIBRARIES})

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

## Two-loop self-energy examples
ADD_EXECUTABLE(tsil_test tsil.cpp )
TARGET_LINK_LIBRARIES(tsil_test ${ROMB_LIBS})

## One-loop E0 function example
ADD_EXECUTABLE(pentagon_test pentagon.cpp )
TARGET_LINK_LIBRARIES(pentagon_test ${ROMB_LIBS})

## One-loop F0 function example
ADD_EXECUTABLE(hexagon_test hexagon.cpp )
TARGET_LINK_LIBRARIES(hexagon_test ${ROMB_LIBS})


## Three-loop bubbles example
ADD_EXECUTABLE(bubbles3loop_test bubbles3loop.cpp )
TARGET_LINK_LIBRARIES(bubbles3loop_test ${ROMB_LIBS})

## Two-loop Double-Box function example
ADD_EXECUTABLE(doubleBox_test doubleBox.cpp )
TARGET_LINK_LIBRARIES(doubleBox_test ${ROMB_LIBS})

## 5loop Bubble  example
ADD_EXECUTABLE(5loopBubble_test 5loopBubble.cpp )
TARGET_LINK_LIBRARIES(5loopBubble_test ${ROMB_LIBS})
