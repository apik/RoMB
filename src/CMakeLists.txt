
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
#ADD_SUBDIRECTORY( ara )

#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ara )

ADD_LIBRARY(romb uf.cpp utils.cpp romb.cpp constracc.cpp mbintegral.cpp collect_square.cpp romb_excompiler.cpp code_generation.cc shift.cpp ppl_interface.cpp )
#TARGET_LINK_LIBRARIES(romb cuba ginac cln ppl)

