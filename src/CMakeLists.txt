message(STATUS "(EC) Generating random_tools Makefile")

file(GLOB PROJECT_INCLUDE_FILES randomfunctions.h randomgeneric.h)
set(PROJECT_C_FILES randomgeneric.c random_gaussian.c random_mt19937.c random_r250.c random_shr3.c random_xsr128.c random_xsr128r.c)

add_library(random STATIC
   ${PROJECT_INCLUDE_FILES}
   ${PROJECT_C_FILES}
)
set_target_properties(random PROPERTIES
   VERSION ${PROJECT_VERSION}
   PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
   INTERFACE_COMPILE_DEFINITIONS HAVE_RANDOM=${PROJECT_VERSION}
   POSITION_INDEPENDENT_CODE ON)
  
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.h
   PRE_BUILD
   COMMAND echo '\#include \"randomgeneric.h\"' >  ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.h
   COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && cat ${PROJECT_C_FILES} | grep -w InTc | sed 's:[\t ]*//[\t ]*!InTc!.*:;:' >> ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.h
   COMMENT "Generating randomfunctions.h"
   DEPENDS ${PROJECT_C_FILES}
   USES_TERMINAL
)
# Add a custom target to regenerate randomfunctions.h
add_custom_target(randomfunctionsH DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.h)
add_dependencies(random randomfunctionsH)
     
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.inc
   PRE_BUILD
   COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && cat ${PROJECT_C_FILES} | grep -w InTf | sed 's:[\t ]*!InTf!.*::' >> ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.inc
   COMMENT "Generating randomfunctions.inc"
   DEPENDS ${PROJECT_C_FILES}
   USES_TERMINAL
)
# Add a custom target to regenerate randomfunctions.inc
add_custom_target(randomfunctionsINC DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.inc)
add_dependencies(random randomfunctionsINC)

# Installation
set(DISTINCLUDES
   ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.h
   ${CMAKE_CURRENT_BINARY_DIR}/randomfunctions.inc
)

install(TARGETS random
   EXPORT ${PROJECT_NAME}-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   PUBLIC_HEADER DESTINATION include
   INCLUDES DESTINATION include
)

install(FILES ${DISTINCLUDES} DESTINATION include)
