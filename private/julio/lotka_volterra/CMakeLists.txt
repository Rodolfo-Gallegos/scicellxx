# Indicate source files and dependencies in the files
SET(SRC_lotka_volterra demo_lotka_volterra.cpp cc_lotka_volterra_odes.cpp)

# Do not create executable but add them to the 'make' list
ADD_EXECUTABLE(julio_demo_lotka_volterra EXCLUDE_FROM_ALL ${SRC_lotka_volterra})

# Indicate linking libraries
SET(LIB_lotka_volterra ${LIB_scicellxx_lib})

# ... and link againts them  
TARGET_LINK_LIBRARIES(julio_demo_lotka_volterra ${LIB_lotka_volterra})

# Check if the output/bin directory exists
IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  # Then create the directory
  FILE(MAKE_DIRECTORY "${bin}")
ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Set directory where to create the executables
set_target_properties( julio_demo_lotka_volterra
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )

