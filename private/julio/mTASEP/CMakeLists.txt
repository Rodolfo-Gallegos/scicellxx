# Indicate source files and dependencies in the files
#SET(SRC_demo_mTASEP cc_mTASEP.cpp demo_mTASEP.cpp)
SET(SRC_demo_mTASEP demo_mTASEP.cpp)

# Do not create executable but add them to the 'make' list
ADD_EXECUTABLE(julio_demo_mTASEP EXCLUDE_FROM_ALL ${SRC_demo_mTASEP})

# Indicate linking libraries
SET(LIB_demo_mTASEP ${LIB_scicellxx_lib})

# ... and link againts them  
#TARGET_LINK_LIBRARIES(julio_demo_mTASEP ${LIB_demo_mTASEP})
# Temporal solution for MPI
# when linking a target use the following (supported in CMake 3.9+)
#target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)

TARGET_LINK_LIBRARIES(julio_demo_mTASEP PUBLIC ${LIB_demo_mTASEP} MPI::MPI_CXX)

# Check if the output/bin directory exists
IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  # Then create the directory
  FILE(MAKE_DIRECTORY "${bin}")
ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Set directory where to create the executables
set_target_properties( julio_demo_mTASEP
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )

