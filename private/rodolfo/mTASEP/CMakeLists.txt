# Indicate source files and dependencies in the files
#SET(SRC_demo_mTASEP_LK cc_mTASEP.cpp demo_mTASEP_LK.cpp)

# Encuentra el kit de herramientas de CUDA
find_package(CUDAToolkit REQUIRED)

if(CUDAToolkit_FOUND)
    message(STATUS "CUDAToolkit encontrado: ${CUDAToolkit_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "CUDAToolkit no encontrado. Asegúrate de tener CUDAToolkit instalado.")
endif()

SET(SRC_demo_mTASEP_LK 
    demo_mTASEP_LK_CUDA.cu)
    # cuda_functions.cu

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G -lcudadevrt -lineinfo")


# Do not create executable but add them to the 'make' list
ADD_EXECUTABLE(rodolfo_demo_mTASEP_LK EXCLUDE_FROM_ALL ${SRC_demo_mTASEP_LK})

# Indicate linking libraries
SET(LIB_demo_mTASEP_LK ${LIB_scicellxx_lib})

# ... and link againts them  
#TARGET_LINK_LIBRARIES(julio_demo_mTASEP_LK ${LIB_demo_mTASEP_LK})
# Temporal solution for MPI
# when linking a target use the following (supported in CMake 3.9+)
#target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)

# Agregar las arquitecturas de GPU para las que se compilará el código CUDA
set_target_properties(rodolfo_demo_mTASEP_LK PROPERTIES CUDA_ARCHITECTURES 86)

# Agregar directorios de inclusión para archivos de encabezado CUDA
include_directories(${CUDAToolkit_INCLUDE_DIRS})

# Enlazar bibliotecas CUDA
target_link_libraries(rodolfo_demo_mTASEP_LK ${CUDAToolkit_LIBRARIES})
TARGET_LINK_LIBRARIES(rodolfo_demo_mTASEP_LK PUBLIC ${LIB_demo_mTASEP_LK})
target_link_libraries(rodolfo_demo_mTASEP_LK PRIVATE CUDA::cudart CUDA::curand)


# Check if the output/bin directory exists
IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  # Then create the directory
  FILE(MAKE_DIRECTORY "${bin}")
ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bin)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-profile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-profile-arcs -ftest-coverage")

# Set directory where to create the executables
set_target_properties( rodolfo_demo_mTASEP_LK
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )

