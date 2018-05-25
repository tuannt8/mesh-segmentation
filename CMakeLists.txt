cmake_minimum_required(VERSION 3.6)
project(dsc_seg)

set(CMAKE_CXX_STANDARD 14)

set(EXEC dsc_seg)

# DSC lib
add_subdirectory(2D-DSC)
include_directories(2D-DSC/src)

# GEL
add_subdirectory(GEL)
include_directories(GEL)

if (APPLE)
FIND_LIBRARY(COCOA_LIBRARY Cocoa)
endif (APPLE)

# Open GL
find_package (OpenGL)
find_package (GLEW REQUIRED)
find_package(GLUT REQUIRED)
find_package(glfw3)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} $(GLEW_INCLUDE_DIRS))

# Lib png
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

# Other lib: CImg
include_directories(other)

# Build executable file
include_directories(segment)
aux_source_directory(segment SRC)
add_executable(${EXEC} ${SRC})

# Link library
target_link_libraries(${EXEC} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} 2D-DSC GEL OpenGL::GLU OpenGL::GL ${GLFW_LIBRARIES} ${PNG_LIBRARY}  GLEW::GLEW ${COCOA_LIBRARY})
