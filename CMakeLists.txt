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
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework OpenGL -framework IOKit")
endif (APPLE)

# Open GL
find_package (OpenGL)
find_package (GLEW REQUIRED)
find_package(GLUT REQUIRED)
find_package(glfw3)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} $(GLEW_INCLUDE_DIRS))



# Other lib: CImg
include_directories(other)

###################################################
# Lib png
# search for pkg-config
include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()
# check for libpng
pkg_check_modules (LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng16 development libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    link_directories (${LIBPNG_LIBRARY_DIRS})
    link_libraries (${LIBPNG_LIBRARIES})
endif ()

# Build executable file
include_directories(segment)
aux_source_directory(segment SRC)
add_executable(${EXEC} ${SRC} ${LIBPNG_LINK_FLAGS})

# Link library
target_link_libraries(${EXEC} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OPENGL_glu_LIBRARY} 2D-DSC GEL ${COCOA_LIBRARY})

