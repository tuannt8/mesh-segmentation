# 2D image segmentation with triangle mesh
![demo](/resources/2D-demo.jpg)
# Run

```shell
dsc_seg hamster.txt
```

# Build note
The project contains some submodules, so clone the project with 
```shell
git clone --recurse-submodules https://github.com/tuannt8/mesh-segmentation
```

Dependencies: `libpng`, `OpenGL`, and `Glut`

Comment out the `PyGEL` in `GEL/CMakeLists.txt` as we dont need it

```cmake
include_directories(./src)
aux_source_directory(./src/PyGEL PYG_SRC_LIST)
add_library(PyGEL SHARED ${PYG_SRC_LIST})
target_link_libraries(PyGEL GEL glfw ${OPENGL_glu_LIBRARY} ${GLFW_LIBRARIES})
```