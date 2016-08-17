//
//  image.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/18/15.
//

#ifndef __DSC__image__
#define __DSC__image__

#include <stdio.h>
//#include <GL/glew.h>
#include "define.h"


#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>
#include <Util/ArgExtracter.h>
#else
#include <GEL/Util/ArgExtracter.h>
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#endif

#include "CImg.h"

typedef unsigned char BYTE;
#define MAX_BYTE 255

/*
 Image 
 Dimension: 
 y ^
   |
  0|------>
          x (width)
 Scale: 0-1
 */

struct intensity_out{
    double total_differ;
    double area;
    int total_pixel;
};

class image: public cimg_library::CImg<BYTE>{
private:
    // Row order storage idx = y * width + x
    Vec2_array gradient_;
    
private:
    void compute_gradient();

public:
    double ** _buffer;

	GLuint tex_ID;
public:
    void load_image(std::string const file_path);
    // Draw in OpenGL coordinate
    void draw_image(int window_width);
    void draw_grad(int window_width);
    
    // double: 0 - 1.0
    double get_intensity(int x, int y);
    
    /**
    // Interpolate intensity to smooth function
    /// Computation base on its 4 neighbor pixel
     */
    double get_intensity_f(double x, double y);
    // get \int (g - ci)^2 d\Omega
    double get_tri_differ_f(Vec2_array tris, double ci);
    // \int g d\Omega
    double get_tri_intensity_f(Vec2_array tris, double * area = nullptr);

    double get_sum_gradient_tri(Vec2_array tris, double * area = nullptr);
    
    void debug_integral(Vec2_array tris);
    
    // Get value on the whole triangle
    template <typename T>
    T get_sum_on_tri(Vec2_array tris, std::function<T(Vec2)> get_v);
    
    double get_sum_on_tri_intensity(Vec2_array tris);
    double get_sum_on_tri_variation(Vec2_array tris);
    double get_sum_on_tri_differ(Vec2_array tris, double ci);
    double ci_temp;
    
    /*
     * Get total variation with filtered
     */
    double get_sum_on_tri_variation(Vec2_array tris, double pixel_gap);
    double smooth_filter(double xi1, double xi2, double xi3, double gap);
    
    /*
     * Edge energy
     * By multiply gradient with normal of the edge
     */
    double get_edge_energy(Vec2 p1, Vec2 p2);
    double get_edge_energy(Vec2 p1, Vec2 p2, double spread);
    
    // total intensity inside a triangle
    void get_tri_intensity(Vec2_array tris, int * total_pixel, double * total_intensity);
    
    // Intensity differ with assumed mean intensity.
    void get_tri_differ(Vec2_array tris, int *total_pixel, double * total_differ, double ci);
    intensity_out get_tri_differ(Vec2_array tris, double ci);
    
    Vec2 size(){return Vec2(width(), height());}
    
    // gradient
    Vec2 grad(int x, int y);
    // Interpolated gradient
    Vec2 grad_f(double x, double y);
    
    void set_gl_texture();
};

#endif /* defined(__DSC__image__) */
