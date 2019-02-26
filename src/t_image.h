//
//  image.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/18/15.
//

#ifndef __DSC__image__
#define __DSC__image__

#include <stdio.h>
#include "define.h"
#include <memory>

#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>
#include <Util/ArgExtracter.h>
#else
#include <GEL/Util/ArgExtracter.h>
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#endif

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
 data stored by row/column? order
 */

struct intensity_out{
    double total_differ;
    double area;
    int total_pixel;
};

#define UPSAMPLE_LEVEL 2
class t_image_upsample;

class t_image
{
    friend t_image_upsample;
    std::shared_ptr<t_image_upsample> m_data_upsample;
    
    std::vector<float> m_data; // data scale [0-1]
    int m_width, m_height;
    
    GLuint tex_ID;
    
private:
    void normalize();
    void set_gl_texture();
    
    // Get value on the whole triangle
    // See https://goo.gl/MZTqWq
    template <typename T>
    T get_sum_on_tri(Vec2_array tris, std::function<T(Vec2)> get_v);
    
    void draw_openGL(double scale = 1);
public:
    void load_image(std::string const file_path); // Using SOIL from GEL
    
    float & operator()(int x, int y);
//    float operator()(float x, float y);
    Vec2i size(){return Vec2i(m_width, m_height);}
    
    void draw_image();
    void draw_upsample_image();
    
    float get_intensity_linear(float x, float y);
    float get_intensity_bilinear(float x, float y);
    float get_intensity_bilinear_upscale(float x, float y);
    float get_intensity_bicubic(float x, float y);
    float get_intensity(int x, int y);
    
    double get_sum_on_tri_intensity(Vec2_array tris);
    double get_sum_on_tri_differ(Vec2_array tris, double ci);
    double get_tri_intensity_f(Vec2_array tris, double * area = nullptr);
};



#endif /* defined(__DSC__image__) */
