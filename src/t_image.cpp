//
//  image.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/18/15.
//

#include "t_image.h"
#include <stdlib.h>
#include <string.h>
#include "helper.h"
#include <math.h>
#include "define.h"
#include <GEL/GLGraphics/SOIL.h>

#include "CImg.h"




int scale_factor = 3;
t_image scaled_image;

inline Vec2 get_coord_barry(Vec2_array & tris, double xi1, double xi2) {
    return tris[0]*xi1 + tris[1]*xi2 + tris[2]*(1- xi1 -xi2);
}

void t_image::load_image(std::string const file_path) // Using SOIL from GEL
{
    try {
        unsigned char* _image = SOIL_load_image(file_path.c_str(),
                                                &m_width, &m_height, nullptr, SOIL_LOAD_L);
        

        
        // resolve jagged pixel
        {
            cimg_library::CImg<unsigned char> img(_image, m_width, m_height, 1,1);
//            img.resize_doubleXY();
            img.resize_tripleXY();
            int n = std::min(m_width, m_height)*0.02;
            img = img.blur_median(n);
            
            // init scale image
            
            scaled_image.m_width = scale_factor * m_width;
            scaled_image.m_height = scale_factor* m_height;
            
            scaled_image.m_data = std::vector<float>(scaled_image.m_width * scaled_image.m_height, 0);
            for (int i = 0; i < scaled_image.m_width; i++) {
                for (int j = 0; j < scaled_image.m_height; j++) {
                    scaled_image(i,j) = (float)img(i,j)/255.;
                }
            }

            scaled_image.set_gl_texture();
        }
        
        
        SOIL_free_image_data(_image);
        
        m_data = std::vector<float>(m_width * m_height, 0);
        for (int i = 0; i < m_data.size(); i++) {
            m_data[i] = (float)_image[i]/255.;
        }
        
        set_gl_texture();
    } catch (std::exception e) {
        std::cerr << "Error " << e.what();
    }
    

}

template <typename T>
T t_image::get_sum_on_tri(Vec2_array tris, std::function<T(Vec2)> get_v){
    double area = helper_t::area(tris);
    int N = std::max(1.0,std::ceil(std::sqrt(2*area)));
    double dxi = 1./(double)N;
    
    double xi1, xi2;
    Vec2 p;
    T sum = T(0.0);
    for (int i1 = 0; i1 < N; i1 ++) {
        
        xi1 = i1*dxi;
        
        for (int i2 = 0; i2 < N - i1; i2++) {
            
            xi2 = i2 * dxi;
            
            if (i1 + i2 == N - 1) {
                p = get_coord_barry(tris, (xi1 + dxi/3.), xi2 + dxi/3.);
                sum += get_v(p)*0.5;
            }else{
                p = get_coord_barry(tris, (xi1 + dxi/2.), xi2 + dxi/2.);
                sum += get_v(p);
            }
        }
    }
    
    sum *= 2*area/(double)(N*N); // |j| = 2A
    
    return sum;
}

float & t_image::operator()(int x, int y)
{

    if(x<0)x=0;
    if(x >= m_width) x= m_width -1;
    

    if(y<0)y=0;
    if(y >= m_height) y= m_height -1;
    
    return m_data[x + y * m_width];
}

void t_image::draw_upsample_image()
{
    scaled_image.draw_openGL(1.0/scale_factor);
}

void t_image::draw_image()
{
    draw_openGL();
}

void t_image::normalize()
{
    float maxv = -INFINITY;
    for (auto & v : m_data) {
        maxv = std::max(v, maxv);
    }
    
    for (auto & v : m_data) {
        v = v/maxv;
    }
}

void t_image::set_gl_texture() {
    
    BYTE* texture_buf = (BYTE*)malloc( m_width* m_height * 3 * sizeof(BYTE));
    BYTE* ptr = texture_buf;
    
    for (int j = 0; j < m_height; ++j) {
        for (int i = 0; i < m_width; ++i) {
            
            BYTE color = (BYTE) (get_intensity_bilinear(i,j) * 255);
            
            *(ptr++) = color;
            *(ptr++) = color;
            *(ptr++) = color;
        }
    }
    
    glPixelStorei ( GL_UNPACK_ALIGNMENT,   1 );
    
    glGenTextures(1, &tex_ID);
    
    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, m_width, m_height, 0, GL_RGB, GL_UNSIGNED_BYTE, texture_buf);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

// Using one point
float t_image::get_intensity_linear(float x, float y)
{
    int x_i = (int)x;
    int y_i = (int)y;
    
    return (*this)(x_i, y_i);
}


// Linear interpolation
// Using for point
float t_image::get_intensity_bilinear(float x, float y)
{
//    CROP_W(x)
//    CROP_H(y)

    int x_i = (int)x;
    int y_i = (int)y;
    double ep_x = x - x_i;
    double ep_y = y - y_i;

    double vdown = (*this)(x_i, y_i)*(1-ep_x) + (*this)(x_i + 1, y_i)*ep_x;
    double vup = (*this)(x_i, y_i+1)*(1-ep_x) + (*this)(x_i + 1, y_i+1)*ep_x;
    double v = vdown * (1 - ep_y) + vup * ep_y;

    return v;
}

// Using 16 point
float get_intensity_bicubic(float x, float y)
{
    int x_i = (int)x;
    int y_i = (int)y;
    double ep_x = x - x_i;
    double ep_y = y - y_i;
    
    return INFINITY;
}

float t_image::get_intensity(int x, int y)
{
    return (*this)(x,y);
}

double t_image::get_sum_on_tri_intensity(Vec2_array tris)
{
    return get_sum_on_tri<double>(tris,
      std::function<double(Vec2)>([this](Vec2 p)
      {
          return this->get_intensity_bilinear(p[0], p[1]);
      }
      ));
}

double t_image::get_sum_on_tri_differ(Vec2_array tris, double ci)
{
    auto ci_temp = ci;
    return get_sum_on_tri<double>(tris,
std::function<double(Vec2)>([this, ci_temp](Vec2 p)
      {
          return std::pow((this->get_intensity_bilinear(p[0], p[1]) - ci_temp), 2);
      }
      ));
}

double t_image::get_tri_intensity_f(Vec2_array tris, double * area)
{
    if(area)
        *area = helper_t::area(tris);
    
    return get_sum_on_tri_intensity(tris);
}

void t_image::draw_openGL(double scale)
{
    double w = m_width;
    double h = m_height;
    
    glPushMatrix();
    glScaled(scale, scale, scale);
    
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);
    
    glTexCoord2f(0.0, 0.0);
    glVertex2f(0.0, 0.0);
    
    glTexCoord2f(1.0, 0.0);
    glVertex2f(w, 0.0);
    
    glTexCoord2f(1.0, 1.0);
    glVertex2f(w, h);
    
    glTexCoord2f(0.0, 1.0);
    glVertex2f(0.0, h);
    
    glEnd();
    
    glDisable(GL_TEXTURE_2D);
    
    //  Draw border
    glColor3f(1, 0, 0);
    glBegin(GL_LINE_LOOP);
    glVertex2f(0.0, 0.0);
    glVertex2f(w, 0.0);
    glVertex2f(w, h);
    glVertex2f(0.0, h);
    glEnd();
    
    glPopMatrix();
}

float t_image::get_intensity_bilinear_upscale(float x, float y)
{
    return scaled_image.get_intensity_bilinear(x*scale_factor, y*scale_factor);
}
