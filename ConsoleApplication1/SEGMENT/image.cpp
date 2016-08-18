//
//  image.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/18/15.
//

#include "image.h"
#include <stdlib.h>
#include <string.h>
#include "helper.h"
#include <math.h>




#define NOISE 10
#define BLUR 1.0

void image::load_image(std::string const file_path){
    try
    {
        load(file_path.c_str());
        this->mirror('y');
        
#ifdef ADD_NOISE
        blur(BLUR); //isotropically
        noise(NOISE);
#endif
        
        set_gl_texture();
//        compute_gradient();
    }
    catch (std::exception e)
    {
        cout << e.what();
    }
}

// Draw in OpenGL coordinate
void image::draw_image(int window_width){

        double w = width();
        double h = height();

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
}

void image::draw_grad(int window_width){
    double pointSize = (double)window_width / this->width();
    glPointSize(pointSize);
    glBegin(GL_LINES);
    for (int i = 0; i < width(); i++) {
        for (int j = 0; j < height(); j++) {
            Vec2 g = grad_f(i, j)*10;
            glColor3f(1.0, 0.0, 0.0);
            
            Vec2 pt((double)i+0.5, (double)j+0.5);
            glVertex2d(pt[0], pt[1]);
            Vec2 top = pt + g;
            glVertex2d(top[0], top[1]);
        }
    }
    glEnd();
    glPointSize(1.0);
}

// 0 - 1.0
double image::get_intensity(int x, int y){
    if (x < 0 or x > width()
        or y < 0 or y > height()) {
        return 0;
    }
    return ((double)(*this)(x, y)) / (double)MAX_BYTE;
}

// 0 - 1.0
double image::get_intensity_f(double x, double y){
    if (x < 0 or x > width()
        or y < 0 or y > height())
    {
        return 0;
    }
    
    int x_i = (int)x;
    int y_i = (int)y;
    double ep_x = x - x_i;
    double ep_y = y - y_i;

    double vdown = get_intensity(x_i, y_i)*(1-ep_x) + get_intensity(x_i + 1, y_i)*ep_x;
    double vup = get_intensity(x_i, y_i+1)*(1-ep_x) + get_intensity(x_i + 1, y_i+1)*ep_x;
    double v = vdown * (1 - ep_y) + vup * ep_y;
    
    return v;
}

Vec2 image::grad_f(double x, double y){
    if (x < 0 or x > width()
        or y < 0 or y > height()) {
        return Vec2(0.0);
    }
    
    int x_i = (int)x;
    int y_i = (int)y;
    double ep_x = x - x_i;
    double ep_y = y - y_i;
    
    auto vdown = grad(x_i, y_i)*(1-ep_x) + grad(x_i + 1, y_i)*ep_x;
    auto vup = grad(x_i, y_i+1)*(1-ep_x) + grad(x_i + 1, y_i+1)*ep_x;
    auto v = vdown * (1 - ep_y) + vup * ep_y;
    
    return v;
}

double image::get_tri_intensity_f(Vec2_array tris, double * area_in){
    double area = helper_t::area(tris);
    int res = (int)std::sqrt(area) + 1;
    double step = 1.0/(double)res;
    
    double inten = 0.0;
    
    for (int i = 0; i < res; i++) {
        double s1 = i*step;
        for (int j = res - i; j > 0; j--) {
            double s2 = j*step;
            
            // Triangle 1
            {
                double ep1 = s1 + step/3.0;
                double ep2 = s2 - step*2.0/3.;
                double ep3 = 1 - ep1 - ep2;
                Vec2 pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                inten += get_intensity_f(pt[0], pt[1]);
            }
            // Triangle 2
            {
                double ep1 = s1 + step*2/3;
                double ep2 = s2 - step*4/3;
                if (ep2 > 0) {
                    double ep3 = 1 - ep1 - ep2;
                    Vec2 pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                    inten += get_intensity_f(pt[0], pt[1]);

                }
            }
            
        }
    }
//    for (double ep1 = step/2.0; ep1 < 1.0; ep1 += step){
//        for(double ep2 = 1 - ep1 - step/2.0; ep2 > 0; ep2 -= step){
//            double ep3 = 1 - ep1 - ep2;
//            Vec2 pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
//            inten += get_intensity_f(pt[0], pt[1]);
//        }
//
//    }
    
    if (area_in) {
        *area_in = area;
    }
    
    return inten * area / (double)(res*res);
}

void image::get_tri_intensity(Vec2_array tris, int * total_pixel, double * total_intensity){
    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for (auto p: tris){
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
    }
    
    int t_pixel = 0;
    double total_inten = 0.0;
    
    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
            if (helper_t::is_point_in_tri(Vec2(i,j), tris)) {
                t_pixel ++;
                total_inten += get_intensity(i, j);
            }
        }
    }
    
    *total_intensity = total_inten;
    *total_pixel = t_pixel;
}

inline Vec2 get_coord_barry(Vec2_array & tris, double xi1, double xi2) {
    return tris[0]*xi1 + tris[1]*xi2 + tris[2]*(1- xi1 -xi2);
}

double image::get_sum_on_tri_differ(Vec2_array tris, double ci){
    ci_temp = ci;
    return get_sum_on_tri<double>(tris,
              std::function<double(Vec2)>([this](Vec2 p)
                  {
                      return std::pow((this->get_intensity_f(p[0], p[1]) - ci_temp), 2);
                  }
                  ));
}

double image::get_sum_on_tri_intensity(Vec2_array tris){
    return get_sum_on_tri<double>(tris,
                std::function<double(Vec2)>([this](Vec2 p)
                    {
                        return this->get_intensity_f(p[0], p[1]);
                    }
                ));
}


double image::get_sum_on_tri_variation(Vec2_array tris){
    return get_sum_on_tri<double>(tris,
          std::function<double(Vec2)>([this](Vec2 p)
                                      {
                                          return this->grad_f(p[0], p[1]).length();
                                      }
                                      ));
}

template <typename T>
T image::get_sum_on_tri(Vec2_array tris, std::function<T(Vec2)> get_v){
    double area = helper_t::area(tris);
    int N = std::ceil(std::sqrt(2*area));
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

intensity_out image::get_tri_differ(Vec2_array tris, double ci){
    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for (auto p: tris){
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
    }
    
    int t_pixel = 0;
    double total_diff = 0.0;
    
    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
            if (helper_t::is_point_in_tri(Vec2(i,j), tris)) {
                t_pixel ++;
                total_diff += (ci - get_intensity(i, j)) * (ci - get_intensity(i, j));
            }
        }
    }
    
    intensity_out out;
    out.total_pixel = t_pixel;
    out.total_differ = total_diff;
    
    return out;
}

double image::smooth_filter(double xi1, double xi2, double xi3, double gap){
    double xi123 = xi1 * xi2 * xi3;
    if (xi123 < gap*(1-gap)*(1-gap)/4) {
        return xi123/( gap*(1-gap)*(1-gap)/4.0 );
    }else{
        return 1;
    }
}

double image::get_edge_energy(Vec2 p1, Vec2 p2){
    double ene = 0;
    Vec2 p12 = p2 - p1;
    Vec2 norm(p12[1], -p12[0]);
    norm.normalize();
    
    double L = p12.length();
    int N = (int)L;
    double dl = L/(double)N;
    for (int i = 0; i < N; i++) {
        auto p = p1 + p12*((i+0.5)*dl);
        ene += std::abs(CGLA::cross(norm, grad_f(p[0], p[1]))) * dl/L;
    }
    
    return ene;
}

double image::get_edge_energy(Vec2 p1, Vec2 p2, double spread){
    double ene = 0;
    Vec2 p12 = p2 - p1;
    Vec2 norm(p12[1], -p12[0]);
    norm.normalize();
    
    double L = p12.length();
    int N = (int)L;
    double dl = L/(double)N;
    for (int i = 0; i < N; i++) {
        auto p = p1 + p12*((i+0.5)*dl);
        auto p1 = p + norm*spread;
        auto p2 = p - norm*spread;
        
        ene += std::abs(CGLA::dot(norm, grad_f(p[0], p[1])));
        ene += std::abs(CGLA::dot(norm, grad_f(p1[0], p1[1])));
        ene += std::abs(CGLA::dot(norm, grad_f(p2[0], p2[1])));
    }
    
    return ene * dl/L;
}

double image::get_sum_on_tri_variation(Vec2_array tris, double pixel_gap){
    
//    Vec2 center = (tris[0] + tris[1] + tris[2])/3.0;
//    for (int i = 0; i < 3; i++) {
//        double  l = (tris[i] - center).length();
//        tris[i] = center + (tris[i] - center)*(l-pixel_gap)/l;
//    }

    Vec2_array new_tris;
    for (int i = 0; i < 3; i++) {
        auto p0 = tris[i];
        auto p1 = tris[(i+1)%3];
        auto p2 = tris[(i+2)%3];
        
        auto p01 = p1 - p0; p01.normalize();
        auto p02 = p2 - p0; p02.normalize();
        auto pd = p01 + p02;
        double sinA = std::abs(CGLA::cross(p01, p02));
        
        Vec2 pn = p0 + pd*(pixel_gap/sinA);
        new_tris.push_back(pn);
    }
    
    return get_sum_on_tri_variation(new_tris);
    
//    double area = helper_t::area(tris);
//    int N = std::ceil(std::sqrt(2*area));
//    double dxi = 1./(double)N;
//    double gap = pixel_gap/std::sqrt(2*area);
//
//    double xi1, xi2;
//    Vec2 p;
//    double sum = 0.0;
//    for (int i1 = 0; i1 < N; i1 ++) {
//
//        xi1 = i1*dxi;
//
//        for (int i2 = 0; i2 < N - i1; i2++) {
//
//            xi2 = i2 * dxi;
//
//            if (i1 + i2 == N - 1) {
//                p = get_coord_barry(tris, (xi1 + dxi/3.), xi2 + dxi/3.);
//                double filterd = smooth_filter(xi1+ dxi/3., xi2 + dxi/3.,
//                                               1 - xi1 - xi2 - 2*dxi/3., gap);
//                sum += grad_f(p[0], p[1]).length()*0.5 * filterd;
//            }else{
//                p = get_coord_barry(tris, (xi1 + dxi/2.), xi2 + dxi/2.);
//                double filterd = smooth_filter(xi1+ dxi/2., xi2 + dxi/2.,
//                                               1 - xi1 - xi2 - dxi, gap);
//                sum += grad_f(p[0], p[1]).length() * filterd;
//            }
//        }
//    }
//
//    sum *= 2*area/(double)(N*N); // |j| = 2A
//
//    return sum;
}

double image::get_sum_gradient_tri(Vec2_array tris, double * area){
    double gE = 0;
    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for (auto p: tris){
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
    }
    
    int t_pixel = 0;
    
    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
            if (helper_t::is_point_in_tri(Vec2(i,j), tris)) {
                t_pixel ++;
                gE += grad(i, j).length();
            }
        }
    }
    
    return gE/t_pixel;
}

void image::get_tri_differ(Vec2_array tris, int *total_pixel, double * total_differ, double ci){
    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for (auto p: tris){
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
    }
    
    int t_pixel = 0;
    double total_diff = 0.0;
    
    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
            if (helper_t::is_point_in_tri(Vec2(i,j), tris)) {
                t_pixel ++;
                total_diff += (ci - get_intensity(i, j)) * (ci - get_intensity(i, j));
            }
        }
    }
    
    *total_differ = total_diff;
    if(total_pixel)
        *total_pixel = t_pixel;
}


void image::compute_gradient(){
    gradient_.resize(width()*height());
    
    cimg_library::CImgList<int> grad_img = get_gradient("xy", 0);
    CImg<int> *gX = grad_img.data(0);
    CImg<int> *gY = grad_img.data(1);
    
    for (int x = 0; x < width(); x++) {
        for (int y = 0; y< height(); y++) {
            double x_g = (*gX)(x, y)/(double)MAX_BYTE;
            double y_g = (*gY)(x, y)/(double)MAX_BYTE;
            gradient_[y * width() + x] = Vec2(x_g,y_g);
        }
    }
  
//    FILE *f = fopen("../grad.txt", "w");
//    for (int x = 0; x < width(); x++) {
//        for (int y = 0; y < height(); y++) {
//            fprintf(f, "%f ", grad(x, y)[0]);
//        }
//        fprintf(f, "\n");
//    }
//    fclose(f);
//    gX->display("gX");
//    gY->display("gY");
//    grad_img.display();
}


Vec2 image::grad(int x, int y){
    if(x < 0 or x >= width()
       or y < 0 or y >= height())
    {
        return Vec2(0.0);
    }
    
    return gradient_[y * width() + x];
}

void image::set_gl_texture() {
    
    BYTE* texture_buf = (BYTE*)malloc( width()* height() * 3 * sizeof(BYTE));
    BYTE* ptr = texture_buf;
    
    for (int j = 0; j < height(); ++j) {
    for (int i = 0; i < width(); ++i) {
        
            BYTE color = (BYTE) (get_intensity_f(i,j) * 255);
        
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
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, width(), height(), 0, GL_RGB, GL_UNSIGNED_BYTE, texture_buf);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

void image::debug_integral(Vec2_array tris){

}

double image::get_tri_differ_f(Vec2_array tris, double ci) {
    double area = helper_t::area(tris);
    int res = (int)std::sqrt(area) + 1;
    double step = 1.0/(double)res;

    double energy = 0.0;
    
    for (int i = 0; i < res; i++) {
        double s1 = i*step;
        for (int j = res - i; j > 0; j--) {
            double s2 = j*step;
            
            // Triangle 1
            {
                double ep1 = s1 + step/3.0;
                double ep2 = s2 - step*2.0/3.;
                double ep3 = 1 - ep1 - ep2;
                Vec2 pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                double g = get_intensity_f(pt[0], pt[1]);
                energy += (g - ci)*(g - ci);
            }
            // Triangle 2
            {
                double ep1 = s1 + step*2/3;
                double ep2 = s2 - step*4/3;
                if (ep2 > 0) {
                    double ep3 = 1 - ep1 - ep2;
                    Vec2 pt = tris[0]*ep1 + tris[1]*ep2 + tris[2]*ep3;
                    double g = get_intensity_f(pt[0], pt[1]);
                    energy += (g - ci)*(g - ci);
                }
            }
            
        }
    }

    return energy * area / (double)(res*res);
}
