////
////  texture_helper.cpp
////  DSC
////
////  Created by Tuan Nguyen Trung on 2/11/15.
////  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
////

//#include "texture_helper.h"
//#include <stdlib.h>
//#include <string.h>
//#include "helper.h"
//#include <math.h>


////GLuint texture_helper::LoadTexture( const char * filename ){
////    std::string filePath = std::string(DATA_PATH) + std::string(filename);
//    return -1;
////    std::string filePath(filename);
    
////    std::cout << filePath << std::endl;

////    GLuint tex = SOIL_load_OGL_texture(filePath.c_str(),
////                                       SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS);

    
////    int width, height, chanel;
////    unsigned char* _image = SOIL_load_image(filePath.c_str(),
////                                        &width, &height, &chanel, SOIL_LOAD_RGB);
////    size_ = Vec2(width, height);
    
////    // Gray image
////    gray_image_ = new double[width*height];
//////    for (int i = 0; i < height; i++) {
//////        for (int j = 0; j < width; j++) {
//////            int idx = i*width + j;
//////            gray_image_[idx] = 0.0;
//////            gray_image_[idx] += _image[idx*chanel];
//////            gray_image_[idx] += _image[idx*chanel + 1];
//////            gray_image_[idx] += _image[idx*chanel + 2];
//////        }
//////    }
    

////    texture.push_back(tex);
////    tex_sizes.push_back(Vec2(width, height));
////    SOIL_free_image_data(_image);
    
////    image_.load(filePath.c_str());
////   // image_.mirror('y');
    
////    image_ = (image_.get_channel(0) + image_.get_channel(1) + image_.get_channel(2)) / 3;
    
////    return tex;
//}

//texture_helper::texture_helper(){
//    LoadTexture(IMAGE_NAME.c_str());
//}

//texture_helper::~texture_helper(){
//    // glDeleteTexture
//    if (gray_image_) {
//        delete gray_image_;
//    }
    
//}

//void texture_helper::map_texture(int tex_id){
//    return;
//    glEnable(GL_TEXTURE_2D);
//    glBindTexture(GL_TEXTURE_2D, texture[tex_id]);
//}

//void texture_helper::end(){
//    glDisable(GL_TEXTURE_2D);
//}

//#pragma mark - Public

//void texture_helper::drawImage(int window_x){
//    double pointSize = (double)window_x / image_.width();
//    glPointSize(pointSize);
//    glBegin(GL_POINTS);
//    for (int i = 0; i < image_.width(); i++) {
//        for (int j = 0; j < image_.height(); j++) {
//            double g =grayd(i, j);
//            glColor3f(g, g, g);
//            glVertex2d((double)i, (double)j);
//        }
//    }
//    glEnd();
//    glPointSize(1.0);
//}

//bool texture_helper::is_tri_intersect_phase(std::vector<Vec2> pts){
    
//    int pixel_cont = 0, phase_count = 0;
//    get_triangle_intensity(pts, pixel_cont, phase_count);
    
//    return  (double)phase_count/pixel_cont > 0.7;
//}

//void texture_helper::get_triangle_intensity(std::vector<Vec2> pts,
//                                            int &pixel_count, int &total_i){
//    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
//    for (auto p: pts){
//        min[0] = std::min(min[0], p[0]);
//        min[1] = std::min(min[1], p[1]);
//        max[0] = std::max(max[0], p[0]);
//        max[1] = std::max(max[1], p[1]);
//    }
    
    
//    //Convert to int
//    total_i = 0.0;
//    pixel_count = 0;
//    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
//        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
//            if (helper_t::is_point_in_tri(Vec2(i,j), pts)) {
//                pixel_count ++;
//                total_i += phase(i, j);
//            }
//        }
//    }
//}

//double texture_helper:: average_intensity(std::vector<Vec2> pts){
//    Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
//    for (auto p: pts){
//        min[0] = std::min(min[0], p[0]);
//        min[1] = std::min(min[1], p[1]);
//        max[0] = std::max(max[0], p[0]);
//        max[1] = std::max(max[1], p[1]);
//    }
    
//    //Convert to int
//    int num = 0;
//    double count = 0;
//    for (int i = floor(min[0]); i < ceil(max[0]); i++) {
//        for (int j = floor(min[1]); j < ceil(max[1]); j++) {
//            if (helper_t::is_point_in_tri(Vec2(i,j), pts)) {
//                num ++;
//                count += phase(i, j);
//            }
//        }
//    }
    
//    return count / (double)num;
//}

//double texture_helper::grayd(int x, int y){
//    return gray(x, y) / 255.0;
//}

//unsigned int texture_helper::gray(int x, int y){
//    int r = image_.height() - y - 1;
//    int c = x;
//    return image_.data()[r * image_.width() + c];
//}

//int texture_helper::phase(double x, double y){
//    return phase((int)std::floor(x), (int)std::floor(y));
//}

//int texture_helper::phase(int x, int y){
//    if (gray(x, y) > 250) {
//        return 0;
//    }else{
//        return 1;
//    }
//}

//Vec2 texture_helper::get_local_norm(dsc_obj &complex, Node_key key, bool outside){
//    Vec2 norm(0.0);
    
//    if(key.get_index() == debug_num[0]){
//        //Debug here
//    }
    
//    int phase = outside? 0 : 1;
    
//    for (auto hew = complex.walker(key); !hew.full_circle(); hew = hew.circulate_vertex_cw()) {
//        auto fid = hew.face();
        
//        if (complex.get_label(fid) == phase) { // Currently take only outward pointer
//            auto vids = complex.get_verts(fid);
            
//            // Get bisector vector
//            Vec2 bisector(0.0);
//            Vec2 e[2];int num = 0;
//            for (auto vid : vids) {
//                if (vid != key) {
//                    Vec2 n = complex.get_pos(vid) - complex.get_pos(key);
//                    n.normalize();
//                    bisector += n;
//                    e[num++] = n;
//                }
//            }
//            bisector.normalize();
//            double cos_angle = sin(std::acos( DSC2D::Util::dot(e[0], e[1]) ));
            
//            auto pts = complex.get_pos(fid);
//            int pixel_cont = 0, phase_count = 0;
//            get_triangle_intensity(pts, pixel_cont, phase_count);
//            double intsity = (outside? (double)phase_count : (double)(pixel_cont - phase_count) )/ pixel_cont;
            
//            norm += bisector * intsity * cos_angle;
//        }
//    }
//    if (norm.length() > 0) {
//        norm.normalize();
        
//        norm = norm * 0.7 + complex.get_normal(key)*(outside? 1:-1);
        
//        norm.normalize();
//        return norm;
//    }
//    else{
//        return complex.get_normal(key)*(outside? 1:-1);
//    }
//}
