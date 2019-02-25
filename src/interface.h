//
//  interface_dsc.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/9/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#pragma once

#include "velocity_function.h"
#include "DSC.h"
#include "log.h"
#include "texture_helper.h"
#include "setting_io.hpp"
#include <memory>
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

#include "dynamics_mul.h"

class interface_dsc{
#pragma mark - Local variable
public:
    std::unique_ptr<DSC2D::DeformableSimplicialComplex> dsc;
    
    // Windows size
    int     WIN_SIZE_X;
    int     WIN_SIZE_Y;
    double  SCALE;
    Vec2i    imageSize;
    std::vector<bool> bDiplay_;
    bool    RUN = false;
    
    // triangle size
    double DISCRETIZATION;
    
    static interface_dsc *instance;
    
    std::unique_ptr<dynamics_mul> dyn_;
    std::unique_ptr<t_image> image_;
    
    int iter = 0;
private:
    int debug_num_[10];
    
#pragma mark - Glut display
public:
    void display();
    
    void animate();
    
    void reshape(int width, int height);
    
    void visible(int v);

    void keyboard(unsigned char key, int x, int y);
    
    void initGL();
    
#pragma mark - class functions
public:
    interface_dsc(int &argc, char** argv);
    virtual ~interface_dsc(){}
    
    static interface_dsc* get_instance(){
        return instance;
    }
  
private:
    void draw();
    void draw_image();
    void draw_coord();
    void update_title();
    void thres_hold_init();
    
private:
    void draw_test();
    void draw_edge_energy();
    void draw_tri_variant();
    void draw_tri_MS_energy();
    
    void load_dsc();
    void back_up_dsc();
    
    void write_triangle_energy();
    
#pragma mark - data
public:
    void init_dsc();
    void threshold_initialization();
    
    void manual_init_dsc();
    void random_init_dsc(int nb_phase);
    
    void init_sqaure_boundary();
    void init_boundary();
    void init_boundary_brain();
    void dynamics_image_seg();
    void circle_init(Vec2 center, double radius, int label);
    
    void export_dsc();
};
