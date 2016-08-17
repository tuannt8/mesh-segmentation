//
//  gl_debug_helper.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/17/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

/*
 This class help draw a rectangle on the screen
 And return the rectangle for debug information
 */
#ifndef __DSC__gl_debug_helper__
#define __DSC__gl_debug_helper__

#include <stdio.h>
#include "define.h"
#include "image.h"

#define LOG_VERTEX  true
#define LOG_EDGE    true
#define LOG_FACE    true

class gl_debug_helper{
private: // Singleton
    gl_debug_helper(){};
    gl_debug_helper(gl_debug_helper const &) = delete;
    void operator = (gl_debug_helper const &) = delete;
    
private:
    bool active_ = false;
    bool drawing_ = false;
    Vec2 left_down_, right_up_; // The rectangle
    
    Vec2 cur_mouse_pos_;
    
    Vec2 window_draw_left_;
    Vec2 scale_gl_over_win_;
    double WIN_HEIGHT_;
    
    Vec2 to_gl_coord(Vec2 win_coord);
    
    dsc_obj * s_dsc_;
public:
    static void set_dsc(dsc_obj *dsc){get_instance().s_dsc_ = dsc;}
    static bool is_debugging(){return get_instance().active_;};
    static void change_state();
    static void begin();
    static void end();
    
    static gl_debug_helper & get_instance();
    
    static void mouseDown(int button, int state, int x, int y);
    static void mouseMove(int x, int y);
    
    static void coord_transform(Vec2 window_draw_left, Vec2 scale, double win_height);
    // Index of element inside the box
    static void print_debug_info(dsc_obj &complex);
    
    // Index of nearest element
    static void print_debug_info_nearest(dsc_obj &complex);

    static void print_image_info(image &img);
    
    static void draw();
    
private:
    static void update_dsc();
};

#endif /* defined(__DSC__gl_debug_helper__) */
