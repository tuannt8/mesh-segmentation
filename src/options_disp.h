//
//  options_disp.h
//  DSC_seg_integral
//
//  Created by Tuan Nguyen Trung on 7/14/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef DSC_seg_integral_options_disp_h
#define DSC_seg_integral_options_disp_h
#include <map>
#include <stdlib.h>
#include <string>

class options_disp{
public:
    
    static std::map<std::string, bool> map_options;
    
    static int width_view;
    static int font_height;
    static int gap;
    
    static int WIN_SIZE_Y;
public:
    static bool get_option(std::string opt, bool default_value = false);
    
    static void draw(int width, int height);
    
    static void mouse_func(int button, int state, int x, int y);
};
#endif
