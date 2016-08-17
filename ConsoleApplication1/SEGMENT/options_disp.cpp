//
//  options_disp.cpp
//  DSC_seg_integral
//
//  Created by Tuan Nguyen Trung on 7/14/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include <stdio.h>
#include "options_disp.h"
#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>
#include <Util/ArgExtracter.h>
#else
#include <GEL/Util/ArgExtracter.h>
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#endif

std::map<std::string, bool> options_disp::map_options;
 int options_disp::width_view = 150;
 int options_disp::font_height = 8;
 int options_disp::gap = 2;
 int options_disp::WIN_SIZE_Y;

bool options_disp::get_option(std::string opt, bool default_value)
{
    if (map_options.find(opt) == map_options.end())
    {
        map_options.insert(std::make_pair(opt, default_value));
    }
    
    return map_options[opt];
}

void print_gl(const double &x, const double &y, const char* str){
    glRasterPos2f(x, y);
    for (const char *c = str; *c != '\0'; c++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }
}

void options_disp::draw(int width, int height)
{
    WIN_SIZE_Y = height;
    glLoadIdentity();

    gluOrtho2D(0, width_view, 0, height);
    glViewport(0, 0, width_view, height);
    
    glPushMatrix();
    
    glColor4f(0.2, 0.2, 0.2, 0.2);
    glBegin(GL_QUADS);{
        glVertex2f(0, 0);
        glVertex2f(width_view, 0);
        glVertex2f(width_view, height);
        glVertex2f(0, height);
    }glEnd();
    
    int i = 0;
    for (auto m : map_options)
    {
        int y = i*(font_height + 2*gap) + gap;
        if (m.second) {
            glColor3f(0, 1, 0);
        }else{
            glColor3f(1, 0, 0);
        }
        print_gl(gap, y, m.first.c_str());
        i ++;
    }
    glPopMatrix();
    
    
}

void options_disp::mouse_func(int button, int state, int x, int y){
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP){
        if (x < width_view) {
            int idx = (WIN_SIZE_Y-y)/(font_height + 2*gap);
            if (idx < map_options.size()) {
                auto iter = map_options.begin();
                for (int i = 0; i<idx; i++) {
                    iter++;
                }
                iter->second = !iter->second;
            }
        }
    }
}
