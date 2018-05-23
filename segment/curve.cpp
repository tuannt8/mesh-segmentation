//
//  curve.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "curve.h"
#ifdef WIN32
#include <GL/glew.h>
#include <GL/glut.h>
#include <GLGraphics/SOIL.h>
#else
#include <GEL/GL/glew.h>
#include <GLUT/glut.h>
#include <GEL/GLGraphics/SOIL.h>
#endif

#include "helper.h"

curve::curve(){
    
}

curve::~curve(){
    
}

#define POS(i) dsc.get_pos(cycle_at(i))

std::vector<Edge_key> curve::get_edge_list(dsc_obj &dsc){
    
    std::vector<Edge_key> edge_list;
    
    for (int i = 0; i < size(); i++) {
        auto v1 = cycle_at(i);
        auto v2 = cycle_at(i+1);
        
        for (auto hew = dsc.walker(v1); !hew.full_circle(); hew = hew.circulate_vertex_cw()) {
            if (hew.vertex() == v2) {
                edge_list.push_back(hew.halfedge());
                break;
            }
        }
    }
    
    assert(edge_list.size() == size());
    
    return edge_list;
}

Vec2 get_edge_norm(dsc_obj &dsc, Edge_key ek){
    auto hew = dsc.walker(ek);
    if (dsc.get_label(hew.face()) == 0) {
        hew = hew.opp();
    }
    
    Vec2 outer_norm = dsc.get_pos(hew.vertex()) - dsc.get_pos(hew.opp().vertex());
    outer_norm = Vec2(outer_norm[1], -outer_norm[0]);
    assert(outer_norm.length() > 0);
    outer_norm.normalize();
    return outer_norm;
}

std::vector<Vec2> curve::get_node_force(dsc_obj &dsc, image &img){
    std::vector<Vec2> node_force(size(), Vec2(0.));
    
    std::vector<Edge_key> edge_list = get_edge_list(dsc);
    
    for (int i = 0; i < size(); i++) {
        auto pt1 = dsc.get_pos(cycle_at(i));
        auto pt2 = dsc.get_pos(cycle_at(i+1));
        
        Vec2 e_norm = get_edge_norm(dsc, edge_list[i]);
        
        int l = (int)(pt2-pt1).length();
        double f1 = 0., f2 = 0.;
        for (int j = 0; j <= l; j++) {
            
            Vec2 cur_p = pt1 + (pt2-pt1)*((double)j/(double)l);
            
            double force_mag = (m_out() - m_in()) *
                    (img.get_intensity_i(cur_p[0], cur_p[1])*2 - m_out() - m_in());
            
            f1 += force_mag * (l - j);
            f2 += force_mag * j;
        }
        
        node_force[i] += e_norm * (f1 / (double)l);
        node_force[(i+1)%size()] += e_norm * (f2 / (double)l);
        
    }
    
    return node_force;
}

void curve::draw(DSC2D::DeformableSimplicialComplex &dsc, DSC2D::vec3 color){
    
    glLineWidth(1.5f);
    glColor3dv(color.get());
    glBegin(GL_LINES);
    
    for (int i = 0; i < size(); i++) {
        glVertex2dv(dsc.get_pos(at(i)).get());
        glVertex2dv(dsc.get_pos(at( (i + 1)%size() )).get());
    }
    
    glEnd();
}

void curve::update_derivative(DSC2D::DeformableSimplicialComplex &dsc){
    
    // Second derivative
    second_derivative_.clear();
    second_derivative_.resize(size());
    for (int i = 0; i < size(); i++) {
        second_derivative_[i] = POS(i-1) - POS(i)*2 + POS(i+1);
    }
    
    // Forth derivative
    forth_derivative_.clear();
    forth_derivative_.resize(size());
    for (int i = 0; i < size(); i++) {
        forth_derivative_[i] = POS(i-2) - POS(i-1)*4 + POS(i)*6 - POS(i+1)*4 + POS(i+2);
    }
}

curve::node_key curve::cycle_at(int idx){

    return (*this)[idx % size()];
}

void curve::update_mean_intensity(dsc_obj &complex, image &img){
    m_in_ = 0.0;
    m_out_ = 0.0;
    int count_in = 0, count_out = 0;
    
    for (auto fi = complex.faces_begin(); fi != complex.faces_end(); fi++) {
        auto pts = complex.get_pos(*fi);
        int count;
        int total_i = img.get_triangle_intensity_count(pts, &count);
        
        if (complex.get_label(*fi) == 0) { // Outside
            m_out_ += total_i;
            count_out += count;
        }else{
            m_in_ += total_i;
            count_in += count;
        }
    }
    
    m_in_ /= count_in;
    m_out_ /= count_out;
}

void curve::split_edge(dsc_obj &dsc, image &img) {
    double inten_threshold = 10;

    for (int i = 0; i < size(); i++){
        Node_key n1 = cycle_at(i);
        Node_key n2 = cycle_at(i+1);

        Edge_key ek = HMesh::InvalidHalfEdgeID;
        for(auto hew = dsc.walker(n1); hew.full_circle(); hew.circulate_vertex_cw()){
            if(hew.opp().vertex() == n2){
                if(dsc.get_label(hew.face()) == 0)
                    ek = hew.opp().halfedge();
                else
                    ek = hew.halfedge();
                break;
            }
        }

#ifdef DEBUG
        assert(ek != HMesh::InvalidHalfEdgeID);
#endif

        int count_out;
        auto pt_out = dsc.get_pos(dsc.walker(ek).face());
        auto pt_in = dsc.get_pos(dsc.walker(ek).opp().face());

        int out = img.get_triangle_intensity_count(pt_out, &count_out);

        double inten_out = out / (double)count_out;
        // We solve this case first
        if (inten_out > inten_threshold){
            dsc.split_edge(ek);
        }
    }
}
