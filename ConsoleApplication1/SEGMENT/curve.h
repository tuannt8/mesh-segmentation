//
//  curve.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC__curve__
#define __DSC__curve__

#include <stdio.h>
#include "define.h"
#include "DSC.h"
#include "texture_helper.h"
#include "image.h"

// Parametric curve
class curve:public std::vector<DSC2D::DeformableSimplicialComplex::node_key>{
public:
    typedef DSC2D::DeformableSimplicialComplex::node_key node_key;
    

private:
    std::vector<Vec2> second_derivative_;
    std::vector<Vec2> forth_derivative_;
    
    double m_in_, m_out_; // Mean intensity iinside and outside
public:
    curve();
    ~curve();
    
    // Draw the curve
    void split_edge(dsc_obj &dsc, image &img);
    void draw(DSC2D::DeformableSimplicialComplex &dsc,
              DSC2D::vec3 color = DSC2D::vec3(1., 0., 0.));
    
    std::vector<Vec2> get_node_force(dsc_obj &dsc, image &img);
    std::vector<Edge_key> get_edge_list(dsc_obj &dsc);
    
    // Compute second and forth derivative
    void update_derivative(DSC2D::DeformableSimplicialComplex &dsc);
    void update_mean_intensity(dsc_obj &complex, image &img);
    
    // Get derivative
    Vec2 derive2(int idx){
        return second_derivative_[idx];
    };
    Vec2 derive4(int idx){
        return forth_derivative_[idx];
    };
    double m_in(){return m_in_;}
    double m_out(){return m_out_;}
    
private:
    node_key cycle_at(int idx);
};

#endif /* defined(__DSC__curve__) */
