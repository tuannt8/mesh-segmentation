//
//  dynamics.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC__dynamics__
#define __DSC__dynamics__

#include <stdio.h>
#include "define.h"
#include "DSC.h"
#include "texture_helper.h"
#include "curve.h"
#include "velocity_function.h"
#include "image.h"

using std::vector;

class dynamics : public DSC2D::VelocityFunc<> {
public:
    dynamics();
    ~dynamics();
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("Image segmentation");
    }
    
    /**
     Computes the motion of each interface vertex and stores the destination in 
     the simplicial complex class.
     */
    virtual void deform(DSC2D::DeformableSimplicialComplex& dsc);
    
public:
    bool update_dsc(DSC2D::DeformableSimplicialComplex &dsc, image &img);
    
    /*
     Split edge with high energy
     */
    void split_edge(dsc_obj &dsc, image &img);
    
#pragma mark - Debug
public:
    void draw_curve(DSC2D::DeformableSimplicialComplex &dsc);
    
#pragma mark - Data
private:
    std::vector<curve> curve_list_;
    dynamics_param d_param_ = g_param;
    
#pragma mark - private
private:
    std::vector<curve> extract_curve(DSC2D::DeformableSimplicialComplex &dsc);
    void compute_internal_force(std::vector<curve> &curve_list,
                                DSC2D::DeformableSimplicialComplex &dsc);
    void compute_external_force(std::vector<curve> &curve_list
                                ,dsc_obj &complex
                                ,image &tex);
    void compute_external_force_edge(std::vector<curve> &curve_list
                                ,dsc_obj &complex
                                ,image &tex);
    void compute_displacement(dsc_obj &dsc);
    
private: // Utility
    curve::node_key get_next_key(curve & cu, std::vector<curve::node_key> pt_keys);
    bool is_in_array(std::vector<curve> &curve_list, Node_key nk);
    bool is_in_array(std::vector<std::vector<Face_key>> &regions, Face_key fk);
    bool is_in_array(std::vector<Face_key> &region, Face_key fk);
    bool is_in_array(std::vector<std::vector<HMesh::HalfEdgeID>> curves_edge, HMesh::HalfEdgeID fk);
    void grow_region(dsc_obj& dsc, std::vector<Face_key> &region, Face_key fk);
    
    HMesh::Walker next_edge(dsc_obj &dsc, HMesh::Walker hew, bool is_CCW);
    HMesh::Walker next_edge_interface(dsc_obj &dsc, HMesh::Walker hew);
    HMesh::HalfEdgeID get_next_edge(dsc_obj &dsc, std::vector<HMesh::HalfEdgeID> newCurve);
};

#endif /* defined(__DSC__dynamics__) */
