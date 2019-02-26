//
//  adapt_mesh.cpp
//  DSC_seg_integral
//
//  Created by Tuan Nguyen Trung on 7/6/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "adapt_mesh.h"
#include "dynamics_mul.h"



adapt_mesh::adapt_mesh(){
}

adapt_mesh::~adapt_mesh(){
    
}

#define PROTECT_BOUND






struct edge_s_e
{
    edge_s_e(Edge_key ekey_, double length_):ekey(ekey_), length(length_){}
    Edge_key ekey;
    double length;
};

void adapt_mesh::remove_needles(DSC2D::DeformableSimplicialComplex &dsc)
{
    // TODO: Check inverted triangles when collapsing the two edge
    dsc_ = &dsc;
    
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (!dsc.mesh->in_use(*fit))
        {
            continue;
        }
        
        HMesh::Walker het = dsc.walker(*fit);
        
        double thres_hold  = 0.7 * dsc.MIN_ANGLE;
        
        if (dsc.min_angle(*fit) < thres_hold) //dsc.MIN_ANGLE
        {
            if( dsc.max_angle(*fit, het) < 120*M_PI/180.)
            {
                dsc.remove_degenerate_needle(*fit);
            }
            else if(dsc.max_angle(*fit, het) > 170*M_PI/180.)
            {
                dsc.remove_degenerate_needle2(*fit);
            }
        }
    }
}

//void adapt_mesh::thinning(DSC2D::DeformableSimplicialComplex &dsc, image &img)
//{
//    dsc_ = & dsc;

    
//    double flip_thres = SPLIT_FACE_COEFFICIENT;
    
////    std::vector<Node_key> to_collapse;
//    for (auto nkey : dsc_->vertices())
//    {
//        if (dsc_->is_interface(nkey)
//            || HMesh::boundary(*dsc_->mesh, nkey))
//        {
//            continue;
//        }
        
//        // Check if the one ring have low variation
//        bool low_var = true;
//        auto smallest = dsc_->walker(nkey);
//        double shortest = INFINITY;
//        for (auto hew = dsc_->walker(nkey); !hew.full_circle(); hew = hew.circulate_vertex_ccw())
//        {
//            auto fkey = hew.face();
//            auto pts = dsc_->get_pos(fkey);
//            double area;
//            double mi = img.get_tri_intensity_f(pts, &area); mi /= area;
//            double e = img.get_tri_differ_f(pts, mi)/ (area + SINGULAR_AREA);
            
//            if (e > flip_thres)
//            {
//                low_var = false;
//            }
            
//            if (dsc_->length(hew.halfedge()) < shortest)
//            {
//                shortest = dsc_->length(hew.halfedge());
//                smallest = hew;
//            }
//        }
        
//        if (low_var)
//        {
//            dsc_->collapse(smallest.halfedge(), true);
//            //dsc.collapse(smallest, 0.0);
//        }
//    }
    

//}

inline bool is_bound(DSC2D::DeformableSimplicialComplex * dsc, HMesh::HalfEdgeID e)
{
  //  return false;
    
    auto hew = dsc->walker(e);
    if (dsc->get_label(hew.face()) == BOUND_FACE
        || dsc->get_label(hew.opp().face()) == BOUND_FACE)
    {
        return true;
    }
    
    return false;
}


bool adapt_mesh::collapse_edge(HMesh::Walker hew)
{
    if (dsc_->is_crossing(hew.opp().vertex()) )
    {
        return false;
    }
    
    Vec2 p0 = dsc_->get_pos(hew.opp().vertex());
    Vec2 p1 = dsc_->get_pos(hew.vertex());
    Vec2 p2 = dsc_->get_pos(dsc_->previous_interface(hew).opp().vertex());
    if (DSC2D::Util::cos_angle(p2, p0, p1) < dsc_->COS_MIN_ANGLE)
    {
        return dsc_->collapse(hew, 0.0);
    }
    
    return false;
}

void adapt_mesh::coarsening_triangles(DSC2D::DeformableSimplicialComplex &dsc)
{

}

void adapt_mesh::coarsening_interface(DSC2D::DeformableSimplicialComplex &dsc, t_image &img)
{
    
}

void adapt_mesh::split_single_edge(Edge_key ekey)
{
    auto hew = dsc_->walker(ekey);
    
    // Add point to adjencent triangle that has small cap
    add_point_if_need(hew);
    add_point_if_need(hew.opp());
}

void adapt_mesh::add_point_if_need(HMesh::Walker hew)
{
    auto p1 = dsc_->get_pos(hew.opp().vertex());
    auto p2 = dsc_->get_pos(hew.vertex());
    auto p0 = dsc_->get_pos(hew.next().vertex());
    
    if(DSC2D::Util::cos_angle(p1, p0, p2) > std::cos(60*180/PI_V1))
    {
        // Add point
        
    }
}







