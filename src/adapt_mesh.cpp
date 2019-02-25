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
        if (
            dsc.min_angle(*fit) < 10*M_PI/180.// dsc.DEG_ANGLE
            && dsc.max_angle(*fit, het) < 120*M_PI/180.
            )
        {
            dsc.remove_degenerate_needle(*fit);
        }
        else if(dsc.max_angle(*fit, het) > 170*M_PI/180.)
        {
            dsc.remove_degenerate_needle2(*fit);
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
void adapt_mesh::collapse_interface(DSC2D::DeformableSimplicialComplex &dsc, t_image &img)
{
    dsc_ = & dsc;

    auto mean_inten_ = g_param.mean_intensity;
    
    for (auto vkey : dsc.vertices())
    {
        if (dsc.is_interface(vkey) && !dsc.is_crossing(vkey))
        {
            std::vector<HMesh::Walker> edges;
            for (auto hew = dsc.walker(vkey); !hew.full_circle(); hew = hew.circulate_vertex_cw())
            {
                if (dsc.is_interface(hew.halfedge()))
                {
                    edges.push_back(hew);
                }
            }
            assert(edges.size()==2);
            
            auto cangle = DSC2D::Util::cos_angle(dsc.get_pos(edges[0].opp().vertex()),
                                         dsc.get_pos(edges[0].vertex()),
                                         dsc.get_pos(edges[1].opp().vertex()));
            if (cangle < dsc.COS_MIN_ANGLE)
            {
                // check energy
                double ev = 0;
                double c0 = mean_inten_[dsc.get_label(edges[0].face())];
                double c1 = mean_inten_[dsc.get_label(edges[0].opp().face())];
                
                auto p0 = dsc.get_pos(edges[0].opp().vertex());
                auto p1 = dsc.get_pos(edges[1].opp().vertex());
                
                double length = (p1 - p0).length();
                int N = (int)length;
                double dl = length/(double)N;
                for (int i = 0; i <= N; i++) {
                    auto p = p0 + (p1 - p0)*(i/(double)N)*dl;
                    double I = img.get_intensity_f(p[0], p[1]);
                    
                    // Normalize force
                    double f = (2*I - c0 - c1) / (c0-c1);
                    
                    ev += std::abs(f)*dl;
                }
                
                ev = ev / (length + SINGULAR_EDGE);
                
                double thres = SPLIT_EDGE_COEFFICIENT*(c0-c1) * (c0-c1);
                if (dsc.bStable[vkey] == 1
                    && ev < thres)
                {
                    if (HMesh::precond_collapse_edge(*dsc.mesh, edges[0].halfedge())
                        && dsc.unsafe_editable(edges[0].halfedge())
                        && (dsc.is_collapsable(edges[0], true) || dsc.is_collapsable(edges[0].opp(), true))
                        )
                    {
                        dsc.collapse(edges[0].halfedge(), 1);
                    }
                }
            }
        }
    }
    
    
}
void adapt_mesh::split_edge(DSC2D::DeformableSimplicialComplex &dsc, t_image &img)
{
    dsc_ = &dsc;
    
    
    std::vector<Edge_key> edges;
    for(auto hei = dsc.halfedges_begin(); hei != dsc.halfedges_end(); ++hei)
    {
        if (dsc.is_interface(*hei)) {
            auto hew = dsc.walker(*hei);
            if(dsc.is_movable(*hei)
               && dsc.get_label(hew.face()) < dsc.get_label(hew.opp().face()))
            {
                edges.push_back(*hei);
            }
        }
    }
    
    auto mean_inten_ = g_param.mean_intensity;
    
    for (auto ekeyp = dsc.halfedges_begin(); ekeyp != dsc.halfedges_end(); ekeyp++){

        auto ekey = *ekeyp;
        auto hew = dsc.walker(ekey);
        
        
        
        if (! dsc.mesh->in_use(*ekeyp)
            || !dsc.is_interface(ekey)
            || hew.vertex() < hew.opp().vertex()
            || hew.face() == HMesh::InvalidFaceID
            || hew.opp().face() == HMesh::InvalidFaceID)
        {
            continue;
        }
        
        double ev = 0;
        double c0 = mean_inten_[dsc.get_label(hew.face())];
        double c1 = mean_inten_[dsc.get_label(hew.opp().face())];
        
        // Loop on the edge
        auto p0 = dsc.get_pos(hew.opp().vertex());
        auto p1 = dsc.get_pos(hew.vertex());

        double length = (p1 - p0).length();
        int N = (int)length;
        double dl = length/(double)N;
        for (int i = 0; i <= N; i++) {
            auto p = p0 + (p1 - p0)*(i/(double)N)*dl;
            double I = img.get_intensity_f(p[0], p[1]);
            
            // Normalize force
            double f = (2*I - c0 - c1) / (c0-c1);
            
            ev += std::abs(f)*dl;
        }
        
        ev = ev / (length + SINGULAR_EDGE);
        
        double thres = SPLIT_EDGE_COEFFICIENT*(c0-c1) * (c0-c1);
        
        if (dsc.bStable[hew.vertex()] == 1
            && dsc.bStable[hew.opp().vertex()] == 1)
        {
            if (ev > thres && length > 3*SMALLEST_SIZE
                && !is_bound(&dsc, ekey)
                ) // High energy. Split
            {
                
                
                dsc.split_adpat_mesh(ekey);
            }
            else // Low energy, consider collapse
            {
                // Only collapse edge on the interface
                // And doesnot reduce mesh quality
                // conflict with face split
                
//                if(dsc.collapse(ekey, true))
//                {
//                }
            }
        }
    }
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







