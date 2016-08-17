//
//  dynamics.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "dynamics.h"
#include "helper.h"

#define IS_IN_ARRAY_LIST(a, b, c) \
{ \
    c = false; \
    for (auto & it : a) { \
        if(c) break; \
        for (auto &n : it) { \
            if (n == b) { \
                c = true; break; \
            } \
        } \
    } \
}



bool dynamics:: update_dsc(DSC2D::DeformableSimplicialComplex &dsc, image &img){
    // 1. Process interface vertices
    curve_list_ = extract_curve(dsc);
    assert(curve_list_.size() > 0);

    // 2. Internal forces
    compute_internal_force(curve_list_, dsc);
    
    // 3. External forces
    compute_external_force_edge(curve_list_, dsc, img);

    // 4. Compute displacement
    compute_displacement(dsc);
    
    // 5. Update DSC
    deform(dsc);
    
//    // Debug
//    // 1. Process interface vertices
    curve_list_ = extract_curve(dsc);
//
//    // 2. Internal forces
    compute_internal_force(curve_list_, dsc);
//
//    // 3. External forces
    compute_external_force_edge(curve_list_, dsc, img);

    
    return  true;
}

void dynamics::deform(DSC2D::DeformableSimplicialComplex& dsc)
{
    auto init_time = std::chrono::system_clock::now();
    update_compute_time(init_time);
    init_time = std::chrono::system_clock::now();
    
    dsc.deform();
    
    update_deform_time(init_time);
}

dynamics::dynamics(): VelocityFunc(0.1, 0.01){
    
}

dynamics::~dynamics(){
    
}

bool dynamics::is_in_array(std::vector<curve> &curve_list, Node_key nk){
    for (auto & c : curve_list) {
        for (auto &n : c) {
            if (n == nk) {
                return true;
            }
        }
    }
    return false;
}

bool dynamics::is_in_array(std::vector<std::vector<Face_key>> &regions, Face_key fk){
    for (auto & fs : regions) {
        for (auto f : fs) {
            if (f == fk) {
                return true;
            }
        }
    }
    return false;
}

bool dynamics::is_in_array(std::vector<std::vector<HMesh::HalfEdgeID>> curves_edge, HMesh::HalfEdgeID fk){
    bool out;
    IS_IN_ARRAY_LIST(curves_edge, fk, out);
    return out;
}

bool dynamics::is_in_array(std::vector<Face_key> &region, Face_key fk){
    for (auto f : region) {
        if (f == fk) {
            return true;
        }
    }
    return false;
}

void dynamics::grow_region(dsc_obj& dsc, std::vector<Face_key> &region, Face_key fk){
    // Check if this triangle existed
    if (is_in_array(region, fk)) {
        return;
    }
    
    region.push_back(fk);
    
    // Get all neightbor share edge
    for (auto hew = dsc.walker(fk); hew.full_circle(); hew = hew.circulate_face_cw()) {
        auto nb = hew.opp().face();
        if (dsc.get_label(nb) != 0) {
            grow_region(dsc, region, nb);
        }
    }
}

HMesh::Walker dynamics::next_edge(dsc_obj &dsc, HMesh::Walker hew, bool is_CCW){
    
    if (is_CCW) {
        return hew.next().opp();
    }else{
        return hew.opp().prev();
    }
}

HMesh::Walker dynamics::next_edge_interface(dsc_obj &dsc, HMesh::Walker hew){
#ifdef DEBUG
    assert(dsc.is_interface(hew.halfedge()));
#endif
    bool isCCW = dsc.get_label(hew.face()) != 0;
    
    while (1) {
        hew = next_edge(dsc, hew, isCCW);
        if (dsc.is_interface(hew.halfedge())) {
            return hew.opp();
        }
    }
}

std::vector<curve> dynamics::extract_curve(DSC2D::DeformableSimplicialComplex &dsc){
    
    std::vector<std::vector<HMesh::HalfEdgeID>> curves_edge;
    std::vector<curve> curve_list;

    for(auto hei = dsc.halfedges_begin(); hei != dsc.halfedges_end(); ++hei)
    {
        auto hew = dsc.walker(*hei);
        if (dsc.is_interface(hew.halfedge())
            && !is_in_array(curves_edge, hew.halfedge())
            && !is_in_array(curves_edge, hew.opp().halfedge()))
        {
            // Start growing this guy
            std::vector<HMesh::HalfEdgeID> newCurve;
            newCurve.push_back(hew.halfedge());
            curve new_curve;
            new_curve.push_back(hew.vertex());
            
            while (1) {
                hew =  next_edge_interface(dsc, hew);
                assert(dsc.is_interface(hew.halfedge()));
                if (hew.halfedge() == newCurve[0]) {
                    break;
                }
                else{
                    newCurve.push_back(hew.halfedge());
                    new_curve.push_back(hew.vertex());
                }
            }
            
            curves_edge.push_back(newCurve);
            curve_list.push_back(new_curve);
        }
    }
    
    return curve_list;
}

curve::node_key dynamics::get_next_key(curve & cu, std::vector<curve::node_key> pt_keys){
    if (cu.size() == 1) {
        return pt_keys[0];
    }else{
        return ( pt_keys[0] == cu[cu.size() - 2] )? pt_keys[1] : pt_keys[0];
    }
    assert(0);
}

void dynamics::draw_curve(DSC2D::DeformableSimplicialComplex &dsc){
    helper_t::autoColor co;
    for (int i = 0; i < curve_list_.size(); i++) {
        curve_list_[i].draw(dsc, co.next());
    }
}

void dynamics::compute_internal_force(std::vector<curve> &curve_list,
                            DSC2D::DeformableSimplicialComplex &dsc){
    for (auto & cu : curve_list) {
        cu.update_derivative(dsc);
        for (int i = 0; i < cu.size(); i++) {
//            if (dsc.is_crossing(cu[i])) {
//                continue; // By pass crossing vertex
//            }
            Vec2 force = cu.derive2(i)*d_param_.alpha - cu.derive4(i)*d_param_.beta;
            dsc.set_node_internal_force(cu[i], force);
        }
    }
}

void dynamics::compute_external_force_edge(std::vector<curve> &curve_list
                                 ,dsc_obj &complex
                                 ,image &img){
    for (auto &cu : curve_list) {
        cu.update_mean_intensity(complex, img);
        
        std::vector<Edge_key> edge_list = cu.get_edge_list(complex);
        std::vector<Vec2> node_force = cu.get_node_force(complex, img);
        
        for (int i = 0; i < node_force.size(); i++) {
            complex.set_node_external_force(cu[i], -node_force[i]*g_param.gamma);
        }
    }
}

void dynamics::compute_external_force(std::vector<curve> &curve_list
                                      ,dsc_obj &complex
                                      ,image &img){
    for(auto &cu : curve_list){
        cu.update_mean_intensity(complex, img);
        
        for (int i = 0; i < cu.size(); i++) {
            
            Vec2 pt = complex.get_pos(cu[i]);
            double inten = img.get_intensity_i(pt[0], pt[1]);
            double scale = (cu.m_out() - cu.m_in())* (inten - cu.m_out() + inten - cu.m_in());
            
            Vec2 norm = complex.get_normal(cu[i]);
            Vec2 force = -norm*scale * d_param_.gamma;
            
            //Vec2 norm = img.get_local_norm(complex, cu[i], scale < 0);
            //Vec2 force = norm*std::abs(scale) * d_param_.gamma;
            
            complex.set_node_external_force(cu[i], force);
        }
    }

}

void dynamics::compute_displacement(dsc_obj &dsc){

    double el = dsc.get_avg_edge_length();
    for (auto ni = dsc.vertices_begin(); ni != dsc.vertices_end(); ni++) {
        Vec2 dis = (dsc.get_node_internal_force(*ni) + dsc.get_node_external_force(*ni))
                        / d_param_.mass;
        
        if (dis.length() > 0.0001*el
       //     && !dsc.is_crossing(*ni)
            ) {
            dsc.set_destination(*ni, dsc.get_pos(*ni) + dis);
        }
        
    }
}

void dynamics::split_edge(dsc_obj &dsc, image &img) {
    std::vector<Edge_key> edges;
    for(auto hei = dsc.halfedges_begin(); hei != dsc.halfedges_end(); ++hei)
    {
        auto hew = dsc.walker(*hei);
        if(dsc.is_movable(*hei) && dsc.get_label(hew.face()) < dsc.get_label(hew.opp().face()))
        {
            edges.push_back(*hei);
        }
    }
    
    int thres_hold = 10;
    
    std::map<Edge_key, bool> edge_used;
    for (auto e : edges)
    {
        if (edge_used[e]) {
            continue;
        }

        
        auto hew = dsc.walker(e);
        if(dsc.get_label(hew.face()) != 0)
            hew = hew.opp();
        
        edge_used.insert(std::pair<Edge_key, bool>(e, true));
        edge_used.insert(std::pair<Edge_key, bool>(hew.opp().halfedge(), true));
        
        int count_out;
        int out = img.get_triangle_intensity_count(dsc.get_pos(hew.face()), &count_out);
        
        if(out/(double)count_out > thres_hold)
        {
            bool success = dsc.split(e);
#ifdef DEBUG
            if(success)
            {
                std::cout << "Split interface" << std::endl;
            }
#endif
        }
    }
    
    curve_list_ = extract_curve(dsc);
}
