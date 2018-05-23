//
//  dynamics_mul.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 3/20/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC__dynamics_mul__
#define __DSC__dynamics_mul__

#include <stdio.h>
#include "define.h"
#include "image.h"

//#define EX_BOUND 1
//#define IN_BOUND 2
//#define IMAGE_GRAD 3
//#define INDEX_VERT 4
//#define FORCE_TEMP 5
//#define V_COUNT 6
//#define E2 7

enum {
    EX_BOUND = 0,
    IN_BOUND,
    IMAGE_GRAD,
    INDEX_VERT,
    FORCE_TEMP,
    V_COUNT,
    E2_GRAD,
    AREA_FORCE,
    STAR_DIFFER,
};

enum {
    FACE_IDX = 0,
    PHASE_PROBABILITY, // Probability of phase
};

// #define STABLE_MOVE 1e-2

class dynamics_mul {

    
public:
    dynamics_mul();
    ~dynamics_mul();
    
    void  update_dsc(dsc_obj &dsc, image &img);
    
    void update_vertex_stable();
    
    void write_energy();
private:
    // temporary variable
    dsc_obj * s_dsc;
    image * s_img;
    
    // Mean intensity
    std::map<int, double> mean_inten_;
    std::map<int, double> alpha_map_;

    // Measure energy change
    FILE *f;
    std::vector<Vec3> data_log;
    
    // Adaptive dt
    double E0_ = 0.0, E1_ = 0.0, dE_0_ = 0., dE2 = 0.;
    std::vector<Vec2> E_grad0_;
    double dt = DT_;
private:
    /*
     Compute on the whole domain
     */
    void update_dsc_implicit(dsc_obj &dsc, image &img);
    void compute_intensity_force_implicit();
    void compute_curvature_force_implicit();
    void compute_image_gradient_force_implicit(std::vector<Vec2> & grad_force);
    void indexing_vertices();
    std::vector<int> get_vert_idx(std::vector<HMesh::VertexID> vids);
    void build_and_solve();
    
    /*
     Update with adaptive mesh
     */
    void update_dsc_with_adaptive_mesh();
private:
    void update_dsc_explicit_whole_domain(dsc_obj &dsc, image &img);
    void update_dsc_area(dsc_obj &dsc, image &img);
    void update_dsc_build_and_solve(dsc_obj &dsc, image &img);
    
private:
    void update_probability(dsc_obj &dsc, image &img);
    
private:
    void update_dsc_explicit(dsc_obj &dsc, image &img);
    void compute_mean_intensity(std::map<int, double> & mean_inten_o);
    void compute_intensity_force();
    void displace_dsc(dsc_obj *obj = nullptr);
    void compute_curvature_force();
    
    void compute_difference();
    
public:
    void compute_mean_intensity(dsc_obj &dsc, image &img);
    
public:
    void displace_dsc_2();
    void debug_optimum_dt();
    void debug_optimum_dt_2();
    
    void optimize_phase();
    
    // Phase relabeling with small variant condition
    void optimize_phase_with_variant();
    
    double furthest_move(Node_key nid, Vec2 direction);
    double energy_change(Node_key nid, Vec2 new_pos);
    
    double star_energy(Node_key nid, Vec2 new_pos);
    // intensity different in node link
    double intensity_energy(Node_key nid, Vec2 new_pos);
    // Interface length in node link
    double curve_length(Node_key nid, Vec2 new_pos);
    
private: public:
    double get_total_energy(); // Oct. 20
    void get_energy(double &e_, double &l_); // 25.01.16
    
    double get_total_energy(dsc_obj *obj, std::map<int, double>  intesity_map);
    
    double energy_gradient_by_moving_distance(dsc_obj *obj, std::map<int, double>  intesity_map);
    double u_gradient(dsc_obj *obj, std::map<int, double>  intensity_map);
    double image_gradient_count(dsc_obj *obj, std::map<int, double>  intensity_map);
    double gradient_length(dsc_obj *obj);
    void get_curvature(dsc_obj *obj, HMesh::Walker hew, double &Kcur, double &Kpre);
    double get_curvature(dsc_obj *obj, HMesh::Walker hew);
    
    double energy_triangle(HMesh::FaceID fid, double c, int new_phase);
    
    double optimal_dt(dsc_obj * clone);
    
    Vec2 get_vertex_norm(dsc_obj *obj, HMesh::Walker hew);
    HMesh::Walker pre_edge(dsc_obj *obj, HMesh::Walker hew);
    HMesh::Walker next_edge(dsc_obj *obj, HMesh::Walker hew);
    
    // Energy within a star of the node
    bool energy_with_location(double &E, Node_key nkey , Vec2 displace, double * real_dis = nullptr);
};

#endif /* defined(__DSC__dynamics_mul__) */
