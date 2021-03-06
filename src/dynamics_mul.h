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
#include "t_image.h"

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
    
    void  update_dsc(dsc_obj &dsc, t_image &img);
    
    void update_vertex_stable();
    
    void write_energy();

    // temporary variable
    dsc_obj * s_dsc;
    t_image * s_img;
    
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
     Update with adaptive mesh
     */
    void update_dsc_with_adaptive_mesh();
    
    Vec2 get_node_displacement(Node_key vid);
    
    HMesh::VertexAttributeVector<Vec2> internal_node_forces;
    HMesh::VertexAttributeVector<Vec2> external_node_forces;
private:
    void update_probability(dsc_obj &dsc, t_image &img);
    
public:
//    void update_dsc_explicit(dsc_obj &dsc, image &img);
    void compute_mean_intensity();
    void compute_intensity_force();
    void displace_dsc(dsc_obj *obj = nullptr);
    void compute_internal_force();
    
    void compute_difference();

private:
    void adapt_triangle();
    void thinning();

    void thinning_interface();

    double get_energy_assume_label(Face_key fid, int assumed_label);
    
public:
//    void compute_mean_intensity(dsc_obj &dsc, image &img);
    
private: public:
    void get_energy(double &e_, double &l_); // 25.01.16
    
    
    double energy_triangle(HMesh::FaceID fid, double c, int new_phase);
    
//    double optimal_dt(dsc_obj * clone);
    
    Vec2 get_vertex_norm(dsc_obj *obj, HMesh::Walker hew);
    HMesh::Walker pre_edge(dsc_obj *obj, HMesh::Walker hew);
    HMesh::Walker next_edge(dsc_obj *obj, HMesh::Walker hew);
    
    // Energy within a star of the node
//    bool energy_with_location(double &E, Node_key nkey , Vec2 displace, double * real_dis = nullptr);
    
    
    int relabel_triangles(); // Purely base on MS energy

    /******************************/
    // Coarsening approach
    /******************************/
    std::map<int,double> get_energy_thres();
//    void subdivide_triangles();
    
    void subdivide_triangles();

    /******************************/
    // Thinning approach
    /******************************/
    void thinning_triangles(bool colapse_all =  false);
    bool collapse_edge(Edge_key ek, Node_key n_to_remove, bool safe = true);
    
};

#endif /* defined(__DSC__dynamics_mul__) */
