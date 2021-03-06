//
//  adapt_mesh.h
//  DSC_seg_integral
//
//  Created by Tuan Nguyen Trung on 7/6/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC_seg_integral__adapt_mesh__
#define __DSC_seg_integral__adapt_mesh__

#include <stdio.h>
#include "DSC.h"
#include "define.h"
#include "t_image.h"

#define BOUND_FACE 0

#define SINGULAR_AREA 0.        // In computation of triangle variation
#define SINGULAR_EDGE (SMALLEST_SIZE / 2.0)


class adapt_mesh{
public:
    adapt_mesh();
    ~adapt_mesh();
    
    // Split edge manually
    void split_edge(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
    
    // Split face mannually
    void collapse_interface(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
    void split_face(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
    
    // remove steiner vertices
    void adapt_triangle(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
    void thinning(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
    
    void remove_needles(DSC2D::DeformableSimplicialComplex &dsc);
    
    // First step: Split face and relabeling. Not working.
//    void split_face_and_relabel(DSC2D::DeformableSimplicialComplex &dsc, image &img);
private:
    DSC2D::DeformableSimplicialComplex *dsc_;
    void split_single_edge(Edge_key ekey);
    void add_point_if_need(HMesh::Walker hew);
    bool collapse_edge(HMesh::Walker hew);
    
public:// From dense to sparse, another approach
    void coarsening_triangles(DSC2D::DeformableSimplicialComplex &dsc);
    void coarsening_interface(DSC2D::DeformableSimplicialComplex &dsc, t_image &img);
};

#endif /* defined(__DSC_seg_integral__adapt_mesh__) */
