//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#pragma once

#include "DSC.h"
#include "velocity_function.h"
#include "define.h"

#include <map>

#include "t_image.h"

#define FORCE_SCALE 10

const static double POINT_SIZE = 0.05;
const static double LINE_WIDTH = 0.03;

const static DSC2D::vec3 BACKGROUND_COLOR = DSC2D::vec3(0.7); // DSC2D::vec3(1.0);
const static DSC2D::vec3 INVISIBLE = DSC2D::vec3(-1.);
const static DSC2D::vec3 DARK_RED = DSC2D::vec3(0.66,0.11,0.15);
const static DSC2D::vec3 RED = DSC2D::vec3(0.96,0.11,0.15);
const static DSC2D::vec3 DARK_BLUE = DSC2D::vec3(0.14,0.16,0.88);
const static DSC2D::vec3 BLUE = DSC2D::vec3(0.45,0.7,0.9);
const static DSC2D::vec3 GREEN = DSC2D::vec3(0.05,1.,0.15);
const static DSC2D::vec3 ORANGE = DSC2D::vec3(0.9,0.4,0.);
const static DSC2D::vec3 BLACK = DSC2D::vec3(0.);
const static DSC2D::vec3 DARK_GRAY = DSC2D::vec3(0.5);
const static DSC2D::vec3 GRAY = DSC2D::vec3(0.8);



/**
 A painter handles all draw functionality using OpenGL.
 */
class Painter {
    
    
public:
    /**
     Saves the current painting to the selected folder.
     */
    static void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1);
    
    static void save_painting_no_overwite(int width, int height, std::string folder = std::string(""));
    
    static void save_painting_dsc(int originx, int originy, int width, int height, std::string folder); 
    /**
     Begins drawing.
     */
    static void begin();
    
    /**
     Ends drawing.
     */
    static void end();
    
    /**
     Draw force
     */
    static void draw_internal_force(const DSC2D::DeformableSimplicialComplex& complex, double scale = 1.0);
    static void draw_external_force(const DSC2D::DeformableSimplicialComplex& complex, double scale = 1.0);
    
    /**
     Draws the simplicial complex.
     */
    static void draw_complex(const DSC2D::DeformableSimplicialComplex& complex);
    
    /**
     Draws the domain.
     */
    static void draw_domain(const DSC2D::DesignDomain& domain, DSC2D::vec3 color = GRAY);
    
    /**
     Draws the vertices with the colors defined by the get_vertex_colors function in the simplicial complex.
     */
    static void draw_vertices(const DSC2D::DeformableSimplicialComplex& complex);
    static void draw_vertices_index(const DSC2D::DeformableSimplicialComplex& complex);
    /**
     Draws the edges with the colors defined by the get_edge_colors function in the simplicial complex.
     */
    static void draw_edges(const DSC2D::DeformableSimplicialComplex& complex);
    static void draw_edges_index(const DSC2D::DeformableSimplicialComplex& complex);
    
    /**
     Draws the faces with the colors defined by the get_face_colors function in the simplicial complex.
     */
    static void draw_faces(const DSC2D::DeformableSimplicialComplex& complex);
    static void draw_faces_index(const DSC2D::DeformableSimplicialComplex& complex);
    // Draw faces with color from phase intensity
    static void draw_faces_intensity(const DSC2D::DeformableSimplicialComplex& complex);
    /**
     Draws the faces with the colors given as input.
     */
    static void draw_faces(const DSC2D::DeformableSimplicialComplex& complex, const HMesh::FaceAttributeVector<DSC2D::vec3> &colors);
    
    /**
     Draw the phase index (Phase label)
     */
    static void draw_face_label(const DSC2D::DeformableSimplicialComplex& complex);
    
    /**
     Draws the faces using the 'jet' color scheme.
     */
    static void draw_faces(const DSC2D::DeformableSimplicialComplex& complex, const HMesh::FaceAttributeVector<double> &values);
    
    /**
     Draws the interface with the color given as input.
     */
    static void draw_interface(const DSC2D::DeformableSimplicialComplex& complex, DSC2D::vec3 color = ORANGE);
    
    /**
     Draws the arrows given as input with the color given as input.
     */
    static void draw_arrows(const DSC2D::DeformableSimplicialComplex& complex, const HMesh::VertexAttributeVector<DSC2D::vec2> &arrows, DSC2D::vec3 color = ORANGE, double scale = 1.0);

    /**
     Draws the lines given as input with the color given as input.
     */
    static void draw_lines(const DSC2D::DeformableSimplicialComplex& complex, const HMesh::VertexAttributeVector<DSC2D::vec2> &lines, DSC2D::vec3 color = GREEN);
    
    // Utility
    static void print_gl(const double &x, const double &y, const char* str);
    
    static void draw_face_energy(const DSC2D::DeformableSimplicialComplex& dsc, t_image & img, std::map<int,double> thres_hold);
};
