//
//  define.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/11/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef DSC_define_h
#define DSC_define_h
#include "DSC.h"
#include <string>



// Data directory
#ifdef WIN32

#define and &&
#define or ||

#define DATA_PATH "./DATA/"
#define LOG_PATH "./LOG/"

#else
#define DATA_PATH "./DATA/"
#define LOG_PATH "./LOG/"
#endif



#define PI_V1 3.14159

//#define USE_SETTING_FILE

#ifdef USE_SETTING_FILE

extern std::string IMAGE_PATH;
extern int DISCRETIZE_RES;
extern float SMALLEST_SIZE;
extern float SPLIT_FACE_COEFFICIENT;
extern float SPLIT_EDGE_COEFFICIENT;
extern float ALPHA;
extern float DT_;

#define STABLE_MOVE 1e-2
#define BETA    1.0 // external force

#else

/*
 * Two phase synthetics image
 */
// Coefficient
#define IMAGE_PATH "test.png"
// Mesh control
#define DISCRETIZE_RES 10.0
#define SMALLEST_SIZE   10
// Adaptive mesh
#define SPLIT_FACE_COEFFICIENT  0.2 // split thres = coe*(1-coe)*(ci-cj)^2
#define SPLIT_EDGE_COEFFICIENT 0.03    // Split thres = coe*(ci-cj)^2
// Mumford
#define ALPHA 1  // internal force
#define BETA    1.0 // external force

#define DT_ 0.2
//#define ADD_NOISE
#define STABLE_MOVE 1e-2

///*
// * Sound
// */
//// Coefficient
//#define IMAGE_PATH "Data/sound_gap.png"
//// Mesh control
//#define DISCRETIZE_RES 20.0
//#define SMALLEST_SIZE   3.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.08 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 3    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 1.0
//#define STABLE_MOVE 1e-2

///*
// * Fuel cell
// */
//// Coefficient
//#define IMAGE_PATH "Data/fuel_cell/small_gap.png"
//// Mesh control
//#define DISCRETIZE_RES 20.0
//#define SMALLEST_SIZE   3.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.6 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 1    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 0.2
//#define ADD_NOISE
//#define STABLE_MOVE 1e-2

///*
// * Hamster -- easy
// */
//// Coefficient
//#define IMAGE_PATH "Data/Gomu/TallHamster_x_135.png"
//// Mesh control
//#define DISCRETIZE_RES 10.0
//#define SMALLEST_SIZE   2.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.08 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 10    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 0.5
//#define STABLE_MOVE 1e-2

///*
// * Carbon filber
// */
//// Coefficient
//#define IMAGE_PATH "../../Large_data/Filber/filber.png"
//// Mesh control
//#define DISCRETIZE_RES 30.0
//#define SMALLEST_SIZE   1
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.1 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 4    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.3  // internal force
//#define BETA    1.0 // external force
//#define DT_ 1
//#define STABLE_MOVE 1e-2

///*
// * Star
// */
//// Coefficient
//#define IMAGE_PATH "Star.png"
//// Mesh control
//#define DISCRETIZE_RES 10.0
//#define SMALLEST_SIZE   15.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.2 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 20    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 0.5
//#define STABLE_MOVE 1e-1
//#define ADD_NOISE

///*
// * Multi phase
// */
//// Coefficient
//#define IMAGE_PATH "multi_700x700.png"
//// Mesh control
//#define DISCRETIZE_RES 10.0
//#define SMALLEST_SIZE   5.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.2 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 2    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 3
//#define STABLE_MOVE 1e-2
//#define ADD_NOISE

///*
// * adapt mesh
// */
//// Coefficient
//#define IMAGE_PATH "adaptive.png"
//// Mesh control
//#define DISCRETIZE_RES 10.0
//#define SMALLEST_SIZE   3.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.2 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 0.3    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.1  // internal force
//#define BETA    1.0 // external force
//#define DT_ 3
//#define STABLE_MOVE 1e-2
//#define ADD_NOISE

///*
// * Cement
// */
//// Coefficient
//#define IMAGE_PATH "Data/cement/sample.png"
//// Mesh control
//#define DISCRETIZE_RES 20.0
//#define SMALLEST_SIZE   4.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.08 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 10    // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.2  // internal force
//#define BETA    1.0 // external force
//#define DT_ 0.5
//#define STABLE_MOVE 1e-2

///*
// * Dental
// */
//// Coefficient
//#define IMAGE_PATH "Data/dental/sample_gap.png"
//// Mesh control
//#define DISCRETIZE_RES 20.0
//#define SMALLEST_SIZE   8.0
//// Adaptive mesh
//#define SPLIT_FACE_COEFFICIENT  0.5 // split thres = coe*(1-coe)*(ci-cj)^2
//#define SPLIT_EDGE_COEFFICIENT 5   // Split thres = coe*(ci-cj)^2
//// Mumford
//#define ALPHA 0.2  // internal force
//#define BETA    1.0 // external force
//#define DT_ 0.5
//#define STABLE_MOVE 1e-2

#endif

using std::vector;
using std::cout;
using std::endl;

/*********************************************************************/
/* Type def
 */
typedef DSC2D::vec2 Vec2;
typedef DSC2D::vec3 Vec3;
typedef DSC2D::DeformableSimplicialComplex dsc_obj;
typedef dsc_obj::node_key Node_key;
typedef dsc_obj::face_key Face_key;
typedef dsc_obj::edge_key Edge_key;
typedef std::vector<Vec2> Vec2_array;

struct dynamics_param{
    dynamics_param(){}
    dynamics_param(double a, double b, double m){alpha = a; beta = b; mass = m;};
    
    double alpha = 1.0; // Curvature
    double beta = 1.0; // Forth derivative. Keep the curve straight
    double mass = 10.0; // Display scale
    
    std::map<int, double> mean_intensity;
    
    std::vector<bool> bDisplay;
} ;

/*******
 Flags
 */

/*********************************************************************/
/* GLobal variable
 */
extern int debug_num[10];

// Image to load
extern std::string IMAGE_NAME;
// Dynamics parameter
extern dynamics_param g_param;


#endif // File protection
