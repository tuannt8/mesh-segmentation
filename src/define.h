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

#include <GEL/CGLA/Vec2i.h>

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

/*********************************************************************/
/* Type def
 */
typedef DSC2D::vec2 Vec2;
typedef CGLA::Vec2i Vec2i;
typedef DSC2D::vec3 Vec3;
typedef DSC2D::DeformableSimplicialComplex dsc_obj;
typedef dsc_obj::node_key Node_key;
typedef dsc_obj::face_key Face_key;
typedef dsc_obj::edge_key Edge_key;
typedef std::vector<Vec2> Vec2_array;


#define PI_V1 3.14159

struct init_circle
{
    Vec2 _center;
    double _radius;
    
    bool is_in_circle(Vec2 pt){ return (pt-_center).length() < _radius; }
    bool is_in_circle(std::vector<Vec2> pts)
    {
        for(auto p: pts){
            if(!is_in_circle(p))
                return false;
        }
        return true;
    }
    
    Vec2 project_to_circle(Vec2 pt)
    {
        return _center + (pt - _center)*( _radius / (pt - _center).length());
    }
};

extern std::string IMAGE_PATH;
extern int DISCRETIZE_RES;
extern float SMALLEST_SIZE;
extern float SPLIT_FACE_COEFFICIENT;
extern float SPLIT_EDGE_COEFFICIENT;
extern float ALPHA;
extern float DT_;
extern bool RELABEL;
extern int NB_PHASE;
extern int ADAPTIVE;
extern std::map<int, std::vector<init_circle>> _circle_inits;

#define STABLE_MOVE 1e-2 // ???

using std::vector;
using std::cout;
using std::endl;



struct dynamics_param{
    dynamics_param(){}
        
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
//extern std::string IMAGE_NAME;
// Dynamics parameter
extern dynamics_param g_param;


#endif // File protection
