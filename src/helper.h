//
//  helper.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC__helper__
#define __DSC__helper__

#include <stdio.h>
#include "define.h"

namespace helper_t {
    // Check point inside triangle
    float sign(Vec2 p1, Vec2 p2, Vec2 p3);
    
    class autoColor{
    public:
        autoColor(){};
        ~autoColor(){};
        int index = 0;
        Vec3 next();
    };

    /**
    * Profiling
    */
    void start_timer();
    void stop_timer();
    double get_time_and_start();

    /**
    * Geometry
    */
    double area(std::vector<Vec2> const & pts);
    bool is_point_in_tri(Vec2 p, std::vector<Vec2> const & pts);
    bool is_point_in_tri(Vec2 p, Vec2 t1, Vec2 t2, Vec2 t3);
    double distance_to_edge(Vec2 p, std::vector<Vec2> const & pts);
}

#endif /* defined(__DSC__helper__) */
