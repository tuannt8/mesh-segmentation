//
//  helper.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/12/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "helper.h"

float helper_t::sign(Vec2 p1, Vec2 p2, Vec2 p3){
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
}
bool helper_t::is_point_in_tri(Vec2 pt, Vec2 t1, Vec2 t2, Vec2 t3){
    bool b1, b2, b3;
    b1 = sign(pt, t1, t2) < 0.0f;
    b2 = sign(pt, t2, t3) < 0.0f;
    b3 = sign(pt, t3, t1) < 0.0f;
    
    return ((b1 == b2) && (b2 == b3));
}

bool helper_t::is_point_in_tri(Vec2 p, std::vector<Vec2> const & pts){
    return is_point_in_tri(p, pts[0], pts[1], pts[2]);
}

double dis_point_to_edge(Vec2 p, Vec2 e1, Vec2 e2){
    // Return minimum distance between line segment vw and point p
    const double l2 = (e1-e2).length();
    if (l2 == 0.0) return (p-e1).length();
    
    const float t = dot(p - e1, e2 - e1) / l2;
    
    if (t < 0.0) return (p - e1).length();
    else if (t > 1.0) return (p - e2).length();
    
    const Vec2 projection = e1 + (e2 - e1)*t;
    
    return (p - projection).length();
}

double helper_t::distance_to_edge(Vec2 p, std::vector<Vec2> const & pts){
    double dis = INFINITY;
    for (int i = 0; i < 3; i++) {
        double dis_to_e = dis_point_to_edge(p, pts[i], pts[(i+1)%3]);
        if(dis > dis_to_e){
            dis = dis_to_e;
        }
    }
    
    return dis;
}

namespace helper_t {
    Vec3 autoColor::next(){
        static std::vector<Vec3> colorList =
        {
            Vec3(1,0,0)     // red
            ,Vec3(0,1,0)    // green
            , Vec3(0,0,1)  // Blue
            , Vec3(1, 1, 0) // yellow
            , Vec3(1, 0, 1) // pink
            , Vec3(0, 1, 1)
            , Vec3(0.3, 0.3, 0.3)
        };
        
        return colorList[index++];
    }
    
//    std::chrono::time_point<std::chrono::system_clock> start_time;
//    std::chrono::time_point<std::chrono::system_clock> stop_time;
    void start_timer(){
//        start_time = std::chrono::system_clock::now();
    }
    void stop_timer(){
//        stop_time = std::chrono::system_clock::now();
    }
    double get_time_and_start(){
//        stop_timer();
//        std::chrono::duration<double> t = (stop_time - start_time);
//        start_timer();
//        return t.count();
        return 0;
    }

    double area(std::vector<Vec2> const &pts) {
        double A = 0.5*cross(pts[1]-pts[0], pts[2]-pts[0]);
        return std::abs(A);
    }
}
