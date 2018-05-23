//
//  polygon_f.h
//  DSC_seg
//
//  Created by Tuan Nguyen Trung on 5/4/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC_seg__polygon_f__
#define __DSC_seg__polygon_f__

#include <stdio.h>
#include "define.h"

// f = PI(x - x_i)
// The abstract class
class polygon_f{
public:
    std::vector<double> x_i;
    
    
public:
    polygon_f();
    ~polygon_f();
    
    double operator ()(double x);
    double at(double x);
    double diff(double x);
};

// XI_i = PI(x - j)/ with j != i
class XI{
    polygon_f numerator_;
    double denuminator_;
    int n_;
    int i_;
public:
    XI(int i, int n);
    ~XI();
    
    double at(double x);
    double diff(double x);
};

class KPhi: public polygon_f{
    polygon_f f_;
    int n_;
public:
    KPhi(int n);
    ~KPhi();
    double at(double x);
    double diff(double x);
};

#endif /* defined(__DSC_seg__polygon_f__) */
