//
//  polygon_f.cpp
//  DSC_seg
//
//  Created by Tuan Nguyen Trung on 5/4/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "polygon_f.h"

double polygon_f::operator ()(double x){
    double f = 1.0;
    for (auto const xi : x_i) {
        f *= (x - xi);
    }
    return f;
}
double polygon_f::at(double x){
    double f = 1.0;
    for (auto const xi : x_i) {
        f *= (x - xi);
    }
    return f;
}
double polygon_f::diff(double x){
    double d = 0.0;
    for (int i = 0; i < x_i.size(); i++) {
        double di = 1.0;
        for (int j = 0; x_i.size(); j++) {
            if (i != j) {
                di *= (x - x_i[j]);
            }
        }
        d += di;
    }
    
    return d;
}

polygon_f::polygon_f(){
    
}

polygon_f::~polygon_f(){
    
}

XI::XI(int index, int n):
i_(index),
n_(n)
{
    std::vector<double> numer;
    for (int i = 0; i < n; i++) {
        if (i != i_) {
            numer.push_back(i);
        }
    }
    numerator_.x_i = numer;
    
    denuminator_ = 1.0;
    for (int i = 0; i < n; i++) {
        if (i != i_) {
            denuminator_ *= (i-i_);
        }
    }
}

XI::~XI(){
    
}

double XI::at(double x){
    return numerator_.at(x)/denuminator_;
}
double XI::diff(double x){
    return numerator_.diff(x)/denuminator_;
}

KPhi::KPhi(int n){
    n_ = n;
    std::vector<double> numer;
    for (int i = 0; i < n; i++) {
        numer.push_back(i);
    }
    f_.x_i = numer;
}
KPhi::~KPhi(){
    
}
double KPhi::at(double x){
    return f_.at(x);
}
double KPhi::diff(double x){
    return f_.diff(x);
}