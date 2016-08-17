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

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices in the normal direction.
 */
class NormalFunc: public DSC2D::VelocityFunc<> {
    
    
public:
    /**
     Creates a velocity function which moves the interface vertices in the normal direction.
     */
    NormalFunc(double velocity, double accuracy): VelocityFunc(velocity/100., accuracy/100.)
    {
        
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("NORMAL MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the destination in the simplicial complex class.
     */
    virtual void deform(DSC2D::DeformableSimplicialComplex& dsc);
};
