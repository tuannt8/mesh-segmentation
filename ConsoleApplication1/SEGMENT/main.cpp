//
//  SEGMENT.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/9/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//


#include "interface.h"
//#include "DSC.h"

int main(int argc, char ** argv) {
    interface_dsc ui(argc, argv);
    std::cout<<"GL ver: " << glGetString(GL_VERSION)<<std::endl;
    glutMainLoop();
    return 0;
}
