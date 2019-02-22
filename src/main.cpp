//
//  SEGMENT.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/9/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//


#include "interface.h"
//#include "DSC.h"

using std::cout;

void print_help()
{
    cout << "===============================================\n"
        << "| > Image segmentation with triangle mesh\n"
        << "| > Arguments\n"
        << "| > dsc_seg setting_file\n"
        << "| >------------------------\n"
        << "| > Short key:\n"
        << "| >  - 0-9: Choose phase to label triangle\n"
        << "===============================================\n"
    ;
}

int main(int argc, char ** argv) {
    
    print_help();
    
    interface_dsc ui(argc, argv);
    glutMainLoop();
    return 0;
}
