//
//  setting_io.cpp
//  DSC_seg_integral_cinema
//
//  Created by Tuan Nguyen Trung on 8/16/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#include "setting_io.hpp"


//#ifdef USE_SETTING_FILE

std::string IMAGE_PATH;
int DISCRETIZE_RES;
float SMALLEST_SIZE;
float SPLIT_FACE_COEFFICIENT;
float SPLIT_EDGE_COEFFICIENT;
float ALPHA;
float DT_;
bool RELABEL = true;
std::map<int , std::vector<init_circle>> _circle_inits;
int NB_PHASE = 0;
int ADAPTIVE = 1;

//#endif

std::string setting_io::_file_path;

void setting_io::setup_parameter()
{
    IMAGE_PATH = get_string(IMAGE_PATH_S);
    DISCRETIZE_RES = get_int(DISCRETIZE_RES_S);
    SMALLEST_SIZE = get_float(SMALLEST_SIZE_S);
    ALPHA = get_float(ALPHA_S);
    DT_ = get_float(DT_S);
    RELABEL = get_int(RELABEL_S);
    NB_PHASE = get_int(NB_PHASE_S);
    ADAPTIVE = get_int(ADAPTIVE_S);
}

