//
//  setting_io.hpp
//  DSC_seg_integral_cinema
//
//  Created by Tuan Nguyen Trung on 8/16/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef setting_io_hpp
#define setting_io_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "define.h"

#ifdef WIN32
#include <ctype.h>
#endif


#define IMAGE_PATH_S                "file_path"

#define SMALLEST_SIZE_S             "shortest_edge"
#define ALPHA_S                     "smooth_coefficient"
#define SPLIT_EDGE_COEFFICIENT_S        "threshold_adapt_edge"
#define SPLIT_FACE_COEFFICIENT_S        "threshold_adapt_face"

#define RELABEL_S               "relabel"

#define DT_S                        "time_step"
#define DISCRETIZE_RES_S           "dsc_init_resolution"

#define NB_PHASE_S            "number_phases_to_segment"
#define ADAPTIVE_S            "adaptive_mesh"

#define REMOVE_SPACE(line) \
    line.erase(std::remove_if( \
                        line.begin(), \
                        line.end(), \
                        [](char x){return isspace(x);}), \
                        line.end());

class setting_io
{
    static std::string _file_path;
    
public:
    setting_io(){}
    
    ~setting_io(){}
    
    static void set_param(std::string file_path)
    {
        setting_io a;
        a.load(file_path);
    }
    
private:
    void load(std::string file_path)
    {
        _file_path = file_path;
        
        try
        {
            std::ifstream infile(file_path);
            std::string line;
            std::string delimiter = ":";
            
            while(std::getline(infile, line))
            {
                auto origin = line;
                                
                // Ignore comment with #
                if (line.size()==0 || line[0] == '#')
                {
                    continue;
                }
                
                line = line.substr(0, line.find('#'));

                
                // Read key and value
                std::string key = line.substr(0, line.find(delimiter));
                std::string value = line.substr(key.size()+1, line.size());
                
                REMOVE_SPACE(key);
                auto tt = key.substr(0, 5);
                
                if (key.substr(0, 5) == "phase")
                {
                    int p = std::atoi(key.substr(5,6).c_str()); // phase
                    int x, y;
                    double r;
                    std::istringstream is(value.c_str());
                    is >> x >> y >> r;
                    
                    if (_circle_inits.find(p) == _circle_inits.end())
                    {
                        std::vector<init_circle> a = {{Vec2(x,y), r}};
                        _circle_inits.insert(std::make_pair(p, a) );
                    }else{
                        init_circle a = {Vec2(x,y), r};
                        _circle_inits[p].push_back(a);
                    }
                }
                else
                {
                    REMOVE_SPACE(value);
                    _setting.insert(std::make_pair(key, value));
                }
            }
            
            // Assign parameter
            setup_parameter();
        }
        catch (std::exception & e)
        {
            throw e;
        }

    }
    
    int get_int(std::string key)
    {
        auto v = get_string(key);
        return ::atoi(v.c_str());
    }
    float get_float(std::string key)
    {
        auto v = get_string(key);
        return ::atof(v.c_str());
    }
    
    std::string get_string(std::string key)
    {
        if(_setting.find(key) == _setting.end())
            throw std::runtime_error("Setting file error");
            
        return _setting[key];
    }
    
private:
    std::map<std::string, std::string> _setting;
    
    void setup_parameter();
    

};

#endif /* setting_io_hpp */
