#pragma once
#include "DSC.h"
#include <map>
#include <sstream>
#include <fstream>

namespace export_dsc {

  typedef DSC2D::DeformableSimplicialComplex dsc_type;

  bool export_fem(dsc_type * dsc, int label, std::string out_path){
    std::stringstream node_str;
    std::map<int, int> vertex_idx_map; // map[dsc_node] = fem_idx
    std::map<int, int> edge_idx; // map[half_edge] = fem_idx
    int fem_node_idx = 1;

    
//    std::stringstream obj_node, obj_face;

    std::stringstream face_str;
    int f_idx = 1;
    for(auto fit : dsc->faces())
    {
        if(dsc->get_label(fit) == label)
        {
            face_str << f_idx++;
//            obj_face << "f ";
    
            // coner nodes
            for(auto hew = dsc->walker(fit); !hew.full_circle(); hew = hew.circulate_face_ccw())
            {
                if(vertex_idx_map[hew.vertex().get_index()] == 0)
                {
                    auto pos = dsc->get_pos(hew.vertex());
                    node_str << fem_node_idx 
                             << ", " << pos[0]
                             << ", " << pos[1]
                             << ", " << 0.0 << std::endl;
                    
//                    obj_node << "v " 
//                             << " " << pos[0]
//                             << " " << pos[1]
//                             << " " << 0.0 << std::endl;
                    
                    vertex_idx_map[hew.vertex().get_index()] = fem_node_idx++;
                }
                    
                face_str << ", " << vertex_idx_map[hew.vertex().get_index()];
//                obj_face <<  " " << vertex_idx_map[hew.vertex().get_index()];
            }
    
            // edge node
            for(auto hew = dsc->walker(fit); !hew.full_circle(); hew = hew.circulate_face_ccw())
            {
                if(edge_idx[hew.halfedge().get_index()] == 0)
                {
                    auto pos = (dsc->get_pos(hew.vertex()) + dsc->get_pos(hew.next().vertex()))*0.5;
                    node_str << fem_node_idx << ", "
                             << pos[0] <<", "
                             << pos[1] <<", "
                             << 0.0 << std::endl;
    
                    edge_idx[hew.halfedge().get_index()] = fem_node_idx;
                    edge_idx[hew.opp().halfedge().get_index()] = fem_node_idx;
                    fem_node_idx++;
                }
                face_str << ", " << edge_idx[hew.halfedge().get_index()] ;
            }
    
            face_str << std::endl;
//            obj_face << std::endl;
        }
    }

    try
    {
        std::ofstream f_node(out_path + "_node.inp");
        if( f_node.is_open() )
        {
            f_node << node_str.str();
        }
        
        std::ofstream f_face(out_path + "_ele.inp");
        if( f_face.is_open() )
        {
            f_face << face_str.str();
        }
        
//        std::ofstream f_obj(out_path + ".obj");
//        if( f_obj.is_open() )
//        {
//            f_obj << obj_node.str() << obj_face.str();
//        }
    }
    catch(std::exception e)
    {
        std::cout << "Error writing to file " << out_path;
        return false;
    }

    return true;
  }
}
