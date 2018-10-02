//
//  dynamics_mul.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 3/20/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "dynamics_mul.h"
#include "helper.h"
#ifdef WIN32
#include <CGLA/Mat3x3f.h>
#else
#include <GEL/CGLA/Mat3x3f.h>
#endif

#include "define.h"
#include "util.h"
#include "polygon_f.h"
#include "adapt_mesh.h"
#include <chrono>

#include "draw.h"

vector<int> state(5,0);

using namespace std;

#define GET_FACE_IDX(fid) (int)s_dsc->get_phase_attr(fid, FACE_IDX][0]

dynamics_mul::dynamics_mul(){
    
}

dynamics_mul::~dynamics_mul(){
    
}


void dynamics_mul::write_energy()
{
    f = fopen("./LOG/MSE.txt", "w");
    
    if(f)
    {
        for (auto d : data_log)
        {
            fprintf(f, "%f %f %f \n", d[0], d[1], d[2]);
        }
        fclose(f);
    }
    
    if(s_dsc)
    cout << "DSC: " << s_dsc->get_no_faces() << " faces; " << s_dsc->get_no_vertices() << " vertices " << endl;
}

void dynamics_mul::update_dsc_with_adaptive_mesh()
{/*cout << "----------------------------" <<endl;
    thinning_triangles(true);
    s_dsc->clean_attributes();
    s_dsc->update_attributes();
    return;*/
    
    auto init_time = std::chrono::system_clock::now();
    
    int nb_displace = 5;


    // Deform the mesh
    displace_dsc();
    
    // Init force vector
    internal_node_forces = HMesh::VertexAttributeVector<Vec2>(s_dsc->get_no_nodes_allocated(), Vec2(0.0));
    external_node_forces = HMesh::VertexAttributeVector<Vec2>(s_dsc->get_no_nodes_allocated(), Vec2(0.0));
    
    // Compute the force
    compute_mean_intensity();
    compute_intensity_force();
    compute_internal_force();


    static int count = 0;

    relabel_triangles();    
    
    if (count++ % nb_displace == 0)
    {
        update_vertex_stable();
        
//        if(count < 50)
        {
            subdivide_triangles();
            cout << relabel_triangles() << "Relabel"; 
        }


        thinning_triangles();
        
        if(ADAPTIVE)
        {
//            adapt_triangle();
//            thinning();
        }

//        thinning_interface();
    }
    
    static double time = 0;
    std::chrono::duration<double> t = std::chrono::system_clock::now() - init_time;
    time += t.count();
    
    
//    double E, l;
//    get_energy(E, l);
    //    data_log.push_back(Vec3(E, l, time));
}

Vec2 dynamics_mul::get_node_displacement(Node_key vid)
{
    if(vid.get_index() > external_node_forces.size())
        return Vec2(0.0);
    
    static double n_dt = dt;
    auto dis = external_node_forces[vid] + internal_node_forces[vid]*ALPHA;
    return dis*n_dt;
}
void dynamics_mul::compute_difference()
{
    for (auto nkey : s_dsc->vertices())
    {
        double total = 0.0;
        int count_v = 0;
        for(auto hew = s_dsc->walker(nkey); !hew.full_circle(); hew = hew.circulate_vertex_cw())
        {
            auto fid = hew.face();
            
            if (fid == HMesh::InvalidFaceID) {
                continue;
            }
            
            double ci = mean_inten_[s_dsc->get_label(fid)];
            auto tris = s_dsc->get_pos(fid);
            
            Vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
            for (auto p: tris){
                min[0] = std::min(min[0], p[0]);
                min[1] = std::min(min[1], p[1]);
                max[0] = std::max(max[0], p[0]);
                max[1] = std::max(max[1], p[1]);
            }
            
            double tri_total = 0.0;
            int count = 0;
            for (int i = floor(min[0]); i < ceil(max[0]); i++) {
                for (int j = floor(min[1]); j < ceil(max[1]); j++) {
                    if (helper_t::is_point_in_tri(Vec2(i,j), tris)) {
                        double I = s_img->get_intensity(i, j);
                        tri_total += (I-ci)*(I-ci);
                        count ++;
                    }
                }
            }
            
            count_v += count;
            total += tri_total;
        } /* for(auto hew) */
        
        s_dsc->set_node_force(nkey, Vec2(total/(double)count_v), STAR_DIFFER);
    } /* for (auto nkey) */
}

void dynamics_mul::adapt_triangle()
{

}

void dynamics_mul::thinning()
{

}

void dynamics_mul::thinning_interface()
{
    for (auto vkey : s_dsc->vertices())
    {
        if (s_dsc->is_interface(vkey) && !s_dsc->is_crossing(vkey))
        {
            // interface
            std::vector<HMesh::Walker> edges;
            for (auto hew = s_dsc->walker(vkey); !hew.full_circle(); hew = hew.circulate_vertex_cw())
            {
                if (s_dsc->is_interface(hew.halfedge()))
                {
                    edges.push_back(hew);
                }
            }
            assert(edges.size()==2); 

            auto cangle = DSC2D::Util::cos_angle(
                        s_dsc->get_pos(edges[0].vertex()),
                        s_dsc->get_pos(edges[0].opp().vertex()),
                        s_dsc->get_pos(edges[1].vertex()));

            auto shortest_edge =
                    s_dsc->length(edges[0].halfedge()) <  s_dsc->length(edges[1].halfedge())?
                        edges[0] : edges[1];

            static double thres = cos(186*M_PI/180.);
            if (cangle < thres)
            {
                collapse_edge(shortest_edge.halfedge(), vkey);
//                // check quality if collapse

//                if(shortest_edge.vertex().get_index() != vkey.get_index())
//                    shortest_edge = shortest_edge.opp();

//                auto survive_vertex = shortest_edge.opp().vertex();
//                auto pos = s_dsc->get_pos(survive_vertex);

//                double min_quality = INFINITY;
//                for (auto hew = s_dsc->walker(vkey); !hew.full_circle();
//                     hew = hew.circulate_vertex_ccw())
//                {
//                    auto v1 = hew.vertex();
//                    auto v2 =  hew.next().vertex();

//                    assert(v1.get_index() != vkey.get_index() && v2.get_index() != vkey.get_index());

//                    if(v1.get_index() != survive_vertex.get_index()
//                       && v2.get_index() != survive_vertex.get_index())
//                    {
//                        min_quality = std::min(DSC2D::Util::min_angle(s_dsc->get_pos(v1), s_dsc->get_pos(v2), pos), min_quality);
//                    }
//                }

//                assert(min_quality < 100);

//                if(min_quality > M_PI * 10./180.)
//                    s_dsc->collapse(shortest_edge, 1.0);
            }
        }
    }
}

double dynamics_mul::get_energy_assume_label(Face_key fid, int assumed_label)
{
    auto pts = s_dsc->get_pos(fid);
    double external = s_img->get_sum_on_tri_differ(pts, mean_inten_[assumed_label]);
    double internal = 0;
    for(auto hew = s_dsc->walker(fid); !hew.full_circle(); hew = hew.circulate_face_ccw())
    {
        if(hew.opp().face() != HMesh::InvalidFaceID
                && s_dsc->get_label(hew.opp().face()) != assumed_label)
            internal += s_dsc->length(hew.halfedge());
    }

    return external + internal * ALPHA;
}

void dynamics_mul::update_dsc(dsc_obj &dsc, image &img){

    s_dsc = &dsc;
    s_img = &img;

    update_dsc_with_adaptive_mesh();
}

void dynamics_mul::update_vertex_stable()
{
    auto obj = s_dsc;
    
    for (auto ni = obj->vertices_begin(); ni != obj->vertices_end(); ni++)
    {
        obj->bStable[*ni] = 1; // default stable
        
        if(obj->is_crossing(*ni))
        {
            obj->bStable[*ni] = 1;
        }
        else if ((obj->is_interface(*ni)))
        {
            Vec2 dis = get_node_displacement(*ni);
            assert(dis.length() != NAN);

//            auto norm = obj->get_normal(*ni);
            double move = dis.length();// DSC2D::Util::dot(dis, norm);
            
            if (move > STABLE_MOVE) // stable
            {
                obj->bStable[*ni] = 0;
            }
        }
    }
    
}


HMesh::Walker dynamics_mul::pre_edge(dsc_obj *obj, HMesh::Walker hew){
    auto pre_e = hew.circulate_vertex_ccw();
    for (; !pre_e.full_circle(); pre_e = pre_e.circulate_vertex_ccw()) {
        if (obj->is_interface(pre_e.halfedge())) {
            pre_e = pre_e.opp();
            break;
        }
    }
    assert(obj->is_interface(pre_e.halfedge()));
    assert( pre_e.vertex() == hew.opp().vertex());
    return pre_e;
}

HMesh::Walker dynamics_mul::next_edge(dsc_obj *obj, HMesh::Walker hew){
    auto next_e = hew.opp().circulate_vertex_cw();
    for (; !next_e.full_circle(); next_e = next_e.circulate_vertex_cw()) {
        if (obj->is_interface(next_e.halfedge())) {            break;
        }
    }
    assert(obj->is_interface(next_e.halfedge()));
    assert( next_e.opp().vertex() == hew.vertex());
    
    return next_e;
}

Vec2 dynamics_mul::get_vertex_norm(dsc_obj *obj, HMesh::Walker hew){
    auto p = obj->get_pos(hew.vertex());
    auto p_opp = obj->get_pos(hew.opp().vertex());
    
    Vec2 N;
    auto next_e = next_edge(obj, hew);
    
    Vec2 Ne = p - p_opp;
    Ne.normalize();
    Ne = Vec2(Ne[1], -Ne[0]);
    
    Vec2 Ne_next = obj->get_pos(next_e.vertex()) - p;
    Ne_next.normalize();
    Ne_next = Vec2(Ne_next[1], -Ne_next[0]);
    
    N = Ne + Ne_next;
    N.normalize();
    return N;
}


int dynamics_mul::relabel_triangles()
{
    int count = 0;
    for(auto fid : s_dsc->faces())
    {
        auto current_label = s_dsc->get_label(fid);
        
        int min_label = -1;
        double min_E = INFINITY;
        for(int ll = 0; ll < mean_inten_.size(); ll++)
        {
            double E = get_energy_assume_label(fid, ll);
            if(min_E > E)
            {
                min_E = E;
                min_label = ll;
            }
        }
        
        if(current_label != min_label)
        {
            s_dsc->set_label(fid, min_label);
            count ++;
        }
    }
    
    return count;
}

std::map<int,double> dynamics_mul::get_energy_thres()
{
    std::map<int,double> e_thres;
    for(auto phase = mean_inten_.begin(); phase != mean_inten_.end(); phase ++)
    {
        auto c = phase->first;
        double mean_differ = INFINITY;
        
        for(int ll = 1; ll < mean_inten_.size(); ll++)
        {
            if(ll != phase->first)
            {
                mean_differ = std::min(mean_differ, 
                  (mean_inten_[ll] - phase->second)*(mean_inten_[ll] - phase->second)
                                       );
            }
        }
        
        e_thres.insert(std::make_pair(phase->first, mean_differ));
    }
    
    double mean_area = SMALLEST_SIZE*SMALLEST_SIZE*0.5;
    double mean_circum = 3*SMALLEST_SIZE;
    for(auto & et : e_thres )
    {
        et.second = et.second * mean_area + mean_circum*ALPHA;
    }
    
    return e_thres;
}

void dynamics_mul::subdivide_triangles()
{
    auto thres_hold = get_energy_thres();
    double min_area = SMALLEST_SIZE*SMALLEST_SIZE*0.5;
    
    HMesh::AttributeVector<int, Face_key> faces_to_split(s_dsc->get_no_faces_allocated(), 0);
    int count = 0;
    for(auto fid : s_dsc->faces())
    {
        if(s_dsc->area(fid) < min_area)
            continue;
        
        auto l = s_dsc->get_label(fid);
        
        // Compute the variance
        auto tris = s_dsc->get_pos(fid);
        auto mean_c = s_img->get_sum_on_tri_intensity(tris) / s_dsc->area(fid);
        auto ext_E = s_img->get_sum_on_tri_differ(tris, mean_c);
        
        if(ext_E > thres_hold[l]) // Consider splitting
        {
            count ++;
            faces_to_split[fid] = 1;
        }
    }
    
    std::cout << count << " faces need to be subdivided" <<endl;
    
    s_dsc->recursive_split(faces_to_split);
}

void dynamics_mul::thinning_triangles(bool colapse_all)
{    
    auto thres_hold = get_energy_thres();
    double min_area = SMALLEST_SIZE*SMALLEST_SIZE*0.5;
    
    // Mark the triangles that are homeogeneous
    HMesh::AttributeVector<int, Face_key> face_can_collapse(s_dsc->get_no_faces_allocated(), 1);
    
    int count = 0;
    
    if(!colapse_all)
    {
        for(auto fid : s_dsc->faces())
        {        
            auto l = s_dsc->get_label(fid);
            
            // Compute the variance
            auto tris = s_dsc->get_pos(fid);
            auto mean_c = s_img->get_sum_on_tri_intensity(tris) / s_dsc->area(fid);
            auto ext_E = s_img->get_sum_on_tri_differ(tris, mean_c);
            
            if(ext_E > thres_hold[l]) // Consider splitting
            {
                face_can_collapse[fid] = 0;
                count ++;
            }
        }
    }
    

    // Collapse interior vertices
    int count_internal = 0, counter_interface = 0, count_boundary = 0;
    int cc = 0;
    for(auto vid : s_dsc->vertices())
    {

        bool can_collaspe = true;
        for(auto hew = s_dsc->walker(vid); !hew.full_circle(); hew = hew.circulate_face_ccw())
        {
            if( hew.face().get_index() <  face_can_collapse.size() 
                    &&  face_can_collapse[hew.face()] == 0)
            {
                can_collaspe = false;
                break;
            }
        }
        
        if(can_collaspe)
        {
           cc++;
           
           if(s_dsc->is_interface(vid))
           {
               if(!s_dsc->is_crossing(vid) 
                  && !HMesh::boundary(*s_dsc->mesh, vid))
               {
                   // only collapse flat edges
                   std::vector<HMesh::Walker> edges;
                   for (auto hew = s_dsc->walker(vid); !hew.full_circle(); hew = hew.circulate_vertex_cw())
                   {
                       if (s_dsc->is_interface(hew.halfedge()))
                       {
                           edges.push_back(hew);
                       }
                   }
                   
                   if(edges.size()!=2){
                       continue;
                   }
                   
                   // Check if it is stable
//                   if(s_dsc->bStable[edges[0].vertex()] == 0
//                           || s_dsc->bStable[edges[0].opp().vertex()] == 0
//                           || s_dsc->bStable[edges[1].vertex()] == 0)
//                   {
//                       continue;
//                   }
                   
                   auto cangle = DSC2D::Util::cos_angle(
                               s_dsc->get_pos(edges[0].vertex()),
                               s_dsc->get_pos(edges[0].opp().vertex()),
                               s_dsc->get_pos(edges[1].vertex()));
                   
                   auto shortest_edge =
                           s_dsc->length(edges[0].halfedge()) <  s_dsc->length(edges[1].halfedge())?
                               edges[0] : edges[1];     
                   
                   static double thres = cos(175*M_PI/180.);
                   
//                   s_dsc->collapse(shortest_edge, true);
                           
                   if(cangle < thres &&
                           collapse_edge(shortest_edge.halfedge(), vid, false))
//                           s_dsc->collapse(shortest_edge, false))
                       counter_interface++;
               }
           }
           else
           {
               if(!HMesh::boundary(*s_dsc->mesh, vid))
               {
                    // Collapse the shortest edge
                    // Not the optimal choice, but we priorize performance
                    double shortest_length = INFINITY;
                    Edge_key shortest_edge;
                    for(auto hew = s_dsc->walker(vid); !hew.full_circle(); hew = hew.circulate_vertex_ccw())
                    {
                        auto length = s_dsc->length(hew.halfedge());
                        if(length < shortest_length)
                        {
                            shortest_length = length;
                            shortest_edge = hew.halfedge();
                        }
                    }
                    
        
                    
                    if(s_dsc->collapse(shortest_edge, true))
                        count_internal++;
               }
               else
               {
                   // Collapse the shortest edge
                   std::vector<HMesh::Walker> edges;
                   for (auto hew = s_dsc->walker(vid); !hew.full_circle(); hew = hew.circulate_vertex_cw())
                   {
                       if (HMesh::boundary(*s_dsc->mesh, hew.halfedge()))
                       {
                           edges.push_back(hew);
                       }
                   }
                   
                   if(edges.size()!=2){
                       continue;
                   }
                   
                   auto cangle = DSC2D::Util::cos_angle(
                               s_dsc->get_pos(edges[0].vertex()),
                               s_dsc->get_pos(edges[0].opp().vertex()),
                               s_dsc->get_pos(edges[1].vertex()));
                   
                   auto shortest_edge =
                           s_dsc->length(edges[0].halfedge()) <  s_dsc->length(edges[1].halfedge())?
                               edges[0] : edges[1];  
                   static double thres = cos(174*M_PI/180.);

                   if(cangle < thres 
                      && collapse_edge(shortest_edge.halfedge(), vid, false))
                       count_boundary++;
               }
           }
        }
    }
    
    cout << "thinning (interial/interface/boundary) " << count_internal << "/"  << counter_interface << "/" << count_boundary << endl;
    
    s_dsc->update_attributes();
    
}

bool dynamics_mul::collapse_edge(Edge_key ek, Node_key n_to_remove, bool safe)
{
    
    auto hew = s_dsc->walker(ek);
    if(hew.vertex().get_index() == n_to_remove.get_index())
        hew = hew.opp();
    
    if (!precond_collapse_edge(*s_dsc->mesh, hew.halfedge()) )
    {
        return false;
    }    
    
    if(!s_dsc->is_collapsable(hew, safe))
    {
        return false;
    }
    
    // Find neighbors
    std::vector<Face_key> e_fids = {hew.face(), hew.opp().face()};
    std::vector<Edge_key> eids0;
    
    for(auto hw = s_dsc->walker(n_to_remove); !hw.full_circle(); hw = hw.circulate_vertex_cw())
    {
        if(hw.face() != e_fids[0] && hw.face() != e_fids[1])
        {
            eids0.push_back(hw.next().halfedge());
        }
    }

    
    // Check the quality after collapsing
    Vec2 p_new = s_dsc->get_pos(hew.vertex());
    double q = s_dsc->min_quality(eids0, s_dsc->get_pos(n_to_remove), p_new);

    static double min_angle = M_PI * 20. / 180.;
    if(q > min_angle)
    {
        auto result = s_dsc->collapse(hew, 0, safe);
//        if(result)
//        {
//
//            cout << "Collapsed " << n_to_remove.get_index() << endl;
//        }
            
        return result;
    }
    
    return false;
}

double dynamics_mul::energy_triangle(HMesh::FaceID fid, double c,  int new_phase){
    auto tris = s_dsc->get_pos(fid);
    
    double ET = s_img->get_sum_on_tri_differ(tris, c);
    
    // length
    auto hew = s_dsc->walker(fid);
    double length = 0.;
    for (int i = 0; i < 3; i++) {
        if(HMesh::boundary(*s_dsc->mesh, hew.halfedge()))
            continue;
        
        int phase_opp = s_dsc->get_label( hew.opp().face() );
        if(phase_opp != new_phase){
            length += s_dsc->length(hew.halfedge());
        }
        hew = hew.next();
    }
    
    return ET + length * ALPHA;
}

#define NB_PHASE 3

void dynamics_mul::update_probability(dsc_obj &dsc, image &img){
    s_img = &img;
    s_dsc = &dsc;
    
    // Init
    static std::vector<XI> XI_f;
    static KPhi k_f(NB_PHASE);
    
    static bool inited = false;
    if (!inited) {
        //
        for (auto fkey : s_dsc->faces()) {
            s_dsc->set_phase_attr(fkey, Vec2(s_dsc->get_label(fkey)), PHASE_PROBABILITY);
        }
        
        for (int i = 0; i < NB_PHASE; i++) {
            XI_f.push_back(XI(i, NB_PHASE));
        }
        
        inited = true;
    }
    
    // Evolution
    // Compute intensity
    
}


//void dynamics_mul::indexing_vertices()
//{
//    int idx = 0;
//    for(auto vid = s_dsc->vertices_begin(); vid != s_dsc->vertices_end(); vid++){
//        if(!HMesh::boundary(*s_dsc->mesh, *vid)){
//            s_dsc->set_node_force(*vid, Vec2(idx), INDEX_VERT);
//            idx++;
//        }
//    }
//    for(auto vid = s_dsc->vertices_begin(); vid != s_dsc->vertices_end(); vid++){
//        if(HMesh::boundary(*s_dsc->mesh, *vid)){
//            s_dsc->set_node_force(*vid, Vec2(idx), INDEX_VERT);
//            idx++;
//        }
//    }
//    assert(idx == s_dsc->get_no_vertices());
    
//    idx = 0;
//    for (auto fid = s_dsc->faces_begin(); fid != s_dsc->faces_end(); fid++) {
//        s_dsc->set_phase_attr(*fid, Vec2(idx++), FACE_IDX);
//    }
//}

//std::vector<int> dynamics_mul::get_vert_idx(std::vector<HMesh::VertexID> vids){
//    std::vector<int> idxs;
//    for (auto v : vids) {
//        idxs.push_back((int)s_dsc->get_node_force(v, INDEX_VERT)[0]);
//    }
    
//    return idxs;
//}

//void dynamics_mul::build_and_solve(){


//}


void dynamics_mul::get_energy(double &e_, double &l_)
{
    // intensity difference
    double E = 0;
    for (auto fkey : s_dsc->faces())
    {
        auto pts = s_dsc->get_pos(fkey);
        double ci = mean_inten_[s_dsc->get_label(fkey)];
        E += s_img->get_tri_differ_f(pts, ci);
    }
    
    // Length
    double L = 0;
    for (auto ekey : s_dsc->halfedges())
    {
        auto hew = s_dsc->walker(ekey);
        if (hew.vertex().get_index() > hew.opp().vertex().get_index()
            && s_dsc->is_interface(hew.halfedge()))
        {
            L += s_dsc->length(ekey);
        }
    }
    
    e_ = E;
    l_ = L;
}

//double dynamics_mul::furthest_move(Node_key nid, Vec2 direction){
//    double max_move = s_dsc->intersection_with_link(nid, s_dsc->get_pos(nid) + direction);
//    return max_move;
//}


void dynamics_mul::compute_internal_force(){
    
    for (auto vkey : s_dsc->vertices())
    {
        if (s_dsc->is_interface(vkey)
            or s_dsc->is_crossing(vkey))
        {
            for(auto w = s_dsc->walker(vkey); !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                if (s_dsc->is_interface(w.halfedge()))
                {
                    auto p12 = s_dsc->get_pos(w.vertex()) - s_dsc->get_pos(w.opp().vertex());

                    if(p12.length() > 0.001)
                    {
                        p12.normalize();
                        internal_node_forces[vkey] += p12;
                    }
                    else
                        cout << "Small edge length" <<endl;
                }
            }
        }
    }
}


inline bool is_bound(DSC2D::DeformableSimplicialComplex * dsc, HMesh::VertexID n)
{
    
    for (auto hew = dsc->walker(n); !hew.full_circle(); hew = hew.circulate_vertex_ccw())
    {
        if (dsc->mesh->in_use(hew.face()) && dsc->get_label(hew.face()) == BOUND_FACE)
        {
            return true;
        }
    }
    return false;
}

inline void crop(Vec2 & p, Vec2 & bound)
{
    for(int i = 0; i < 2; i++)
    {
        if(p[i] < 0)
            p[i] = 0;
        if(p[i] > bound[i])
            p[i] = bound[i];
    }
}

inline void proximity_snap(Vec2 const & origin, Vec2 & des, Vec2 const & bound)
{
    static double eps = SMALLEST_SIZE * 0.1;
    for(int i = 0; i < 2; i++)
    {
        if(origin[i] < eps )
            des[i] = 0;
        if( bound[i] - origin[i] < eps)
            des[i] = bound[i];
    }
}

void dynamics_mul::displace_dsc(dsc_obj *obj){
    if (!obj) {
        obj = s_dsc;
    }

    auto bound = s_img->size();

    for (auto ni = obj->vertices_begin(); ni != obj->vertices_end(); ni++)
    {
        obj->bStable[*ni] = 1;
        
        if ((obj->is_interface(*ni) or obj->is_crossing(*ni)))
        {
            Vec2 dis = get_node_displacement(*ni);
            assert(dis.length() != NAN);
            
            // Crop destination to avoid boundary corruption
            auto origin = obj->get_pos(*ni);
            auto des = origin + dis;
            
//            crop(des, bound);
//            proximity_snap(origin, des, bound);

            obj->set_destination(*ni, des);
        }
    }
    
    obj->deform(ADAPTIVE);
}

void dynamics_mul::compute_mean_intensity()
{
    std::map<int, double> mean_inten_o;
    
    std::map<int, double> total_area;
    
    mean_inten_o.clear();
    
    
    for (auto fid : s_dsc->faces()) {

        double area = 0.0;
        
        auto tris = s_dsc->get_pos(fid);
        double total_inten = s_img->get_tri_intensity_f(tris, &area);
        
        
        int phase = s_dsc->get_label(fid);
        if (mean_inten_o.find(phase) != mean_inten_o.end())
        {//Existed
            mean_inten_o[phase] += total_inten;
            total_area[phase] += area;
        }
        else
        {
            total_area.insert(std::make_pair(phase, area));
            mean_inten_o.insert(std::make_pair(phase, total_inten));
        }
    }
    
    
    for (auto mit = mean_inten_o.begin(); mit != mean_inten_o.end(); mit++)
    {
        mit->second /= (double)total_area[mit->first];
    }
    
//    mean_inten_o[BOUND_FACE] = -0.2;

    cout << "Mean intensity: ";
    for(auto m : mean_inten_o) cout << m.second << " ";
    cout << endl;

    g_param.mean_intensity = mean_inten_o;
    mean_inten_ = mean_inten_o;
}

void dynamics_mul::compute_intensity_force(){
    
    for (auto eit = s_dsc->halfedges_begin(); eit != s_dsc->halfedges_end(); eit++) {
        auto hew = s_dsc->walker(*eit);
        
        if( s_dsc->is_interface(*eit) and
                hew.halfedge() < hew.opp().halfedge()) // Avoid retouch the edge
        {
            int phase0 = s_dsc->get_label(hew.face());
            int phase1 = s_dsc->get_label(hew.opp().face());
            
            double c0 = mean_inten_[phase0];
            double c1 = mean_inten_[phase1];
            
            // Loop on the edge
            auto p0 = s_dsc->get_pos(hew.opp().vertex());
            auto p1 = s_dsc->get_pos(hew.vertex());
            
            double length = (p1 - p0).length();
            
            double f0 = 0.0, f1 = 0.0;
            
            // Integrate on the edge
            int N = ceil(length);

            double dl = (double)length / (double)N;
            for (int i = 0; i < N; i++) {
                auto p = p0 + (p1 - p0)*((double)i / (double)N);
                double I = s_img->get_intensity_f(p[0], p[1]);
                
                double f = 0.0;
                // Same coefficient
//                f = (2*I - c0 - c1) / (c0-c1) * dl;
                
                // No normalization
                f = (2*I - c0 - c1) * (c0-c1) * dl;
                
                assert(!DSC2D::Util::isnan(f));
                
                // Barry Centric coordinate
                f0 += f* (p-p1).length() / length;
                f1 += f* (p-p0).length() / length;
            }
            
            // Normal vector from phase0 to phase1
            Vec2 L01 = p1 - p0;
            Vec2 N01(L01[1], -L01[0]);
//            auto other_n = s_dsc->get_pos(hew.next().vertex());
//            N01 = N01 * DSC2D::Util::dot(N01, other_n - p0);
            N01.normalize();

            
            Vec2 f_x0 = N01*f0;
            Vec2 f_x1 = N01*f1;
            
            assert(f0 != NAN);
            assert(f1 != NAN);
            
            external_node_forces[hew.opp().vertex()] += f_x0;
            external_node_forces[hew.vertex()] += f_x1;
        }
    }
    
}



