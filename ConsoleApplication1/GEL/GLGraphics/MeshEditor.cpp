//
//  MeshEditor.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 09/10/13.
//
//
#include <thread>
#if !defined (WIN32)
#include <stdarg.h>
#include <unistd.h>
#endif
#include <GL/glew.h>
#include <functional>
#include "MeshEditor.h"
#include <string>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

#include <GLGraphics/Console.h>
#include <GLGraphics/glsl_shader.h>
#include <GLGraphics/ShadowBuffer.h>

#include <CGLA/CGLA.h>

#include <HMesh/HMesh.h>

#include <Util/Timer.h>

#include "VisObj.h"

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Util;

namespace GLGraphics {
    
    namespace {
        
        bool wantshelp(const std::vector<std::string> & args)
        {
            if(args.size() == 0)
                return false;
            
            string str = args[0];
            
            if(str=="help" || str=="HELP" || str=="Help" || str=="?")
                return true;
            
            return false;
        }

// Commented out since apparently unused:
//        bool wantshelp(MeshEditor* me, const std::vector<std::string> & args)
//        {
//            if(args.size() == 0)
//                return false;
//            
//            string str = args[0];
//            
//            if(str=="help" || str=="HELP" || str=="Help" || str=="?")
//                return true;
//            
//            return false;
//        }
        
        /// Function that aligns two meshes.
        void console_align(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: align <dest> <src>");
                me->printf("This function aligns dest mesh with src");
                me->printf("In practice the GLViewController of src is copied to dst.");
                me->printf("both arguments are mandatory and must be numbers between 1 and 9.");
                me->printf("Note that results might be unexpexted if the meshes are not on the same scale");
            }
            
            int dest = 0;
            
            if(args.size()>0){
                istringstream a0(args[0]);
                a0 >> dest;
                --dest;
                
                if(dest <0 || dest>8)
                {
                    me->printf("dest mesh out of range (1-9)");
                    return;
                }
            }
            else
            {
                me->printf("neither source nor destination mesh?!");
                return;
            }
            
            int src = 0;
            if(args.size()>1){
                istringstream a1(args[1]);
                a1 >> src;
                --src;
                
                if(src <0 || src>8)
                {
                    me->printf("src mesh out of range (1-9)");
                    return;
                }
            }
            else
            {
                me->printf("no src mesh?");
                return;
            }
            me->align(src,dest);
        }
        
        
        
        
        void console_ridge_lines(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: ridge_lines");
                return;
            }
            
            me->save_active_mesh();
            
            Manifold& mani = me->active_mesh();
            
            VertexAttributeVector<Mat3x3d> curvature_tensors(mani.allocated_vertices());
            VertexAttributeVector<Vec3d> min_curv_direction(mani.allocated_vertices());
            VertexAttributeVector<Vec3d> max_curv_direction(mani.allocated_vertices());
            VertexAttributeVector<Vec2d> curvature(mani.allocated_vertices());
            
            
            
            curvature_paraboloids(mani,
                                  min_curv_direction,
                                  max_curv_direction,
                                  curvature);
            
            for(auto vid : mani.vertices())
            {
                Vec3d max_curv_dir = normalize(max_curv_direction[vid]);
                Vec3d min_curv_dir = normalize(min_curv_direction[vid]);
                double vid_min_pc = curvature[vid][0];
                double vid_max_pc = curvature[vid][1];
                bool ridge = true;
                bool ravine = true;
                Walker w = mani.walker(vid);
                Vec3d r(0);
                for(; !w.full_circle();w = w.circulate_vertex_ccw())
                {
                    Vec3d e = (mani.pos(w.vertex()) - mani.pos(vid));
                    
                    if(abs(dot(min_curv_dir,e)) > abs(dot(max_curv_dir,e)))
                    {
                        if(curvature[w.vertex()][0]<vid_min_pc+20)
                            ravine = false;
                        
                    }
                    else
                    {
                        if(curvature[w.vertex()][1]>vid_max_pc-20)
                            ridge = false;
                    }
                }
                DebugRenderer::vertex_colors[vid] = Vec3f(ridge,ravine,0.0);
            }
            for(auto fid : mani.faces())
                DebugRenderer::face_colors[fid] = Vec3f(.3,.3,.6);
            for(auto hid : mani.halfedges()) {
                
                Walker w = mani.walker(hid);
                Vec3f c0 = DebugRenderer::vertex_colors[w.opp().vertex()];
                Vec3f c1 = DebugRenderer::vertex_colors[w.vertex()];
                
                DebugRenderer::edge_colors[hid] = (c0==c1) ? c0 : Vec3f(0.1,0.1,0.3);
                
            }
        }

        void transform_mesh(Manifold& mani, const Mat4x4d& m)
        {
            for(VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
                mani.pos(*vid) = m.mul_3D_point(mani.pos(*vid));
        }
        
        void console_scale(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.scale sx sy sz");
                me->printf("Note: If only sx is provided, uniform scaling is applied");
                return;
            }
            
            Vec3d s;
            if(args.size() > 0){
                istringstream a0(args[0]);
                double scale;
                a0 >> scale;
                s = Vec3d(scale);
            }
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> s[1];
            }
            if(args.size() > 2){
                istringstream a0(args[2]);
                a0 >> s[2];
            }
            
            me->save_active_mesh();
            transform_mesh(me->active_mesh(),scaling_Mat4x4d(s));
        }

        void console_rotate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.rotate axis_x axis_y axis_z angle");
                return;
            }
            
            Vec3d a;
            double angle;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> a[0];
            }
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> a[1];
            }
            if(args.size() > 2){
                istringstream a0(args[2]);
                a0 >> a[2];
            }
            if(args.size() > 3){
                istringstream a0(args[3]);
                a0 >> angle;
            }
            
            me->save_active_mesh();
            Quatd q;
            q.make_rot(angle, a);
            transform_mesh(me->active_mesh(),q.get_Mat4x4d());
        }

        
        void console_translate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: transform.translate tx ty tz");
                me->printf("Note: recenters if no arguments are provided.");
                return;
            }
            
            Vec3d t;
            
            if(args.size()==0)
            {
                float rad;
                bsphere(me->active_mesh(), t, rad);
            }
            else {
                if(args.size() > 0){
                    istringstream a0(args[0]);
                    a0 >> t[0];
                }
                if(args.size() > 1){
                    istringstream a0(args[0]);
                    a0 >> t[1];
                }
                if(args.size() > 2){
                    istringstream a0(args[0]);
                    a0 >> t[2];
                }
            }
            
            me->save_active_mesh();
            transform_mesh(me->active_mesh(),translation_Mat4x4d(-t));
        }

        
        void console_merge_1_ring(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.merge_1_ring");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_vertex_selection();
            for(auto v: sel)
                if(m.in_use(v))
                    m.merge_one_ring(v);
        }
        
        
        void console_flip_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.flip_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h) && precond_flip_edge(m, h))
                    m.flip_edge(h);
        }

        void console_collapse_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.collapse_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h) && precond_collapse_edge(m, h))
                    m.collapse_edge(h,true);
        }

        void console_dissolve_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.dissolve_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h))
                    m.merge_faces(m.walker(h).face(), h);
        }
        
        void console_split_edge(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.split_edge");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_halfedge_selection();
            for(auto h: sel)
                if(m.in_use(h))
                    m.split_edge(h);
        }
        
        void console_stellate_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.stellate_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            for(auto f: sel)
                if(m.in_use(f))
                    m.split_face_by_vertex(f);
        }
        
        
        void console_triangulate_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.triangulate_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            for(auto f: sel)
                if(m.in_use(f))
                    triangulate_face_by_edge_split(m, f);
        }
        
        void console_bridge_faces(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.bridge_faces");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto sel = me->get_face_selection();
            if(sel.size() != 2) {
                me->printf("You must select exactly two faces");
                return;
            }

            int i=0;
            FaceID f0, f1;
            for(auto f: sel)
            {
                int n = no_edges(m, f);
                if(i==0) {
                    f0 = f;
                    i+= n;
                }
                else {
                    i-= n;
                    f1 = f;
                }
            }
            if(i!= 0){
                me->printf("The selected faces must have same number of edges");
                return;
            }
            vector<VertexID> loop0;
            circulate_face_ccw(m, f0, [&](VertexID v) {
                loop0.push_back(v);
            });
            
            vector<VertexID> loop1;
            circulate_face_ccw(m, f1, [&](VertexID v) {
                loop1.push_back(v);
            });
            
            vector<pair<VertexID, VertexID> > connections;
            
            size_t L0= loop0.size();
            size_t L1= loop1.size();
            
            assert(L0==L1);
            
            size_t L = L0;
            
            float min_len = FLT_MAX;
            int j_off_min_len = -1;
            for(int j_off = 0; j_off < L; ++j_off)
            {
                float len = 0;
                for(int i=0;i<L;++i)
                    len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
                if(len < min_len)
                {
                    j_off_min_len = j_off;
                    min_len = len;
                }
            }
            for(int i=0;i<L;++i)
                connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
            // Merge the two one rings producing two faces.
            
            // Bridge the just created faces.
            vector<HalfEdgeID> newhalfedges = m.bridge_faces(f0, f1, connections);
        }


        void console_delete(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.merge_1_ring");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto vsel = me->get_vertex_selection();
            auto hsel = me->get_halfedge_selection();
            auto fsel = me->get_face_selection();
            for(auto v: vsel)
                if(m.in_use(v))
                    m.remove_vertex(v);
            for(auto h: hsel)
                if(m.in_use(h))
                    m.remove_edge(h);

            for(auto f: fsel)
                if(m.in_use(f))
                    m.remove_face(f);
        }
        
        void console_split_face(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: edit.selected.split_face");
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            auto vsel = me->get_vertex_selection();
            auto fsel = me->get_face_selection();
            for(auto f: fsel)
                if(m.in_use(f))
                {
                    VertexID v0=InvalidVertexID,v1=InvalidVertexID;
                    circulate_face_ccw(m, f, [&](VertexID v){
                        if(vsel.count(v))
                        {
                            if(v0==InvalidVertexID)
                                v0 = v;
                            else
                                v1 = v;
                        
                        }
                    
                    });
                    m.split_face_by_edge(f, v0, v1);
                }
        }
        


        
        
        void console_refit_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("display.refit_trackball");
                return;
            }
            me->refit();
        }

        
        void console_test(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: test");
                return;
            }
            
            me->save_active_mesh();
            me->active_mesh().slit_edges(me->get_vertex_selection());
        }
        
        void console_save(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: save <name->x3d|name->obj> ");
                
                return;
            }
            const string& file_name = args[0];
            if(args.size() == 1){
                if(file_name.substr(file_name.length()-4,file_name.length())==".obj"){
                    obj_save(file_name, me->active_mesh());
                    
                    return;
                }
                else if(file_name.substr(file_name.length()-4,file_name.length())==".off"){
                    off_save(file_name, me->active_mesh());
                    
                    return;
                }
                else if(file_name.substr(file_name.length()-4,file_name.length())==".x3d"){
                    x3d_save(file_name, me->active_mesh());
                    
                    return;
                }
                me->printf("unknown format");
                return;
            }
            me->printf("usage: save <name->x3d|name->obj> ");
        }
        
        
        void console_refine_edges(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: refine.split_edges <length>");
                me->printf("splits edges longer than <length>; default is 0.5 times average length");
                return;
            }
            
            me->save_active_mesh();
            
            float thresh = 0.5f;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> thresh;
            }
            
            float avg_length = average_edge_length(me->active_mesh());
            
            refine_edges(me->active_mesh(), thresh * avg_length);
            
            return;
            
        }
        
        void console_refine_faces(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: refine.split_faces ");
                me->printf("usage:  Takes no arguments. Inserts a vertex at the centre of each face.");
                
                return;
            }
            me->save_active_mesh();
            
            triangulate_by_vertex_face_split(me->active_mesh());
            
            return;
            
        }
        
        void console_cc_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.catmull_clark ");
                me->printf("Does one step of Catmull-Clark subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            cc_split(me->active_mesh(),me->active_mesh());
            cc_smooth(me->active_mesh());
            
            return;
        }
        
        void console_root_cc_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.catmull_clark ");
                me->printf("Does one step of Catmull-Clark subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            rootCC_subdivide(me->active_mesh(),me->active_mesh());
            return;
        }

        
        void console_loop_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.loop");
                me->printf("Does one step of Loop subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            loop_split(me->active_mesh(),me->active_mesh());
            loop_smooth(me->active_mesh());
            
            return;
        }
        
        void console_stitch(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.stitch <rad>");
                me->printf("Stitches faces");
                
                return;
            }
            double r = 0.001;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> r;
            }
            
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            Vec3d c;
            float rad;
            bsphere(m, c, rad);
            stitch_mesh(me->active_mesh(), r * rad);
            return;
        }
        
        void console_remove_duplicates(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.remove_duplicates <rad>");
                me->printf("Removes duplicate vertices and incident faces");
                
                return;
            }
            double r = 0.001;
            
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> r;
            }
            
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            Vec3d c;
            float rad;
            bsphere(m, c, rad);
            remove_duplicates(me->active_mesh(), r * rad);
            return;
        }
        
        void console_compact(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.compact");
                me->printf("Removes unreferenced vertices");
                
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            m.cleanup();
            return;
        }


        
        void console_remove_val2(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: cleanup.remove_val2");
                me->printf("Removes valence 2 vertices");
                
                return;
            }
            me->save_active_mesh();
            Manifold& m = me->active_mesh();
            remove_valence_two_vertices(m);
            return;
        }

        
        void console_root3_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.root3");
                me->printf("Does one step of sqrt(3) subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            root3_subdivide(me->active_mesh(),me->active_mesh());
            
            return;
        }
        
        
        void console_doosabin_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.doo_sabin ");
                me->printf("Does one step of Doo-Sabin Subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            cc_split(me->active_mesh(),me->active_mesh());
            dual(me->active_mesh());
            
            return;
        }
        
        void console_butterfly_subdivide(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: subdivide.butterfly ");
                me->printf("Does one step of Modified Butterfly Subdivision");
                
                return;
            }
            me->save_active_mesh();
            
            butterfly_subdivide(me->active_mesh(),me->active_mesh());
            
            return;
        }
        
        void console_dual(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: dual ");
                me->printf("Produces the dual by converting each face to a vertex placed at the barycenter.");
                return;
            }
            me->save_active_mesh();
            
            dual(me->active_mesh());
            
            return;
        }
        
        
        void console_minimize_curvature(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.minimize_curvature <anneal>");
                me->printf("Flip edges to minimize mean curvature.");
                me->printf("If anneal is true, simulated annealing (slow) is used rather than a greedy scheme");
                return;
            }
            me->save_active_mesh();
            
            bool anneal=false;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> anneal;
            }
            
            minimize_curvature(me->active_mesh(), anneal);
            me->post_create_display_list();
            return;
        }
        
        void console_minimize_dihedral(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.minimize_dihedral <iter> <anneal> <use_alpha> <gamma> ");
                me->printf("Flip edges to minimize dihedral angles.");
                me->printf("Iter is the max number of iterations. anneal tells us whether to use ");
                me->printf("simulated annealing and not greedy optimization. use_alpha (default=true) ");
                me->printf("means to use angle and not cosine of anglegamma (default=4) is the power ");
                me->printf("to which we raise the dihedral angle");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1000;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> iter;
            }
            
            bool anneal = false;
            if(args.size() > 1)
            {
                istringstream a0(args[1]);
                a0 >> anneal;
            }
            
            bool use_alpha = true;
            if(args.size() > 2)
            {
                istringstream a0(args[2]);
                a0 >> use_alpha;
            }
            
            float gamma = 4.0f;
            if(args.size() > 3)
            {
                istringstream a0(args[3]);
                a0 >> gamma;
            }
            
            
            minimize_dihedral_angle(me->active_mesh(), iter, anneal, use_alpha, gamma);
            return;
        }
        
        void console_maximize_min_angle(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.maximize_min_angle <thresh> <anneal>");
                me->printf("Flip edges to maximize min angle - to make mesh more Delaunay.");
                me->printf("If the dot product of the normals between adjacent faces < thresh");
                me->printf("no flip will be made. anneal selects simulated annealing rather ");
                me->printf("nthan greedy optimization.");
                return;
            }
            me->save_active_mesh();
            
            float thresh = 0.0f;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> thresh;
            }
            bool anneal = false;
            if(args.size() > 1)
            {
                istringstream a0(args[1]);
                a0 >> anneal;
            }
            maximize_min_angle(me->active_mesh(),thresh,anneal);
            return;
        }
        
        
        void console_optimize_valency(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: optimize.valency <anneal> ");
                me->printf("Optimizes valency for triangle meshes. Anneal selects simulated annealing rather than greedy optim.");
                return;
            }
            me->save_active_mesh();
            
            bool anneal = false;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> anneal;
            }
            optimize_valency(me->active_mesh(), anneal);
            return;
        }
        
//            void console_analyze(MeshEditor* me, const std::vector<std::string> & args)
//            {
//                if(wantshelp(args))
//                {
//                    me->printf("usage:  harmonics.analyze");
//                    me->printf("Creates the Laplace Beltrami operator for the mesh and finds all eigensolutions.");
//                    me->printf("It also projects the vertices onto the eigenvectors - thus transforming the mesh");
//                    me->printf("to this basis.");
//                    me->printf("Note that this will stall the computer for a large mesh - as long as we use Lapack.");
//                    return;
//                }
//                me->harmonics_analyze_mesh();
//                return;
//            }
//        
//        
//        void console_partial_reconstruct(MeshEditor* me, const std::vector<std::string> & args)
//        {
//            if(args.size() != 3)
//                me->printf("usage: haramonics.partial_reconstruct <e0> <e1> <s>");
//            
//            if(wantshelp(args)) {
//                me->printf("Reconstruct from projections onto eigenvectors. The two first arguments indicate");
//                me->printf("the eigenvector interval that we reconstruct from. The last argument is the ");
//                me->printf("scaling factor. Thus, for a vertex, v, the formula for computing the position, p, is:");
//                me->printf("for (i=e0; i<=e1;++i) p += proj[i] * Q[i][v] * s;");
//                me->printf("where proj[i] is the 3D vector containing the x, y, and z projections of the mesh onto");
//                me->printf("eigenvector i. Q[i][v] is the v'th coordinate of the i'th eigenvector.");
//                me->printf("Note that if vertex coordinates are not first reset, the result is probably unexpected.");
//            }
//            me->save_active_mesh();
//            
//            if(args.size() != 3)
//                return;
//            
//            int E0,E1;
//            float scale;
//            istringstream a0(args[0]);
//            a0 >> E0;
//            istringstream a1(args[1]);
//            a1 >> E1;
//            istringstream a2(args[2]);
//            a2 >> scale;
//            me->harmonics_partial_reconstruct(E0,E1,scale);
//            return;
//        }
//        
//        void console_reset_shape(MeshEditor* me, const std::vector<std::string> & args)
//        {
//            if(wantshelp(args))
//            {
//                me->printf("usage: harmonics.reset_shape ");
//                me->printf("Simply sets all vertices to 0,0,0. Call this before doing partial_reconstruct");
//                me->printf("unless you know what you are doing.");
//                return;
//            }
//            me->save_active_mesh();
//            me->harmonics_reset_shape();
//            return;
//        }
        
        
        void console_close_holes(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: cleanup.close_holes");
                me->printf("This function closes holes. It simply follows the loop of halfvectors which");
                me->printf("enclose the hole and add a face to which they all point.");
                return;
            }
            me->save_active_mesh();
            
            close_holes(me->active_mesh());
            return;
        }
        
        void console_reload(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  load <file>");
                me->printf("(Re)loads the current file if no argument is given, but");
                me->printf("if an argument is given, then that becomes the current file");
                return;
            }
            me->save_active_mesh();
            
            if(!me->reload_active_from_file(args.size() > 0 ? args[0]:""))
                me->printf("failed to load");
            
            return;
        }
        
        
        void console_add_mesh(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  add_mesh <file>");
                me->printf("Loads the file but without clearing the mesh. Thus, the loaded mesh is added to the");
                me->printf("current model.");
                return;
            }
            me->save_active_mesh();
            
            if(!me->add_to_active_from_file(args.size() > 0 ? args[0]:""))
                me->printf("failed to load");
            
            return;
        }
        
        void console_valid(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  validity");
                me->printf("Tests validity of Manifold");
                return;
            }
            if(valid(me->active_mesh()))
                me->printf("Mesh is valid");
            else
                me->printf("Mesh is invalid - check console output");
            return;
        }
        
        void console_Dijkstra(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  Dijkstra");
                return;
            }
            
            Manifold& m = me->active_mesh();
            
            VertexID source_vertex = *me->get_vertex_selection().begin();
            
            auto d_out = Dijkstra(m, source_vertex);
            auto& dist = d_out.dist;
            auto& pred = d_out.pred;
            
            
            auto& scalar_field = me->active_visobj().get_scalar_field_attrib_vector();
            double max_dist = 0;
            for (auto v: m.vertices())
                if(dist[v] > max_dist)
                    max_dist = dist[v];
            for(auto vid : m.vertices())
                scalar_field[vid] = 1-dist[vid]/max_dist;
            
            
            VertexAttributeVector<int> main_branch(m.allocated_vertices(),-1);
            queue<VertexID> vertex_queue;
            for(auto v: m.vertices())
                if(v != source_vertex)
                    vertex_queue.push(v);
            
            int main_branch_no = 1;
            
            while(!vertex_queue.empty()) {
            
                VertexID v = vertex_queue.front();
                vertex_queue.pop();
                
                if(main_branch[v]==-1)
                {
                    vector<VertexID> path;
                    
                    do {
                        path.push_back(v);
                        v = pred[v];
                    } while( v != source_vertex &&
                            main_branch[v] == -1);
                    
                    if(v == source_vertex)
                    {
                        for(auto w : path)
                            main_branch[w] = main_branch_no;
                        ++main_branch_no;
                    }
                    else {
                        int current_branch = main_branch[v];
                        for(auto w : path)
                            main_branch[w] = current_branch;
                    }
                }
            
            }

// Commented out since apparently unused:
//            auto int_to_col = [](int x, int xmax) {
//                Vec3f vec(((x*2)%xmax)/float(xmax),
//                          1-((x*5)%xmax)/float(xmax),
//                          ((x*11)%xmax)/float(xmax));
//                vec /= vec.max_coord();
//                return vec;
//            };
            
            auto tip = [&](VertexID tip_candidate) {
                bool is_tip = true;
                circulate_vertex_ccw(m, tip_candidate, [&](VertexID v){
                    if(pred[v] == tip_candidate)
                        is_tip = false;
                });
                return is_tip;
            };
            
            auto grad = [&](VertexID tip_v) {
                VertexID v2 = pred[tip_v];
                return normalize(m.pos(tip_v)-m.pos(v2));
            };
            
            double shortest_loop_len = DBL_MAX;
            VertexID loop_v0, loop_v1;
            for(auto h: m.halfedges()) {
                Walker w = m.walker(h);
                VertexID l0 = w.vertex();
                VertexID l1 = w.opp().vertex();
                if(tip(l0) && tip(l1) &&
                   main_branch[l0] != main_branch[l1])
                {
                    double s = dot(grad(l0),grad(l1));
                    if(s<-0.75) {
                        double loop_len = dist[l0]+dist[l1];
                        if(loop_len< shortest_loop_len)
                        {
                            shortest_loop_len = loop_len;
                            loop_v0 = l0;
                            loop_v1 = l1;
                        }
                    }
                }
                DebugRenderer::edge_colors[h] = Vec3f(0);
            }
            
            if(shortest_loop_len<DBL_MAX)
            {
                while(loop_v0 != source_vertex) {
                    circulate_vertex_ccw(m, loop_v0, [&](Walker w) {
                        VertexID vn = w.vertex();
                        if ((vn == pred[loop_v0])) {
                            DebugRenderer::edge_colors[w.halfedge()] = Vec3f(1);
                            DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                            
                        }
                    });
                    loop_v0 = pred[loop_v0];
                }
                while(loop_v1 != source_vertex) {
                    circulate_vertex_ccw(m, loop_v1, [&](Walker w) {
                        VertexID vn = w.vertex();
                        if ((vn == pred[loop_v1])) {
                            DebugRenderer::edge_colors[w.halfedge()] = Vec3f(1);
                            DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                            
                        }
                    });
                    loop_v1 = pred[loop_v1];
                }
            }
            
            
            
            for(auto v: m.vertices())
                DebugRenderer::vertex_colors[v] = Vec3f(0);
            for(auto f: m.faces())
                DebugRenderer::face_colors[f] = Vec3f(0.3);
            
            return;
        }

        const Vec3f& get_color(int i)
        {
            static Vec3f ctable[100000];
            static bool was_here;
            gel_srand(0);
            if(!was_here)
            {
                was_here = true;
                ctable[0] = Vec3f(0);
                for(int j=1;j<100000;++j)
                    ctable[j] = Vec3f(0.3)+0.7*normalize(Vec3f(gel_rand(),gel_rand(),gel_rand()));
                ctable[3] = Vec3f(1,0,0);
                ctable[4] = Vec3f(0,1,0);
                ctable[5] = Vec3f(0,0,1);
                ctable[6] = Vec3f(1,0,1);
            }
            return ctable[i%100000];
        }

        void console_info(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage:  info");
                me->printf("Provides information about mesh.");
                return;
            }
            Vec3d p0, p7;
            bbox(me->active_mesh(), p0, p7);
            stringstream bbox_corners;
            bbox_corners << p0 << " - " << p7 << endl;
            me->printf("Bounding box corners : %s", bbox_corners.str().c_str());
            map<int,int> val_hist;
            
            Manifold& m = me->active_mesh();
            
            double avg_len = 0;
            for(HalfEdgeID h: m.halfedges()) {
                DebugRenderer::edge_colors[h] = Vec3f(0.3);
                avg_len += length(m,h);
            }
            avg_len /= m.no_halfedges();
            me->printf("Avg. edge length: %f", avg_len);
            for(VertexID v: m.vertices())
            {
                int val = valency(m,v);
                DebugRenderer::vertex_colors[v] = get_color(val);
                ++val_hist[val];
                
                if(val != 4)
                    circulate_vertex_ccw(m, v, static_cast<std::function<void(HalfEdgeID)>>([&](HalfEdgeID h){
                        Walker w = m.walker(h);
                        DebugRenderer::edge_colors[h] = Vec3f(1);
                        DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                        while(valency(m, w.vertex())==4) {
                            w = w.next().opp().next();
                            DebugRenderer::edge_colors[w.halfedge()] = Vec3f(1);
                            DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1);
                        }
                    }));
            }
            map<int, int> ngon_hist;
            for(FaceID f: m.faces()) {
                int ne = no_edges(m, f);
                ++ngon_hist[ne];
                DebugRenderer::face_colors[f] = 0.7*get_color(ne);
            }
            
            me->printf("Valency histogram");
            for(map<int,int>::iterator iter = val_hist.begin(); iter != val_hist.end(); ++iter)
            {
                me->printf("%d, %d", iter->first, iter->second);
            }
            
            me->printf("Ngon histogram");
            for(map<int,int>::iterator iter = ngon_hist.begin(); iter != ngon_hist.end(); ++iter)
            {
                me->printf("%d, %d", iter->first, iter->second);
            }

            
            me->printf("Mesh contains %d faces", me->active_mesh().no_faces());
            me->printf("Mesh contains %d halfedges", me->active_mesh().no_halfedges());
            me->printf("Mesh contains %d vertices", me->active_mesh().no_vertices());
            return;
        }
        
        
        void console_simplify(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: simplify <fraction> ");
                me->printf("Performs Garland Heckbert (quadric based) mesh simplification.");
                me->printf("The only argument is the fraction of vertices to keep.");
                return;
            }
            me->save_active_mesh();
            
            float keep_fraction;
            if(args.size() == 0)
            {
                me->printf("you must specify fraction of vertices to keep");
                return;
            }
            istringstream a0(args[0]);
            a0 >> keep_fraction;

            bool optimal_positions = true;
            if(args.size() == 2)
            {
                istringstream a1(args[1]);
                a1 >> optimal_positions;
            }
            
            Vec3d p0, p7;
            Manifold m = me->active_mesh();
            bbox(m, p0, p7);
            Vec3d d = p7-p0;
            float s = d.max_coord();
            Vec3d pcentre = (p7+p0)/2.0;
            for(VertexID v: m.vertices())
            {
                m.pos(v) -= pcentre;
                m.pos(v) /= s;
            }
            cout << "Timing the Garland Heckbert (quadric based) mesh simplication..." << endl;
            Timer timer;
            timer.start();
            
            //simplify
            quadric_simplify(me->active_mesh(),keep_fraction,0.0001f,optimal_positions);
            
            cout << "Simplification complete, process time: " << timer.get_secs() << " seconds" << endl;
            
            for(VertexID v: m.vertices())
            {
                m.pos(v) *= s;
                m.pos(v) += pcentre;
            }
            return;
        }
        
        void console_vertex_noise(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args))
            {
                me->printf("usage: noise.perturb_vertices <amplitude>");
                me->printf("adds a random vector to each vertex. A random vector in the unit cube is generated and");
                me->printf("to ensure an isotropic distribution, vectors outside the unit ball are discarded.");
                me->printf("The vector is multiplied by the average edge length and then by the amplitude specified.");
                me->printf("If no amplitude is specified, the default (0.5) is used.");
                return;
            }
            me->save_active_mesh();
            
            float avg_length = average_edge_length(me->active_mesh());
            
            float noise_amplitude = 0.5f;
            if(args.size() > 0) {
                istringstream a0(args[0]);
                a0 >> noise_amplitude;
            }
            
            gel_srand(0);
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi){
                Vec3d v;
                do{
                    v = Vec3d(gel_rand(),gel_rand(),gel_rand());
                    v /= (float)(GEL_RAND_MAX);
                    v -= Vec3d(0.5);
                    v *= 2.0;
                }
                while(sqr_length(v) > 1.0);
                
                v *= noise_amplitude;
                v *= avg_length;
                me->active_mesh().pos(*vi) += v;
            }
            return;
        }
        
        void console_perpendicular_vertex_noise(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage: noise.perturb_vertices_perpendicular <amplitude>");
                me->printf("adds the normal times a random scalar times amplitude times");
                me->printf("times average edge length to the vertex. (default amplitude=0.5)");
                return;
            }
            me->save_active_mesh();
            
            float avg_length = average_edge_length(me->active_mesh());
            
            float noise_amplitude = 0.5;
            if(args.size() > 0)
            {
                istringstream a0(args[0]);
                a0 >> noise_amplitude;
            }
            
            VertexAttributeVector<Vec3d> normals(me->active_mesh().allocated_vertices());
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi)
                normals[*vi] = normal(me->active_mesh(), *vi);
            
            gel_srand(0);
            for(VertexIDIterator vi = me->active_mesh().vertices_begin(); vi != me->active_mesh().vertices_end(); ++vi)
            {
                float rval = 0.5-gel_rand() / float(GEL_RAND_MAX);
                me->active_mesh().pos(*vi) += normals[*vi]*rval*noise_amplitude*avg_length*2.0;
            }
            return;
        }
        
        void console_noisy_flips(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage:  noise.perturb_topology <iter>");
                me->printf("Perform random flips. iter (default=1) is the number of iterations.");
                me->printf("mostly for making nasty synthetic test cases.");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            
            randomize_mesh(me->active_mesh(),  iter);
            return;
        }
        
        void console_laplacian_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  smooth.laplacian <weight> <iter>");
                me->printf("Perform Laplacian smoothing. weight is the scaling factor for the Laplacian.");
                me->printf("default weight = 1.0. Default number of iterations = 1");
                return;
            }
            me->save_active_mesh();
            
            float t=1.0;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> t;
            }
            int iter = 1;
            if(args.size()>1){
                istringstream a0(args[1]);
                a0 >> iter;
            }
            Util::Timer tim;
            tim.start();
            /// Simple laplacian smoothing with an optional weight.
            laplacian_smooth(me->active_mesh(), t, iter);
            cout << "It took "<< tim.get_secs();
            return;
        }
        
        
        void console_mean_curvature_smooth(MeshEditor* me, const std::vector<std::string> & args){
            if(wantshelp(args)) {
                me->printf("usage:  smooth.mean_curvature <weight> <iter>");
                me->printf("Perform mean curvature smoothing. weight is the scaling factor for the");
                me->printf("mean curvature vector which has been normalized by dividing by edge lengths");
                me->printf("this allows for larger steps as suggested by Desbrun et al.");
                me->printf("default weight = 1.0. Default number of iterations = 1");
                return;
            }
            me->save_active_mesh();
            
            double t=1.0;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> t;
            }
            int iter=1;
            if(args.size() > 1){
                istringstream a0(args[1]);
                a0 >> iter;
            }
            VertexAttributeVector<Vec3d> new_pos(me->active_mesh().allocated_vertices());
            for(int j = 0; j < iter; ++j){
                for(VertexIDIterator v = me->active_mesh().vertices_begin(); v != me->active_mesh().vertices_end(); ++v) {
                    Vec3d m;
                    double w_sum;
                    unnormalized_mean_curvature_normal(me->active_mesh(), *v, m, w_sum);
                    new_pos[*v] = Vec3d(me->active_mesh().pos(*v))  + (t * m/w_sum);
                }
                for(VertexIDIterator v = me->active_mesh().vertices_begin(); v != me->active_mesh().vertices_end(); ++v)
                    me->active_mesh().pos(*v) = new_pos[*v];
            }
            return;
        }
        
        void console_taubin_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage:  smooth.taubin <iter>");
                me->printf("Perform Taubin smoothing. iter (default=1) is the number of iterations.");
                return;
            }
            me->save_active_mesh();
            
            int iter = 1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            /// Taubin smoothing is similar to laplacian smoothing but reduces shrinkage
            taubin_smooth(me->active_mesh(),  iter);
            
            return;
        }
        
        void console_fvm_anisotropic_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: smooth.fuzzy_vector_median <iter>");
                me->printf("Smooth normals using fuzzy vector median smoothing. iter (default=1) is the number of iterations");
                me->printf("This function does a very good job of preserving sharp edges.");
                return;
            }
            me->save_active_mesh();
            
            int iter=1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            // Fuzzy vector median smoothing is effective when it comes to preserving sharp edges.
            anisotropic_smooth(me->active_mesh(),  iter, FVM_NORMAL_SMOOTH);
            
            return;
        }
        
        void console_bilateral_anisotropic_smooth(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: smooth.fuzzy_vector_median <iter>");
                me->printf("Smooth normals using fuzzy vector median smoothing. iter (default=1) is the number of iterations");
                me->printf("This function does a very good job of preserving sharp edges.");
                return;
            }
            me->save_active_mesh();
            
            int iter=1;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> iter;
            }
            
            anisotropic_smooth(me->active_mesh(),  iter, BILATERAL_NORMAL_SMOOTH);
            
            return;
        }
        
        void console_triangulate(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  triangulate");
                me->printf("This function triangulates all non triangular faces of the mesh.");
                me->printf("you may want to call it after hole closing. For a polygon it simply connects");
                me->printf("the two closest vertices in a recursive manner until only triangles remain");
                return;
            }
            me->save_active_mesh();
            
            shortest_edge_triangulate(me->active_mesh());
            me->active_mesh().cleanup();
            valid(me->active_mesh());
            return;
        }
        
        
        void console_remove_caps(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  cleanup.remove_caps thresh");
                me->printf("Remove caps (triangles with one very big angle). The thresh argument is the fraction of PI to");
                me->printf("use as threshold for big angle. Default is 0.85. Caps are removed by flipping.");
                return;
            }
            me->save_active_mesh();
            
            float t = 0.85f;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> t;
            }
            remove_caps(me->active_mesh(), static_cast<float>(M_PI) *t);
            me->active_mesh().cleanup();
            
            return;
        }
        
        void console_remove_needles(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: cleanup.remove_needles <thresh>");
                me->printf("Removes very short edges by collapse. thresh is multiplied by the average edge length");
                me->printf("to get the length shorter than which we collapse. Default = 0.1");
                return;
            }
            me->save_active_mesh();
            
            float thresh = 0.1f;
            if(args.size() > 0){
                istringstream a0(args[0]);
                a0 >> thresh;
            }
            float avg_length = average_edge_length(me->active_mesh());
            remove_needles(me->active_mesh(), thresh * avg_length);
            me->active_mesh().cleanup();
            
            return;
        }
        
        void console_flip_orientation(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)){
                me->printf("usage: cleanup.flip_orientation");
                me->printf("reorients all faces - flipping normal direction");
                return;
            }
            me->save_active_mesh();
            flip_orientation(me->active_mesh());
            return;
        }

        
        void console_undo(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  undo");
                me->printf("This function undoes one operation. Repeated undo does nothing");
                return;
            }
            me->restore_active_mesh();
            return;
        }
        
        void console_save_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  display.save_trackball");
                me->printf("This function saves the trackball to disk");
                return;
            }
            me->save_ball();
            return;
        }

        void console_load_trackball(MeshEditor* me, const std::vector<std::string> & args)
        {
            if(wantshelp(args)) {
                me->printf("usage:  display.load_trackball");
                me->printf("This function loads the trackball from disk");
                return;
            }
            me->load_ball();
            return;
        }

        
    }
    
    
    void MeshEditor::register_console_function(const std::string& name,
                                               const std::function<void(MeshEditor*, const std::vector<std::string>&)>& con_fun,
                                               const std::string& help_txt)
    {
        std::function<void (const std::vector<std::string>&)> f = bind(con_fun, this, placeholders::_1);
        theConsole.reg_cmdN(name, f, help_txt);
    }
    
    void MeshEditor::keyparse(unsigned short key){
        //toggle console with ESC
        if (key == 27)
        {
            console_visible = !console_visible;
            return;
        }
        
        if (console_visible)
        {
            static Timer tim;
            if(key==13)
                tim.start();
            theConsole.keyboard(key);
            if(key == 13)
            {
                active_visobj().post_create_display_list();
                double t = tim.get_secs();
                printf("%f seconds",t);
            }
            return;
        }
        else {
            
            switch(key) {
                    case 'R':
                {
                    theConsole.key_up();
                    theConsole.keyboard(13);
                    active_visobj().post_create_display_list();
                }
                    break;
                case 'q': exit(0);
                case '\033':
                    console_visible = false;
                    break;
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                {
                    int w[4];
                    glGetIntegerv(GL_VIEWPORT, w);
                    active = key - '1';
                    active_view_control().reshape(w[2],w[3]);
                }
                    break;
                case 'f': display_smooth_shading = !display_smooth_shading; break;
                case 'w':
                    display_render_mode = "wire"; break;
                case 'n':
                    display_render_mode = "normal"; break;
                case 'i':
                    display_render_mode = "isophotes"; break;
                case 'r':
                    display_render_mode = "reflection"; break;
//                case 'h':
//                    display_render_mode = "harmonics"; break;
                case 't':
                    display_render_mode = "toon"; break;
                case 'g':
                    display_render_mode = "glazed"; break;
                case 'a':
                    display_render_mode = "ambient_occlusion"; break;
                case 'c':
                    display_render_mode = "copper"; break;
                case 's':
                    display_render_mode = "scalar"; break;
                case 'd':
                    display_render_mode = "debug"; break;
                case 'C':
                    display_render_mode = "curvature_lines"; break;
                case 'M':
                    display_render_mode = "mean_curvature"; break;
                case 'G':
                    display_render_mode = "gaussian_curvature"; break;
                case ' ':
                    active_visobj().clear_vertex_selection();
                    active_visobj().clear_face_selection();
                    active_visobj().clear_halfedge_selection();
                    break;
                
            }
            
            if(key != '\033') post_create_display_list();
        }
        
    }
    
    void MeshEditor::printf(const char* format, ...)
    {
        //format text
        char buffer[1024];
        va_list args;
        va_start(args, format);
        vsprintf(buffer, format, args);
        va_end(args);
        theConsole.print(buffer);
    }

    void MeshEditor::key_up(){theConsole.key_up();}
    void MeshEditor::key_down(){theConsole.key_down();}
    void MeshEditor::key_left(){theConsole.key_left();}
    void MeshEditor::key_right(){theConsole.key_right();}
    void MeshEditor::key_home(){theConsole.key_home();}
    void MeshEditor::key_end(){theConsole.key_end();}
    
    
    
    void MeshEditor::grab_ball(TrackBallAction action, const CGLA::Vec2i& pos){
        active_view_control().grab_ball(action, pos);
    }
    void MeshEditor::roll_ball(const CGLA::Vec2i& pos){
        active_view_control().roll_ball(pos);
    }
    void MeshEditor::release_ball(){
        active_view_control().release_ball();
    }
    bool MeshEditor::try_spinning_ball(){
        return active_view_control().try_spin();
    }
    
    void MeshEditor::save_ball() {
        ofstream ofs("trackball.bin", ios_base::binary);
        active_view_control().save(ofs);
        
    }
    void MeshEditor::load_ball() {
        ifstream ifs("trackball.bin", ios_base::binary);
        active_view_control().load(ifs);
    
    }

    
    bool MeshEditor::grab_mesh(const CGLA::Vec2i& pos)
    {
        if(depth_pick(pos[0], pos[1], depth))
        {
            dragging = true;
            mouse_x = pos[0];
            mouse_y = pos[1];
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Manifold& m = active_mesh();
            active_visobj().save_old();
            Vec3d c;
            float r;
            bsphere(m, c, r);
            for(auto vid : m.vertices())
            {
                double l = sqr_length(p0-m.pos(vid));
                weight_vector[vid] = exp(-l/(brush_size*r*r));
            }
            return true;
        }
        return false;
    }
    
    bool MeshEditor::drag_mesh(const CGLA::Vec2i& pos)
    {
        auto deform_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, Vec3d& v)
        {
            const Manifold& m_old = active_visobj().mesh_old();
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m_old.vertices())
                new_pos[vid] = m_old.pos(vid) + weight_vector[vid] * v;
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto inflate_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * (0.25 * laplacian(m, vid) + (r*0.0025)*normal(m, vid));
            }
            m.positions_attribute_vector() = new_pos;
        };
        
        
        auto smooth_mesh = [&](Manifold& m, VertexAttributeVector<float>& wv, const Vec3d& p0)
        {
            Vec3d c;
            float r;
            bsphere(m, c, r);
            VertexAttributeVector<Vec3d> new_pos;
            for(auto vid : m.vertices()) {
                Vec3d p = m.pos(vid);
                double l = sqr_length(p0-p);
                double wgt = exp(-l/(brush_size*r*r));
                new_pos[vid] = m.pos(vid) +
                wgt * 0.5 * laplacian(m, vid);
            }
            m.positions_attribute_vector() = new_pos;
        };

        
        if(dragging)
        {
            Vec3d p0 = screen2world(mouse_x, mouse_y, depth);
            Vec3d p1 = screen2world(pos[0], pos[1], depth);
            Vec3d v = p1-p0;
            Manifold& m = active_mesh();
            auto vset = active_visobj().get_vertex_selection();
            auto hset = active_visobj().get_halfedge_selection();
            auto fset = active_visobj().get_face_selection();
            
            if(!vset.empty())
                for(auto vid : vset)
                    m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
            else if(!hset.empty())
                for(auto h : hset) {
                    Walker w = m.walker(h);
                    auto vid0 = w.vertex();
                    auto vid1 = w.opp().vertex();
                    m.pos(vid0) = active_visobj().mesh_old().pos(vid0) + v;
                    m.pos(vid1) = active_visobj().mesh_old().pos(vid1) + v;
                }
            else if(!fset.empty())
                for(auto f: fset)
                {
                    circulate_face_ccw(m, f, [&](VertexID vid){
                        m.pos(vid) = active_visobj().mesh_old().pos(vid) + v;
                    });
                }
            else {
                if(string(brush_type) == "smooth")
                    smooth_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "inflate")
                    inflate_mesh(m, weight_vector, p1);
                else if(string(brush_type) == "deform")
                    deform_mesh(m, weight_vector, v);
            }
            
            post_create_display_list();
            return true;
        }
        return false;
    }
    
    void MeshEditor::release_mesh()
    {
        dragging = false;
    }
    
    
    
    
//    mutex parallel_work;
    void MeshEditor::display(int scale){

        // Clear screen.
        glClearColor(1, 1, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        //        parallel_work.lock();
        //        active_visobj().post_create_display_list();
        
        
        // Display object.
        active_visobj().display(display_render_mode, theConsole, display_smooth_shading, display_gamma);
        
        // Draw console
        if(console_visible)
        {
            glUseProgram(0);
            theConsole.display(scale);
        }
        
        // Static variable controlling whether we render shadow at all.
        static Console::variable<int> shadow_enable(0);
        shadow_enable.reg(theConsole, "display.shadow.enable", "");
        if(shadow_enable) {
            // Static variables that control the shadow display created and linked to console
            static Console::variable<float> zenith(1.571);
            static Console::variable<float> azimuth(0);
            static Console::variable<float> shadow_alpha(0.3);
            azimuth.reg(theConsole, "display.shadow.azimuth", "");
            zenith.reg(theConsole, "display.shadow.zenith", "");
            shadow_alpha.reg(theConsole, "display.shadow.shadow_alpha", "");

            
            // Shadow buffer (really a frame buffer object)
            static ShadowBuffer sb;
            
            // String containing fragment program for shadow rendering
            static string shadow_shdr_fp = "#version 120\n"
            "    uniform sampler2DShadow shadow_map;\n"
            "    uniform mat4 Mat;\n"
            "   uniform float shadow_alpha;\n"
            "    varying vec4 ep;\n"
            "    void main()\n"
            "    {\n"
            "        vec4 light_pos =  Mat * ep;\n"
            "        light_pos.z = max(0.001,min(0.999,light_pos.z-0.003));\n"
            "        if(light_pos.x <=0 || light_pos.x >=1 || light_pos.y <=0 || light_pos.y>=1) gl_FragColor = vec4(1);"
            "   else      gl_FragColor= vec4(1.0-shadow_alpha)+shadow_alpha*0.25*(shadow2D(shadow_map, light_pos.xyz+vec3(0.001,-0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz-vec3(0.001,0.001,0))+shadow2D(shadow_map, light_pos.xyz+vec3(-0.001,0.001,0)));\n"
            "    }\n";
            
            // Shader program for shadow rendering is compiled and linked.
            static GLuint prog = 0;
            if(!prog)
            {
                GLuint vp = create_glsl_shader(GL_VERTEX_SHADER,"#version 120\n varying vec4 ep; void main(){ep = gl_Vertex; gl_Position = ftransform();}");
                GLuint fp = create_glsl_shader(GL_FRAGMENT_SHADER, shadow_shdr_fp);
                prog = glCreateProgram();
                
                glAttachShader(prog, vp);
                glAttachShader(prog, fp);
                glLinkProgram(prog);
            }
            
            
            // Setup OpenGL state for lighting - used when rendering 3D object casting shadow
            Vec4f lpos = Vec4f(cos(azimuth)*cos(zenith),sin(zenith),sin(azimuth)*cos(zenith),0);
            glLightfv(GL_LIGHT0, GL_POSITION, lpos.get());
            Vec4f mamb(.8,.8,.8,0);
            Vec4f lamb(.4,.4,.5,0);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mamb.get());
            glLightfv(GL_LIGHT0, GL_AMBIENT, lamb.get());

            // Get old viewport dimensions. Setup rendering to FBO and set viewport
            // dimensions for rendering shadows.
            GLint viewp[4];
            glGetIntegerv(GL_VIEWPORT, viewp);
            sb.enable();
            glViewport(0, 0, 1024, 1024);
            
            // Setup object transformations using old school GL.
            Mat4x4f m;
            float r = active_visobj().get_bsphere_radius();
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            glOrtho(-2*r, 2*r, -2*r, 2*r, -2*r, 2*r);
            gluLookAt(0.1*r*cos(azimuth)*cos(zenith),0.1*r*sin(zenith),0.1*r*sin(azimuth)*cos(zenith), 0,0,0, 0,1,0);
            
            // Copy the transformation matrix to user code.
            glGetFloatv(GL_PROJECTION_MATRIX, m.get());
            
            // Draw the object in light space.
            draw(active_visobj().mesh());
            
            // Restore transformation matrices.
            glPopMatrix();
            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();
            
            // Restore the usual on screen framebuffer.
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);
            glViewport(0, 0, viewp[2], viewp[3]);
            
            // Use shadow: bind the shadow texture to unit 0.
            sb.bind_textures(0);
            
            // Use the shadow rendering shader program.
            glUseProgram(prog);
            active_visobj().view_control().set_gl_modelview();
            glUniform1i(glGetUniformLocation(prog, "shadow_map"), 0);
            glUniform1f(glGetUniformLocation(prog, "shadow_alpha"), shadow_alpha);
            m = translation_Mat4x4f(Vec3f(0.5)) * scaling_Mat4x4f(Vec3f(0.5)) * transpose(m);
            glUniformMatrix4fv(glGetUniformLocation(prog, "Mat"), 1, 1, m.get());
            
            
            // Setup blending such that the shadow is alpha blended with model.
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ZERO, GL_SRC_COLOR);

            // Draw ground plane for shadows.
            Vec3d p0, p7;
            bbox(active_visobj().mesh(), p0, p7);
            glBegin(GL_QUADS);
            glVertex3f(-100*r, p0[1],-100*r);
            glVertex3f(-100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],100*r);
            glVertex3f(100*r, p0[1],-100*r);
            glEnd();
            
            // Draw model again ... just to add shadow.
            draw(active_visobj().mesh());
            
            // Disable blending and shader program.
            glDisable(GL_BLEND);
            glUseProgram(0);
        }
        //        parallel_work.unlock();
        //        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    void MeshEditor::reshape(int w, int h) {
        for(VisObj& v : vo)
            v.view_control().reshape(w, h);
    }

    void MeshEditor::init() {
        glewInit();
        
        GLint vp[4];
        glGetIntegerv(GL_VIEWPORT, vp);
        for(VisObj& vis_obj : vo)
            vis_obj.view_control().reshape(vp[2], vp[3]);
        
//        glEnable(GL_CULL_FACE);
//        glCullFace(GL_BACK);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glClearColor(1,1,0, 0.f);
        glColor4f(1.0f, 1.0f, 1.0f, 0.f);
        float material[4] = {1,1,1,1};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material);
        glEnable(GL_DEPTH_TEST);
        
        register_console_function("simplify", console_simplify,"");
        
        register_console_function("ridge_lines", console_ridge_lines,"");
        
        register_console_function("smooth.mean_curvature", console_mean_curvature_smooth,"");
        register_console_function("smooth.laplacian", console_laplacian_smooth,"");
        register_console_function("smooth.taubin", console_taubin_smooth,"");
        register_console_function("smooth.fuzzy_vector_median_anisotropic", console_fvm_anisotropic_smooth ,"");
        register_console_function("smooth.bilateral_anisotropic", console_bilateral_anisotropic_smooth ,"");
        
        register_console_function("optimize.valency", console_optimize_valency,"");
        register_console_function("optimize.minimize_dihedral_angles", console_minimize_dihedral,"");
        register_console_function("optimize.minimize_curvature", console_minimize_curvature,"");
        register_console_function("optimize.maximize_min_angle", console_maximize_min_angle,"");
        register_console_function("load_mesh", console_reload,"");
        register_console_function("add_mesh", console_add_mesh,"");
        
        register_console_function("cleanup.stitch", console_stitch,"");
        register_console_function("cleanup.remove_duplicates", console_remove_duplicates,"");
        register_console_function("cleanup.remove_val2", console_remove_val2, "");
        register_console_function("cleanup.flip_orientation", console_flip_orientation,"");
        register_console_function("cleanup.remove_caps", console_remove_caps,"");
        register_console_function("cleanup.remove_needles", console_remove_needles,"");
        register_console_function("cleanup.close_holes", console_close_holes,"");
        register_console_function("cleanup.compact", console_compact,"");

        register_console_function("triangulate", console_triangulate,"");
        register_console_function("refine.split_edges", console_refine_edges,"");
        register_console_function("refine.split_faces", console_refine_faces,"");
        register_console_function("subdivide.catmull_clark", console_cc_subdivide,"");
        register_console_function("subdivide.rootcc", console_root_cc_subdivide,"");
        register_console_function("subdivide.loop", console_loop_subdivide,"");
        register_console_function("subdivide.root3", console_root3_subdivide,"");
        register_console_function("subdivide.doo_sabin", console_doosabin_subdivide,"");
        register_console_function("subdivide.butterfly", console_butterfly_subdivide,"");
        register_console_function("save_mesh", console_save,"");
        register_console_function("noise.perturb_vertices", console_vertex_noise,"");
        register_console_function("noise.perturb_vertices_perpendicular", console_perpendicular_vertex_noise,"");
        register_console_function("noise.perturb_topology", console_noisy_flips,"");
        
        register_console_function("dual", console_dual,"");
//        register_console_function("flatten", console_flatten,"");
        
        register_console_function("align", console_align,"");
        register_console_function("undo", console_undo,"");
        
        register_console_function("validity", console_valid,"");
        register_console_function("info", console_info,"");
        
//        register_console_function("harmonics.reset_shape", console_reset_shape, "");
//        register_console_function("harmonics.analyze", console_analyze, "");
//        register_console_function("harmonics.partial_reconstruct", console_partial_reconstruct,"");

        register_console_function("Dijkstra", console_Dijkstra,"");
        
        register_console_function("display.refit_trackball", console_refit_trackball, "Resets trackball");
        register_console_function("display.save_trackball", console_save_trackball, "Saves trackball to disk");
        register_console_function("display.load_trackball", console_load_trackball, "Load trackball to disk");
        
        register_console_function("transform.scale", console_scale, "Scale mesh");
        register_console_function("transform.rotate", console_rotate, "Rotate mesh");
        register_console_function("transform.translate", console_translate, "Translate mesh");
        
        register_console_function("edit.selected.merge_1_ring", console_merge_1_ring, "Merge 1-ring of selected vertices");
        
        register_console_function("edit.selected.split_edge", console_split_edge , "");
        register_console_function("edit.selected.collapse_edge", console_collapse_edge, "");
        register_console_function("edit.selected.flip_edge", console_flip_edge, "");
        register_console_function("edit.selected.dissolve_edge", console_dissolve_edge, "");
        register_console_function("edit.selected.stellate_face", console_stellate_face, "");
        register_console_function("edit.selected.split_face", console_split_face, "");
        register_console_function("edit.selected.delete", console_delete, "");
        register_console_function("edit.selected.triangulate_face", console_triangulate_face, "");
        register_console_function("edit.selected.bridge_faces", console_bridge_faces, "");

        register_console_function("test", console_test, "Test some shit");
        
        selection_mode.reg(theConsole, "selection.mode", "The selection mode. 0 = vertex, 1 = halfedge, 2 = face");
        active.reg(theConsole, "active_mesh", "The active mesh");
        display_render_mode.reg(theConsole, "display.render_mode", "Display render mode");
        brush_size.reg(theConsole, "brush.size", "Size of brush used for editing");
        brush_type.reg(theConsole, "brush.type", "smooth, deform, or inflate");
        
        display_smooth_shading.reg(theConsole, "display.smooth_shading", "1 for smooth shading 0 for flat");
        display_gamma.reg(theConsole, "display.gamma", "The gamma setting for the display");

        theConsole.print("Welcome to MeshEdit");
        theConsole.newline();
    }
    
    bool MeshEditor::add_file(const std::string& str)
    {
        while (active_mesh().no_vertices()>0 && active<NO_MESHES)
            active  = active + 1;
        if(active == NO_MESHES)
            active = 0;
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }

    bool MeshEditor::reload_active_from_file(const std::string& str)
    {
        if(active_visobj().reload(str)) {
            active_visobj().post_create_display_list();
            return true;
        }
        return false;
    }
    
    bool MeshEditor::add_to_active_from_file(const std::string& str)
    {
        if(active_visobj().add_mesh(str)) {
            active_visobj().post_create_display_list();
            return  true;
        }
        return false;
    }

    
}
