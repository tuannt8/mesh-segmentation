//
//  interface_dsc.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/9/15.
//  Copyright (c) 2015 Asger Nyman Christiansen. All rights reserved.
//
#include "interface.h"

#include "trializer.h"
#include "object_generator.h"
#include "draw.h"

#include "dyn_integral.h"

#include "gl_debug_helper.h"
#include "dynamics_edge.h"
#include "adapt_mesh.h"

#include "options_disp.h"

#include "profile.h"
#include "define.h"

int ox, oy, w,h;

void _check_gl_error(const char *file, int line)
{
    GLenum err (glGetError());
    
    while(err!=GL_NO_ERROR) {
        std::string error;
        
        switch(err) {
            case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
            case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
            case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
            case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
        }
        
        std::cerr << "GL_" << error.c_str() <<" - "<<file<<":"<<line<<std::endl;
        err=glGetError();
    }
}

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)


void display_(){
    interface_dsc::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    interface_dsc::get_instance()->keyboard(key, x, y);
}

void reshape_(int width, int height){
    interface_dsc::get_instance()->reshape(width, height);
}

void visible_(int v){
    interface_dsc::get_instance()->visible(v);
}

void animate_(){
    interface_dsc::get_instance()->animate();
}

void glutMouseFunc_(int button, int state, int x, int y){
    gl_debug_helper::mouseDown(button, state, x, y);
    options_disp::mouse_func(button, state, x, y);
};;


void glutMotion_(int x, int y){
    gl_debug_helper::mouseMove(x, y);
};

interface_dsc* interface_dsc::instance = NULL;

#pragma mark - Glut display
void interface_dsc::display(){
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
    

    
    draw();
    update_title();
    
    check_gl_error();
    
    if (RUN) {
        dynamics_image_seg();
        glutPostRedisplay();
    }
}

void interface_dsc::animate(){
//    if (RUN) {
//        dynamics_image_seg();
//        glutPostRedisplay();
//    }
//    else{
//        sleep(0.7);
//    }
}

void interface_dsc::reshape(int width, int height){
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    
    double real_width = width - options_disp::width_view;
    if(dsc)
    {
        double image_ratio = imageSize[1] / imageSize[0];
        double gl_ratio = (double)WIN_SIZE_Y / real_width;
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluOrtho2D(0, imageSize[0], 0, imageSize[1]);
//        gluOrtho2D(-DISCRETIZE_RES, imageSize[0]-DISCRETIZE_RES, -DISCRETIZE_RES, imageSize[1]-DISCRETIZE_RES);
        
        double lx = (gl_ratio < image_ratio)? WIN_SIZE_Y/image_ratio : real_width;
        double ly = (gl_ratio < image_ratio)? WIN_SIZE_Y : real_width*image_ratio;


        glViewport(options_disp::width_view + (real_width - lx)/2, (WIN_SIZE_Y - ly)/2, lx, ly);
      //  glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
        

        ox = options_disp::width_view + (real_width - lx)/2 + DISCRETIZE_RES / (imageSize[0] / lx);
        oy = (WIN_SIZE_Y - ly)/2 + DISCRETIZE_RES / (imageSize[1] / ly);
        w = lx - 2*DISCRETIZE_RES / (imageSize[0] / lx) + 1;
        h = ly - 2*DISCRETIZE_RES / (imageSize[1] / ly) + 1;

        gl_debug_helper::coord_transform(Vec2(ox,oy),
                                         Vec2(imageSize[0] / lx, imageSize[1] / ly), WIN_SIZE_Y);
    }
}

void interface_dsc::visible(int v){
//    if(v==GLUT_VISIBLE){
//        glutIdleFunc(animate_);
//    }
//    else
//        glutIdleFunc(0);
}



void interface_dsc::keyboard(unsigned char key, int x, int y){
    
    switch (key) {
        case 'd':
            gl_debug_helper::change_state();
            break;
        case 'e':
            export_dsc();
            break;
        case 'r':
        {
            profile::close();
        }
            break;
        case ' ':
        {
            RUN = !RUN;
            update_title();
            //dynamics_image_seg();
            break;
        }

        case '\t':
            Painter::save_painting_dsc(ox, oy, w, h, "./LOG");
            break;
        case 'i':
            dyn_->write_energy();
            break;
        case 'f': // Flipping phase
        {
            adapt_mesh am;
//            am.split_face(*dsc, *image_);
            am.split_face_and_relabel(*dsc, *image_);
        }
            break;
        case 's': // Split edge
        {
            static int c = 1;
            std::ostringstream s;
            s << "DATA/tues/im_" << c++ << ".png";
            image_->load_image(s.str());
        }
            break;
        case 'b': // Split edge
        {
            back_up_dsc();
        }
            break;
        case 'l': // Split edge
        {
            load_dsc();
        }
            break;
//        case 'e': // Split edge
//        {
//            export_mesh();
//        }
//            break;
        case 'w': // Split edge
        {
            write_triangle_energy();
        }
            break;
        case 'u':
            std::cout << g_param.alpha << " - New alpha: ";
            std::cin >> g_param.alpha;
            std::cout << g_param.beta << " - New beta: ";
            std::cin >> g_param.beta;
            break;
        default:
            break;
    }
    
    int dis = (int)key - 48;
    if (dis <= 10 and dis >= 0) {
//        bDiplay_[dis] = ! bDiplay_[dis];
        gl_debug_helper::_label_idx = dis;
    }
    g_param.bDisplay = bDiplay_;
    
    glutPostRedisplay();
}



using namespace DSC2D;
void interface_dsc::load_dsc()
{
    std::ostringstream os;
    os << LOG_PATH << "dsc.dsc";
    
    std::ifstream myfile(os.str().c_str());
    
    if (myfile.is_open()) {
        int nb_vertice, nb_face;
        myfile >> nb_vertice;
        myfile >> nb_face;
        
        std::vector<double> points(3*nb_vertice, 0);
        for (int i = 0; i < nb_vertice; i++) {
            myfile >> points[3*i];
            myfile >> points[3*i+1];
        }
        std::vector<int> faces; faces.resize(3*nb_face);
        std::vector<int> labels; labels.resize(nb_face);
        for (int i = 0; i < nb_face; i++) {
            myfile >> faces[3*i];
            myfile >> faces[3*i + 1];
            myfile >> faces[3*i + 2];
            
            myfile >> labels[i];
        }
        
        // Init DSC
        double width = imageSize[0];
        double height = imageSize[1];
        
        DISCRETIZATION = (double) height / (double)DISCRETIZE_RES;
        
        width -= 2*DISCRETIZATION;
        height -= 2*DISCRETIZATION;
        
        DesignDomain *domain = new DesignDomain(DesignDomain::RECTANGLE, width, height, 0 /* DISCRETIZATION */);
        
        dsc = std::unique_ptr<DeformableSimplicialComplex>(
                                                           new DeformableSimplicialComplex(DISCRETIZATION, points, faces, domain));
        
        // set label
        int idx = 0;
        for (auto fiter = dsc->faces_begin(); fiter != dsc->faces_end(); fiter++)
        {
            dsc->update_attributes(*fiter, labels[idx++]);
        }
            
            
#ifdef TUAN_MULTI_RES
        dsc->img = &*image_;
#endif
        gl_debug_helper::set_dsc(&(*dsc));
        
        myfile.close();
    }else{
        std::cout << "Fail to load dsc mesh \n";
    }
}

void interface_dsc::export_mesh()
{
    std::ostringstream os;
    os << LOG_PATH << "mesh.dsc";
    
    std::ofstream myfile(os.str());
    if (myfile.is_open()) {
        myfile << dsc->get_no_vertices() << " " << dsc->get_no_faces() << "\n";
        // Write vertices
        std::map<int,int> index_map;
        int idx = 0;
        for (auto vkey : dsc->vertices())
        {
            index_map.insert(std::make_pair(vkey.get_index(), idx++));
            auto p = dsc->get_pos(vkey);
            myfile << p[0] << " " << p[1] << "\n";
        }
        
        // write face
        for (auto fkey : dsc->faces())
        {
            auto verts = dsc->get_verts(fkey);
            myfile << index_map[(int)verts[0].get_index()] << " "
            << index_map[(int)verts[1].get_index()] << " "
            << index_map[(int)verts[2].get_index()] << " "
            << dsc->get_label(fkey) << "\n";
        }
        
        myfile.close();
    }
}

void interface_dsc::back_up_dsc()
{
    std::ostringstream os;
    os << LOG_PATH << "dsc.dsc";
    
    dsc->save(os.str().c_str());
}

void interface_dsc::initGL(){
    glutInitWindowSize(WIN_SIZE_X,WIN_SIZE_Y);
#ifdef WIN32
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE);
#else
    glutInitDisplayString("rgba double samples=16");
#endif
    glutCreateWindow("");
    
    glEnable(GL_MULTISAMPLE);
    
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    glutMouseFunc(glutMouseFunc_);
    glutMotionFunc(glutMotion_);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

#pragma mark - interface_dsc
void interface_dsc::draw()
{
    Painter::begin();
    
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
    
    
    if (options_disp::get_option("Image", true)) {
        image_->draw_image(WIN_SIZE_X);
    }
    else
        glColor3f(0.4, 0.4, 0.4);
    
    if (options_disp::get_option("DSC faces", true) and dsc) {
        Painter::draw_faces(*dsc);
    }
    
//    draw_test();
//    
//    if(options_disp::get_option("Image gradient", false)){
//        image_->draw_grad(WIN_SIZE_X);
//    }
//    
    if (options_disp::get_option("Face intensity", false) and dsc) {
        if(g_param.mean_intensity.size()==0)
            dyn_->compute_mean_intensity(*dsc, *image_);
        
        Painter::draw_faces_intensity(*dsc);
    }
    
    if (options_disp::get_option("Edge and vertices ", true) and dsc) {
        glLineWidth(1.0);
        Painter::draw_edges(*dsc);
        glColor3f(1, 0.0, 0.0);
        glPointSize(1.0);
    //    Painter::draw_vertices(*dsc);
    }
    
//    if(options_disp::get_option("Phase index", false)){
//        Painter::draw_face_label(*dsc);
//    }
//    
    if(options_disp::get_option("Triangle variation", false)){
        draw_tri_variant();
    }
    
    if(options_disp::get_option("Triangle MS energy", false)){
        draw_tri_MS_energy();
    }
//
//    
//    if(options_disp::get_option("Edge energy", false)){
//        draw_edge_energy();
//    }
//    
//    if(options_disp::get_option("Vertices index", false)){
//        glColor3f(1, 0, 0);
//        Painter::draw_vertices_index(*dsc);
//    }
    
//    if(options_disp::get_option("Edge index", false)){
//        glColor3f(1, 0, 0);
//        // Painter::draw_vertices_index(*dsc);
//        // Painter::draw_faces_index(*dsc);
//        Painter::draw_edges_index(*dsc);
//    }
    
//    if (options_disp::get_option("Face index")){
//        Painter::draw_faces_index(*dsc);
//    }
//    
//    if (options_disp::get_option("Node external force")) {
//        Painter::draw_external_force(*dsc);
//    }
    
    // Debug
//    auto fid = dsc->faces_begin();
//    auto pts = dsc->get_pos(*fid);
//    image_->debug_integral(pts);
    
    
    gl_debug_helper::draw();
    
    options_disp::draw(WIN_SIZE_X, WIN_SIZE_Y);
    
    Painter::end();
}

void interface_dsc::draw_tri_MS_energy(){
    // Energy base on mumfordshah energy
    dyn_->compute_mean_intensity(*dsc, *image_);
    auto mean_inten = g_param.mean_intensity;
    for (auto fkey : dsc->faces()){
        auto pts = dsc->get_pos(fkey);
//        /*
//         Shrink the triangle
//         */
//        {
//            auto center = (pts[0] + pts[1] + pts[2]) / 3.0;
//            double aa = 1 - 0.1;
//            pts[0] = center + (pts[0] - center)*aa;
//            pts[1] = center + (pts[1] - center)*aa;
//            pts[2] = center + (pts[2] - center)*aa;
//        }
//        /**/
        double area = dsc->area(fkey);
        double mi = mean_inten[dsc->get_label(fkey)];
        double e = image_->get_tri_differ_f(pts, mi)/ (area + SINGULAR_AREA);
        
        auto center = (pts[0] + pts[1] + pts[2])/3.0;
        
        std::ostringstream is;
        is.precision(3);
        is << e;
        Painter::print_gl(center[0], center[1], is.str().c_str());
    }
}

void interface_dsc::draw_tri_variant(){

        for (auto fkey : dsc->faces()){
            auto pts = dsc->get_pos(fkey);
//            /*
//             Shrink the triangle
//             */
//            {
//            auto center = (pts[0] + pts[1] + pts[2]) / 3.0;
//            double aa = 1 - 0.1;
//            pts[0] = center + (pts[0] - center)*aa;
//            pts[1] = center + (pts[1] - center)*aa;
//            pts[2] = center + (pts[2] - center)*aa;
//            }
//            /**/
            double area;
            double mi = image_->get_tri_intensity_f(pts, &area); mi /= area;
            double e = image_->get_tri_differ_f(pts, mi)/ (area + SINGULAR_AREA);
            
            auto center = (pts[0] + pts[1] + pts[2])/3.0;
            
            std::ostringstream is;
            is.precision(3);
            is << e;
            Painter::print_gl(center[0], center[1], is.str().c_str());
        }
}

void interface_dsc::write_triangle_energy()
{
    std::ofstream f_variant("variant.txt");
    std::ofstream f_ms("MS.txt");
    
    dyn_->compute_mean_intensity(*dsc, *image_);
    auto mean_inten = g_param.mean_intensity;
    
    for (auto fkey : dsc->faces()){
        
        if (dsc->get_label(fkey) == BOUND_FACE)
        {
            continue;
        }
        
        auto pts = dsc->get_pos(fkey);
        
        // variency
        {
            double area;
            double mi = image_->get_tri_intensity_f(pts, &area); mi /= area;
            double e = image_->get_tri_differ_f(pts, mi)/ (area + SINGULAR_AREA);
            
            f_variant << e << std::endl;
        }
        
        // Mumford-Shah energy
        {
            double area = dsc->area(fkey);
            double mi = mean_inten[dsc->get_label(fkey)];
            double e = image_->get_tri_differ_f(pts, mi)/ (area + SINGULAR_AREA);
            f_ms << e << std::endl;
        }
    }
    
    f_variant.close();
    f_ms.close();
}

void interface_dsc::draw_edge_energy(){
    
    // By total force energy
    if (1) {
        
        glColor3f(0, 0, 0);
        
        std::vector<Edge_key> edges;
        for(auto hei = dsc->halfedges_begin(); hei != dsc->halfedges_end(); ++hei)
        {
            if (dsc->is_interface(*hei)) {
                auto hew = dsc->walker(*hei);
                if(dsc->is_movable(*hei)
                   and dsc->get_label(hew.face()) < dsc->get_label(hew.opp().face()))
                {
                    edges.push_back(*hei);
                }
            }
        }
        
        auto mean_inten_ = g_param.mean_intensity;
        
        for (auto ekey : edges){
            auto hew = dsc->walker(ekey);
            
            double ev = 0;
            double c0 = mean_inten_[dsc->get_label(hew.face())];
            double c1 = mean_inten_[dsc->get_label(hew.opp().face())];
            
            // Loop on the edge
            auto p0 = dsc->get_pos(hew.opp().vertex());
            auto p1 = dsc->get_pos(hew.vertex());
            
            
                double length = (p1 - p0).length();
                int N = (int)length;
                double dl = length/(double)N;
                for (int i = 0; i <= N; i++) {
                    auto p = p0 + (p1 - p0)*(i/(double)N)*dl;
                    double I = image_->get_intensity_f(p[0], p[1]);
                    
                    // Normalize force
                    double f = (2*I - c0 - c1) / (c0-c1);
                    
                    ev += std::abs(f)*dl;
                }
            
            
            ev = ev / ((double)length + 3);
            
            // draw
            auto tris = dsc->get_pos(ekey);
            auto center = (tris[0] + tris[1])/2.0;
            std::ostringstream is;
            is << ev;
            Painter::print_gl(center[0], center[1], is.str().c_str());
        }
    }
}


void interface_dsc::draw_test(){
    
    if (options_disp::get_option("Face energy", false)) {
        if (g_param.mean_intensity.size() == 0) {
            return;
        }
        
        HMesh::FaceAttributeVector<Vec3> intensity(dsc->get_no_faces(), Vec3(0.0));
        glColor3f(0, 0, 0);
        for (auto fkey : dsc->faces())
        {
            auto tris = dsc->get_pos(fkey);
            auto area = dsc->area(fkey);
            double ci = g_param.mean_intensity[dsc->get_label(fkey)];
            
            double sum = image_->get_sum_on_tri_differ(tris, ci);
            //        double sum = image_->get_sum_on_tri_variation(tris, 3);
            
            if(sum < 0.001) sum = 0;
            
            auto center = (tris[0] + tris[1] + tris[2])/3.0;
            std::ostringstream is;
            is << sum/area;
            Painter::print_gl(center[0], center[1], is.str().c_str());
        }
    }

//    if (options_disp::get_option("Edge energy", false)) {
//        for (auto ekey : dsc->halfedges()){
//            auto hew = dsc->walker(ekey);
//            if (dsc->is_interface_dsc(ekey) and
//                hew.vertex().get_index() > hew.opp().vertex().get_index())
//            {
//                auto pts = dsc->get_pos(ekey);
//                double energy = std::abs(image_->get_edge_energy(pts[0], pts[1], 1));
//                auto c = (pts[0] + pts[1])/2;
//                std::ostringstream str;
//                str << energy;
//                Painter::print_gl(c[0], c[1], str.str().c_str());
//            }
//        }
//    }

}

std::vector<DSC2D::vec2> get_quad(double minx, double miny, double maxx, double maxy){
    std::vector<DSC2D::vec2> quad_tex;
    quad_tex.push_back(DSC2D::vec2(minx, maxy));
    quad_tex.push_back(DSC2D::vec2(maxx, maxy));
    quad_tex.push_back(DSC2D::vec2(maxx, miny));
    quad_tex.push_back(DSC2D::vec2(minx, miny));
    
    return quad_tex;
}

void interface_dsc::draw_coord(){
    DSC2D::DesignDomain const * domain = dsc->get_design_domain();
    std::vector<DSC2D::vec2> corners = domain->get_corners();
    DSC2D::vec2 min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for(auto p : corners){
        if (min[0] > p[0])
            min[0] = p[0];
        if (min[1] > p[1])
            min[1] = p[1];
        if (max[0] < p[0])
            max[0] = p[0];
        if (max[1] < p[1])
            max[1] = p[1];
    }
    
    DSC2D::vec2 center = (min + max)/2.0;
    DSC2D::vec2 length = (max - min);
    double l = std::min(length[0], length[1]);
    DSC2D::vec2 lx = center; lx[0] += l/2.0;
    DSC2D::vec2 ly = center; ly[1] += l/2.0;
    
    
    glLineWidth(2.0);
    glBegin(GL_LINES);
    
    glColor3f(1, 0, 0);
    glVertex2dv(center.get());
    glVertex2dv(lx.get());
    
    glColor3f(0, 1, 0);
    glVertex2dv(center.get());
    glVertex2dv(ly.get());
    
    glEnd();
    glLineWidth(1.0);
}

void interface_dsc::draw_image(){

    DSC2D::DesignDomain const * domain = dsc->get_design_domain();
    std::vector<DSC2D::vec2> corners = domain->get_corners();

    std::vector<DSC2D::vec2> quad_v = get_quad(0, 0, imageSize[0], imageSize[1]);
    std::vector<DSC2D::vec2> quad_tex;// = get_quad(0.0, 0.0, 1.0, 1.0);

    quad_tex.push_back(DSC2D::vec2(1, 0));
    quad_tex.push_back(DSC2D::vec2(1, 1));
    quad_tex.push_back(DSC2D::vec2(0, 1));
    quad_tex.push_back(DSC2D::vec2(0, 0));
    
    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);
    for (int i = 0; i < 4; i++) {
        glVertex2dv((GLdouble*)quad_v[i].get());
        glTexCoord2dv((GLdouble*)quad_tex[i].get());
    }
    glEnd();
    
    glLineWidth(2.0);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (int i = 0; i < 4; i++) {
        glVertex2dv((GLdouble*)quad_v[i].get());
        glVertex2dv((GLdouble*)quad_v[(i+1)%4].get());
    }
    glEnd();
}

void interface_dsc::update_title()
{
    std::ostringstream oss;
    oss << "2D DSC\t";
    if (RUN) {
        oss << " running iteration " << iter;
    }

    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

interface_dsc::interface_dsc(int &argc, char** argv){
    
    
    // 1. Parser arguments
    if (argc < 2)
    {
        std::cout << "Invalid setting file" << std::endl;
        exit(1412);
    }
    try
    {
        setting_io::set_param(argv[1]);
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
    }
    
    
    // 2.
    
    instance = this;
    WIN_SIZE_X = 900;
    WIN_SIZE_Y = 600;
    
    bDiplay_.resize(10);
    std::fill(bDiplay_.begin(), bDiplay_.end(), true);
    
    glutInit(&argc, argv);
    initGL();
    
    dyn_ = std::unique_ptr<dynamics_mul>(new dynamics_mul);
    dsc = nullptr;
    
    image_ = std::unique_ptr<image>(new image);
    image_->load_image(IMAGE_PATH);
    
    
    imageSize = Vec2(image_->width(), image_->height());
    // Boundary
//    imageSize = Vec2(image_->width() + 2*DISCRETIZE_RES, image_->height() + 2*DISCRETIZE_RES);

    check_gl_error();
    
    init_dsc();
    
//    threshold_initialization();
    
    gl_debug_helper::set_dsc(&(*dsc));
    
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
    display();
    
    
}

#pragma mark - Data
using namespace DSC2D;

bool is_boundary(DeformableSimplicialComplex & dsc, Face_key fkey){
    for(auto hew = dsc.walker(fkey); !hew.full_circle(); hew = hew.circulate_face_ccw()){
        if (HMesh::boundary(*dsc.mesh, hew.vertex())) {
            return true;
        }
    }
    return false;
}


void interface_dsc::init_dsc(){
    
    
    double width = imageSize[0];
    double height = imageSize[1];
    
    DISCRETIZATION = DISCRETIZE_RES;
    
    // Boundary gap
    width -= 2*DISCRETIZATION;
    height -= 2*DISCRETIZATION;
    
    std::vector<real> points;
    std::vector<int> faces;
    Trializer::trialize(width, height, DISCRETIZATION, points, faces);
    
    // Boundary gap
//    // Offset the mesh
//    for (auto & p:points)
//    {
//        p -= DISCRETIZE_RES;
//    }
//    //
    
    // Boundary gap
    width += 2*DISCRETIZATION;
    height += 2*DISCRETIZATION;
    
    DesignDomain *domain = new DesignDomain(DesignDomain::RECTANGLE, width, height, 0 /*,  DISCRETIZATION */);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex>(
                            new DeformableSimplicialComplex(DISCRETIZATION, points, faces, domain));
    
    if (ADAPTIVE == 1)
    {
        dsc->set_smallest_feature_size(SMALLEST_SIZE);
    }else
    {
        dsc->set_uniform_smallest_feature(SMALLEST_SIZE);
    }
    
#ifdef TUAN_MULTI_RES
    dsc->img = &*image_;
#endif
    
    // Boundary gap
//    // Label all margin triangle
//    for (auto fkey : dsc->faces())
//    {
//        if (is_boundary(*dsc, fkey))
//        {
//            dsc->update_attributes(fkey, 100);
//        }
//    }
//    dsc->clean_attributes();
    
    //
//    dsc->deform();
    // Initialize if need
//    manual_init_dsc();
    
    random_init_dsc(NB_PHASE);
    
    printf("Average edge length: %f ; # faces: %d\n", dsc->get_avg_edge_length(), dsc->get_no_faces());
}

void interface_dsc::export_dsc()
{
    dsc->clean_attributes();
    
    std::map<unsigned int, unsigned int> node_idx_map;
    std::vector<Vec2> vertices;
    int idx = 1;
    for (auto nkey : dsc->vertices())
    {
//        if (!HMesh::boundary(*dsc->mesh, nkey))
        {
            node_idx_map.insert(std::make_pair(nkey.get_index(), idx));
            vertices.push_back(dsc->get_pos(nkey));
            idx++;
        }
    }
    
    std::ofstream f("mesh.obj");
    if (f.is_open())
    {
        // vertices
        for (auto v:vertices)
        {
            f << "v " << v[0] << " " << v[1] << " 0" << endl;
        }
        
        // triangle
        for (auto fkey : dsc->faces())
        {
//            if (dsc->get_label(fkey) != BOUND_FACE)
            {
                auto verts = dsc->get_verts(fkey);
                f << "f " << node_idx_map[verts[0].get_index()] << " "
                << node_idx_map[verts[1].get_index()] << " "
                << node_idx_map[verts[2].get_index()] << endl;
            }
        }
        
        // label
        for (auto fkey : dsc->faces())
        {
//            if (dsc->get_label(fkey) != BOUND_FACE)
            {
                f << "l " << dsc->get_label(fkey) << endl;
            }
        }
        
        f.close();
    }
    else{
        cout<<"Fail to write file mesh.obj" << endl;
    }
    
}

void interface_dsc::random_init_dsc(int nb_phase)
{
    if (NB_PHASE <= 1)
    {
        return;
    }
    
    // Relabel
    for (auto tri : dsc->faces())
    {
        if (dsc->get_label(tri) != BOUND_FACE)
        {
            int idx = rand()%nb_phase;
            dsc->update_attributes(tri, idx);
        }
    }
    
    dsc->clean_attributes();
    
    cout << "Random initialization with " << NB_PHASE << " phases \n";
}

void interface_dsc::manual_init_dsc()
{
    
    // Relabel
    for (auto tri : dsc->faces())
    {
        auto pts = dsc->get_pos(tri);
        // check if it belong to any circle
        for (auto & cc : _circle_inits)
        {
            auto circles = cc.second;
            for (int j = 0; j < circles.size(); j++)
            {
                auto & circle = circles[j];
                if(circle.is_in_circle(pts)
                   && dsc->get_label(tri) != BOUND_FACE){
                    dsc->update_attributes(tri, cc.first); // Phase 0 is reserved for background
                }
            }
        }
    }
    
    // Deform
    for (auto v:dsc->vertices())
    {
        if (dsc->is_interface(v))
        {
            auto pt = dsc->get_pos(v);
            // check if it belongs to any circle
            for (auto & cc : _circle_inits)
            {
                auto circles = cc.second;
                for (int j = 0; j < circles.size(); j++)
                {
                    auto & circle = circles[j];
                    if(circle.is_in_circle(pt)){
                        auto pc = circle.project_to_circle(pt);
                        dsc->set_destination(v, pc);
                    }
                }
            }
        }
    }
    dsc->deform();
}

void interface_dsc::threshold_initialization()
{
    double c[2] = {0.6, 0.85};
    for (auto fkey : dsc->faces())
    {
        
        auto pts = dsc->get_pos(fkey);
        auto sumIntensity = image_->get_sum_on_tri_intensity(pts);
        double area = dsc->area(fkey);
        
        double average = sumIntensity / area;
        
        int new_label = 0;
        if (average > 0.8)
        {
            new_label = 0;
        }else
        {
            new_label = 1;
        }
        
        dsc->set_label(fkey, new_label);
    }
    
    dsc->clean_attributes();
}

void interface_dsc::thres_hold_init(){
    for (auto fid = dsc->faces_begin(); fid != dsc->faces_end(); fid++) {
        auto tris = dsc->get_pos(*fid);
        int totalPixel = 0;
        double total_inten = 0.0;
        image_->get_tri_intensity(tris, & totalPixel, & total_inten);
        double mean_c = total_inten / (double)totalPixel;
        
        if (mean_c < 0.2) {
            dsc->set_label(*fid, 0);
        }
        else if (mean_c < 0.7)
            dsc->set_label(*fid, 1);
        else
            dsc->set_label(*fid, 3);
    }
}

void interface_dsc::init_sqaure_boundary(){
    Vec2 s = imageSize;// tex->get_image_size();
    double left = 0.2;
    double right = 0.6;
    ObjectGenerator::create_square(*dsc, vec2(left*s[0], left*s[1]), vec2(right*s[0], right*s[1]), 1);
}

void interface_dsc::init_boundary(){
    std::vector<DSC2D::DeformableSimplicialComplex::face_key> faceKeys;
    for (auto p = dsc->faces_begin(); p != dsc->faces_end(); p++) {
        // Compute average intensity inside triangle
        auto pts = dsc->get_pos(*p);
        int count;
//        if(image_->get_triangle_intensity_count(pts, &count) > 0.1*count*MAX_BYTE ){
//            faceKeys.push_back(*p);
//        }
    }
    
    ObjectGenerator::label_tris(*dsc, faceKeys, 1);
}


void interface_dsc::dynamics_image_seg(){
    // Old approach
    // Edge-based force
    
    profile t("Total time");
    
    dyn_->update_dsc(*dsc, *image_);
    iter ++;
    // Virtual displacement
    // Compute energy change with assumed movement
//    static dyn_integral dyn;
//    dyn.update_dsc(*dsc, *image_);

//    static dynamics_edge dyn;
//    dyn.update_dsc(*dsc, *image_);
}

int closest(std::vector<double> & array, double v){
    
    int idx = -1;
    double smallest = INFINITY;
    
    for (int i = 0; i < array.size(); i++) {
        double dis = std::abs(v - array[i]);
        if (dis < smallest) {
            smallest = dis;
            idx = i;
        }
    }
    
    return idx;
}

void interface_dsc::init_boundary_brain(){
    // Threshold initialization
    std::vector<double> threshold = {0.0, 70./255., 77./255.};
    
    std::vector<vector<Face_key>> labels_list;
    labels_list.resize(threshold.size());
    
    for (auto fkey : dsc->faces()) {
        auto pts = dsc->get_pos(fkey);
        double area;
        double inten = image_->get_tri_intensity_f(pts, &area);
        
        double mean_i = inten / area;
        
        int idx = closest(threshold, mean_i);
        labels_list[idx].push_back(fkey);
    }
    
    
    for (int i = 0; i < labels_list.size(); i++) {
        ObjectGenerator::label_tris(*dsc, labels_list[i], i);
    }
}
