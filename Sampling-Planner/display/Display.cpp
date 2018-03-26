

#include "Display.hpp"

std::vector<std::vector<float>> color_table;
//bool do_triangulation(false);
//std::vector<Point2d> triangles;

extern FileParameter *f_prm;
extern Planner *planner;

/*
 * CONSTRUCTOR
 *
 * Define data members
 */
Display::Display(QWidget* parent):
    QOpenGLWidget(parent) {}

/*
 * DESTRUCTOR
 *
 * Destroy shader program
 */
Display::~Display()
{}

void Display::initializeGL () {
    glClearColor(0.5f, 0.5f, 0.5f, 0.5f);
    glClearDepth(1.0f);

    color_table.push_back({0,0,0}); // black
    color_table.push_back({1,1,1}); // white
    color_table.push_back({1,0,0}); // red
    color_table.push_back({0,1,0}); // green
    color_table.push_back({0,0,1}); // blue
    color_table.push_back({1,1,0}); // yellow
    color_table.push_back({0,1,1}); // cyan
    color_table.push_back({1,0,1}); // magenta
    for(int i=0;i<1000;++i)
        color_table.push_back({(float)drand48(), (float)drand48(), (float)drand48()});
}


/*
 * PAINT GL
 *
 * Regenerates (if necessary) and draws scene to
 * dispaly screen.
 */
void Display::paintGL () {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    renderScene();
    glFlush();
}


/*
 * RESIZE GL
 *
 * Sets size of viewport.
 * If 'Display' is the wooden frame of a painting,
 * the viewport is the canvas.
 */
void Display::resizeGL (int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height); // set origin to bottom left corner
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void Display::renderScene () {

    drawPolygons(planner->getPolygons());

    if(dynamic_cast<PRM*>(planner) != NULL && f_prm->show_prm_graph) {
        PRM * prm=dynamic_cast<PRM*>(planner);
        drawGraph(prm, prm->getGraph(), color_table[3], color_table[2]);
        update();
    }
    if (dynamic_cast<RRT*>(planner) != NULL && f_prm->show_rrt_graph) {
        RRT * rrt=dynamic_cast<RRT*>(planner);
        drawTree(rrt, rrt->getTree());
        update();
    }

    if (!f_prm->no_path) {
        drawPath(planner, planner->getPath());

//        if (f_prm->replay_anim) {
//            f_prm->replay_anim = false;

//            f_prm->show_anim = true;
//            f_prm->pause_anim = false;
//        }

//        if(f_prm->show_anim){
//            usleep((99-animationSpeed)*animationSpeedScale);

//            const PATH& path = planner->getPath();

//            if(path_index < 0) path_index = 0;
//            if(path_index >= path.size()) path_index = path.size()-1;

//            if (f_prm->pause_anim) {
//                drawRobot(planner->getRobot(), planner->to_physical(path[path_index]));
//            }
//            else {
//                if(path_index < planner->getPath().size()-1){
//                    drawRobot(planner->getRobot(), planner->to_physical(path[path_index]));
//                    path_index++;
//                } else if(path_index+1 == planner->getPath().size()){
//                    drawRobot(planner->getRobot(), planner->to_physical(path[path_index]));

//                    showAnim=false;
//                    path_index = 0;
//                }
//            }
//            update();
//        }
//        if (f_prm->show_trace) {
//            const PATH& path = planner->getPath();
//            for(int i=0;i<planner->getPath().size();++i){
//                drawRobot(planner->getRobot(), planner->toPhysical(path[i]));
//            }
//        }
    }

    drawRobot(planner->getRobot(), planner->toPhysical(planner->getStart()), color_table[6]);
    drawRobot(planner->getRobot(), planner->toPhysical(planner->getGoal()), color_table[7]);
}

void Display::drawPath (Planner* planner, const PATH& path) {
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(3);
    glBegin(GL_LINE_STRIP);
    for (int i=0;i<path.size();++i) {
        if (path[i].dim_t == 2) glVertex2d(path[i].t[0], path[i].t[1]);
        if (path[i].dim_t == 3) glVertex3d(path[i].t[0], path[i].t[1], path[i].t[2]);
    }
    glEnd();
}

void Display::drawRobot (Robot robot, const Config& cfg, std::vector<float> clr) {
    robot.setConfig(cfg);
    if (robot.name.compare("disc") == 0)
        drawCircle(robot, clr);
}

void Display::drawLink (Point2d a, Point2d b, std::vector<float> clr) {
//    glColor3fv(clr.rgb);
//    glLineWidth(thickness*2.0f);
//    glBegin(GL_LINE_STRIP);
//        glVertex2f(a.x, a.y);
//        glVertex2f(b.x, b.y);
//    glEnd();
//    glLineWidth(1);

//    // remember to do the exchange of the vector
//    double vec_x = b[0]-a[0];
//    double vec_y = b[1]-a[1];
//    double norm2 = sqrt(vec_x*vec_x+vec_y*vec_y);
//    vec_x = -(b[1]-a[1])/norm2;
//    vec_y = (b[0]-a[0])/norm2;

//    Vector2d p[4];

//    int dx[4] = {1,-1,1,-1};
//    int dy[4] = {1,-1,1,-1};

//    thickness *= 0.5;
//    for(int i=0;i<2;++i){
//        p[i]  = Vector2d(a[0]+dx[i]*thickness*vec_x, a[1]+dy[i]*thickness*vec_y);
//    }
//    for(int i=2;i<4;++i){
//        p[i]  = Vector2d(b[0]+dx[i]*thickness*vec_x, b[1]+dy[i]*thickness*vec_y);
//    }
//    thickness *= 2;


//    qsort(p, 4, sizeof(Vector2d), cmp);

//    if(p[0][1] > p[3][1]){
//        drawTriangle(p[0], p[3], p[1], clr);
//        drawTriangle(p[0], p[2], p[3], clr);
//    }
//    else{
//        drawTriangle(p[0], p[3], p[2], clr);
//        drawTriangle(p[0], p[1], p[3], clr);
//    }
}
void Display::drawTriangle (Point2d a, Point2d b, Point2d c, std::vector<float> clr) {
    glColor3f(clr[0], clr[1], clr[2]);
    glBegin(GL_TRIANGLES);
        glVertex2f( a.X(), a.Y() );
        glVertex2f( b.X(), b.Y() );
        glVertex2f( c.X(), c.Y() );
    glEnd();
}
void Display::drawCircle(Robot robot, std::vector<float> clr){
    glColor3f(clr[0], clr[1], clr[2]);
    glBegin(GL_TRIANGLE_FAN);
    for (int ii = 0; ii < 360; ii++) {
        GLfloat theta = M_PI*GLfloat(ii)/180.0f;
        GLfloat x = robot.R*0.5f*cosf(theta);
        GLfloat y = robot.R*0.5f*sinf(theta);
        glVertex2f(x + robot.cfg.t[0], y + robot.cfg.t[1]);
    }
    glEnd();
}

void Display::drawPolygons(std::vector<Polygon2d>& objs) {

    glLineWidth(2);
    for (unsigned i=0;i<objs.size();++i) {
        glColor3f(0, 0, 1);
        glBegin(GL_LINE_LOOP);
        std::vector<Point2d> pts = objs[i].points();
        for (unsigned j=0;j<pts.size();++j) {
            glVertex2d(pts[j].X(), pts[j].Y());
        }
        glEnd();
    }

//    if(showFilledObstacles){
//        if(!doTriangulation){
//            doTriangulation = true;
//            triangles.clear();
//            for(RIT it = objs.begin(); it != objs.end(); ++it) {
//                Vector2dVector a;
//                size_t size = it->getSize();
//                int bb = -1;
//                if((int)size == 4){
//                    ply_vertex* pp1 = it->operator [](0);
//                    const double *p1 = pp1->getPos().get();
//                    ply_vertex* pp2 = it->operator [](1);
//                    const double *p2 = pp2->getPos().get();
//                    ply_vertex* pp3 = it->operator [](2);
//                    const double *p3 = pp3->getPos().get();
//                    ply_vertex* pp4 = it->operator [](3);
//                    const double *p4 = pp4->getPos().get();
//                    if((int)p1[0] == 0 && (int)p1[1] == 0 &&
//                       (int)p2[0] == 0 && (int)p2[1] == 512 &&
//                       (int)p3[0] == 512 && (int)p3[1] == 512 &&
//                       (int)p4[0] == 512 && (int)p4[1] == 0){
//                        bb = 1;
//                    }
//                }
//                if(bb > 0) continue;
//                for(size_t i=0;i<size-1;i++) {
//                    ply_vertex* p = it->operator [](i);
//                    const double *pp = p->getPos().get();
//                    a.push_back(m_Vector2d(pp[0], pp[1]));
//                }
//                // allocate an STL vector to hold the answer.
//                Vector2dVector result;
//                //  Invoke the triangulator to triangulate this polygon.
//                m_Triangulate::Process(a,result);
//                // print out the results.
//                int tcount = result.size()/3;

//                for (int i=0; i<tcount; i++){
//                    const m_Vector2d &p1 = result[i*3+0];
//                    const m_Vector2d &p2 = result[i*3+1];
//                    const m_Vector2d &p3 = result[i*3+2];
//                    triangles.push_back(p1);
//                    triangles.push_back(p2);
//                    triangles.push_back(p3);
//                }
//            }
//        }
//        else{
//            for(int tri=0;tri<(int)triangles.size();tri+=3){
//                glColor3f(0.8,0.8,0.8);
//                glBegin(GL_TRIANGLES);
//                    glVertex2d(triangles[tri].GetX(), triangles[tri].GetY());
//                    glVertex2d(triangles[tri+1].GetX(), triangles[tri+1].GetY());
//                    glVertex2d(triangles[tri+2].GetX(), triangles[tri+2].GetY());
//                glEnd();
//            }
//        }
//    }
}

void Display::drawTree(RRT* planner, const RRT_Tree& tree){
    glLineWidth(0.5);
    glColor3f(0,0,0);
    glBegin(GL_LINES);
    for(unsigned i=0;i<tree.size();++i){
        const RRT_NODE* node = tree[i];
        const RRT_NODE* parent = tree[i]->parent;
        if(!parent) continue;

        Config ncfg=planner->toPhysical(node->cfg);
        Config pcfg=planner->toPhysical(parent->cfg);
        glVertex2d(ncfg.t[0], ncfg.t[1]);
        glVertex2d(pcfg.t[0], pcfg.t[1]);

//        glVertex2d(node->cfg.t[0], node->cfg.t[1]);
//        glVertex2d(parent->cfg.t[0], parent->cfg.t[1]);
    }
    glEnd();
}

void Display::drawGraph(PRM* planner, const UndirectedGraph* graph, std::vector<float> vert_clr, std::vector<float> edge_clr) {

    glDisable(GL_LIGHTING);

    glLineWidth(0.5);
    glColor3f(vert_clr[0], vert_clr[1], vert_clr[2]);
    ///////////////////////////////////////////////////////////////////////////
    //draw nodes
    glPointSize(2);
    glBegin(GL_POINTS);
    vertex_iterator u, v, nxt;
    boost::tie(u, v) = boost::vertices(*graph);
    for (nxt=u; nxt!=v; ++nxt) {
        Config cfg = planner->toPhysical((*graph)[*nxt].cfg);
        glVertex2d(cfg.t[0], cfg.t[1]);
    }
    glEnd();

    ///////////////////////////////////////////////////////////////////////////
    //draw edges
    std::pair<edge_iterator, edge_iterator> ei = boost::edges(*graph);
    glColor3f(edge_clr[0], edge_clr[1], edge_clr[2]);
    glPushAttrib(GL_CURRENT_BIT);
    glBegin(GL_LINES);
    for (edge_iterator it=ei.first; it!=ei.second; ++it) {
        unsigned src = boost::source(*it,*graph);
        unsigned tgt = boost::target(*it,*graph);
        Config cfg1 = planner->toPhysical((*graph)[src].cfg);
        Config cfg2 = planner->toPhysical((*graph)[tgt].cfg);
        glVertex2d(cfg1.t[0], cfg1.t[1]);
        glVertex2d(cfg2.t[0], cfg2.t[1]);
    }
    glEnd();
}