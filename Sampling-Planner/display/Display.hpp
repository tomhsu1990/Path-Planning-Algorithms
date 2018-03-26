

#ifndef Display_H_
#define Display_H_

// Qt & OpenGL
#include <QOpenGLWidget>
#include <QGLWidget>
#include <QtOpenGL>
#include <glu.h>

// Standard Library
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <ctime>

// Custom
#include "config/FileParameter.hpp"
#include "planner/RRT.hpp"
#include "planner/PRM.hpp"
#include "utils/Triangulate.hpp"

class Display : public QOpenGLWidget
{
    Q_OBJECT
public:
    Display(QWidget* parent = 0);
    virtual ~Display();

protected:
    // Essential Functions Inherited from QOpenGLWidget
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

private:
    void renderScene();

    /*********************************
     * Functions used to paint scene *
     *********************************/
    void drawPath(Planner* planner, const PATH& path);
    void drawRobot(Robot robot, const Config& cfg, std::vector<float> clr);

    void drawLink(Point2d a, Point2d b, std::vector<float> clr);
    void drawTriangle(Point2d a, Point2d b, Point2d c, std::vector<float> clr);
    void drawCircle(Robot robot, std::vector<float> clr);
    void drawPolygons(std::vector<Polygon2d>& objs);

    void drawTree(RRT* planner, const RRT_Tree& tree);
    void drawGraph(PRM* planner, const UndirectedGraph* graph, std::vector<float> vert_clr, std::vector<float> edge_clr);
};

#endif // Display_H_
