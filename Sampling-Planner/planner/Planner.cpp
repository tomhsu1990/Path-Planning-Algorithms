/*
 *  Planner.cpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "Planner.hpp"

Planner::Planner (const Env &env, double tr, double rr, const Config& start, const Config& goal) {
    m_env.dim = env.dim;
    for(int i=0;i<m_env.dim;++i)
        m_env.space.push_back(env.space[i]);
    env_TR=tr;
    env_RR=rr;

    m_start = start;
    m_goal = goal;

    m_path.clear();
    m_path_index = 0;

    if (env.dim == 2) {
        Polygon2d bb;
        //add a bounding box
        bb.addPoint(0, 0);
        bb.addPoint(0, env.space[0] + 0);
        bb.addPoint(env.space[1] + 0, env.space[0] + 0);
        bb.addPoint(env.space[1] + 0, 0);
        bb.addPoint(0, 0);
        bb.setOrientation();
        this->addPolygon(bb);
    }
}
Planner::~Planner () {fprintf(stderr, "in ~\n");}

Config Planner::toPhysical (const Config& cfg) {
    std::vector<double> ref;
    for (int i=0;i<m_env.dim;++i)
        ref.push_back(m_env.space[i]/2);
    return cfg.toPhysical(m_env.space, ref);
}

Config Planner::toParametric (const Config& cfg) {
    std::vector<double> ref;
    for (int i=0;i<m_env.dim;++i)
        ref.push_back(m_env.space[i]/2);
    return cfg.toParamtric(m_env.space, ref);
}

//check whether a cfg is valid or not
bool Planner::isValid (const Config& cfg) {
    Config phy_cfg=toPhysical(cfg);
    m_robot.setConfig(phy_cfg);

    std::vector<double> robot_ang;
    std::vector<Point2d> robot_joint;
    Point2d robot_base(m_robot.cfg.t[0], m_robot.cfg.t[1]);
    for (int i=0;i<m_robot.cfg.dim_r;++i)
        robot_ang.push_back(m_robot.cfg.r[i]*M_PI/180.0f);
    robot_joint.push_back(robot_base);
    for (unsigned i=1;i<=m_robot.L.size();++i)
        robot_joint.push_back(Point2d(robot_joint[i-1].X()+m_robot.L[i-1]*cos(robot_ang[i-1]),
                                      robot_joint[i-1].Y()+m_robot.L[i-1]*sin(robot_ang[i-1])));

    // can use flann to speed up the search, ex: radius search
    for (int i=0;i<m_env.polys.size();++i) {
        Polygon2d poly(m_env.polys[i]);
        //different robot has different collision detection rules
        // disc
        if (m_robot.name.compare("disc") == 0) {

            //check cfg base
            if (poly.orientation() > 0 && !poly.inside(robot_base)) return false;
            if (poly.orientation() < 0 && poly.inside(robot_base))  return false;

            std::vector<Point2d> pts = poly.points();
            for (int j=0;j<pts.size()-1;++j) {
                Segment2d poly_seg(pts[j], pts[j+1]);
                if (poly_seg.distance(robot_base) <= m_robot.R) {
                    return false;
                }
            }
        }
        // link
        if (m_robot.name.compare("link") == 0) {

            //check cfg base
            for (unsigned j=0;j<robot_joint.size();++j) {
                if (poly.orientation() > 0 && !poly.inside(robot_joint[j])) return false;
                if (poly.orientation() < 0 && poly.inside(robot_joint[j]))  return false;
            }

            std::vector<Point2d> pts = poly.points();
            for (int j=0;j<pts.size()-1;++j) {
                Segment2d poly_seg(pts[j], pts[j+1]);
                Line2d poly_line(poly_seg.toLine());
                for (unsigned k=0;k<robot_joint.size()-1;++k) {
                    Segment2d robot_seg(robot_joint[k], robot_joint[k+1]);
                    Line2d robot_line(robot_seg.toLine());
                    // polygon's segment has an overlap with a robot link
                    if (poly_seg.contains(robot_seg.startPt()) || poly_seg.contains(robot_seg.stopPt()) ||
                        robot_seg.contains(poly_seg.startPt()) || robot_seg.contains(poly_seg.stopPt())) {
                        return false;
                    }
                    // check the distance between polygon's segment and the robot link
                    if (poly_seg.distance(robot_seg.startPt()) <= m_robot.Thickness*0.5f ||
                        poly_seg.distance(robot_seg.stopPt()) <= m_robot.Thickness*0.5f){
                        return false;
                    }
                    // check if they have an intersection
                    if (poly_line.intersects(robot_line) == 0) {
                        Point2d *intersection_pt = (Point2d *)poly_line.intersection(robot_line);
                        if (poly_seg.distance(*intersection_pt) <= m_robot.Thickness*0.5f &&
                            robot_seg.distance(*intersection_pt) <= m_robot.Thickness*0.5f){
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

//connect from c1 to c2
bool Planner::isValid (const Config& c1, const Config& c2) {
    Config dir = (c2 - c1); //general moving direction
    int steps = (int) ceil( std::max((dir.normT()/env_TR), (dir.normR()/env_RR)) );
    Config step = dir/steps;
    // this linear check needs to be improved !!!
    // modify it into binary check...
    for (double i=1;i<=steps;++i) {
        Config now_cfg = c1 + step*i;
        if (!this->isValid(now_cfg)) {
            return false;
        }
    }
    return true;
}

//connect from c1 to c2 and record c3 as a witness
bool Planner::isValid (const Config& c1, const Config& c2, Config& c3) {
    Config dir = (c2 - c1); //general moving direction
    int steps = (int) ceil( std::max((dir.normT()/env_TR), (dir.normR()/env_RR)) );
    Config step = dir/steps;

    for (int i=1;i<=steps;++i) {
        Config now_cfg = c1 + step*i;
        if (!this->isValid(now_cfg)) {
            c3 = c1 + step*(i-1);
            return false;
        }
    }
    return true;
}
