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

    // can use flann to speed up the search, ex: radius search
    for (int i=0;i<m_env.polys.size();++i) {
        Polygon2d poly(m_env.polys[i]);
        Point2d robot_base(m_robot.cfg.t[0], m_robot.cfg.t[1]);
        //check cfg base
        if (poly.orientation() > 0) {
            if (!poly.inside(robot_base)) {
                return false;
            }
        }
        else {
            if (poly.inside(robot_base)) {
                return false;
            }
        }

        //different robot has different collision detection rules
        // disc
        if (m_robot.name.compare("disc") == 0) {
            std::vector<Point2d> pts = poly.points();
            for (int j=0;j<pts.size()-1;++j) {
                Segment2d seg(pts[j], pts[j+1]);
                if (seg.distance(robot_base) <= m_robot.R) {
                    return false;
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
