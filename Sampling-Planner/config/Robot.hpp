/*
 *  Robot.hpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#ifndef ROBOT_H_
#define ROBOT_H_

#include "Config.hpp"

class Robot {
public:
    Robot(){
        R=Thickness=0;
        L.clear();
    }
    ~Robot(){}

    void init (Robot r, const Config &c) {
        (*this) = r;
        this->setConfig(c);
    }
    void init (double r, const Config& c) {
        R=r;
        this->setConfig(c);
    }
    void init (std::vector<double> l, double thickness, const Config& c) {
        for (unsigned i=0;i<l.size();++i)
            L.push_back(l[i]);
        Thickness = thickness;
        this->setConfig(c);
    }

    void setConfig (const Config& c) {
        this->cfg = c;
    }

    std::string name;
    // disc
    double R;
    // link
    std::vector<double> L;
    double Thickness;
    // triangle
    // polygon

    Config cfg;
    std::vector<Pose> p;
};

#endif // ROBOT_H_
