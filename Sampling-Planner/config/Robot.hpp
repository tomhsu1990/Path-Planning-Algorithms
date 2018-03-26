

#ifndef ROBOT_H_
#define ROBOT_H_

#include "Config.hpp"

class Robot {
public:
    Robot(){
        R=L1=L2=Thickness=0;
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
    void init (double l1, double thickness, const Config& c) {
        L1=l1;
        Thickness = thickness;
        this->setConfig(c);
    }
    void init (double l1, double l2, double thickness, const Config& c) {
        L1=l1;
        L2=l2;
        Thickness = thickness;
        this->setConfig(c);
    }
//    void init (const Config& c) {
//        this->setConfig(c);
//    }
//    void init (double thickness, const Config& c) {
//        Thickness = thickness;
//        this->setConfig(c);
//    }

    void setConfig (const Config& c) {
        this->cfg = c;
    }

    std::string name;
    // disc
    double R;
    // rod
    double L1, L2;
    double Thickness;
    // triangle
    // polygon

    Config cfg;
    std::vector<Pose> p;
};

#endif // ROBOT_H_
