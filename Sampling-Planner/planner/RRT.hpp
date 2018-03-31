/*
 * RRT.hpp
 *
 *  Created on: Mar. 23 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef RRT_H_
#define RRT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

#include "config/Config.hpp"
#include "Planner.hpp"
#include "utils/Timer.hpp"

struct RRT_NODE {
    Config cfg;
    RRT_NODE* parent;
    RRT_NODE(const Config& c, RRT_NODE* parent = nullptr) {
        this->cfg = c;
        this->parent = parent;
    }
};

typedef std::vector<RRT_NODE*> RRT_Tree;

class RRT : public Planner {
public:

    RRT(const Env &env, double tr, double rr, const Config& start, const Config& goal,
        unsigned int max_sample, float expand_step, float bias, float close_to_goal);
    virtual ~RRT();

    virtual bool findPath();
    RRT_Tree& getTree() { return m_tree; }

protected:
    bool moveTo (const Config& nearest_cfg, const Config& rand_cfg, Config& out_cfg);
    RRT_NODE* nearest (const Config& cfg);

    RRT_Tree m_tree;
    std::shared_ptr< flann::Index< flann::L2<double> > > index_;

    unsigned int m_max_sample;
    float m_expand_step;
    float m_goal_bias;
    float m_close_to_goal;

private:
    void clean () {
        if (index_ != nullptr) {
            for (size_t ii = 0; ii < index_->size(); ++ii) {
                double* point = index_->getPoint(ii);
                delete[] point;
            }
        }
    }

    void addPoint (RRT_NODE* node) ;
};

#endif // RRT_H_
