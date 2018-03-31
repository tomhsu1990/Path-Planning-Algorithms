/*
 * PRM.hpp
 *
 *  Created on: Mar. 23 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef PRM_H_
#define PRM_H_

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>

#include "config/Config.hpp"
#include "Planner.hpp"

enum GRAPH_TYPE {
    FREE=0, OBST
};

class PRM : public Planner {
public:

    PRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
        unsigned int n_sample, unsigned int k_connection);
    virtual ~PRM ();

    virtual bool findPath ();

    UndirectedGraph& getGraph() { return m_graph[FREE]; }

protected:

    virtual void sample();
    virtual void connect();

    int connect2Map (const Config& cfg, std::vector<int> &cc);
    bool findPathV1V2 (int v1, int v2, PATH& path);

    UndirectedGraph m_graph[2];
    std::shared_ptr<flann::Index< flann::L2<double>>> index_[2];

    unsigned int m_n_sample;
    unsigned int m_k_closest;

private:
    void clean () {
        for (int i=0;i<2;++i) {
            index_[i].reset();
            index_[i] = nullptr;
            m_graph[i].clear();
        }
    }
    void addPoint (Config cfg, UndirectedGraph &m_graph, std::shared_ptr< flann::Index< flann::L2<double> > > &index_);
};

#endif // PRM_H_
