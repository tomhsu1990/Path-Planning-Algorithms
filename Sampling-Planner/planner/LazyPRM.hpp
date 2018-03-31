/*
 * LazyPRM.hpp
 *
 *  Created on: Mar. 30 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef LAZYPRM_H_
#define LAZYPRM_H_

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>

#include "config/Config.hpp"
#include "PRM.hpp"

class LazyPRM : public PRM {
public:

    LazyPRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
        unsigned int n_sample, unsigned int k_connection);
    virtual ~LazyPRM ();

    virtual bool findPath ();

    UndirectedGraph* getGraph() { return m_graph_free; }

protected:

    virtual void lazyConnect();

    int connect2Map (const Config& cfg, std::vector<int> &cc);
    bool findPathV1V2 (int v1, int v2, PATH& path);

    UndirectedGraph *m_graph_free;
    std::vector<vertex_descriptor> cfgs;
    std::shared_ptr< flann::Index< flann::L2<double> > > index_;

    unsigned int m_n_sample;
    unsigned int m_k_closest;

    bool m_skip_same_cc;
    bool m_create_good_loops;

private:
    void addPoint (Config cfg) ;
};

#endif // PRM_H_
