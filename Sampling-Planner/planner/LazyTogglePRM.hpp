/*
 * LazyTogglePRM.hpp
 *
 *  Created on: Mar. 31 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef LAZYTOGGLEPRM_H_
#define LAZYTOGGLEPRM_H_

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>

#include "config/Config.hpp"
#include "LazyPRM.hpp"
#include "TogglePRM.hpp"

class LazyTogglePRM : public LazyPRM, public TogglePRM {
public:

    LazyTogglePRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
        unsigned int n_sample, unsigned int k_connection);
    virtual ~LazyTogglePRM ();

    virtual bool findPath ();

    UndirectedGraph& getObstGraph() { return m_graph[OBST]; }

protected:

    bool pathValidation (int v1, int v2, PATH& path, std::vector<Config> &witness);
    void witnessProcessing (std::vector<Config> &witness);

    unsigned int m_k_closest_free, m_k_closest_obst;
};

#endif // LAZYTOGGLEPRM_H_
