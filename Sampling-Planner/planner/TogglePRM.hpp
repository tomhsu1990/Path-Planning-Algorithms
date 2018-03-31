/*
 * TogglePRM.hpp
 *
 *  Created on: Mar. 30 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef TOGGLEPRM_H_
#define TOGGLEPRM_H_

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

class TogglePRM : public PRM {
public:

    TogglePRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
        unsigned int n_sample, unsigned int k_connection);
    virtual ~TogglePRM ();

    virtual bool findPath ();

protected:
    virtual void toggleConnect();

    int connect2Map (const Config& cfg, std::vector<int> &cc);
    bool findPathV1V2 (int v1, int v2, PATH& path);
};

#endif // TOGGLEPRM_H_
