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

protected:

    virtual void lazyConnect();

    bool findPathV1V2 (int v1, int v2, PATH& path);
};

#endif // PRM_H_
