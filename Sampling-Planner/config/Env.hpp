/*
 *  Env.hpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#ifndef ENV_H_
#define ENV_H_

#include <vector>

#include "CORE/geom2d/point2d.h"
#include "CORE/geom2d/line2d.h"
#include "CORE/geom2d/segment2d.h"
#include "CORE/geom2d/polygon2d.h"

class Env {
public:
    int dim;
    std::vector<double> space;

    std::vector<Point2d> pts;
    std::vector<Segment2d> segs;
    std::vector<Polygon2d> polys;
};

#endif // ENV_H_
