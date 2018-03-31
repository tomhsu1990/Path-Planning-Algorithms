/*
 * Planner.hpp
 *
 *  Created on: Mar. 23 2018
 *      Author: Ching-Hsiang Hsu
 */

#ifndef PLANNER_H_
#define PLANNER_H_

#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Core>
#include "flann/flann.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

#include "config/Config.hpp"
#include "config/Robot.hpp"
#include "config/Env.hpp"
#include "utils/Timer.hpp"

// Create a struct to hold properties for each vertex
struct VertexProperties {
    Config cfg;
};
// Create a struct to hold properties for each edge
//struct EdgeProperties { double weight; };

typedef std::vector<Config> PATH;
typedef boost::property<boost::edge_weight_t, double> EdgeWeight;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, VertexProperties, EdgeWeight> UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<UndirectedGraph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<UndirectedGraph>::edge_iterator edge_iterator;

class Planner {
public:
    //workspace, translational resolution, rotational resolution...
    Planner(const Env &env, double tr, double rr, const Config& start, const Config& goal);

    virtual bool findPath () = 0;
    const PATH& getPath () const { return this->m_path; }

    void addPolygon (const Polygon2d& poly) { this->m_env.polys.push_back(poly); }
    std::vector<Polygon2d>& getPolygons () { return this->m_env.polys; }


    void setStart (const Config &cfg) {
        assert(cfg.ws);
        m_start = cfg;
    }
    void setGoal (const Config &cfg) {
        assert(cfg.ws);
        m_goal = cfg;
    }

    const Config& getStart () const { return this->m_start; }
    const Config& getGoal () const { return this->m_goal; }

    Robot& getRobot () { return this->m_robot; }

    Config toParametric (const Config& cfg);
    Config toPhysical (const Config& cfg);


protected:

    bool isValid(const Config& cfg);
    bool isValid(const Config& c1, const Config& c2);
    bool isValid (const Config& c1, const Config& c2, Config& c3);

    Config m_start;
    Config m_goal;

    Robot m_robot;
    Env m_env;

    bool m_found;
	
    double env_TR; // Translational RESOLUTION
    double env_RR; // Rotational RESOLUTION (deg)

	PATH m_path;
};

#endif // PLANNER_H_
