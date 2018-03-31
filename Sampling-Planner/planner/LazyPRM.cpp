/*
 *  LazyPRM.cpp
 *
 *  Created on: Mar. 30, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "LazyPRM.hpp"

LazyPRM::LazyPRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
         unsigned int n_sample, unsigned int k_connection)
: PRM(env,tr,rr,start,goal,n_sample,k_connection) {
}

LazyPRM::~LazyPRM()
{}

bool LazyPRM::findPath () {
    if (!isValid(m_start) || !isValid(m_goal)) return false;

    sample();
    lazyConnect();

    std::vector<int> component(boost::num_vertices(m_graph[FREE]));
    int num = boost::connected_components(m_graph[FREE], &component[0]);
    std::unordered_map<int, std::map<double, int>> CC_start, CC_goal;


    vertex_iterator u, v, nxt;
    boost::tie(u, v) = boost::vertices(m_graph[FREE]);
    for (nxt=u; nxt!=v; ++nxt) {
        Config cfg = m_graph[FREE][*nxt].cfg;
        CC_start[component[*nxt]].insert({cfg.distance(m_start), *nxt});
    }
    for (nxt=u; nxt!=v; ++nxt) {
        Config cfg = m_graph[FREE][*nxt].cfg;
        CC_goal[component[*nxt]].insert({cfg.distance(m_goal), *nxt});
    }

    for (auto it=CC_start.begin();it!=CC_start.end();) {
        int v1 = it->second.begin()->second;
        auto it2=CC_goal.find(it->first);
        int v2 = it2->second.begin()->second;

        PATH mypath;
        if (!findPathV1V2(v1,v2,mypath)) {
            if (it->second.size() == 0) ++it;
            continue;
        }

        m_path.push_back(toPhysical(m_start));
        m_path.insert(m_path.end(), mypath.begin(), mypath.end());
        m_path.push_back(toPhysical(m_goal));

        return m_found = true;
    }

    return false;
}

void LazyPRM::lazyConnect () {
    //for each cfg in free graph find k closest
    vertex_iterator u, v, nxt;
    boost::tie(u, v) = boost::vertices(m_graph[FREE]);
    for (nxt=u; nxt!=v; ++nxt) {
        Config cfg1 = m_graph[FREE][*nxt].cfg;
        // Convert the input point to the FLANN format.
        flann::Matrix<double> flann_query(new double[m_env.dim], 1, m_env.dim);
        for (int j=0;j<m_env.dim;++j)
            flann_query[0][j] = cfg1.t[j];
        // Search the kd tree for the nearest neighbor to the query.
        std::vector< std::vector<int> > query_match_indices;
        std::vector< std::vector<double> > query_distances;
        int num_neighbors_found = index_[FREE]->knnSearch(
                                flann_query, query_match_indices, query_distances, m_k_closest+1,
                                flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED) /* no approx */);
        for (int j=1;j<num_neighbors_found;++j) {
            if (boost::edge(*nxt, query_match_indices[0][j],m_graph[FREE]).second) continue;
            Config cfg2 = m_graph[FREE][query_match_indices[0][j]].cfg;
            boost::add_edge(*nxt, query_match_indices[0][j], EdgeProperties(cfg1.distance(cfg2),0), m_graph[FREE]);
        }
    }
}

bool LazyPRM::findPathV1V2 (int v1, int v2, PATH& path) {
    if (v1 == v2) return true;
    std::vector<vertex_descriptor> predecessors(boost::num_vertices(m_graph[FREE]));
    std::vector<double> distances(boost::num_vertices(m_graph[FREE]));

    // Dijkstra's Shortest Paths
    boost::dijkstra_shortest_paths(m_graph[FREE], v1,
                                   boost::weight_map(boost::get(&EdgeProperties::weight,m_graph[FREE]))
                                   .distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index,m_graph[FREE])))
                                   .predecessor_map(boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index,m_graph[FREE])))
                                   );
    if (predecessors.size() == 0) return false;

    int current = v2;
    while (current != v1) {
        path.push_back(toPhysical(m_graph[FREE][current].cfg));
        if (!isValid(m_graph[FREE][current].cfg, m_graph[FREE][predecessors[current]].cfg)) {
            boost::remove_edge(current, predecessors[current], m_graph[FREE]);
            return false;
        }
        else {
            std::pair<edge_descriptor, bool> edge_pair = boost::edge(current, predecessors[current], m_graph[FREE]);
            m_graph[FREE][edge_pair.first].status = 1;
        }
        current = predecessors[current];
    }
    path.push_back(toPhysical(m_graph[FREE][v1].cfg));

    std::reverse(path.begin(), path.end());
    return true;
}

