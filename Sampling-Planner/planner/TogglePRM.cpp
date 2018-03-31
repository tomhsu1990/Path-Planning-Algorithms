/*
 *  TogglePRM.cpp
 *
 *  Created on: Mar. 23, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "TogglePRM.hpp"

TogglePRM::TogglePRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
         unsigned int n_sample, unsigned int k_connection)
: PRM(env,tr,rr,start,goal, n_sample, k_connection) {
}

TogglePRM::~TogglePRM()
{}

bool TogglePRM::findPath () {
    if (!isValid(m_start) || !isValid(m_goal)) return false;

    sample();
    toggleConnect();

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

    for (auto it=CC_start.begin();it!=CC_start.end();++it) {
        std::vector<int> cc_start;
        for (auto len=it->second.begin();len!=it->second.end()&&cc_start.size()<m_k_closest;++len) {
            cc_start.push_back(len->second);
        }
        int v1 = connect2Map(m_start, cc_start);
        if (v1<0) continue;

        auto it2=CC_goal.find(it->first);
        std::vector<int> cc_goal;
        for (auto len=it2->second.begin();len!=it2->second.end()&&cc_goal.size()<m_k_closest;++len) {
            cc_goal.push_back(len->second);
        }
        int v2 = connect2Map(m_goal, cc_goal);
        if (v2<0) continue;

        PATH mypath;
        if (!findPathV1V2(v1,v2,mypath)) continue;

        m_path.push_back(toPhysical(m_start));
        m_path.insert(m_path.end(), mypath.begin(), mypath.end());
        m_path.push_back(toPhysical(m_goal));

        return m_found = true;
    }

    return false;
}

void TogglePRM::toggleConnect () {
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
                                flann_query, query_match_indices, query_distances, m_k_closest,
                                flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED) /* no approx */);
        for (int j=0;j<num_neighbors_found;++j) {
            Config cfg2 = m_graph[FREE][query_match_indices[0][j]].cfg;
            if (boost::edge(*nxt, query_match_indices[0][j], m_graph[FREE]).second) continue;
            if (isValid(cfg1,cfg2) || isValid(cfg2,cfg1)) {
                boost::add_edge(*nxt, query_match_indices[0][j], cfg1.distance(cfg2), m_graph[FREE]);
            }
        }
    }
}

//connect cfg to rmap
int TogglePRM::connect2Map (const Config& cfg, std::vector<int> &cc) {
    for (int i=0;i<cc.size();++i) {
        Config cfg2 = m_graph[FREE][cc[i]].cfg;
        if (isValid(cfg,cfg2) || isValid(cfg2,cfg))
            return cc[i];
    }
    return -1;
}

bool TogglePRM::findPathV1V2 (int v1, int v2, PATH& path) {
    if (v1 == v2) return true;
    std::vector<int> predecessors(boost::num_vertices(m_graph[FREE]));
    std::vector<double> distances(boost::num_vertices(m_graph[FREE]));

    // Dijkstra's Shortest Paths
    boost::dijkstra_shortest_paths(m_graph[FREE], v1,
                                   boost::predecessor_map(&predecessors[0]).distance_map(&distances[0]));
    if (predecessors.size() == 0) return false;

    int current = v2;
    while (current != v1) {
        path.push_back(toPhysical(m_graph[FREE][current].cfg));
        current = predecessors[current];
    }
    path.push_back(toPhysical(m_graph[FREE][v1].cfg));

    std::reverse(path.begin(), path.end());
    return true;
}
