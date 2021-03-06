/*
 *  LazyTogglePRM.cpp
 *
 *  Created on: Mar. 31, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "LazyTogglePRM.hpp"

LazyTogglePRM::LazyTogglePRM(const Env &env, double tr, double rr, const Config& start, const Config& goal,
         unsigned int n_sample, unsigned int k_connection)
: LazyPRM(env,tr,rr,start,goal,n_sample,k_connection),
  TogglePRM(env,tr,rr,start,goal,n_sample,k_connection),
  PRM(env,tr,rr,start,goal,n_sample,k_connection) {
}

LazyTogglePRM::~LazyTogglePRM()
{}

bool LazyTogglePRM::findPath () {
    if (!isValid(m_start) || !isValid(m_goal)) return false;

    sample();

    while (true) {
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

        // path validation
        std::vector<Config> witness;
        for (auto it=CC_start.begin();it!=CC_start.end();) {
            int v1 = it->second.begin()->second;
            auto it2=CC_goal.find(it->first);
            int v2 = it2->second.begin()->second;

            PATH mypath;
            if (!pathValidation(v1,v2,mypath,witness)) {
                if (it->second.size() == 0) ++it;
                continue;
            }

            m_path.push_back(toPhysical(m_start));
            m_path.insert(m_path.end(), mypath.begin(), mypath.end());
            m_path.push_back(toPhysical(m_goal));

            return m_found = true;
        }

        witnessProcessing(witness);
    }

    return false;
}

bool LazyTogglePRM::pathValidation (int v1, int v2, PATH& path, std::vector<Config> &witness) {
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

    int cnt(0), current = v2;
    while (current != v1) {
        path.push_back(toPhysical(m_graph[FREE][current].cfg));
        Config cfg;
        if (!isValid(m_graph[FREE][current].cfg, m_graph[FREE][predecessors[current]].cfg, cfg)) {
            boost::remove_edge(current, predecessors[current], m_graph[FREE]);
            witness.push_back(cfg);
            ++cnt;
        }
        else {
            std::pair<edge_descriptor, bool> edge_pair = boost::edge(current, predecessors[current], m_graph[FREE]);
            m_graph[FREE][edge_pair.first].status = 1;
        }
        current = predecessors[current];
    }
    path.push_back(toPhysical(m_graph[FREE][v1].cfg));

    std::reverse(path.begin(), path.end());
    return cnt == 0;
}

void LazyTogglePRM::witnessProcessing (std::vector<Config> &witness) {
    for (unsigned i=0;i<witness.size();++i)
        addPoint(witness[i], m_graph[FREE], index_[FREE]);
}
