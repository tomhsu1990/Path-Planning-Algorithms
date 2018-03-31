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
: Planner(env,tr,rr,start,goal), m_graph_free(NULL) {
    m_n_sample=n_sample;
    m_k_closest=k_connection;
    m_skip_same_cc=true;
    m_create_good_loops=false;
}

PRM::~PRM()
{}

void PRM::addPoint (Config cfg) {
    vertex_descriptor v = boost::add_vertex(*m_graph_free);
    (*m_graph_free)[v].idx = cfgs.size()-1;
    (*m_graph_free)[v].cfg = cfg;
    cfgs.push_back(v);
    flann::Matrix<double> flann_point(new double[m_env.dim], 1, m_env.dim);
    for (int i=0;i<m_env.dim;++i) {
        flann_point[0][i] = cfg.t[i];
    }
    if (index_ == NULL) {
        index_.reset(new flann::Index< flann::L2<double> >(
                         flann_point, flann::KDTreeIndexParams(1)));
        index_->buildIndex();
        return ;
    }
    index_->addPoints(flann_point, 2);
}

bool PRM::findPath () {
    if (!isValid(m_start) || !isValid(m_goal)) return false;

    sample();
    connect();

    std::vector<int> component(boost::num_vertices(*m_graph_free));
    int num = boost::connected_components(*m_graph_free, &component[0]);
    std::unordered_map<int, std::map<double, int>> CC_start;
    for (int i=0;i<component.size();++i) {
        Config cfg = (*m_graph_free)[cfgs[i]].cfg;
        CC_start[component[i]].insert({cfg.distance(m_start), i});
    }
    std::unordered_map<int, std::map<double, int>> CC_goal;
    for (int i=0;i<component.size();++i) {
        Config cfg = (*m_graph_free)[cfgs[i]].cfg;
        CC_goal[component[i]].insert({cfg.distance(m_goal), i});
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

void PRM::sample () {
    m_graph_free = new UndirectedGraph(0);
    for (int i=0;i<m_n_sample;++i) {
        Config cfg = Config::randomCfg(m_start.dim_t, m_start.dim_r);
        if (!isValid(cfg)) {
            continue;
        }
        addPoint(cfg);
    }
}

void PRM::connect () {
    //for each cfg find k closest
    for (int i=0;i<cfgs.size();++i) {
        Config cfg1 = (*m_graph_free)[cfgs[i]].cfg;
        // Convert the input point to the FLANN format.
        flann::Matrix<double> flann_query(new double[m_env.dim], 1, m_env.dim);
        for (int j=0;j<m_env.dim;++j)
            flann_query[0][j] = cfg1.t[j];

        // Search the kd tree for the nearest neighbor to the query.
        std::vector< std::vector<int> > query_match_indices;
        std::vector< std::vector<double> > query_distances;
        int num_neighbors_found = index_->knnSearch(
                                flann_query, query_match_indices, query_distances, m_k_closest,
                                flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED) /* no approx */);

        for (int j=0;j<num_neighbors_found;++j) {
            Config cfg2 = (*m_graph_free)[cfgs[query_match_indices[0][j]]].cfg;
            if (boost::edge(cfgs[i], cfgs[query_match_indices[0][j]], *m_graph_free).second) continue;
            if (isValid(cfg1,cfg2) || isValid(cfg2,cfg1)) {
                boost::add_edge(cfgs[i], cfgs[query_match_indices[0][j]], cfg1.distance(cfg2), *m_graph_free);
            }
        }
    }
}

//connect cfg to rmap
int PRM::connect2Map (const Config& cfg, std::vector<int> &cc) {
    for (int i=0;i<cc.size();++i) {
        Config cfg2 = (*m_graph_free)[cfgs[cc[i]]].cfg;
        if (isValid(cfg,cfg2) || isValid(cfg2,cfg))
            return cc[i];
    }
    return -1;
}

bool PRM::findPathV1V2 (int v1, int v2, PATH& path) {
    if (v1 == v2) return true;
    std::vector<int> predecessors(boost::num_vertices(*m_graph_free));
    std::vector<double> distances( boost::num_vertices(*m_graph_free));

    // Dijkstra's Shortest Paths
    boost::dijkstra_shortest_paths(*m_graph_free, v1,
                                   boost::predecessor_map(&predecessors[0]).distance_map(&distances[0]));
    if (predecessors.size() == 0) return false;

    int current = v2;
    while (current != v1) {
        path.push_back(toPhysical((*m_graph_free)[cfgs[current]].cfg));
        current = predecessors[current];
    }
    path.push_back(toPhysical((*m_graph_free)[cfgs[v1]].cfg));

    std::reverse(path.begin(), path.end());
    return true;
}
