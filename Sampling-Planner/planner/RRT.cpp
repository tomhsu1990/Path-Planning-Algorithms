/*
 * RRT.cpp
 *
 *  Created on: Mar. 23 2018
 *      Author: Ching-Hsiang Hsu
 */

#include "RRT.hpp"

RRT::RRT (const Env &env, double tr, double rr, const Config& start, const Config& goal,
         unsigned int max_sample, float expand_step, float bias, float close_to_goal)
    : Planner(env,tr,rr,start,goal), index_(NULL) {
    m_max_sample=max_sample;
    m_expand_step=expand_step;
    m_goal_bias=bias;
    m_close_to_goal=close_to_goal;
}

RRT::~RRT () {
    if (index_ != nullptr) {
        for (size_t ii = 0; ii < index_->size(); ++ii) {
            double* point = index_->getPoint(ii);
            delete[] point;
        }
    }
}

void RRT::addPoint (RRT_NODE* node) {
    m_tree.push_back(node);
    flann::Matrix<double> flann_point(new double[m_env.dim], 1, m_env.dim);
    for (int i=0;i<m_env.dim;++i) {
        flann_point[0][i] = node->cfg.t[i];
    }
    if (index_ == NULL) {
        index_.reset(new flann::Index< flann::L2<double> >(
                         flann_point, flann::KDTreeIndexParams(1)));
        index_->buildIndex();
        return ;
    }
    index_->addPoints(flann_point, 2);
}

bool RRT::findPath () {
    if (!isValid(m_start) || !isValid(m_goal)) return false;

    RRT_NODE* start_node = new RRT_NODE(m_start);
    addPoint(start_node);

    double min_dist = (m_start-m_goal).norm();

	this->m_found = false;

    for (int i=0;i<m_max_sample && !m_found;i++) {
        Config rand_cfg;

        //bias to the goal...
        if (drand48()<this->m_goal_bias)
            rand_cfg=m_goal;
        else rand_cfg=Config::randomCfg(this->m_start.dim_t, this->m_start.dim_r);

        RRT_NODE* nearest_node = this->nearest(rand_cfg);
        Config out_cfg;

        if (this->moveTo(nearest_node->cfg, rand_cfg, out_cfg)) {
            RRT_NODE* new_node = new RRT_NODE(out_cfg, nearest_node);
            addPoint(new_node);

            double dist = (out_cfg-m_goal).norm();
            if (dist < min_dist) {
                min_dist = dist;
            }
            if (dist <= m_close_to_goal) {
                RRT_NODE* goal_node = new RRT_NODE(m_goal, new_node);
                addPoint(goal_node);

                m_found = true;
            }
        }
    }

    if (m_found) {
        RRT_NODE* node = this->m_tree.back();
        while (node) {
            m_path.push_back(toPhysical(node->cfg));
            node = node->parent;
        }
        std::reverse(m_path.begin(), m_path.end());
    }

	return m_found;
}

bool RRT::moveTo (const Config& nearest_cfg, const Config& rand_cfg, Config& out_cfg) {
    Config dir = (rand_cfg - nearest_cfg); //general moving direction
    double dist = dir.norm();              // distance form nearest to rand
    double exp_dist=std::min(dist,(double)m_expand_step); // move in dir for dist
    dir=dir.normalize()*exp_dist;

    int steps = (int) ceil( std::max((dir.normT()/env_TR), (dir.normR()/env_RR)) );
    Config step = dir/steps;

    for (int i=1;i<=steps;i++) {
        Config now_cfg = nearest_cfg + step*i;

        if (!this->isValid(now_cfg)) {
		    return false;
	    }

        if (i==steps || ((now_cfg-m_goal).normT() < env_TR && (now_cfg-m_goal).normR() < env_RR)) {
		    out_cfg = now_cfg;
		    return true;
	    }
    }

    return false;
}

RRT_NODE* RRT::nearest (const Config& cfg) {
    // Convert the input point to the FLANN format.
    flann::Matrix<double> flann_query(new double[m_env.dim], 1, m_env.dim);
    for (int i=0;i<m_env.dim;++i)
        flann_query[0][i] = cfg.t[i];

    // Search the kd tree for the nearest neighbor to the query.
    std::vector< std::vector<int> > query_match_indices;
    std::vector< std::vector<double> > query_distances;

    int num_neighbors_found = index_->knnSearch(
        flann_query, query_match_indices, query_distances, 1,
        flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED) /* no approx */);

    return m_tree[query_match_indices[0][0]];

//    // linear search
//    int idx(0);
//    double min_dist(m_tree[0]->cfg.distance(cfg));
//    for (int i=1;i<m_tree.size();++i) {
//        double dist(m_tree[i]->cfg.distance(cfg));
//        if (min_dist > dist) {
//            min_dist = dist;
//            idx = i;
//        }
//    }
//    return m_tree[idx];
}
