

#include "ANNHelper.hpp"

int ANNHelper::findClosestPt(const Config &cfg){
    for(int i=0;i<dim;++i)
        query_pt[i] = cfg.t[i];
	// find closest point to compute the error distance
	kd_tree->annkSearch(				// search
		query_pt,						// query point
		k,								// number of near neighbors
		nn_idx,							// nearest neighbors (returned)
		dists,							// distance (returned)
		eps);							// error bound
    return nn_idx[0];
}

void ANNHelper::initializeKdTree () {
	kd_tree = new ANNkd_tree(			// build search structure
				data_pts,				// the data points
				n_pts,					// number of points
				dim);					// dimension of space
}

void ANNHelper::addPt (Config cfg) {
    //kd_tree->
}
