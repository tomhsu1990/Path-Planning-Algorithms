/*
 *  FileParameter.hpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#ifndef FILEPARAMETER_H
#define FILEPARAMETER_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <QDir>
#include <QApplication>
#include <QSurfaceFormat>

#include "CORE/geom2d/point2d.h"
#include "CORE/geom2d/line2d.h"
#include "CORE/geom2d/segment2d.h"
#include "CORE/geom2d/polygon2d.h"
#include "Config.hpp"
#include "Robot.hpp"
#include "planner/RRT.hpp"
#include "planner/PRM.hpp"
#include "planner/LazyPRM.hpp"
#include "planner/TogglePRM.hpp"
#include "planner/LazyTogglePRM.hpp"

class FileParameter {
public:
    FileParameter ();
    ~FileParameter ();

    void checkPwd ();
    void parseExampleList ();
    void parseExampleFile ();
    void parseMapFile (Planner *planner);

    //working environment parameters
    std::string working_dir, input_dir;
    std::string cfg_name, file_name;
    std::vector<std::string> cfg_name_list;
    int num_cfg;

    // working space;
    Env env;

    // configuration parameters
    int dim_t, dim_r;
    Config start, goal;

    // robot paramters
    Robot robot;

    // environment(map) parameters
    float env_scale;
    std::pair<double, double> env_delta;

    // running parameters
    int seed;
    double timeout;
    double env_TR;			// TRANSLATIONAL RESOLUTION
    double env_RR;			// Rotational RESOLUTION (deg)

    // planner parameters
    std::string method;
    unsigned int sample_size;
    unsigned int max_sample_size;
    //rrt parameters
    double rrt_step_size;
    double rrt_bias;
    double rrt_close_to_goal;
    //lazy toggle prm parameters
    unsigned int prm_closest_free_k;
    unsigned int prm_closest_obst_k;

    // display parameters
    std::pair<int, int> window_pos; // Position of Window
    // animation
    bool no_path;
    bool show_anim;
    bool pause_anim;
    bool replay_anim;
    bool show_rrt_graph;
    bool show_prm_graph;
    bool show_prm_graph_mixed;
    bool show_prm_graph_free;
    bool show_prm_graph_obst;
    bool show_prm_graph_edge;
    int animation_speed;
    int animation_speed_scale;
    bool show_trace;
    bool show_filled_obstacles;

    // timing parameters
    double elapsed_time, elapsed_CPU_time;
};

#endif // FILEPARAMETER_H
