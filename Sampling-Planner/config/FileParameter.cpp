#include "FileParameter.hpp"

FileParameter::FileParameter(){
    cfg_name = "rand-tri-100.cfg";
    input_dir = "inputs";
    file_name = "simple.txt";
    num_cfg = 0;

    env_scale = 1.0;
    env_delta = {0,0};

    seed = 101;
    timeout = 60;

    env.dim = 2;
    for(int i=0;i<env.dim;++i)
        env.space.push_back(512);
    env_TR = 0.1;
    env_RR = 1;

    method = "prm";
    sample_size = 100;
    max_sample_size = 1e4;

    rrt_step_size = 0.01;
    rrt_bias = 0.1;
    rrt_close_to_goal = 0.01;

    prm_closest_free_k = 10;
    prm_closest_obst_k = 2;

    window_pos = {50, 50};
    show_anim = true;
    pause_anim = replay_anim = false;
    show_prm_graph = true;
    show_prm_graph_mixed = show_prm_graph_free = show_prm_graph_obst = show_prm_graph_edge = false;
    show_rrt_graph = false;

    animation_speed = 99;
    animation_speed_scale = 5000;

    show_trace = false;
    show_filled_obstacles = true;
}

FileParameter::~FileParameter(){
}

void FileParameter::checkPwd () {
    bool found = false;
    QString project_name = QCoreApplication::applicationName();
    working_dir = QDir::currentPath().toStdString();

    // Test if the build directory is "sampling-planner". If so,
    // the path to the current working directory will
    // include "sampling-planner"
    unsigned long desired_dir = working_dir.rfind("/"+project_name.toStdString()+"/");
    if (desired_dir != std::string::npos) {
        working_dir = working_dir.substr(0, desired_dir+project_name.size()+1);

        // Set current working directory to "sampling-planner"
        QDir::setCurrent(working_dir.c_str());
        found = true;
    }

    // Test if a build directory (/build-sampling-planner-...) was created. This directory
    // will reside in the same directory as "sampling-planner"
    if (!found && (desired_dir = working_dir.rfind("/build/sampling-planner")) != std::string::npos) {
        QDir dir(working_dir.substr(0, desired_dir).c_str());
        if (dir.exists("sampling-planner.pro")) { // Test if /sampling-planner exists
            working_dir = working_dir.substr(0, desired_dir);

            // Set current working directory to "sampling-planner"
            QDir::setCurrent(working_dir.c_str());
            found = true;
        }
    }

    // /sampling-planner could not be found
    if (!found) {
        std::cerr << std::endl << "!! WARNING !!\n"
        << "The program may not work correctly or at all because the folder "
           "containing the program's files cannot be found.\n"
           "Make sure that the program is inside of a folder named \"sampling-planner\".\n";
    }
}

void FileParameter::parseExampleList () {
    std::string s;
    QDir cfg_dir(std::string(working_dir + "/" + input_dir).c_str());
    QStringList file_list = cfg_dir.entryList();
    while (!file_list.empty()){
        s = file_list.front().toStdString();
        int size = s.size();
        if (size > 4 && s[size-1] == 'g' && s[size-2] == 'f' && s[size-3] == 'c' && s[size-4] == '.'){
            cfg_name_list.push_back(s);
            ++num_cfg;
        }
        file_list.pop_front();
    }
}

void FileParameter::parseExampleFile () {
    char tmp[256];
    FILE *fptr = fopen(std::string(input_dir+"/"+cfg_name).c_str(), "r");
    if (fptr == NULL) return ;

    while (fgets(tmp, 256, fptr) != NULL){
        char *sptr = strtok(tmp, "=: \t\n");

        if (sptr == NULL) {
            continue;
        }
        // comments
        if (strcmp(sptr, "#") == 0) {
            continue;
        }

        if (strcmp(sptr, "dimension") == 0) {
            sptr = strtok(NULL, "=: #\t");
            env.dim = atoi(sptr);
            for (int i=0;i<env.dim;++i) {
                sptr = strtok(NULL, "=: #\t");
                env.space.push_back(atof(sptr));
            }
            sptr = strtok(NULL, "=: #\t");
            dim_t = atoi(sptr);
            sptr = strtok(NULL, "=: #\t");
            dim_r = atoi(sptr);
            start.init(dim_t, dim_r, false);
            goal.init(dim_t, dim_r, false);
        }

        // start configuration
        if (strcmp(sptr, "start") == 0) {
            for (int i=0;i<dim_t;++i) {
                sptr = strtok(NULL, "=: #\t(),");
                start.t[i] = atof(sptr);
            }
            for (int i=0;i<dim_r;++i) {
                sptr = strtok(NULL, "=: #\t(),");
                start.r[i] = atof(sptr);
            }
        }

        // goal configuration
        if (strcmp(sptr, "goal") == 0) {
            for (int i=0;i<dim_t;++i) {
                sptr = strtok(NULL, "=: #\t(),");
                goal.t[i] = atof(sptr);
            }
            for (int i=0;i<dim_r;++i) {
                sptr = strtok(NULL, "=: #\t(),");
                goal.r[i] = atof(sptr);
            }
        }

        // robot
        if (strcmp(sptr, "robot") == 0) {
            sptr = strtok(NULL, "=: #\t");
            robot.name = sptr;

            if (robot.name.compare("disc") == 0) {
                sptr = strtok(NULL, "=: #\t");
                robot.R = atof(sptr);
            }
            if (robot.name.compare("1link") == 0) {
                sptr = strtok(NULL, "=: #\t");
                robot.L1 = atof(sptr);
                sptr = strtok(NULL, "=: #\t");
                robot.Thickness = atof(sptr);
            }
            if (robot.name.compare("2links") == 0) {
                sptr = strtok(NULL, "=: #\t");
                robot.L1 = atof(sptr);
                sptr = strtok(NULL, "=: #\t");
                robot.L2 = atof(sptr);
                sptr = strtok(NULL, "=: #\t");
                robot.Thickness = atof(sptr);
            }
        }

        if (strcmp(sptr, "input_dir") == 0) {
            sptr = strtok(NULL, "=: #\t");
            input_dir = sptr;
        }
        if (strcmp(sptr, "file_name") == 0) {
            sptr = strtok(NULL, "=: #\t");
            file_name = sptr;
        }

        // windows position
        if (strcmp(sptr, "windows_pos") == 0) {
            sptr = strtok(NULL, "=: \t");
            window_pos.first = atoi(sptr);
            sptr = strtok(NULL, "=: \t");
            window_pos.second = atoi(sptr);
        }

        if (strcmp(sptr, "method") == 0) {
            sptr = strtok(NULL, "=: \t\n");
            method = sptr;
        }

        if (strcmp(sptr, "seed") == 0) {
            sptr = strtok(NULL, "=: \t");
            seed = atoi(sptr);
        }
        if (strcmp(sptr, "timeout") == 0) {
            sptr = strtok(NULL, "=: \t");
            timeout = atof(sptr);
        }

        // environment delta
        if (strcmp(sptr, "env_delta") == 0) {
            sptr = strtok(NULL, "=: \t");
            env_delta.first = atof(sptr);
            sptr = strtok(NULL, "=: \t");
            env_delta.second = atof(sptr);
            sptr = strtok(NULL, "=: \t");
            env_scale = atof(sptr);
        }

        if (strcmp(sptr, "sample_size") == 0) {
            sptr = strtok(NULL, "=: \t");
            sample_size = atoi(sptr);
        }
        if (strcmp(sptr, "max_sample_size") == 0) {
            sptr = strtok(NULL, "=: \t");
            max_sample_size = atoi(sptr);
        }

        if (strcmp(sptr, "prm_closest_free_k") == 0) {
            sptr = strtok(NULL, "=: \t");
            prm_closest_free_k = atoi(sptr);
        }
        if (strcmp(sptr, "prm_closest_free_k") == 0) {
            sptr = strtok(NULL, "=: \t");
            prm_closest_free_k = atoi(sptr);
        }
        if (strcmp(sptr, "prm_closest_obst_k") == 0) {
            sptr = strtok(NULL, "=: \t");
            prm_closest_obst_k = atoi(sptr);
        }

        if (strcmp(sptr, "rrt_step_size") == 0) {
            sptr = strtok(NULL, "=: \t");
            rrt_step_size = atof(sptr);
        }
        if (strcmp(sptr, "rrt_bias") == 0) {
            sptr = strtok(NULL, "=: \t");
            rrt_bias = atof(sptr);
        }
        if (strcmp(sptr, "rrt_close_to_goal") == 0) {
            sptr = strtok(NULL, "=: \t");
            rrt_close_to_goal = atof(sptr);
        }

        if (strcmp(sptr, "show_prm_graph") == 0) {
            sptr = strtok(NULL, "=: \t");
            show_prm_graph = atoi(sptr);
        }
        if (strcmp(sptr, "show_rrt_graph") == 0) {
            sptr = strtok(NULL, "=: \t");
            show_rrt_graph = atoi(sptr);
        }

        if (strcmp(sptr, "animation_speed") == 0) {
            sptr = strtok(NULL, "=: \t");
            animation_speed = atoi(sptr);
        }
        if (strcmp(sptr, "animation_speed_scale") == 0) {
            sptr = strtok(NULL, "=: \t");
            animation_speed_scale = atoi(sptr);
        }
    }
}

void FileParameter::parseMapFile (Planner **planner) {

    free((*planner));
    if (!method.compare("prm") || !method.compare("PRM") || !method.compare("Prm")) {
        (*planner) = new PRM(env, env_TR, env_RR, start, goal,
                             max_sample_size,prm_closest_free_k);
    }
    else if (!method.compare("rrt") || !method.compare("RRT") || !method.compare("Rrt")) {
        (*planner) = new RRT(env, env_TR, env_RR, start, goal,
                             max_sample_size, rrt_step_size,rrt_bias,rrt_close_to_goal);
    }
    if ((*planner) == NULL) return;

    std::string full_name = input_dir+"/"+file_name;
    std::ifstream fin;
    fin.open(full_name.c_str());

    char buf[1024];
    int num_pts = -1;
    int num_polys = -1;
    std::vector<Point2d> pts;
    while (fin.getline(buf, 1024)) {
        std::string lline(buf);
        if(buf[0] == '#' || lline.size() == 0) continue;

        std::stringstream ss(lline);
        if (num_pts < 0) ss>>num_pts;
        else if (num_pts != pts.size()) {
            double x, y;
            ss>>x>>y;
            pts.push_back(Point2d(x, y));
        }
        else if (num_polys < 0) {
            ss>>num_polys;
        }
        else {
            Polygon2d poly;
            int idx;
            while (ss>>idx) {
                --idx;
                poly.addPoint(pts[idx].X()*env_scale+env_delta.first,
                              pts[idx].Y()*env_scale+env_delta.second);
            }
            poly.setOrientation();
            (*planner)->addPolygon(poly);
        }
    }
}
