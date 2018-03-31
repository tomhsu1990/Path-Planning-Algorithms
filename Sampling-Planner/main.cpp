/*
 *  main.cpp
 *
 *  Created on: Mar. 22, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "main.hpp"

#define mw_out (*window)
MainWindow *window;

void run();

int main(int argc, char** argv) {

    f_prm = new FileParameter();
    srand(f_prm->seed);

    QApplication app(argc, argv);
    f_prm->checkPwd();
    f_prm->parseExampleList();
    f_prm->parseExampleFile();
    f_prm->parseMapFile(&planner);

    window = new MainWindow();
    run();
    window->show();

    return app.exec();
}

int run_count = 1;
void run () {
    mw_out<<"Run No. "<<run_count++<<"\n";
    mw_out<<"- method = "<<f_prm->method<<"\n";
    mw_out<<"- start = \t";
    for (int i=0;i<f_prm->dim_t;++i)
        mw_out<<f_prm->start.t[i]<<" ";
    mw_out<<"\t";
    for (int i=0;i<f_prm->dim_r;++i)
        mw_out<<f_prm->start.r[i]<<" ";
    mw_out<<"\n";
    mw_out<<"- goal = \t";
    for (int i=0;i<f_prm->dim_t;++i)
        mw_out<<f_prm->goal.t[i]<<" ";
    mw_out<<"\t";
    for (int i=0;i<f_prm->dim_r;++i)
        mw_out<<f_prm->goal.r[i]<<" ";
    mw_out<<"\n";
    mw_out<<"- planner max sample size = "<<(int)f_prm->max_sample_size<<"\n";

    if (!f_prm->method.compare("prm") || !f_prm->method.compare("PRM") || !f_prm->method.compare("Prm")) {
        mw_out<<"- prm connection k = "<<(int)f_prm->prm_closest_free_k<<"\n";
    }
    else if (!f_prm->method.compare("lazyprm") || !f_prm->method.compare("LAZYPRM") || !f_prm->method.compare("LazyPRM") || !f_prm->method.compare("LazyPrm")) {
        mw_out<<"- lazy prm connection k = "<<(int)f_prm->prm_closest_free_k<<"\n";
    }
    else if (!f_prm->method.compare("toggleprm") || !f_prm->method.compare("TOGGLEPRM") || !f_prm->method.compare("TogglePRM") || !f_prm->method.compare("TogglePrm")) {
        mw_out<<"- toggle prm connection k = "<<(int)f_prm->prm_closest_free_k<<"\n";
    }
    else if (!f_prm->method.compare("rrt") || !f_prm->method.compare("RRT") || !f_prm->method.compare("Rrt")) {
        mw_out<<"- rrt step size = "<<f_prm->rrt_step_size<<"\n";
        mw_out<<"- rrt goal bias = "<<f_prm->rrt_bias<<"\n";
        mw_out<<"- rrt close to goal = "<<f_prm->rrt_close_to_goal<<"\n";
    }
    mw_out<<"- seed = "<<f_prm->seed<<"\n";

    Timer t; t.start();
    planner->getRobot().init(f_prm->robot, f_prm->start);
    planner->setStart(planner->toParametric(f_prm->start));
    planner->setGoal(planner->toParametric(f_prm->goal));
    f_prm->no_path = !planner->findPath();
    f_prm->elapsed_time = t.getElapsedMilliseconds();
    f_prm->elapsed_CPU_time = t.getElapsedCPUMilliseconds();
    window->reportTime(run_count);

    if (!f_prm->no_path) {
        std::vector<Config> path = planner->getPath();
        mw_out<<"\n! Path found "<<"length = "<<(int)path.size()<<")\n";
    }
    else {
        mw_out<<"\n! Path not found\n";
    }
    mw_out<<"\nend\n\n";
}
