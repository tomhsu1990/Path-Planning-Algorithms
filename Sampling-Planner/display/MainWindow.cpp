/*
 *  MainWindow.cpp
 *
 *  Created on: Mar. 24, 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#include "MainWindow.hpp"

extern FileParameter *f_prm;
extern Planner *planner;

extern void run();

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    for (int i=0;i<f_prm->num_cfg;++i) {
        ui->example_list->addItem(f_prm->cfg_name_list[i].c_str());
    }
    for (int i=0;i<f_prm->num_cfg;++i) {
        if (f_prm->cfg_name_list[i].compare(f_prm->cfg_name.c_str()) == 0) {
            ui->example_list->setCurrentIndex(i);
        }
    }
    ui->input_file->setText(QString::fromStdString(f_prm->file_name.substr(0,f_prm->file_name.size()-4)));

    ui->start_x->setValue(f_prm->start.t[0]);
    ui->start_y->setValue(f_prm->start.t[1]);

    ui->goal_x->setValue(f_prm->goal.t[0]);
    ui->goal_y->setValue(f_prm->goal.t[1]);

    ui->robot_name->setText(QString::fromStdString(f_prm->robot.name));
    if (f_prm->robot.name.compare("disc") == 0) {
        ui->R->setValue(f_prm->robot.R);
    }
    else if (f_prm->robot.name.compare("1link") == 0) {

    }
    else if (f_prm->robot.name.compare("2links") == 0) {

    }

    if (f_prm->method.compare("prm") || f_prm->method.compare("PRM") || f_prm->method.compare("Prm")) {
        ui->prm->setChecked(true);
    }
    else if (f_prm->method.compare("rrt") || f_prm->method.compare("RRT") || f_prm->method.compare("Rrt")) {
        ui->rrt->setChecked(true);
    }

    ui->max_sample->setValue(f_prm->max_sample_size);

    ui->prm_closest_free_k->setValue(f_prm->prm_closest_free_k);
    ui->prm_closest_obst_k->setValue(f_prm->prm_closest_obst_k);

    ui->rrt_step_size->setValue(f_prm->rrt_step_size);
    ui->rrt_bias->setValue(f_prm->rrt_bias);
    ui->rrt_close_to_goal->setValue(f_prm->rrt_close_to_goal);

    ui->prm_graph->setChecked(f_prm->show_prm_graph);
    ui->prm_graph_mixed->setChecked(f_prm->show_prm_graph_mixed);
    ui->prm_graph_free->setChecked(f_prm->show_prm_graph_free);
    ui->prm_graph_obst->setChecked(f_prm->show_prm_graph_obst);
    ui->prm_graph_edge->setChecked(f_prm->show_prm_graph_edge);

    ui->rrt_graph->setChecked(f_prm->show_rrt_graph);

    ui->animation_speed->setValue(f_prm->animation_speed);
    ui->seed->setValue(f_prm->seed);
    srand(f_prm->seed);

    ui->timeout->setValue(f_prm->timeout);

    ui->text_output_time->setText(QString::fromStdString(""));
}

MainWindow::~MainWindow() {
    delete ui;
}

//================================//
//     Display Text to Window     //
//================================//
/* Scrolls the text screen to the bottom and prints text */
MainWindow& MainWindow::operator<< (const std::string& str) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(str.c_str());
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

MainWindow& MainWindow::operator<< (const char* text) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(text);
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

MainWindow& MainWindow::operator<< (int i) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(QString::number(i));
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

MainWindow& MainWindow::operator<< (long l) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(QString::number(l));
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

MainWindow& MainWindow::operator<< (float f) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(QString::number(f));
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

MainWindow& MainWindow::operator<< (double d) {
    ui->textOutput->moveCursor (QTextCursor::End);
    ui->textOutput->insertPlainText(QString::number(d));
    ui->textOutput->moveCursor (QTextCursor::End);
    return *this;
}

void MainWindow::on_run_clicked() {

    if (ui->example_list->currentText().toStdString().compare(f_prm->cfg_name) == 0) {
        f_prm->file_name=ui->input_file->text().toStdString()+".txt";

        f_prm->start.t[0] = ui->start_x->value();
        f_prm->start.t[1] = ui->start_y->value();

        f_prm->goal.t[0] = ui->goal_x->value();
        f_prm->goal.t[1] = ui->goal_y->value();

        f_prm->robot.name = ui->robot_name->text().toStdString();
        if (f_prm->robot.name.compare("disc") == 0) {
            f_prm->robot.R = ui->R->value();
        }
        else if (f_prm->robot.name.compare("1link") == 0) {

        }
        else if (f_prm->robot.name.compare("2links") == 0) {

        }


        f_prm->max_sample_size=ui->max_sample->value();

        f_prm->prm_closest_free_k=ui->prm_closest_free_k->value();
        f_prm->prm_closest_obst_k=ui->prm_closest_obst_k->value();

        f_prm->rrt_step_size=ui->rrt_step_size->value();
        f_prm->rrt_bias=ui->rrt_bias->value();
        f_prm->rrt_close_to_goal=ui->rrt_close_to_goal->value();

        f_prm->show_prm_graph=ui->prm_graph->isChecked();
        f_prm->show_rrt_graph=ui->rrt_graph->isChecked();

        // 4/13/2016
        // only initialize the seed when it is changed
        int new_seed = ui->seed->value();
        if (f_prm->seed != new_seed) {
            f_prm->seed = new_seed;
            srand(f_prm->seed);
        }

        f_prm->timeout = ui->timeout->value();
    }
    else {
        f_prm->cfg_name = ui->example_list->currentText().toStdString();
        f_prm->parseExampleFile();

        ui->input_file->setText(QString::fromStdString(f_prm->file_name.substr(0,f_prm->file_name.length()-4)));

        ui->start_x->setValue(f_prm->start.t[0]);
        ui->start_y->setValue(f_prm->start.t[1]);

        ui->goal_x->setValue(f_prm->goal.t[0]);
        ui->goal_y->setValue(f_prm->goal.t[1]);

        ui->robot_name->setText(QString::fromStdString(f_prm->robot.name));
        if (f_prm->robot.name.compare("disc") == 0) {
            ui->R->setValue(f_prm->robot.R);
        }
        else if (f_prm->robot.name.compare("1link") == 0) {

        }
        else if (f_prm->robot.name.compare("2links") == 0) {

        }

        if (f_prm->method.compare("prm") || f_prm->method.compare("PRM") || f_prm->method.compare("Prm")) {
            ui->prm->setChecked(true);
        }
        else if (f_prm->method.compare("rrt") || f_prm->method.compare("RRT") || f_prm->method.compare("Rrt")) {
            ui->rrt->setChecked(true);
        }
        ui->max_sample->setValue(f_prm->max_sample_size);

        ui->prm_closest_free_k->setValue(f_prm->prm_closest_free_k);
        ui->prm_closest_obst_k->setValue(f_prm->prm_closest_obst_k);

        ui->rrt_step_size->setValue(f_prm->rrt_step_size);
        ui->rrt_bias->setValue(f_prm->rrt_bias);
        ui->rrt_close_to_goal->setValue(f_prm->rrt_close_to_goal);

        ui->prm_graph->setChecked(f_prm->show_prm_graph);
        ui->rrt_graph->setChecked(f_prm->show_rrt_graph);

        int old_seed = ui->seed->value();
        if (f_prm->seed != old_seed) {
            f_prm->seed = old_seed;
            ui->seed->setValue(f_prm->seed);
            srand(f_prm->seed);
        }

        ui->timeout->setValue(f_prm->timeout);
    }

    f_prm->show_anim = true;
    f_prm->pause_anim = false;

    f_prm->parseMapFile(&planner);
    run();

    ui->openGLWidget->update();
}

void MainWindow::reportTime(int run) {
    std::string time_info;
    char tmp_buff[200];
    sprintf(tmp_buff, "Run #%d\n", run);
    time_info.append(tmp_buff);
    sprintf(tmp_buff, "Elapsed time:\n%lf (ms)\n", f_prm->elapsed_time);
    time_info.append(tmp_buff);
    sprintf(tmp_buff, "Elapsed CPU time:\n%lf (ms)\n", f_prm->elapsed_CPU_time);
    time_info.append(tmp_buff);
    ui->text_output_time->append(QString::fromStdString(time_info));
}

void MainWindow::on_prm_clicked() {
    f_prm->method = "prm";
}

void MainWindow::on_rrt_clicked() {
    f_prm->method = "rrt";
}

void MainWindow::on_toggle_prm_clicked() {
    f_prm->method = "toggle";
}

void MainWindow::on_lazy_toggle_prm_clicked() {
    f_prm->method = "lazytoggle";
}

void MainWindow::on_prm_graph_clicked(){
    f_prm->show_prm_graph = ui->prm_graph->isChecked();
    this->update();
}
void MainWindow::on_prm_graph_mixed_clicked(){
    f_prm->show_prm_graph_mixed = ui->prm_graph_mixed->isChecked();
    this->update();
}
void MainWindow::on_prm_graph_free_clicked(){
    f_prm->show_prm_graph_free = ui->prm_graph_free->isChecked();
    this->update();
}
void MainWindow::on_prm_graph_obst_clicked(){
    f_prm->show_prm_graph_obst = ui->prm_graph_obst->isChecked();
    this->update();
}
void MainWindow::on_prm_graph_edge_clicked(){
    f_prm->show_prm_graph_edge = ui->prm_graph_edge->isChecked();
    this->update();
}

void MainWindow::on_rrt_graph_clicked(){
    f_prm->show_rrt_graph = ui->rrt_graph->isChecked();
    this->update();
}

void MainWindow::on_exit_clicked() {
    this->close();
}

void MainWindow::on_show_clicked() {
    f_prm->show_anim = true;
    this->update();
}

void MainWindow::on_pause_clicked() {
    f_prm->pause_anim = !f_prm->pause_anim;
    this->update();
}

void MainWindow::on_replay_clicked() {
    f_prm->replay_anim = true;
    this->update();
}

void MainWindow::on_show_trace_clicked() {
    f_prm->show_trace = !f_prm->show_trace;
    this->update();
}
void MainWindow::on_show_filled_obstacles_clicked() {
    f_prm->show_filled_obstacles = !f_prm->show_filled_obstacles;
    this->update();
}

void MainWindow::on_animation_speed_valueChanged(int value) {
    f_prm->animation_speed = value;
}
