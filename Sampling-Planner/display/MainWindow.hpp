#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <string>

#include <QMainWindow>
#include <QMouseEvent>

#include "config/FileParameter.hpp"
#include "planner/Planner.hpp"
#include "Display.hpp"

#include "build/ui_MainWindow.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void reportTime(int run);

    // Print text
    virtual MainWindow& operator<<(const std::string&);
    virtual MainWindow& operator<<(const char*);
    virtual MainWindow& operator<<(int);
    virtual MainWindow& operator<<(long);
    virtual MainWindow& operator<<(float);
    virtual MainWindow& operator<<(double);

private slots:

    void on_run_clicked();
    void on_exit_clicked();

    void on_prm_clicked();
    void on_rrt_clicked();
    void on_toggle_clicked();
    void on_lazytoggle_clicked();

    void on_prm_graph_clicked();
    void on_prm_graph_mixed_clicked();
    void on_prm_graph_free_clicked();
    void on_prm_graph_obst_clicked();
    void on_prm_graph_edge_clicked();
    void on_rrt_graph_clicked();

    void on_show_clicked();
    void on_pause_clicked();
    void on_replay_clicked();

    void on_trace_clicked();
    void on_showFilledObstacles_clicked();

    void on_animationSpeed_valueChanged(int value);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
