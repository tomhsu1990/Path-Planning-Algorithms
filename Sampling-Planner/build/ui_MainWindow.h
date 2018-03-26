/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.5.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <display/Display.hpp>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QTextBrowser *textOutput;
    QPushButton *run;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *inputFile;
    QLabel *label_2;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QLabel *label_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_4;
    QDoubleSpinBox *aX;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_5;
    QDoubleSpinBox *aY;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *verticalLayout_2;
    QLabel *label_8;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_9;
    QDoubleSpinBox *bX;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_10;
    QDoubleSpinBox *bY;
    QWidget *verticalLayoutWidget_3;
    QVBoxLayout *verticalLayout_3;
    QLabel *label_13;
    QHBoxLayout *horizontalLayout_11;
    QLabel *label_15;
    QLineEdit *inputFile_2;
    QHBoxLayout *horizontalLayout_10;
    QLabel *label_14;
    QDoubleSpinBox *l1;
    QWidget *horizontalLayoutWidget_9;
    QHBoxLayout *horizontalLayout_13;
    QLabel *label_17;
    QSpinBox *random;
    QWidget *verticalLayoutWidget_4;
    QVBoxLayout *verticalLayout_4;
    QRadioButton *prm;
    QRadioButton *toggle;
    QRadioButton *lazytoggle;
    QRadioButton *rrt;
    QPushButton *exit;
    QFrame *line;
    QFrame *line_2;
    QFrame *line_3;
    QFrame *line_4;
    Display *openGLWidget;
    QPushButton *show;
    QWidget *verticalLayoutWidget_6;
    QVBoxLayout *verticalLayout_6;
    QLabel *label_18;
    QSlider *animationSpeed;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout_16;
    QLabel *label_20;
    QComboBox *comboBox;
    QPushButton *pause;
    QPushButton *replay;
    QWidget *verticalLayoutWidget_5;
    QVBoxLayout *verticalLayout_5;
    QLabel *label_16;
    QHBoxLayout *horizontalLayout_12;
    QLabel *label_19;
    QDoubleSpinBox *prm_closest_free_k;
    QHBoxLayout *horizontalLayout_24;
    QLabel *label_34;
    QDoubleSpinBox *prm_closest_obst_k;
    QWidget *verticalLayoutWidget_8;
    QVBoxLayout *verticalLayout_8;
    QLabel *label_25;
    QHBoxLayout *horizontalLayout_18;
    QLabel *label_26;
    QDoubleSpinBox *rrt_step_size;
    QHBoxLayout *horizontalLayout_21;
    QLabel *label_29;
    QDoubleSpinBox *rrt_bias;
    QHBoxLayout *horizontalLayout_19;
    QLabel *label_27;
    QDoubleSpinBox *rrt_close_to_goal;
    QWidget *layoutWidget_2;
    QHBoxLayout *horizontalLayout_22;
    QLabel *label_30;
    QDoubleSpinBox *max_sample;
    QCheckBox *prm_graph;
    QCheckBox *rrt_graph;
    QTextBrowser *textOutputTime;
    QCheckBox *prm_graph_free;
    QCheckBox *prm_graph_obst;
    QWidget *layoutWidget_3;
    QHBoxLayout *horizontalLayout_23;
    QLabel *label_31;
    QDoubleSpinBox *timeout;
    QLabel *label_33;
    QCheckBox *prm_graph_edge;
    QPushButton *trace;
    QPushButton *showFilledObstacles;
    QCheckBox *prm_graph_mixed;
    QMenuBar *menuBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(1117, 725);
        MainWindow->setMaximumSize(QSize(1200, 750));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        textOutput = new QTextBrowser(centralWidget);
        textOutput->setObjectName(QStringLiteral("textOutput"));
        textOutput->setGeometry(QRect(730, 450, 381, 261));
        run = new QPushButton(centralWidget);
        run->setObjectName(QStringLiteral("run"));
        run->setGeometry(QRect(1024, 410, 91, 32));
        horizontalLayoutWidget = new QWidget(centralWidget);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(740, 30, 231, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        inputFile = new QLineEdit(horizontalLayoutWidget);
        inputFile->setObjectName(QStringLiteral("inputFile"));

        horizontalLayout->addWidget(inputFile);

        label_2 = new QLabel(horizontalLayoutWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout->addWidget(label_2);

        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QStringLiteral("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(940, 70, 171, 160));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(verticalLayoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        verticalLayout->addWidget(label_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_4 = new QLabel(verticalLayoutWidget);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_2->addWidget(label_4);

        aX = new QDoubleSpinBox(verticalLayoutWidget);
        aX->setObjectName(QStringLiteral("aX"));
        aX->setDecimals(1);
        aX->setMaximum(1e+08);

        horizontalLayout_2->addWidget(aX);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_5 = new QLabel(verticalLayoutWidget);
        label_5->setObjectName(QStringLiteral("label_5"));

        horizontalLayout_3->addWidget(label_5);

        aY = new QDoubleSpinBox(verticalLayoutWidget);
        aY->setObjectName(QStringLiteral("aY"));
        aY->setDecimals(1);
        aY->setMaximum(1e+08);

        horizontalLayout_3->addWidget(aY);


        verticalLayout->addLayout(horizontalLayout_3);

        verticalLayoutWidget_2 = new QWidget(centralWidget);
        verticalLayoutWidget_2->setObjectName(QStringLiteral("verticalLayoutWidget_2"));
        verticalLayoutWidget_2->setGeometry(QRect(940, 240, 171, 160));
        verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_8 = new QLabel(verticalLayoutWidget_2);
        label_8->setObjectName(QStringLiteral("label_8"));

        verticalLayout_2->addWidget(label_8);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_9 = new QLabel(verticalLayoutWidget_2);
        label_9->setObjectName(QStringLiteral("label_9"));

        horizontalLayout_6->addWidget(label_9);

        bX = new QDoubleSpinBox(verticalLayoutWidget_2);
        bX->setObjectName(QStringLiteral("bX"));
        bX->setDecimals(1);
        bX->setMaximum(1e+08);

        horizontalLayout_6->addWidget(bX);


        verticalLayout_2->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        label_10 = new QLabel(verticalLayoutWidget_2);
        label_10->setObjectName(QStringLiteral("label_10"));

        horizontalLayout_7->addWidget(label_10);

        bY = new QDoubleSpinBox(verticalLayoutWidget_2);
        bY->setObjectName(QStringLiteral("bY"));
        bY->setDecimals(1);
        bY->setMaximum(1e+08);

        horizontalLayout_7->addWidget(bY);


        verticalLayout_2->addLayout(horizontalLayout_7);

        verticalLayoutWidget_3 = new QWidget(centralWidget);
        verticalLayoutWidget_3->setObjectName(QStringLiteral("verticalLayoutWidget_3"));
        verticalLayoutWidget_3->setGeometry(QRect(740, 70, 181, 151));
        verticalLayout_3 = new QVBoxLayout(verticalLayoutWidget_3);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_13 = new QLabel(verticalLayoutWidget_3);
        label_13->setObjectName(QStringLiteral("label_13"));

        verticalLayout_3->addWidget(label_13);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setSpacing(6);
        horizontalLayout_11->setObjectName(QStringLiteral("horizontalLayout_11"));
        label_15 = new QLabel(verticalLayoutWidget_3);
        label_15->setObjectName(QStringLiteral("label_15"));

        horizontalLayout_11->addWidget(label_15);

        inputFile_2 = new QLineEdit(verticalLayoutWidget_3);
        inputFile_2->setObjectName(QStringLiteral("inputFile_2"));

        horizontalLayout_11->addWidget(inputFile_2);


        verticalLayout_3->addLayout(horizontalLayout_11);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        label_14 = new QLabel(verticalLayoutWidget_3);
        label_14->setObjectName(QStringLiteral("label_14"));

        horizontalLayout_10->addWidget(label_14);

        l1 = new QDoubleSpinBox(verticalLayoutWidget_3);
        l1->setObjectName(QStringLiteral("l1"));
        l1->setDecimals(0);
        l1->setMaximum(1000);

        horizontalLayout_10->addWidget(l1);


        verticalLayout_3->addLayout(horizontalLayout_10);

        horizontalLayoutWidget_9 = new QWidget(centralWidget);
        horizontalLayoutWidget_9->setObjectName(QStringLiteral("horizontalLayoutWidget_9"));
        horizontalLayoutWidget_9->setGeometry(QRect(1000, 30, 111, 31));
        horizontalLayout_13 = new QHBoxLayout(horizontalLayoutWidget_9);
        horizontalLayout_13->setSpacing(6);
        horizontalLayout_13->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_13->setObjectName(QStringLiteral("horizontalLayout_13"));
        horizontalLayout_13->setContentsMargins(0, 0, 0, 0);
        label_17 = new QLabel(horizontalLayoutWidget_9);
        label_17->setObjectName(QStringLiteral("label_17"));

        horizontalLayout_13->addWidget(label_17);

        random = new QSpinBox(horizontalLayoutWidget_9);
        random->setObjectName(QStringLiteral("random"));

        horizontalLayout_13->addWidget(random);

        verticalLayoutWidget_4 = new QWidget(centralWidget);
        verticalLayoutWidget_4->setObjectName(QStringLiteral("verticalLayoutWidget_4"));
        verticalLayoutWidget_4->setGeometry(QRect(740, 230, 181, 111));
        verticalLayout_4 = new QVBoxLayout(verticalLayoutWidget_4);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        verticalLayout_4->setContentsMargins(0, 0, 0, 0);
        prm = new QRadioButton(verticalLayoutWidget_4);
        prm->setObjectName(QStringLiteral("prm"));

        verticalLayout_4->addWidget(prm);

        toggle = new QRadioButton(verticalLayoutWidget_4);
        toggle->setObjectName(QStringLiteral("toggle"));

        verticalLayout_4->addWidget(toggle);

        lazytoggle = new QRadioButton(verticalLayoutWidget_4);
        lazytoggle->setObjectName(QStringLiteral("lazytoggle"));

        verticalLayout_4->addWidget(lazytoggle);

        rrt = new QRadioButton(verticalLayoutWidget_4);
        rrt->setObjectName(QStringLiteral("rrt"));

        verticalLayout_4->addWidget(rrt);

        exit = new QPushButton(centralWidget);
        exit->setObjectName(QStringLiteral("exit"));
        exit->setGeometry(QRect(940, 410, 91, 32));
        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setGeometry(QRect(920, 70, 16, 371));
        line->setFrameShadow(QFrame::Sunken);
        line->setLineWidth(2);
        line->setFrameShape(QFrame::VLine);
        line_2 = new QFrame(centralWidget);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setGeometry(QRect(740, 60, 371, 16));
        line_2->setLineWidth(2);
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        line_3 = new QFrame(centralWidget);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setGeometry(QRect(740, 220, 181, 16));
        line_3->setLineWidth(2);
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);
        line_4 = new QFrame(centralWidget);
        line_4->setObjectName(QStringLiteral("line_4"));
        line_4->setGeometry(QRect(940, 230, 171, 16));
        line_4->setLineWidth(2);
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);
        openGLWidget = new Display(centralWidget);
        openGLWidget->setObjectName(QStringLiteral("openGLWidget"));
        openGLWidget->setGeometry(QRect(30, 10, 512, 512));
        openGLWidget->setMaximumSize(QSize(512, 512));
        show = new QPushButton(centralWidget);
        show->setObjectName(QStringLiteral("show"));
        show->setGeometry(QRect(730, 350, 61, 32));
        verticalLayoutWidget_6 = new QWidget(centralWidget);
        verticalLayoutWidget_6->setObjectName(QStringLiteral("verticalLayoutWidget_6"));
        verticalLayoutWidget_6->setGeometry(QRect(750, 390, 160, 51));
        verticalLayout_6 = new QVBoxLayout(verticalLayoutWidget_6);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        verticalLayout_6->setContentsMargins(0, 0, 0, 0);
        label_18 = new QLabel(verticalLayoutWidget_6);
        label_18->setObjectName(QStringLiteral("label_18"));

        verticalLayout_6->addWidget(label_18);

        animationSpeed = new QSlider(verticalLayoutWidget_6);
        animationSpeed->setObjectName(QStringLiteral("animationSpeed"));
        animationSpeed->setValue(90);
        animationSpeed->setSliderPosition(90);
        animationSpeed->setOrientation(Qt::Horizontal);

        verticalLayout_6->addWidget(animationSpeed);

        horizontalLayoutWidget_2 = new QWidget(centralWidget);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(740, 0, 371, 31));
        horizontalLayout_16 = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout_16->setSpacing(6);
        horizontalLayout_16->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_16->setObjectName(QStringLiteral("horizontalLayout_16"));
        horizontalLayout_16->setContentsMargins(0, 0, 0, 0);
        label_20 = new QLabel(horizontalLayoutWidget_2);
        label_20->setObjectName(QStringLiteral("label_20"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label_20->sizePolicy().hasHeightForWidth());
        label_20->setSizePolicy(sizePolicy);

        horizontalLayout_16->addWidget(label_20);

        comboBox = new QComboBox(horizontalLayoutWidget_2);
        comboBox->setObjectName(QStringLiteral("comboBox"));

        horizontalLayout_16->addWidget(comboBox);

        pause = new QPushButton(centralWidget);
        pause->setObjectName(QStringLiteral("pause"));
        pause->setGeometry(QRect(790, 350, 61, 32));
        replay = new QPushButton(centralWidget);
        replay->setObjectName(QStringLiteral("replay"));
        replay->setGeometry(QRect(850, 350, 71, 32));
        verticalLayoutWidget_5 = new QWidget(centralWidget);
        verticalLayoutWidget_5->setObjectName(QStringLiteral("verticalLayoutWidget_5"));
        verticalLayoutWidget_5->setGeometry(QRect(260, 540, 201, 101));
        verticalLayout_5 = new QVBoxLayout(verticalLayoutWidget_5);
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setContentsMargins(11, 11, 11, 11);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        verticalLayout_5->setContentsMargins(0, 0, 0, 0);
        label_16 = new QLabel(verticalLayoutWidget_5);
        label_16->setObjectName(QStringLiteral("label_16"));

        verticalLayout_5->addWidget(label_16);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setObjectName(QStringLiteral("horizontalLayout_12"));
        label_19 = new QLabel(verticalLayoutWidget_5);
        label_19->setObjectName(QStringLiteral("label_19"));

        horizontalLayout_12->addWidget(label_19);

        prm_closest_free_k = new QDoubleSpinBox(verticalLayoutWidget_5);
        prm_closest_free_k->setObjectName(QStringLiteral("prm_closest_free_k"));
        QSizePolicy sizePolicy1(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(prm_closest_free_k->sizePolicy().hasHeightForWidth());
        prm_closest_free_k->setSizePolicy(sizePolicy1);
        prm_closest_free_k->setDecimals(0);
        prm_closest_free_k->setMaximum(1e+06);

        horizontalLayout_12->addWidget(prm_closest_free_k);


        verticalLayout_5->addLayout(horizontalLayout_12);

        horizontalLayout_24 = new QHBoxLayout();
        horizontalLayout_24->setSpacing(6);
        horizontalLayout_24->setObjectName(QStringLiteral("horizontalLayout_24"));
        label_34 = new QLabel(verticalLayoutWidget_5);
        label_34->setObjectName(QStringLiteral("label_34"));

        horizontalLayout_24->addWidget(label_34);

        prm_closest_obst_k = new QDoubleSpinBox(verticalLayoutWidget_5);
        prm_closest_obst_k->setObjectName(QStringLiteral("prm_closest_obst_k"));
        sizePolicy1.setHeightForWidth(prm_closest_obst_k->sizePolicy().hasHeightForWidth());
        prm_closest_obst_k->setSizePolicy(sizePolicy1);
        prm_closest_obst_k->setDecimals(0);
        prm_closest_obst_k->setMaximum(1e+06);

        horizontalLayout_24->addWidget(prm_closest_obst_k);


        verticalLayout_5->addLayout(horizontalLayout_24);

        verticalLayoutWidget_8 = new QWidget(centralWidget);
        verticalLayoutWidget_8->setObjectName(QStringLiteral("verticalLayoutWidget_8"));
        verticalLayoutWidget_8->setGeometry(QRect(480, 540, 221, 124));
        verticalLayout_8 = new QVBoxLayout(verticalLayoutWidget_8);
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setContentsMargins(11, 11, 11, 11);
        verticalLayout_8->setObjectName(QStringLiteral("verticalLayout_8"));
        verticalLayout_8->setContentsMargins(0, 0, 0, 0);
        label_25 = new QLabel(verticalLayoutWidget_8);
        label_25->setObjectName(QStringLiteral("label_25"));

        verticalLayout_8->addWidget(label_25);

        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setSpacing(6);
        horizontalLayout_18->setObjectName(QStringLiteral("horizontalLayout_18"));
        label_26 = new QLabel(verticalLayoutWidget_8);
        label_26->setObjectName(QStringLiteral("label_26"));

        horizontalLayout_18->addWidget(label_26);

        rrt_step_size = new QDoubleSpinBox(verticalLayoutWidget_8);
        rrt_step_size->setObjectName(QStringLiteral("rrt_step_size"));
        rrt_step_size->setDecimals(4);
        rrt_step_size->setMaximum(1000);

        horizontalLayout_18->addWidget(rrt_step_size);


        verticalLayout_8->addLayout(horizontalLayout_18);

        horizontalLayout_21 = new QHBoxLayout();
        horizontalLayout_21->setSpacing(6);
        horizontalLayout_21->setObjectName(QStringLiteral("horizontalLayout_21"));
        label_29 = new QLabel(verticalLayoutWidget_8);
        label_29->setObjectName(QStringLiteral("label_29"));

        horizontalLayout_21->addWidget(label_29);

        rrt_bias = new QDoubleSpinBox(verticalLayoutWidget_8);
        rrt_bias->setObjectName(QStringLiteral("rrt_bias"));
        rrt_bias->setDecimals(4);
        rrt_bias->setMaximum(1000);

        horizontalLayout_21->addWidget(rrt_bias);


        verticalLayout_8->addLayout(horizontalLayout_21);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setSpacing(6);
        horizontalLayout_19->setObjectName(QStringLiteral("horizontalLayout_19"));
        label_27 = new QLabel(verticalLayoutWidget_8);
        label_27->setObjectName(QStringLiteral("label_27"));

        horizontalLayout_19->addWidget(label_27);

        rrt_close_to_goal = new QDoubleSpinBox(verticalLayoutWidget_8);
        rrt_close_to_goal->setObjectName(QStringLiteral("rrt_close_to_goal"));
        rrt_close_to_goal->setDecimals(4);
        rrt_close_to_goal->setMaximum(1000);

        horizontalLayout_19->addWidget(rrt_close_to_goal);


        verticalLayout_8->addLayout(horizontalLayout_19);

        layoutWidget_2 = new QWidget(centralWidget);
        layoutWidget_2->setObjectName(QStringLiteral("layoutWidget_2"));
        layoutWidget_2->setGeometry(QRect(30, 540, 214, 36));
        horizontalLayout_22 = new QHBoxLayout(layoutWidget_2);
        horizontalLayout_22->setSpacing(6);
        horizontalLayout_22->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_22->setObjectName(QStringLiteral("horizontalLayout_22"));
        horizontalLayout_22->setContentsMargins(0, 0, 0, 0);
        label_30 = new QLabel(layoutWidget_2);
        label_30->setObjectName(QStringLiteral("label_30"));
        QSizePolicy sizePolicy2(QSizePolicy::Maximum, QSizePolicy::Minimum);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_30->sizePolicy().hasHeightForWidth());
        label_30->setSizePolicy(sizePolicy2);

        horizontalLayout_22->addWidget(label_30);

        max_sample = new QDoubleSpinBox(layoutWidget_2);
        max_sample->setObjectName(QStringLiteral("max_sample"));
        QSizePolicy sizePolicy3(QSizePolicy::MinimumExpanding, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(max_sample->sizePolicy().hasHeightForWidth());
        max_sample->setSizePolicy(sizePolicy3);
        max_sample->setDecimals(0);
        max_sample->setMaximum(1e+07);

        horizontalLayout_22->addWidget(max_sample);

        prm_graph = new QCheckBox(centralWidget);
        prm_graph->setObjectName(QStringLiteral("prm_graph"));
        prm_graph->setGeometry(QRect(560, 410, 111, 20));
        rrt_graph = new QCheckBox(centralWidget);
        rrt_graph->setObjectName(QStringLiteral("rrt_graph"));
        rrt_graph->setGeometry(QRect(560, 490, 111, 20));
        textOutputTime = new QTextBrowser(centralWidget);
        textOutputTime->setObjectName(QStringLiteral("textOutputTime"));
        textOutputTime->setGeometry(QRect(550, 10, 171, 271));
        prm_graph_free = new QCheckBox(centralWidget);
        prm_graph_free->setObjectName(QStringLiteral("prm_graph_free"));
        prm_graph_free->setGeometry(QRect(650, 430, 61, 20));
        prm_graph_obst = new QCheckBox(centralWidget);
        prm_graph_obst->setObjectName(QStringLiteral("prm_graph_obst"));
        prm_graph_obst->setGeometry(QRect(650, 450, 61, 20));
        layoutWidget_3 = new QWidget(centralWidget);
        layoutWidget_3->setObjectName(QStringLiteral("layoutWidget_3"));
        layoutWidget_3->setGeometry(QRect(560, 360, 188, 36));
        horizontalLayout_23 = new QHBoxLayout(layoutWidget_3);
        horizontalLayout_23->setSpacing(6);
        horizontalLayout_23->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_23->setObjectName(QStringLiteral("horizontalLayout_23"));
        horizontalLayout_23->setContentsMargins(0, 0, 0, 0);
        label_31 = new QLabel(layoutWidget_3);
        label_31->setObjectName(QStringLiteral("label_31"));
        sizePolicy2.setHeightForWidth(label_31->sizePolicy().hasHeightForWidth());
        label_31->setSizePolicy(sizePolicy2);

        horizontalLayout_23->addWidget(label_31);

        timeout = new QDoubleSpinBox(layoutWidget_3);
        timeout->setObjectName(QStringLiteral("timeout"));
        sizePolicy3.setHeightForWidth(timeout->sizePolicy().hasHeightForWidth());
        timeout->setSizePolicy(sizePolicy3);
        timeout->setDecimals(0);
        timeout->setMaximum(1e+07);

        horizontalLayout_23->addWidget(timeout);

        label_33 = new QLabel(layoutWidget_3);
        label_33->setObjectName(QStringLiteral("label_33"));

        horizontalLayout_23->addWidget(label_33);

        prm_graph_edge = new QCheckBox(centralWidget);
        prm_graph_edge->setObjectName(QStringLiteral("prm_graph_edge"));
        prm_graph_edge->setGeometry(QRect(650, 470, 61, 20));
        trace = new QPushButton(centralWidget);
        trace->setObjectName(QStringLiteral("trace"));
        trace->setGeometry(QRect(670, 320, 61, 32));
        showFilledObstacles = new QPushButton(centralWidget);
        showFilledObstacles->setObjectName(QStringLiteral("showFilledObstacles"));
        showFilledObstacles->setGeometry(QRect(550, 320, 121, 32));
        prm_graph_mixed = new QCheckBox(centralWidget);
        prm_graph_mixed->setObjectName(QStringLiteral("prm_graph_mixed"));
        prm_graph_mixed->setGeometry(QRect(650, 410, 61, 20));
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1117, 22));
        MainWindow->setMenuBar(menuBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Sampling-based Path Planning", 0));
        run->setText(QApplication::translate("MainWindow", "Run", 0));
        label->setText(QApplication::translate("MainWindow", "<html><head/><body><p>Input File Name</p></body></html>", 0));
        label_2->setText(QApplication::translate("MainWindow", "<html><head/><body><p>.txt</p></body></html>", 0));
        label_3->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">Start Configuration</p></body></html>", 0));
        label_4->setText(QApplication::translate("MainWindow", "<html><head/><body><p>x:</p></body></html>", 0));
        label_5->setText(QApplication::translate("MainWindow", "<html><head/><body><p>y:</p></body></html>", 0));
        label_8->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">Goal Configuration</p></body></html>", 0));
        label_9->setText(QApplication::translate("MainWindow", "<html><head/><body><p>x:</p></body></html>", 0));
        label_10->setText(QApplication::translate("MainWindow", "<html><head/><body><p>y:</p></body></html>", 0));
        label_13->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">Robot Specs: </p></body></html>", 0));
        label_15->setText(QApplication::translate("MainWindow", "<html><head/><body><p>robot name: </p></body></html>", 0));
        label_14->setText(QApplication::translate("MainWindow", "<html><head/><body><p>R:</p></body></html>", 0));
        label_17->setText(QApplication::translate("MainWindow", "<html><head/><body><p>Seed:</p></body></html>", 0));
        prm->setText(QApplication::translate("MainWindow", "PRM", 0));
        toggle->setText(QApplication::translate("MainWindow", "Toggle PRM", 0));
        lazytoggle->setText(QApplication::translate("MainWindow", "Lazy Toggle PRM", 0));
        rrt->setText(QApplication::translate("MainWindow", "RRT", 0));
        exit->setText(QApplication::translate("MainWindow", "Exit", 0));
        show->setText(QApplication::translate("MainWindow", "Show", 0));
        label_18->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">Animation Speed</p></body></html>", 0));
        label_20->setText(QApplication::translate("MainWindow", "<html><head/><body><p>Example Name </p></body></html>", 0));
        pause->setText(QApplication::translate("MainWindow", "Pause", 0));
        replay->setText(QApplication::translate("MainWindow", "Replay", 0));
        label_16->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">PRM(Lazy, Toggle):</p></body></html>", 0));
        label_19->setText(QApplication::translate("MainWindow", "<html><head/><body><p>closest k (free):</p></body></html>", 0));
        label_34->setText(QApplication::translate("MainWindow", "<html><head/><body><p>closest k (obst):</p></body></html>", 0));
        label_25->setText(QApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">RRT:</p></body></html>", 0));
        label_26->setText(QApplication::translate("MainWindow", "<html><head/><body><p>step size:</p></body></html>", 0));
        label_29->setText(QApplication::translate("MainWindow", "<html><head/><body><p>bias:</p></body></html>", 0));
        label_27->setText(QApplication::translate("MainWindow", "<html><head/><body><p>close to goal:</p></body></html>", 0));
        label_30->setText(QApplication::translate("MainWindow", "<html><head/><body><p># sample:</p></body></html>", 0));
        prm_graph->setText(QApplication::translate("MainWindow", "prm graph", 0));
        rrt_graph->setText(QApplication::translate("MainWindow", "rrt graph", 0));
        prm_graph_free->setText(QApplication::translate("MainWindow", "free", 0));
        prm_graph_obst->setText(QApplication::translate("MainWindow", "obst", 0));
        label_31->setText(QApplication::translate("MainWindow", "<html><head/><body><p>timeout:</p></body></html>", 0));
        label_33->setText(QApplication::translate("MainWindow", "<html><head/><body><p>sec</p></body></html>", 0));
        prm_graph_edge->setText(QApplication::translate("MainWindow", "edge", 0));
        trace->setText(QApplication::translate("MainWindow", "trace", 0));
        showFilledObstacles->setText(QApplication::translate("MainWindow", "show obstacles", 0));
        prm_graph_mixed->setText(QApplication::translate("MainWindow", "mixed", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
