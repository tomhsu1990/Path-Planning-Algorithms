#
#  sampling-planner.pro
#
#  Created on: Mar. 22, 2018
#     Author: Ching-Hsiang Hsu

QT       += \
core \
gui \
opengl

greaterThan(QT_MAJOR_VERSION, 5): QT += widgets

TARGET = sampling-planner

TEMPLATE = app

CONFIG += c++11

SOURCES += \
    main.cpp \
    config/FileParameter.cpp \
    config/Config.cpp \
    config/Robot.cpp \
    config/Env.cpp \
    planner/Planner.cpp \
    planner/RRT.cpp \
    planner/PRM.cpp \
    planner/LazyPRM.cpp \
    planner/TogglePRM.cpp \
    planner/LazyTogglePRM.cpp \
    display/MainWindow.cpp \
    display/Display.cpp \
    utils/Triangulate.cpp

HEADERS  += \
    main.hpp \
    config/FileParameter.hpp \
    config/Config.hpp \
    config/Robot.hpp \
    config/Env.hpp \
    planner/Planner.hpp \
    planner/RRT.hpp \
    planner/PRM.hpp \
    planner/LazyPRM.hpp \
    planner/TogglePRM.hpp \
    planner/LazyTogglePRM.hpp \
    display/MainWindow.hpp \
    display/Display.hpp \
    utils/Triangulate.hpp \
    utils/Timer.hpp

# CORE geometry library source code
SOURCES += \
    ../lib/CORE/MpfrIO.cpp \
    ../lib/CORE/linearAlgebra.cpp \
    ../lib/CORE/geom2d/point2d.cpp \
    ../lib/CORE/geom2d/line2d.cpp \
    ../lib/CORE/geom2d/circle2d.cpp \
    ../lib/CORE/geom2d/segment2d.cpp \
    ../lib/CORE/geom2d/triangle2d.cpp \
    ../lib/CORE/geom2d/polygon2d.cpp

HEADERS  += \
    ../lib/CORE/CORE.h \
    ../lib/CORE/linearAlgebra.h \
    ../lib/CORE/geom2d/point2d.h \
    ../lib/CORE/geom2d/line2d.h \
    ../lib/CORE/geom2d/circle2d.h \
    ../lib/CORE/geom2d/segment2d.h \
    ../lib/CORE/geom2d/triangle2d.h \
    ../lib/CORE/geom2d/polygon2d.h

FORMS += \
    display/MainWindow.ui

# Eigen
INCLUDEPATH += /usr/local/include/eigen3/

# ANN library
#unix: LIBS += -L$$PWD/../lib/ann_1.1.2/lib/ -lANN
#INCLUDEPATH += $$PWD/../lib/ann_1.1.2/include
#DEPENDPATH += $$PWD/../lib/ann_1.1.2/include
#unix: PRE_TARGETDEPS += $$PWD/../lib/ann_1.1.2/lib/libANN.a

# Flann library
unix: LIBS += -L$$PWD/../lib/flann_1.8.4/build/lib/ -lflann.1.8.4
INCLUDEPATH += $$PWD/../lib/flann_1.8.4/src/cpp/
DEPENDPATH += $$PWD/../lib/flann_1.8.4/src/cpp/


# Boost C++ graph library
unix: PRE_TARGETDEPS += $$PWD/../lib/boost_1_66_0/bin.v2/libs/graph/build/darwin-darwin-4.2.1/release/threadapi-pthread/threading-multi/libboost_graph.dylib
INCLUDEPATH += $$PWD/../lib/boost_1_66_0/
DEPENDPATH += $$PWD/../lib/boost_1_66_0/

# CORE geometry library
INCLUDEPATH += \
    /usr/local/Cellar/gmp/6.1.1/include/ \
    /usr/local/Cellar/mpfr/3.1.5/include/ \
    ../lib/
LIBS += \
    -L/usr/local/Cellar/gmp/6.1.1/lib/ -lgmp \
    -L/usr/local/Cellar/mpfr/3.1.5/lib/ -lmpfr
