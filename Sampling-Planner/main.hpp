/*
 * main.h
 *
 *  Created on: Mar. 22 2018
 *      Author: Ching-Hsiang Hsu
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>

using namespace std;

#include "config/FileParameter.hpp"
#include "utils/Timer.hpp"
#include "planner/Planner.hpp"
#include "display/MainWindow.hpp"

FILE *fpt_debug;
Planner *planner;
FileParameter *f_prm;

#endif // MAIN_H_
