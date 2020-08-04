#pragma once
//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1   
//# Copyright (C) 2010 Wen Huang
//# 
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details. 
//# http://www.gnu.org/copyleft/gpl.html 
//##########################################################################

// zdcommunity.h
//    community dection method
//          August/2020
//                    by zdeng

#include "wstring.h"
#include "wfile.h"
#include <cstring>
#include "TreeOPE.h"
#include "warray.cpp"
#include <iostream>
#undef max
#undef min
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#else
#include <cctype>
#endif
#include <bitset>
#include <cmath>
//#include <valarray>
#include "Sparse_matrix.h"
#include "slicer.h"
#include "greedy_louvain.h"
#include "wdef.h"
#include "queue"
#include "wmatrix.h"

String make_stdname(String s, std::map<String, String> &paras);

template<class T>
void print_comm_array(Matrix<T> &arr, int n, File &output, bool arr_is_covariance, double highfreq, double lowfreq, int &covariance_freeid_size, int &covariance_nonfree_id_size, int *covariance_freeid, int *covariance_nonfree_id);

int read_conf(char* filename, int* &conf, int* &sign);

void create_resolution(double lp, double ln, int nb_layers, int* sign, double* &lambda);

bool community_detection_automatically(Matrix<double> &mat, map<String, String> &paras);

