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

#include "wstring.hpp"
#include <cstring>
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
// #include "Sparse_matrix.hpp"
#include "SpecMat.hpp"
// #include "slicer.h"
#include "greedy_louvain.h"
// #include "wdef.h"
#include "queue"
// #include "wmatrix.h"
#include "header_info.hpp"

void check_self_covariance(SpecMat::LowerTri<PRECISION>& adjacency, double highfrequence, double lowfrequence, int* covariance_freeid, int* covariance_nonfree_id);

Graph *SymmetricOneSliceAdjaceny2Graph(SpecMat::LowerTri<PRECISION> &adjacency, int *node_id, int node_size, double weight_threshold);

// int read_conf(char* filename, int* &conf, int* &sign);

void create_resolution(double lp, double ln, int nb_layers, int* sign, double* &lambda);

bool community_detection_automatically(SpecMat::LowerTri<PRECISION> &mat, map<String, String> &paras);

