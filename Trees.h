
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

// Trees.h
// difinition of a Trees class
//             by whuang gzhou

#ifndef TREES_H
#define TREES_H

#include "wstring.h"
#include "wfile.h"
#include <cstring>
#include "TreeOPE.h"
#include "warray.cpp"
#include "hash.hh"
#include "hashfunc.hh"
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
#include "rspr.hpp"
#ifndef COMMAND_LINE_VERSION
    #include <QString>
    #include <QStringList>
    #include <QProgressDialog>
    #include <QProgressBar>
#endif
//#include "pBar.h";

enum Treefileformat{NEXUS, NEWICK};
enum ConsensusTree{MAJORITYTREE, STRICTTREE};

extern std::map<String, void ***> StrToDist;
//extern std::map<String, double ***> StrToDoubleDist;
//extern std::map<String, int ***> StrToIntDist;

class Trees{
public:
    Trees();
    ~Trees();
    void destructor();
    void deletetrees();
    void deleteBipartitionMatrix();
    void deleteConsensustree();

    void initialTrees(string fname);
    void Settreeroottype(bool isrt);
    void Settreeweighttype(bool iswt);
    void ReadTrees(); // newick nexus
    void WriteTrees(string &outfile, Treefileformat tf); // newick nexus
    void WriteConsensusTree(string &outfile, Treefileformat tf); // newick nexus
    string WriteTreesFilename(string fname, string type);
    string WriteConsensusTreeFilename(string fname, string type);

    void Printf(int *idx, int length);

    NEWICKTREE * &operator[](const int idx){return treeset[idx];};               // pick the (index + 1)-th element

//    const NEWICKTREE * &operator[](const int idx) const{return treeset[idx];};   // pick the (index + 1)-th element

    void Compute_Hash();

    void Compute_Bipart_Matrix();
    void Compute_Bipart_Matrix(std::map<String, String>& paras);

    string make_Bipart_Matrix_name(string fname, String format);
	string make_Bipart_Matrix_name(string fname);

#ifdef COMMAND_LINE_VERSION
    void Compute_Bipart_Covariance();

    // by hash table
    bool Compute_RF_dist_by_hash(bool ISWEIGHTED);

    bool Compute_Matching_dist();
    bool Compute_SPR_dist();

    bool compute_community_automatically(String str_matrix, int modelType, string highfreq, string lowfreq);

    bool compute_community_manually(String str_matrix, int modelType, Array<double> param1, Array<double> param2, string highfreq, string lowfreq);
#else
    void Compute_Bipart_Covariance();

    // by hash table
    bool Compute_RF_dist_by_hash(bool ISWEIGHTED);

    bool Compute_Matching_dist();
    bool Compute_SPR_dist();

    bool compute_community_automatically(String str_matrix, int modelType, string highfreq, string lowfreq);

    bool compute_community_manually(String str_matrix, int modelType, Array<double> param1, Array<double> param2, string highfreq, string lowfreq);
#endif

    // by bipartmatrix
    void Compute_Unweight_RF_dist_by_bipartmatrix(); // TODO...

    void pttree(struct Ptree *treeA, int node);
    int tree_mmdis(struct Ptree *tree1, struct Ptree *tree2, int num_leaf);
    void compute_matrix(int *r, int range, struct Ptree *tree1, struct Ptree *tree2);
    int **array_to_matrix (int* m, int rows, int cols);

    void Compute_Geo_dist();

    void Compute_Affinity_dist(String str_matrix, int type);

    void Community_Detection();

    unsigned int Get_n_trees(){return n_trees;}

    int Get_treecov_size(){return treecov_size;}

    void OutputBipartitionMatrix(std::ostream &output, SparseMatrixOutputType smtype = RCVLIST);

    void compute_numofbipart();

    void load_distfile(string fname);

    void load_coordinatefile(string stdfname);

    void load_affinityfile(string fname);

    void load_covariancefile(string fname);

    void delete_matrix(String str_matrix);

    string make_DISToutput_name(String str_matrix);

    void print_matrix(String str_matrix, string outfile);

	void print_matrix2(String str_matrix, string outfile);

    bool bipartmatrixIsexisting();

    bool covarianceMatrixIsexisting();

    bool treesAreexisting();

    bool consensusTreeIsexisting();

    bool compute_community_fixedlambda(String str_matrix, int modelType, double lambdapos, double lambdaneg, string highfreq, string lowfreq);

    bool compute_consensus_tree(ConsensusTree type, const char *listname = NULL);

//    void Compute_Cumulative(double** bipartFreq, int* bipartFreqIdx, int splits, int increments);

//    void Compute_Slide(double** bipartFreq, int* bipartFreqIdx, int splits, int increments);

#ifndef COMMAND_LINE_VERSION
    QString get_treefilename_without_path();
#endif

    string Print_selected_indices();

    string Print_selected_trees(Treefileformat tf); // newick nexus

    string WriteSelectedTreesFilename(string type); // newick nexus
    void WriteSelectedTrees(string &outfile, Treefileformat tf); // newick nexus


    const NEWICKTREE *get_tree(int idx);
    const NEWICKTREE *get_contree(int idx)
    {
        if(idx < consensustrees.get_length() && idx >= 0)
            return consensustrees[idx];
        return NULL;
    };

    double **GetdistURF(){return dist_URF;}
    double **GetdistRF(){return dist_RF;}
    int **GetdistMatch(){return dist_match;}
    int **GetdistSPR(){return dist_SPR;}
    double **GetdistGeo(){return dist_geo;}
    double **GetdistFile(){return dist_file;}
    double **GetcoordFile(){return coord_file;}
    int Get_filedistsize(){return file_distsize;}
    int Get_fileAffinitysize(){return affinityfile_size;}
    int Get_filedcoordinatesize(){return file_coordinatesize;}
    int Get_filedcoordinatedim(){return file_coordinatedim;}
    bool Get_isrooted(){return isrooted;}
    bool Get_isweighted(){return isweighted;}
    const LabelMap *Get_labelmap(){return &leaveslabelsmaps;}
    Array<int> *getidxlist(){return &idxlist;}
	unsigned int *get_bipartcount() { return bipart_count; };

    map<unsigned long long, Array<char> *> hash2bitstr;
    Array<int> selected_trees;
    void Get_community_info(double ** &info, int &length){info = com_info; length = com_info_col;}

    string set_commfilename(string str){commfilename = str; return commfilename;}

private:

    //any functions that may be used.
    template<class T>
    void delete_double_array(T ** (&arr), int n);
    template<class T>
    void delete_double_array(T *** arr, int n);
    template<class T>
    void print_double_array(T *** arr, int n, string outfile);
	template<class T>
	void print_double_array2(T *** arr, int n, string outfile);
    template<class T>
    void print_coordinate_matrix(T *** arr, int n, int m, string outfile);

    Treefileformat CheckTreeFormat();
    void Sort(unsigned long long *matrix_hv,
                           unsigned int *matrix_treeIdx,
                           double *matrix_weight, int &idx);  // unacceptably slow......

    double **Vec_multiply(const double* Vec1, const double* Vec2, int Unique_idx);

    //functions used in community detection.
    string create_temp_name(String str_matrx);
    void print_community_file(String Str_matrix, string outfile, double highfreq, double lowfreq);
    template<class T>
    void print_comm_array(T *** arr, int n, string outfile, bool arr_is_covariance, double highfreq, double lowfreq);
    string create_out_name(String str_matrix);
    string create_node_name(String str_matrix);
    string create_conf_name(String str_matrix);
    int read_conf(char* filename, int* &conf, int* &sign);
    void create_resolution(double lp, double ln, int nb_layers, int* sign, double* &lambda);
    string create_comm_name(String str_matrix);

    // functions for Matching distance
    void print_bipartitionofonetree(NEWICKNODE*currentnode, bool isrooted, int depth);
    void Get_bipartitionofonetree(NEWICKNODE*currentnode, bool isrooted, int depth, Array<Array<char> > &bitstrofatree, int &idx);

    LabelMap leaveslabelsmaps;
    bool isrooted;
    bool isweighted;
    string treesfilename;
    Treefileformat treeformat;
    unsigned n_trees;
    NEWICKTREE **treeset;
    HashRFMap vec_hashrf;

    Array<NEWICKTREE *> consensustrees;
    Array<string> consensuslist;

    // sparse representation of bipartition matrix
    SparseMatrix *sbipartmatrix;
    map<long, Array<char> > MapHashBitsString;

//    double **bipartcovariance;
    string commfilename;
    int treecov_size;
    double **treecov;
    int filecov_size;
    double **filecov;
    int *numberofbipartition;
    unsigned int *bipart_count;

    double **dist_URF;
    double **dist_RF;
    int **dist_match;
    int **dist_SPR;
    double **dist_geo;

    double **affi_Recip_URF;
    double **affi_Recip_RF;
    double **affi_Recip_match;
    double **affi_Recip_SPR;
    double **affi_Recip_geo;

    double **affi_Exp_URF;
    double **affi_Exp_RF;
    double **affi_Exp_match;
    double **affi_Exp_SPR;
    double **affi_Exp_geo;

    int affinityfile_size;
    double **fileaffinity;

    int file_distsize;
    double **dist_file;
    double **affi_Recip_file;
    double **affi_Exp_file;
    int file_coordinatesize;
    int file_coordinatedim;
    double **coord_file;
    Array<int> idxlist;

    //--- to add storages for communities
    int *covariance_freeid;
    int covariance_freeid_size;
    int *covariance_nonfree_id;
    int covariance_nonfree_id_size;
    double **com_info;
    int com_info_col;


//    double** bipartFreq;
//    int* bipartFreqIdx;
};

#endif

