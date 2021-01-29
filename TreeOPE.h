
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

// TreeOPE.h
// definition of a TreeOPE class
//             by whuang gzhou

#ifndef TREEOPE_H
#define TREEOPE_H

#include "wstring.h"
//#include <cstring>
//#include <stdlib.h>
//#include <stdio.h>
#undef max
#undef min
#include <sstream>
#include "label-map.hh"
#include "hash.hh"
#include "warray.h"
#include "stdlib.h"

extern "C" {
#include "hungarian.h"
}

struct Ptree
{
    int leaf_number;
    int *parent;
    int *lchild;
    int *rchild;
    int **edge; //vector representation for bipatitions
};

typedef struct newicknode
{
    int Nchildren;              /* number of children (0 for leaves) */
    char *label;                /* node label, can be null */
    double weight;              /* node weight */
    struct newicknode **child;  /* list of children */
    unsigned long long hv1;
    unsigned long long hv2;
    Array<char> *bitstr;
    struct newicknode *parent;
} NEWICKNODE;

typedef struct
{
    NEWICKNODE *root;
} NEWICKTREE;

//static NEWICKTREE *TREEDEBUG = NULL;

class TreeOPE{
public:

    static int newickmain(int argc, char **argv);

    static NEWICKTREE *loadnewicktree(char *fname, int *error);
    static NEWICKTREE *loadnewicktree2(FILE *fp, int *error);
    static NEWICKTREE *floadnewicktree(FILE *fp, int *error);
    static void killnewicktree(NEWICKTREE *tree);
    static void killnewicktreeBitstr(NEWICKTREE *tree);
    static char *makenewicklabel(const char *str);
    static void printnewicktree(const NEWICKTREE *tree);
    static void printnewicknode(NEWICKNODE *node, int indent);

    static void killnoder(NEWICKNODE *node);
    static void killnoderBitstr(NEWICKNODE *node);
    static NEWICKNODE *loadnode(FILE *fp, int *error);
    static int addchild(NEWICKNODE *parent, NEWICKNODE *child);
    static NEWICKNODE *loadleaf(FILE *fp, int *error);
    static char *readlabelandweight(FILE *fp, double *weight, int *error);
    static char *readlabel(FILE *fp);
    static char *readquotedlabel(FILE *fp);
    static char *readplainlabel(FILE *fp);
    static void skipspace(FILE *fp);
    static char *mystrdup(const char *str);
    static int mystrcount(const char *str, int ch);

    static NEWICKTREE *parsetree(char *str, int *error, NEWICKTREE *testtree);
    static NEWICKNODE *parsenode(char **str, int *error);
    static NEWICKNODE *parseleaf(char **str, int *error);
    static char *parselabelandweight(char **str, double *weight, int *error);
    static char *parselabel(char **str);
    static char *parsequotedlabel(char **str);
    static char *parseplainlabel(char **str);
    static void skipstrspace(char *str);
    static void printTree_nex(NEWICKNODE *node, int length, char *gzbuff, int depth, int N);
    static void printTree_new(NEWICKNODE *node, char **str, char *gzbuff, int depth, int N);


    static void dfs_compute_hash(NEWICKNODE* startNode,
                                LabelMap &lm,
                                HashRFMap &vec_hashrf,
                                unsigned treeIdx,
                                unsigned &numBitstr,
                                unsigned long long m1,
                                unsigned long long m2,
                                bool WEIGHTED,
                                unsigned int NUM_Taxa,
                                map<unsigned long long, Array<char> *> &hash2bitstr,
                                int numofbipartition,
                                std::ostream& file_collusion,
                                int &collusion_cnt);

    static void GetTaxaLabels(NEWICKNODE *node, LabelMap &lm);

    static void bipart(NEWICKNODE * const startnode, unsigned int &treeIdx,
                             unsigned long long *matrix_hv,
                             unsigned int *matrix_treeIdx,
                             double *matrix_weight, int &idx, int depth, bool isrooted);

    static NEWICKNODE *findleaf(std::string leafname, NEWICKNODE *currentnode, NEWICKNODE *parent, int *icpt);
    static void normailzedTree(NEWICKNODE *lrpt, NEWICKTREE *newickTree, int indexchild);
    static void normalizedNode(NEWICKNODE *currentnode, NEWICKNODE *parent, double currentweight);
    static void Label_strint(NEWICKNODE *node, LabelMap &lm);

    static bool newick2lcbb(const NEWICKTREE *nwtree, int num_leaves, struct Ptree *tree);
    static bool newick2ptree(NEWICKNODE *node, struct Ptree *tree, bool &currentnode, int &vidx, int &curidx);
    static int sumofdegree(NEWICKNODE *node, bool isrooted, int depth = -1);
    static void obtainbipartcount(NEWICKTREE *nwtree, bool isrooted, map<unsigned long long, unsigned long long> &bipcount);

    static void bipartcount(NEWICKNODE *node, bool isrooted, map<unsigned long long, unsigned long long> &bipcount, int depth = 0);

    static bool buildconsensustree(NEWICKTREE *&tree, const Array<double> &confreq, const Array<unsigned long long> &conhash, Array<char> *contreebitstr, unsigned int conbipnum, int bitlength);

    static bool Addbipart(NEWICKNODE* startNode, double freq, unsigned long long hash, Array<char> &bitstr, int NumTaxa, bool &iscontained);
};

#endif
