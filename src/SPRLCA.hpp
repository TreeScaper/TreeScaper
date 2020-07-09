#ifndef SPRLCA_hpp
#define SPRLCA_hpp

#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include "SPRNode.hpp"

int mylog2(int val);
class SPRNode;

class LCA{
private:
    SPRNode *SPR_tree;
    std::vector<int> E;		// preorder numbers of euler tour
    std::vector<int> L;		// levels of euler tour
    std::vector<int> H;		// first occurence of a preorder number in E
    std::vector<int> T;    // real preorder to internal preorder mapping
    std::vector<SPRNode *> N;	// preorder to SPRNode mapping
    std::vector<std::vector<int> > RMQ;	// precomputed RMQ values
    
public:
    LCA(SPRNode *SPR_tree);
    LCA();

    void euler_tour(SPRNode *node, int SPR_depth);
    void precompute_rmq();
    int rmq(int i, int j);
    SPRNode *get_lca(SPRNode *a, SPRNode *b);
    SPRNode *get_SPR_tree();

};
#endif /* SPRLCA_hpp */
