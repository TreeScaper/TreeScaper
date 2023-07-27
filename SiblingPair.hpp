#ifndef SiblingPair_hpp
#define SiblingPair_hpp

#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <set>
#include "SPRNode.hpp"

class SiblingPair{
public:
    SPRNode *a;
    SPRNode *c;
    int key;
    int key2;
    
    SiblingPair();
    SiblingPair(SPRNode *A, SPRNode *C);
    bool operator< (const SiblingPair &sp) const;
};

void find_sibling_pairs_set_hlpr(SPRNode *n,
                                 std::set<SiblingPair> *sibling_pairs);
void append_sibling_pairs_set(SPRNode *n,std::set<SiblingPair> *sibling_pairs);
std::set<SiblingPair> *find_sibling_pairs_set(SPRNode *n);
std::set<SiblingPair> *find_sibling_pairs_set(Forest *f);



#endif /* SiblingPair_hpp */
