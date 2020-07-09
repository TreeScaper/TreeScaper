
#include "SiblingPair.hpp"

SiblingPair::SiblingPair()
{
    a = NULL;
    c = NULL;
    key = -1;
    key2 = -1;
}

SiblingPair::SiblingPair(SPRNode *A, SPRNode *C)
{
    a = A;
    c = C;
    key = a->get_preorder_number();
    if (c->get_preorder_number() < key)
        key = c->get_preorder_number();
}

bool SiblingPair::operator< (const SiblingPair &sp) const
{
    return key < sp.key;
}

void find_sibling_pairs_set_hlpr(SPRNode *n,
                                 std::set<SiblingPair> *sibling_pairs) {
    SPRNode *lchild = n->lchild();
    SPRNode *rchild = n->rchild();
    bool lchild_leaf = false;
    bool rchild_leaf = false;
    if (lchild != NULL) {
        if (lchild->is_leaf())
            lchild_leaf = true;
        else
            find_sibling_pairs_set_hlpr(lchild,sibling_pairs);
    }
    if (rchild != NULL) {
        if (rchild->is_leaf())
            rchild_leaf = true;
        else
            find_sibling_pairs_set_hlpr(rchild,sibling_pairs);
    }
    if (lchild_leaf && rchild_leaf) {
        sibling_pairs->insert(SiblingPair(lchild,rchild));
        //lchild->add_to_sibling_pairs(sibling_pairs, 1);
        //rchild->add_to_sibling_pairs(sibling_pairs, 2);
    }
}

// find the sibling pairs in this SPRNode's subtree
void append_sibling_pairs_set(SPRNode *n,std::set<SiblingPair> *sibling_pairs) {
    find_sibling_pairs_set_hlpr(n,sibling_pairs);
}

// find the sibling pairs in this SPRNode's subtree
std::set<SiblingPair> *find_sibling_pairs_set(SPRNode *n) {
    std::set<SiblingPair> *sibling_pairs = new std::set<SiblingPair>();
    find_sibling_pairs_set_hlpr(n,sibling_pairs);
    return sibling_pairs;
}

// return a set of the sibling pairs
std::set<SiblingPair> *find_sibling_pairs_set(Forest *f) {
    std::set<SiblingPair> *sibling_pairs = new std::set<SiblingPair>();
    for(int i = 0; i < f->num_components(); i++) {
        SPRNode *component = f->get_component(i);
        append_sibling_pairs_set(component,sibling_pairs);
    }
    return sibling_pairs;
}


