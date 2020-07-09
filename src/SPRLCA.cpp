#include "SPRLCA.hpp"

#ifdef __GNUC__
#define clz(x) __builtin_clz(x)
#define ctz(x) __builtin_ctz(x)
#else
static int popcnt( int x )
{
    x -= ((x >> 1) & 0x55555555);
    x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
    x = (((x >> 4) + x) & 0x0f0f0f0f);
    x += (x >> 8);
    x += (x >> 16);
    return x & 0x0000003f;
}
static int clz( int x )
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return 32 - popcnt(x);
}
static int ctz( int x )
{
    return popcnt((x & -x) - 1);
}

#endif

int mylog2(int val)
{
    if (val <= 0)
        return -1;
    return 31 - clz(val);
    int ret = -1;
    while (val > 0) {
        val >>= 1;
        ret++;
    }
    std::cout << ret << std::endl;
    std::cout << std::endl;
    return ret;
}

LCA::LCA(SPRNode *SPR_tree) {
    this->SPR_tree = SPR_tree;
    T = std::vector<int>();
    if (SPR_tree->get_preorder_number() == -1)
        SPR_tree->preorder_number();
    euler_tour(SPR_tree, 0);
    precompute_rmq();
}

LCA::LCA() {
    this->SPR_tree = NULL;
}

void LCA::euler_tour(SPRNode *node, int SPR_depth) {
    // First visit
    int preorder_number = (int) N.size();
    int euler_number = (int) E.size();
    N.push_back(node);
    if (T.size() <= node->get_preorder_number())
        T.resize(node->get_preorder_number()+1,-1);
    T[node->get_preorder_number()] = preorder_number;
    
    //cout << preorder_number << "\t";
    //SPRNode->print_subSPR_tree();
    H.push_back(euler_number);
    L.push_back(SPR_depth);
    E.push_back(preorder_number);
    
    // TODO: check that this is correct for multifurcating SPR_trees
    
    std::list<SPRNode *>::const_iterator c;
    for(c = node->get_children().begin(); c != node->get_children().end();
        c++) {
        euler_tour(*c, SPR_depth+1);
        // Middle/Last visit
        L.push_back(SPR_depth);
        E.push_back(preorder_number);
    }
}

void LCA::precompute_rmq() {
    // E gives queries of length 1
    RMQ.push_back(E);
    for(int length = 1; length < E.size(); length *= 2) {
        std::vector<int> V = std::vector<int>();
        for(int start = 0; start < E.size() - length; start++) {
            V.push_back(rmq(start, start + length));
        }
        RMQ.push_back(V);
    }
}

int LCA::rmq(int i, int j) {
    if (i == j)
        return E[i];
    int length = j - i - 1;
    length = (mylog2(length));
    int rmq1 = RMQ[length+1][i];
    int rmq2;
    if (length >= 0)
        rmq2 = RMQ[length+1][j - (1 << (length))];
    else
        rmq2 = RMQ[length+1][j];
    if (rmq1 < rmq2)
        return rmq1;
    return rmq2;
}

SPRNode *LCA::get_lca(SPRNode *a, SPRNode *b) {
    int preorder_a = T[a->get_preorder_number()];
    int preorder_b = T[b->get_preorder_number()];
    int lca_index;
    if (preorder_a <= preorder_b)
        lca_index = rmq(H[preorder_a], H[preorder_b]);
    else
        lca_index = rmq(H[preorder_b], H[preorder_a]);
    return N[lca_index];
}

SPRNode *LCA::get_SPR_tree() {
    return SPR_tree;
}


