#pragma once

#include <fstream>
#include <iostream>

#include "array.hpp"
#include "sparse.hpp"
#include "wstring.hpp"
#include "zdtree.hpp"
#include "header_info.hpp"

int get_file_lines(std::fstream &fs);

int get_file_cols(std::fstream &fs);

String make_stdname(String s, std::map<String, String> &paras);

String get_path(String fname);

bool read_subset(String fname, Array<size_t> &subset);

void read_paras_from_csv(String fname, map<String, String> &paras, bool allow_insert);

map<String, String> read_paras(int argc, char *argv[], int key_size, String *default_paras, String *options);

void binary_print_smat(std::ostream &output, SparseMatrix *sb2t_mat);

SparseMatrix *binary_read_smat(std::istream &input);

void binary_print_lowertri(std::ostream &output, SpecMat::LowerTri<PRECISION> *ltmat);

SpecMat::LowerTri<PRECISION> *binary_read_lowertri(std::istream &input);

void print_lowertri(std::ostream &output, SpecMat::LowerTri<PRECISION> *ltmat);

SpecMat::LowerTri<PRECISION> *read_lowertri(std::istream &input);

void binary_print_taxon(std::ostream &output, TaxonList *taxa);

TaxonList *binary_read_taxon(std::istream &input);

template <class T>
void binary_print_bitstr(std::ostream &output, BitString<T> &bitstr)
{
    size_t l, b_size;
    l = bitstr.get_length();
    b_size = bitstr.get_bitsize();

    T *v = bitstr.get_vec();

    output << l << b_size;
    for (size_t i = 0; i < l; i++)
        output << v[i];
}

template <class T>
BitString<T> *binary_read_bitstr(std::istream &input, T dummy)
{
    size_t l, b_size;
    input >> l >> b_size;

    T *v = new T[l];
    for (size_t i = 0; i < l; i++)
        input >> v[i];

    BitString<T> *Ans = new BitString<T>(l, b_size, v);
    return Ans;
}

template <class T>
void binary_print_bipart(std::ostream &output, Bipartition<T> *Bipart)
{
    T hb, inv_h;
    size_t n;
    bool iss;

    Array2D<int> *h2i;
    Array<int> *h2i_col;
    int *h2i_col_vec;
    size_t h2i_size, h2i_col_size;

    hb = Bipart->get_hash_bound();
    inv_h = Bipart->get_invariant_hashing();
    n = Bipart->get_size();
    h2i = Bipart->get_hash_table();
    h2i_size = h2i->get_size();

    output << hb << inv_h;
    output << h2i_size;
    for (size_t i = 0; i < h2i_size; i++)
    {
        h2i_col = h2i->get_c_ptr(i);
        h2i_col_size = (h2i_col == nullptr) ? 0 : h2i_col->get_size();
        h2i_col_vec = (h2i_col == nullptr) ? nullptr : h2i_col->get_vec();
        output << h2i_col_size;
        for (size_t j = 0; j < h2i_col_size; j++)
        {
            output << h2i_col_vec[j];
        }
    }

    output << n;
    for (size_t i = 0; i < n; i++)
        binary_print_bitstr(output, &Bipart[i]);

    iss = Bipart->get_issorted();
    output << iss;
}

template <class T>
Bipartition<T> *binary_read_bipart(std::istream &input, T dummy, TaxonList *Taxa)
{
    T hb, inv_h;
    size_t n;
    int ind;
    bool iss;

    size_t h2i_size, h2i_col_size;

    input >> hb >> inv_h >> h2i_size;
    Array2D<int> *h2i = new Array2D<int>(hb);
    for (size_t i = 0; i < hb; i++)
    {
        input >> h2i_col_size;
        if (h2i_col_size != 0)
            h2i->set_c(i, new Array<int>(0, 50));
        for (size_t j = 0; j < h2i_col_size; j++)
        {
            input >> ind;
            h2i->push(ind, i);
        }
    }

    Bipartition<T> *Ans = new Bipartition<T>();
    Ans->set_taxon_list(Taxa);
    Ans->set_hash_bound(hb);
    Ans->set_invariant(inv_h);
    Ans->set_hash_table(h2i);

    input >> n;
    for (size_t i = 0; i < n; i++)
        Ans->push_bs_only(*binary_read_bitstr(input, dummy));

    return Ans;
}
