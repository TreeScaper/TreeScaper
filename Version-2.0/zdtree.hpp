#pragma once

#include <stdlib.h>
#include <cstring>
#include <map>
#include <fstream>
#include <iostream>
#include "zdarray.hpp"
#include "zdarray.cpp"
#include "Sparse_matrix.hpp"

using std::string;
using std::map;
using std::ostream;


// This is a simple translation from Newwick formatted string to an array of tree.
// It contains rows of integers. Each row represents all index of child-nodes of 
// some internal node except for the first integer, which is the row index of the
// parent node. 

template <class T>
T GetPrime(T topNum, unsigned from)
{
	T primeNum = 0;
	T candidate = 0;

	if (topNum <= 100)
		candidate = 2;
	else
		candidate = topNum;

	while (candidate <= topNum + from) {
		T trialDivisor = 2;
		int prime = 1;

		while (trialDivisor * trialDivisor <= candidate) {
			if (candidate % trialDivisor == 0) {
				prime = 0;
				break;
			}
			trialDivisor++;
		}
		if (prime) primeNum = candidate;

		candidate++;
	}

	return primeNum;
};

class TaxonList {
public:
	TaxonList() : size(0), bitstr_size(0){};
	int size;
	int bitstr_size;
	Array<string> Ind2Taxon;
	map<string, int> Taxon2Ind;
	Array<int> IndB2IndA;
	
	void push(string item, bool flag_str_format = true);
	template <class T>
	void set_bitstr_size(T ele) { 
		bitstr_size = (size + (-size) % (8 * sizeof(T))) / (8 * sizeof(T));
	};
	int get_bitstr_size() const { return bitstr_size; };
	size_t ReadTaxa(std::string fname);
};


template <class T>
class BitString {
	// This is a bitstring representing the bipartition. It must normalized with 0 at the first bit.
	// T must be unsigned.
public:
	friend void swap(BitString<T>& lhs, BitString<T>& rhs) {
		using std::swap;
		swap(lhs.length, rhs.length);
		swap(lhs.bit_size, rhs.bit_size);
		swap(lhs.vec, rhs.vec);
	}

	BitString() : length(0),vec(nullptr) {};
	BitString(int n, int b_size) : length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr) { memset(vec, 0, length * sizeof(T)); };
	BitString(int n, int b_size, int index) : length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr) {
		if (n > 0) {
			memset(vec, 0, length * sizeof(T));
			int pos_byte = index / (8 * sizeof(T));
			int pos_bit = index % (8 * sizeof(T));
			T bit = 1;
			vec[pos_byte] |= bit << pos_bit;
		}
	}
	BitString(int n, int b_size, T* src) :length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr) {
		memcpy(vec, src, length * sizeof(T));
	}
	BitString(const BitString<T>& src): length(src.length), bit_size(src.bit_size), vec(length > 0 ? new T[length] : nullptr) {
		memcpy(vec, src.vec, length * sizeof(T));
	}


	~BitString() { delete[] vec; };

	BitString<T> operator|(const BitString<T>& rhs) {
		T* temp = new T[length];
		for (int i = 0; i < length; i++)
			temp[i] = vec[i] | rhs.vec[i];
		return BitString<T>(length, temp);
	}

	BitString<T>& operator|=(const BitString<T>& rhs) {
		for (int i = 0; i < length; i++)
			this->vec[i] |= rhs.vec[i];
		return *this;
	}

	BitString<T>& operator=(const BitString<T>& rhs) {
		BitString<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	bool operator==(const BitString<T>& rhs) {
		if (length != rhs.length || bit_size != rhs.bit_size) {
			return false;
		}
		T temp = vec[0] ^ rhs.vec[0];

		if (temp & (T)1) {
			for (int i = 0; i < length - 1; i++) {
				if (temp != (T)-1)
					return false;
				temp = vec[i + 1] ^ rhs.vec[i + 1];
			}
			size_t remainder = 8 * length * sizeof(T) - bit_size;
			T mask = ((T)-1 >> remainder);
			if ((temp & mask) != mask)
				return false;
			return true;
		}
		else {
			for (int i = 0; i < length; i++) {
				if (vec[i] != rhs.vec[i])
					return false;
			}
			return true;
		}
	}

	bool operator<(const BitString<T>& rhs) { // This only works for normalized bitstring
		T mask;

		for (int i = 0; i < length; i++) {
			mask = vec[i] ^ rhs.vec[i];			// mask for all different bits
			if (mask == 0)						// no different bit found
				continue;
			else{
				mask = (mask & -mask);			// mask for the rightmost bit
				if (vec[i] & mask)
					return false;
			}
		}
		if (mask == 0)
			return false;
		else
			return true;
	}

	friend bool operator<(const BitString<T>& lhs, const BitString<T>& rhs) {
		T mask = 0;

		for (int i = 0; i < lhs.length; i++) {
			mask = lhs.vec[i] ^ rhs.vec[i];			// mask for all different bits
			if (mask == 0)							// no different bit found
				continue;
			else {
				mask = (mask & (~mask + 1));		// mask for the rightmost bit
				if (lhs.vec[i] & mask)
					return false;
			}
		}

		if (mask == 0)
			return false;
		else
			return true;
	}

	bool operator!=(const BitString<T>& rhs) {
		if (length != rhs.length) {
			return true;
		}
		else {
			for (int i = 0; i < length; i++) {
				if (vec[i] != rhs.vec[i])
					return true;
			}
			return false;
		}
	}

	T& operator[](size_t i) { return vec[i]; };

	T operator()(size_t i) const { return vec[i]; };

	int size() const { return length; };

	BitString<T>& normalized() {
		if (!this->is_normalized()) { // First bit is 1
				for (int i = 0; i < length; i++)
					vec[i] = ~vec[i];
				size_t remainder = 8 * length * sizeof(T) - bit_size;
				T temp = -1;
				temp = (temp >> remainder);
				vec[length - 1] = vec[length - 1] & temp;
		}
		return *this;
	}

	bool is_normalized() const {
		//assert(this->length != 0);
		return !((bool) (vec[0] & (T) 1));
	}

	void print_BitString(ostream &fout) {
		bool bit;
		size_t n = sizeof(T) * 8;
		T one = 1, temp = 0;
		for (int i = 0; i < length - 1; i++) {
			temp = vec[i];
			for (int j = 0; j < n; j++) {
				bit = temp & one;
				temp = temp >> 1;
				fout << bit;
			}
		}
		temp = vec[length - 1];
		n = bit_size - (length - 1) * sizeof(T) * 8;
		for (int j = 0; j < n; j++) {
			bit = temp & one;
			temp = temp >> 1;
			fout << bit;
		}
	}

	T* get_vec() { return vec; };

private:
	size_t length;
	size_t bit_size;
	T* vec;
};


template <class T>
class Bipartition {
public:
	Bipartition() : Taxa(nullptr), hash_bound(0), 
		Id2BitString(Array<BitString<T> >(0, 1000)), 
		Hash2Id(nullptr), Id2Tree(0, 500) {};
	Bipartition(TaxonList* taxa, int n_tree, T hb) : Taxa(taxa), hash_bound(hb),
		Id2BitString(Array<BitString<T> >(0, 1000)),
		Hash2Id(new Array<int>*[hash_bound]), Id2Tree(0, 500) {
		this->set_invariant();
		for (int i = 0; i < Taxa->size; i++)
			this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i), 0);
	};

	Bipartition(TaxonList* taxa, int n_tree) : Taxa(taxa), Id2BitString(Array<BitString<T> >(0, 1000)),
		Hash2Id(nullptr), Id2Tree(0, 500) {
		this->set_hash_bound(n_tree);
		this->set_invariant();
		
		for (int i = 0; i < Taxa->size; i++)
			this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i).normalized(), 0);
	};

	bool is_empty() { return Id2BitString.is_empty(); };
	

	int push(BitString<T>& bs, int tree_id, int index_hash = -1) {
		if (index_hash == -1)
			index_hash = hashing(bs);
		int index = -1;
		if (Hash2Id[index_hash] == nullptr) {
			Hash2Id[index_hash] = new Array<int>(0, 10);
			Hash2Id[index_hash]->push(Id2BitString.get_size());
			Id2BitString.push(bs);
			Id2BitString.back().normalized();

			Id2Tree.push(Array<int>(0, 10));
			Id2Tree.back().push(tree_id);

			return Id2BitString.get_size() - 1;
		}
		else {
			for (int i = 0; i < Hash2Id[index_hash]->get_size(); i++)
			{
				index = Hash2Id[index_hash]->operator[](i);
				if (Id2BitString[index] == bs){
					Id2Tree[index].push(tree_id);
					return index;
				}
			}
			Hash2Id[index_hash]->push(Id2BitString.get_size());
			Id2BitString.push(bs);
			Id2BitString.back().normalized();

			Id2Tree.push(Array<int>(0, 10));
			Id2Tree.back().push(tree_id);

			return Id2BitString.get_size() - 1;
		}
	};

	// Add BitString to the list and return the id assigned to it. If it has encountered before, do nothing.


	T hashing(const BitString<T>& bs) const {
		T result = 0;
		for (int i = 0; i < bs.size(); i++) {
			result += bs(i) % hash_bound;
			result %= hash_bound;
		}
		if (bs.is_normalized())		//Is there a better way to access the top bit?
			return result;
		else
			return (result > invariant_hashing ? hash_bound + invariant_hashing - result : invariant_hashing - result);
	}

	void set_hash_bound(int n_tree) {
		T p = 0;
		unsigned int mul = 1;
		T top = (n_tree * (2 * Taxa->size - 3)) / 10;
		do {
			unsigned from = 100 * mul;
			p = GetPrime(top, from);
			++mul;
		} while (p == 0);
		hash_bound = (p < 101 ? 101 : p);
		Hash2Id = new Array<int> * [hash_bound];

		for (int i = 0; i < hash_bound; i++)
			Hash2Id[i] = new Array<int>(0, 10);		
	}

	void set_invariant() {
		//compute the hash value of the sum of normalized bitstring and its complement.
		// example: taxa size: 6, bit size of T: 4
		// bitstring:	0010 0001
		// complement:	1101 0010
		// sum:			1111 0011
		T result = 0, comp = -1;
		int bitstr_size = Taxa->bitstr_size;
		int remainder = bitstr_size * (8 * sizeof(T)) - Taxa->size;

		for (int i = 0; i < bitstr_size - 1; i++)
		{
			result += comp % hash_bound;
			result %= hash_bound;
		}
		comp = comp >> remainder; //fill zeros in the last container.
		result += comp % hash_bound;
		result %= hash_bound;
		invariant_hashing = result;
	}

	int get_leaf_size() const { return Taxa->size; };
	int get_bitstr_size() const { return Taxa->bitstr_size; };
	int get_size() const { return Id2BitString.get_size(); };
	T get_hash_bound() const { return hash_bound; };
	T get_invariant_hashing() const { return invariant_hashing; };
	Array<Array<int> >* get_Id2Tree_ptr() { return &(Id2Tree); };


	void print_summary(int n_trees) {
		int max_collision_size = 0, collision_cnt_5 = 0, collision_cnt_10 = 0, empty_cnt = 0;
		int cur_size;
		for (int i = 0; i < hash_bound; i++) {
			if (Hash2Id[i] != nullptr) {
				cur_size = Hash2Id[i]->get_size();
				if (cur_size > max_collision_size)
					max_collision_size = cur_size;
				if (cur_size >= 5)
					collision_cnt_5++;
				if (cur_size >= 10)
					collision_cnt_10++;
				if (cur_size == 0)
					empty_cnt++;
				//std::cout << "Collusion beam # " << i << ":\t";
				//for (int j = 0; j < cur_size; j++)
				//	std::cout << Hash2Id[i]->operator[](j) <<"  ";
				//std::cout << '\n';
			}
			else {
				empty_cnt++;
			}
		}
		std::cout << "-------------Bipartition info summary------------------\n";
		std::cout << "\tDistinct bipartition number: " << Id2BitString.get_size() << ",\n";
		std::cout << "\tHash container size: " << hash_bound << ",\n";
		std::cout << "\tEmpty container: " << empty_cnt << ",\n";
		std::cout << "\tMaximum Collision beam size: " << max_collision_size << ",\n";
		std::cout << "\tCount of Collision beams over 5: " << collision_cnt_5 << ",\n";
		std::cout << "\tCount of Collision beams over 10: " << collision_cnt_10 << ".\n";
		std::cout << "-------------Bipartition info summary------------------\n";
	}

	void print_Bipart(ostream& fout) {
		size_t n = Id2BitString.get_size();
		//fout << Taxa->size << '\t' << n << '\n';
		for (int i = 0; i < n; i++){
			fout << i << " ";
			Id2BitString[i].print_BitString(fout);
			fout << ' ' << Id2Tree[i].get_size() << '\n';
		}
	}

	void print_Bipart(ostream& fout, size_t i) {
		Id2BitString[i].print_BitString(fout);
	}

	BitString<T>& operator[](size_t i) { return Id2BitString[i]; };


private:
	TaxonList* Taxa;
	T hash_bound;
	Array<BitString<T> > Id2BitString;
	Array<int>** Hash2Id;
	T invariant_hashing;
	map<BitString<T>, int> BitStr2Id;
	Array<Array<int> > Id2Tree;
};

class TreeArray {
	// This class use edges forms as internal data structure of trees. It translates Newick form to edge form 
	// by labelling all internal node with integers. The labelling tends to assign small integer to nodes that
	// are close to leaf. The leaf nodes are labelled with (0, ..., l-1) determined by the taxa normalization.
	// Note that the labelling depends on the Newick form, which is not unique, i.e., different Newick form
	// of a same tree gives different edge array. 
public:
	TreeArray() : size(0), levels(Array<Array<int> >(0, 100)), weights(Array<double>(0, 100)), edges(nullptr) {};
	~TreeArray() {};

	TreeArray(string& tree, TaxonList& taxon_list, Array<int>& active_rows, Array<int>& level_flag, int flag_label);

	void print() {
		for (int i = 0; i < levels.get_size(); i++) {
			for (int j = 0; j < levels[i].get_size(); j++)
				std::cout << levels[i][j] << '\t';
			std::cout << '\n';
		}


		for (int i = 0; i < size; i++)
			std::cout << edges[0][i] << '\t' << edges[1][i] << '\t' << weights[i] << '\n';
		std::cout << "\n\n\n";
	}


	void label_internal_node(Array<int>& active_rows, Array<int>& unlabeled);
	void label_internal_node(Array<int>& active_rows, Array<int>& unlabeled, Array<Array<double> >& w_temp);

	template <class T>
	void compute_bitstring(Bipartition<T>& Bipart, int tree_id) {
		bitstrs = new int[size];
		T* index_hash = new T[size];
		memset(index_hash, 0, size * sizeof(T));
		BitString<T>* bs = new BitString<T>[size];
		T iv_h = Bipart.get_invariant_hashing();
		T h_b = Bipart.get_hash_bound();
		int l_s = Bipart.get_leaf_size();


		for (int i = 0; i < size; i++) {
			if (i < l_s) {
				bs[i] = BitString<T>(Bipart.get_bitstr_size(), l_s, i);
				index_hash[i] = Bipart.hashing(bs[i]);
			}
			else
				bs[i] = BitString<T>(Bipart.get_bitstr_size(), l_s);
		}
		// Note that bs[0] = 000000000...000001 = 1 is not normalized.
		// Therefore hashing() will return the hash of its complement's hash.
		// However, we need 1 mod M = 1 here.

		index_hash[0] = 1;

		for (int i = 0; i < size; i++) {
			if (edges[1][i] < l_s) {
				bitstrs[i] = edges[1][i];
			}
			else{
				T hv = index_hash[edges[1][i]];
				if (bs[edges[1][i]].is_normalized())
					bitstrs[i] = Bipart.push(bs[edges[1][i]], tree_id, hv);
				// use the computed hashing value for normalized bitstring generated from bitwise OR.
				else
					bitstrs[i] = Bipart.push(bs[edges[1][i]], tree_id, (hv > iv_h ? h_b + iv_h - hv : iv_h - hv));
				// use the hashing of complement bitstring.
			}
			
			if (edges[0][i] < size) {
				bs[edges[0][i]] |= bs[edges[1][i]];
				index_hash[edges[0][i]] += index_hash[edges[1][i]];
				index_hash[edges[0][i]] = index_hash[edges[0][i]] % Bipart.get_hash_bound();
			}
		}
	}

	void release_edge_form() {
		delete edges[0];
		delete edges[1];
		delete[] edges;
		edges = nullptr;

		weights.release();
	}
	
	int find_bipart_i(int idx) {
		for (int i = 0; i < size; i++)
			if (bitstrs[i] == idx)
				return i;
		return -1;
	}


	void check_levels() {
		for (int i = 0; i < levels.get_size(); i++) {
			std::cout << levels[i] << '\n';
		}
	}

	void sort_by_bitstrs(bool flag_include_edges = false) {
		// Bubble sort implemented here.
		size_t n = size;
		bool swapped = false;
		bool isweighted = weights.get_size();

		while (n > 1) {
			size_t new_n = 0;
			for (int i = 1; i < n; i++) {
				if (bitstrs[i - 1] > bitstrs[i]) {
					// Swapped
					using std::swap;
					swap(bitstrs[i - 1], bitstrs[i]);
					if (isweighted)
						swap(weights[i - 1], weights[i]);
					if (flag_include_edges) {
						swap(edges[0][i - 1], edges[0][i]);
						swap(edges[1][i - 1], edges[1][i]);
					}
					new_n = i;
				}
			}
			n = new_n;
		}
	}

	int get_size() const { return size; };

	int* get_bitstr() { return bitstrs; };

	double* get_weight() { return weights.get_vec(); };

	bool isweighted() const { return weights.get_size(); }
	

private:
	int size;
	Array<Array<int> > levels;
	Array<double> weights;
	int** edges;
	int* bitstrs;
};

template <class T>
class TreeSet {
public:
	TreeSet(size_t N) : size(N), trees_newick(Array<TreeArray*>(0, N)),
		Bipart(nullptr), Taxa(nullptr) {};
	TreeSet() : size(0), trees_newick(Array<TreeArray*>(0, 100)),
		Bipart(nullptr), Taxa(nullptr) {};
	TreeSet(Bipartition<T>* bipart, TaxonList* taxa): size(0), trees_newick(Array<TreeArray*>(0, 100)),
		Bipart(bipart), Taxa(taxa) {};
	~TreeSet() {
		for (int i = 0; i < size; i++)
			trees_newick[i]->~TreeArray();
		trees_newick.release();
	}

	void push(TreeArray* tree) { trees_newick.push(tree); size++; };
	TreeArray* operator[](size_t i) { return trees_newick[i]; };

	void ReadTree(std::string fname, size_t pos, int flag_label) {
		Array<int> active_levels(0, 100);
		Array<int> unlabeled(0, 100);

		std::ifstream fin;
		fin.open(fname);
		string temp;

		fin.seekg(pos);

		while (!fin.eof()) {
			temp.clear();
			std::getline(fin, temp, ';');
			if (temp.find('(') != std::string::npos) {
				this->push(new TreeArray(temp, *Taxa, active_levels, unlabeled, flag_label));
			}
		}
	};

	void compute_Bipart() {
		if (Bipart->is_empty()) {
			for (int i = 0; i < Taxa->size; i++) {
				BitString<T> temp = BitString<T>(Taxa->bitstr_size, Taxa->size, i);
				Bipart->push(temp, 0);
			}
		}
		for (int i = 0; i < this->size; i++) {
			this->trees_newick[i]->compute_bitstring(*(this->Bipart), i);
			//this->trees_newick[i]->print_bipart();
		}
	};

	int get_size() { return trees_newick.get_size(); };

	int get_all_bipart() {
		int result = 0;
		for (int i = 0; i < size; i++) 
			result += trees_newick[i]->get_size();
		return result;
	}

	int get_unique_bipart() { return Bipart->get_size(); }

	Bipartition<T>* get_bipart_ptr() { return Bipart; };

	void release_tree() {
		for (int i = 0; i < size; i++)
			trees_newick[i]->release_edge_form();
	}

	void check_bipart(int bipart_idx) {
		std::cout << "Finding bipartition\n";
		Bipart->operator[](bipart_idx).print_BitString(std::cout);
		TreeArray* tree_cur = nullptr;
		for (int i = 0; i < size; i++) {
			tree_cur = trees_newick[i];
			if (tree_cur->find_bipart_i(bipart_idx) != -1)
				std::cout << "Bipartition found in tree " << i << '\n';
		}
	}

private:
	int size;
	Array<TreeArray*> trees_newick;
	Bipartition<T>* Bipart;
	TaxonList* Taxa;
};

template <class T>
class TreeObjects {
	// This class includes essential objects of trees that require further compuatations.

public:
	TreeObjects() : Trees(nullptr), sbipart_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
		n_trees(0), all_bipart(0), unique_bipart(0), isweighted(false) {};

	TreeObjects(TreeSet<T>* trees) : Trees(trees), sbipart_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
		n_trees(trees->get_size()), all_bipart(trees->get_all_bipart()), unique_bipart(trees->get_unique_bipart()),
		isweighted(trees->operator[](0)->isweighted()){};

	void Compute_Bipart_Matrix() {
		for (int i = 0; i < n_trees; i++)
			Trees->operator[](i)->sort_by_bitstrs();

		if (sbipart_mat != nullptr) {
			sbipart_mat->~SparseMatrix();
			sbipart_mat = nullptr;
		}

		double* Vals = new double[all_bipart];
		int* RowInds = new int[all_bipart];
		int* ColPtr = new int[n_trees + 1];

		int index = 0;

		for (int i = 0; i < n_trees; i++) {
			ColPtr[i] = index;
			int n_bipart_in_tree_i = Trees->operator[](i)->get_size();
			if (n_bipart_in_tree_i) {
				memcpy(RowInds + index, Trees->operator[](i)->get_bitstr(), n_bipart_in_tree_i * sizeof(int));
				if (isweighted)
					memcpy(Vals + index, Trees->operator[](i)->get_weight(), n_bipart_in_tree_i * sizeof(double));
				index += n_bipart_in_tree_i;
			}
		}
		ColPtr[n_trees] = index;

		if (!isweighted)
			for (int i = 0; i < unique_bipart * n_trees; i++)
				Vals[i] = 1.0;

		sbipart_mat = new SparseMatrix(unique_bipart, n_trees, Vals, RowInds, ColPtr);

	}

	void Compute_Covariance_Matrix() {};

	void Compute_Distance_Matrix() {};

	bool Compute_RF_Distance(bool ISWEIGHTED) {
		bool bUbid = false; // for counting the number pf unique bipartitions
		unsigned long uBID = 0;
		int index_j, index_k;

		Bipartition<T>* Bipart = Trees->get_bipart_ptr();
		Array<Array<int> >* Id2Tree = Bipart->get_Id2Tree_ptr();

		if (ISWEIGHTED == false)
		{
			dis_mat = new double* [n_trees];
			for (int i = 0; i < n_trees; i++) {
				dis_mat[i] = new double[n_trees];
				memset(dis_mat[i], 0, n_trees * sizeof(double));
			}

			for (int i = 0; i < unique_bipart; i++) {
				int tree_size = Id2Tree->operator()(i).get_size();
				for (int j = 0; j < tree_size; j++)									//There may be better way to do this.
				{
					index_j = Id2Tree->operator()(i)(j);
					for (int k = 0; k < tree_size; k++)
					{
						index_k = Id2Tree->operator()(i)(k);
						if (j == k)
							continue;
						else
							dis_mat[index_j][index_k] += 1;
					}
				}
			}

			//for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti)
			//{
			//	unsigned int sizeVec = vec_hashrf._hashtab2[hti].size();
			//	if (sizeVec)
			//	{
			//		uBID += sizeVec;
			//		if (!bUbid) //This is no longer necessary
			//		{
			//			for (unsigned int i = 0; i < sizeVec; ++i)
			//			{
			//				unsigned int sizeTreeIdx = vec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
			//				if (sizeTreeIdx > 1)
			//				{
			//					for (unsigned int j = 0; j < sizeTreeIdx; ++j)
			//					{
			//						for (unsigned int k = 0; k < sizeTreeIdx; ++k)
			//						{
			//							if (j == k)
			//								continue;
			//							else
			//								dist_URF[vec_hashrf._hashtab2[hti][i]._vec_treeidx[j]][vec_hashrf._hashtab2[hti][i]._vec_treeidx[k]] += 1;
			//						}
			//					}
			//				}
			//			}
			//		}
			//	}
			//}

			for (int i = 0; i < n_trees; i++)
			{
				int n_bipart_i = Trees->operator[](i)->get_size();
				for (int j = 0; j < i; j++)
				{
					int n_bipart_j = Trees->operator[](j)->get_size();
					dis_mat[i][j] = (double)((n_bipart_i + n_bipart_j - 2 * dis_mat[i][j]) / 2);
					dis_mat[j][i] = dis_mat[i][j];
				}
			}
		}
		//else if (ISWEIGHTED == true)
		//{
		//	//vec_hashrf._hashtab.resize(vec_hashrf._hashtab2.size());

		//	//for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti) // For every hash_value_1 appear, pointed by hti
		//	//{
		//	//	unsigned int sizeLinkedList = vec_hashrf._hashtab2[hti].size();
		//	//	if (sizeLinkedList > 0)
		//	//	{
		//	//		for (unsigned int i1 = 0; i1 < sizeLinkedList; ++i1) // For every bipartition has hti as hv1, pointed by i1
		//	//		{
		//	//			unsigned int bidi = vec_hashrf._hashtab2[hti][i1]._vec_treeidx.size();
		//	//			for (unsigned int i2 = 0; i2 < bidi; ++i2) // For every trees that obtain bipartition i1
		//	//			{
		//	//				BUCKET_STRUCT_T bk;
		//	//				bk._hv2 = vec_hashrf._hashtab2[hti][i1]._hv2;
		//	//				bk._t_i = vec_hashrf._hashtab2[hti][i1]._vec_treeidx[i2];
		//	//				bk._dist = vec_hashrf._hashtab2[hti][i1]._vec_dist[i2];
		//	//				vec_hashrf._hashtab[hti].push_back(bk);
		//	//			}// Flatten the storage level as some of the instance are stored are stored in a vector started with i1
		//	//		}
		//	//	}
		//	//}
		//	//vec_hashrf._hashtab2.clear();

		//	// Move all records from _hashtab2 to _hashtab. Records has key (hv1) and beam with values (hv2, tree_id, weight)

		//	for (unsigned int hti = 0; hti < vec_hashrf._hashtab.size(); ++hti) // For hv1 appear, pointed by hti
		//	{
		//		unsigned int sizeLinkedList = vec_hashrf._hashtab[hti].size();
		//		if (sizeLinkedList > 1)
		//		{
		//			vector<unsigned long> vec_hv2;
		//			vector<unsigned long>::iterator itr_vec_hv2;
		//			for (unsigned int i = 0; i < sizeLinkedList; ++i) // For every bipartition with hv1 = hti, pointed by i
		//			{
		//				unsigned long hv2 = vec_hashrf._hashtab[hti][i]._hv2;
		//				if (vec_hv2.empty())
		//					vec_hv2.push_back(hv2);
		//				else
		//				{
		//					itr_vec_hv2 = find(vec_hv2.begin(), vec_hv2.end(), hv2);
		//					if (itr_vec_hv2 == vec_hv2.end())
		//						vec_hv2.push_back(hv2);
		//				}
		//			}// Collect unique hv2 values under hv1 = hti

		//			vector<vector<float> > vec_dist(vec_hv2.size(), vector<float>(n_trees, 0));
		//			//vec_dist is a matrix of size_vec_hv2 * n_trees.

		//			// Set the distance array with distance at proper tree index
		//			for (unsigned int i = 0; i < sizeLinkedList; ++i) // For every bipartition with hv1 = hti, pointed by i
		//			{
		//				for (unsigned int j = 0; j < vec_hv2.size(); ++j) // For every unique bipartition with hv1 = hti, pointed by j
		//				{
		//					if (vec_hashrf._hashtab[hti][i]._hv2 == vec_hv2[j])
		//						vec_dist[j][vec_hashrf._hashtab[hti][i]._t_i] = vec_hashrf._hashtab[hti][i]._dist;
		//				}// store the weight of i-th tree's root at vec_dist[j][T_i], where T_i is the unrooted tree i lies in.
		//			}
		//			// Update floatsim matrix using vec_dist
		//			for (unsigned int i = 0; i < vec_dist.size(); ++i) // For every unique bipartition appear, pointed by i
		//			{
		//				for (unsigned int j = 0; j < vec_dist[i].size(); ++j) // For every trees, pointed by j
		//				{
		//					for (unsigned int k = 0; k < vec_dist[i].size(); ++k) // For every trees, pointed by k
		//					{
		//						if (j == k)
		//							continue;
		//						else
		//							dist_RF[j][k] += (vec_dist[i][j] - vec_dist[i][k] > 0) ? (vec_dist[i][j] - vec_dist[i][k]) : (vec_dist[i][k] - vec_dist[i][j]);
		//					}
		//					// Accumulate absolute difference of the roots weight's dist_RF[i][j] where they are instance(tree) of bipartition i at tree j and k
		//				}
		//			}
		//			vec_hv2.clear();
		//			vec_dist.clear();
		//		}
		//		else if (sizeLinkedList == 1)
		//		{
		//			// propagate the dist value to other tree's distance
		//			for (unsigned int i = 0; i < n_trees; ++i)
		//			{
		//				if (i == vec_hashrf._hashtab[hti][0]._t_i)
		//					continue;
		//				else
		//				{
		//					dist_RF[i][vec_hashrf._hashtab[hti][0]._t_i] += vec_hashrf._hashtab[hti][0]._dist;
		//					dist_RF[vec_hashrf._hashtab[hti][0]._t_i][i] += vec_hashrf._hashtab[hti][0]._dist;
		//				}
		//			}
		//		}
		//	}
		//	for (int i = 0; i < n_trees; i++)
		//	{
		//		for (int j = 0; j < i; j++)
		//		{
		//			dist_RF[i][j] = (dist_RF[i][j] + dist_RF[j][i]) / 4;
		//			dist_RF[j][i] = dist_RF[i][j];
		//		}
		//	}
		//}

		// Temperate sanity check for trees with missing taxa.

		//for (int i = 0; i < n_trees; i++) {
		//	if (Get_num_leaves(treeset[i]->root) != leaveslabelsmaps.size()) {
		//		cout << "Warning! Tree with missing taxa detected! Tree ID: " << i + 1 << ". The associated distances are replaced by -1.\n";
		//		if (!ISWEIGHTED) {
		//			for (int j = 0; j < n_trees; j++) {
		//				dist_URF[i][j] = -1;
		//				dist_URF[j][i] = -1;
		//			}
		//		}
		//		else {
		//			for (int j = 0; j < n_trees; j++) {
		//				dist_RF[i][j] = -1;
		//				dist_RF[j][i] = -1;
		//			}
		//		}
		//	}
		//}

		return true;
	}


private:
	TreeSet<T>* Trees;
	SparseMatrix* sbipart_mat;
	double** cov_mat;
	double** dis_mat;

	int n_trees;
	int all_bipart;
	int	unique_bipart;

	bool isweighted;
};

