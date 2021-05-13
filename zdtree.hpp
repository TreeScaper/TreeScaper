#pragma once

#include <stdlib.h>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>

#include "array.hpp"
#include "sparse.hpp"

using std::string;
using std::map;
using std::ostream;



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
	void dummy_leaf() {
		// dummy_leaf must be inserted as the first leaf.
		assert(size == 0);
		string taxon("dummy_leaf");
		Ind2Taxon.push(taxon);
		Taxon2Ind.insert(std::pair<string, int>(taxon, Ind2Taxon.get_size() - 1));
		size++;
	}
	size_t ReadTaxa(std::string fname);
};


template <class T>
class BitString {
	// This is a bitstring representing the bipartition. It must normalized with 0 at the first bit.
	// T must be unsigned.
private:
	size_t length;
	size_t bit_size;
	T* vec;
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
};


template <class T>
class Bipartition {

private:
	TaxonList* Taxa;
	T hash_bound;
	Array<BitString<T> > Id2BitString;
	//darray_dc<int>* Hash2Id;
	Array2D<int>* Hash2Id;
	T invariant_hashing;
	map<BitString<T>, int> BitStr2Id;

	bool issorted;
public:
	Bipartition() : Taxa(nullptr), hash_bound(0), 
		Id2BitString(Array<BitString<T> >(0, 1000)), 
		Hash2Id(nullptr), /*Id2Tree(0, 500),*/ issorted(false) {};
	Bipartition(TaxonList* taxa, int n_tree, T hb) : Taxa(taxa), hash_bound(hb),
		Id2BitString(Array<BitString<T> >(0, 1000)),
		Hash2Id(new Array2D<int>(hash_bound)), issorted(false) {
		Hash2Id->resize(hash_bound);
		this->set_invariant();
		for (int i = 0; i < Taxa->size; i++)
			this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i).normalized());
	};

	Bipartition(TaxonList* taxa, int n_tree) : Taxa(taxa), Id2BitString(Array<BitString<T> >(0, 1000)),
		Hash2Id(nullptr), issorted(false) {
		this->set_hash_bound(n_tree);
		this->set_invariant();
		
		for (int i = 0; i < Taxa->size; i++)
			this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i).normalized());
	};

	bool is_empty() { return Id2BitString.is_empty(); };
	

	int push(BitString<T>& bs, int index_hash = -1) {

		if (index_hash == -1)
			index_hash = hashing(bs);
		int index = -1;
		if (Hash2Id->get_c_ptr(index_hash) == nullptr) {

			Hash2Id->set_c(index_hash, new Array<int>(0, 50));
			Hash2Id->push(Id2BitString.get_size(), index_hash);
			Id2BitString.push(bs);
			Id2BitString.back().normalized();

			return Id2BitString.get_size() - 1;
		}
		else {
			size_t collison_size = Hash2Id->get_c_ptr(index_hash)->get_size();
			for (int i = 0; i < collison_size; i++)
			{
				index = Hash2Id->ele(index_hash, i);
				if (Id2BitString[index] == bs)
					return index;
			}
			Hash2Id->push(Id2BitString.get_size(), index_hash);
			Id2BitString.push(bs);
			Id2BitString.back().normalized();

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
		Hash2Id = new Array2D<int>(hash_bound);	
		Hash2Id->resize(hash_bound);
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
	//darray_dc<int>* get_Id2Tree_ptr() { return &(Id2Tree); };


	// void print_summary(int n_trees) {
	// 	// int max_collision_size = 0, collision_cnt_5 = 0, collision_cnt_10 = 0, empty_cnt = 0;
	// 	// int cur_size;
	// 	// for (int i = 0; i < hash_bound; i++) {
	// 	// 	if (!Hash2Id->get_c_ptr(i)->is_empty()) {
	// 	// 		cur_size = Hash2Id->get_c_ptr(i)->get_size();
	// 	// 		if (cur_size > max_collision_size)
	// 	// 			max_collision_size = cur_size;
	// 	// 		if (cur_size >= 5)
	// 	// 			collision_cnt_5++;
	// 	// 		if (cur_size >= 10)
	// 	// 			collision_cnt_10++;
	// 	// 		if (cur_size == 0)
	// 	// 			empty_cnt++;
	// 	// 		//std::cout << "Collusion beam # " << i << ":\t";
	// 	// 		//for (int j = 0; j < cur_size; j++)
	// 	// 		//	std::cout << Hash2Id[i]->operator[](j) <<"  ";
	// 	// 		//std::cout << '\n';
	// 	// 	}
	// 	// 	else {
	// 	// 		empty_cnt++;
	// 	// 	}
	// 	// }
	// 	std::cout << "-------------Bipartition info summary------------------\n";
	// 	std::cout << "\tDistinct bipartition number: " << Id2BitString.get_size() << ",\n";
	// 	std::cout << "\tHash container size: " << hash_bound << ",\n";
	// 	// std::cout << "\tEmpty container: " << empty_cnt << ",\n";
	// 	// std::cout << "\tMaximum Collision beam size: " << max_collision_size << ",\n";
	// 	// std::cout << "\tCount of Collision beams over 5: " << collision_cnt_5 << ",\n";
	// 	// std::cout << "\tCount of Collision beams over 10: " << collision_cnt_10 << ".\n";
	// 	std::cout << "-------------Bipartition info summary------------------\n";
	// }

	void print_Bipart(ostream& fout, size_t i) {
		Id2BitString[i].print_BitString(fout);
	}

	BitString<T>& operator[](size_t i) { 
		std::cout << "returning bipartition " << i << "\n";
		return Id2BitString[i]; 
		};
};


class TreeArray {
	// This class use edges forms as internal data structure of trees. It translates Newick form to edge form 
	// by labelling all internal node with integers. The labelling tends to assign small integer to nodes that
	// are close to leaf. The leaf nodes are labelled with (0, ..., l-1) determined by the taxa normalization.
	// Note that the labelling depends on the Newick form, which is not unique, i.e., different Newick form
	// of a same tree gives different edge array. 

private:
	int size;
	//darray<darray<int> > levels;
	Array<Array<int> >* levels;
	Array<double>* weights; // weights will be kept and used in generalizing B2T matrix.
	int** edges;
	Array<size_t>* t2b;
	bool issorted;

public:
	TreeArray() : size(0), levels(nullptr), weights(nullptr),
		t2b(nullptr), edges(nullptr), issorted(false) {};
	~TreeArray() {
		this->release_edge_form();
		delete levels;
		delete weights;
		delete t2b;
	};
	TreeArray(string& tree, TaxonList& taxon_list, Array<int>& active_rows, Array<int>& level_flag, int flag_label, bool ISROOTED, bool ISWEIGHTED);
	//TreeArray(string& tree, TaxonList& taxon_list, Array<int>& active_rows, Array<int>& level_flag, Array2D<double>& w_temp, int flag_label, bool ISROOTED);

	//void print_level(std::ostream& output) { levels.print(output); };

	void label_internal_node(Array<int>& active_rows, Array<int>& unlabeled);
	void label_internal_node(Array<int>& active_rows, Array<int>& unlabeled, Array<Array<double> >& w_temp);

	template <class T>
	void compute_bitstring(Bipartition<T>* Bipart) {

		//std::cout << "------------------Bitstring computation started----------------------\n";
		t2b = new Array<size_t>(size);
		t2b->resize(size, (size_t)-1);
		T* index_hash = new T[size];
		memset(index_hash, 0, size * sizeof(T));

		//darray<int>* t2b_ptr = TreeId2Bipart->get_c_ptr(tree_id);

		//if (!t2b_ptr->is_empty())
		//	t2b_ptr->resize(size, -1);

		//t2b_ptr->resize(size, -1);

		BitString<T>* bs = new BitString<T>[size];
		T iv_h = Bipart->get_invariant_hashing();
		T h_b = Bipart->get_hash_bound();
		int bitstr_size = Bipart->get_bitstr_size();
		int l_s = Bipart->get_leaf_size();


		for (int i = 0; i < size; i++) {
			if (i < l_s) {
				bs[i] = BitString<T>(bitstr_size, l_s, i);
				index_hash[i] = Bipart->hashing(bs[i]);
			}
			else
				bs[i] = BitString<T>(bitstr_size, l_s);
		}

		index_hash[0] = 1;
		// Note that bs[0] = 000000000...000001 = 1 is not normalized.
		// Therefore hashing() will return the hash of its complement's hash.
		// However, we need 1 mod M = 1 here.


		for (int i = 0; i < size; i++) {
			if (edges[1][i] < l_s) 
				t2b->ele(i) = edges[1][i];
			else {
				T hv = index_hash[edges[1][i]];
				if (bs[edges[1][i]].is_normalized()) 
					t2b->ele(i) = Bipart->push(bs[edges[1][i]], hv);
				// use the computed hashing value for normalized bitstring generated from bitwise OR.
				else 
					t2b->ele(i) = Bipart->push(bs[edges[1][i]], (hv > iv_h ? h_b + iv_h - hv : iv_h - hv));
				// use the hashing of complement bitstring.
			}

			if (edges[0][i] < size) {
				bs[edges[0][i]] |= bs[edges[1][i]];
				index_hash[edges[0][i]] += index_hash[edges[1][i]];
				index_hash[edges[0][i]] = index_hash[edges[0][i]] % h_b;
			}
		}

		// weights.release();
		// Weights record have been transferred to a table referenced by bipartition, which helps inquiry for computing distance.
		// By keeping this weights, we can keep fast inquiry from tree, but then it keeps two copies of the weights info.

	}

	void release_edge_form() {
		delete edges[0];
		delete edges[1];
		delete[] edges;
		edges = nullptr;

		//weights.release();
		//weights is needed for B2T matrix
	}

	void release_level() {
		delete levels;
		levels = nullptr;
	}

	int find_bipart_i(int idx) {
		auto it_size = t2b->get_size();
		for (int i = 0; i < it_size; i++)
			if (t2b->ele(i) == idx)
				return i;
		return -1;
	}

	/*
	//void sort_by_bitstrs(bool flag_include_edges = false) {
	//	// Bubble sort implemented here.
	//	if (issorted)
	//		return;

	//	size_t n = size;
	//	bool swapped = false;
	//	bool isweighted = !weights->is_empty();

	//	while (n > 1) {
	//		size_t new_n = 0;
	//		int i = 1;
	//		for (auto it = t2b->begin(); i < n; it++, i++) {
	//			if (*it > *(it + 1)) {
	//				// Swapped
	//				using std::swap;
	//				swap(*it, *(it + 1));
	//				if (isweighted)
	//					swap(weights->ele(i - 1), weights->ele(i));
	//				if (flag_include_edges) {
	//					swap(edges[0][i - 1], edges[0][i]);
	//					swap(edges[1][i - 1], edges[1][i]);
	//				}
	//				new_n = i;
	//			}
	//		}
	//		n = new_n;
	//	}
	//	issorted = true;
	//}
	*/

	int get_size() const { return size; };

	Array<size_t>* get_t2b() { return t2b; };

	Array<size_t>* pop_t2b() { Array<size_t>* b_ptr = t2b; t2b = nullptr; return b_ptr; };

	Array<double>* get_weight() { assert(weights != nullptr);  return weights; };

	Array<double>* pop_weight() { Array<double>* w_ptr = weights; weights = nullptr; return w_ptr; };
};

template <class T>
class TreeSet {
private:
	int size;
	Array<TreeArray*> trees;
	Bipartition<T>* Bipart;
	TaxonList* Taxa;
public:
	TreeSet(size_t N) : size(0), trees(Array<TreeArray*>(0, N)),
		Bipart(nullptr), Taxa(nullptr) {};
	TreeSet() : size(0), trees(Array<TreeArray*>(0, 100)),
		Bipart(nullptr), Taxa(nullptr) {};
	TreeSet(Bipartition<T>* bipart, TaxonList* taxa): size(0), trees(Array<TreeArray*>(0, 100)),
		Bipart(bipart), Taxa(taxa) {};
	~TreeSet() {
		for (int i = 0; i < size; i++)
			delete trees[i];
	}

	void push(TreeArray* tree) { trees.push(tree); size++; };

	TreeArray* ele_ptr(size_t i) { return trees[i]; };

	TreeArray& ele(size_t i) { return *(trees[i]); };

	TreeArray& operator[](size_t i) { return this->ele(i); };

	void ReadTree(std::string fname, size_t pos, int flag_label, bool flag_rooted, bool flag_weighted) {
		size_t leaf_size = Taxa->size;
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
					this->push(new TreeArray(temp, *Taxa, active_levels, unlabeled, flag_label, flag_rooted, flag_weighted));
			}
		}

		fin.close();
	};

	int get_size() { return trees.get_size(); };

	int get_all_bipart_size() {
		int result = 0;
		for (int i = 0; i < size; i++) 
			result += trees[i]->get_size();
		return result;
	}

	int get_unique_bipart_size() { return Bipart->get_size(); }

	Bipartition<T>* get_bipart_ptr() { return Bipart; };

	TaxonList* get_taxa_ptr() { return Taxa; };

	void release_tree_edge_form() {
		for (int i = 0; i < size; i++)
			trees[i]->release_edge_form();
	}

	void release_tree(){
		for (int i = 0; i < size; i++)
			delete trees[i];
		size = 0;
	}

};
