#pragma once

#include <stdlib.h>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>

#include "array.hpp"
#include "sparse.hpp"

using std::map;
using std::ostream;
using std::string;

template <class T>
T GetPrime(T topNum, unsigned from)
{
	T primeNum = 0;
	T candidate = 0;

	if (topNum <= 100)
		candidate = 2;
	else
		candidate = topNum;

	while (candidate <= topNum + from)
	{
		T trialDivisor = 2;
		int prime = 1;

		while (trialDivisor * trialDivisor <= candidate)
		{
			if (candidate % trialDivisor == 0)
			{
				prime = 0;
				break;
			}
			trialDivisor++;
		}
		if (prime)
			primeNum = candidate;

		candidate++;
	}

	return primeNum;
};

class TaxonList
{
public:
	TaxonList() : size(0), bitstr_size(0){};
	int size;
	int bitstr_size;
	Array<string> Ind2Taxon;
	map<string, int> Taxon2Ind;
	Array<int> IndB2IndA;

	void push(string item, bool flag_str_format = true);
	void insert(string item);
	bool cmp_taxa(string &tree, bool *barray);
	bool report_missing_taxon(size_t tree_id, string &tree, bool *barray, std::ostream &out);

	template <class T>
	void set_bitstr_size(T ele)
	{
		bitstr_size = (size + (-size) % (8 * sizeof(T))) / (8 * sizeof(T));
	};
	int get_bitstr_size() const { return bitstr_size; };
	void dummy_leaf()
	{
		// dummy_leaf must be inserted as the first leaf.
		assert(size == 0);
		string taxon("dummy_leaf");
		Ind2Taxon.push(taxon);
		Taxon2Ind.insert(std::pair<string, int>(taxon, Ind2Taxon.get_size() - 1));
		size++;
	}
	size_t ReadTaxa(std::string fname);
	size_t ReadTaxa(std::string fname, bool same_leaf);
	bool ScanTaxa(std::string fname, size_t tree_pos, std::string outname);
	string operator()(int i) const { return Ind2Taxon(i); };
	string get_Taxon(int i) const { return Ind2Taxon(i); };
};

template <class T>
class BitString
{
	// This is a bitstring representing the bipartition. It must normalized with 0 at the first bit.
	// T must be unsigned.
private:
	size_t length;
	size_t bit_size; // These two terms can be better managed.
	T *vec;

public:
	friend void swap(BitString<T> &lhs, BitString<T> &rhs)
	{
		using std::swap;
		swap(lhs.length, rhs.length);
		swap(lhs.bit_size, rhs.bit_size);
		swap(lhs.vec, rhs.vec);
	}

	BitString() : length(0), vec(nullptr){};
	BitString(int n, int b_size) : length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr) { memset(vec, 0, length * sizeof(T)); };
	BitString(int n, int b_size, int index) : length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr)
	{
		if (n > 0)
		{
			if (index >= 0)
			{
				memset(vec, 0, length * sizeof(T));
				int pos_byte = index / (8 * sizeof(T));
				int pos_bit = index % (8 * sizeof(T));
				T bit = 1;
				vec[pos_byte] |= bit << pos_bit;
			}
			else // This can be done more efficiently
			{
				throw(1);
				// memset(vec, 0, length * sizeof(T));
				// T bit = 1;
				// for (int i = 0; i < b_size; i++)
				// {
				// 	int pos_byte = i / (8 * sizeof(T));
				// 	int pos_bit = i % (8 * sizeof(T));
				// 	vec[pos_byte] |= bit << pos_bit;
				// }
			}
		}
	}
	BitString(int n, int b_size, T *src) : length(n), bit_size(b_size), vec(n > 0 ? new T[n] : nullptr)
	{
		memcpy(vec, src, length * sizeof(T));
	}
	BitString(const BitString<T> &src) : length(src.length), bit_size(src.bit_size), vec(length > 0 ? new T[length] : nullptr)
	{
		memcpy(vec, src.vec, length * sizeof(T));
	}

	~BitString() { delete[] vec; };

	BitString<T> operator|(const BitString<T> &rhs)
	{
		T *temp = new T[length];
		for (int i = 0; i < length; i++)
			temp[i] = vec[i] | rhs.vec[i];
		return BitString<T>(length, bit_size, temp);
	}

	BitString<T> &operator|=(const BitString<T> &rhs)
	{
		for (int i = 0; i < length; i++)
			this->vec[i] |= rhs.vec[i];
		return *this;
	}

	BitString<T> operator%(const BitString<T> &rhs)
	{
		T *temp = new T[length];
		for (int i = 0; i < length; i++)
			temp[i] = vec[i] - rhs.vec[i];
		return BitString<T>(length, bit_size, temp);
	}

	friend BitString<T> operator%(const BitString<T> &lhs, const BitString<T> &rhs)
	{
		T *temp = new T[lhs.length];
		for (int i = 0; i < lhs.length; i++)
			temp[i] = lhs.vec[i] - rhs.vec[i];
		return BitString<T>(lhs.length, lhs.bit_size, temp);
	}

	BitString<T> &operator=(const BitString<T> &rhs)
	{
		BitString<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	bool operator==(const BitString<T> &rhs)
	{
		if (length != rhs.length || bit_size != rhs.bit_size)
		{
			return false;
		}
		T temp = vec[0] ^ rhs.vec[0];

		if (temp & (T)1)
		{
			for (int i = 0; i < length - 1; i++)
			{
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
		else
		{
			for (int i = 0; i < length; i++)
			{
				if (vec[i] != rhs.vec[i])
					return false;
			}
			return true;
		}
	}

	bool operator<(const BitString<T> &rhs)
	{ // This only works for normalized bitstring
		T mask;

		for (int i = 0; i < length; i++)
		{
			mask = vec[i] ^ rhs.vec[i]; // mask for all different bits
			if (mask == 0)				// no different bit found
				continue;
			else
			{
				mask = (mask & -mask); // mask for the rightmost bit
				if (vec[i] & mask)
					return false;
			}
		}
		if (mask == 0)
			return false;
		else
			return true;
	}

	bool operator>(const BitString<T> &rhs)
	{
		for (int i = 0; i < length; i++)
		{
			if (vec[i] > rhs.vec[i])
				return true;
			else if (vec[i] < rhs.vec[i])
				return false;
		}
		return false;
	}
	
	friend int bitstr_cmp(const BitString<T> &lhs, const BitString<T> &rhs)
	{
		for (int i = 0; i < lhs.length; i++)
			if (lhs.vec[i] > rhs.vec[i])
				return 1;
			else if(lhs.vec[i] < rhs.vec[i])
				return -1;
		return 0;
	}

	friend bool operator<(const BitString<T> &lhs, const BitString<T> &rhs)
	{
		T mask = 0;

		for (int i = 0; i < lhs.length; i++)
		{
			mask = lhs.vec[i] ^ rhs.vec[i]; // mask for all different bits
			if (mask == 0)					// no different bit found
				continue;
			else
			{
				mask = (mask & (~mask + 1)); // mask for the rightmost bit
				if (lhs.vec[i] & mask)
					return false;
			}
		}

		if (mask == 0)
			return false;
		else
			return true;
	}

	friend bool operator>(const BitString<T> &lhs, const BitString<T> &rhs) {return bitstr_cmp(lhs, rhs) == 1;};

	bool operator!=(const BitString<T> &rhs)
	{
		if (length != rhs.length)
		{
			return true;
		}
		else
		{
			for (int i = 0; i < length; i++)
			{
				if (vec[i] != rhs.vec[i])
					return true;
			}
			return false;
		}
	}

	BitString<T> &set_bit(int index)
	{
		int pos_byte = index / (8 * sizeof(T));
		int pos_bit = index % (8 * sizeof(T));
		T bit = 1;
		this->vec[pos_byte] |= (bit << pos_bit);
		return *this;
	};

	T &operator[](size_t i) { return vec[i]; };

	T operator()(size_t i) const { return vec[i]; };

	int size() const { return length; };

	BitString<T> &normalized()
	{
		if (!this->is_normalized())
		{ // First bit is 1
			for (int i = 0; i < length; i++)
				vec[i] = ~vec[i];
			size_t remainder = 8 * length * sizeof(T) - bit_size;
			T temp = -1;
			temp = (temp >> remainder);
			vec[length - 1] = vec[length - 1] & temp;
		}
		return *this;
	}

	bool is_normalized() const
	{
		// assert(this->length != 0);
		return !((bool)(vec[0] & (T)1));
	}

	friend bool is_subset(const BitString<T> &lhs, const BitString<T> &rhs)
	{
		// This function assumes compatible bitstring, i.e., there could be 1 out of the following 3 scenarios
		// 		1: lhs \subsetneq rhs, 				it returns true;
		// 		2: lhs \subsetneq rhs,				it returns false;
		//		3: lhs \cup rhs = \emptyset, 		it returns false;

		// Therefore, it performs check on lhs & rhs == lhs only.
		// For more stable check on arbitrary bitstrings, use check_subset instead, which is slower but includes sanity check.

		for (int i = 0; i < lhs.length; i++)
			if ((lhs.vec[i] & rhs.vec[i]) != lhs.vec[i])
				return false;
		return true;
	}

	int popcount() const
	{
		int population = 0;
		for (int i = 0; i < this->length; i++)
			population += __builtin_popcount(this->vec[i]);
		return population;
	}

	void print_BitString(ostream &fout)
	{
		bool bit;
		size_t n = sizeof(T) * 8;
		T one = 1, temp = 0;
		for (int i = 0; i < length - 1; i++)
		{
			temp = vec[i];
			for (int j = 0; j < n; j++)
			{
				bit = temp & one;
				temp = temp >> 1;
				fout << bit;
			}
		}
		temp = vec[length - 1];
		n = bit_size - (length - 1) * sizeof(T) * 8;
		for (int j = 0; j < n; j++)
		{
			bit = temp & one;
			temp = temp >> 1;
			fout << bit;
		}
	}

	T *get_vec() { return vec; };
	size_t get_length() { return length; };
	size_t get_bitsize() { return bit_size; };

	int get_taxon_id() const
	{
		int pos_bit = 0;
		int pos_byte = 0;

		for (int i = 0; i < length; i++)
		{
			if (vec[i] != 0)
			{
				pos_byte = i;
				while (vec[i] >> (pos_bit + 1) != 0)
					pos_bit += 1;
				break;
			}
		}

		return 8 * sizeof(T) * pos_byte + pos_bit;
	}

	friend int get_taxon_id(const BitString<T> &bitstr, const BitString<T> &leafset)
	{
		int pop = bitstr.popcount();

		if (pop == 1)
			return bitstr.get_taxon_id();
		else
		{
			auto comp_bitstr = leafset % bitstr;
			pop = comp_bitstr.popcount();
			if (pop == 1)
				return comp_bitstr.get_taxon_id();
			else
				return -1;
		}
	}
};

template <class T>
class Bipartition
{

private:
	TaxonList *Taxa;
	T hash_bound;
	Array<BitString<T>> Id2BitString;
	Array<BitString<T>> Id2FullLeaf;
	// darray_dc<int>* Hash2Id;
	Array2D<int> *Hash2Id;
	T invariant_hashing;
	// map<BitString<T>, int> BitStr2Id;

	bool issorted;

public:
	Bipartition() : Taxa(nullptr), hash_bound(0),
					Id2BitString(Array<BitString<T>>(0, 1000)),
					Id2FullLeaf(Array<BitString<T>>(0, 1000)),
					Hash2Id(nullptr), /*Id2Tree(0, 500),*/ issorted(false){};
	Bipartition(TaxonList *taxa, int n_tree, T hb) : Taxa(taxa), hash_bound(hb),
													 Id2BitString(Array<BitString<T>>(0, 1000)),
													 Id2FullLeaf(Array<BitString<T>>(0, 1000)),
													 Hash2Id(new Array2D<int>(hash_bound)), issorted(false)
	{
		Hash2Id->resize(hash_bound);
		this->set_invariant();
		for (int i = 0; i < Taxa->size; i++)
			this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i).normalized());
	};

	Bipartition(TaxonList *taxa, int n_tree, bool sameleaf) : Taxa(taxa), Id2BitString(Array<BitString<T>>(0, 1000)),
											   				  Hash2Id(nullptr), issorted(false)
	{
		this->set_hash_bound(n_tree);
		this->set_invariant();

		for (int i = 0; i < Taxa->size; i++)
			if (sameleaf)
				this->push(BitString<T>(Taxa->bitstr_size, Taxa->size, i).normalized());
			else
				this->push_full(BitString<T>(Taxa->bitstr_size, Taxa->size, i));
	};

	bool is_empty() { return Id2BitString.is_empty(); };

	int push(BitString<T> &bs, int index_hash = -1)
	{

		if (index_hash == -1)
			index_hash = hashing(bs);
		int index = -1;
		if (Hash2Id->get_c_ptr(index_hash) == nullptr)
		{

			Hash2Id->set_c(index_hash, new Array<int>(0, 50));
			Hash2Id->push(Id2BitString.get_size(), index_hash);
			Id2BitString.push(bs);
			Id2BitString.back().normalized();

			return Id2BitString.get_size() - 1;
		}
		else
		{
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

	int push_full(const BitString<T> &bs, int index_hash = -1)
	{

		if (index_hash == -1)
			index_hash = hashing_mod(bs);
		int index = -1;
		if (Hash2Id->get_c_ptr(index_hash) == nullptr)
		{

			Hash2Id->set_c(index_hash, new Array<int>(0, 50));
			Hash2Id->push(Id2BitString.get_size(), index_hash);
			Id2BitString.push(bs);

			return Id2BitString.get_size() - 1;
		}
		else
		{
			size_t collison_size = Hash2Id->get_c_ptr(index_hash)->get_size();
			for (int i = 0; i < collison_size; i++)
			{
				index = Hash2Id->ele(index_hash, i);
				if (Id2BitString[index] == bs)
					return index;
			}
			Hash2Id->push(Id2BitString.get_size(), index_hash);
			Id2BitString.push(bs);

			return Id2BitString.get_size() - 1;
		}
	};

	// Add BitString to the list and return the id assigned to it. If it has encountered before, do nothing.

	void push_bs_only(BitString<T> &bs) { Id2BitString.push(bs); };

	T hashing(const BitString<T> &bs) const
	{
		T result = 0;
		for (int i = 0; i < bs.size(); i++)
		{
			result += bs(i) % hash_bound;
			result %= hash_bound;
		}
		if (bs.is_normalized()) // Is there a better way to access the top bit?
			return result;
		else
			return (result > invariant_hashing ? hash_bound + invariant_hashing - result : invariant_hashing - result);
	}

	T hashing_mod(const BitString<T> &bs) const
	{
		T result = 0;
		for (int i = 0; i < bs.size(); i++)
		{
			result += bs(i) % hash_bound;
			result %= hash_bound;
		}
		return result;
	}

	void set_hash_bound(int n_tree)
	{
		T p = 0;
		unsigned int mul = 1;
		T top = (n_tree * (2 * Taxa->size - 3)) / 10;
		do
		{
			unsigned from = 100 * mul;
			p = GetPrime(top, from);
			++mul;
		} while (p == 0);
		hash_bound = (p < 101 ? 101 : p);
		Hash2Id = new Array2D<int>(hash_bound);
		Hash2Id->resize(hash_bound);
	}

	void set_hash_table(Array2D<int> *h2i)
	{
		if (Hash2Id != nullptr)
			delete Hash2Id;
		Hash2Id = h2i;
	}

	void set_issorted(bool iss)
	{
		issorted = iss;
	};

	void set_invariant()
	{
		// compute the hash value of the sum of normalized bitstring and its complement.
		//  example: taxa size: 6, bit size of T: 4
		//  bitstring:	0010 0001
		//  complement:	1101 0010
		//  sum:			1111 0011
		T result = 0, comp = -1;
		int bitstr_size = Taxa->bitstr_size;
		int remainder = bitstr_size * (8 * sizeof(T)) - Taxa->size;

		for (int i = 0; i < bitstr_size - 1; i++)
		{
			result += comp % hash_bound;
			result %= hash_bound;
		}
		comp = comp >> remainder; // fill zeros in the last container.
		result += comp % hash_bound;
		result %= hash_bound;
		invariant_hashing = result;
	}

	void set_taxon_list(TaxonList *taxa)
	{

		if (Taxa != nullptr)
			delete Taxa;
		Taxa = taxa;
	}

	int get_leaf_size() const { return Taxa->size; };
	int get_bitstr_size() const { return Taxa->bitstr_size; };
	int get_size() const { return Id2BitString.get_size(); };
	T get_hash_bound() const { return hash_bound; };
	T get_invariant_hashing() const { return invariant_hashing; };
	Array2D<int> *get_hash_table() const { return Hash2Id; };
	bool get_issorted() const { return issorted; };
	const BitString<T> &get_BitString(int i) const { return Id2BitString(i); };
	// darray_dc<int>* get_Id2Tree_ptr() { return &(Id2Tree); };

	// void print_Bipart(ostream& fout, int n_tree) {
	//	size_t n = Id2BitString.get_size();
	//	fout << Taxa->size << '\t' << n << '\n';
	//	for (int i = 0; i < n; i++){
	//		fout << i << " ";
	//		Id2BitString[i].print_BitString(fout);
	//		if (Id2Tree[i][0] == -1)
	//			fout << ' ' <<  n_tree << '\n';
	//		else
	//			fout << ' ' << Id2Tree[i].get_size() << '\n';
	//	}
	// }

	void print_Bipart(ostream &fout, size_t i)
	{
		Id2BitString[i].print_BitString(fout);
	}

	BitString<T> &operator[](size_t i) { return Id2BitString[i]; };

	// BitString<T> &operator()(size_t i) const { return Id2BitString[i]; };
};

class TreeArray
{
	// This class use edges forms as internal data structure of trees. It translates Newick form to edge form
	// by labelling all internal node with integers. The labelling tends to assign small integer to nodes that
	// are close to leaf. The leaf nodes are labelled with (0, ..., l-1) determined by the taxa normalization.
	// Note that the labelling depends on the Newick form, which is not unique, i.e., different Newick form
	// of a same tree gives different edge array.

private:
	int size;
	int leafsize;
	int maxleaflabel;
	// darray<darray<int> > levels;
	Array<Array<int>> *levels;
	Array<PRECISION> *weights; // weights will be kept and used in generalizing B2T matrix.
	int **edges;
	Array<size_t> *t2b;
	Array<size_t> *t2b_upper;

	bool issorted;
	bool sameleaf;

public:
	TreeArray() : size(0), leafsize(0), maxleaflabel(-1),
				  levels(nullptr), weights(nullptr), edges(nullptr),
				  t2b(nullptr), t2b_upper(nullptr),
				  issorted(false), sameleaf(true){};
	~TreeArray()
	{
		this->release_edge_form();
		if (levels)
			delete levels;
		if (weights)
			delete weights;
		if (t2b)
			delete t2b;
		if (t2b_upper)
			delete t2b_upper;
	};
	TreeArray(string &tree, TaxonList &taxon_list, Array<int> &active_rows, Array<int> &level_flag, int flag_label, bool ISROOTED, bool ISWEIGHTED);
	// TreeArray(string& tree, TaxonList& taxon_list, Array<int>& active_rows, Array<int>& level_flag, Array2D<PRECISION>& w_temp, int flag_label, bool ISROOTED);

	// void print_level(std::ostream& output) { levels.print(output); };

	void label_internal_node(Array<int> &active_rows, Array<int> &unlabeled, int leaf_size);
	void label_internal_node(Array<int> &active_rows, Array<int> &unlabeled, Array<Array<PRECISION>> &w_temp, int leaf_size);
	void set_same_leaf(bool sl) {sameleaf = sl;};
	// template <class T>
	// void compute_bitstring(Bipartition<T> *Bipart)
	// {

	// 	// std::cout << "------------------Bitstring computation started----------------------\n";
	// 	t2b = new Array<size_t>(size);
	// 	t2b->resize(size, (size_t)-1);
	// 	T *index_hash = new T[size];
	// 	memset(index_hash, 0, size * sizeof(T));

	// 	// darray<int>* t2b_ptr = TreeId2Bipart->get_c_ptr(tree_id);

	// 	// if (!t2b_ptr->is_empty())
	// 	//	t2b_ptr->resize(size, -1);

	// 	// t2b_ptr->resize(size, -1);

	// 	BitString<T> *bs = new BitString<T>[size];
	// 	T iv_h = Bipart->get_invariant_hashing();
	// 	T h_b = Bipart->get_hash_bound();
	// 	int bitstr_size = Bipart->get_bitstr_size();
	// 	int l_s = Bipart->get_leaf_size();

	// 	for (int i = 0; i < size; i++)
	// 	{
	// 		if (i < l_s)
	// 		{
	// 			bs[i] = BitString<T>(bitstr_size, l_s, i);
	// 			index_hash[i] = Bipart->hashing(bs[i]);
	// 		}
	// 		else
	// 			bs[i] = BitString<T>(bitstr_size, l_s);
	// 	}

	// 	index_hash[0] = 1;
	// 	// Note that bs[0] = 000000000...000001 = 1 is not normalized.
	// 	// Therefore hashing() will return the hash of its complement's hash.
	// 	// However, we need 1 mod M = 1 here.

	// 	for (int i = 0; i < size; i++)
	// 	{
	// 		if (edges[1][i] < l_s)
	// 			t2b->ele(i) = edges[1][i];
	// 		else
	// 		{
	// 			T hv = index_hash[edges[1][i]];
	// 			if (bs[edges[1][i]].is_normalized())
	// 				t2b->ele(i) = Bipart->push(bs[edges[1][i]], hv);
	// 			// use the computed hashing value for normalized bitstring generated from bitwise OR.
	// 			else
	// 				t2b->ele(i) = Bipart->push(bs[edges[1][i]], (hv > iv_h ? h_b + iv_h - hv : iv_h - hv));
	// 			// use the hashing of complement bitstring.
	// 		}

	// 		if (edges[0][i] < size)
	// 		{
	// 			bs[edges[0][i]] |= bs[edges[1][i]];
	// 			index_hash[edges[0][i]] += index_hash[edges[1][i]];
	// 			index_hash[edges[0][i]] = index_hash[edges[0][i]] % h_b;
	// 		}
	// 	}

	// 	// weights.release();
	// 	// Weights record have been transferred to a table referenced by bipartition, which helps inquiry for computing distance.
	// 	// By keeping this weights, we can keep fast inquiry from tree, but then it keeps two copies of the weights info.
	// }

	template <class T>
	BitString<T> *compute_bitstring(Bipartition<T> *Bipart)
	{
		if (sameleaf)
			return compute_bitstring_sameleaf(Bipart);
		else
			return compute_bitstring_diffleaf(Bipart);
	}

	template <class T>
	BitString<T> *compute_bitstring_sameleaf(Bipartition<T> *Bipart)
	{
		// std::cout << "------------------Bitstring computation started----------------------\n";
		t2b = new Array<size_t>(size);
		t2b->resize(size, (size_t)-1);


		// darray<int>* t2b_ptr = TreeId2Bipart->get_c_ptr(tree_id);

		// if (!t2b_ptr->is_empty())
		//	t2b_ptr->resize(size, -1);

		// t2b_ptr->resize(size, -1);

		T iv_h = Bipart->get_invariant_hashing();
		T h_b = Bipart->get_hash_bound();
		int bitstr_size = Bipart->get_bitstr_size();
		int all_leaf_size = Bipart->get_leaf_size();

		// BitString<T> *bs = new BitString<T>[size];
		// T *index_hash = new T[size];
		// memset(index_hash, 0, size * sizeof(T));
		// When the node labels are not sequential, it is not ideal to use simple array for reverse reference is required.
		// They are sequential only when all trees share the same leafset that is identical to the provided taxon list.

		std::map<int, BitString<T>> bs;
		std::map<int, T> index_hash;


		for (int i = 0; i < size; i++)
		{
			if (edges[1][i] <= maxleaflabel)
			{
				bs[edges[1][i]] = BitString<T>(bitstr_size, all_leaf_size, edges[1][i]);
				index_hash[edges[1][i]] = Bipart->hashing(bs[edges[1][i]]);

				// bs[edges[1][i]].print_BitString(std::cout);
				// std::cout << "\n";
			}
			else
				bs[edges[1][i]] = BitString<T>(bitstr_size, all_leaf_size);
		}

		for (int i = 0; i < size; i++)
		{
			if (edges[1][i] <= maxleaflabel)
				t2b->ele(i) = edges[1][i];
			else
			{
				T hv = index_hash[edges[1][i]];
				if (bs[edges[1][i]].is_normalized())
					t2b->ele(i) = Bipart->push(bs[edges[1][i]], hv);
				// use the computed hashing value for normalized bitstring generated from bitwise OR.
				else
					t2b->ele(i) = Bipart->push(bs[edges[1][i]], (hv > iv_h ? h_b + iv_h - hv : iv_h - hv));
				// use the hashing of complement bitstring.
			}

			if (edges[0][i] < size)
			{
				bs[edges[0][i]] |= bs[edges[1][i]];
				index_hash[edges[0][i]] += index_hash[edges[1][i]];
				index_hash[edges[0][i]] = index_hash[edges[0][i]] % h_b;

				// std::cout << "\n";
				// bs[edges[1][i]].print_BitString(std::cout);
				// std::cout << "\n>>\n";
				// bs[edges[0][i]].print_BitString(std::cout);
				// std::cout << "\n";
			}
		}

		// delete[] index_hash;
		// delete[] bs;

		// weights.release();
		// Weights record have been transferred to a table referenced by bipartition, which helps inquiry for computing distance.
		// By keeping this weights, we can keep fast inquiry from tree, but then it keeps two copies of the weights info.

		return nullptr;
	}

	template <class T>
	BitString<T> *compute_bitstring_diffleaf(Bipartition<T> *Bipart)
	{

		// std::cout << "------------------Bitstring computation started----------------------\n";
		t2b = new Array<size_t>(size);
		t2b_upper = new Array<size_t>(size);

		t2b->resize(size, (size_t)-1);
		t2b_upper->resize(size, (size_t)-1);

		// darray<int>* t2b_ptr = TreeId2Bipart->get_c_ptr(tree_id);

		// if (!t2b_ptr->is_empty())
		//	t2b_ptr->resize(size, -1);

		// t2b_ptr->resize(size, -1);

		int bitstr_size = Bipart->get_bitstr_size();
		int all_leaf_size = Bipart->get_leaf_size();

		// BitString<T> *bs = new BitString<T>[all_leaf_size + size - leafsize];
		std::map<int, BitString<T>> bs;

		auto leafset = new BitString<T>(bitstr_size, all_leaf_size);

		for (int i = 0; i < size; i++)
		{
			if (edges[1][i] <= maxleaflabel)
			{
				bs[edges[1][i]] = BitString<T>(bitstr_size, all_leaf_size, edges[1][i]);
				leafset->set_bit(edges[1][i]);
			}
			else
				bs[edges[1][i]] = BitString<T>(bitstr_size, all_leaf_size);
		}
		// leafset->print_BitString(std::cout);
		// std::cout<<"\n";

		for (int i = 0; i < size; i++)
		{
			if (edges[1][i] <= maxleaflabel)
			{
				auto cbs = (*leafset) % bs[edges[1][i]];

				// std::cout << "Original bitstring\n";
				// bs[edges[1][i]].print_BitString(std::cout);
				// std::cout << "\nComplementary bitstring\n";
				// cbs.print_BitString(std::cout);
				// std::cout << '\n';

				if (bitstr_cmp(cbs, bs[edges[1][i]]) == 1)
				{
					t2b->ele(i) = edges[1][i];
					t2b_upper->ele(i) = Bipart->push_full(cbs);
				}
				else
				{
					t2b_upper->ele(i) = edges[1][i];
					t2b->ele(i) = Bipart->push_full(cbs);
				}
			}
			else
			{
				auto cbs = (*leafset) % bs[edges[1][i]];

				// std::cout << "Original bitstring\n";
				// bs[edges[1][i]].print_BitString(std::cout);
				// std::cout << "\nComplementary bitstring\n";
				// cbs.print_BitString(std::cout);
				// std::cout << '\n';

				if (bitstr_cmp(cbs, bs[edges[1][i]]) == 1)
				{
					t2b->ele(i) = Bipart->push_full(bs[edges[1][i]]);
					t2b_upper->ele(i) = Bipart->push_full(cbs);
				}
				else
				{
					t2b_upper->ele(i) = Bipart->push_full(bs[edges[1][i]]);
					t2b->ele(i) = Bipart->push_full(cbs);
				}

				
				// std::cout << "\n";
				// (*Bipart)[t2b->ele(i)].print_BitString(std::cout);
				// std::cout << "\n";
				// (*Bipart)[t2b_upper->ele(i)].print_BitString(std::cout);
				// std::cout << "\n";
			}

			if (edges[0][i] < size)
				bs[edges[0][i]] |= bs[edges[1][i]];
		}

		// weights.release();
		// Weights record have been transferred to a table referenced by bipartition, which helps inquiry for computing distance.
		// By keeping this weights, we can keep fast inquiry from tree, but then it keeps two copies of the weights info.

		return leafset;
	}

	template <class T>
	BitString<T> *compute_bitstring_leafset(Bipartition<T> *Bipart)
	{

		int bitstr_size = Bipart->get_bitstr_size();
		int all_leaf_size = Bipart->get_leaf_size();

		auto leafset = new BitString<T>(bitstr_size, all_leaf_size);

		for (int i = 0; i < size; i++)
			if (edges[1][i] < all_leaf_size)
				leafset->set_bit(edges[1][i]);
		return leafset;
	}


	void release_edge_form()
	{
		delete edges[0];
		delete edges[1];
		delete[] edges;
		edges = nullptr;

		// weights.release();
		// weights is needed for B2T matrix
	}

	void release_level()
	{
		delete levels;
		levels = nullptr;
	}

	int find_bipart_i(int idx)
	{
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

	Array<size_t> *get_t2b() { return t2b; };

	Array<size_t> *pop_t2b()
	{
		Array<size_t> *b_ptr = t2b;
		t2b = nullptr;
		return b_ptr;
	};

	Array<size_t> *pop_t2b_upper()
	{
		Array<size_t> *b_ptr = t2b_upper;
		t2b_upper = nullptr;
		return b_ptr;
	};

	Array<PRECISION> *get_weight()
	{
		assert(weights != nullptr);
		return weights;
	};

	Array<PRECISION> *pop_weight()
	{
		Array<PRECISION> *w_ptr = weights;
		weights = nullptr;
		return w_ptr;
	};
};

template <class T>
class TreeSet
{
private:
	int size;
	Array<TreeArray *> trees;
	Bipartition<T> *Bipart;
	TaxonList *Taxa;

	BitString<T> **LeafSets;
	BitString<T> *Uniform_LeafSet;

	bool sameleaf;

public:
	TreeSet(size_t N) : size(0), trees(Array<TreeArray *>(0, N)),
						LeafSets(nullptr), Uniform_LeafSet(nullptr),
						Bipart(nullptr), Taxa(nullptr), sameleaf(true){};

	TreeSet() : size(0), trees(Array<TreeArray *>(0, 100)),
				LeafSets(nullptr), Uniform_LeafSet(nullptr),
				Bipart(nullptr), Taxa(nullptr), sameleaf(true){};

	TreeSet(Bipartition<T> *bipart, TaxonList *taxa) : size(0), trees(Array<TreeArray *>(0, 100)),
													   LeafSets(nullptr), Uniform_LeafSet(nullptr),
													   Bipart(bipart), Taxa(taxa), sameleaf(true){};

	TreeSet(Bipartition<T> *bipart, TaxonList *taxa, bool sl) : size(0), trees(Array<TreeArray *>(0, 100)),
																LeafSets(nullptr), Uniform_LeafSet(nullptr),
																Bipart(bipart), Taxa(taxa), sameleaf(sl){};

	~TreeSet()
	{
		for (int i = 0; i < size; i++)
			delete trees[i];
		if (Uniform_LeafSet)
			delete Uniform_LeafSet;
		if (LeafSets)
		{
			for (int i = 0; i < size; i++)
				delete LeafSets[i];
			delete[] LeafSets;
		}
	}

	void push(TreeArray *tree)
	{
		trees.push(tree);
		size++;
	};

	TreeArray *ele_ptr(size_t i) { return trees[i]; };

	TreeArray &ele(size_t i) { return *(trees[i]); };

	TreeArray &ele(int i) { return *(trees[i]); };

	TreeArray &operator[](size_t i) { return this->ele(i); };

	TreeArray &operator[](int i) { return this->ele(i); };


	void ReadTree(std::string fname, size_t pos, int flag_label, bool flag_rooted, bool flag_weighted)
	{
		size_t leaf_size = Taxa->size;
		Array<int> active_levels(0, 100);
		Array<int> unlabeled(0, 100);

		std::ifstream fin;
		fin.open(fname);
		string temp;

		fin.seekg(pos);

		while (!fin.eof())
		{
			temp.clear();
			std::getline(fin, temp, ';');
			if (temp.find('(') != std::string::npos)
			{
				auto tree = new TreeArray(temp, *Taxa, active_levels, unlabeled, flag_label, flag_rooted, flag_weighted);
				tree->set_same_leaf(sameleaf);
				this->push(tree);
			}
		}

		fin.close();

		// if (sameleaf)
		// {
		// 	Uniform_LeafSet = new BitString<T>(Bipart->get_bitstr_size(), Bipart->get_leaf_size());
		// 	auto tree1 = *trees(0);
		// 	auto t2b1 = *(tree1.get_t2b());
		// 	for(int i = 0; i < tree1.get_size(); i++)
		// 	{
		// 		(*Uniform_LeafSet) |= (*Bipart)(t2b1(i));
		// 	}
		// }
	};

	void init_leafsets()
	{
		assert(LeafSets == nullptr);
		LeafSets = new BitString<T> *[this->size];
	}

	void set_leafsets(int ind, BitString<T> *leafset)
	{
		LeafSets[ind] = leafset;
	}

	void set_uniform_leafset(BitString<T> *leafset) { Uniform_LeafSet = leafset; };

	BitString<T>* compute_leafset(int i) { return trees[i]->compute_bitstring_leafset(Bipart); };

	BitString<T>* compute_uniform_leafset() { Uniform_LeafSet = trees[0]->compute_bitstring_leafset(Bipart); return Uniform_LeafSet; };

	bool issameleaf() const { return sameleaf; };

	int get_size() { return trees.get_size(); };

	int get_all_bipart_size()
	{
		int result = 0;
		for (int i = 0; i < size; i++)
			result += trees[i]->get_size();
		return result;
	}

	int get_unique_bipart_size() { return Bipart->get_size(); }

	Bipartition<T> *get_bipart_ptr() { return Bipart; };

	TaxonList *get_taxa_ptr() { return Taxa; };

	auto get_leafsets() const { return LeafSets; };

	BitString<T> *get_uniform_leafset() const { return Uniform_LeafSet; };

	void release_tree_edge_form()
	{
		for (int i = 0; i < size; i++)
			trees[i]->release_edge_form();
	}

	void release_tree()
	{
		for (int i = 0; i < size; i++)
			delete trees[i];
		trees.release();
		size = 0;
	}
};

template <class T>
class NewickEdge
{
	int label;
	PRECISION weight;
	BitString<T> leafstr;
	Array<NewickEdge *> child;

public:
	NewickEdge(const BitString<T> &leafset) : label(-1), weight(0.0), leafstr(BitString<T>(leafset)), child(Array<NewickEdge *>(0, 10)){};
	NewickEdge(const BitString<T> &leafstr, PRECISION w, int l) : label(l), weight(w), leafstr(BitString<T>(leafstr)), child(Array<NewickEdge *>(0, 10)){};

	NewickEdge(BitString<T> *leafset) : label(-1), weight(0.0), leafstr(BitString<T>(*leafset)), child(Array<NewickEdge *>(0, 10)){};
	NewickEdge(BitString<T> *leafstr, PRECISION w, int l) : label(l), weight(w), leafstr(BitString<T>(*leafstr)), child(Array<NewickEdge *>(0, 10)){};
	~NewickEdge()
	{
		for (int i = 0; i < child.get_size(); i++)
		{
			child[i]->~NewickEdge();
			delete child[i];
		}
	};

	BitString<T> &get_leafstr() { return leafstr; };

	void insert_edge(NewickEdge<T> *edge)
	{
		this->child.push(edge);
	}

	void newick_push(NewickEdge<T> *edge)
	{
		auto bitstr = edge->get_leafstr();

		// std::cout << "Current Edge\t";
		// this->leafstr.print_BitString(std::cout);
		// std::cout << "\nIncoming Edge\t";
		// edge->get_leafstr().print_BitString(std::cout);
		// std::cout << "\n\n";

		for (int i = 0; i < child.get_size(); i++)
		{
			if (is_subset(this->child[i]->get_leafstr(), bitstr))
			{
				auto sub_edge = this->child.pop(i);
				edge->insert_edge(sub_edge);
				i--;
			}
		}

		for (int i = 0; i < child.get_size(); i++)
		{
			if (is_subset(bitstr, this->child[i]->get_leafstr()))
			{
				this->child[i]->newick_push(edge);
				return;
			}
		}
		this->child.push(edge);
		return;
	}

	void newick_push(const BitString<T> &leafstr, PRECISION weight, int label) { this->newick_push(new NewickEdge(leafstr, weight, label)); };

	void Build_Newick_Tree(Bipartition<T> *Bipart, int *bipart_id, PRECISION *weights, int edge_num)
	{
		// This function assumes the current NewickEdge being the dummy edge with full leaf set.
		int label = get_taxon_id((*Bipart)[bipart_id[0]], this->leafstr);

		if (label == 0)
		{
			this->insert_edge(new NewickEdge((*Bipart)[bipart_id[0]], weights[0], -1));
			this->insert_edge(new NewickEdge(this->leafstr % (*Bipart)[bipart_id[0]], weights[0], label));
		}
		else if(label != -1)
		{
			this->insert_edge(new NewickEdge((*Bipart)[bipart_id[0]], weights[0], label));
			this->insert_edge(new NewickEdge(this->leafstr % (*Bipart)[bipart_id[0]], weights[0], -1));
		}
		else
		{
			this->insert_edge(new NewickEdge((*Bipart)[bipart_id[0]], weights[0], -1));
			this->insert_edge(new NewickEdge(this->leafstr % (*Bipart)[bipart_id[0]], weights[0], -1));
		}

		for (int i = 1; i < edge_num; i++)
		{
			label = get_taxon_id((*Bipart)[bipart_id[i]], this->leafstr);
			if (is_subset((*Bipart)[bipart_id[i]], child[0]->get_leafstr()))
				child[0]->newick_push((*Bipart)[bipart_id[i]], weights[i], label);
			else if (is_subset((*Bipart)[bipart_id[i]], child[1]->get_leafstr()))
				child[1]->newick_push((*Bipart)[bipart_id[i]], weights[i], label);
			else
			{
				auto comp_leafstr = this->leafstr % (*Bipart)[bipart_id[i]];
				if (is_subset(comp_leafstr, child[0]->get_leafstr()))
					child[0]->newick_push(comp_leafstr, weights[i], label);
				else
					child[1]->newick_push(comp_leafstr, weights[i], label);
			}
		}
	}

	void Print_Newick_Edge(std::ostream &out, const TaxonList &taxa)
	{
		// This function assumes the current NewickEdge being the actural edge with subset of leaf nodes.
		if (this->label != -1)
			out << taxa(this->label) << ":" << this->weight;
		else
		{
			out << "(";
			this->child[0]->Print_Newick_Edge(out, taxa);
			for (int i = 1; i < this->child.get_size(); i++)
			{
				out << ",";
				this->child[i]->Print_Newick_Edge(out, taxa);
			}
			out << "):" << this->weight;
		}
	}

	void Print_Newick_Tree(std::ostream &out, const TaxonList &taxa)
	{
		// This function assumes the current NewickEdge being the dummy edge with full leaf set.
		out << "(";
		this->child[0]->Print_Newick_Edge(out, taxa);
		this->child[1]->Print_Newick_Edge(out, taxa);
		out << "):0.0;\n";
	};
};

// template <class T>
// class TreeDirectedEdge
// {
// public:
// 	int num_next_edges;
// 	int num_population;
// 	PRECISION weight;
// 	BitString<T> bitstr;
// 	TreeDirectedEdge **next_edges;

// 	TreeDirectedEdge() : num_next_edges(0), num_population(0), bitstr(), next_edges(nullptr){};
// 	TreeDirectedEdge(const BitString<T> &bs, PRECISION w)
// 	{
// 		weight = w;
// 		bitstr = bs;

// 		num_population = bs.popcount();
// 		num_next_edges = 0;

// 		if (num_population == 1)
// 			next_edges = nullptr;
// 		else
// 		{
// 			next_edges = new TreeDirectedEdge *[num_population];
// 			for (int i = 0; i < num_population; i++)
// 				next_edges[i] = nullptr;
// 		}
// 	}
// 	~TreeDirectedEdge()
// 	{
// 		for (int i = 0; i < num_next_edges; i++)
// 			if (next_edges[i] != nullptr)
// 			{
// 				std::cout << "Error in deconstructing TreeDirectedEdge object: non-empty adjoint edge at position " << i << ".\n";
// 				throw(1);
// 			}

// 		delete[] this->next_edges;
// 	}
// };

// template <class T>
// class TreeLinkedList
// {
// public:
// 	TreeDirectedEdge<T> *full_taxa;

// 	TreeLinkedList() : full_taxa(nullptr){};
// 	TreeLinkedList(const BitString<T> &bs, PRECISION w) : full_taxa(new TreeDirectedEdge<T>(bs, w)){};

// 	TreeLinkedList(const BitString<T> &full, BitString<T> **bitstrs, PRECISION *weights, int num_edges) : full_taxa(new TreeDirectedEdge<T>(full, 0.0))
// 	{
// 		for (int i = 0; i < num_edges; i++)
// 			insert_edge(*(bitstrs[i]), weights[i]);
// 	}
// 	~TreeLinkedList()
// 	{
// 		release_edge(full_taxa);
// 	}

// 	TreeLinkedList(const BitString<T> &full, int *bitstrs, PRECISION *weights, int num_edges, Bipartition<T> *Bipart) : full_taxa(new TreeDirectedEdge<T>(full, 0.0))
// 	{
// 		for (int i = 0; i < num_edges; i++)
// 			insert_edge((*Bipart)[bitstrs[i]], weights[i]);
// 	}

// 	void release_edge(TreeDirectedEdge<T> *cur_edge)
// 	{

// 		for (int i = cur_edge->num_next_edges - 1; i >= 0; i++)
// 		{
// 			release_edge(cur_edge->next_edges[i]);
// 			cur_edge->next_edges[i] = nullptr;
// 			cur_edge->num_next_edges = cur_edge->num_next_edges - 1;
// 		}

// 		delete cur_edge;
// 	}

// 	TreeDirectedEdge<T> *locate_edge(TreeDirectedEdge<T> *new_edge, TreeDirectedEdge<T> *cur_edge)
// 	{
// 		if (is_subset(new_edge->bitstr, cur_edge->bitstr))
// 		{
// 			if (cur_edge->num_next_edges == 0)
// 				return cur_edge;
// 			else
// 			{
// 				auto found_edge = nullptr;
// 				for (int i = 0; i < cur_edge->num_next_edges; i++)
// 				{
// 					found_edge = locate_edge(new_edge, cur_edge->next_edges[i]);
// 					if (found_edge != nullptr)
// 						return found_edge;
// 				}
// 				return nullptr;
// 			}
// 		}
// 		else
// 			return nullptr;
// 	}

// 	void insert_edge(const BitString<T> &bs, PRECISION w)
// 	{
// 		auto e = new TreeDirectedEdge<T>(bs, w);
// 		auto c = new TreeDirectedEdge<T>(full_taxa->bitstr % bs, w);

// 		if (full_taxa->num_next_edges == 0)
// 		{
// 			full_taxa->next_edges[0] = e;
// 			full_taxa->next_edges[1] = c;
// 			full_taxa->num_next_edges == 2;
// 		}
// 		else
// 		{
// 			auto found_edge = locate_edge(e, full_taxa);
// 			auto new_edge = e;

// 			if (found_edge == nullptr)
// 			{
// 				found_edge = locate_edge(c, full_taxa);
// 				new_edge = c;
// 				delete e;
// 			}
// 			else
// 				delete c;

// 			for (int i = 0; i < found_edge->num_next_edges; i++)
// 			{
// 				if (is_subset(found_edge->next_edges[i]->bitstr, new_edge->bitstr))
// 				{
// 					new_edge->next_edges[new_edge->num_next_edges++] = found_edge->next_edges[i];
// 					for (int j = i + 1; j < found_edge->num_next_edges; j++)
// 						found_edge->next_edges[j - 1] = found_edge->next_edges[j];
// 					found_edge->num_next_edges = found_edge->num_next_edges - 1;
// 					i--;
// 				}
// 			}
// 			found_edge->next_edges[found_edge->num_next_edges++] = new_edge;
// 		}
// 	}

// 	void print_NewickStr(ostream &out, TreeDirectedEdge<T> *cur_edge, const TaxonList &taxa)
// 	{
// 		if (cur_edge->num_population == 1)
// 		{
// 			int taxon_id = cur_edge->bitstr.get_taxon_id();
// 			out << taxa.Ind2Taxon(taxon_id) << ':' << cur_edge->weight;
// 		}
// 		else
// 		{
// 			out << '(';
// 			for (int i = 0; i < cur_edge->num_next_edges - 1; i++)
// 			{
// 				print_NewickStr(out, cur_edge->next_edges[i], taxa);
// 				out << ',';
// 			}
// 			print_NewickStr(out, cur_edge->next_edges[cur_edge->num_next_edges - 1], taxa);
// 			out << "):" << cur_edge->weight;
// 		}
// 	}

// 	void print_NewickStr(ostream &out, const TaxonList &taxa)
// 	{
// 		out << '(';
// 		print_NewickStr(out, full_taxa->next_edges[0], taxa);
// 		out << ',';
// 		print_NewickStr(out, full_taxa->next_edges[0], taxa);
// 		out << "):" << full_taxa->weight << '\n';
// 	}
// };