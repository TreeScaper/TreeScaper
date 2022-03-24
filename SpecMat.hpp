#pragma once
#include <iostream>
#include <fstream>
//#include<string>
#include <cstring>
#include <cassert>
#include "array.hpp"

namespace SpecMat
{
	template <class T>
	class LowerTri
	{
	protected:
		unsigned dim;
		T *data;
		T **row;
	public:
		LowerTri()
		{
			dim = 0;
			data = nullptr;
			row = nullptr;
		}
		LowerTri(unsigned d) : dim(d)
		{
			data = new T[((1 + dim) * dim) / 2];
			memset(data, 0, ((1 + dim) * dim) / 2 * sizeof(T));
			row = nullptr;
		}

		LowerTri(unsigned d, bool row_ptr) : LowerTri(d)
		{
			if (row_ptr)
				this->form_row_ptr();
		}
		//LowerTri(const mat::matrix &src):LowerTri(src.getRow()) {
		//	for (unsigned i = 0; i < dim; i++) {
		//		memcpy(data + map(i), src(i), (i + 1) * sizeof(PRECISION));
		//	}
		//}
		LowerTri(const LowerTri<T> &src) : LowerTri(src.dim)
		{
			data = new T[((1 + dim) * dim) / 2];
			int s = ((1 + dim) * dim) / 2;
			for (auto i = 0; i < s; i++)
				data[i] = src.data[i];
			if (src.row)
				this->form_row_ptr();
		}
		LowerTri(unsigned d, T *vec) : dim(d), data(vec){};
		~LowerTri() 
		{ 
			if (row)
				delete[] row;
			delete[] data; 
		};

		LowerTri(unsigned d, std::iostream& file) : LowerTri(d)
		{
			char* line = new char[d * 10];
			this->form_row_ptr();
			for(int i = 0; i < d; i++){
				for(int j = 0; j <= i; j++)
					file >> this->row[i][j];
				file.getline(line, d * 10);
			}
			delete line;
		}
		void form_row_ptr(){
			if (row)
				return;
			row = new PRECISION*[dim];
			for(auto i = 0; i < dim; i++)
				row[i] = data + map(i);
		}
		



		// Default fast access though extra row. Return error if row not initialized.
		T *operator[](unsigned i) { 
			return row[i]; 
			};
		T *ele(unsigned i) const { 
			return row[i]; 
			};
		T &operator()(unsigned i, unsigned j) { 
			// assert(i >= j);
			return row[i][j];
		};
		T &ele(unsigned i, unsigned j) const { 
			// assert(i >= j);
			return row[i][j];
		};

		// Slow access that recompute start of each row 
		T *ele_raw(unsigned i) { return data + map(i); };
		T *ele_raw_const(unsigned i) const { return data + map(i); };
		T *ele_raw(unsigned i, unsigned j) { return data[map(i, j)]; };
		T *ele_raw_const(unsigned i, unsigned j) const { return data[map(i, j)]; };



		LowerTri<T> & operator+=(const LowerTri<T> &rhs)
		{
			int size = ((1 + dim) * dim) / 2;
			for(int i = 0; i < size; i++)
				this->data[i] += rhs.data[i];
			return *this;
		}

		LowerTri<T> & operator-=(const LowerTri<T> &rhs)
		{
			int size = ((1 + dim) * dim) / 2;
			for(int i = 0; i < size; i++)
				this->data[i] -= rhs.data[i];
			return *this;
		}

		const LowerTri<T> operator+(const LowerTri<T> &rhs) const
		{
			return LowerTri<T>(*this) += rhs;
		}

		const LowerTri<T> operator-(const LowerTri<T> &rhs) const
		{
			return LowerTri<T>(*this) -= rhs;
		}

		void swap(LowerTri<T> &rhs)
		{
			using std::swap;
			swap(this->data, rhs.data);
			swap(this->dim, rhs.dim);
		}
		LowerTri<T> &operator=(const LowerTri &rhs)
		{
			SpecMat::LowerTri<T> temp(rhs);
			swap(temp);
			return (*this);
		}

		unsigned dimension() const { return dim; };
		unsigned map(unsigned i) const { return ((1 + i) * i) / 2; };
		unsigned map(unsigned i, unsigned j) const
		{
			//assert(j <= i);
			return (i > j) ? ((1 + i) * i) / 2 + j : ((1 + j) * j) / 2 + i;
		};

		struct Row_Iterator 
		{
    		using iterator_category = std::forward_iterator_tag;
    		using difference_type   = std::ptrdiff_t;
			using self_type = Row_Iterator;
			using value_type        = T;
    		using pointer           = T*;  // or also value_type*
    		using reference         = T&;  // or also value_type&

			Row_Iterator(pointer ptr, unsigned r, unsigned c) : m_ptr(ptr), row(r), col(c){};


			reference operator*() const { return *m_ptr; }
    		pointer operator->() { return m_ptr; }

    		// Prefix increment
    		self_type& operator++() {
				m_ptr += (row > col) ? 1 : (1 + col);
				col++;
				return *this; 
			}  

    		// Postfix increment
    		self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }

    		friend bool operator== (const self_type& a, const self_type& b) { return a.row == b.row && a.col == b.col; };
    		friend bool operator!= (const self_type& a, const self_type& b) { return a.row != b.row || a.col != b.col; };

			unsigned get_col() { return col; };

		private:

    		pointer m_ptr;
			unsigned row;
			unsigned col;
		};

		Row_Iterator row_begin(unsigned r) { return Row_Iterator((*this)[r], r, 0); };
		Row_Iterator row_end(unsigned r) { return Row_Iterator(nullptr, r, dim); };

		PRECISION norm_Frobenius_square()
		{
			PRECISION temp = 0;
			unsigned length = ((1 + dim) * dim) / 2;
			for (unsigned i = 0; i < length; i++)
				temp += data[i] * data[i];
			return temp;
		}

		void print(std::ostream &output = std::cout)
		{
			unsigned ind = 0;
			for (unsigned i = 0; i < dim; i++)
			{
				for (unsigned j = 0; j <= i; j++)
				{
					output << data[ind] << '\t';
					ind++;
				}
				output << '\n';
			}
			output << '\n';
		}

		T *get_vec() { return data; };

		friend PRECISION diff_norm(const LowerTri<PRECISION> &lhs, const LowerTri<PRECISION> &rhs, MAT_NORM_TYPE mnt)
		{
			PRECISION ans = 0, dif = 0;
			unsigned len = (lhs.dim * (lhs.dim + 1)) / 2;
			switch (mnt)
			{
			case M_FOR_NORM:
				for (auto i = 0; i < len; i++)
				{
					dif = lhs.data[i] - rhs.data[i];
					ans += dif * dif;
				}
				ans = sqrt(ans);
				break;
			default:
				std::cout << "Norm function not found!\n";
			}
			return ans;
		}

		PRECISION norm(MAT_NORM_TYPE mnt)
		{
			PRECISION ans = 0;
			unsigned len = (dim * (dim + 1)) / 2;
			switch (mnt)
			{
			case M_FOR_NORM:
				for (auto i = 0; i < len; i++)
					ans += data[i] * data[i];
				ans = sqrt(ans);
				break;
			default:
				std::cout << "Norm function not found!\n";
			}
			return ans;
		}
	};

	// template <class T>
	// class LowerTri
	// {
	// protected:
	// 	unsigned dim;
	// 	PRECISION *data;

	// public:
	// 	LowerTri()
	// 	{
	// 		dim = 0;
	// 		data = nullptr;
	// 	}
	// 	LowerTri(unsigned d) : dim(d)
	// 	{
	// 		data = new T[((1 + dim) * dim) / 2];
	// 		memset(data, 0, ((1 + dim) * dim) / 2 * sizeof(T));
	// 	}
	// 	//LowerTri(const mat::matrix &src):LowerTri(src.getRow()) {
	// 	//	for (unsigned i = 0; i < dim; i++) {
	// 	//		memcpy(data + map(i), src(i), (i + 1) * sizeof(PRECISION));
	// 	//	}
	// 	//}
	// 	LowerTri(const LowerTri<T> &src) : dim(src.dimension())
	// 	{
	// 		data = new T[((1 + dim) * dim) / 2];
	// 		int s = ((1 + dim) * dim) / 2;
	// 		for (auto i = 0; i < s; i++)
	// 			data[i] = src.data[i];
	// 	}
	// 	LowerTri(unsigned d, T *vec) : dim(d), data(vec){};
	// 	~LowerTri() { delete[] data; };

	// 	T *operator[](unsigned i) { return data + map(i); };
	// 	//PRECISION* operator()(unsigned i) const { return data + map(i); };
	// 	T &operator()(unsigned i, unsigned j) { return data[map(i, j)]; };
	// 	T &ele(unsigned i, unsigned j) const { return data[map(i, j)]; };
	// 	T *ele(unsigned i) const { return data + map(i); };
	// 	LowerTri<T> operator+(const LowerTri<T> &rhs) const
	// 	{
	// 		LowerTri<T> result(rhs);
	// 		for (unsigned i = 0; i < dim; i++)
	// 		{
	// 			for (unsigned j = 0; j <= i; j++)
	// 				result[i][j] = result[i][j] + data[this->map(i, j)];
	// 		}
	// 		return result;
	// 	}
	// 	LowerTri<T> operator-(const LowerTri<T> &rhs) const
	// 	{
	// 		LowerTri<T> result(rhs);
	// 		for (unsigned i = 0; i < dim; i++)
	// 		{
	// 			for (unsigned j = 0; j <= i; j++)
	// 				result[i][j] = data[this->map(i, j)] - result[i][j];
	// 		}
	// 		return result;
	// 	}
	// 	void swap(LowerTri<T> &rhs)
	// 	{
	// 		using std::swap;
	// 		swap(this->data, rhs.data);
	// 		swap(this->dim, rhs.dim);
	// 	}
	// 	LowerTri<T> &operator=(const LowerTri &rhs)
	// 	{
	// 		SpecMat::LowerTri<T> temp(rhs);
	// 		swap(temp);
	// 		return (*this);
	// 	}

	// 	unsigned dimension() const { return dim; };
	// 	unsigned map(unsigned i) const { return ((1 + i) * i) / 2; };
	// 	unsigned map(unsigned i, unsigned j) const
	// 	{
	// 		//assert(j <= i);
	// 		return (i > j) ? ((1 + i) * i) / 2 + j : ((1 + j) * j) / 2 + i;
	// 	};
	// 	PRECISION norm_Frobenius_square()
	// 	{
	// 		PRECISION temp = 0;
	// 		unsigned length = ((1 + dim) * dim) / 2;
	// 		for (unsigned i = 0; i < length; i++)
	// 			temp += data[i] * data[i];
	// 		return temp;
	// 	}

	// 	void print(std::ostream &output = std::cout)
	// 	{
	// 		unsigned ind = 0;
	// 		for (unsigned i = 0; i < dim; i++)
	// 		{
	// 			for (unsigned j = 0; j <= i; j++)
	// 			{
	// 				output << data[ind] << '\t';
	// 				ind++;
	// 			}
	// 			output << '\n';
	// 		}
	// 		output << '\n';
	// 	}

	// 	T *get_vec() { return data; };

	// 	friend PRECISION diff_norm(const LowerTri<PRECISION> &lhs, const LowerTri<PRECISION> &rhs, MAT_NORM_TYPE mnt)
	// 	{
	// 		PRECISION ans = 0, dif = 0;
	// 		unsigned len = (lhs.dim * (lhs.dim + 1)) / 2;
	// 		switch (mnt)
	// 		{
	// 		case M_FOR_NORM:
	// 			for (auto i = 0; i < len; i++)
	// 			{
	// 				dif = lhs.data[i] - rhs.data[i];
	// 				ans += dif * dif;
	// 			}
	// 			ans = sqrt(ans);
	// 			break;
	// 		default:
	// 			std::cout << "Norm function not found!\n";
	// 		}
	// 		return ans;
	// 	}

	// 	PRECISION norm(MAT_NORM_TYPE mnt)
	// 	{
	// 		PRECISION ans = 0;
	// 		unsigned len = (dim * (dim + 1)) / 2;
	// 		switch (mnt)
	// 		{
	// 		case M_FOR_NORM:
	// 			for (auto i = 0; i < len; i++)
	// 				ans += data[i] * data[i];
	// 			ans = sqrt(ans);
	// 			break;
	// 		default:
	// 			std::cout << "Norm function not found!\n";
	// 		}
	// 		return ans;
	// 	}
	// };

	class UpperTri
	{
	protected:
		unsigned dim;
		PRECISION *data;

	public:
		UpperTri()
		{
			dim = 0;
			data = new PRECISION[((1 + dim) * dim) / 2];
			memset(data, 0, ((1 + dim) * dim) / 2 * sizeof(PRECISION));
		}
		UpperTri(unsigned d) : dim(d)
		{
			data = new PRECISION[((1 + dim) * dim) / 2];
			memset(data, 0, ((1 + dim) * dim) / 2 * sizeof(PRECISION));
		}
		//UpperTri(const mat::matrix &src) :dim(src.getRow()) {
		//	for (unsigned i = 0; i < dim; i++) {
		//		memcpy(data + map(i), src(i), (dim - i) * sizeof(PRECISION));
		//	}
		//}
		UpperTri(const UpperTri &src) : UpperTri(src.dimension())
		{
			memcpy(data, src(0), (((1 + dim) * dim) / 2) * sizeof(PRECISION));
		}

		unsigned map(unsigned i) const { return ((dim + dim + 1 - i) * i) / 2; };
		PRECISION *operator[](unsigned i) { return data + map(i); };
		PRECISION *operator()(unsigned i) const { return data + map(i); };
		PRECISION operator()(unsigned i, unsigned j) const { return data[map(i) + j]; };
		void swap(UpperTri &rhs)
		{
			using std::swap;
			swap(this->data, rhs.data);
			swap(this->dim, rhs.dim);
		}
		UpperTri &operator=(UpperTri &rhs)
		{
			swap(rhs);
			return (*this);
		}

		unsigned dimension() const { return dim; };

		void print()
		{
			for (unsigned i = 0; i < dim; i++)
			{
				for (unsigned j = 0; j <= dim - i; j++)
				{
					std::cout << data[map(i) + j] << '\t';
				}
				std::cout << '\n';
			}
			std::cout << '\n';
		}
	};

	//class Banded_d {
	//protected:
	//	unsigned dim;
	//	unsigned upper;
	//	unsigned lower;
	//	unsigned len;
	//	PRECISION* data;
	//public:
	//	Banded_d() :dim(5), upper(0), lower(0) {
	//		len = dim + ((dim - 1 + dim - upper)*upper) / 2 + ((dim - 1 + dim - lower)*lower) / 2;
	//		data = new PRECISION[len]; memset(data, 0, len * sizeof(PRECISION));
	//	};
	//	Banded_d(unsigned d, unsigned u, unsigned l) :dim(d), upper(u), lower(l) {
	//		len = dim + ((dim - 1 + dim - upper)*upper) / 2 + ((dim - 1 + dim - lower)*lower) / 2;
	//		data = new PRECISION[len]; memset(data, 0, len * sizeof(PRECISION));
	//	}
	//	Banded_d(unsigned d, unsigned u, unsigned l, unsigned length) :dim(d), upper(u), lower(l), len(length) {
	//		data = new PRECISION[len]; memset(data, 0, len * sizeof(PRECISION));
	//	}
	//	Banded_d(PRECISION *v, unsigned d, unsigned u, unsigned l) :dim(d), upper(u), lower(l) {
	//		len = dim + ((dim - 1 + dim - upper)*upper) / 2 + ((dim - 1 + dim - lower)*lower) / 2;
	//		data = new PRECISION[len]; memcpy(data, v, len * sizeof(PRECISION));
	//	}
	//	Banded_d(PRECISION **v, unsigned d, unsigned u, unsigned l) :dim(d), upper(u), lower(l) {
	//		len = dim + ((dim - 1 + dim - upper)*upper) / 2 + ((dim - 1 + dim - lower)*lower) / 2;
	//		data = new PRECISION[len];
	//		PRECISION* current = data;
	//		for (unsigned i = 0; i < u + l + 1; i++) {
	//			if(i < u){
	//				memcpy(current, v[i], (dim - u + i) * sizeof(PRECISION));
	//				current = current + (dim - u + i);
	//			}
	//			else if (i == u) {
	//				memcpy(current, v[i], dim * sizeof(PRECISION));
	//				current = current + dim;
	//			}
	//			else {
	//				memcpy(current, v[i], (dim - i + u) * sizeof(PRECISION));
	//				current = current + (dim - i + u);
	//			}
	//		}
	//	}
	//	Banded_d(const Banded_d &src) :dim(src.dim), upper(src.upper), lower(src.lower), len(src.len) { data = new PRECISION[len]; memcpy(data, src.start(), len * sizeof(PRECISION)); };
	//	~Banded_d() { delete[] data; };

	//	PRECISION* operator[](int i) {
	//		if (i >= 0)
	//			return data + ((2 * dim - upper - i - 1)*(upper - i)) / 2;
	//		else
	//			return data + ((2 * dim - upper - 1)*(upper)) / 2 - ((2 * dim + i + 1)*i) / 2;
	//	}
	//	void swap(Banded_d &rhs) {
	//		using std::swap;
	//		swap(this->data, rhs.data);
	//		swap(this->dim, rhs.dim);
	//		swap(this->len, rhs.len);
	//		swap(this->upper, rhs.upper);
	//		swap(this->lower, rhs.lower);
	//	}
	//	Banded_d& operator=(Banded_d & rhs) { swap(rhs); return (*this); }

	//	unsigned dimension() const { return dim; };
	//	PRECISION* start() const { return data; };

	//};

	//class UpperBiD_d: public Banded_d {
	//public:
	//	UpperBiD_d() : Banded_d(5,1,0,9) {};
	//	UpperBiD_d(unsigned d) :Banded_d(d, 1, 0, 2 * d + 1) {};
	//	UpperBiD_d(PRECISION* v, unsigned d) :Banded_d(d, 1, 0, 2 * d + 1) { memcpy(data, v, len * sizeof(PRECISION)); };
	//	UpperBiD_d(PRECISION* v_u, PRECISION*v_d, unsigned d) :Banded_d(d, 1, 0, 2 * d + 1) {
	//		memcpy(data, v_u, (dim - 1) * sizeof(PRECISION));
	//	};
	//	UpperBiD_d(const UpperBiD_d &src) :Banded_d(src) {};

	//	//using Banded_d::operator[];
	//	PRECISION* operator[](unsigned i) {
	//		if (i == 1)
	//			return data;
	//		else if (i == 0)
	//			return data + (dim - 1);
	//		else {
	//			std::cout << "Error: Overflow in UpperBid_d. Try to access position not on main or first super diagonal.\n";
	//			exit(1);
	//		}
	//	}
	//	void swap(UpperBiD_d &rhs) {
	//		using std::swap;
	//		swap(this->data, rhs.data);
	//		swap(this->dim, rhs.dim);
	//		swap(this->len, rhs.len);
	//		swap(this->upper, rhs.upper);
	//		swap(this->lower, rhs.lower);
	//	}
	//	UpperBiD_d& operator=(UpperBiD_d & rhs) { swap(rhs); return (*this); }

	//};

	//class Toeplitz {
	//protected:
	//	unsigned dim;
	//	unsigned length;
	//	PRECISION* data;
	//	//Important: data store the first row of the extended Circulant matrix
	//	// 3 by 3 Toeplitz with data = {3, 4, 5, 1, 2}:
	//	// % 3 4 5 | 1 2 %
	//	// % 2 3 4 | 5 1 %
	//	// % 1 2 3 | 4 5 %
	//	// % - - - - - - %
	//	// % 5 1 2 | 3 4 %
	//	// % 4 5 1 | 2 3 %
	//public:
	//	Toeplitz() {
	//		dim = 5;
	//		length = dim * 2 - 1;
	//		data = new PRECISION[length];
	//		memset(data, 0, length * sizeof(PRECISION));
	//	}
	//	Toeplitz(unsigned d) : dim(d), length(2 * dim - 1) {
	//		data = new PRECISION[length];
	//		memset(data, 0, length * sizeof(PRECISION));
	//	}
	//	Toeplitz(PRECISION* src, unsigned d) : dim(d), length(2 * d - 1) {
	//		data = new PRECISION[length];
	//		memcpy(data, src, length * sizeof(PRECISION));
	//	}
	//	Toeplitz(const vec::vector &src, unsigned d) : Toeplitz(src.start(), d) {};
	//	Toeplitz(const Toeplitz &src) : dim(src.dimension()), length(2 * src.dimension() - 1) {
	//		data = new PRECISION[length];
	//		memcpy(data, src.start(), length * sizeof(PRECISION));
	//	}
	//	virtual ~Toeplitz() { delete[] data; };

	//	void swap(Toeplitz &rhs) {
	//		using std::swap;
	//		swap(this->data, rhs.data);
	//		swap(this->dim, rhs.dim);
	//		swap(this->length,rhs.length);
	//	}
	//	Toeplitz& operator=(Toeplitz & rhs) { swap(rhs); return (*this); }

	//	bool SymmetricQ()const {
	//		for (unsigned i = 0; i < dim; i++) {
	//			if (data[i + 1] != data[length - 1 - i])
	//				return false;
	//		}
	//		return true;
	//	}
	//	virtual mat::matrix dense() const {
	//		mat::matrix result(dim, dim);
	//		for (unsigned i = 0; i < dim; i++) {
	//			for (unsigned j = 0; j < dim; j++)
	//				result[i][j] = data[(length + j - i) % length];
	//		}
	//		return result;
	//	};
	//	vec::vector FactoredForm() const {
	//		PRECISION *temp = new PRECISION[dim];
	//		temp[0] = sqrt(data[0]);
	//		for (unsigned i = 1; i < dim; i++)
	//			temp[i] = data[i] / temp[0];
	//		return vec::vector(temp, dim);
	//	};
	//	unsigned dimension() const { return dim; };
	//	PRECISION* start() const { return data; };
	//};

	//class Circulant :public Toeplitz {
	//public:
	//	Circulant() {
	//		dim = 5;
	//		length = dim;
	//		data = new PRECISION[length];
	//		memset(data, 0, length * sizeof(PRECISION));
	//	}
	//	Circulant(unsigned d) {
	//		dim = d;
	//		length = dim;
	//		data = new PRECISION[length];
	//		memset(data, 0, length * sizeof(PRECISION));
	//	}
	//	Circulant(PRECISION* src, unsigned d) : Circulant(d) { memcpy(data, src, length * sizeof(PRECISION)); }
	//	Circulant(const vec::vector &src) : Circulant(src.start(), src.length()) {};
	//};

	//class HypRotation {
	//public:
	//	unsigned dim;
	//	unsigned i;
	//	unsigned j;
	//	PRECISION ch;
	//	PRECISION sh;

	//	HypRotation() :dim(2), i(0), j(1), ch(1), sh(0) {};
	//	HypRotation(PRECISION rho) :dim(2), i(0), j(1), ch(cosh(rho)), sh(sinh(rho)) {};
	//	HypRotation(PRECISION c, PRECISION s) :dim(2), i(0), j(1), ch(c), sh(s) {};
	//	HypRotation(unsigned d, unsigned i1, unsigned i2, PRECISION rho) :dim(d), i(i1), j(i2), ch(cosh(rho)), sh(sinh(rho)) {};
	//	HypRotation(unsigned d, unsigned i1, unsigned i2, PRECISION c, PRECISION s) :dim(d), i(i1), j(i2), ch(c), sh(s) {};
	//	HypRotation(HypRotation &src) :dim(src.dim), i(src.i), j(src.j), ch(src.ch), sh(src.sh) {};
	//	~HypRotation() {}

	//	void swap(HypRotation &rhs) {
	//		using std::swap;
	//		swap(this->dim, rhs.dim);
	//		swap(this->i, rhs.i);
	//		swap(this->j, rhs.j);
	//		swap(this->ch, rhs.ch);
	//		swap(this->sh, rhs.sh);
	//	}
	//	HypRotation& operator=(HypRotation & rhs) { swap(rhs); return (*this); }

	//	void action(PRECISION *v1, PRECISION *v2, unsigned n) {
	//		PRECISION* temp = new PRECISION[n];
	//		memcpy(temp, v1, n * sizeof(PRECISION));
	//		for (unsigned i = 0; i < n; i++) {
	//			v1[i] = v1[i] * ch + v2[i] * sh;
	//			v2[i] = temp[i] * sh + v2[i] * ch;
	//		}
	//	}

	//	void action(mat::matrix& src, unsigned start, unsigned M_LEN, char c) {
	//		PRECISION* temp = new PRECISION[M_LEN];
	//		if (c == 'L')
	//			memcpy(temp, src(i) + start, M_LEN * sizeof(PRECISION));
	//		else if (c == 'R') {
	//			for (unsigned k = 0; k < M_LEN; k++)
	//				temp[k] = src(k + start, i);
	//		}
	//		else {
	//			std::cout << "Error: the control variable in action() of HypRotation is incorrect.\n";
	//			exit(1);
	//		}
	//		for (unsigned k = 0; k < M_LEN; k++) {
	//			if (c == 'L') {
	//				src[i][k] = ch * src[i][k] + sh * src[j][k];
	//				src[j][k] = sh * temp[k] + ch * src[j][k];
	//			}
	//			else if (c == 'R') {
	//				src[k][i] = ch * src[k][i] + sh * src[k][j];
	//				src[k][j] = ch * src[k][j] + sh * temp[k];
	//			}
	//		}
	//	};

	//	void action(mat::matrix& src, unsigned start, char c) {
	//		if (c == 'L')
	//			action(src, start, src.getCol() - start, c);
	//		else if (c == 'R')
	//			action(src, start, src.getRow() - start, c);
	//		else {
	//			std::cout << "Error: the control variable is incorrect in action() of HypRotation.\n";
	//			exit(1);
	//		}
	//	}

	//	void action(vec::vector &src, unsigned start) {
	//		PRECISION temp = 0;
	//		unsigned ii = i + start;
	//		unsigned jj = j + start;
	//		temp = src[ii];
	//		src[ii] = src[ii] * ch + src[jj] * sh;
	//		src[jj] = temp * sh + src[jj] * ch;
	//	}
	//	void action(vec::vector &src) { action(src, 0); };
	//	void action_transpose(vec::vector &src, unsigned start) {
	//		PRECISION temp = 0;
	//		unsigned ii = i + start;
	//		unsigned jj = j + start;
	//		temp = src[ii];
	//		src[ii] = src[ii] * ch + src[jj] * sh;
	//		src[jj] = src[jj] * ch + temp * sh;
	//	}
	//	void action_transpose(vec::vector &src) { action(src, 0); };
	//};

	//class GivRotation {
	//public:
	//	unsigned dim;
	//	unsigned i;
	//	unsigned j;
	//	PRECISION co;
	//	PRECISION si;

	//	GivRotation() :dim(2), i(0), j(1), co(1), si(0) {};
	//	GivRotation(PRECISION rho) :dim(2), i(0), j(1), co(cos(rho)), si(sin(rho)) {};
	//	GivRotation(PRECISION c, PRECISION s) :dim(2), i(0), j(1), co(c), si(s) {};
	//	GivRotation(unsigned d, unsigned i1, unsigned i2, PRECISION rho) :dim(d), i(i1), j(i2), co(cos(rho)), si(sin(rho)) {};
	//	GivRotation(unsigned d, unsigned i1, unsigned i2, PRECISION c, PRECISION s) :dim(d), i(i1), j(i2), co(c), si(s) {};
	//	GivRotation(GivRotation &src) :dim(src.dim), i(src.i), j(src.j), co(src.co), si(src.si) {};
	//	~GivRotation() {}

	//	void swap(GivRotation &rhs) {
	//		using std::swap;
	//		swap(this->dim, rhs.dim);
	//		swap(this->i, rhs.i);
	//		swap(this->j, rhs.j);
	//		swap(this->co, rhs.co);
	//		swap(this->si, rhs.si);
	//	}
	//	GivRotation& operator=(GivRotation & rhs) { swap(rhs); return (*this); }

	//	void action(PRECISION *v1, PRECISION *v2, unsigned n) {
	//		PRECISION* temp = new PRECISION[n];
	//		memcpy(temp, v1, n * sizeof(PRECISION));
	//		for (unsigned i = 0; i < n; i++) {
	//			v1[i] = v1[i] * co - v2[i] * si;
	//			v2[i] = temp[i] * si + v2[i] * co;
	//		}
	//	}

	//	void action(mat::matrix& src, unsigned start,unsigned M_LEN,char c) {
	//		PRECISION* temp = new PRECISION[M_LEN];
	//		if (c == 'L')
	//			memcpy(temp, src(i) + start, M_LEN * sizeof(PRECISION));
	//		else if (c == 'R') {
	//			for (unsigned k = 0; k < M_LEN; k++)
	//				temp[k] = src(k + start, i);
	//		}
	//		else {
	//			std::cout << "Error: the control variable in action() of GivRotation is incorrect.\n";
	//			exit(1);
	//		}
	//		for (unsigned k = 0; k < M_LEN; k++) {
	//			if(c == 'L'){
	//				src[i][k] = co*src[i][k] - si*src[j][k];
	//				src[j][k] = si*temp[k] + co*src[j][k];
	//			}
	//			else if (c == 'R') {
	//				src[k][i] = co*src[k][i] + si*src[k][j];
	//				src[k][j] = co*src[k][j] - si*temp[k];
	//			}
	//		}
	//	};

	//	void action(mat::matrix& src, unsigned start, char c) {
	//		if (c == 'L')
	//			action(src, start, src.getCol() - start, c);
	//		else if (c == 'R')
	//			action(src, start, src.getRow() - start, c);
	//		else {
	//			std::cout << "Error: the control variable is incorrect in action() of GivRotation.\n";
	//			exit(1);
	//		}
	//	}

	//	void action(vec::vector &src, unsigned start) {
	//		PRECISION temp = 0;
	//		unsigned ii = i + start;
	//		unsigned jj = j + start;
	//		temp = src[ii];
	//		src[ii] = src[ii] * co - src[jj] * si;
	//		src[jj] = temp*si + src[jj] * co;
	//	}
	//	void action(vec::vector &src) { action(src, 0); };
	//	void action_transpose(vec::vector &src, unsigned start){
	//		PRECISION temp = 0;
	//		unsigned ii = i + start;
	//		unsigned jj = j + start;
	//		temp = src[ii];
	//		src[ii] = src[ii] * co + src[jj] * si;
	//		src[jj] = src[jj] * co - temp * si;
	//}
	//	void action_transpose(vec::vector &src) { action(src, 0); };
	//};

	//class Permutation {
	//public:
	//	unsigned i;
	//	unsigned j;
	//	unsigned dim;
	//	Permutation() { dim = 2; i = 0; j = 1; };
	//	Permutation(unsigned d, unsigned ii, unsigned jj) :dim(d), i(ii), j(jj) {};
	//	~Permutation() {};

	//	void action_l(mat::matrix &src) {
	//		unsigned d = src.getCol();
	//		PRECISION *temp = new PRECISION[d];
	//		memcpy(temp, src[i], d * sizeof(PRECISION));
	//		memcpy(src[i], src[j], d * sizeof(PRECISION));
	//		memcpy(src[j], temp, d * sizeof(PRECISION));
	//	}
	//	void action_r(mat::matrix &src) {
	//		unsigned d = src.getRow();
	//		PRECISION *temp = new PRECISION[d];
	//		for (unsigned k = 0; k < d; k++)
	//			temp[k] = src[k][i];
	//		for (unsigned k = 0; k < d; k++)
	//			src[k][i] = src[k][j];
	//		for (unsigned k = 0; k < d; k++)
	//			src[k][j] = temp[k];

	//	}

	//	void action(vec::vector &src) {
	//		PRECISION temp = src[i];
	//		src[i] = src[j];
	//		src[j] = temp;
	//	}
	//};

	//class HouseHolder {
	//private:
	//	unsigned len;
	//	PRECISION *data;
	//	PRECISION alpha;
	//public:
	//	HouseHolder() : len(1), alpha(0) { data = new PRECISION[len]; memset(data, 0, len * sizeof(PRECISION)); };
	//	HouseHolder(unsigned l) : len(l), alpha(0) { data = new PRECISION[len]; memset(data, 0, len * sizeof(PRECISION)); };
	//	HouseHolder(PRECISION *d, unsigned l, PRECISION a) : len(l), alpha(a) { data = new PRECISION[len]; memcpy(data, d, len * sizeof(PRECISION));};
	//	HouseHolder(PRECISION *d, unsigned l) : HouseHolder(d, l, -2.0 / vec::norm_pow(d, len, 2)) {};
	//	HouseHolder(const vec::vector &src) : HouseHolder(src.start(), src.length()) {};
	//	HouseHolder(const HouseHolder &src) : len(src.len), alpha(src.alpha) { data = new PRECISION[len]; memcpy(data, src.v(), len * sizeof(PRECISION)); };
	//	~HouseHolder() { delete[] data; };

	//	void swap(HouseHolder &rhs) {
	//		using std::swap;
	//		swap(this->len, rhs.len);
	//		swap(this->data, rhs.data);
	//		swap(this->alpha, rhs.alpha);
	//	}
	//	HouseHolder& operator=(HouseHolder & rhs) { swap(rhs); return (*this); }

	//	PRECISION* v() const { return data; };
	//	unsigned length() const { return len; };

	//	void action(PRECISION *vv) {
	//		PRECISION *temp = new PRECISION[len];
	//		PRECISION scalar = 0;
	//		memcpy(temp, vv, len * sizeof(PRECISION));
	//		scalar = alpha*vec::dot(data, vv, len);
	//		for (unsigned i = 0; i < len; i++)
	//			vv[i] += scalar*data[i];
	//	};
	//	void action(vec::vector &vv,unsigned i_start) {
	//		action(vv.start() + i_start);
	//	}
	//	void action(vec::vector &vv) { action(vv.start()); };
	//	void action(mat::matrix &M, unsigned i_start, unsigned j_start, unsigned M_LEN, char c) {
	//		PRECISION *scalars = new PRECISION[M_LEN];
	//		memset(scalars, 0, M_LEN * sizeof(PRECISION));
	//		if (c == 'L') {
	//			for (unsigned j = 0; j < M_LEN; j++) {
	//				for (unsigned i = 0; i < len; i++)
	//					scalars[j] += data[i] * M(i + i_start, j + j_start);
	//				scalars[j] *= alpha;
	//			}
	//		}
	//		else if (c == 'R') {
	//			for (unsigned i = 0; i < M_LEN; i++) {
	//				for (unsigned j = 0; j < len; j++)
	//					scalars[i] += M(i + i_start, j + j_start)*data[j];
	//				scalars[i] *= alpha;
	//			}
	//		}
	//		else {
	//				std::cout << "Error: the control variable in action() of HH is incorrect.\n";
	//				exit(1);
	//		}
	//
	//		if (c == 'L') {
	//			for (unsigned i = 0; i < len; i++) {
	//				for (unsigned j = 0; j < M_LEN; j++)
	//					M[i + i_start][j + j_start] = M[i + i_start][j + j_start] + scalars[j] * data[i];
	//			}
	//		}
	//		else if (c == 'R') {
	//			for (unsigned i = 0; i < M_LEN; i++) {
	//				for(unsigned j=0;j<len;j++)
	//					M[i + i_start][j + j_start] = M[i + i_start][j + j_start] + scalars[i] * data[j];
	//			}
	//		}
	//		else {
	//			std::cout << "Error: the control variable in action() of HH is incorrect.\n";
	//			exit(1);
	//		}
	//		delete[] scalars;
	//	}
	//	void action(mat::matrix &M, unsigned i_start, unsigned j_start, char c) {
	//		if (c == 'L')
	//			action(M, i_start, j_start, M.getCol() - j_start, c);
	//		else if (c == 'R')
	//			action(M, i_start, j_start, M.getRow() - i_start, c);
	//		else {
	//			std::cout << "Error: the control variable in action() of HH is incorrect.\n";
	//			exit(1);
	//		}
	//	};
	//};
}
