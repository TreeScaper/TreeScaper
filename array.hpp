#pragma once

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstring>

extern "C"
{
#include "f2c.h"
#include "clapack.h"
}
#undef abs
#undef max
#undef min

#define ARRAY_CONTAINER 100

#define PRECISION double

enum MAT_NORM_TYPE
{
	M_FOR_NORM
};

template <class T>
class Array
{
private:
	friend std::ostream &operator<<(std::ostream &output, const Array<T> &arr) // overload <<
	{
		output << "{" << arr.size << ": ";

		for (int i = 0; i < arr.size; i++)
			output << arr.vec[i] << ", ";
		output << "}";

		return output; // output form, for example {length: a1, a2, a3, }
	};

	int size;
	int container;
	int max_size;
	T *vec;

public:
	friend void swap(Array<T> &lhs, Array<T> &rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.container, rhs.container);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.vec, rhs.vec);
	}

	Array() : size(0), container(0), max_size(0), vec(nullptr){};

	Array(int n, int c)
		: size(n),
		  container(c),
		  max_size(size ? ((size - 1) / container + 1) * container : container),
		  vec(max_size ? new T[max_size] : nullptr){};

	Array(int n) : Array(n, n){};

	Array(int n, int c, const T *arr)
		: Array(n, c)
	{
		for (int i = 0; i < size; i++)
			vec[i] = arr[i];
	};

	Array(int n, const T *arr) : Array(n, n, arr){};

	Array(int n, int c, const T &ele) : Array(n, c)
	{
		for (int i = 0; i < size; i++)
			vec[i] = ele;
	}

	// Array(int n, const T& ele) : Array(n, n, ele) {};

	// Copy-constructor.
	Array(const Array<T> &arr) : Array(arr.size, arr.container)
	{
		for (int i = 0; i < size; i++)
			vec[i] = arr.vec[i];
	}

	~Array() { if (vec) delete[] vec; };

	friend Array<T> join(const Array<T> &lhs, const Array<T> &rhs)
	{
		T *data = new T[lhs.size + rhs.size];
		for (auto i = 0; i < lhs.size; i++)
			data[i] = lhs(i);
		for (auto i = 0; i < rhs.size; i++)
			data[i + lhs.size] = rhs(i);
		return Array<T>(lhs.size + rhs.size, data);
	}

	friend Array<T> duplicate(int n, const Array<T> &src, bool entrywise = false)
	{
		T *data = new T[n * src.size];
		if (entrywise)
			for (auto i = 0; i < src.size; i++)
				for (auto j = 0; j < n; j++)
					data[i * n + j] = src(i);
		else
			for (auto i = 0; i < n; i++)
				for (auto j = 0; j < src.size; j++)
					data[i * src.size + j] = src(j);
		return Array<T>(n * src.size, data);
	}

	void push(const T &ele)
	{
		if (size == max_size)
		{
			if (container == 0)
				container = ARRAY_CONTAINER;
			max_size = max_size ? (size / container + 1) * container : container;
			T *temp = new T[max_size];
			for (int i = 0; i < size; i++)
				temp[i] = vec[i];
			// memcpy(temp, vec, length * sizeof(T)); // Cause copy constructor and destructor issues.
			delete[] vec;
			vec = temp;
		}
		vec[size] = ele;
		size++;
	}

	void push(const T &ele, size_t pos)
	{
		assert(pos <= size);
		this->push(ele);
		if (pos == size - 1)
			return;
		else
		{
			// for (int i = pos; i < size - 1; i++)
			//	memcpy(vec + i + 1, vec + i, sizeof(T));
			memmove(vec + pos + 1, vec + pos, (size - pos - 1) * sizeof(T));
			vec[pos] = ele;
		}
	}

	int order_push(const T &ele, bool ISINCREASING)
	{
		int pos = this->binary_search(ele, ISINCREASING, true);
		this->push(ele, pos);
		return pos;
	}

	Array<T> &operator=(const Array<T> &rhs)
	{
		Array<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T &operator[](const int index) { return vec[index]; };

	T &ele(const int index) { return vec[index]; };

	T operator()(const int index) const { return vec[index]; };

	T &front() { return vec[0]; };

	T &back() { return vec[size - 1]; };

	T pop()
	{
		assert(size > 0);
		size--;
		return vec[size];
	};

	T pop(int index)
	{
		assert(index < size);
		auto val = vec[index];
		for (int i = index + 1; i < size; i++)
			vec[i - 1] = vec[i];
		size--;
		return val;
	}

	T *get_vec() const { return vec; };

	T *pop_vec()
	{
		T *temp = vec;
		vec = nullptr;
		size = 0;
		max_size = 0;
		return temp;
	}

	void erase() { size = 0; };

	void clean() { size = 0; };

	void release()
	{
		size = 0;
		max_size = 0;
		delete[] vec;
		vec = nullptr;
	}

	void set_ele(const T &ele)
	{
		for (int i = 0; i < size; i++)
			vec[i] = ele;
	}

	void set_ele(const T &ele, int start, int end)
	{
		assert((0 <= start) && (start <= end) && (end < size));
		for (int i = start; i < end + 1; i++)
			vec[i] = ele;
	}

	void set_ele(const T &ele, int i)
	{
		assert((0 <= i) && (i < size));
		vec[i] = ele;
	}

	void set_container(int c)
	{
		assert(container == 0);
		container = c;
	};

	void resize(int n)
	{
		if (n > max_size)
		{
			if (container == 0)
				container = n;
			max_size = ((n + 1) / container) * container;
			T *temp = new T[max_size];
			for (int i = 0; i < size; i++)
				temp[i] = vec[i];
			delete[] vec;
			vec = temp;
			size = n;
		}
		else
			size = n;
	};

	void resize_uinit(int n)
	{
		if (n > max_size)
		{
			if (container == 0)
				container = n;
			max_size = ((n + 1) / container) * container;
			T *temp = new T[max_size];
			delete[] vec;
			vec = temp;
			size = n;
		}
		else
			size = n;
	}

	void resize(int n, const T &ele)
	{
		if (n > max_size)
		{
			if (container == 0)
				container = n;
			max_size = ((n + 1) / container) * container;
			T *temp = new T[max_size];
			for (int i = 0; i < n; i++)
				temp[i] = ele;
			delete[] vec;
			vec = temp;
			size = n;
		}
		else
		{
			size = n;
			this->set_ele(ele);
		}
	};

	Array<T> &align()
	{
		assert(size > 0);
		Array<T> temp(size, vec);
		swap(*this, temp);
		return *this;
	}

	int get_size() const { return size; };

	bool is_empty() { return (size == 0); };

	void bubble_sort()
	{
		int n = size;
		bool swapped = false;

		while (n > 1)
		{
			int new_n = 0;
			int i = 1;
			for (int i = 1; i < n; i++)
			{
				if (vec[i - 1] > vec[i])
				{
					// Swapped
					using std::swap;
					swap(vec[i - 1], vec[i]);
					new_n = i;
				}
			}
			n = new_n;
		}
	}

	template <class S>
	void bubble_sort(Array<S> &arr_adj)
	{
		assert(size == arr_adj.get_size());
		int n = size;
		bool swapped = false;

		while (n > 1)
		{
			int new_n = 0;
			int i = 1;
			for (int i = 1; i < n; i++)
			{
				if (vec[i - 1] > vec[i])
				{
					// Swapped
					using std::swap;
					swap(vec[i - 1], vec[i]);
					swap(arr_adj[i - 1], arr_adj[i]);
					new_n = i;
				}
			}
			n = new_n;
		}
	}

	int binary_search(const T &src, bool ISINCREASING, bool GET_NEXT_IF_NOT_FOUND = false)
	{
		if (size == 0)
		{
			if (GET_NEXT_IF_NOT_FOUND)
				return 0;
			else
				return -1;
		}
		else if (size == 1)
		{
			if (vec[0] == src)
				return 0;
			else if ((vec[0] < src) ^ ISINCREASING)
			{
				if (GET_NEXT_IF_NOT_FOUND)
					return 0;
				else
					return -1;
			}
			else
			{
				if (GET_NEXT_IF_NOT_FOUND)
					return size;
				else
					return -1;
			}
		}

		int a = 0, b = size - 1;

		if ((vec[a] < src) ^ ISINCREASING)
		{
			if (GET_NEXT_IF_NOT_FOUND)
				return 0;
			else
				return -1;
		}
		else if (!((vec[b] < src) ^ ISINCREASING))
		{
			if (GET_NEXT_IF_NOT_FOUND)
				return size;
			else
				return -1;
		}
		if (vec[a] == src)
			return a;
		else if (vec[b] == src)
			return b;

		int len = b - a + 1;
		int ind = a + (len / 2);

		while (len > 2)
		{
			if (vec[ind] == src)
				return ind;
			else if ((vec[ind] < src) ^ ISINCREASING)
			{
				b = ind;
				len = b - a + 1;
				ind = a + (len / 2);
			}
			else
			{
				a = ind;
				len = b - a + 1;
				ind = a + (len / 2);
			}
		}

		if (GET_NEXT_IF_NOT_FOUND)
			return b;
		else
			return -1;
	}
};

template <class T>
class Array2D
{
	// This is a 2D array implementation. It contains Array< Array<T>*>
	// recording pointer to a container implemented as Array<T>.
private:
	Array<Array<T> *> ptr_container;
	size_t size;
	size_t max_size;

public:
	Array2D() : size(0), max_size(0), ptr_container(Array<Array<T> *>()){};
	Array2D(int c) : size(0), max_size(c), ptr_container(Array<Array<T> *>(c, c, (Array<T> *)nullptr)) { ptr_container.resize(0); };
	Array2D(int c, int e) : size(c), max_size(c), ptr_container(Array<Array<T> *>(c))
	{
		for (int i = 0; i < size; i++)
			ptr_container.set_ele(new Array<T>(0, e), i);
	};
	Array2D(int c, int e, const T &ele) : size(c), max_size(c), ptr_container(Array<Array<T> *>(c))
	{
		for (int i = 0; i < size; i++)
			ptr_container[i] = new Array<T>(e, e, ele);
	};
	Array2D(const Array2D &src) : Array2D(src.max_size)
	{
		for (int i = 0; i < max_size; i++)
			if (src.ptr_container(i) != nullptr)
				ptr_container[i] = new Array<T>(*src.ptr_container(i));
	}

	~Array2D()
	{
		for (int i = 0; i < max_size; i++)
			delete ptr_container[i];
	}

	friend void swap(Array2D<T> &lhs, Array2D<T> &rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.ptr_container, rhs.ptr_container);
	}

	Array2D<T> &operator=(const Array2D<T> &rhs)
	{
		Array2D<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	Array<T> &operator[](int i) { return *ptr_container[i]; };

	Array<T> &operator()(int i) const { return *ptr_container[i]; };

	Array<T> *get_c_ptr(int i) { return ptr_container[i]; };

	size_t get_max_size() { return max_size; };

	size_t get_size() { return size; };

	T &ele(int i, int j) { return ptr_container[i]->ele(j); };

	void set_max_size(int c)
	{
		assert(max_size <= c);
		max_size = c;
	}

	void set_c(int i, Array<T> *src_ptr)
	{
		assert((0 <= i) && (i < max_size));
		if (ptr_container[i] != nullptr)
		{
			delete ptr_container[i];
			ptr_container[i] = src_ptr;
		}
		else
		{
			ptr_container[i] = src_ptr;
			size++;
		}
	}

	void push_c(Array<T> *src_ptr)
	{
		assert(size < max_size);
		ptr_container.push(src_ptr);
		size++;
	}

	void push(const T &src, int i)
	{
		size_t ind = (i >= 0) ? i : (size + i);
		Array<T> *c_ptr = ptr_container[ind];
		if (c_ptr == nullptr)
		{
			c_ptr = new Array<T>();
			c_ptr->push(src);
			ptr_container[ind] = c_ptr;
		}
		else
			c_ptr->push(src);
	}

	void push(const T &src)
	{
		this->push(src, -1);
	}

	Array<T> *pop_c(int i)
	{
		assert((0 <= i) && (i < max_size));
		Array<T> *temp = (ptr_container[i] == nullptr) ? nullptr : ptr_container[i];
		if (temp != nullptr)
			size--;
		ptr_container[i] = nullptr;
		return temp;
		// if (ptr_container[i] == nullptr)
		//	return nullptr;
		// else
		//	return ptr_container[i];
	}

	void earse()
	{
		for (int i = 0; i < max_size; i++)
			if (ptr_container[i] != nullptr)
				ptr_container[i]->earse();
	}

	void earse(int i)
	{
		if (ptr_container[i] != nullptr)
			ptr_container[i]->earse();
	}

	void clean()
	{
		for (int i = 0; i < max_size; i++)
		{
			if (ptr_container[i] != nullptr)
			{
				delete ptr_container[i];
				ptr_container.set_ele(nullptr, i);
			}
		}
		size = 0;
	}

	void resize(int n)
	{
		ptr_container.resize(n);
		size = n;
		max_size = n;
	}
};

// template<class... Types>
// Variadic templates ?

template <class T1, class T2>
class Array2D_tuple
{
	// This one stores arrays of tuple (T1, T2), where T1 is the key type and T2 is the value type.
	typedef Array2D_tuple self_type;
	typedef T1 key_type;
	typedef T2 val_type;
	typedef Array2D<T1> key_arr_type;
	typedef Array2D<T2> val_arr_type;

private:
	key_arr_type key_;
	val_arr_type val_;
	size_t size;
	size_t max_size;

public:
	Array2D_tuple() : size(0), max_size(0), key_(key_arr_type()), val_(val_arr_type()){};
	Array2D_tuple(int c) : size(0), max_size(c), key_(key_arr_type(c)), val_(val_arr_type(c)){};
	Array2D_tuple(int c, int e) : size(c), max_size(c), key_(key_arr_type(c, e)), val_(val_arr_type(c, e)){};
	Array2D_tuple(int c, int e, const key_type &k, const val_type &v) : size(c),
																		max_size(c),
																		key_(key_arr_type(c, e, k)),
																		val_(val_arr_type(c, e, v)){};
	Array2D_tuple(const self_type &src) : size(src.size),
										  max_size(src.max_size),
										  key_(key_arr_type(src.key_)),
										  val_(val_arr_type(src.val_)){};
	~Array2D_tuple(){};

	friend void swap(Array2D_tuple<T1, T2> &lhs, Array2D_tuple<T1, T2> &rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.key_, rhs.key_);
		swap(lhs.val_, rhs.val_);
	}

	self_type &operator=(const self_type &rhs)
	{
		self_type temp(rhs);
		swap(*this, temp);
		return *this;
	}

	void init_c(int i, int e)
	{
		assert((0 <= i) && (i < max_size));
		key_.set_c(i, new Array<key_type>(0, e));
		val_.set_c(i, new Array<val_type>(0, e));
	}

	void init_c(int i, Array<key_type> *src_key_ptr, Array<val_type> *src_val_ptr)
	{
		assert((0 <= i) && (i < max_size));
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.set_c(i, src_key_ptr);
		val_.set_c(i, src_val_ptr);
	}

	key_arr_type &key() { return key_; };
	val_arr_type &val() { return val_; };

	Array<key_type> *get_key_ptr(int i) { return key_.get_c_ptr(i); };
	Array<val_type> *get_val_ptr(int i) { return val_.get_c_ptr(i); };

	key_type &ele_key(int i, int j) { return key_.ele(i, j); };
	val_type &ele_val(int i, int j) { return val_.ele(i, j); };

	void set_max_size(int c)
	{
		assert(max_size <= c);
		max_size = c;
	}

	size_t get_max_size() { return max_size; }

	size_t get_size() { return size; };

	void set_c(int i, Array<key_type> *src_key_ptr, Array<val_type> *src_val_ptr)
	{
		assert((0 <= i) && (i < max_size));
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.set_c(i, src_key_ptr);
		val_.set_c(i, src_val_ptr);
	}

	void push_c(Array<key_type> *src_key_ptr, Array<val_type> *src_val_ptr)
	{
		assert(key_.get_size() < max_size);
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.push_c(src_key_ptr);
		val_.push_c(src_val_ptr);
		size++;
	}

	Array<key_type> *pop_key_c(int i)
	{
		assert((0 <= i) && (i < max_size));
		if (val_.get_c_ptr(i) == nullptr)
			size--;
		return key_.pop_c(i);
	}

	Array<val_type> *pop_val_c(int i)
	{
		assert((0 <= i) && (i < max_size));
		if (key_.get_c_ptr(i) == nullptr)
			size--;
		return val_.pop_c(i);
	}

	void sort_by_key(int i) { key_[i].bubble_sort(val_[i]); };

	void sort_by_val(int i) { val_[i].bubble_sort(key_[i]); };

	void earse()
	{
		for (int i = 0; i < max_size; i++)
		{
			if (key_.get_c_ptr(i) != nullptr)
			{
				key_.get_c_ptr(i)->erase();
				val_.get_c_ptr(i)->erase();
			}
		}
	}

	void earse(int i)
	{
		if (key_.get_c_ptr(i) != nullptr)
		{
			key_.get_c_ptr(i)->erase();
			val_.get_c_ptr(i)->erase();
		}
	}

	void resize(size_t n)
	{
		key_.resize(n);
		val_.resize(n);
		size = n;
		max_size = n;
	}
};

template <class T>
class Matrix
{
	typedef Matrix<T> self_type;

private:
	size_t row;
	size_t col;
	T **vec;

public:
	Matrix<T>(size_t r, size_t c) : row(r), col(c), vec((r > 0) ? new T *[row] : nullptr)
	{
		if ((row > 0) && (col > 0))
		{
			for (int i = 0; i < row; i++)
				vec[i] = new T[col];
		}
	};
	Matrix<T>() : Matrix(0, 0){};
	Matrix<T>(size_t r, size_t c, const T &ele) : Matrix(r, c)
	{
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vec[i][j] = ele;
	};
	Matrix<T>(int r, int c, T **src) : row(r), col(c), vec(src){};
	Matrix<T>(const Matrix &src) : Matrix(src.row, src.col)
	{
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vec[i][j] = src.vec[i][j];
	}
	~Matrix<T>()
	{
		for (int i = 0; i < row; i++)
			delete[] vec[i];
		delete vec;
	}

	friend void swap(self_type &lhs, self_type &rhs)
	{
		using std::swap;
		swap(lhs.row, rhs.row);
		swap(lhs.col, rhs.col);
		swap(lhs.vec, rhs.vec);
	}

	self_type &operator=(const self_type &rhs)
	{
		self_type temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T *operator[](int i) { return vec[i]; };

	T &operator()(int i, int j) { return vec[i][j]; };

	T ele(int i, int j) const { return vec[i][j]; };

	size_t get_row() const { return row; };

	size_t get_col() const { return col; };

	T **get_vec() { return vec; };

	void set_row(int i, T *src)
	{
		assert((0 <= i) && (i < row));
		delete vec[i];
		vec[i] = src;
	}

	void resize(size_t r, size_t c)
	{
		if (r != row || c != col)
		{
			if (vec != nullptr)
			{
				for (int i = 0; i < col; i++)
					delete[] vec[i];
				delete[] vec;
			}
			row = r;
			col = c;
			vec = (r > 0) ? new T *[row] : nullptr;
			if ((row > 0) && (col > 0))
			{
				for (int i = 0; i < row; i++)
					vec[i] = new T[col];
			}
		}
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vec[i][j] = 0;
	}

	void print(std::ostream &output) const
	{
		if (row == 0 || col == 0)
			return;
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
				output << vec[i][j] << '\t';
			output << '\n';
		}
	}

	friend void shift(Matrix<T> &mat, const T &s)
	{
		for (auto i = 0; i < mat.row; i++)
			for (auto j = 0; j < mat.col; j++)
				mat.vec[i][j] += s;
	}

	friend void shift(Matrix<T> &mat, const Matrix<T> &displace, const PRECISION scalar)
	{
		for (auto i = 0; i < mat.row; i++)
			for (auto j = 0; j < mat.col; j++)
				mat.vec[i][j] += scalar * displace.vec[i][j];
	}

	friend void shift(Matrix<T> &mat, const Matrix<T> &displace)
	{
		for (auto i = 0; i < mat.row; i++)
			for (auto j = 0; j < mat.col; j++)
				mat.vec[i][j] += displace.vec[i][j];
	}

	friend void shift(Matrix<T> &mat, const Array<T> &displace, char Broadcast, const PRECISION scalar)
	{
		if (Broadcast == 'r' || Broadcast == 'R')
			for (auto i = 0; i < mat.row; i++)
				for (auto j = 0; j < mat.col; j++)
					mat.vec[i][j] += scalar * displace[j];
		else if (Broadcast == 'c' || Broadcast == 'C')
			for (auto i = 0; i < mat.row; i++)
				for (auto j = 0; j < mat.col; j++)
					mat.vec[i][j] += scalar * displace[i];
		else
		{
			std::cout << "Error! Broadcast pattern not found! Use `r' / `R' for Broadcast_By_Row or `c' / `C' for Broadcast_By_Column.\n";
			throw(1);
		}
	}

	friend void shift(Matrix<T> &mat, const Array<T> &displace, char Broadcast)
	{
		if (Broadcast == 'r' || Broadcast == 'R')
			for (auto i = 0; i < mat.row; i++)
				for (auto j = 0; j < mat.col; j++)
					mat.vec[i][j] += displace(j);
		else if (Broadcast == 'c' || Broadcast == 'C')
			for (auto i = 0; i < mat.row; i++)
				for (auto j = 0; j < mat.col; j++)
					mat.vec[i][j] += displace(i);
		else
		{
			std::cout << "Error! Broadcast pattern not found! Use `r' / `R' for Broadcast_By_Row or `c' / `C' for Broadcast_By_Column.\n";
			throw(1);
		}
	}

	friend void row_sum(Matrix<T> &mat, Array<T> &r_s)
	{
		for (auto j = 0; j < mat.col; j++)
			r_s[j] = 0;
		for (auto i = 0; i < mat.row; i++)
			for (auto j = 0; j < mat.col; j++)
				r_s[j] += mat.vec[i][j];
	}

	friend void row_mean(Matrix<T> &mat, Array<T> &r_m)
	{
		row_sum(mat, r_m);
		for (auto j = 0; j < mat.col; j++)
			r_m[j] /= mat.row;
	}

	friend void centralize(Matrix<T> &mat)
	{
		Array<T> c(mat.col);
		row_mean(mat, c);
		shift(mat, c, 'R');
	}

	// friend bool inverse(Matrix<T>& mat) { return true; };

	friend bool inverse(Matrix<T> &mat)
	{
		assert(mat.row == mat.col);
		unsigned row = mat.row, col = mat.col;

		real *V = new real[row * col];
		real *Vn = new real[row * col];
		integer *ipiv = new integer[row * col];
		integer M = row;
		integer lda = row;
		integer lwork = row;
		integer INFO;

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				V[i * col + j] = mat.vec[i][j];

		sgetrf_(&M, &M, V, &M, ipiv, &INFO);
		if (INFO != 0)
		{
			std::cout << "warning: step1 fail.\n";
			delete[] V;
			delete[] Vn;
			delete[] ipiv;
			return false;
		}

		sgetri_(&M, V, &lda, ipiv, Vn, &lwork, &INFO);

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				mat.vec[i][j] = V[i * col + j];

		delete[] V;
		delete[] Vn;
		delete[] ipiv;
		if (INFO == 0)
		{
			return true;
		}
		else
		{
			std::cout << "warning: step2 fail.\n";
			return false;
		}
	};

	PRECISION norm2(MAT_NORM_TYPE mnt)
	{
		PRECISION ans = 0;
		switch (mnt)
		{
		case M_FOR_NORM:
			for (auto i = 0; i < row; i++)
				for (auto j = 0; j < col; j++)
					ans += vec[i][j] * vec[i][j];
			return ans;
		default:
			std::cout << "Error! Matrix norm type not found!\n";
			throw(1);
		}
	};

	PRECISION norm(MAT_NORM_TYPE mnt)
	{
		PRECISION ans = 0;
		switch (mnt)
		{
		case M_FOR_NORM:
			return sqrt(this->norm2(mnt));
		default:
			std::cout << "Error! Matrix norm type not found!\n";
			throw(1);
		}
	};
};