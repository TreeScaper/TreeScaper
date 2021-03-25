#pragma once

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cmath>

#define ARRAY_CONTAINER 100

enum MatrixOutputType { FULL, LOWERTRI };


template<class T>
class Array {
private:
	friend std::ostream& operator<<(std::ostream& output, const Array<T>& arr)           // overload <<
	{
		output << "{" << arr.size << ": ";

		for (int i = 0; i < arr.size; i++)
			output << arr.vec[i] << ", ";
		output << "}";

		return output;                                                         // output form, for example {length: a1, a2, a3, }
	};

	int size;
	int container;
	int max_size;
	T* vec;
public:

	friend void swap(Array<T>& lhs, Array<T>& rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.container, rhs.container);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.vec, rhs.vec);
	}

	Array() : size(0), container(0), max_size(0), vec(nullptr) {};

	Array(int n, int c)
		: size(n),
		container(c),
		max_size(size ? ((size - 1) / container + 1) * container : container),
		vec(max_size ? new T[max_size] : nullptr) {};

	Array(int n) : Array(n, n) {};
	
	Array(int n, int c, const T* arr)
		: Array(n, c)
	{
		for (int i = 0; i < size; i++)
			vec[i] = arr[i];
	};


	Array(int n, const T* arr) : Array(n, n, arr) {};

	Array(int n, int c, const T& ele) : Array(n, c){
		for (int i = 0; i < size; i++)
			vec[i] = ele;
	}

	// Array(int n, const T& ele) : Array(n, n, ele) {};

	//Copy-constructor.
	Array(const Array<T>& arr) : Array(arr.size, arr.container){
		for (int i = 0; i < size; i++)
			vec[i] = arr.vec[i];
	}

	~Array() { delete[] vec; };


	void push(const T& ele) {
		if (size == max_size) {
			if (container == 0)
				container = ARRAY_CONTAINER;
			max_size = max_size ? (size / container + 1) * container : container;
			T* temp = new T[max_size];
			for (int i = 0; i < size; i++) 
				temp[i] = vec[i];
			//memcpy(temp, vec, length * sizeof(T)); // Cause copy constructor and destructor issues.
			delete[] vec;
			vec = temp;
		}
		vec[size] = ele;
		size++;
	}

	Array<T>& operator=(const Array<T>& rhs)
	{
		Array<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T& operator[](const int index) { return vec[index]; };

	T& ele(const int index) { return vec[index]; };

	T operator()(const int index) const { return vec[index]; };

	T& front() { return vec[0]; };

	T& back() { return vec[size - 1]; };

	T pop() { return vec[size - 1]; size--; };

	T* get_vec() const { return vec; };


	void erase() { size = 0; };

	void clean() { size = 0; };

	void set_ele(const T& ele) {
		for (int i = 0; i < size; i++)
			vec[i] = ele;
	}

	void set_ele(const T& ele, int start, int end) {
		assert((0 <= start) && (start <= end) && (end < size));
		for (int i = start; i < end + 1; i++)
			vec[i] = ele;
	}

	void set_ele(const T& ele, int i) {
		assert((0 <= i) && (i < size));
		vec[i] = ele;
	}

	void set_container(int c) { assert(container == 0); container = c; };

	void resize(int n) { 
		if (n > max_size) {
			max_size = (n / container + 1) * container;
			T* temp = new T[max_size];
			for (int i = 0; i < size; i++)
				temp[i] = vec[i];
			delete[] vec;
			vec = temp;
			size = n;
		}
		else
			size = n;
	};

	void resize_uinit(int n) {
		if (n > max_size) {
			max_size = (n / container + 1) * container;
			T* temp = new T[max_size];
			delete[] vec;
			vec = temp;
			size = n;
		}
		else
			size = n;
	}

	void resize(int n, const T& ele) {
		if (n > max_size) {
			max_size = (n / container + 1) * container;
			T* temp = new T[max_size];
			for (int i = 0; i < n; i++)
				temp[i] = ele;
			delete[] vec;
			vec = temp;
			size = n;
		}
		else {
			size = n;
			this->set_ele(ele);
		}
	};

	//void align() { resize(this->length, this->length); };

	int get_size() const { return size; };

	bool is_empty() { return (size == 0); };

	void bubble_sort() {
		int n = size;
		bool swapped = false;

		while (n > 1) {
			int new_n = 0;
			int i = 1;
			for (int i = 1; i < n; i++) {
				if (vec[i - 1] > vec[i]) {
					// Swapped
					using std::swap;
					swap(vec[i - 1], vec[i]);
					new_n = i;
				}
			}
			n = new_n;
		}
	}

	template<class S>
	void bubble_sort(Array<S>& arr_adj) {
		assert(size == arr_adj.get_size());
		int n = size;
		bool swapped = false;

		while (n > 1) {
			int new_n = 0;
			int i = 1;
			for (int i = 1; i < n; i++) {
				if (vec[i - 1] > vec[i]) {
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
};

template<class T>
class Array2D {
	// This is a 2D array implementation. It contains Array< Array<T>*> 
	// recording pointer to a container implemented as Array<T>.
private:
	Array<Array<T>*> ptr_container;
	size_t size;
	size_t max_size;

public:
	
	Array2D() : size(0), max_size(0), ptr_container(Array<Array<T>* >()) {};
	Array2D(int c) : size(0), max_size(c), ptr_container(Array<Array<T>* >(c, c, (Array<T>*)nullptr)) { ptr_container.resize(0); };
	Array2D(int c, int e) : size(c), max_size(c), ptr_container(Array<Array<T>* >(c)) {
		for (int i = 0; i < size; i++)
			ptr_container.set_ele(new Array<T>(0, e), i);
	};
	Array2D(int c, int e, const T& ele) : size(c), max_size(c), ptr_container(Array<Array<T>* >(c)) {
		for (int i = 0; i < size; i++)
			ptr_container[i] = new Array<T>(e, e, ele);
	};
	Array2D(const Array2D& src) : Array2D(src.max_size) {
		for (int i = 0; i < max_size; i++)
			if (src.ptr_container(i) != nullptr)
				ptr_container[i] = new Array<T>(*src.ptr_container(i));
	}

	~Array2D() {
		for (int i = 0; i < max_size; i++)
			delete ptr_container[i];
	}

	friend void swap(Array2D<T>& lhs, Array2D<T>& rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.ptr_container, rhs.ptr_container);
	}

	Array2D<T>& operator=(const Array2D<T>& rhs)
	{
		Array2D<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	Array<T>& operator[](int i) { return *ptr_container[i]; };

	Array<T>& operator()(int i) const { return *ptr_container[i]; };

	Array<T>* get_c_ptr(int i) { return ptr_container[i]; };

	size_t get_max_size() { return max_size; };

	size_t get_size() { return size; };

	T& ele(int i, int j) { return ptr_container[i]->ele(j); };

	void set_max_size(int c) { assert(max_size <= c); max_size = c; }

	void set_c(int i, Array<T>* src_ptr) {
		assert((0 <= i) && (i < max_size));
		if (ptr_container[i] != nullptr) {
			delete ptr_container[i];
			ptr_container[i] = src_ptr;
		}
		else {
			ptr_container[i] = src_ptr;
			size++;
		}
		
	}

	void push_c(Array<T>* src_ptr) {
		assert(size < max_size);
		ptr_container.push(src_ptr);
		size++;
	}

	void push(const T& src, int i) {
		size_t ind = (i >= 0) ? i : (size + i);
		Array<T>* c_ptr = ptr_container[ind];
		if (c_ptr == nullptr) {
			c_ptr = new Array<T>();
			c_ptr->push(src);
			ptr_container[ind] = c_ptr;
		}
		else
			c_ptr->push(src);
	}

	void push(const T& src) {
		this->push(src, -1);
	}

	Array<T>* pop_c(int i) {
		assert((0 <= i) && (i < max_size));
		Array<T>* temp = (ptr_container[i] == nullptr) ? nullptr : ptr_container[i];
		if (temp != nullptr)
			size--;
		ptr_container[i] = nullptr;
		return temp;
		//if (ptr_container[i] == nullptr)
		//	return nullptr;
		//else
		//	return ptr_container[i];
	}

	void earse() {
		for (int i = 0; i < max_size; i++)
			if (ptr_container[i] != nullptr)
				ptr_container[i]->earse();
	}

	void earse(int i) {
		if (ptr_container[i] != nullptr) 
			ptr_container[i]->earse();
	}

	void clean() {
		for (int i = 0; i < max_size; i++) {
			if (ptr_container[i] != nullptr) {
				delete ptr_container[i];
				ptr_container.set_ele(nullptr, i);
			}
		}
		size = 0;
	}

	void resize(int n) {
		ptr_container.resize(n);
		size = n;
		max_size = n;
	}
};

// template<class... Types>
// Variadic templates ?

template<class T1, class T2>
class Array2D_tuple {
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
	Array2D_tuple() : size(0), max_size(0), key_(key_arr_type()), val_(val_arr_type()) {};
	Array2D_tuple(int c) : size(0), max_size(c), key_(key_arr_type(c)), val_(val_arr_type(c)) {};
	Array2D_tuple(int c, int e) : size(c), max_size(c), key_(key_arr_type(c, e)), val_(val_arr_type(c, e)) {};
	Array2D_tuple (int c, int e, const key_type& k, const val_type& v) : 
		size(c),
		max_size(c),
		key_(key_arr_type(c, e, k)), 
		val_(val_arr_type(c, e, v)) {};
	Array2D_tuple(const self_type& src) : 
		size(src.size),
		max_size(src.max_size),
		key_(key_arr_type(src.key_)), 
		val_(val_arr_type(src.val_)) {};
	~Array2D_tuple() {};

	friend void swap(Array2D_tuple<T1, T2>& lhs, Array2D_tuple<T1, T2>& rhs)
	{
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.key_, rhs.key_);
		swap(lhs.val_, rhs.val_);
	}

	self_type& operator=(const self_type& rhs)
	{
		self_type temp(rhs);
		swap(*this, temp);
		return *this;
	}

	void init_c(int i, int e) {
		assert((0 <= i) && (i < max_size));
		key_.set_c(i, new Array<key_type>(0, e));
		val_.set_c(i, new Array<val_type>(0, e));
	}

	void init_c(int i, Array<key_type>* src_key_ptr, Array<val_type>* src_val_ptr) {
		assert((0 <= i) && (i < max_size));
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.set_c(i, src_key_ptr);
		val_.set_c(i, src_val_ptr);
	}

	key_arr_type& key() { return key_; };
	val_arr_type& val() { return val_; };

	Array<key_type>* get_key_ptr(int i) { return key_.get_c_ptr(i); };
	Array<val_type>* get_val_ptr(int i) { return val_.get_c_ptr(i); };

	key_type& ele_key(int i, int j) { return key_.ele(i, j); };
	val_type& ele_val(int i, int j) { return val_.ele(i, j); };

	void set_max_size(int c) { assert(max_size <= c); max_size = c; }

	size_t get_max_size() { return max_size; }

	size_t get_size() { return size; };

	void set_c(int i, Array<key_type>* src_key_ptr, Array<val_type>* src_val_ptr) {
		assert((0 <= i) && (i < max_size));
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.set_c(i, src_key_ptr);
		val_.set_c(i, src_val_ptr);
	}

	void push_c(Array<key_type>* src_key_ptr, Array<val_type>* src_val_ptr) {
		assert(key_.get_size() < max_size);
		assert(src_key_ptr->get_size() == src_val_ptr->get_size());
		key_.push_c(src_key_ptr);
		val_.push_c(src_val_ptr);
		size++;
	}

	Array<key_type>* pop_key_c(int i) {
		assert((0 <= i) && (i < max_size));
		if (val_[i] == nullptr)
			size--;
		return key_.pop_c(i);
	}

	Array<val_type>* pop_val_c(int i) {
		assert((0 <= i) && (i < max_size));
		if (key_[i] == nullptr)
			size--;
		return val_.pop_c(i);
	}

	void sort_by_key(int i) { key_[i].bubble_sort(val_[i]); };

	void sort_by_val(int i) { val_[i].bubble_sort(key_[i]); };
	
	void earse() {
		for (int i = 0; i < max_size; i++) {
			if (key_.get_c_ptr(i) != nullptr) {
				key_.get_c_ptr(i)->erase();
				val_.get_c_ptr(i)->erase();
			}
		}
	}

	void earse(int i) {
		if (key_.get_c_ptr(i) != nullptr) {
			key_.get_c_ptr(i)->erase();
			val_.get_c_ptr(i)->erase();
		}
	}
	
	void resize(size_t n) {
		key_.resize(n);
		val_.resize(n);
		size = n;
		max_size = n;
	}
};

template<class T>
class Matrix {
	typedef Matrix<T> self_type;
private:
	size_t row;
	size_t col;
	T** vec;

public:
	Matrix(size_t r, size_t c) : row(r), col(c), vec((row > 0) ? new T*[row] : nullptr) {
		if ((row > 0) && (col > 0)) {
			for (int i = 0; i < row; i++)
				vec[i] = new T[col];
		}
	};
	Matrix() : Matrix(0, 0) {};
	Matrix(size_t r, size_t c, const T& ele) : Matrix(r, c) {
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vec[i][j] = ele;
	};
	Matrix(int r, int c, T** src) : row(r), col(c), vec(src) {};
	Matrix(const Matrix& src) : Matrix(src.row, src.col) {
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vec[i][j] = src.vec[i][j];
	}
	~Matrix() {
		for (int i = 0; i < row; i++)
			delete[] vec[i];
		delete vec;
	}

	friend void swap(self_type& lhs, self_type& rhs)
	{
		using std::swap;
		swap(lhs.row, rhs.row);
		swap(lhs.col, rhs.col);
		swap(lhs.vec, rhs.vec);
	}

	self_type& operator=(const self_type& rhs)
	{
		self_type temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T* operator[](int i) { return vec[i]; };

	T& operator()(int i, int j) { return vec[i][j]; };

	size_t get_row() { return row; };

	size_t get_col() { return col; };

	size_t get_vec() { return vec; };

	void set_row(int i, T* src) {
		assert((0 <= i) && (i < row));
		delete vec[i];
		vec[i] = src;
	}

	void print(std::ostream &output, MatrixOutputType motype){
		if (motype == FULL){
			for (int i = 0; i < row; i++){
				for (int j = 0; j < col - 1; j++)
					output << vec[i][j] << '\t';
				output << vec[i][col - 1] << '\n';
			}
		}
		else{
			for (int i = 0; i < row; i++){
				for (int j = 0; j < i; j++)
					output << vec[i][j] << '\t';
				output << vec[i][i] << '\n';
			}
		}
	}
};