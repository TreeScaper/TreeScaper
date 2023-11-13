#pragma once

#include <iomanip>
#include <iostream>
#include <cassert>
#define ADDRESS_CONTAINER 10
#define DYNAMIC_ARRAY_CONTAINER 100

#ifndef INTEGER_TYPE
#define INTEGER_TYPE int
#endif

class address_list {
public:
	address_list() : vec(nullptr), size(0), max_size(0) {};
	address_list(INTEGER_TYPE n) : size(0), max_size(n) { 
		vec = new void* [n]; 
		for (int i = 0; i < n; i++)
			vec[i] = nullptr;
	};
	address_list(const address_list& src) : size(src.size), max_size(src.max_size) {
		vec = new void* [max_size];
		for (int i = 0; i < max_size; i++)
			vec[i] = src.vec[i];
	};

	~address_list() {
		// assert(this->is_empty());
		// Only release the address list, the actual vectors recorded are not released.
		// To release vectors in record, use release_vector().

		assert(this->is_empty());
   		delete[] vec;
	}

	bool is_empty() {
		for (int i = 0; i < max_size; i++) {
			if (vec[i] != nullptr)
				return false;
		}
		return true;
	}


	template<class S>
	void release_vector(S& junk) {
		for (int i = 0; i < max_size; i++) {
			if (vec[i] != nullptr) {
				S* temp = (S*)vec[i];
				delete[] temp;
				vec[i] = nullptr;
			}
		}
	}

	template<class S>
	void release_vector(INTEGER_TYPE j, S& junk) {
		if (vec[j] != nullptr) {
			S* temp = (S*)vec[j];
			delete[] temp;
			vec[j] = nullptr;
		}
	}

	friend void swap(address_list& lhs, address_list& rhs){
		using std::swap;
		swap(lhs.size, rhs.size);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.vec, rhs.vec);
	}

	template<class T>
	void push(T* ele) {
		assert(size <= max_size);
		if (size == max_size) {
			if (size == 0)
				max_size = ADDRESS_CONTAINER;
			else
				max_size = (size / ADDRESS_CONTAINER + 1) * ADDRESS_CONTAINER;
			void** temp = new void*[max_size];
			for (int i = 0; i < size; i++)
				temp[i] = vec[i];
			for (int i = size; i < max_size; i++)
				temp[i] = nullptr;
			//memcpy(temp, vec, length * sizeof(T)); // Cause copy constructor and destructor issues.
			delete[] vec;
			vec = temp;
		}
		vec[size] = (void*) ele;
		size++;
	}

	void* pop() {
		assert(size > 0);
		void* ptr = vec[size - 1];
		vec[size - 1] = nullptr;
		size--;
		return ptr;
	}

	void* pop(INTEGER_TYPE i) {
		assert((i >= 0) || (i < size));
		if (i == (size - 1))
			return pop();
		else {
			void* ptr = vec[i];
			for (int j = i + 1; j < size; j++)
				vec[j - 1] = vec[j];
			vec[size - 1] = nullptr;
			size--;
			return ptr;
		}
	}

	INTEGER_TYPE get_size() const { return size; };

	INTEGER_TYPE get_max_size() const { return max_size; };

	void* operator[](INTEGER_TYPE i) { return vec[i]; };

	void* ele(INTEGER_TYPE i) { return vec[i]; };

	template<class T>
	void set_ele(INTEGER_TYPE i, T* ptr) { vec[i] = (void*)ptr; };

	void set_ele(INTEGER_TYPE i, void* ptr) { vec[i] = ptr; };

	void set_size(INTEGER_TYPE i) { assert(i <= max_size); size = i; }

	void* back() { return vec[size - 1]; };

private:
	void** vec;
	INTEGER_TYPE size;
	INTEGER_TYPE max_size;
};

template <class T>
class darray {
public:
	darray() : container(DYNAMIC_ARRAY_CONTAINER), size(0), address(new address_list), end_ptr(nullptr) {};
	darray(INTEGER_TYPE c) : container(c), size(0), address(new address_list()), end_ptr(nullptr) {};
	darray(INTEGER_TYPE n, T* src_ptr) : container(n), size(n), address(new address_list(1)), end_ptr(src_ptr + (n - 1)) {
		// Construction from shallow copy of given array.
		address->set_ele(0, src_ptr);
		src_ptr = nullptr; 
		// Avoid releasing array outside darray class, drop the reference kept in src_ptr.
	};
	darray(const darray<T>& src) : container(src.container), size(src.size) {
		if (size == 0) {
			address = new address_list;
			end_ptr = nullptr;
		}
		else{
			INTEGER_TYPE q = (size - 1) / container + 1; // All containers
			INTEGER_TYPE r = size - (q - 1) * container; // Vaild entries number of the last container
			address = new address_list(q);
			T* dest_vec = nullptr;
			T* src_vec = nullptr;
			for (INTEGER_TYPE i = 0; i < q; i++) {
				 dest_vec = new T[container];
				 src_vec = src.get_vec_const(i);
				 for (INTEGER_TYPE j = 0; j < container; j++)
					 dest_vec[j] = src_vec[j];
				 address->set_ele(i, dest_vec);
			}
			end_ptr = dest_vec + (r - 1);
		}
	}
	darray(INTEGER_TYPE c, INTEGER_TYPE s, const address_list& adds, T* endptr = nullptr) : container(c), size(s), address(adds), end_ptr(endptr) {
		// This is a construction by shallow copy. The argument is usually sent from a dynamic array that is about to be released.
		if (size > 0 && end_ptr == nullptr) {
			INTEGER_TYPE q = (size - 1) / container;	// Container indices
			INTEGER_TYPE r = size - q * container;	// Vaild entries number of the last container
			end_ptr = (T*)address->ele(q);
			end_ptr += (r - 1);
		}
	}
	~darray() {
		//std::cout << "Releasing darray with container\t" << container << ", size\t" << size << ".\n";
		this->release();

		
		delete address;
	}

	void shal_cpy(darray<T>& dest) {
		// This is a shallow copy that kept current array and overwrite dest.
		// shal_cpy must be followed with a deconstructor of the current array.
		darray<T> temp(container, size, address, end_ptr);
		swap(dest, temp);
	}

	class iterator {
	public:
		typedef typename darray<T>::iterator self_type;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;

		darray<T>::iterator() : cur_block(-1), cur_index(-1), container(0), address(nullptr), ptr_(nullptr) {};
		darray<T>::iterator(INTEGER_TYPE b, INTEGER_TYPE i, INTEGER_TYPE c, address_list* a, pointer ptr) : cur_block(b), cur_index(i), container(c), address(a), ptr_(ptr) {};
		//darray<T>::iterator(self_type&& src) : iterator() {
		//	swap(*this, src);
		//	std::cout << "this iterator:\t" << cur_block << '\t' << cur_index << '\n';
		//	std::cout << "src iterator:\t" << src.cur_block << '\t' << src.cur_index << '\n';
		//};
		darray<T>::iterator(const self_type& src) : cur_block(src.cur_block), cur_index(src.cur_index), container(src.container), address(src.address), ptr_(src.ptr_) {};

		//darray<T>::iterator(INTEGER_TYPE b, INTEGER_TYPE i, INTEGER_TYPE c, address_list* a, pointer ptr) {
		//	cur_block = b;
		//	cur_index = i;
		//	container = c;
		//	address = a;
		//	ptr_ = ptr;
		//};
		//darray<T>::iterator(const self_type& src) {
		//	cur_block = src.cur_block;
		//	cur_index = src.cur_index;
		//	container = src.container;
		//	address = src.address;
		//	ptr_ = src.ptr_;
		//};

		friend void swap(self_type& lhs, self_type& rhs) {
			using std::swap;
			swap(lhs.cur_block, rhs.cur_block);
			swap(lhs.cur_index, rhs.cur_index);
			swap(lhs.container, rhs.container);
			swap(lhs.address, rhs.address);
			swap(lhs.ptr_, rhs.ptr_);
		}

		self_type& operator=(const self_type& rhs) {
			self_type temp = self_type(rhs);
			swap(*this, temp);
			return *this;
		}


		self_type operator++(int) {
			if (cur_index < (container - 1)) { cur_index++; ptr_++; }
			else { 
				cur_index = 0;
				ptr_ = (cur_block++ == (address->get_size() - 1)) ?  nullptr : (pointer)address->ele(cur_block); }
			return *this;
		}
		self_type operator++() {
			self_type i = *this;
			if (cur_index < (container - 1)) { cur_index++; ptr_++; }
			else { 
				cur_index = 0;
				(cur_block == (address->get_size() - 1)) ? ptr_ = nullptr : ptr_ = (pointer)address->ele(cur_block++); }
			return i;
		}
		self_type operator--(int) {
			if (cur_index != 0) { cur_index--; ptr_--; }
			else {
				ptr_ = (cur_block-- == 0) ? nullptr : ((pointer)address->ele(cur_block) + container);
				(ptr_ == nullptr) ? cur_index = -1 : cur_index = container - 1;
			}
			return *this;
		}
		self_type operator--() {
			self_type i = *this;
			if (cur_index != 0) { cur_index--; ptr_--; }
			else {
				(--cur_block == 0) ? ptr_ = nullptr : ptr_ = ((pointer)address->ele(cur_block) + container);
				(ptr_ == nullptr) ? cur_index = -1 : cur_index = container - 1;
			}
			return i;
		}

		self_type operator+(const int& rhs) {
			self_type it_new(*this);
			if (rhs > 0) {
				for (int i = 0; i < rhs; i++)
					it_new++;
			}
			else if (rhs < 0) {
				for (int i = 0; i < -rhs; i++)
					it_new--;
			}
			return it_new;
		}
		self_type operator-(const int& rhs) {
			return (*this + (-rhs));
		}

		reference operator*() { return *ptr_; }
		pointer operator->() { return ptr_; }
		bool operator==(const self_type& rhs) { return ((cur_block == rhs.cur_block) && (cur_index == rhs.cur_index)); }
		bool operator!=(const self_type& rhs) { return ((cur_block != rhs.cur_block) || (cur_index != rhs.cur_index)); }
		self_type& get_ref() { return *this; }
	private:
		INTEGER_TYPE cur_block;
		INTEGER_TYPE cur_index;
		INTEGER_TYPE container;
		address_list* address;
		pointer ptr_;
	};

	class const_iterator {
	public:
		typedef typename darray<T>::const_iterator self_type;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;

		const_iterator() : cur_block(-1), cur_index(-1), container(0), address(nullptr), ptr_(nullptr) {};
		const_iterator(INTEGER_TYPE b, INTEGER_TYPE i, INTEGER_TYPE c, address_list* a, pointer ptr) : cur_block(b), cur_index(i), container(c), address(a), ptr_(ptr) {};
		const_iterator(const self_type& src) : cur_block(src.cur_block), cur_index(src.cur_index), container(src.container), address(src.address), ptr_(src.ptr_) {};

		friend void swap(self_type& lhs, self_type& rhs) {
			using std::swap;
			swap(lhs.cur_block, rhs.cur_block);
			swap(lhs.cur_index, rhs.cur_index);
			swap(lhs.container, rhs.container);
			swap(lhs.address, rhs.address);
			swap(lhs.ptr_, rhs.ptr_);
		}

		self_type& operator=(const self_type& rhs) {
			self_type temp = self_type(rhs);
			swap(*this, temp);
			return *this;
		}

		self_type operator++(int) {
			if (cur_index < (container - 1)) { cur_index++; ptr_++; }
			else { cur_index = 0; (cur_block == (address->get_size() - 1)) ? ptr_ = nullptr : ptr_ = address->ele(cur_block++); }
			return *this;
		}
		self_type operator++() {
			self_type i = *this;
			if (cur_index < (container - 1)) { cur_index++; ptr_++; }
			else { cur_index = 0; (cur_block == (address->get_size() - 1)) ? ptr_ = nullptr : ptr_ = address->ele(cur_block++); }
			return i;
		}
		self_type operator--(int) {
			if (cur_index != 0) { cur_index--; ptr_--; }
			else {
				(--cur_block == 0) ? ptr_ = nullptr : ptr_ = ((pointer)address->ele(cur_block) + container);
				(ptr_ == nullptr) ? cur_index = -1 : cur_index = container - 1;
			}
			return *this;
		}
		self_type operator--() {
			self_type i = *this;
			if (cur_index != 0) { cur_index--; ptr_--; }
			else { 
				(--cur_block == 0) ? ptr_ = nullptr : ptr_ = ((pointer)address->ele(cur_block) + container); 
				(ptr_ == nullptr) ? cur_index = -1 : cur_index = container - 1;
			}
			return i;
		}

		self_type operator+(const int& rhs) {
			self_type it_new(*this);
			if (rhs > 0) {
				for (int i = 0; i < rhs; i++)
					it_new++;
			}
			else if (rhs < 0) {
				for (int i = 0; i < -rhs; i++)
					it_new--;
			}
			return it_new;
		}
		self_type operator-(const int& rhs) {
			return (*this + (-rhs));
		}

		reference operator*() { return *ptr_; }
		pointer operator->() { return ptr_; }
		bool operator==(const self_type& rhs) { return ((cur_block == rhs.cur_block) && (cur_index == rhs.cur_index)); }
		bool operator!=(const self_type& rhs) { return ((cur_block != rhs.cur_block) || (cur_index != rhs.cur_index)); }
	private:
		INTEGER_TYPE cur_block;
		INTEGER_TYPE cur_index;
		INTEGER_TYPE container;
		address_list* address;
		pointer ptr_;
	};

	typedef typename darray<T>::iterator it_type;
	typedef typename darray<T>::const_iterator const_it_type;

	it_type begin() {
		return (size == 0) ?
			it_type(0, 0, container, address, nullptr).get_ref() :
			it_type(0, 0, container, address, (T*)address->ele(0));
	};

	it_type rbegin() { return it_type(-1, -1, container, address, nullptr); };

	it_type end() { return it_type(size / container, size % container, container, address, nullptr); };

	it_type rend() {
		return (size == 0) ?
			this->rbegin() :
			it_type((size - 1) / container, (size - 1) % container, container, address, end_ptr);
	};


	const_it_type cbegin() const {
		return (size == 0) ?
			const_it_type(0, 0, container, address, nullptr) :
			const_it_type(0, 0, container, address, (T*)address->ele(0));
	};

	const_it_type crbegin() const { return const_it_type(-1, -1, container, address, nullptr); };

	const_it_type cend() const { return const_it_type(size / container, size % container, container, address, nullptr); };

	const_it_type crend() const {
		return (size == 0) ?
			this->crbegin() :
			const_it_type((size - 1) / container, (size - 1) % container, container, address, end_ptr);
	};

	void push(const T& ele) {
		if (size % container == 0) {
			if (size < address->get_size() * container){
				end_ptr = (T*)address->ele(size / container);
			}
			else {
				//if (size)
				//	std::cout << "Check point!\n";

				T* vec = new T[container];
				end_ptr = vec;
				address->push(vec);
			}
		}
		else
			end_ptr++;
		assert(end_ptr != nullptr);
		*end_ptr = ele;
		size++;
	}

	friend void swap(darray<T>& lhs, darray<T>& rhs) {
		using std::swap;
		swap(lhs.container, rhs.container);
		swap(lhs.size, rhs.size);
		swap(*lhs.address, *rhs.address);
		swap(lhs.end_ptr, rhs.end_ptr);
	}

	darray<T>& operator=(const darray<T>& rhs) {
		darray<T> temp = darray<T>(rhs);
		swap(*this, rhs);
		return *this;
	}

	T& ele(INTEGER_TYPE index) {
		// start from 0
		return *ele_ptr(index);
	};

	T* ele_ptr(INTEGER_TYPE index) {
		// start from 0
		assert(index < size);
		INTEGER_TYPE q = index / container;
		INTEGER_TYPE r = index - q * container;
		T* found = (T*)address->ele(q);
		found += r;
		return found;
	};

	T& operator[](INTEGER_TYPE index) { return this->ele(index); };

	T* get_vec(INTEGER_TYPE index = 0) { return (T*)address->ele(index); };

	T* get_vec_const(INTEGER_TYPE index = 0) const { return (T*)address->ele(index); };

	address_list* get_address() { return address; };

	INTEGER_TYPE get_container_size() { return container; };

	INTEGER_TYPE get_size() const { return (end_ptr == nullptr) ? 0 : size; };

	T* get_end_ptr() { return end_ptr; };

	void resize(INTEGER_TYPE n) {
		INTEGER_TYPE max_size = address->get_size() * container;
		if (n > max_size) {
			INTEGER_TYPE r = n - max_size; // r > 0
			INTEGER_TYPE q = (r - 1) / container + 1; // q > 0
			r = n - (address->get_size() + q - 1) * container;
			for (INTEGER_TYPE i = 0; i < q; i++)
				address->push(new T[container]);
			end_ptr = (T*)address->back();
			end_ptr += (r - 1);
			size = n;
		}
		else {
			end_ptr = (n != 0) ? &(this->ele(n - 1)) : nullptr;
			size = n;
		}
	}

	void set_value(const T& ele) {
		assert(size != 0);
		INTEGER_TYPE q = (size - 1) / container + 1;
		T* cur_vec;
		for (int i = 0; i < q; i++) {
			cur_vec = (T*)address->ele(i);
			for (int j = 0; j < container; j++)
				cur_vec[j] = ele;
		}
	}

	void resize(INTEGER_TYPE n, const T& ele) {
		this->resize(n);
		this->set_value(ele);
	}

	void erase() { this->resize(0); };

	void release_ptr() {
		// This is used for releasing darray<S*> with pointer to an object S.
		auto it_end = this->end();
		for (auto it = this->begin(); it != it_end; it++) 
			delete *it;
		this->release();
	}

	void release_arr() {
		// This is used for releasing darray<S*> with pointer to an array of S.
		auto it_end = this->end();
		for (auto it = this->begin(); it != it_end; it++) 
			delete[] *it;
		this->release();
	}

	//template<class S>
	//void release(S& junk) {
	//	// This is used for releasing darray<S> with basic data type S.
	//	end_ptr = nullptr;
	//	address->release_vector(junk);
	//	container = DYNAMIC_ARRAY_CONTAINER;
	//	size = 0;
	//}

	void release() {
		// This is used for releasing darray<T> with basic data type S.
		while (!address->is_empty()) {
			T* ptr = (T*)address->pop();
			delete ptr;
		}
		size = 0;
		end_ptr = nullptr;
	}

	bool is_empty() { return (this == nullptr) ? true : (size == 0); };

	T* back_ptr() { return end_ptr; };

	T& back() { return *this->back_ptr(); };

	T* front_ptr() { return this->ele_ptr(0); };

	T& front() { return *(this->front_ptr()); };

private:
	INTEGER_TYPE container;
	INTEGER_TYPE size;
	address_list* address;
	T* end_ptr;

};

template <class T>
class darray_dc {
public:
	darray_dc() : c_size(0), e_size(0), address(new address_list()) {};
	darray_dc(INTEGER_TYPE n) : c_size(0), e_size(0), address(new address_list(n)) {};
	darray_dc(INTEGER_TYPE n, INTEGER_TYPE c) : darray_dc(n) { this->initial_container(c); c_size = n; };
	darray_dc(INTEGER_TYPE n, INTEGER_TYPE c, T& ele) : c_size(n), e_size(n* c), address(new address_list(n)) {
		for (INTEGER_TYPE i = 0; i < c_size; i++) {
			darray<T>* temp = new darray<T>(c);
			temp->resize(c, ele);
			address->set_size(n);
			address->set_ele(i, temp);
		}
	};
	darray_dc(const darray_dc<T>& src) : c_size(src.c_size), e_size(src.e_size) {
		this->address = new address_list(*src.address);
		darray<T>* temp_src;
		for (INTEGER_TYPE i = 0; i < c_size; i++) {
			temp_src = (darray<T>*) src.address->ele(i);
			if (temp_src != nullptr) {
				auto temp = new darray<T>(*temp_src); //deep copy.
				this->address->set_ele(i, temp);
			}
			else
				this->address->set_ele(i, (T*)nullptr);
		}
	}
	~darray_dc() {
		//this->release();
		//std::cout << "Releasing darray_dc with container\t" << c_size << ", entry\t" << e_size << ".\n";
		this->earse_drop_c();
		delete this->address;
	}

	void initial_container(INTEGER_TYPE c = DYNAMIC_ARRAY_CONTAINER) {
		for (INTEGER_TYPE i = 0; i < c_size; i++)
			address->set_ele(i, new darray<T>(c));
	}

	class iterator {
	public:
		typedef typename darray_dc<T>::iterator self_type;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;
		typedef darray<T> c_type;
		typedef darray<T>* c_ptr;
		typedef typename darray<T>::iterator c_iterator;

		darray_dc<T>::iterator() : cur_block(-1), cur_index(-1), cur_iterator(typename darray<T>::iterator().get_ref()), address(nullptr), cur_block_ptr(nullptr) {};
		darray_dc<T>::iterator(INTEGER_TYPE b, INTEGER_TYPE i, address_list* a, c_iterator it, c_ptr cur_b_p) : cur_block(b), cur_index(i), address(a), cur_iterator(it), cur_block_ptr(cur_b_p) {};
		darray_dc<T>::iterator(const self_type& src) : cur_block(src.cur_block), cur_index(src.cur_index), address(src.address), cur_iterator(src.cur_iterator), cur_block_ptr(src.cur_block_ptr) {};

		friend void swap(self_type& lhs, self_type& rhs) {
			using std::swap;
			swap(lhs.cur_block, rhs.cur_block);
			swap(lhs.cur_index, rhs.cur_index);
			swap(lhs.address, rhs.address);
			swap(lhs.cur_iterator, rhs.cur_iterator);
		}

		self_type& operator=(const self_type& rhs) {
			self_type temp = self_type(rhs);
			swap(*this, temp);
			return *this;
		}


		self_type operator++(int) {
			if (cur_index < (cur_block_ptr->get_size() - 1)) { cur_index++; cur_iterator++; }
			else {
				if (cur_block == (address->get_size() - 1)) {
					cur_iterator = cur_block_ptr->end();
					cur_index++;
				}
				else {
					cur_block++;
					cur_index = 0;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					while (cur_block_ptr->is_empty() && cur_block < (address->get_size() - 1)) {
						cur_block++;
						cur_block_ptr = (c_ptr)address->ele(cur_block);
					}
					cur_iterator = cur_block_ptr->begin();
				}
			}
			return *this;
		}
		self_type operator++() {
			self_type i = *this;
			if (cur_index < (cur_block_ptr->get_size() - 1)) { cur_index++; cur_iterator++; }
			else {
				if (cur_block == (address->get_size() - 1)) {
					cur_iterator = cur_block_ptr->end();
					cur_index++;
				}
				else {
					cur_block++;
					cur_index = 0;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					while (cur_block_ptr->is_empty() && cur_block < (address->get_size() - 1)) {
						cur_block++;
						cur_block_ptr = (c_ptr)address->ele(cur_block);
					}
					cur_iterator = cur_block_ptr->begin();

				}
			}
			return i;
		}
		self_type operator--(int) {
			if (cur_index != 0) { cur_index--; cur_iterator--; }
			else {
				if (cur_block == 0) {
					cur_iterator = cur_block_ptr->rbegin();
					cur_index = (INTEGER_TYPE)-1;
					cur_block = (INTEGER_TYPE)-1;
					cur_block_ptr = nullptr;
				}
				else {
					cur_block--;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					cur_index = cur_block_ptr->get_size() - 1;
					cur_iterator = cur_block_ptr->rend();
				}
			}
			return *this;
		}
		self_type operator--() {
			self_type i = *this;
			if (cur_index != 0) { cur_index--; cur_iterator--; }
			else {
				if (cur_block == 0) {
					cur_iterator = cur_block_ptr->rbegin();
					cur_index = (INTEGER_TYPE)-1;
					cur_block = (INTEGER_TYPE)-1;
					cur_block_ptr = nullptr;
				}
				else {
					cur_block--;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					cur_index = cur_block_ptr->get_size() - 1;
					cur_iterator = cur_block_ptr->rend();
				}
			}
			return i;
		}
		reference operator*() { return *cur_iterator; };
		pointer operator->() { return cur_iterator.operator->(); };
		bool operator==(const self_type& rhs) { return ((cur_block == rhs.cur_block) && (cur_index == rhs.cur_index)); }
		bool operator!=(const self_type& rhs) { return ((cur_block != rhs.cur_block) || (cur_index != rhs.cur_index)); }

		INTEGER_TYPE get_block_ind() { return cur_block; };
		INTEGER_TYPE get_entry_ind() { return cur_index; };
	private:
		INTEGER_TYPE cur_block;
		INTEGER_TYPE cur_index;
		address_list* address;
		c_iterator cur_iterator;
		c_ptr cur_block_ptr;
	};

	class const_iterator {
	public:
		typedef typename darray_dc<T>::const_iterator self_type;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;
		typedef darray<T> c_type;
		typedef darray<T>* c_ptr;
		typedef typename darray<T>::const_iterator& c_iterator;

		const_iterator() : cur_block(0), cur_index(0), address(nullptr), cur_iterator(darray<T>::iterator()), cur_block_ptr(nullptr) {};
		const_iterator(INTEGER_TYPE b, INTEGER_TYPE i, address_list* a, c_iterator it, c_ptr cur_b_p) : cur_block(b), cur_index(i), address(a), cur_iterator(it), cur_block_ptr(cur_b_p) {};
		const_iterator(const self_type& src) : cur_block(src.cur_block), cur_index(src.cur_index), address(src.address), cur_iterator(src.cur_iterator), cur_block_ptr(src.cur_block_ptr) {};

		friend void swap(self_type& lhs, self_type& rhs) {
			using std::swap;
			swap(lhs.cur_block, rhs.cur_block);
			swap(lhs.cur_index, rhs.cur_index);
			swap(lhs.address, rhs.address);
			swap(lhs.cur_iterator, rhs.cur_iterator);
		}

		self_type& operator=(const self_type& rhs) {
			self_type temp = self_type(rhs);
			swap(*this, temp);
			return *this;
		}


		self_type operator++(int) {
			if (cur_index < (cur_block_ptr->get_size() - 1)) { cur_index++; cur_iterator++; }
			else {
				if (cur_block == (address->get_size() - 1)) {
					cur_iterator = cur_block_ptr->end();
					cur_index++;
				}
				else {
					cur_block++;
					cur_index = 0;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					while (cur_block_ptr->is_empty() && cur_block < (address->get_size() - 1)) {
						cur_block++;
						cur_block_ptr = (c_ptr)address->ele(cur_block);
					}
					cur_iterator = cur_block_ptr->begin();
				}
			}
			return *this;
		}
		self_type operator++() {
			self_type i = *this;
			if (cur_index < (cur_block_ptr->get_size() - 1)) { cur_index++; cur_iterator++; }
			else {
				if (cur_block == (address->get_size() - 1)) {
					cur_iterator = cur_block_ptr->end();
					cur_index++;
				}
				else {
					cur_block++;
					cur_index = 0;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					while (cur_block_ptr->is_empty() && cur_block < (address->get_size() - 1)) {
						cur_block++;
						cur_block_ptr = (c_ptr)address->ele(cur_block);
					}
					cur_iterator = cur_block_ptr->begin();
				}
			}
			return i;
		}
		self_type operator--(int) {
			if (cur_index != 0) { cur_index--; cur_iterator--; }
			else {
				if (cur_block == 0) {
					cur_iterator = cur_block_ptr->rbegin();
					cur_index = (INTEGER_TYPE)-1;
					cur_block = (INTEGER_TYPE)-1;
					cur_block_ptr = nullptr;
				}
				else {
					cur_block--;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					cur_index = cur_block_ptr->get_size() - 1;
					cur_iterator = cur_block_ptr->rend();
				}
			}
			return *this;
		}
		self_type operator--() {
			self_type i = *this;
			if (cur_index != 0) { cur_index--; cur_iterator--; }
			else {
				if (cur_block == 0) {
					cur_iterator = cur_block_ptr->rbegin();
					cur_index = (INTEGER_TYPE)-1;
					cur_block = (INTEGER_TYPE)-1;
					cur_block_ptr = nullptr;
				}
				else {
					cur_block--;
					cur_block_ptr = (c_ptr)address->ele(cur_block);
					cur_index = cur_block_ptr->get_size() - 1;
					cur_iterator = cur_block_ptr->rend();
				}
			}
			return i;
		}
		reference operator*() { return *cur_iterator; };
		//pointer operator->() { return cur_iterator.operator->(); };
		bool operator==(const self_type& rhs) { return ((cur_block == rhs.cur_block) && (cur_index == rhs.cur_index)); }
		bool operator!=(const self_type& rhs) { return ((cur_block != rhs.cur_block) || (cur_index != rhs.cur_index)); }

		INTEGER_TYPE get_block_ind() { return cur_block; };
		INTEGER_TYPE get_entry_ind() { return cur_index; };
	private:
		INTEGER_TYPE cur_block;
		INTEGER_TYPE cur_index;
		address_list* address;
		c_iterator cur_iterator;
		c_ptr cur_block_ptr;
	};

	typedef typename darray_dc<T>::iterator it_type;
	typedef typename darray_dc<T>::const_iterator const_it_type;

	it_type begin() {
		return (e_size == 0) ?  it_type() : it_type(0, 0, address, ((darray<T>*)address->ele(0))->begin(), (darray<T>*)address->ele(0));
		//if (e_size == 0) {
		//	return it_type();
		//}	
		//else{
		//	darray<T>* ptr_begin = (darray<T>*)address->ele(0);
		//	darray_dc<T>::iterator it(0, 0, address, ptr_begin->begin(), ptr_begin);
		//	return it;
		//}
	};

	it_type rbegin() { return it_type(); };

	it_type end() {
		if (e_size == 0)
			return this->begin();
		else {
			darray<T>* ptr_back = (darray<T>*)address->back();
			return it_type(address->get_size() - 1, ptr_back->get_size(), address, ptr_back->end(), ptr_back);
		}
	}

	it_type rend() {
		if (e_size == 0)
			return this->rbegin();
		else {
			darray<T>* ptr_back = (darray<T>*)address->back();
			return it_type(address->get_size() - 1, ptr_back->get_size() - 1, address, ptr_back->rend(), ptr_back);
		}
	}

	const_it_type cbegin() const {
		return (e_size == 0) ?  const_it_type() : const_it_type(0, 0, address, ((darray<T>*)address->ele(0))->cbegin(), (darray<T>*)address->ele(0));
	};

	const_it_type crbegin() const { return const_it_type(); };

	const_it_type cend() const {
		if (e_size == 0)
			return this->cbegin();
		else {
			darray<T>* ptr_back = (darray<T>*)address->back();
			return const_it_type(address->get_size() - 1, ptr_back->get_size(), address, ptr_back->cend(), ptr_back);
		}
	}

	const_it_type crend() const {
		if (e_size == 0)
			return this->crbegin();
		else {
			darray<T>* ptr_back = (darray<T>*)address->back();
			return const_it_type(address->get_size() - 1, ptr_back->get_size() - 1, address, ptr_back->crend(), ptr_back);
		}
	}

	void push(const T& ele, int c_index) {
		darray<T>* c_ptr = nullptr;
		INTEGER_TYPE ind = (c_index >= 0) ? c_index : c_size + c_index;

		if (this->get_c_ptr(ind) == nullptr)
			this->set_c(ind, new darray<T>());
		
		this->get_c_ptr(ind)->push(ele);
		e_size++;
	}

	void push(const T& ele) { push(ele, c_size - 1); e_size++; };

	void set_c(INTEGER_TYPE i, darray<T>* src_container) {
		assert(i < address->get_max_size());
		INTEGER_TYPE size = 0;
		if (address->ele(i) != nullptr) {
			size = ((darray<T>*)address->ele(i))->get_size();
			if (size)
				std::cout << "Releasing non-empty container in darray_dc by set_container\n";
			delete (darray<T>*)address->ele(i);
		}
		else
			c_size++;
		e_size -= size;
		size = src_container->get_size();
		e_size += size;
		address->set_ele(i, src_container);
	}

	void push_c(darray<T>* c) {
		address->push(c);
		e_size += c->get_size();
		c_size++;
	}

	darray<T>* pop_c() {
		assert(c_size > 0);
		darray<T>* ptr = (darray<T>*)address->pop();
		e_size -= ptr->get_size();
		c_size--;
	}

	darray<T>* pop_c(INTEGER_TYPE i) {
		assert((i >= 0) && (i < c_size));
		darray<T>* ptr = (darray<T>*)address->pop(i);
		e_size -= ptr->get_size();
		c_size--;
	}

	friend void swap(darray_dc<T>& lhs, darray_dc<T>& rhs) {
		using std::swap;
		swap(lhs.c_size, rhs.c_size);
		swap(lhs.e_size, rhs.e_size);
		swap(*lhs.address, *rhs.address);
	}

	void earse_block(INTEGER_TYPE i) {
		assert(i < c_size);
		darray<T>* cur_block = (darray<T>*)address->ele(i);
		if (cur_block->get_size())
			cur_block->resize(0);
	}

	darray_dc<T>& operator=(const darray_dc<T>& rhs) {
		darray_dc<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T* ele_ptr(INTEGER_TYPE c, INTEGER_TYPE i) { auto c_ptr = get_c_ptr(c); return c_ptr->ele_ptr(i); };

	T& ele(INTEGER_TYPE c, INTEGER_TYPE i) { return *ele_ptr(c, i); };
	
	T& ele(INTEGER_TYPE n) {
		// This is not suggested to use for heavy computation. Dynamic container makes the inquiry without
		// container address hard to locate.
		assert(n < e_size);
		INTEGER_TYPE accumulate = 0;
		for (INTEGER_TYPE i = 0; i < address->get_size(); i++) {
			auto c_ptr = get_c_ptr(i);
			if (n < (accumulate + c_ptr->get_size())) // found.
				return c_ptr->ele(n - accumulate);
			else
				accumulate += c_ptr->get_size();
		}
		throw(1);
	}

	T& operator()(INTEGER_TYPE c, INTEGER_TYPE i) { return this->ele(c, i); };


	// This inquiry is not recommended, especially for repeated inequires. Since the indices convertion is nontrivial 
	// for dynamic container, which is called in every inquiry.
	T& operator()(INTEGER_TYPE n) { return this->ele(n); };

	darray<T>* get_c_ptr(INTEGER_TYPE c) { return (darray<T>*)address->ele(c); };

	darray<T>& get_c(INTEGER_TYPE c) { return *(this->get_c_ptr(c)); };

	darray<T>& operator[](INTEGER_TYPE c) { return get_c(c); };

	darray<T>* back_c_ptr() { return this->get_c_ptr(c_size - 1); };

	darray<T>& back_c() { return *(this->back_c_ptr()); };

	T* back_ptr() { return this->back_c_ptr()->back_ptr(); }

	T* back() { return &(this->back_ptr()); }

	INTEGER_TYPE get_c_size() { return c_size; };

	INTEGER_TYPE get_size() { return e_size; };

	bool is_empty() { return (e_size == 0); };

	void print(std::ostream& output, char c_breaker = '\n') {
		for (INTEGER_TYPE i = 0; i < c_size; i++) {
			darray<T>* c_ptr = this->get_c_ptr(i);
			auto it_end = c_ptr->end();
			INTEGER_TYPE size = c_ptr->get_size(), cur_ind = 0;
			for (auto it = c_ptr->begin(); it != it_end; it++, cur_ind++){
				if (cur_ind == size - 1)
					output << *it << c_breaker;
				else
					output << *it << '\t';
			}
		}

	}

	//void release() {
	//	INTEGER_TYPE n = address->get_max_size();
	//	darray<T>* ptr;
	//	T junk;
	//	for (INTEGER_TYPE i = 0; i < n; i++) 
	//		delete this->get_c_ptr(i);
	//		
	//	c_size = 0;
	//	e_size = 0;
	//	delete address;
	//	address = nullptr;
	//}

	void earse_drop_c() {
		while (!address->is_empty()) {
			darray<T>* ptr = (darray<T>*)address->pop();
			delete ptr;
		}
		c_size = 0;
		e_size = 0;
	}

	void earse_keep_c() {
		for (int i = 0; i < address->get_size(); i++) {
			darray<T>* ptr = (darray<T>*)address->ele(i);
			if (ptr != nullptr)
				if (!ptr->is_empty())
					ptr->resize(0);
		}
		c_size = 0;
		e_size = 0;
	}




private:
	INTEGER_TYPE c_size;
	INTEGER_TYPE e_size;
	address_list* address;

};


//#include "zd_dynamic_array.cpp"