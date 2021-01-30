#pragma once

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cmath>


template<class T>
class Array {
private:
	friend std::ostream& operator<<(std::ostream& output, const Array<T>& arr)           // overload <<
	{
		output << "{" << arr.length << ": ";

		for (int i = 0; i < arr.length; i++)
			output << arr.vec[i] << ", ";
		output << "}";

		return output;                                                         // output form, for example {length: a1, a2, a3, }
	};

	int length;
	int container_size;
	int max_size;
	T* vec;
public:

	friend void swap(Array<T>& lhs, Array<T>& rhs)
	{
		using std::swap;
		swap(lhs.length, rhs.length);
		swap(lhs.container_size, rhs.container_size);
		swap(lhs.max_size, rhs.max_size);
		swap(lhs.vec, rhs.vec);
	}

	Array(const T*arr, int size, int container = 100)
		: length(size),
		container_size(container),
		max_size(length + (-length) % container_size),
		vec(max_size ? new T[max_size] : nullptr)
	{
		for (int i = 0; i < length; i++)
			vec[i] = arr[i];
		memset(vec + length, 0, (max_size - length) * sizeof(T));
	};

	//Copy-constructor.
	Array(const Array<T>& arr)
		: length(arr.length),
		container_size(arr.container_size),
		max_size(arr.max_size),
		vec(max_size ? new T[max_size] : nullptr)
	{
		for (int i = 0; i < max_size; i++)
			vec[i] = arr(i);
	}

	////Copy-constructor (different type).
	//template <class S>
	//Array(const Array<S>& arr)
	//	: length(arr.length),
	//	container_size(arr.container_size),
	//	max_size(arr.max_size),
	//	vec(max_size ? new T[max_size] : nullptr)
	//{
	//	for(int i = 0; i < max_size; i++)
	//		vec[i] = (T) arr.vec[i]
	//}

	Array(int size, int container)
		: length(size),
		container_size(container),
		max_size(length + (-length) % container_size),
		vec(max_size ? new T[max_size] : nullptr)
	{
		memset(vec, 0, max_size * sizeof(T));
	}

	Array(int size, int container, T ele)
		: length(size),
		container_size(container),
		max_size(length + (-length) % container_size),
		vec(max_size ? new T[max_size] : nullptr)
	{
		for (int i = 0; i < size; i++)
			vec[i] = ele;
	}

	Array(int size) : length(size), container_size(size), max_size(size), vec(length ? new T[length] : nullptr) {};

	Array() : length(0), container_size(10), max_size(0), vec(nullptr) {};

	void push(const T& ele) {
		if (length == max_size) {
			max_size += container_size;
			T* temp = new T[max_size];
			for (int i = 0; i < length; i++) 
				temp[i] = vec[i];
			//memcpy(temp, vec, length * sizeof(T)); // Cause copy constructor and destructor issues.
			delete[] vec;
			vec = temp;
		}
		vec[length] = ele;
		length++;
	}

	Array<T>& operator=(const Array<T>& rhs)
	{
		Array<T> temp(rhs);
		swap(*this, temp);
		return *this;
	}

	T& operator[](const int index) { return vec[index]; };

	T operator()(const int index) const { return vec[index]; };

	T& front() { return vec[0]; };

	T& back() { return vec[length - 1]; };

	T pop() { return vec[length - 1]; length--; };

	T* get_vec() const { return vec; };


	void erase() { length = 0; };

	void release() { delete[] vec; vec = nullptr; length = 0; max_size = container_size; }

	void resize(int size, int container, bool flag_reset = false);

	void align() { resize(this->length, this->length); };

	int get_size() const { return length; };

	bool is_empty() { return (length == 0); };
};