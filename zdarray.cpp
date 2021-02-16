#include "zdarray.hpp"


template <class T>
void Array<T>::resize(int size, int container, bool flag_set0) {
	size_t max = ((size + (container - 1)) / container) * container;
	if (max == max_size) {
		length = size;
		container_size = container;
		if (flag_set0)
			memset(vec, 0, length * sizeof(T));
	}
	else {
		T* temp = new T[max];
		if (flag_set0)
			memset(temp, 0, size * sizeof(T));
		else {
			int loop = (size < length ? size : length);
			for (int i = 0; i < loop; i++)
				vec[i] = temp[i];
		}
			
		length = size;
		container_size = container;
		max_size = max;
		delete[] vec;
		vec = temp;
	}
}