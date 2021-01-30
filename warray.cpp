
//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1   
//# Copyright (C) 2010 Wen Huang
//# 
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details. 
//# http://www.gnu.org/copyleft/gpl.html 
//##########################################################################

// warray.cpp
// Member function definitions for class warray.h
//              by whuang Nov/11/2008

#ifndef WARRAY_CPP
#define WARRAY_CPP

#undef max
#undef min

#define mymin(a,b) ((a)<(b)?(a):(b))
#define mymax(a,b) ((a)>(b)?(a):(b))

#include "warray.h"

// assignment
template<class T>
const Array<T> &Array<T>::operator=(const Array<T> &right)
{
	if(&right != this)
	{
		T *left = vec;
		if(length > 0)
        {
			delete [] left;
            left = NULL;
        }

        length = right.length;
        if(length > 0)
            vec = new T[length];
        else
            vec = NULL;

		for(int i = 0; i < length; i++)
		{
			vec[i] = right.vec[i];
		}
	}

	return *this;
}

// concatenation
template<class T>
const Array<T> &Array<T>::operator+=(const Array<T> &right)
{
//	int left_length = length;
//	T *left = vec;
//	length += right.length;
//	vec = new T[length];
	
	Array<T> result(length + right.length);

	//vec.resize(length);

	for(int i = 0; i < length; i++)                    // assignment every elements to the left elements
	{
		result[i] = vec[i];
	}
	
	for(int i = length; i < length + right.length; i++)
	{
		result[i] = right.vec[i - length];
	}

	return (*this) = result;
}

// also leave the elements which the right don't have
template<class T>
const Array<T> &Array<T>::operator-=(const Array<T> &right)
{
	Array<T> result;
	result.length = 0;

	for(int i = 0; i < length; i++)                             // check every element of left Array
	{
		if(! (right.include(vec[i])))                           // if i-th element of left Array
		{                                                       // isn't included in the right, Then it should be left for the result.
			result.resize(result.length + 1);
			//--result.vec.resize(result.length);
			result.vec[result.length - 1] = vec[i];
		}
	}
	return (*this) = result;
}

// check if the Array include value
template<class T>
bool Array<T>::include(const T value) const
{
	for(int i = 0; i < length; i++)
	{
		if(vec[i] == value)
		return true;
	}
	return false;
}

// template<class T>
// bool Array<T>::operator==(const Array<T> &right) const
// {
// 	if(length != right.length)                                     // if they don't have the same size, they are not same things.
// 		return false;

// 	for(int i = 0; i < length; i++)                                // compare every element
// 		if(vec[i] != right.vec[i])
// 			return false;

// 	return true;
// }

template<class T>
bool Array<T>::operator<(const Array<T> &right) const
{
	if((*this) == right)                                            // same things, return false
		return false;

    int i, size = mymin(right.length, length);                          // get minimum length of them

	for(i = 0; i < size; i++)                                   // if the first greater element belong to left Array, Then left Array is greater 
		if(vec[i] > right.vec[i])
			return false;

	if(length > right.length)
	{
		for(i = 0; i < right.length; i++)
		{
			if((*this)[i] != right[i])
				break;
		}
		if(i == right.length)
			return false;
	}

//	if(length > right.length && (*this)(0, (right.length) - 1) == right)   // the same length elements are unique., then the longer Array should be greater.
//		return false;

	return true;                                                    // otherwise
}

// pick some element
template<class T>
Array<T> Array<T>::operator()(const int index, const int end)       // pick the part of Array according to the index.
{
	assert(index >=0 && end < length);

	Array<T> result;
	result.resize(end - index + 1);
	for(int i = index; i <= end; i++)
		result.vec[i - index] = vec[i];

	return result;
}

template<class T>
void Array<T>::add(T element)
{
        (*this).resize((int) ((*this).get_length()) + 1);
	(*this)[(*this).get_length() - 1] = element;
}

template<class T>
void Array<T>::resize(int size, int flag)                          // change the size of the Array.
{
	/*
	Array<T> result(size);
	int min = (length > size) ? size : length;
	for(int i = 0; i < min; i++)
	{
		result.vec[i] = vec[i];
	}
*/
    T *org_vec = vec;

    if(size > 0)
        vec = new T[size];
    else
        vec = NULL;

    int min = (length > size) ? size : length;
//        std::cout << "min: " << min << std::endl;//---
	for(int i = 0; i < min; i++)
	{
		vec[i] = org_vec[i];
    }

	if(length > 0)
		delete [] org_vec;
	length = size;

    flag = 1;
}

template<class T>
bool Array<T>::bitstrXOR(const Array<T> &left, const Array<T> &right, int &result)
{
    if(left.get_length() != right.get_length())
    {
        std::cout << "Warning(warray.cpp): Left bitstring and right bitstring are not the same length!" << std::endl;
        return false;
    }

    Array<T> resultarr(left.get_length());

    for (int i = 0; i < left.get_length(); i++)
        resultarr.vec[i] = left.vec[i]^right.vec[i];

    result = 0;

    int bitnum = 8 * sizeof(T);
    T BITS = 0;
    for(int i = 0; i < left.get_length(); i++)
    {
        BITS = 1;
        for(int j = 0; j < bitnum; j++)
        {
            if(BITS & resultarr.vec[i])
                result++;
            BITS *= 2;
        }
    }

    return true;
}

template<class T>
bool Array<T>::onebitstrXOR(const Array<T> &left, int &result)
{
    if(left.get_length() == 0)
    {
        std::cout << "Warning(warray.cpp): Bitstring cannot have length zero!" << std::endl;
        return false;
    }

    result = 0;
    int bitnum = 8 * sizeof(T);
    T BITS = 0;
    for(int i = 0; i < left.get_length(); i++)
    {
        BITS = 1;
        for(int j = 0; j < bitnum; j++)
        {
            if(BITS & left.vec[i])
                result++;
            BITS *= 2;
        }
    }

    return true;
}

/*
template<class T>
bool Array<T>::ORbitOPE(const Array<T> &left, const Array<T> &right)
{
    if(left.get_length() != right.get_length())
    {
        std::cout << "Warning(warray.cpp): Left bitstring and right bitstring are not the same length!" << std::endl;
        return false;
    }
    if(length > 0)
        delete [] vec;
    length = left.get_length();

    vec = new T [length];

    for(int i = 0; i < length; i++)
        vec[i] = left.vec[i] | right.vec[i];

    return true;
}*/
/*
template<class T>
bool Array<T>::SetBitArray(int idx)
{
    if(idx > sizeof(T) * length)
    {
        std::cout << "Warning(warray.cpp): Wrong index!" << std::endl;
        return false;
    }

    int quotient = (int) idx / sizeof(T);
    int remainder = idx - sizeof(T) * quotient;
    T BITS = 0;
    if(remainder >= 1)
        BITS = 1;
    for(int i = 1; i < remainder; i++)
        BITS *= 2;

    vec[quotient] = vec[quotient] | BITS;

    return true;
}*/
/*
template<class T>
bool Array<T>::SetBitArray(int idx)
{
    if(idx > sizeof(T) * length)
    {
        std::cout << "Warning(warray.cpp): Wrong index!" << std::endl;
        return false;
    }

    int quotient = (int) idx / sizeof(T);
    int remainder = idx - sizeof(T) * quotient;
    T BITS = 0;
    if(remainder >= 1)
        BITS = 1;
    for(int i = 1; i < remainder; i++)
        BITS *= 2;

    vec[quotient] = vec[quotient] | BITS;
    return true;
};*/

#endif
