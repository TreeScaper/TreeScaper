
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

//warray.h
// definition of a array class
//             by whuang Nov/11/2008

#ifndef WARRAY_H
#define WARRAY_H

#include <iomanip>
#include <vector>
#include <iostream>
#include <cassert>
#include "wfile.h"

template<class T>
class Array{

    friend std::istream &operator>>(std::istream &input, Array<T> &arr)                  // overload >>
	{
		return input;                                                          // do nothing
	};

    friend std::ostream &operator<<(std::ostream &output, const Array<T> &arr)           // overload <<
	{
		output << "{" << arr.length << ": ";

		for(int i = 0; i < arr.length; i++)
			output << arr.vec[i] << ", ";
		output << "}";

		return output;                                                         // output form, for example {length: a1, a2, a3, }
	};

	friend Array operator+(const Array<T> &left, const Array<T> &right)
	{
		Array<T> result;
		result = left;
		result += right;
		return result;
	};                                        // concatenation

	friend Array operator-(const Array<T> &left, const Array<T> &right)
	{
		Array<T> result;
		result = left;
		result -= right;
		return result;
	};                                         // leave the elements which the right don't have.


public:

	template <class T2, int N>
	Array(const T2 (&arr)[N])                                                   // intial and copy a nomal array to a Array
	{
		length = SIZE <N> ::cnt;
		vec = new T2[length];
		for(int i = 0; i < length; i++)
		{
			vec[i] = *(arr + i);
		}
	};

	template <class T2>
	Array(Array<T2> &arr)
	{
		length = arr.length;
		vec = new T2[length];

		for(int i = 0; i < length; i++)
		{
			vec[i] = arr.vec[i];
		}
	}

	Array(const Array &arr)
	{
		length = arr.length;
		vec = new T[length];

		for(int i = 0; i < length; i++)
		{
			vec[i] = arr.vec[i];
		}
	}

	template<class S>
	Array(int size, S arr)
	{
        length = size;

        if(length > 0)
            vec = new S[length];
        else
            vec = NULL;

		for(int i = 0; i < length; i++)
		{
			vec[i] = arr;
		}
	};

    template<class S>
    Array(int size, S* arr)
    {
        length = size;

        if (length > 0)
            vec = new S[length];
        else
            vec = NULL;

        memcpy(vec, arr, sizeof(S) * size);
    };

	Array(int size)
	{
		length = size;
        vec = new T[length];
        for(int i = 0; i < length; i++)
            vec[i] = NULL;
	};

	~Array()
	{
		if(length > 0)
		{
			delete [] vec;
        }
	};

	Array(){length = 0; };                                                     // construct a dafault Array

	const Array &operator=(const Array<T> &);                                  // assignment

	const Array &operator+=(const Array<T> &);                                 // concatenation

	const Array &operator-=(const Array<T> &);                                 // also leave the elements which the right don't have
	
	bool include(const T) const;                                               // check if this Array include the element

	bool operator==(const Array<T> &) const;                                         // compare

    friend operator==(const Array<T> &lhs, const Array<T> &rhs){
        if(lhs.get_length() != rhs.get_length())                                     // if they don't have the same size, they are not same things.
		    return false;

        length = lhs.get_length();

	    for(int i = 0; i < length; i++)                                // compare every element
		    if(lhs.vec[i] != rhs.vec[i])
			    return false;

	    return true;
    }

	bool operator<(const Array<T> &) const;

	Array operator()(const int, const int);

	bool operator!=(const Array<T> &right) const{return ! ((*this) == right);};

	bool operator>(const Array<T> &right) const{return right < (*this);};

	bool operator<=(const Array<T> &right) const{return ! (right < (*this));};

	bool operator>=(const Array<T> &right) const{return ! (right > (*this));};

	bool is_empty(){return length == 0;};                                      // check if this Array empty;

	T &operator[](const int index){return vec[index];};                               // pick the (index + 1)-th element

	const T &operator[](const int index) const{return vec[index];};                   // pick the (index + 1)-th element

	operator T *() const{return vec;};                                                // translate to array.

	int get_length() const {return length;};                                          // get length of Array

	void add(T element);

	void resize(int, int = 0);                                                          // reset the length of a Array

    static bool bitAcontainB(const Array<T> &A, const Array<T> &B, int BITSlength, bool &tf)
    {
        // if the 1 bit in A contains all 1 bit in B, tf is true and return true,
        // if the 1 bit in B contains all 1 bit in A, tf is false and return true,
        // otherwise return false.
        int bitnum = sizeof(T) * 8;
        int quotient = (int) ((BITSlength - 1) / bitnum);
        int remainder = BITSlength - bitnum * quotient;
        int flag = 0; // 1 means A should contain B, -1 means B should contain A.
        int bitlength = 0;
        tf = false;
        T BITS = 0;
        T BITA, BITB;

        for(int i = 0; i <= quotient; i++)
        {
            if(i < quotient)
                bitlength = bitnum;
            else
                bitlength = remainder;

            BITS = 1;
            for(int j = 0; j < bitlength; j++)
            {
                BITA = BITS & A[i];
                BITB = BITS & B[i];
                if(flag == 0 && BITA && ! BITB)
                    flag = 1;
                if(flag == 0 && !BITA && BITB)
                    flag = -1;
                if(flag == -1 && BITA && ! BITB)
                    return false;
                if(flag == 1 && !BITA && BITB)
                    return false;

                BITS *= 2;
            }
        }
        if(flag == 1)
            tf = true;
        else if(flag == -1)
            tf = false;

        return true;
    };

    void Flipbit(int BITSlength)
    {
        int bitnum = sizeof(T) * 8;
        int quotient = (int) ((BITSlength - 1) / bitnum);
        int remainder = BITSlength - bitnum * quotient;
        T BITS = 1;

        for(int i = 0; i < quotient; i++)
            vec[i] = ~(vec[i]);

        for(int i = 0; i < remainder; i++)
        {
            if(vec[quotient] & BITS)
                vec[quotient] = vec[quotient] & (~BITS);
            else
                vec[quotient] = vec[quotient] | BITS;
            BITS *= 2;
        }
    };

        bool ORbitOPE(const Array<T> &left, const Array<T> &right)
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
                vec[i] = NULL;

            for(int i = 0; i < length; i++)
                vec[i] = left.vec[i] | right.vec[i];

//            std::cout << "left:";
//            left.printbits(12);
//            std::cout << std::endl;
//            std::cout << "right:";
//            right.printbits(12);
//            std::cout << std::endl;
//            std::cout << "OR:";
//            this->printbits(12);
//            std::cout << std::endl;

            return true;
        };

        bool SetBitArray(int idx)
        {
            int bitnum = sizeof(T) * 8;
            if(idx > bitnum * length)
            {
                std::cout << "Warning(warray.cpp): Wrong index! idx:" << idx << ", nbits:" << bitnum * length << std::endl;
                return false;
            }

            idx--;
            int quotient = (int) idx / bitnum;
            int remainder = idx - bitnum * quotient;
            T BITS = 0;
            if(remainder >= 0)
                BITS = 1;
            for(int i = 0; i < remainder; i++)
                BITS *= 2;

            vec[quotient] = vec[quotient] | BITS;

            return true;
        };

        void printbits(int outputlength) const
        {
            int bitnum = 8 * sizeof(T);
            T BITS = 0;
            int idx = 0;
            for(int i = 0; i < length; i++)
            {
                BITS = 1;
//std::cout << std::endl; std::cout << "i:" << i << ":" << (int) vec[i] << std::endl;//---
                for(int j = 0; j < bitnum; j++)
                {
                    if(BITS & vec[i])
                        std::cout << 1;
                    else
                        std::cout << 0;
                    BITS *= 2;
                    idx++;
                    if(idx == outputlength)
                        break;
                }
                if(idx == outputlength)
                    break;
            }
        }
        void printbits(int outputlength)
        {
            int bitnum = 8 * sizeof(T);
            T BITS = 0;
            int idx = 0;
            for(int i = 0; i < length; i++)
            {
                BITS = 1;
//std::cout << std::endl; std::cout << "i:" << i << ":" << (int) vec[i] << std::endl;//---
                for(int j = 0; j < bitnum; j++)
                {
                    if(BITS & vec[i])
                        std::cout << 1;
                    else
                        std::cout << 0;
                    BITS *= 2;
                    idx++;
                    if(idx == outputlength)
                        break;
                }
                if(idx == outputlength)
                    break;
            }
        }
        void printbits(int outputlength, File &fout)
        {
            int bitnum = 8 * sizeof(T);
            T BITS = 0;
            int idx = 0;
            for (int i = 0; i < length; i++)
            {
                BITS = 1;
                //std::cout << std::endl; std::cout << "i:" << i << ":" << (int) vec[i] << std::endl;//---
                for (int j = 0; j < bitnum; j++)
                {
                    if (BITS & vec[i])
                        fout << 1;
                    else
                        fout << 0;
                    BITS *= 2;
                    idx++;
                    if (idx == outputlength)
                        break;
                }
                if (idx == outputlength)
                    break;
            }
        }
        void printbits(int outputlength, std::ostream &fout)
        {
            int bitnum = 8 * sizeof(T);
            T BITS = 0;
            int idx = 0;
            for (int i = 0; i < length; i++)
            {
                BITS = 1;
                //std::cout << std::endl; std::cout << "i:" << i << ":" << (int) vec[i] << std::endl;//---
                for (int j = 0; j < bitnum; j++)
                {
                    if (BITS & vec[i])
                        fout << 1;
                    else
                        fout << 0;
                    BITS *= 2;
                    idx++;
                    if (idx == outputlength)
                        break;
                }
                if (idx == outputlength)
                    break;
            }
        }

	template<class S>
	void resize(int size, S arr)
	{
		int sub_length = 0;
		length = size;
        this->resize(length);
		sub_length = arr.get_length();

		for(int i = 0; i < length; i++)
			(vec[i]).resize(sub_length, arr[0]);
	};

    static bool bitstrXOR(const Array<T> &left, const Array<T> &right, int &result);

    static bool onebitstrXOR(const Array<T> &left, int &result);

private:

	template <int N>
	struct SIZE{                                                               // as a tool to store the length of a array
		static const int cnt = N;
	};

    T *vec;
    //std::vector<T> vec;
	int length;
};

#endif
