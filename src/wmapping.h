
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

// wmapping.h
//     a kind of storage
//           by whuang Jan/03/2008

#ifndef WMAPPING_H
#define WMAPPING_H

#include <map>
#include <iostream>
#include "warray.cpp"

using namespace std;

template<class T, class S>
class Mapping{

	friend ostream &operator<<(ostream &output, Mapping<T, S> right)
    {
		output << "[";
        Array<T> keys = right.Keys();
//		class Mix key, value;
		T key;
		S value;
		output << right.get_length() << ":\n";
		int l = keys.get_length();
		for(int i = 0; i < l; i++)
        {
			key = keys[i];
			value = right[key];
			if(i < l - 1)
				output << "  " << key << " : " << value << ",\n";
			else
				output << "  " << key << " : " << value << "\n";
		}

		output << "]";

		return output;
	};

	friend Mapping operator+(Mapping<T, S> left, Mapping<T, S> right)
	{
		Mapping<T, S> result;
		map<T, S> result_map;
		typename map<T, S>::iterator iter;
		for(iter = (left.mapping).begin(); iter != (left.mapping).end(); iter++)
		{
			result_map[iter->first] = iter->second;
		}

		for(iter = (right.mapping).begin(); iter != (right.mapping).end(); iter++)
		{
			result_map[iter->first] = iter->second;
		}

		result.mapping = result_map;

		return result;
	};
	
	friend Mapping operator-(Mapping<T, S> left, Mapping<T, S> right)
	{
		Mapping<T, S> result;

		Array<T> arr = left.Keys() - right.Keys();

		for(int i = 0; i < arr.get_length(); i++)
			(result.mapping)[arr[i]] = (left.mapping)[arr[i]];

		return result;
	};

public:
	Mapping(){};

	const Mapping &operator=(const Mapping<T, S> &right);

	const Mapping &operator+=(const Mapping<T, S> &right);

	const Mapping &operator-=(const Mapping<T, S> right);

	bool operator<(const Mapping<T, S> &right) const;

	bool operator>(const Mapping<T, S> &right) const;

	bool operator<=(const Mapping<T, S> &right) const;

	bool operator>=(const Mapping<T, S> &right) const;

	bool operator!=(const Mapping<T, S> &right) const;

	bool operator==(const Mapping<T, S> &right) const;

	Array<T> Keys() const;

	Array<S> Values() const;

	S &operator[](T key);

	bool Include_Key(T key);

	int get_length() const;

private:
	map<T, S> mapping;
};

#endif
