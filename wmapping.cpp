
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

// wmapping.cpp
//      member functions
//           by whuang Jan/04/2008

#ifndef WMAPPING_CPP
#define WMAPPING_CPP

#include "wmapping.h"

template<class T, class S>
const Mapping<T, S> &Mapping<T, S>::operator=(const Mapping<T, S> &right)
{
    (*this).mapping = right.mapping;
	return (*this);
};

template<class T, class S>
const Mapping<T, S> &Mapping<T, S>::operator+=(const Mapping<T, S> &right)
{
	(*this) = (*this) + right;
	return (*this);
};

template<class T, class S>
const Mapping<T, S> &Mapping<T, S>::operator-=(const Mapping<T, S> right)
{
	(*this) = (*this) - right;
	return (*this);
};

template<class T, class S>
bool Mapping<T, S>::operator<(const Mapping<T, S> &right) const
{
	return (*this).Keys() < right.Keys();
}

template<class T, class S>
bool Mapping<T, S>::operator>(const Mapping<T, S> &right) const
{
	return right.Keys() < (*this).Keys();
}

template<class T, class S>
bool Mapping<T, S>::operator<=(const Mapping<T, S> &right) const
{
	return ! (right.Keys() < (*this));
}

template<class T, class S>
bool Mapping<T, S>::operator>=(const Mapping<T, S> &right) const
{
	return ! ((*this).Keys() < right.Keys());
}

template<class T, class S>
bool Mapping<T, S>::operator!=(const Mapping<T, S> &right) const
{
	return ! ((*this) == right);
}

template<class T, class S>
bool Mapping<T, S>::operator==(const Mapping<T, S> &right) const
{
	return ! ((*this) < right || (*this) > right);
}

template<class T, class S>
Array<T> Mapping<T, S>::Keys() const
{
    Array<T> arr;
    map<T, S> m(mapping);
    typename map<T, S>::iterator iter = m.begin();
    for(iter = (m).begin(); iter != (m).end(); iter++)
        arr.add(iter->first);

	return arr;
};

template<class T, class S>
Array<S> Mapping<T, S>::Values() const
{
	Array<S> arr;
	map<T, S> m = mapping;
	typename map<T, S>::iterator iter;
	for(iter = (m).begin(); iter != (m).end(); iter++)
		arr.add(iter->second);

	return arr;
};

template<class T, class S>
S &Mapping<T, S>::operator[](T key)
{
	cout << "";
	return mapping[key];
};

template<class T, class S>
bool Mapping<T, S>::Include_Key(T key)
{
	typename map<T, S>::iterator iter;
	iter = mapping.find(key);
	if(iter == mapping.end())
		return false;
	return true;
}

template<class T, class S>
int Mapping<T, S>::get_length() const
{
	return mapping.size();
};

#endif
