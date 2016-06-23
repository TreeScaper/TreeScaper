
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

// wmix.h
// definition of a mix data type
//                by whuang Dec/05/2008

#ifndef WMIX_H
#define WMIX_H

#include <iostream>
#include <cassert>
#include "wstring.h"
#include "warray.cpp"
#include "wmapping.cpp"

using namespace std;

class Mix{
	friend ostream &operator<<(ostream &output, const Mix &mix)
	{
		if(mix.data_type == INT)
		{
			int *ptr = (int *) mix.point;
			output << (*ptr);
		} else
		if(mix.data_type == DOUBLE)
		{
			double *ptr = (double *) mix.point;
			output << (*ptr);
		} else
		if(mix.data_type == STRING)
		{
			String *ptr = (String *) mix.point;
			output << (*ptr);
		} else
		if(mix.data_type == ARRAY)
		{
			Array<Mix> *ptr = (Array<Mix> *) mix.point;
			output << (*ptr);
		} else
		if(mix.data_type == MAPPING)
		{
			Mapping<Mix, Mix> *ptr = (Mapping<Mix, Mix> *) mix.point;
			output << (*ptr);
		}

		return output;
	};

public:
	enum types{NONE, INT, DOUBLE, STRING, ARRAY, MAPPING};

	Mix(){data_type = NONE;};

	Mix(const Mix &data){(*this) = data;};

	Mix(int data);

	Mix(double data);

	Mix(String data);

	Mix(char *data);

	Mix(Array<Mix> data);

	Mix(Mapping<Mix, Mix> data);

	const Mix &operator=(const Mix &right);

	int operator=(int right);

	Mix operator=(double right);

	const Mix &operator=(const String &right);

	Mix operator=(char *right);

	const Mix &operator=(const Array<Mix> &right);

	const Mix &operator=(const Mapping<Mix, Mix> &right);

	bool operator<(const Mix &right) const;

	bool operator>(const Mix &right) const;

	bool operator<=(const Mix &right) const;

	bool operator>=(const Mix &right) const;

	bool operator!=(const Mix &right) const;

	bool operator==(const Mix &right) const;

	void release_data();

	~Mix();

	operator int() const;

	operator double() const;

	operator String() const;

	operator char *() const;

	operator Array<Mix>() const;

	operator Mapping<Mix, Mix>() const;

	int get_type();

private:
	void *point;

	types data_type;	//enum types{NONE, INT, DOUBLE, STRING, ARRAY, MAPPING};
};

#endif
