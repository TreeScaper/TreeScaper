
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

// wmix.cpp
// difinition of member function of mix.h
//                by whuang Dec/5th/2008

#ifndef WMIX_CPP
#define WMIX_CPP

#include "wmix.h"

Mix::Mix(int data)
{
	int *ptr = new int;
	*ptr = data;
	point = (void *) ptr;
	data_type = INT;
};

Mix::Mix(double data)
{
	double *ptr = new double;
	*ptr = data;
	point = (void *) ptr;
	data_type = DOUBLE;
}

Mix::Mix(String data)
{
	String *ptr = new String;
	*ptr = data;
	point = (void *) ptr;
	data_type = STRING;
}

Mix::Mix(char *data)
{
	String *ptr = new String;
	*ptr = data;
	point = (void *) ptr;
	data_type = STRING;
}

Mix::Mix(Array<Mix> data)
{
	Array<Mix> *ptr = new Array<Mix>;
	*ptr = data;
	point = (void *) ptr;
	data_type = ARRAY;
}

Mix::Mix(Mapping<Mix, Mix> data)
{
	Mapping<Mix, Mix> *ptr = new Mapping<Mix, Mix>;
	*ptr = data;
	point = (void *) ptr;
	data_type = MAPPING;
}

const Mix &Mix::operator=(const Mix &right)
{
	(*this).release_data();
	if(right.data_type == NONE)
		data_type = right.data_type;
	else
	if(right.data_type == INT)
	{
		int *ptr = new int;
		*ptr = *((int *) right.point);
		point = (void *) ptr;
		data_type = right.data_type;
	} else
	if(right.data_type == DOUBLE)
	{
		double *ptr = new double;
		*ptr = *((double *) right.point);
		point = (void *) ptr;
		data_type = right.data_type;
	} else
	if(right.data_type == STRING)
	{
		String *ptr = new String;
		*ptr = (*((String *) right.point));
		point = (void *) ptr;
		data_type = right.data_type;
	} else
	if(right.data_type == ARRAY)
	{
		Array<Mix> *ptr = new Array<Mix>;
		*ptr = *((Array<Mix> *) right.point);
		point = (void *) ptr;
		data_type = right.data_type;
	} else
	if(right.data_type == MAPPING)
	{
		Mapping<Mix, Mix> *ptr = new Mapping<Mix, Mix>;
		*ptr = *((Mapping<Mix, Mix> *) right.point);
		point = (void *) ptr;
		data_type = right.data_type;
	}

	return (*this);
}

int Mix::operator=(int right)
{
	(*this).release_data();
	int *ptr = new int;
	*ptr = right;
	point = (void *) ptr;
	data_type = INT;
	return (*this);
}

Mix Mix::operator=(double right)
{
	(*this).release_data();
	double *ptr = new double;
	*ptr = right;
	point = (void *) ptr;
	data_type = DOUBLE;
	return (*this);
}

const Mix &Mix::operator=(const String &right)
{
	(*this).release_data();
	String *ptr = new String;
	*ptr = right;
	point = (void *) ptr;
	data_type = STRING;
	return (*this);
}

Mix Mix::operator=(char *right)
{
	(*this).release_data();
	String *ptr = new String;
	*ptr = right;
	point = (void *) ptr;
	data_type = STRING;
	return (*this);
}

const Mix &Mix::operator=(const Array<Mix> &right)
{
	(*this).release_data();
	Array<Mix> *ptr = new Array<Mix>;
	*ptr = right;
	point = (void *) ptr;
	data_type = ARRAY;
	return (*this);
}

const Mix &Mix::operator=(const Mapping<Mix, Mix> &right)
{
	(*this).release_data();
	Mapping<Mix, Mix> *ptr = new Mapping<Mix, Mix>;
	*ptr = right;
	point = (void *) ptr;
	data_type = MAPPING;
	return (*this);
}

bool Mix::operator<(const Mix &right) const
{
	if(data_type == INT)
	{
		int *ptr = (int *) point;
		if(right.data_type == INT)
		{
			int *ptr1 = (int *) right.point;
			return (*ptr) < (*ptr1);
		} else
		if(right.data_type == DOUBLE)
		{
			double *ptr1 = (double *) right.point;
			return (double) (*ptr) < (*ptr1);
		}
	} else
	if(data_type == DOUBLE)
	{
		double *ptr = (double *) point;
		if(right.data_type == INT)
		{
			int *ptr1 = (int *) right.point;
			return (*ptr) < (double) (*ptr1);
		} else
		if(right.data_type == DOUBLE)
		{
			double *ptr1 = (double *) right.point;
			return (*ptr) < (*ptr1);
		}
	} else
	if(data_type == STRING && right.data_type == STRING)
	{
		String *ptr = (String *) point;
		String *ptr1 = (String *) right.point;
		return (*ptr) < (*ptr1);
	} else
	if(data_type == ARRAY && right.data_type == ARRAY)
	{
		Array<Mix> *ptr = (Array<Mix> *) point;
		Array<Mix> *ptr1 = (Array<Mix> *) right.point;
		return (*ptr) < (*ptr1);
	} else
	if(data_type == MAPPING && right.data_type == MAPPING)
	{
		Mapping<Mix, Mix> *ptr = (Mapping<Mix, Mix> *) point;
		Mapping<Mix, Mix> *ptr1 = (Mapping<Mix, Mix> *) right.point;
		return (*ptr) < (*ptr1);
	}

	return data_type < right.data_type;
}

bool Mix::operator>(const Mix &right) const
{
	return right < (*this);
}

bool Mix::operator<=(const Mix &right) const
{
	return !((*this) > right);
}

bool Mix::operator>=(const Mix &right) const
{
	return !((*this) < right);
}

bool Mix::operator!=(const Mix &right) const
{
	return ((*this) < right || right < (*this));
}

bool Mix::operator==(const Mix &right) const
{
	return !((*this) != right);
}

void Mix::release_data()
{
	if(data_type == INT)
	{
		int *ptr;
		ptr = (int *) point;
		point = NULL;
        if(ptr != NULL)
    		delete ptr;
		data_type = NONE;
	} else
	if(data_type == DOUBLE)
	{
		double *ptr;
		ptr = (double *) point;
		point = NULL;
        if(ptr != NULL)
		    delete ptr;
		data_type = NONE;
	} else
	if(data_type == STRING)
	{
		String *ptr;
		ptr = (String *) point;
		point = NULL;
		delete ptr;
		data_type = NONE;
	} else
	if(data_type == ARRAY)
	{
		Array<Mix> *ptr;
		ptr = (Array<Mix> *) point;
		point = NULL;
		delete ptr;
		data_type = NONE;
	} else
	if(data_type == MAPPING)
	{
		Mapping<Mix, Mix> *ptr;
		ptr = (Mapping<Mix, Mix> *) point;
		point = NULL;
		delete ptr;
		data_type = NONE;
	}
}

Mix::~Mix()
{
	if(data_type == INT)
	{
		int *ptr = (int *) point;
		point = NULL;
		delete ptr;
	} else
	if(data_type == DOUBLE)
	{
		double *ptr = (double *) point;
		point = NULL;
		delete ptr;
	} else
	if(data_type == STRING)
	{
		String *ptr = (String *) point;
		point = NULL;
		delete ptr;
	} else
	if(data_type == ARRAY)
	{
		Array<Mix> *ptr = (Array<Mix> *) point;
		point = NULL;
		delete ptr;
	} else
	if(data_type == MAPPING)
	{
		Mapping<Mix, Mix> *ptr = (Mapping<Mix, Mix> *) point;
		point = NULL;
		delete ptr;
	}
};

Mix::operator int() const
{
	assert(data_type == INT || data_type == DOUBLE);
	if(data_type == INT)
	{
		int *ptr = (int *) point;
		return *ptr;
	} else
	{
		double *ptr = (double *) point;
		return (int) *ptr;
	}
};

Mix::operator double() const
{
	assert(data_type == DOUBLE || data_type == INT);
	if(data_type == DOUBLE)
	{
		double *ptr = (double *) point;
		return *ptr;
	} else
	{
		int *ptr = (int *) point;
		return (double) *ptr;
	}
};

Mix::operator String() const
{
	assert(data_type == STRING);
	String *ptr = (String *) point;
	return *ptr;
};

Mix::operator char *() const
{
	assert(data_type == STRING);
	String *ptr = (String *) point;
	return (char *) *ptr;
}

Mix::operator Array<Mix>() const
{
	assert(data_type == ARRAY);
	Array<Mix> *ptr = (Array<Mix> *) point;
	return *ptr;
}

Mix::operator Mapping<Mix, Mix>() const
{
	assert(data_type == MAPPING);
	Mapping<Mix, Mix> *ptr = (Mapping<Mix, Mix> *) point;
	return *ptr;
}

int Mix::get_type()
{
	assert(data_type >=0 && data_type <= 5);

	return data_type;
}

#endif
