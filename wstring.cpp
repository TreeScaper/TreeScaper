
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

// wstring.cpp
// Member function definitions for class String
//              by whuang Nov/03/2008

#ifndef WSTRING_CPP
#define WSTRING_CPP

#include <iomanip>
#include <iostream>
#include <string>
#include "wstring.h"
#include <cassert>

String EMPTY_STRING;

// constructor: convert char * to String
String::String(const char *s) : length(strlen(s))
{
    set_String(s);            // call utility function
}

// copy constructor
String::String(const String &copy) : length(copy.length)
{
	set_String(copy.str_ptr);           // call utility function
}

// destructor
String::~String()
{
	delete [] str_ptr;                  // reclaim String
}

// = operator, avoids self assignment
const String &String::operator=(const String &right)
{
	if(&right != this)                 // avoid self assignment
	{
        delete [] str_ptr;             // prevents memory leak
        length = right.length;         // new String length
        set_String( right.str_ptr);    // call utility function
	}
	else
        std::cout << "warning: attempted assignment of a String to itself";

	return *this;                      // enable cascaded assignments
}

// concatenate right operand to this object and
// store in this object.
const String &String::operator+=(const String &right)
{
	char *temp_ptr = str_ptr;        // hold to be able to delete
	length += right.length;          // new String length
	str_ptr = new char[length + 1];  // create space
	assert(str_ptr != 0);            // terminate if memory not allocated
	strcpy(str_ptr, temp_ptr);       // left part of new String
	strcat(str_ptr, right.str_ptr);  // right part of new String
	delete [] temp_ptr;              // reclaim old space
	return *this;                    // enable cascaded calls
}

// Is this String empty?
bool String::is_empty() const
{
	return length == 0;
}

// Is this String equal to right String?
bool String::operator==(const String &right) const
{
	return strcmp(str_ptr, right.str_ptr) == 0;
}

// Is this String less than right String?
bool String::operator<(const String &right) const
{
	return strcmp(str_ptr, right.str_ptr) < 0;
}

// return a reference to a character in a String as an value.
char &String::operator[](int subscript)
{
	// first test for subscript out of range
    assert(subscript >= 0 && subscript <= length);

	return str_ptr[subscript];            // creates value
}

// return a reference to a character in a String as an value.
const char &String::operator[](int subscript) const
{
	// first test for subscript out of range
	assert(subscript >= 0 && subscript < length);

	return str_ptr[subscript];            // creates value
}

// return a subString beginning at index and
// of length sublength as a reference to a String object.

// zd_coded: abandoned the interpretation of the parameters
// so the second interger always represents the length of the substring.
// - 5/21/20
String String::operator()(int index, int sub_length)
{
	//ensure index is in range and subString length >= 0
	assert( index >= 0 && index < length && sub_length >=0 );
	String str;
	String *sub_ptr = &str;       //empty String
	assert(sub_ptr != 0);               // ensure new String allocated

	// determin length of subString
	//if( ( sub_length == 0) || ( index + sub_length > length ) )
		//sub_ptr->length = length - index + 1;
	//else
		//sub_ptr->length = sub_length;

  assert((sub_length != 0) && (index + sub_length <= length + 1));
  // throw when the sub_string went beyond the string
  sub_ptr->length = sub_length;
	// allocate memory for subString
	delete [] sub_ptr->str_ptr;         // delete character array from object
	sub_ptr->str_ptr = new char[sub_ptr->length + 1];
	assert(sub_ptr->str_ptr != 0);      // ensure space allocated

	// copy subString into new String
	strncpy(sub_ptr->str_ptr, &str_ptr[index], sub_ptr->length);
	sub_ptr->str_ptr[sub_ptr->length] = '\0';       // terminate String

	return *sub_ptr;             // return new String
}

// return a substring before the input char
String String::before(char s)                        // return a subString
{
    for(int i = 0; i < this->get_length(); i++)
    {
        if((*this)[i] == s)
        {
            return (*this)(0, i);
        }
    }
    return (*this);
}


int String::get_length() const
{
	return length;
}

void String::add(char element)
{
	char temp[] = " ";
	temp[0] = element;
	(*this) += temp;
}

// utility function to be called by constructors and
// assignment operator.
void String::set_String( const char *String2)
{
	str_ptr = new char[length + 1];    // allocate storage
	assert( str_ptr != 0);             // terminate if memory not allocated
	strcpy( str_ptr, String2);         // copy literal to object
}

#endif
