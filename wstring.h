
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

// wstring.h
// Definition of a String class
//           by whuang Nov/03/2008

#ifndef WSTRING_H
#define WSTRING_H

#include <iomanip>
#include <cstring>
#include <iostream>

class String {
    friend std::istream &operator>>(std::istream &input, String &s)
	{
        char temp[1000];        // buffer to store input

        input >> std::setw(1000) >> temp;
        s = temp;                 // use String class assignment operator
		return input;             // enables cascading
	};

    friend std::ostream &operator<<(std::ostream &output, const String &s)
	{
		output << s.str_ptr;
		return output;        // enables cascading
	};

public:
    String( const char * s = "");                 // constuct and default empty String
	String( const String &);                    // copy constructor
    ~String();                                  // destructor
	const String &operator=(const String &);    // assignment
	const String &operator+=(const String &);   // concatenation

	bool is_empty() const;                      // is String empty
	bool operator==(const String &) const;      // test s1 ?= s2
	bool operator<(const String &) const;       // test s1 ?< s2

	// test s1 ?!= s2
	bool operator!=(const String & right) const{return !(*this == right);};

	// test s1 ?> s2
	bool operator>(const String &right) const{return (right < *this);};

	// test s1 ?<= s2
	bool operator<=(const String &right) const{return !(right < *this);};

	// test s1 ?>= s2
	bool operator>=(const String &right) const{return !(*this < right);};

	char &operator[](int);						// subscript operator
    const char &operator[](int) const;			// subscript operator
    String operator()(int, int);				// return a subString
    String before(char);                        // return a subString
	operator char *() const{return str_ptr;};   // translate to array.
	int get_length() const;						// return String length
	void add(char);

private:
	int length;                             // String length
	char *str_ptr;                          // pointer to start of String

	void set_String(const char *);          // utility function
};

extern String EMPTY_STRING;

#endif

