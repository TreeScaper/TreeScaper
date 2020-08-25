
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

// wfile.h
//    file operation
//          March/11/2010
//                    by whuang

#ifndef WFILE_H
#define WFILE_H

#include <iostream>
#include <fstream>
#include "wstring.h"
#include <map>

class File{
	template<class T>
	friend File &operator>>(File &input, T &right)                  // overload >>
	{
		input.fhandle >> right;
		return input;                                                          // 
	};

    template<class T>
    friend std::fstream &operator<<(File &output, const T &right)           // overload <<
    {
        output.fhandle << right;
        return output.fhandle;                                                         // output form, for example {length: a1, a2, a3, }
    };

public:
	File(){};

	File(String name)
	{
        fname = name;
        fhandle.open(name);
	}

	~File()
	{
		//fname.~String();
		if(fhandle.is_open())
			fhandle.close();
	};

	void clear();

	bool open(String name);

	void close();

	bool is_open();

	bool is_end();

	void seek(int);

	int lines();

	int cols();

	bool clean();

	String get_filename() { return fname; };

	void getline(char* s, int n) {
		fhandle.getline(s, n);
	};

	String prefix_name();

	String postfix_name();

    String prefix_name_lastof();

    String postfix_name_lastof();

	int end_header();
	// Return the end position of header information (return 0 if no header information existed).

	void insert_header(String ** info, int lines);
	// Insert header information stored in info.

	bool check_header(String content);
	// Return true if "content" presents in header informaion of the file.

	int load_header(String** info);
	// Return lines of header information and store them in info.
	
private:
    std::fstream fhandle;
	String fname;
};


class Header_info {
public:
	Header_info() {};

	Header_info(File &input);

	Header_info(String* item, String* content, int length);

	~Header_info() {};

	const Header_info &operator=(const Header_info &rhs) {
		(*this).list = rhs.list;
		return (*this);
	}

	int size() { return list.size(); };

	void insert(String item, String content) {
		list[item] = content;
	};

	bool count(String rhs) {
		return list.count(rhs);
	}

	String operator[](String key) {
		return list[key];
	}

	friend std::istream &operator>>(std::istream &input, Header_info &info);

	friend std::ostream &operator<<(std::ostream &output, Header_info &info);

	friend File &operator>>(File &input, Header_info &info);

	friend File &operator<<(File &output, Header_info &info);
private:
	std::map<String, String> list;

};
#endif
