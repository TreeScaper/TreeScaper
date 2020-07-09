
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

	String prefix_name();

	String postfix_name();

    String prefix_name_lastof();

    String postfix_name_lastof();

private:
    std::fstream fhandle;
	String fname;
};

#endif
