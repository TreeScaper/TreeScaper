
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

#ifndef WFILE_CPP
#define WFILE_CPP

#include "wfile.h"

void File::clear()
{
	if(fhandle.eof())
		fhandle.clear();
};

bool File::open(String name)
{
	fname = name;
	fhandle.open(name);

	if(fhandle.is_open())
		return true;
	return false;
};

void File::close()
{
	if(fhandle.is_open())
		fhandle.close();
}

bool File::is_open()
{
	return fhandle.is_open();
};

bool File::is_end()
{
	return fhandle.eof();
};

void File::seek(int n)
{
	if(fhandle.eof())
		fhandle.clear();
    fhandle.seekg(n, std::ios_base::beg);
};

int File::lines()
{
	int n = 0;
	char str[100000];
	while(fhandle.getline(str, 100000))
	{
		if(strlen(str) != 0)
			n = n + 1;
	}

	return n;
};

int File::cols()
{
	int n = 0;
	char str[100000] = "";
	bool pre_is_table = true;
	fhandle.getline(str, 100000);
	for(int i = 0; i < strlen(str); i++)
	{
		if(str[i] == '\t')
		{
			if(!pre_is_table)
				n++;
			pre_is_table = true;
		} else
		{
			pre_is_table = false;
		}
	}
	if(!pre_is_table)
		n++;
	return n;
};

bool File::clean()
{
	(*this).close();
    std::ofstream fout(fname);
	fout.close();
	(*this).open(fname);
	(*this).seek(0);

	return true;
}

String File::prefix_name()
{
    int i = 0;
    while(fname[i] != '.' && fname[i] != '\0')
    {
        i++;
    }

    return fname(0, i);
};

String File::postfix_name()
{
	int i = 0;
	while(fname[i] != '.' && fname[i] != '\0')
		i++;

	if(i == fname.get_length())
		return "";

	return fname(i + 1, fname.get_length() - i - 1);
};

String File::prefix_name_lastof()
{
    int i = fname.get_length();
    while(fname[i] != '.')
        i--;

    return fname(0, i);
};

String File::postfix_name_lastof()
{
    int i = fname.get_length();
    while(fname[i] != '.')
        i--;

    if(i == 0)
        return "";

    return fname(i + 1, fname.get_length() - i - 1);
};

#endif
