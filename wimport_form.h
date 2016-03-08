
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

// wimport_form.h
//      import a form
//               by whuang

#ifndef WIMPORT_FORM_H
#define WIMPORT_FORM_H

#include <vector>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "warray.cpp"
#include "wstring.h"
#include "wmix.h"
#include "wmapping.cpp"

class wImport_Form{
public:
	Array<Mapping<Mix, Mix> > Read(String filename, Mapping<Mix, Mix> mdata);

    Array<Mapping<Mix, Mix> > Read(String filename);

	bool Modify_Line(const String &line, String &result);

	String Delete_Remark(String line);

	bool String_Array(const String &line, Array<Mix> &result);

	bool String_Element(String &str, String &element, String &rest);

	bool To_Mix(String &element, Mix &result);

	bool Array_Element(String &line, String &element, String &rest);

	bool Mapping_Element(String &line, String &element, String &rest, bool &is_end_mapping);

	bool Filter(char letter);               //filter useless letter

private:
};

extern wImport_Form IMPORT_FORM;
extern wImport_Form *LPIMPORT_FORM;

#endif
