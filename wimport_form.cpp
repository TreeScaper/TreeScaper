
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

// wimport_form.cpp
// member functions
//                by whuang Jan/04/2008

#include "wimport_form.h"

wImport_Form IMPORT_FORM;
wImport_Form *LPIMPORT_FORM = &IMPORT_FORM;

Array<Mapping<Mix, Mix> > wImport_Form::Read(String filename, Mapping<Mix, Mix> mdata)
{
	Array<Mapping<Mix, Mix> > arr, result;

    arr = Read(filename);

	for(int i = 0; i < arr.get_length(); i++)
		if((arr[i]).Include_Key((mdata.Keys())[0]) && arr[i][mdata.Keys()[0]] == mdata.Values()[0])
				result.add(arr[i]);

	return result;
}

Array<Mapping<Mix, Mix> > wImport_Form::Read(String filename)
{
    Array<Mapping<Mix, Mix> > result;
	Mix key, value;
	char wmsg[50];

	Array<String> lines;
	String line;

	Array<Array<Mix> > mix_lines;
	Array<Mix> mix_line;

	char element;

	fstream file(filename, ios::out | ios::in);

	if(!file)
	{
		sprintf(wmsg, "fail to open %s file", (char *) filename);
		cout << "{ ERROR! file : wimport_form.cpp / Read(String filename)\n"
			<< wmsg
			<< "\n}\n";
	}

	while(file.get(element))
	{
		if(element != '\n' && element != 13)
			line.add(element);

		if((element == '\n' || element == 13) && line != EMPTY_STRING)
		{
			lines.add(line);
			line = EMPTY_STRING;
		}
	}

	for(int i = 0; i < lines.get_length(); i++)
	{
		line = lines[i];
		if(! Modify_Line(lines[i], line))
		{
			cout << "{ ERROR! file : wimport_form.cpp/Read(String filename)\n"
				<< "error function : Modify_Line\n"
				<< "parameter1 : \n"
				<< lines[i]
				<< "\nparameter2 : \n"
				<< line
				<< "\n}\n";
		}

		line = Delete_Remark(line);
		if(line == EMPTY_STRING)
			continue;

		if(! String_Array(line, mix_line))
		{
			cout << "{ ERROR! file : wimport_form.cpp/Read(String filename)\n"
				<< "error function : String_Array\n"
				<< "parameter1 : \n"
				<< line
				<< "\nparameter2 : \n"
				<< mix_line
				<< "\n}\n";
		}
		mix_lines.add(mix_line);
	}
	// mix_lines to result mapping array
	if(mix_lines.get_length() > 0)
		result.resize(mix_lines.get_length() - 1);

	for(int j = 0; j < mix_lines.get_length() - 1; j++)
	{
		for(int i = 0; i < (mix_lines[0]).get_length(); i++)
		{
			key = mix_lines[0][i];
			result[j][key] = (mix_lines[j + 1][i]);
        }
    }
	return result;
}

bool wImport_Form::Modify_Line(const String &line, String &result)
{
	String modify_line;
	int quotation_num = 0;
	bool in_one_part = false;

	for(int i = 0; i < line.get_length(); i++)
	{
                if(line[i] == ',')   // if element equal to ','
                {
			if(! in_one_part)  // if this element is not a part
				modify_line.add(' ');  // convert ',' to space
			else
                                modify_line.add(line[i]);
                }

		while(i < line.get_length() && line[i] == '"')
		{
			quotation_num++;   // count the number of '"'
			i++;
		}

		if(quotation_num)
		{
			i--;  // go back to the nearest element which is '"'

			if(quotation_num % 2 == 1)
				in_one_part = ! in_one_part;

			for(int j = 0; j < (int) quotation_num / 2; j++)
				modify_line.add('"');

			quotation_num = 0;  // reset
		}

		if(! Filter(line[i]) && line[i] != ',' && line[i] != '"')
			modify_line.add(line[i]);
	}

	if(in_one_part)
	{
		cout << "{ ERROR! file : wimport_form.cpp / Modify_Line(const String &line, String &result)\n"
			<< "\" can't match\n"
			<< "}\n";
		return false;
	}

	result = modify_line;
	return true;
}

String wImport_Form::Delete_Remark(String line)
{
	for(int i = 0; i < line.get_length() - 1; i++)
		if(line[i] == '/' && line[i + 1] == '/')
			return line(0, i);
	return line;
}

bool wImport_Form::String_Array(const String &line, Array<Mix> &result)
{
	Array<Mix> arr;
	Mix m;

	String str, element, rest;
	str = line;
	do{
		String_Element(str, element, rest);
		if(! To_Mix(element, m))
		{
			cout << "{ ERROR! file : wimport_form.cpp/String_Array(const String &line, Array<Mix> &result)\n"
				<< "error function : To_Mix\n"
				<< "parameter1 : \n"
				<< element
				<< "\nparameter2 : \n"
				<< m
				<< "\n}\n";
			return false;
		}
		arr.add(m);
		str = rest;
	} while(rest != EMPTY_STRING);

	result = arr;
	return true;
}

bool wImport_Form::String_Element(String &str, String &element, String &rest)
{
	element = EMPTY_STRING;
	rest = EMPTY_STRING;
        bool in_one_part = false;
	int i = 0;
	int begin, end, rest_begin;
	bool is_str = false, is_arr = false, is_map = false;

	// looking for begin
	for(i = 0; i < str.get_length(); i++)
		if(str[i] != ' ')
		{
			if(str[i] == '"')
				is_str = true;
			if(str[i] == '{')
				is_arr = true;
			if(str[i] == '[')
				is_map = true;
			begin = i;
			break;
		}

	if(i == str.get_length())
		return true;

	// looking for end
	for(i = begin; i < str.get_length(); i++)
	{
		if(! is_str && ! is_arr && ! is_map && str[i] == ' ')
		{
			end = i;
			break;
		}
		if(str[i] == '}')
			is_arr = false;
		if(str[i] == '"')
			is_str = false;
		if(str[i] == ']')
			is_map = false;
	}

	if(i == str.get_length())
	{
		if(! is_arr && ! is_str && ! is_map)
		{
			end = str.get_length();
			element = str(begin, end - begin);
			return true;
		} else
		{
			cout << "{ ERROR! file : wimport_form.cpp / String_Element(String &str, String &element, String &rest)\n"
				<< "fail to generate regular string element\n"
				<< "}\n";
			return false;
		}
	}

	element = str(begin, end - begin);

	// looking for rest_begin
	for(i = end; i < str.get_length(); i++)
		if(str[i] != ' ')
		{
			rest_begin = i;
			break;
		}

	if(i == str.get_length())
		return true;

	rest = str(rest_begin, str.get_length() - rest_begin);
	return true;
}

bool wImport_Form::To_Mix(String &element, Mix &result)
{
	int length = element.get_length();

	// if element is string
	if(element[0] == '"' && element[length - 1] == '"')
	{
		result = element(1, length - 2);
		return true;            // get string
	}

	// if element is number
	Array<int> left, right;
	bool pass_radix_point = false;
	bool is_number = true;
	bool positive_number = true;

	if(element[0] == '-')
	{
		positive_number = false;
		element = element(1, length - 1);
		length = length - 1;
	}

	for(int i = 0; i < length; i++)
	{
    if(((element[i] < '0' || element[i] > '9') && element[i] != '.') || (element[i] == '.' && pass_radix_point))
		{
			is_number = false;  // element isn't a number
			break;
		}

		if(element[i] == '.')
		{
			pass_radix_point = true;
			continue;
		}

		if(! pass_radix_point)
			left.add(element[i] - '0');
		else
			right.add(element[i] - '0');
	}

	if(is_number && !pass_radix_point)
	{
		int integer_num = 0;
		for(int i = 0; i < left.get_length(); i++)
			integer_num = integer_num * 10 + left[i];

		if(positive_number)
			result = integer_num;
		else
			result = -integer_num;
		return true;                             // get int
	} else
	if(is_number && pass_radix_point)
	{
		double double_num = 0;
		for(int i = 0; i < left.get_length(); i++)
			double_num = double_num * 10 + left[i];

		for(int i = 0; i < right.get_length(); i++)
			double_num = double_num + right[i] * pow(0.1, i + 1);

		if(positive_number)
			result = double_num;
		else
			result = -double_num;
		return true;                                // get double
	}

	// if element is array
	Array<Mix> arr;
	String line, pre_element, rest;
	line = element;
	Mix m;

	if(line[0] == '{' && line[length - 1] == '}')
	{
		line = line(1, line.get_length() - 2);
		do{
			if(! Array_Element(line, pre_element, rest))
			{
				cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
					<< "error function : Array_Element\n"
					<< "parameter1 : \n"
					<< line
					<< "\nparameter2 : \n"
					<< pre_element
					<< "\nparameter3 : \n"
					<< rest
					<< "\n}\n";
				return false;
			}
			if(pre_element == EMPTY_STRING && rest != EMPTY_STRING)
			{
				cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
					<< "generate array : pre_element == EMPTY_STRING && rest != EMPTY_STRING \n"
					<< "}\n";
				return false;
			}

			if(pre_element != EMPTY_STRING)
			{
				if(! To_Mix(pre_element, m))
				{
					cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
						<< "error function : To_Mix\n"
						<< "parameter1 : \n"
						<< pre_element
						<< "\nparameter2 : \n"
						<< m
						<< "\n}\n";
					return false;
				}
				arr.add(m);

				line = rest;
			}
		}while(rest != EMPTY_STRING);

		result = arr;
		return true;
	}

	// if element is mapping
	Mapping<Mix, Mix> mapp;
	bool is_end_mapping = true;
		//  All have declare before
		//	String line, pre_element, rest;
		//	line = element;
	Mix key, value;

	if(line[0] == '[' && line[length - 1] == ']')
	{
		line = line(1, line.get_length() - 2);

		do{
			if(! Mapping_Element(line, pre_element, rest, is_end_mapping))
			{
				cout << "{ ERROR! file : wimport_form.cpp / To_Mix(String &element, Mix &result)\n"
					<< " error function : Mapping_Element\n"
					<< "parameter1 : \n"
					<< line
					<< "\nparameter2 : \n"
					<< pre_element
					<< "\nparameter3 : \n"
					<< rest
					<< "\nparameter4 : \n"
					<< is_end_mapping
					<< "\n}\n";
				return false;
			}

			if(pre_element == EMPTY_STRING && rest != EMPTY_STRING)	
			{
				;
				cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
					<< "generate mapping : pre_element == EMPTY_STRING && rest != EMPTY_STRING \n"
					<< "}\n";
				return false;
			}

			if(pre_element != EMPTY_STRING)
			{
				if(! To_Mix(pre_element, m))
				{
					;
					cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
						<< "error function : To_Mix\n"
						<< "parameter1 : \n"
						<< pre_element
						<< "\nparameter2 : \n"
						<< m
						<< "\n}\n";
					return false;
				}
				
				if(is_end_mapping)
					value = m;
				else
					key = m;

				if(is_end_mapping)
					mapp[key] = value;

				line = rest;
			}

			if(rest == EMPTY_STRING && is_end_mapping == false)
			{
				;
				cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
					<< "only key found\n"
					<< "}\n";
				return false;
			}

		} while(rest != EMPTY_STRING);

		result = mapp;
		return true;
	}

	cout << "{ ERROR! file : wimport_form.cpp/To_Mix(String &element, Mix &result)\n"
		<< "error string form, can't be to Mix\n"
		<< "}\n";
	return false;
}

bool wImport_Form::Array_Element(String &line, String &element, String &rest)
{
	element = EMPTY_STRING;
	rest = EMPTY_STRING;
	bool is_in_array = false, is_in_mapping = false;
	int begin = -1, end = -1, rest_begin = -1;
	int i = 0;

	// looking for begin
	for(i = 0; i < line.get_length(); i++)
	{
		if(line[i] == ' ')
			continue;
		else
		{
			if(line[i] == '{')
				is_in_array = true;
			if(line[i] == '[')
				is_in_mapping = true;
			begin = i;
			break;
		}
	}

	// line only contain space
	if(i == line.get_length())
		return true;

	// looking for end
	for(i = begin;  i < line.get_length(); i++)
	{
		if(is_in_array)
		{
			if(line[i] != '}')
				continue;
			else
			{
				i++;
				end = i;
				is_in_array = false;
				break;
			}
		} else
		if(is_in_mapping)
		{
			if(line[i] != ']')
				continue;
			else
			{
				i++;
				end = i;
				is_in_mapping = false;
				break;
			}
		} else
		{
			if(line[i] != ' ' && line[i] != ',')
			{
				end = -1;
				continue;
			}
			if(line[i] == ' ' && end == -1)
				end = i;
			if(line[i] == ',')
			{
				if(end != -1)
					break;
				else
				{
					end = i;
					break;
				}
			}
		}
	}

	if(is_in_array && i == line.get_length())
	{
		cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
			<< "{} can't match\n"
			<< "}\n";
		return false;
	}

	if(is_in_mapping && i == line.get_length())
	{
		cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
			<< "[] can't match\n"
			<< "}\n";
		return false;
	}

	if(end == -1 && i == line.get_length())
		end = line.get_length();

	element = line(begin, end - begin);
	if(i == line.get_length())
		return true;

	// looking for rest_begin
	for(i = end + 1; i < line.get_length(); i++)
		if(line[i] != ' ' && line[i] != ',')
		{
			rest_begin = i;
			rest = line(i, line.get_length() - i);
			return true;
		}

	// otherwise false
	cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
		<< "wrong array form\n"
		<< "}\n";
	return false;
}

bool wImport_Form::Mapping_Element(String &line, String &element, String &rest, bool &is_end_mapping)
{
	element = EMPTY_STRING;
	rest = EMPTY_STRING;
	bool is_in_array = false, is_in_mapping = false;
	int begin = -1, end = -1, rest_begin = -1;
	int i = 0;
	char end_char;

	// looking for begin
	for(i = 0; i < line.get_length(); i++)
	{
		if(line[i] == ' ')
			continue;
		else
		{
			if(line[i] == '{')
				is_in_array = true;
			if(line[i] == '[')
				is_in_mapping = true;
			begin = i;
			break;
		}
	}

	// line only contain space
	if(i == line.get_length())
	{
		is_end_mapping = ! is_end_mapping;
		return true;
	}

	// looking for end
	for(i = begin;  i < line.get_length(); i++)
	{
		if(is_in_array)
		{
			if(line[i] != '}')
				continue;
			else
			{
				i++;
				end = i;
				is_in_array = false;
				break;
			}
		} else
		if(is_in_mapping)
		{
			if(line[i] != ']')
				continue;
			else
			{
				i++;
				end = i;
				is_in_mapping = false;
				break;
			}
		} else
		{
			if(is_end_mapping)
				end_char = ':';
			else
				end_char = ',';
			if(line[i] != ' ' && line[i] != end_char)
			{
				end = -1;
				continue;
			}
			if(line[i] == ' ' && end == -1)
				end = i;
			if(line[i] == end_char)
			{
				if(end != -1)
					break;
				else
				{
					end = i;
					break;
				}
			}
		}
	}

	if(is_in_array && i == line.get_length())
	{
		cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
			<< "{} can't match\n"
			<< "}\n";
		return false;
	}

	if(is_in_mapping && i == line.get_length())
	{
		cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
			<< "[] can't match\n"
			<< "}\n";
		return false;
	}

	if(end == -1 && i == line.get_length())
		end = line.get_length();

	element = line(begin, end - begin);
	if(i == line.get_length())
	{
		is_end_mapping = ! is_end_mapping;
		return true;
	}

	// looking for rest_begin
	for(i = end + 1; i < line.get_length(); i++)
	{
		if(is_end_mapping)
			end_char = ':';
		else
			end_char = ',';

		if(line[i] != ' ' && line[i] != end_char)
		{
			rest_begin = i;
			rest = line(i, line.get_length() - i);
			is_end_mapping = ! is_end_mapping;
			return true;
		}
	}

	// otherwise false
	cout << "{ ERROR! file : wimport_form.cpp/Array_Element(String &line, String &element, String &rest)\n"
		<< "wrong mapping form\n"
		<< "}\n";
	return false;
}

bool wImport_Form::Filter(char letter)
{
	return false;
}
