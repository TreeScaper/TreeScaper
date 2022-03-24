#pragma once

#include <iostream>
#include <fstream>
#include "wstring.hpp"
#include <map>

class Header_info {
public:
	Header_info() {};

	Header_info(std::fstream &fin);

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
private:
	std::map<String, String> list;

};

namespace Header {

	int end_header(std::fstream& fin);
	// Return the end position of header information (return 0 if no header information existed).

	void insert_header(std::fstream& fin, String** info, int lines);
	// Insert header information stored in info.

	bool check_header(std::fstream& fin, String content);
	// Return true if "content" presents in header informaion of the file.

	int load_header(std::fstream& fin, String** info);
	// Return lines of header information and store them in info.
}