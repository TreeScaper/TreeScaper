#include "header_info.hpp"

int Header::end_header(std::fstream& fin) {
	fin.seekg(0, fin.beg);
	int pos = 0;
	char c;
	fin >> c;
	if (c != '<')
		return 0;
	else {
		while (c != '>' && !fin.eof()) 
			fin.get(c);
		if (fin.eof()) {
			std::cout << "Error! Wrong header information format.\n";
			throw(1);
		}
		else
			pos = fin.tellp();
	}

	return pos;

}

void Header::insert_header(std::fstream& fin, String** info, int lines) {
	if (Header::end_header(fin) != 0) {
		std::cout << "Error! Attempt to insert header information to a file that already have one.\n";
		throw(1);
	}
	fin.seekg(0, fin.beg);
	fin << "<\n\n";
	for (int i = 0; i < lines; i++) 
		fin << info[i][0] << ":" << info[i][1] << '\n';
	fin << "\n>\n";

}

bool Header::check_header(std::fstream& fin, String content) {
	fin.seekg(0, fin.beg);
	char temp[1000];
	std::string stemp;
	std::string scontent((char*)content);
	while (temp[0] != '>' && !fin.eof()) {
		fin.getline(temp, 1000);
		stemp = std::string(temp);
		if (stemp.find(scontent) != std::string::npos) {
			fin.seekg(0, fin.beg);
			return true;
		}
	}
	fin.seekg(0, fin.beg);
	return false;
}

int Header::load_header(std::fstream& fin, String** info) {
	if (Header::end_header(fin) == 0)
		return 0;
	fin.seekg(0, fin.beg);
	int lines = 0;
	int pos = 0;
	char temp[1000];
	std::string stemp;
	while (temp[0] != '>' && fin.eof()) {
		fin.getline(temp, 1000);
		stemp = temp;
		if (stemp.length() > 2) {
			pos = stemp.find(':');
			info[lines][0] = (stemp.substr(0, pos - 1)).c_str();
			info[lines][1] = (stemp.substr(pos + 1, stemp.length() - 1)).c_str();
			lines++;
		}
	}
	return lines;
}

Header_info::Header_info(std::fstream& fin) {
	std::map<String, String> list;
	int header_pos = Header::end_header(fin);
	if (header_pos == 0) 
		//std::cout << "No header was found in " << input.get_filename() << ".\n";
		std::cout << "No header was found.\n";
	else {
		fin.seekg(0, fin.beg);
		char temp[1000];
		std::string stemp;
		String item;
		String content;
		int colon_pos;
		while (temp[0] != '>' && !fin.eof()) {
			fin.getline(temp, 1000);
			stemp = std::string(temp);
			colon_pos = stemp.find(':');
			if (colon_pos != std::string::npos) {
				item = (stemp.substr(0, colon_pos)).c_str();
				content = (stemp.substr(colon_pos + 1, stemp.length())).c_str();
				list[item] = content;
			}
		}
	}
}

Header_info::Header_info(String* item, String* content, int length) {
	for (int i = 0; i < length; i++) {
		(*this).list[item[i]] = content[i];
	}
}

std::istream &operator>>(std::istream &input, Header_info &info) {
	input.seekg(0, std::ios::beg);
	char temp[1000];
	input.getline(temp, 1000);
	if(temp[0]!='<')
		std::cout << "No header was found in the input file.\n";
	else {
		std::string stemp;
		String item;
		String content;
		int colon_pos;
		while (temp[0] != '>' && !input.eof()) {
			input.getline(temp, 1000);
			stemp = std::string(temp);
			colon_pos = stemp.find(':');
			if (colon_pos != std::string::npos) {
				item = (stemp.substr(0, colon_pos)).c_str();
				content = (stemp.substr(colon_pos + 1, stemp.length())).c_str();
				info.insert(item, content);
			}
		}
	}
	return input;

}

std::ostream &operator<<(std::ostream &output, Header_info &info){
	output << "<\n\n";
	if (info.size() != 0) {
		for (auto it = info.list.begin(); it != info.list.end(); it++)
			output << it->first << ':' << it->second << '\n';
	}
	output << "\n>\n";
	return output;
}


