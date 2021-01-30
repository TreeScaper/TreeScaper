
#include <iostream>
#include <fstream>
#include <string>
#include "zdtree.hpp"

using std::string;



char read_sequence(string& tree, size_t& index, string& temp) {
	size_t found = tree.find_first_of(",): ", index);
	if (found == std::string::npos) {
		temp = tree.substr(index);
		index = tree.back();
		return '\0';
	}
	else {
		temp = tree.substr(index, found - index);
		index = found;
		return tree[index];
	}
	
	
}

void TaxonList::push(string item, bool flag_str_format) { // Assume to have format "Index Taxon" separated by space.
	int i = 0;
	while (isspace(item[i])) {
		i++;
	}

	//Skip space.

	if (i == item.size()) {
		return;
	}
	else if (flag_str_format){
		size_t found = item.find(' ', i);
		
		if (found != std::string::npos) {
			size_t end_pos = item.find_first_of(",;.\n\t ", found + 1);
			string taxon = item.substr(found + 1, end_pos - found - 1);
			Ind2Taxon.push(taxon);
			Taxon2Ind.insert(std::pair<string, int>(taxon, Ind2Taxon.get_size() - 1));
			size++;
		}
	}
	else{
		Ind2Taxon.push(item);
		Taxon2Ind.insert(std::pair<string, int>(item, Ind2Taxon.get_size() - 1));
		size++;
	}
}



size_t TaxonList::ReadTaxa(std::string fname) {
	std::ifstream fin;
	fin.open(fname);
	if (!fin.is_open()){
		std::cout << "File `" << fname << "' not opened.\n";
		throw(1);
	}

	size_t pos = std::ios_base::beg;
	bool flag_internal = false;

	string temp;
	std::getline(fin, temp);
	if (temp[0] != '#' || temp.find(string("NEWICK")) != std::string::npos) {
		fin.seekg(std::ios_base::beg);
		while (temp.find_first_of("(") == std::string::npos) {
			pos = fin.tellg();
			std::getline(fin, temp);
		}
		string taxon;

		size_t i = temp.find_first_of('(');
		char ch = temp[i];
		while (ch != ';' && ch != '\0') {
			taxon.clear();
			if (isspace(ch) || ch == ',' || ch == '(' || ch == ')'){	//skip spaces and commas and parentheses
				if (ch == ')')
					flag_internal = true;
				else if (ch == '(' || ch == ',')
					flag_internal = false;
				ch = temp[++i];											
			}
			else if (ch == ':') {										//skip weights
				i++;
				ch = read_sequence(temp, i, taxon);
				flag_internal = false;
			}
			else {														//read taxa
				ch = read_sequence(temp, i, taxon);
				if (!flag_internal) {
					this->push(taxon, false);
				}
				flag_internal = false;
			}
		}
	}
	else if (temp.find(string("NEXUS")) != std::string::npos) {
		while (!fin.eof()) {
			temp.clear();
			std::getline(fin, temp);
			if (temp.find(string("Translate")) != std::string::npos) {
				break;
			}
		}
		if (fin.eof()) {
			std::cout << "Error! Cannot find taxon list in Nexus format.\n";
			throw(1);
		}
		while (!fin.eof()) {
			temp.clear();
			std::getline(fin, temp);
			this->push(temp);
			if (temp.find(';') != std::string::npos)
				break;
		}
		pos = fin.tellg();
	}
	fin.close();

	return pos;
}


TreeArray::TreeArray(string& tree, TaxonList& taxon_list, Array<int> &active_levels, 
	Array<int> &unlabeled, int flag_label) {
	int leaf_size = taxon_list.size;
	bool flag_internal = false;
	char ch;
	levels = Array<Array<int> >(0, 2 * leaf_size - 3);
	Array<Array<double> > w_temp = Array<Array<double> >(0, 2 * leaf_size - 3);

	active_levels.erase();
	unlabeled.erase();

	int current_ind = -1;
	string temp;

	//Build dummy level.

	levels.push(Array<int>(0, 10));		
	unlabeled.push(0);
	levels.back().push(current_ind);	//parent level.
	levels.back().push(-1);				//parent node.
	
	w_temp.push(Array<double>(0, 10));  //build the similar array for weights.
										//Note that weights array do not have the first two columns for locating parent node.

	levels.back().push(-1);				//dummy root.
	unlabeled.back()--;
	current_ind++;

	//Build first level.

	levels.push(Array<int>(0, 10));							
	unlabeled.push(0);
	levels.back().push(current_ind);						//parent level.
	levels.back().push(levels[current_ind].get_size() - 1);		//parent node.

	w_temp.push(Array<double>(0, 10));

	current_ind++;

	//Start scanning/


	size_t i = tree.find_first_of('(') + 1;
	ch = tree[i];

	while(ch != ';' && ch != '\0'){
		temp.clear();
		if (isspace(ch))	
			ch = tree[++i];											//skip spaces
		else if (ch == ',') {
			ch = tree[++i];
			flag_internal = false;
		}
		else if (ch == '(') {
			levels[current_ind].push(-1);
			unlabeled[current_ind]--;								//accumulate unlabeled point.
			levels.push(Array<int>(0, 10));							//add level.
			w_temp.push(Array<double>(0, 10));
			unlabeled.push(0);
			levels.back().push(current_ind);						//parent level.
			levels.back().push(levels[current_ind].get_size() - 1);		//parent node.
			
			current_ind = levels.get_size() - 1;						//go to the new level in the back.
			ch = tree[++i];
			flag_internal = false;
		}
		else if (ch == ')') {
			current_ind = levels[current_ind].front();
			if (current_ind < 0)
				break;
			ch = tree[++i];

			flag_internal = true;
		}
		else {
			if (ch == ':') {										
				i++;
				ch = read_sequence(tree, i, temp);
				w_temp[current_ind].push(atof(temp.c_str()));
				//weights.push(atof(temp.get_vec()));

				flag_internal = false;
			}
			else {
				ch = read_sequence(tree, i, temp);
				if (!flag_internal) {
					if (flag_label == 0)			//labeled with integer
						levels[current_ind].push(atoi(temp.c_str()) - 1);
					else if (flag_label == 1)		//labeled with different normalization
						levels[current_ind].push(taxon_list.IndB2IndA[atoi(temp.c_str())]);
					else if (flag_label == 2) {		// labeled with taxon
						string taxon(temp.c_str());
						auto it = taxon_list.Taxon2Ind.find(taxon);
						if (it != taxon_list.Taxon2Ind.end())
							levels[current_ind].push(it->second);
						else {
							std::cout << "Error: Attempt to find taxon: ``" << taxon << "'' that is not in the leaf set.\n";
							throw(1);
						}
					}
					else {
						std::cout << "Error! Wrong flag for leaf's label.\n";
						throw(1);
					}

					active_levels.push(current_ind);
				}
				flag_internal = false;
			}
			
		}
	}

	
	size = 0;

	for (int i = 1; i < levels.get_size(); i++)
		size += levels[i].get_size() - 2;

	//When first level has only 2 nodes, their parent node can be compressed.
	//Therefore the first level becomes 2 points connected by 1 edge.
	if (levels[1].get_size() == 4)
		size--;

	edges = new int* [2];
	edges[0] = new int[size];	// node far-from-leaf
	edges[1] = new int[size];	// node close-to-leaf



	if (w_temp.get_size()==0)
		this->label_internal_node(active_levels, unlabeled);
	else
		this->label_internal_node(active_levels, unlabeled, w_temp);


	for (int i = 0; i < levels.get_size(); i++)
		levels[i].release();
	levels.release();

}

void TreeArray::label_internal_node(Array<int>& active_levels, Array<int>& unlabeled) {
	int node_index = active_levels.get_size();
	int edge_index = 0;
	int cur_level, parent_level;
	// start labels of internal node from the number of leaves, N.
	// Note that leaves are labeled by 0 to N - 1.

	for (int i = 0; active_levels[i] != 0; i++) {
		cur_level = active_levels[i];
		parent_level = levels[cur_level][0];
		if (unlabeled[cur_level] == 0) {											// all nodes on this level is labeled
			for (int j = 2; j < levels[cur_level].get_size(); j++) {					// add edge
				edges[0][edge_index] = node_index;
				edges[1][edge_index++] = levels[cur_level][j];
			}
			levels[parent_level][levels[cur_level][1]] = node_index++;				// label the parent node.
			unlabeled[parent_level]++;
			unlabeled[cur_level]++;													// mark current level done labeling.
			active_levels.push(parent_level);										// make parent level active.
		}
	}

	if (levels[1].get_size() == 4)														//	merge dummy node,
		edges[0][size - 1] = (edges[1][size - 1] == levels[1][2] ? levels[1][3] : levels[1][2]);
}

void TreeArray::label_internal_node(Array<int>& active_levels, Array<int>& unlabeled, Array<Array<double> > &w_temp) {
	int node_index = active_levels.get_size(); 
	int edge_index = 0;
	int cur_level, parent_level;
	// start labels of internal node from the number of leaves, N.
	// Note that leaves are labeled by 0 to N - 1.
	
	for (int i = 0; active_levels[i] != 0; i++) {
		cur_level = active_levels[i];
		parent_level = levels[cur_level][0];
		if (unlabeled[cur_level] == 0) {											// all nodes on this level is labeled
			for (int j = 2; j < levels[cur_level].get_size(); j++) {					// add edge and weight
				edges[0][edge_index] = node_index;									// far-from-leaf node (parent)
				edges[1][edge_index++] = levels[cur_level][j];						// close-to-leaf node (child)
				weights.push(w_temp[cur_level][j - 2]);								// add weights attached to the close-to-leaf node.
			}
			levels[parent_level][levels[cur_level][1]] = node_index++;				// label the parent node.
			unlabeled[parent_level]++;
			unlabeled[cur_level]++;													// mark current level done labeling.
			active_levels.push(parent_level);										// make the parent level active.
		}
	}

	if (levels[1].get_size() == 4)														//	merge dummy node,
		edges[0][size - 1] = (edges[1][size - 1] == levels[1][2] ? levels[1][3] : levels[1][2]);
}

//template <class T>
//int Bipartition<T>::push(const BitString<T>& bs, int index_hash) {
//	if (index_hash == -1)
//		index_hash = hashing(bs);
//	int index = -1;
//	if (Hash2Id[index_hash].is_empty()) {
//		Hash2Id[index_hash].push(Id2BitString.size());
//		Id2BitString.push(bs);
//		return Id2BitString.size() - 1;
//	}
//	else {
//		for (int i = 0; i < Hash2Id[index_hash].size(); i++) {
//			index = Hash2Id[index_hash][i];
//			if (Id2BitString[index] == bs)
//				return index;
//		}
//		Hash2Id[index_hash].push(Id2BitString.size());
//		Id2BitString.push(bs);
//		return Id2BitString.size() - 1;
//	}
//}

//template <class T>
//void TreeSet<T>::ReadTree(std::string fname, size_t pos) {
//	Array<int> active_levels(0, 100);
//	Array<int> unlabeled(0, 100);
//
//	std::ifstream fin;
//	fin.open(fname);
//	string temp;
//	
//	fin.seekg(pos);
//
//	while (!fin.eof()) {
//		temp.clear();
//		std::getline(fin, temp, ';');
//		if (temp.find('(') != std::string::npos){
//			this->push(new TreeArray(temp, *Taxa, active_levels, unlabeled, 2));
//		}
//	}
//
//	Bipart->set_hash_bound(this->size);
//	Bipart->set_invariant();
//}

//template <class T>
//void TreeSet<T>::compute_Bipart() {
//	if (Bipart->is_empty()) {
//		for (int i = 0; i < Taxa->size; i++)
//			Bipart->push(BitString<T>(Taxa->bitstr_size, i).normalized(Taxa->size));
//	}
//	for (int i = 0; i < this->size; i++) {
//		this->trees_newick[i]->compute_bitstring(*(this->Bipart));
//		this->trees_newick[i]->print_bipart();
//	}
//}
