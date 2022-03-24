
#include <iostream>
#include <fstream>
#include <string>
#include "zdtree.hpp"

using std::string;

char read_sequence(string &tree, size_t &index, string &temp)
{
	size_t found = tree.find_first_of(",): ", index);
	if (found == std::string::npos)
	{
		temp = tree.substr(index);
		index = tree.back();
		return '\0';
	}
	else
	{
		temp = tree.substr(index, found - index);
		index = found;
		return tree[index];
	}
}

void TaxonList::push(string item, bool flag_str_format)
{ // Assume to have format "Index Taxon" separated by space.
	int i = 0;
	while (isspace(item[i]))
	{
		i++;
	}

	//Skip space.

	if (i == item.size())
	{
		return;
	}
	else if (flag_str_format)
	{
		size_t found = item.find(' ', i);

		if (found != std::string::npos)
		{
			size_t end_pos = item.find_first_of(",;.\n\t ", found + 1);
			string taxon = item.substr(found + 1, end_pos - found - 1);
			Ind2Taxon.push(taxon);
			Taxon2Ind.insert(std::pair<string, int>(taxon, Ind2Taxon.get_size() - 1));
			size++;
		}
	}
	else
	{
		Ind2Taxon.push(item);
		Taxon2Ind.insert(std::pair<string, int>(item, Ind2Taxon.get_size() - 1));
		size++;
	}
}

size_t TaxonList::ReadTaxa(std::string fname)
{
	std::ifstream fin;
	fin.open(fname);
	if (!fin.is_open())
	{
		std::cout << "File `" << fname << "' not opened.\n";
		throw(1);
	}

	size_t pos = std::ios_base::beg;
	bool flag_internal = false;

	string temp;
	std::getline(fin, temp);
	if (temp[0] != '#' || temp.find(string("NEWICK")) != std::string::npos)
	{
		fin.seekg(std::ios_base::beg);
		while (temp.find_first_of("(") == std::string::npos)
		{
			pos = fin.tellg();
			std::getline(fin, temp);
		}
		string taxon;

		size_t i = temp.find_first_of('(');
		char ch = temp[i];
		while (ch != ';' && ch != '\0')
		{
			taxon.clear();
			if (isspace(ch) || ch == ',' || ch == '(' || ch == ')')
			{ //skip spaces and commas and parentheses
				if (ch == ')')
					flag_internal = true;
				else if (ch == '(' || ch == ',')
					flag_internal = false;
				ch = temp[++i];
			}
			else if (ch == ':')
			{ //skip weights
				i++;
				ch = read_sequence(temp, i, taxon);
				flag_internal = false;
			}
			else
			{ //read taxa
				ch = read_sequence(temp, i, taxon);
				if (!flag_internal)
				{
					this->push(taxon, false);
				}
				flag_internal = false;
			}
		}
	}
	else if (temp.find(string("NEXUS")) != std::string::npos)
	{
		while (!fin.eof())
		{
			temp.clear();
			std::getline(fin, temp);
			if (temp.find(string("Translate")) != std::string::npos)
			{
				break;
			}
		}
		if (fin.eof())
		{
			std::cout << "Error! Cannot find taxon list in Nexus format.\n";
			throw(1);
		}
		while (!fin.eof())
		{
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

//TreeArray::TreeArray(string& tree, TaxonList& taxon_list, Array<int>& active_levels,
//	Array<int>& unlabeled, Array2D<double>& w_temp,
//	int flag_label, bool ISROOTED) {
//	issorted = false;
//	int leaf_size = taxon_list.size;
//	bool flag_internal = false;
//	char ch;
//
//	levels = new Array2D<int>(2 * leaf_size - 3);
//	w_temp.clean();
//	active_levels.clean();
//	unlabeled.clean();
//
//	int current_ind = -1;
//	string temp;
//
//	//darray<int>* levels_back_ptr;
//	//darray<int>* levels_cur_ptr;
//
//	//Build dummy level.
//
//	levels->push_c(new Array<int>(0, 10));
//	unlabeled.push(0);
//	//levels_back_ptr = levels.back_container();
//	levels->push(current_ind);			//parent level.
//	levels->push(-1);					//parent node.
//
//	w_temp.push_c(new Array<double>(0, 10));			//build the similar array for weights.
//														//Note that weights array do not have the first two columns for locating parent node.
//
//	levels->push(-1);					//dummy root.
//	unlabeled.back()--;
//	current_ind++;
//
//	//Build first level.
//
//	levels->push_c(new Array<int>(0, 10));
//	unlabeled.push(0);
//	levels->push(current_ind);								//parent level.
//	levels->push(levels->get_c_ptr(current_ind)->get_size() - 1);		//parent node.
//
//	w_temp.push_c(new Array<double>(0, 10));
//
//	current_ind++;
//
//	//Start scanning/
//
//
//	size_t i = tree.find_first_of('(') + 1;
//	ch = tree[i];
//
//	while (ch != ';' && ch != '\0') {
//		temp.clear();
//		if (isspace(ch))
//			ch = tree[++i];											//skip spaces
//		else if (ch == ',') {
//			ch = tree[++i];
//			flag_internal = false;
//		}
//		else if (ch == '(') {
//			levels->push(-1, current_ind);
//			levels->push_c(new Array<int>(0, 10));						//add level.
//			w_temp.push_c(new Array<double>(0, 10));
//			unlabeled.push(0);
//			unlabeled[current_ind]--;								//accumulate unlabeled point.
//
//
//			//levels_back_ptr = levels.back_c_ptr();
//
//			levels->push(current_ind);											//parent level.
//			levels->push(levels->get_c_ptr(current_ind)->get_size() - 1);		// parent node.
//
//			current_ind = levels->get_size() - 1;			//go to the new level in the back.
//			ch = tree[++i];
//			flag_internal = false;
//		}
//		else if (ch == ')') {
//			current_ind = levels->ele(current_ind, 0);
//			if (current_ind < 0)
//				break;
//			ch = tree[++i];
//
//			flag_internal = true;
//		}
//		else {
//			if (ch == ':') {
//				i++;
//				ch = read_sequence(tree, i, temp);
//				w_temp.push(atof(temp.c_str()), current_ind);
//				flag_internal = false;
//			}
//			else {
//				ch = read_sequence(tree, i, temp);
//				if (!flag_internal) {
//					if (flag_label == 0)			//labeled with integer
//						levels->push(ISROOTED ? atoi(temp.c_str()) : (atoi(temp.c_str()) - 1), current_ind);
//					else if (flag_label == 1)		//labeled with different normalization
//						levels->push(taxon_list.IndB2IndA[atoi(temp.c_str())], current_ind);
//					else if (flag_label == 2) {		// labeled with taxon
//						string taxon(temp.c_str());
//						auto it = taxon_list.Taxon2Ind.find(taxon);
//						if (it != taxon_list.Taxon2Ind.end())
//							levels->push(it->second, current_ind);
//						else {
//							std::cout << "Error: Attempt to find taxon: ``" << taxon << "'' that is not in the leaf set.\n";
//							throw(1);
//						}
//					}
//					else {
//						std::cout << "Error! Wrong flag for leaf's label.\n";
//						throw(1);
//					}
//
//					active_levels.push(current_ind);
//				}
//				flag_internal = false;
//			}
//
//		}
//	}
//
//	// Insert dummy leaf for rooted trees
//
//	if (ISROOTED) {
//		levels->push(0, 1);
//		w_temp.push(1.0, 1);
//		active_levels.push(1);
//	}
//
//	// Get edge size(bipartition size).
//
//	size = 0;
//
//	for (int i = 1; i < levels->get_size(); i++)
//		size += levels->ele(i).get_size() - 2;
//
//	//When first level has only 2 nodes, their parent node can be compressed.
//	//Therefore the first level becomes 2 points connected by 1 edge.
//	if (levels[1].get_size() == 4)
//		size--;
//
//	edges = new int* [2];
//	edges[0] = new int[size];	// node far-from-leaf
//	edges[1] = new int[size];	// node close-to-leaf
//
//
//
//	this->label_internal_node(active_levels, unlabeled, w_temp);
//
//	this->release_level();
//}

TreeArray::TreeArray(string &tree, TaxonList &taxon_list, Array<int> &active_levels,
					 Array<int> &unlabeled, int flag_label, bool ISROOTED, bool ISWEIGHTED)
{

	weights = nullptr;
	t2b = nullptr;
	edges = nullptr;

	issorted = false;
	int leaf_size = taxon_list.size;
	bool flag_internal = false;
	char ch;

	levels = new Array<Array<int>>(0, 2 * leaf_size - 3);

	Array<Array<double>> w_temp = Array<Array<double>>(0, 2 * leaf_size - 3);
	active_levels.clean();
	unlabeled.clean();

	int current_ind = -1;
	string temp;

	//darray<int>* levels_back_ptr;
	//darray<int>* levels_cur_ptr;

	//Build dummy level.

	levels->push(Array<int>(0, 10));
	unlabeled.push(0);
	//levels_back_ptr = levels.back_container();
	levels->back().push(current_ind); //parent level.
	levels->back().push(-1);		  //parent node.

	levels->back().push(-1); //dummy root.
	unlabeled.back()--;
	current_ind++;

	w_temp.push(Array<double>(0, 10)); //build the similar array for weights.
									   //Note that weights array do not have the first two columns for locating parent node.

	//Build first level.

	levels->push(Array<int>(0, 10));
	unlabeled.push(0);
	levels->back().push(current_ind);							  //parent level.
	levels->back().push(levels->ele(current_ind).get_size() - 1); //parent node.

	w_temp.push(Array<double>(0, 10));

	current_ind++;

	//Start scanning/

	size_t i = tree.find_first_of('(') + 1;
	ch = tree[i];

	while (ch != ';' && ch != '\0')
	{
		temp.clear();
		if (isspace(ch))
			ch = tree[++i]; //skip spaces
		else if (ch == ',')
		{
			ch = tree[++i];
			flag_internal = false;
		}
		else if (ch == '(')
		{
			levels->ele(current_ind).push(-1);
			levels->push(Array<int>(0, 10)); //add level.
			w_temp.push(Array<double>(0, 10));
			unlabeled.push(0);
			unlabeled[current_ind]--; //accumulate unlabeled point.

			//levels_back_ptr = levels.back_c_ptr();

			levels->back().push(current_ind);							  //parent level.
			levels->back().push(levels->ele(current_ind).get_size() - 1); // parent node.

			current_ind = levels->get_size() - 1; //go to the new level in the back.
			ch = tree[++i];
			flag_internal = false;
		}
		else if (ch == ')')
		{
			current_ind = levels->ele(current_ind)[0];
			if (current_ind < 0)
				break;
			ch = tree[++i];

			flag_internal = true;
		}
		else
		{
			if (ch == ':')
			{
				i++;
				ch = read_sequence(tree, i, temp);
				w_temp[current_ind].push(atof(temp.c_str()));
				flag_internal = false;
			}
			else
			{
				ch = read_sequence(tree, i, temp);
				if (!flag_internal)
				{
					if (flag_label == 0) //labeled with integer
						levels->ele(current_ind).push(ISROOTED ? atoi(temp.c_str()) : (atoi(temp.c_str()) - 1));
					else if (flag_label == 1) //labeled with different normalization
						levels->ele(current_ind).push(taxon_list.IndB2IndA[atoi(temp.c_str())]);
					else if (flag_label == 2)
					{ // labeled with taxon
						string taxon(temp.c_str());
						auto it = taxon_list.Taxon2Ind.find(taxon);
						if (it != taxon_list.Taxon2Ind.end())
							levels->ele(current_ind).push(it->second);
						else
						{
							std::cout << "Error: Attempt to find taxon: ``" << taxon << "'' that is not in the leaf set.\n";
							throw(1);
						}
					}
					else
					{
						std::cout << "Error! Wrong flag for leaf's label.\n";
						throw(1);
					}

					active_levels.push(current_ind);
				}
				flag_internal = false;
			}
		}
	}

	// Insert dummy leaf for rooted trees

	if (ISROOTED)
	{
		levels->ele(1).push(0);
		active_levels.push(1);
		if (ISWEIGHTED)
			w_temp[1].push(1.0);
	}

	// Get edge size(bipartition size).

	size = 0;

	for (int i = 1; i < levels->get_size(); i++)
		size += levels->ele(i).get_size() - 2;

	//When first level has only 2 nodes, their parent node can be compressed.
	//Therefore the first level becomes 2 points connected by 1 edge.
	if (levels->ele(1).get_size() == 4)
		size--;

	edges = new int *[2];
	edges[0] = new int[size]; // node far-from-leaf
	edges[1] = new int[size]; // node close-to-leaf

	if (ISWEIGHTED)
		this->label_internal_node(active_levels, unlabeled, w_temp);
	else
		this->label_internal_node(active_levels, unlabeled);

	this->release_level();
}

void TreeArray::label_internal_node(Array<int> &active_levels, Array<int> &unlabeled)
{
	int node_index = active_levels.get_size();
	int edge_index = 0;
	int cur_level, parent_level;
	// start labels of internal node from the number of leaves, N.
	// Note that leaves are labeled by 0 to N - 1.

	for (int i = 0; active_levels[i] != 0; i++)
	{
		cur_level = active_levels[i];
		parent_level = levels->ele(cur_level)[0];
		if (unlabeled[cur_level] == 0)
		{ // all nodes on this level is labeled
			int level_size = levels->ele(cur_level).get_size();
			for (int j = 2; j < level_size; j++)
			{ // add edge
				edges[0][edge_index] = node_index;
				edges[1][edge_index++] = levels->ele(cur_level)[j];
			}
			levels->ele(parent_level)[levels->ele(cur_level)[1]] = node_index++; // label the parent node.
			unlabeled[parent_level]++;
			unlabeled[cur_level]++;			  // mark current level done labeling.
			active_levels.push(parent_level); // make parent level active.
		}
	}

	if (levels[1].get_size() == 4) //	merge dummy node,
		edges[0][size - 1] = (edges[1][size - 1] == (*levels)[1][2] ? (*levels)[1][3] : (*levels)[1][2]);
}

void TreeArray::label_internal_node(Array<int> &active_levels, Array<int> &unlabeled, Array<Array<double>> &w_temp)
{
	int node_index = active_levels.get_size();
	assert(weights == nullptr);
	weights = new Array<double>(0, this->size);
	int edge_index = 0;
	int cur_level, parent_level;
	// start labels of internal node from the number of leaves, N.
	// Note that leaves are labeled by 0 to N - 1.

	for (int i = 0; active_levels[i] != 0; i++)
	{
		cur_level = active_levels[i];
		parent_level = levels->ele(cur_level)[0];
		if (unlabeled[cur_level] == 0)
		{ // all nodes on this level is labeled
			int level_size = levels->ele(cur_level).get_size();
			for (int j = 2; j < level_size; j++)
			{														// add edge and weight
				edges[0][edge_index] = node_index;					// far-from-leaf node (parent)
				edges[1][edge_index++] = levels->ele(cur_level)[j]; // close-to-leaf node (child)
				weights->push(w_temp[cur_level][j - 2]);			// add weights attached to the close-to-leaf node.
			}
			levels->ele(parent_level)[levels->ele(cur_level)[1]] = node_index++; // label the parent node.
			unlabeled[parent_level]++;
			unlabeled[cur_level]++;			  // mark current level done labeling.
			active_levels.push(parent_level); // mark the parent level active.
		}
	}

	(*weights).resize(size);

	if (levels[1].get_size() == 4) //	merge dummy node,
		edges[0][size - 1] = (edges[1][size - 1] == (*levels)[1][2] ? (*levels)[1][3] : (*levels)[1][2]);
}
