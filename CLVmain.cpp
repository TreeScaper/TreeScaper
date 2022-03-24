#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "array.hpp"
#include "zdtree.hpp"
#include "zdtreeobj.hpp"
#include "wstring.hpp"
#include "header_info.hpp"
#include "zdfileio.hpp"
#include "greedy_louvain.h"
#include "zdcommunity.hpp"
#include "NLDR.hpp"

using namespace std;
//using namespace FILEIO;

typedef unsigned char u_8;
typedef unsigned int u_32;
typedef unsigned long long u_64;

#define n_tree_key_option 14
#define n_reload_key_option 12
#define n_nldr_key_option 11
#define n_comm_key_option 18
#define n_adj_key_option 7
#define n_task_key_option 9

template <class Container_type>
TreeSetObjects<Container_type> *trees_driver(map<String, String> &paras, Container_type dummy);

template <class Container_type>
bool reload_driver(map<String, String> &paras, Container_type dummy);

Matrix<PRECISION>* nldr_driver(SpecMat::LowerTri<PRECISION>* Dis,map<String, String> &paras);

bool adj_driver(SpecMat::LowerTri<PRECISION>* Adj, map<String, String> &paras);

bool comm_driver(SpecMat::LowerTri<PRECISION>* Adj, map<String, String> &paras);

template <class Container_type>
bool task_driver(map<String, String> &paras, Container_type dummy);

int main(int argc, char *argv[])
{

	if (argc > 1 && (String)argv[1] == (String) "-trees")
	{
		String default_paras[n_tree_key_option] = {"", "postfix", "taxon", "64", "0", "1", "Cov", "1", "0.01", "none", "none", "0", "none", "0"};
		String options[n_tree_key_option] = {"-f", "-post", "-tm", "-bit", "-r", "-w", "-output", "-bipart-hf", "-bipart-lf", "-sb", "-st", "-saveobj", "-trees-key", "-lessprint"};

		map<String, String> paras = read_paras(argc, argv, n_tree_key_option, default_paras, options);

		if (paras["-trees-key"] != (String) "none")
		{
			read_paras_from_csv(paras["-trees-key"], paras, true);
		}

		if (paras["-bit"] == (String) "8")
			trees_driver(paras, (u_8)1);
		else if (paras["-bit"] == (String) "32")
			trees_driver(paras, (u_32)1);
		else if (paras["-bit"] == (String) "64")
			trees_driver(paras, (u_64)1);
		else
			std::cout << "Bitstring container size not supported. Please choose value from {8, 32, 64} for key `-bit'.\n";
	}
	else if (argc > 1 && (String)argv[1] == (String) "-reload")
	{
		String default_paras[n_reload_key_option] = {"", "obj", "postfix", "64", "0", "1", "COV", "1", "0.01", "none", "none", "none"};
		String options[n_reload_key_option] = {"-f", "-ft", "-post", "-bit", "-r", "-w", "-output", "-hf", "-lf", "-sb", "-st", "-reload-key"};

		map<String, String> paras = read_paras(argc, argv, n_reload_key_option, default_paras, options);

		if (paras["-key"] != (String) "none")
			read_paras_from_csv(paras["-reload-key"], paras, true);

		if (paras["-bit"] == (String) "8")
			reload_driver(paras, (u_8)1);
		else if (paras["-bit"] == (String) "32")
			reload_driver(paras, (u_32)1);
		else if (paras["-bit"] == (String) "64")
			reload_driver(paras, (u_64)1);
		else
			std::cout << "Bitstring container size not supported. Please choose value from {8, 32, 64} for key `-bit'.\n";
	}
	else if (argc > 1 && (String)argv[1] == (String) "-adj")
	{
		String default_paras[n_adj_key_option] = { "", "-1", "DIS", "EXP", "time", "none", "0"};
		String options[n_adj_key_option] = { "-f", "-adj-size", "-adj-type", "-am", "-post", "-adj-key", "-lessprint"};
		
		map<String, String> paras = read_paras(argc, argv, n_adj_key_option, default_paras, options);

		if (paras["-key"] != (String) "none")
		{
			std::cout << "Reading parameters from " << paras["-adj-key"] << endl;
			read_paras_from_csv(paras["-adj-key"], paras, true);
		}
		
		adj_driver(nullptr, paras);
	}
	else if (argc > 1 && (String)argv[1] == (String) "-nldr")
	{
		String default_paras[n_nldr_key_option] = {"", "COR", "postfix", "-1", "2", "CCA", "STOCHASTIC", "RAND", "1", "nldr_parameters.csv", "0"};
		String options[n_nldr_key_option] = {"-f", "-nldr-type", "-post", "-nldr-size", "-nldr-dim", "-nldr-cost", "-nldr-algo", "-nldr-init", "-nldr-seed", "-nldr-key", "-lessprint"};

		map<String, String> paras = read_paras(argc, argv, n_nldr_key_option, default_paras, options);

		if (paras["-key"] != (String) "none")
		{
			std::cout << "Reading parameters from " << paras["-nldr-key"] << endl;
			read_paras_from_csv(paras["-nldr-key"], paras, true);
		}

		nldr_driver(nullptr, paras);
	}
	else if (argc > 1 && (String)argv[1] == (String) "-comm")
	{
		String default_paras[n_comm_key_option] = { "", "BIPART", "-1", "CNM", "1", "0", "1",
			"0", "1", "0", "1", "0", "1", "0", "AUTO", "comm_parameters.csv", "postfix", "0"};
		String options[n_comm_key_option] = { "-f", "-comm-type", "-comm-size", "-cm", "-lp", "-lps", "-lpe",
			"-lpiv", "-ln", "-lns", "-lne", "-lniv", "-comm-hf", "-comm-lf", "-lm", "-comm-key", "-post", "-lessprint"};

		map<String, String> paras = read_paras(argc, argv, n_comm_key_option, default_paras, options);

		if (paras["-comm-key"] != (String) "none")
		{
			std::cout << "Reading parameters from " << paras["-comm-key"] << endl;
			read_paras_from_csv(paras["-comm-key"], paras, true);
		}


		comm_driver(nullptr, paras);
	}
	else if (argc > 1 && (String)argv[1] == (String) "-task")
	{
		String default_paras[n_task_key_option] = {"", "postfix", "none", "none", "none", "none", "none", "0", "64"};
		String options[n_task_key_option] = {"-f", "-post", "-key", "-trees", "-adj", "-comm", "-nldr", "-lessprint", "-bit"};

		map<String, String> paras = read_paras(argc, argv, n_task_key_option, default_paras, options);


		if (paras["-bit"] == (String) "8")
			task_driver(paras, (u_8)1);
		else if (paras["-bit"] == (String) "32")
			task_driver(paras, (u_32)1);
		else if (paras["-bit"] == (String) "64")
			task_driver(paras, (u_64)1);
		else
			std::cout << "Bitstring container size not supported. Please choose value from {8, 32, 64} for key `-bit'.\n";
	}
	return 0;
}

template <class Container_type>
TreeSetObjects<Container_type> *trees_driver(map<String, String> &paras, Container_type dummy)
{
	std::clock_t start, end;
	int flag_label;
	bool rooted = (paras["-r"] == (String) "1");
	bool weighted = (paras["-w"] == (String) "1");
	bool less_print = (paras["-lessprint"] == (String) "1");

	std::cout << (rooted ? "rooted" : "unrooted") << "\t" << (weighted ? "weighted" : "unweighted") << '\n';

	if (paras["-tm"] == (String) "TAXON")
		flag_label = 2;
	else
		flag_label = 0;
	// flag_label:	0:	labelled with integer
	//				1:	labelled with integer from different normalization
	//				2:	labelled with taxon name
	int n_tree = 100;
	std::string fname((char *)paras["-f"]);

	std::cout << "Reading taxa...\n";
	start = std::clock();
	auto Taxa_ptr = new TaxonList();
	size_t tree_pos = Taxa_ptr->ReadTaxa(fname);

	std::cout << "Reading trees...\n";
	Taxa_ptr->set_bitstr_size((Container_type)1);
	// Bipartition<Container_type> Bipart(&Taxa, n_tree);
	auto Bipart_ptr = new Bipartition<Container_type>(Taxa_ptr, n_tree);
	std::cout << "Bipartition initialized...\n";
	TreeSet<Container_type> trees(Bipart_ptr, Taxa_ptr);
	trees.ReadTree(fname, tree_pos, flag_label, rooted, weighted);
	end = std::clock();
	cout << "Read tree time(s):\t\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

	auto treesobj_ptr = new TreeSetObjects<Container_type>(&(trees), rooted, weighted);

	std::cout << "Computing bipartition...\n";
	start = std::clock();
	treesobj_ptr->Compute_Bipart();
	end = std::clock();
	cout << "Compute Bipartition time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

	std::cout << "Computing Bipart2Tree sparse matrix...\n";
	start = std::clock();
	treesobj_ptr->Compute_Bipart_Matrix();
	end = std::clock();
	cout << "Compute Bipart2Tree sparse matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

	ofstream fout;
	// Print Bipartition List
	String outname_Bipartcnt("Bipartition_Count");
	outname_Bipartcnt.make_stdname(paras);
	Header_info Header_Bipartcnt;
	Header_Bipartcnt.insert("created", time_stamp());
	Header_Bipartcnt.insert("output type", "Bipartition count");
	Header_Bipartcnt.insert("size", to_string(treesobj_ptr->get_unique_bipart_num()));
	Header_Bipartcnt.insert("format", "bipartition id, bitstring, appear times");
	Header_Bipartcnt.insert("source", paras["-f"]);

	fout.open((char *)outname_Bipartcnt);
	fout << Header_Bipartcnt;
	treesobj_ptr->print_Bipart_List(fout);
	fout.close();
	cout << "Sucessfully printed bipartition list info in Bipartition_count_" << paras["-post"] << ".out file.\n\n";

	// Print Bipartition2Tree Sparse Matrix
	String outname_Bipart("Bipartition");
	outname_Bipart.make_stdname(paras);
	Header_info Header_Bipart;
	Header_Bipart.insert("created", time_stamp());
	Header_Bipart.insert("output type", "Bipartition-to-Tree sparse matrix");
	string bxt_temp = to_string(treesobj_ptr->get_unique_bipart_num()) + (string) " x " + to_string(treesobj_ptr->get_tree_set_size());
	Header_Bipart.insert("size", bxt_temp);
	Header_Bipart.insert("format", "bipartition id, tree id, values");
	Header_Bipart.insert("source", paras["-f"]);

	fout.open((char *)outname_Bipart);
	fout << Header_Bipart;
	treesobj_ptr->print_Bipart2Tree_Matrix(fout, RCVLIST);
	fout.close();
	cout << "Sucessfully printed bipartition-to-tree sparse matrix in Bipartition_" << (char *)paras["-post"] << ".out file.\n\n";

	if (paras["-output"] == (String) "COV")
	{
		Array<size_t> subset, id_mapping;
		double hf = atof(paras["-bipart-hf"]), lf = atof(paras["-bipart-lf"]);
		if (paras["-sb"] != (String) "none")
		{
			if (!read_subset(paras["-sb"], subset))
			{
				std::cout << "Fail to open file for subset of bipartitions. Using default frequency setting (lf, hf) = \t(" << lf << ", " << hf << ").\n";
				treesobj_ptr->bipart_frequency_check(lf, hf, subset);
			}
		}
		else
			treesobj_ptr->bipart_frequency_check(lf, hf, subset);

		std::cout << "Computing covariance matrix...\n";
		start = std::clock();
		treesobj_ptr->Compute_Covariance_Matrix(subset, id_mapping);
		end = std::clock();
		cout << "Compute covariance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

		// Print Covariance Matrix
		String outname_ID("Bipart_ID_for_Cov");
		outname_ID.make_stdname(paras);
		Header_info Header_ID;
		Header_ID.insert("created", time_stamp());
		Header_ID.insert("output type", "Bipartition id of the covariance matrix");
		Header_ID.insert("size", to_string(subset.get_size()));
		Header_ID.insert("source", outname_Bipartcnt);

		fout.open((char *)outname_ID);
		fout << Header_ID;
		for (int i = 0; i < subset.get_size(); i++)
			fout << subset[i] << '\n';
		fout.close();
		cout << "Sucessfully printed bipartition id of sub-covariance matrix in Bipart_ID_for_Cov_" << (char *)paras["-post"] << ".out file.\n\n";

		// Print Covariance Matrix
		if (!less_print)
		{
			String outname_Cova("Covariance");
			outname_Cova.make_stdname(paras);
			Header_info Header_Cova;
			Header_Cova.insert("created", time_stamp());
			Header_Cova.insert("output type", "Covariance matrix");
			string cov_size_temp = to_string(subset.get_size()) + (string) " x " + to_string(subset.get_size());
			Header_Cova.insert("size", cov_size_temp);
			Header_Cova.insert("bipartition id", outname_ID);
			Header_Cova.insert("source", paras["-f"]);

			fout.open((char *)outname_Cova);
			fout << Header_Cova;
			treesobj_ptr->print_Covariance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed covariance matrix in Covariance_" << (char *)paras["-post"] << ".out file.\n\n";
		}
	}
	else if (paras["-output"] == (String) "DIS")
	{
		Array<size_t> subset, id_mapping;
		bool flag_subtree = false;
		if (paras["-st"] != (String) "none")
		{
			if (!read_subset(paras["-st"], subset))
			{
				std::cout << "Fail to open file for subset of tree. Compute the full tree set.\n";
				flag_subtree = false;
			}
			else
				flag_subtree = true;
		}

		String outname_ID("Tree_ID_for_Dist");
		outname_ID.make_stdname(paras);

		if (flag_subtree)
		{
			Header_info Header_ID;
			Header_ID.insert("created", time_stamp());
			Header_ID.insert("output type", "Tree id of the distance matrix");
			Header_ID.insert("size", to_string(subset.get_size()));
			Header_ID.insert("source", paras["-f"]);
			fout.open((char *)outname_ID);
			fout << Header_ID;
			for (int i = 0; i < subset.get_size(); i++)
				fout << subset[i] << '\n';
			fout.close();
			cout << "Sucessfully printed tree id of sub-distance matrix in Tree_ID_for_Dis" << (char *)paras["-post"] << ".out file.\n\n";
		}

		std::cout << "Computing RF-distance matrix...\n";
		start = std::clock();
		if (flag_subtree)
			treesobj_ptr->Compute_RF_Distance_Matrix(subset, id_mapping);
		else
			treesobj_ptr->Compute_RF_Distance_Matrix();
		end = std::clock();
		cout << "Compute RF-distance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

		// Print Distance Matrix
		String outname_Dist("Distance");
		outname_Dist.make_stdname(paras);

		Header_info Header_Dist;
		Header_Dist.insert("created", time_stamp());
		Header_Dist.insert("output type", "Distance matrix");
		string dist_size_temp;
		if (flag_subtree)
			dist_size_temp = to_string(subset.get_size()) + (string) " x " + to_string(subset.get_size());
		else
			dist_size_temp = to_string(treesobj_ptr->get_tree_set_size()) + (string) " x " + to_string(treesobj_ptr->get_tree_set_size());
		Header_Dist.insert("size", dist_size_temp);
		if (flag_subtree)
			Header_Dist.insert("Tree id", outname_ID);
		Header_Dist.insert("source", paras["-f"]);

		fout.open((char *)outname_Dist);
		fout << Header_Dist;
		treesobj_ptr->print_Distance_Matrix(fout);
		fout.close();
		cout << "Sucessfully printed distance matrix in Distance_" << (char *)paras["-post"] << ".out file.\n\n";
	}
	else if (paras["-output"] == (String) "BOTH")
	{

	}
	treesobj_ptr->print_summary();
	if (paras["-saveobj"] == (String) "1")
	{
		String outname_Obj("Object_bin");
		outname_Obj.make_stdname(paras);
		treesobj_ptr->save_treeobj(outname_Obj, true, false, true);
	}
	treesobj_ptr->pop_trees();
	trees.release_tree();

	return treesobj_ptr;
}

template <class Container_type>
bool reload_driver(map<String, String> &paras, Container_type dummy)
{
	std::clock_t start, end;
	bool rooted = (paras["-r"] == (String) "1");
	bool weighted = (paras["-w"] == (String) "1");
	bool filetype = (paras["-ft"] == (String) "obj");
	bool tree_s = (paras["-st"] != (String) "none");
	bool bipart_s = (paras["-sb"] != (String) "none");

	std::cout << "reloading " << (rooted ? "rooted" : "unrooted") << "\t" << (weighted ? "weighted" : "unweighted") << (filetype ? "Tree objects\n" : "Treeset\n");

	if (filetype)
	{
		TreeSetObjects<Container_type> treesobj;
		treesobj.read_treeobj(paras["-f"]);

		Array<size_t> subset, id_mapping;

		if (paras["-output"] == (String) "COV")
		{
			double hf = atof(paras["-hf"]), lf = atof(paras["-lf"]);
			if (bipart_s)
			{
				if (!read_subset(paras["-sb"], subset))
				{
					std::cout << "Fail to open file for subset of bipartitions. Using default frequency setting (lf, hf) = \t(" << lf << ", " << hf << ").\n";
					treesobj.bipart_frequency_check(lf, hf, subset);
				}
			}
			else
			{
				std::cout << "Using default frequency setting (lf, hf) = \t(" << lf << ", " << hf << ") for computing covariance matrix.\n";
				treesobj.bipart_frequency_check(lf, hf, subset);
			}

			std::cout << "Computing covariance matrix...\n";
			start = std::clock();
			treesobj.Compute_Covariance_Matrix(subset, id_mapping);
			end = std::clock();
			cout << "Compute covariance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_ID("Bipart_ID_for_Cov");
			outname_ID.make_stdname(paras);
			Header_info Header_ID;
			Header_ID.insert("created", time_stamp());
			Header_ID.insert("output type", "Bipartition id of the covariance matrix");
			Header_ID.insert("size", to_string(subset.get_size()));

			std::ofstream fout;
			fout.open((char *)outname_ID);
			fout << Header_ID;
			for (int i = 0; i < subset.get_size(); i++)
				fout << subset[i] << '\n';
			fout.close();
			cout << "Sucessfully printed bipartition id of sub-covariance matrix in Bipart_ID_for_Cov" << (char *)paras["-post"] << ".out file.\n\n";

			// Print Covariance Matrix
			String outname_Cova("Covariance");
			outname_Cova.make_stdname(paras);
			Header_info Header_Cova;
			Header_Cova.insert("created", time_stamp());
			Header_Cova.insert("output type", "Covariance matrix");
			string cov_size_temp = to_string(subset.get_size()) + (string) " x " + to_string(subset.get_size());
			Header_Cova.insert("size", cov_size_temp);
			Header_Cova.insert("bipartition id", outname_ID);
			Header_Cova.insert("source", paras["-f"]);

			fout.open((char *)outname_Cova);
			fout << Header_Cova;
			treesobj.print_Covariance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed covariance matrix in Covariance_" << (char *)paras["-post"] << ".out file.\n\n";
		}
		else if (paras["-output"] == (String) "DIS")
		{
			if (tree_s)
			{
				if (!read_subset(paras["-st"], subset))
				{
					std::cout << "Fail to open file for subset of tree. Compute the full tree set.\n";
					tree_s = false;
				}
				else
					tree_s = false;
			}

			std::cout << "Computing RF-distance matrix...\n";
			start = std::clock();
			std::cout << subset.get_size() << '\n';

			if (subset.get_size() != 0)
				tree_s = true;
			if (tree_s)
				treesobj.Compute_RF_Distance_Matrix(subset, id_mapping);
			else
				treesobj.Compute_RF_Distance_Matrix();

			// treesobj.Compute_RF_Distance_Matrix(subset, id_mapping);
			end = std::clock();
			cout << "Compute RF-distance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			String outname_Dist("Distance");
			outname_Dist.make_stdname(paras);

			Header_info Header_Dist;
			Header_Dist.insert("created", time_stamp());
			Header_Dist.insert("output type", "Distance matrix");
			string dist_size_temp;
			if (tree_s)
				dist_size_temp = to_string(subset.get_size()) + (string) " x " + to_string(subset.get_size());
			else
				dist_size_temp = to_string(treesobj.get_tree_set_size()) + (string) " x " + to_string(treesobj.get_tree_set_size());
			Header_Dist.insert("size", dist_size_temp);
			if (tree_s)
				Header_Dist.insert("Tree id", paras["-st"]);
			Header_Dist.insert("source", paras["-f"]);

			std::ofstream fout;
			fout.open((char *)outname_Dist);
			fout << Header_Dist;
			treesobj.print_Distance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed distance matrix in Distance_" << (char *)paras["-post"] << ".out file.\n\n";
		}

		treesobj.print_summary();
	}

	return true;
}

Matrix<PRECISION>* nldr_driver(SpecMat::LowerTri<PRECISION>* Dis, map<String, String> &paras)
{
	Header_info info;
	int size = 0, file_size = 0;
	int dim = 0;
	int seed = atoi((char *)paras["-nldr-seed"]);
    std::default_random_engine eng{seed};
	bool less_print = (paras["-lessprint"] == (String) "1");

	if (Dis == nullptr) // Adjacency should be read from file.
	{
		std::fstream fin;
		fin.open((char *)paras["-f"]);
		fin >> info;
		int pos_header = Header::end_header(fin);
		fin.seekg(pos_header, ios::beg);

		int size = atoi(paras["-nldr-size"]);
		int file_size = info.count("size") ? atoi(info["size"]) : get_file_lines(fin);
		
		if (size < 0)
			size = file_size;
		else if (size != file_size)
		{
			std::cout << "Error: Matrix size indicated in parameters, " << file_size << " , is different than the dimension found from the input file, " << size << ".";
			return;
		}
		paras["-nldr-size"] = to_string(size);

		Dis = new SpecMat::LowerTri<PRECISION>(size, fin);

		if (info.count("distance_type") && !paras.count("-dm")) {
			std::string dm = (char*) info["distance_type"];
			if (dm.find("Unweighted Robinson-Foulds") != std::string::npos || dm.find("URF") != std::string::npos)
				paras["-dm"] = (String) "URF";
			else if (dm.find("Robinson-Foulds") != std::string::npos || dm.find("RF") != std::string::npos)
				paras["-dm"] = (String) "RF";
			else if (dm.find("Matching") != std::string::npos || dm.find("Mat") != std::string::npos)
				paras["-dm"] = (String) "Mat";
			else if (dm.find("SPR") != std::string::npos)
				paras["-dm"] = (String) "SPR";
			else
				paras["-dm"] = (String) "";
		}
	}
	else
	{
		size = Dis->dimension();
	}

	dim = atoi((char *)paras["-nldr-dim"]);

	auto CX = new Matrix<PRECISION>(size, dim);
	if (paras["-nldr-init"] == (String) "RAND")
	{
		PRECISION r_start, r_stop;
		r_start = atof(paras["RANDOM_START"]);
		r_stop = atof(paras["RANDOM_STOP"]);

    	uniform_real_distribution<PRECISION> rand(r_start, r_stop);

    	for (auto i = 0; i < size; i++)
        	for (auto k = 0; k < dim; k++)
            	(*CX)(i, k) = rand(eng);
	}
	else
	{
		std::cout << "Initialization not supported yet.\n";
		throw(1);
	}

    SpecMat::LowerTri<PRECISION> DX(size, true);
	NLDR::Cor2Dis(*CX, DX);
	DX.form_row_ptr();

	
	uniform_int_distribution<> rand_uniform(0, size - 1);



	unsigned MaxIter = 0;
    unsigned MaxTime = -1;
    PRECISION ErrTol = 0;
    PRECISION DifTol = -1.0;

	PRECISION lambda0 = 12000;
    PRECISION lambdan = 1;

	PRECISION *lambda = nullptr;
	PRECISION *alpha = nullptr;
	unsigned *s_ind = nullptr;

	PRECISION metropolis_T = 0;


	if (paras["-nldr-cost"] == (String) "KRUSKAL1" && paras["-nldr-algo"] == (String) "STOCHASTIC")
	{
		MaxIter = atoi((char *) paras["KRU_STO_MAXITER"]);
		MaxTime = atoi((char *) paras["KRU_STO_MAXTIME"]);
		ErrTol = atof((char *)paras["KRU_STO_ERRTOL"]);

		s_ind = new unsigned[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
	        s_ind[i] = rand_uniform(eng);

    	alpha = new PRECISION[MaxIter];
		PRECISION alpha0 = atof(paras["KRU_STO_ALPHA0"]), alphan = atof(paras["KRU_STO_ALPHAn"]);;
		for (auto i = 0; i < MaxIter; i++)
        	alpha[i] = alpha0 * pow(alphan, (double)i / (MaxIter - 1));
		

		NLDR::STOCHASTIC(Dis, size, dim, CX, &DX, KRUSKAL1, MaxIter, MaxTime, ErrTol, DifTol, alpha, nullptr, s_ind);
		delete[] lambda;
		delete[] s_ind;
	}
	else if(paras["-nldr-cost"] == (String) "KRUSKAL1" && paras["-nldr-algo"] == (String) "MAJORIZATION")
	{
		MaxIter = atoi((char *) paras["KRU_MAJ_MAXITER"]);
		MaxTime = atoi((char *) paras["KRU_MAJ_MAXTIME"]);
		ErrTol = atof((char *)paras["KRU_MAJ_ERRTOL"]);

		NLDR::MAJORIZATION(Dis, size, dim, CX, &DX, KRUSKAL1, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if(paras["-nldr-cost"] == (String) "KRUSKAL1" && paras["-nldr-algo"] == (String) "METROPOLIS")
	{
		MaxIter = atoi((char *)paras["KRU_MET_MAXITER"]);
		MaxTime = atoi((char *)paras["KRU_MET_MAXTIME"]);
		ErrTol = atof((char *)paras["KRU_MET_ERRTOL"]);

		metropolis_T = atof((char *)paras["KRU_MET_T"]);

		NLDR::METROPOLIS(Dis, size, dim, CX, &DX, KRUSKAL1, MaxIter, MaxTime, ErrTol, DifTol, nullptr, &eng, metropolis_T);
	}
	else if(paras["-nldr-cost"] == (String) "KRUSKAL1" && paras["-nldr-algo"] == (String) "GAUSS_SEIDEL")
	{
		MaxIter = atoi((char *)paras["KRU_GAU_MAXITER"]);
		MaxTime = atoi((char *)paras["KRU_GAU_MAXTIME"]);
		ErrTol = atof((char *)paras["KRU_GAU_ERRTOL"]);

		NLDR::GAUSS_SEIDEL(Dis, size, dim, CX, &DX, KRUSKAL1, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if (paras["-nldr-cost"] == (String) "NORMALIZED" && paras["-nldr-algo"] == (String) "STOCHASTIC")
	{
		MaxIter = atoi((char *) paras["NOR_STO_MAXITER"]);
		MaxTime = atoi((char *) paras["NOR_STO_MAXTIME"]);
		ErrTol = atof((char *)paras["NOR_STO_ERRTOL"]);

		s_ind = new unsigned int[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
	        s_ind[i] = rand_uniform(eng);

		alpha = new PRECISION[MaxIter];
		PRECISION alpha0 = atof(paras["NOR_STO_ALPHA0"]), alphan = atof(paras["NOR_STO_ALPHAn"]);;
		for (auto i = 0; i < MaxIter; i++)
        	alpha[i] = alpha0 * pow(alphan, (double)i / (MaxIter - 1));

		NLDR::STOCHASTIC(Dis, size, dim, CX, &DX, NORMALIZED, MaxIter, MaxTime, ErrTol, DifTol, alpha, nullptr, s_ind);
		delete[] lambda;
		delete[] s_ind;
	}
	else if(paras["-nldr-cost"] == (String) "NORMALIZED" && paras["-nldr-algo"] == (String) "MAJORIZATION")
	{
		MaxIter = atoi((char *) paras["NOR_MAJ_MAXITER"]);
		MaxTime = atoi((char *) paras["NOR_MAJ_MAXTIME"]);
		ErrTol = atof((char *)paras["NOR_MAJ_ERRTOL"]);

		NLDR::MAJORIZATION(Dis, size, dim, CX, &DX, NORMALIZED, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if(paras["-nldr-cost"] == (String) "NORMALIZED" && paras["-nldr-algo"] == (String) "METROPOLIS")
	{
		MaxIter = atoi((char *)paras["NOR_MET_MAXITER"]);
		MaxTime = atoi((char *)paras["NOR_MET_MAXTIME"]);
		ErrTol = atof((char *)paras["NOR_MET_ERRTOL"]);

		metropolis_T = atof((char *)paras["NOR_MET_T"]);

		NLDR::METROPOLIS(Dis, size, dim, CX, &DX, NORMALIZED, MaxIter, MaxTime, ErrTol, DifTol, nullptr, &eng, metropolis_T);
	}
	else if(paras["-nldr-cost"] == (String) "NORMALIZED" && paras["-nldr-algo"] == (String) "GAUSS_SEIDEL")
	{
		MaxIter = atoi((char *)paras["NOR_GAU_MAXITER"]);
		MaxTime = atoi((char *)paras["NOR_GAU_MAXTIME"]);
		ErrTol = atof((char *)paras["NOR_GAU_ERRTOL"]);

		NLDR::GAUSS_SEIDEL(Dis, size, dim, CX, &DX, NORMALIZED, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if (paras["-nldr-cost"] == (String) "SAMMON" && paras["-nldr-algo"] == (String) "STOCHASTIC")
	{
		MaxIter = atoi((char *) paras["NLM_STO_MAXITER"]);
		MaxTime = atoi((char *) paras["NLM_STO_MAXTIME"]);
		ErrTol = atof((char *)paras["NLM_STO_ERRTOL"]);

		s_ind = new unsigned int[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
	        s_ind[i] = rand_uniform(eng);

		alpha = new PRECISION[MaxIter];
		PRECISION alpha0 = atof(paras["NLM_STO_ALPHA0"]), alphan = atof(paras["NLM_STO_ALPHAn"]);;
		for (auto i = 0; i < MaxIter; i++)
        	alpha[i] = alpha0 * pow(alphan, (double)i / (MaxIter - 1));

		NLDR::STOCHASTIC(Dis, size, dim, CX, &DX, SAMMON, MaxIter, MaxTime, ErrTol, DifTol, alpha, nullptr, s_ind);
		delete[] lambda;
		delete[] s_ind;
	}
	else if(paras["-nldr-cost"] == (String) "SAMMON" && paras["-nldr-algo"] == (String) "MAJORIZATION")
	{
		MaxIter = atoi((char *) paras["NLM_MAJ_MAXITER"]);
		MaxTime = atoi((char *) paras["NLM_MAJ_MAXTIME"]);
		ErrTol = atof((char *)paras["NLM_MAJ_ERRTOL"]);

		NLDR::MAJORIZATION(Dis, size, dim, CX, &DX, SAMMON, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if(paras["-nldr-cost"] == (String) "SAMMON" && paras["-nldr-algo"] == (String) "METROPOLIS")
	{
		MaxIter = atoi((char *)paras["NLM_MET_MAXITER"]);
		MaxTime = atoi((char *)paras["NLM_MET_MAXTIME"]);
		ErrTol = atof((char *)paras["NLM_MET_ERRTOL"]);

		metropolis_T = atof((char *)paras["NLM_MET_T"]);

		NLDR::METROPOLIS(Dis, size, dim, CX, &DX, SAMMON, MaxIter, MaxTime, ErrTol, DifTol, nullptr, &eng, metropolis_T);
	}
	else if(paras["-nldr-cost"] == (String) "SAMMON" && paras["-nldr-algo"] == (String) "GAUSS_SEIDEL")
	{
		MaxIter = atoi((char *)paras["NLM_GAU_MAXITER"]);
		MaxTime = atoi((char *)paras["NLM_GAU_MAXTIME"]);
		ErrTol = atof((char *)paras["NLM_GAU_ERRTOL"]);

		NLDR::GAUSS_SEIDEL(Dis, size, dim, CX, &DX, SAMMON, MaxIter, MaxTime, ErrTol, DifTol, nullptr);
	}
	else if (paras["-nldr-cost"] == (String) "CCA" && paras["-nldr-algo"] == (String) "STOCHASTIC")
	{
		MaxIter = atoi((char *) paras["CCA_STO_MAXITER"]);
		MaxTime = atoi((char *) paras["CCA_STO_MAXTIME"]);
		ErrTol = atof((char *)paras["CCA_STO_ERRTOL"]);

		lambda0 = atof((char *) paras["CCA_STO_ALPHA0"]);
		lambdan = atof((char *) paras["CCA_STO_ALPHAn"]);

		lambda = new PRECISION[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
        	lambda[i] = lambda0 * pow((lambdan / lambda0), (double)i / (MaxIter - 1));

		s_ind = new unsigned int[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
	        s_ind[i] = rand_uniform(eng);

		alpha = new PRECISION[MaxIter];
		PRECISION alpha0 = atof(paras["CCA_STO_ALPHA0"]), alphan = atof(paras["CCA_STO_ALPHAn"]);;
		for (auto i = 0; i < MaxIter; i++)
        	alpha[i] = alpha0 * pow(alphan, (double)i / (MaxIter - 1));

		NLDR::STOCHASTIC(Dis, size, dim, CX, &DX, CCA, MaxIter, MaxTime, ErrTol, DifTol, alpha, lambda, s_ind);
		delete[] lambda;
		delete[] s_ind;
	}
	else if(paras["-nldr-cost"] == (String) "CCA" && paras["-nldr-algo"] == (String) "MAJORIZATION")
	{
		MaxIter = atoi((char *) paras["CCA_MAJ_MAXITER"]);
		MaxTime = atoi((char *) paras["CCA_MAJ_MAXTIME"]);
		ErrTol = atof((char *)paras["CCA_MAJ_ERRTOL"]);

		lambda0 = atof((char *) paras["CCA_MAJ_ALPHA0"]);
		lambdan = atof((char *) paras["CCA_MAJ_ALPHAn"]);
		lambda = new PRECISION[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
        	lambda[i] = lambda0 * pow((lambdan / lambda0), (double)i / (MaxIter - 1));

		NLDR::MAJORIZATION(Dis, size, dim, CX, &DX, CCA, MaxIter, MaxTime, ErrTol, DifTol, lambda);
		delete[] lambda;
	}
	else if(paras["-nldr-cost"] == (String) "CCA" && paras["-nldr-algo"] == (String) "METROPOLIS")
	{
		MaxIter = atoi((char *)paras["CCA_MET_MAXITER"]);
		MaxTime = atoi((char *)paras["CCA_MET_MAXTIME"]);
		ErrTol = atof((char *)paras["CCA_MET_ERRTOL"]);

		metropolis_T = atof((char *)paras["CCA_MET_T"]);

		lambda0 = atof((char *) paras["CCA_MET_ALPHA0"]);
		lambdan = atof((char *) paras["CCA_MET_ALPHAn"]);
		lambda = new PRECISION[MaxIter];
		for (auto i = 0; i < MaxIter; i++)
        	lambda[i] = lambda0 * pow((lambdan / lambda0), (double)i / (MaxIter - 1));

		NLDR::METROPOLIS(Dis, size, dim, CX, &DX, CCA, MaxIter, MaxTime, ErrTol, DifTol, lambda, &eng, metropolis_T);
		delete[] lambda;
	}
	else if(paras["-nldr-cost"] == (String) "CCA" && paras["-nldr-algo"] == (String) "GAUSS_SEIDEL")
	{
		MaxIter = atoi((char *)paras["CCA_GAU_MAXITER"]);
		MaxTime = atoi((char *)paras["CCA_GAU_MAXTIME"]);
		ErrTol = atof((char *)paras["CCA_GAU_ERRTOL"]);

		lambda0 = atof((char *) paras["CCA_GAU_ALPHA0"]);
		lambdan = atof((char *) paras["CCA_GAU_ALPHAn"]);
		lambda = new PRECISION[MaxIter];
		auto sMaxIter = MaxIter * size;
		for (auto i = 0; i < sMaxIter; i++)
        	lambda[i] = lambda0 * pow((lambdan / lambda0), (double)i / (sMaxIter - 1));

		NLDR::GAUSS_SEIDEL(Dis, size, dim, CX, &DX, CCA, MaxIter, MaxTime, ErrTol, DifTol, lambda);
		delete[] lambda;
	}

	(*Dis) = DX;

	info.insert("created", time_stamp());
	info.insert("output_type", "Coordinates matrix");
	info.insert("size", to_string(size));
	info.insert("dimension", to_string(dim));
	info.insert("source", paras["-f"]);

	String outname_cor("Coordinate");
	outname_cor.make_stdname(paras);
	std::ofstream file_Cor;
	file_Cor.open((char*)outname_cor);
	file_Cor << info;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j <= dim; j++)
			file_Cor << (*CX)(i, j) << "\t";
		file_Cor << std::endl;
	}

	if (!less_print)
	{
		info.insert("created", time_stamp());
		info.insert("output_type", "Distance matrix on Euclidean space");
		info.insert("size", to_string(size));
		info.insert("dimension", to_string(dim));
		info.insert("source", paras["-f"]);
		String outname_nldr("Distance_NLDR");
		outname_nldr.make_stdname(paras);
		std::ofstream file_NLDR;
		file_NLDR.open((char*)outname_nldr);
		file_NLDR << info;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j <= size; j++)
				file_NLDR << DX(i, j) << "\t";
			file_NLDR << std::endl;
		}
	}

	if (lambda != nullptr)
		delete[] lambda;
	if (alpha != nullptr)
		delete[] alpha;
	if (s_ind != nullptr)
		delete[] s_ind;
	
	return CX;
}

bool adj_driver(SpecMat::LowerTri<PRECISION>* Adj, map<String, String> &paras)
{
	Header_info info;
	int size = atoi(paras["-adj-size"]);
	bool less_print = (paras["-lessprint"] == (String) "1");

	if (Adj == nullptr) // Adjacency should be read from file.
	{
		std::fstream fin;
		fin.open((char *)paras["-f"]);
		fin >> info;
		int pos_header = Header::end_header(fin);
		fin.seekg(pos_header, ios::beg);
		
		int file_size = info.count("size") ? atoi(info["size"]) : get_file_lines(fin);
		
		if (size < 0)
			size = file_size;
		else if (size != file_size)
		{
			std::cout << "Error: Dimension indicates in parameters, " << size << " , is different than the dimension found from the input file, " << file_size << ".";
			return;
		}
		paras["-adj-size"] = to_string(size);

		Adj = new SpecMat::LowerTri<PRECISION>(size, fin);
	}

	Adj->form_row_ptr();

	double eps = 1000000;
	double ratio = 0.1;


	info.insert("created", time_stamp());
	info.insert("output_type", "Affinity matrix");
	info.insert("size", to_string(size));
	info.insert("source", paras["-f"]);


	if (paras["-am"] == (String) "REC")
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j <= i; j++)
				if ((*Adj)(i, j) > 0 && eps > (*Adj)(i, j))
					eps = (*Adj)(i, j);
		eps = eps * ratio;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j <= i; j++)
				(*Adj)(i, j) = 1.0 / (eps + (*Adj)(i, j));
		}
	}
	else if ((paras["-am"] == (String) "EXP")) {
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j <= i; j++)
				(*Adj)(i, j) = exp(-(*Adj)(i, j));
		}
	}

	if (!less_print)
	{
		String outname_aff("Affinity");
		outname_aff.make_stdname(paras);
		std::ofstream file_Aff;
		file_Aff.open((char*)outname_aff);
		file_Aff << info;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j <= i; j++)
				file_Aff << (*Adj)(i, j) << "\t";
			file_Aff << "" << std::endl;
		}
	}


	return true;
}

bool comm_driver(SpecMat::LowerTri<PRECISION>* Adj, map<String, String> &paras)
{
	Header_info info;
	if (Adj == nullptr) // Adjacency should be read from file.
	{
		std::fstream fin;
		fin.open((char *)paras["-f"]);
		fin >> info;
		int pos_header = Header::end_header(fin);
		fin.seekg(pos_header, ios::beg);

		int size = atoi(paras["-comm-size"]);
		int file_size = info.count("size") ? atoi(info["size"]) : get_file_lines(fin);
		
		if (size < 0)
			size = file_size;
		else if (size != file_size)
		{
			std::cout << "Error: Number of nodes indicated in parameters, " << size << " , is different than the number of nodes found from the input file, " << file_size << ".";
			return;
		}
		paras["-comm-size"] = to_string(size);

		Adj = new SpecMat::LowerTri<PRECISION>(size, fin);

		if (info.count("distance_type") && !paras.count("-dm")) {
			std::string dm = (char*) info["distance_type"];
			if (dm.find("Unweighted Robinson-Foulds") != std::string::npos || dm.find("URF") != std::string::npos)
				paras["-dm"] = "URF";
			else if (dm.find("Robinson-Foulds") != std::string::npos || dm.find("RF") != std::string::npos)
				paras["-dm"] = "RF";
			else if (dm.find("Matching") != std::string::npos || dm.find("Mat") != std::string::npos)
				paras["-dm"] = "Mat";
			else if (dm.find("SPR") != std::string::npos)
				paras["-dm"] = "SPR";
			else
				paras["-dm"] = "";
		}
	}

	

	if (paras["-comm-type"] == (String) "") {
		if (info.count("node_type")) {
			paras["-comm-type"] = info["node_type"];
		}
		else {
			paras["-comm-type"] = "Unknown_node";
		}
	}

	if (paras["-lm"] == (String) "MANU") {
		// if (community_detection_manually(mat, paras)) {
		// 	cout << "Successfully detected communities of adjacency matrix /" << paras["-f"] << "/ by model: " << paras["-cm"] << ".\n";
		// 	cout << "Lambdas are chosen manually." << endl;
		// }
			cout << "Not supported yet" << endl;
	} else if (paras["-lm"] == (String) "AUTO") {
		if (community_detection_automatically(*Adj, paras)) {
			cout << "Successfully detected communities of adjacency matrix /" << paras["-f"] << "/ by model: " << paras["-cm"] << " with high freq. bound:" << paras["-comm-hf"] << ", low freq. bound:" << paras["-comm-lf"] << "!" << endl;
			cout << "Lambdas are chosen automatically." << endl;
		}
	}
	else {
		cout << "Error! Set parameter -lm only have options: auto, manu\n";
	}
	
	return true;
}


template <class Container_type>
bool task_driver(map<String, String> &paras, Container_type dummy)
{
	if (paras["-trees"] == (String) "none")
	{
		return false;
	}

	bool one_csv = false;
	if (paras["-key"] != (String) "none")
	{
		std::cout << "Loading all parameters from " << paras["-key"] << ".\n";
		read_paras_from_csv(paras["-key"], paras, true);
		one_csv = true;
	}

	if (!one_csv)
	{
		std::cout << "Initiate -trees module.\nLoading parameters from " << paras["-trees"] << ".\n";
		read_paras_from_csv(paras["-trees"], paras, true);
	}
	auto treeobj_ptr = trees_driver(paras, dummy);

	SpecMat::LowerTri<PRECISION>* Adj = nullptr;
	SpecMat::LowerTri<PRECISION>* Dis_NLDR = nullptr;

	if (paras["-output"] == (String) "DIS")
		Adj = new SpecMat::LowerTri<PRECISION>(*treeobj_ptr->get_dis_mat());
	else if (paras["-output"] == (String) "COV")
		Adj = new SpecMat::LowerTri<PRECISION>(*treeobj_ptr->get_cov_mat());

	if (paras["-adj"] != (String) "none")
	{
		if (!one_csv)
		{
			std::cout << "Initiate -adj module.\nLoading parameters from " << paras["-adj"] << ".\n";
			read_paras_from_csv(paras["-adj"], paras, true);
		}
		adj_driver(Adj, paras);
	}

	if (paras["-comm"] != (String) "none")
	{
		if (!one_csv)
		{
			std::cout << "Initiate -comm module.\nLoading parameters from " << paras["-comm"] << ".\n";
			read_paras_from_csv(paras["-comm"], paras, true);
		}
		comm_driver(Adj, paras);
	}

	if (paras["-nldr"] != (String) "none")
	{
		if (!one_csv)
		{
			std::cout << "Initiate -nldr module.\nLoading parameters from " << paras["-nldr"] << ".\n";
			read_paras_from_csv(paras["-nldr"], paras, true);
		}
		Dis_NLDR = new SpecMat::LowerTri<PRECISION>(treeobj_ptr->get_dis_mat());
		auto Cordinates = nldr_driver(Dis_NLDR, paras);
	}
}
