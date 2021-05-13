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
#include "version.hpp"
#include "cra.hpp"

using namespace std;

typedef unsigned char u_8;
typedef unsigned int u_32;
typedef unsigned long long u_64;

String make_stdname(String s, std::map<String, String> &paras) {
	String Ans = paras["-path"];
	Ans += s;
	if (paras["-post"] != String("none")) {
		Ans += "_";
		// if (paras["-post"] != String("time"))
		// 	Ans += paras["-post"];
		// else
		// 	Ans += "_post"
		Ans += paras["-post"];
	}
	Ans += ".out";
	return Ans;
}

String get_path(String fname) {
	String path = fname;
	std::string temp = (char*)path;
	size_t loc_slash = temp.find_last_of('/');
	if (loc_slash != string::npos)
		path = path(0, loc_slash + 1);
	else
		path = "";
	return path;
}

bool trees_driver(map<String, String> paras);

map<String, String> read_paras(int argc, char* argv[], int key_size, String* default_paras, String* options){
	cout << "loading all parameters" << endl;
	for(int i = 1; i < argc; i++){
		for(int j = 0; j < key_size; j++){
			if((String) argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-'){
				i++;
                default_paras[j] = argv[i];
                break;
			}
		}
	}
	map<String, String> paras;
	for(int i = 0; i < key_size; i++)
        paras[options[i]] = default_paras[i];
	paras["-path"] = get_path(paras[String("-f")]);
	return paras;
}

int main(int argc, char* argv[]) {
	if (argc == 1) {
		cerr << "No arguments supplied." << endl;
		return 1;
	}
	if ((String) argv[1] == (String) "-inference") {
		String default_paras[] = {"", "", ""};
		String options[] = {"-f", "-p", "-cl"};

		map<String, String> paras = read_paras(argc, argv, sizeof(options)/sizeof(String), default_paras, options);

		string inputfile = string((char*)paras["-f"]);
		string paramfile = string((char*)paras["-p"]);
		cra_log_level = string((char*)paras["-cl"]) == "debug" ? DEBUG : NONE;

		CRAHandle crahandle;
		if (crahandle.submit_jobs(inputfile, paramfile) == false) {
			return 1;
		}
		return 0;

    } else if(argc > 1 && (String) argv[1] == (String) "-trees") {
        String default_paras[] = {"", 		"postfix", 	"taxon", 	"64", 	"0",	"1",	"Cov"};
        String options[] =       {"-f", "	-post", 	"-tm", 		"-bit", "-r",	"-w", 	"-o"};
        
		map<String, String> paras = read_paras(argc, argv, sizeof(options)/sizeof(String), default_paras, options);
        
        trees_driver(paras);
	} else if ((String) argv[1] == (String) "-version" || (String) argv[1] == (String) "-v") {
		cout << program_version << endl;
	}

	return 0;
}

bool trees_driver(map<String, String> paras){
	std::clock_t start, end;
	int flag_type;
	int flag_label;
	bool rooted = (paras["-r"] == (String) "1");
	bool weighted = (paras["-w"] == (String) "1");

	std::cout << (rooted ? "rooted" : "unrooted") << "\t" << (weighted ? "weighted" : "unweighted") << '\n';

	if (paras["-tm"] == (String) "taxon")
		flag_label = 2;
	else
		flag_label = 0;

	if (paras["-bit"] == (String) "8")
		flag_type = 0;
	else if (paras["-bit"] == (String) "32")
		flag_type = 1;
	else if (paras["-bit"] == (String) "64")
		flag_type = 2;
	else{
		std::cout << "Bitstring container size not supported. Please choose value from {8, 32, 64} for key `-bit'.\n";
		return false;
	}
	// flag_label:	0:	labelled with integer
	//				1:	labelled with integer from different normalization
	//				2:	labelled with taxon name
	int n_tree = 100;
	std::string fname((char *) paras["-f"]);

	std::cout << "Reading taxa...\n" << endl;
	start = std::clock();

	TaxonList Taxa;	
	size_t tree_pos = Taxa.ReadTaxa(fname);


	std::cout << "Reading trees...\n"; 


	
	if (flag_type == 0) {
		Taxa.set_bitstr_size((u_8)1);
		Bipartition<u_8> Bipart_8(&Taxa, n_tree);
		std::cout << "Bipartition initialized...\n"; 
		TreeSet<u_8> trees_8(&Bipart_8, &Taxa);	
		trees_8.ReadTree(fname, tree_pos, flag_label, rooted, weighted);
		end = std::clock();
		cout << "Read tree time(s):\t\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		TreeSetObjects<u_8> treesobj_8(&(trees_8), rooted, weighted);
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		treesobj_8.Compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(s):\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		cout << "Distinct bipartition number:\t" << treesobj_8.get_unique_bipart_num() << "\n\n";

		ofstream fout;
		// // Print Bipartition List
		// String outname_Bipartcnt("Bipartition_Count");
		// outname_Bipartcnt.make_stdname(paras);
		// fout.open((char *) outname_Bipartcnt);
		// treesobj_8.print_Bipart_List(fout);
		// fout.close();
		// cout << "Sucessfully printed bipartition list info in Bipartition_count_" << paras["-post"] << ".out file.\n\n";


		std::cout << "Computing Bipart2Tree sparse matrix...\n";
		start = std::clock();
		treesobj_8.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Compute Bipart2Tree sparse matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

		// Print Bipartition2Tree Sparse Matrix
		String outname_Bipart("Bipartition");
		outname_Bipart.make_stdname(paras);
		fout.open((char *) outname_Bipart);
		treesobj_8.print_Bipart2Tree_Matrix(fout, RCVLIST);
		fout.close();
		cout << "Sucessfully printed bipartition-to-tree sparse matrix in Bipartition_" << (char*)paras["-post"] << ".out file.\n\n";



		if (paras["-o"] == (String) "Cov" || paras["-o"] == (String) "cov" || paras["-o"] == (String) "covariance"){
			std::cout << "Computing covariance matrix...\n";
			start = std::clock();
			treesobj_8.Compute_Covariance_Matrix();
			end = std::clock();
			cout << "Compute covariance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Cova("Covariance");
			outname_Cova.make_stdname(paras);
			fout.open((char *) outname_Cova);
			treesobj_8.print_Covariance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed covariance matrix in Covariance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		else if (paras["-o"] == (String) "Dis" || paras["-o"] == (String) "dis" || paras["-o"] == (String) "distance"){
			std::cout << "Computing RF-distance matrix...\n";
			start = std::clock();
			treesobj_8.Compute_RF_Distance_Matrix();
			end = std::clock();
			cout << "Compute RF-distance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Dist("Distance");
			outname_Dist.make_stdname(paras);
			fout.open((char *) outname_Dist);
			treesobj_8.print_Distance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed distance matrix in Distance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		treesobj_8.print_summary();
		trees_8.release_tree();
	}
	else if (flag_type == 1) {
		Taxa.set_bitstr_size((u_32)1);
		Bipartition<u_32> Bipart_32(&Taxa, n_tree);
		std::cout << "Bipartition initialized...\n"; 
		TreeSet<u_32> trees_32(&Bipart_32, &Taxa);	
		trees_32.ReadTree(fname, tree_pos, flag_label, rooted, weighted);
		end = std::clock();
		cout << "Read tree time(s):\t\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		TreeSetObjects<u_32> treesobj_32(&(trees_32), rooted, weighted);
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		treesobj_32.Compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(s):\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		cout << "Distinct bipartition number:\t" << treesobj_32.get_unique_bipart_num() << "\n\n";

		ofstream fout;
		// // Print Bipartition List
		// String outname_Bipartcnt("Bipartition_Count");
		// outname_Bipartcnt.make_stdname(paras);
		// fout.open((char *) outname_Bipartcnt);
		// treesobj_32.print_Bipart_List(fout);
		// fout.close();
		// cout << "Sucessfully printed bipartition list info in Bipartition_count_" << paras["-post"] << ".out file.\n\n";


		std::cout << "Computing Bipart2Tree sparse matrix...\n";
		start = std::clock();
		treesobj_32.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Compute Bipart2Tree sparse matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

		// Print Bipartition2Tree Sparse Matrix
		String outname_Bipart("Bipartition");
		outname_Bipart.make_stdname(paras);
		fout.open((char *) outname_Bipart);
		treesobj_32.print_Bipart2Tree_Matrix(fout, RCVLIST);
		fout.close();
		cout << "Sucessfully printed bipartition-to-tree sparse matrix in Bipartition_" << (char*)paras["-post"] << ".out file.\n\n";



		if (paras["-o"] == (String) "Cov" || paras["-o"] == (String) "cov" || paras["-o"] == (String) "covariance"){
			std::cout << "Computing covariance matrix...\n";
			start = std::clock();
			treesobj_32.Compute_Covariance_Matrix();
			end = std::clock();
			cout << "Compute covariance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Cova("Covariance");
			outname_Cova.make_stdname(paras);
			fout.open((char *) outname_Cova);
			treesobj_32.print_Covariance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed covariance matrix in Covariance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		else if (paras["-o"] == (String) "Dis" || paras["-o"] == (String) "dis" || paras["-o"] == (String) "distance"){
			std::cout << "Computing RF-distance matrix...\n";
			start = std::clock();
			treesobj_32.Compute_RF_Distance_Matrix();
			end = std::clock();
			cout << "Compute RF-distance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Dist("Distance");
			outname_Dist.make_stdname(paras);
			fout.open((char *) outname_Dist);
			treesobj_32.print_Distance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed distance matrix in Distance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		treesobj_32.print_summary();
		trees_32.release_tree();
	}
	else {
		Taxa.set_bitstr_size((u_64)1);
		Bipartition<u_64> Bipart_64(&Taxa, n_tree);
		std::cout << "Bipartition initialized...\n"; 
		TreeSet<u_64> trees_64(&Bipart_64, &Taxa);	
		trees_64.ReadTree(fname, tree_pos, flag_label, rooted, weighted);
		end = std::clock();
		cout << "Read tree time(s):\t\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		TreeSetObjects<u_64> treesobj_64(&(trees_64), rooted, weighted);
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		treesobj_64.Compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(s):\t" << (end - start)/ (double) CLOCKS_PER_SEC << "\n";
		cout << "Distinct bipartition number:\t" << treesobj_64.get_unique_bipart_num() << "\n\n";

		ofstream fout;
		// Print Bipartition List
		// String outname_Bipartcnt("Bipartition_Count");
		// outname_Bipartcnt.make_stdname(paras);
		// fout.open((char *) outname_Bipartcnt);
		// treesobj_64.print_Bipart_List(fout);
		// fout.close();
		// cout << "Sucessfully printed bipartition list info in Bipartition_count_" << paras["-post"] << ".out file.\n\n";


		std::cout << "Computing Bipart2Tree sparse matrix...\n";
		start = std::clock();
		treesobj_64.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Compute Bipart2Tree sparse matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

		// Print Bipartition2Tree Sparse Matrix
		String outname_Bipart("Bipartition");
		outname_Bipart.make_stdname(paras);
		fout.open((char *) outname_Bipart);
		treesobj_64.print_Bipart2Tree_Matrix(fout, RCVLIST);
		fout.close();
		cout << "Sucessfully printed bipartition-to-tree sparse matrix in Bipartition_" << (char*)paras["-post"] << ".out file.\n\n";



		if (paras["-o"] == (String) "Cov" || paras["-o"] == (String) "cov" || paras["-o"] == (String) "covariance"){
			std::cout << "Computing covariance matrix...\n";
			start = std::clock();
			treesobj_64.Compute_Covariance_Matrix();
			end = std::clock();
			cout << "Compute covariance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Cova("Covariance");
			outname_Cova.make_stdname(paras);
			fout.open((char *) outname_Cova);
			treesobj_64.print_Covariance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed covariance matrix in Covariance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		else if (paras["-o"] == (String) "Dis" || paras["-o"] == (String) "dis" || paras["-o"] == (String) "distance"){
			std::cout << "Computing RF-distance matrix...\n";
			start = std::clock();
			treesobj_64.Compute_RF_Distance_Matrix();
			end = std::clock();
			cout << "Compute RF-distance matrix time(s):\t" << (end - start) / (double)CLOCKS_PER_SEC << "\n";

			// Print Covariance Matrix
			String outname_Dist("Distance");
			outname_Dist.make_stdname(paras);
			fout.open((char *) outname_Dist);
			treesobj_64.print_Distance_Matrix(fout);
			fout.close();
			cout << "Sucessfully printed distance matrix in Distance_" << (char*)paras["-post"] << ".out file.\n\n";
		}
		treesobj_64.print_summary();
		trees_64.release_tree();
	}
	return true;

}
