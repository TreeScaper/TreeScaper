#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "zdarray.hpp"
#include "zdtree.hpp"
#include "wstring.hpp"

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


	if(argc > 1 && (String) argv[1] == (String) "-trees")
    {
        String default_paras[4] = {"", "postfix", "taxon", "64"};
        String options[4] =       {"-f", "-post", "-tm", "-bit"};
        
        map<String, String> paras = read_paras(argc, argv, 4, default_paras, options);
        
        trees_driver(paras);
    } 


	
	return 0;
}

bool trees_driver(map<String, String> paras){
	std::clock_t start, end;
	int flag_type;
	int flag_label;

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

	std::cout << "Reading taxa...\n"; 
	start = std::clock();

	TaxonList Taxa;	
	size_t tree_pos = Taxa.ReadTaxa(fname);


	std::cout << "Reading trees...\n"; 


	
	if (flag_type == 0) {
		Taxa.set_bitstr_size((u_8)1);
		Bipartition<u_8> Bipart_8(&Taxa, n_tree);
		TreeSet<u_8> trees_8(&Bipart_8, &Taxa);
		trees_8.ReadTree(fname, tree_pos, flag_label);
		end = std::clock();
		cout << "Read tree time(ms):\t\t" << (end - start) / 1000.0 << "\n";
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		trees_8.compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(ms):\t" << (end - start) / 1000.0 << "\n";
		Bipart_8.print_summary(trees_8.get_size());
		std::cout << "Computing bipartition matrix...\n"; 
		start = std::clock();
		TreeObjects<u_8> treesobj_8(&(trees_8));
		treesobj_8.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Form Bipartition Matrix time(ms):\t" << (end - start) / 1000.0 << "\n";
		String outname_Bipartcnt("Bipartition_Count");
		outname_Bipartcnt.make_stdname(paras);
		ofstream fout;
		fout.open((char *) outname_Bipartcnt);
		Bipart_8.print_Bipart(fout);
		fout.close();
		trees_8.release_tree();
	}
	else if (flag_type == 1) {
		Taxa.set_bitstr_size((u_32)1);
		Bipartition<u_32> Bipart_32(&Taxa, n_tree);
		TreeSet<u_32> trees_32(&Bipart_32, &Taxa);
		trees_32.ReadTree(fname, tree_pos, flag_label);
		end = std::clock();
		cout << "Read tree time(ms):\t\t" << (end - start) / 1000.0 << "\n";
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		trees_32.compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(ms):\t" << (end - start) / 1000.0 << "\n";
		Bipart_32.print_summary(trees_32.get_size());
		std::cout << "Computing bipartition matrix...\n"; 
		start = std::clock();
		TreeObjects<u_32> treesobj_32(&(trees_32));
		treesobj_32.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Form Bipartition Matrix time(ms):\t" << (end - start) / 1000.0 << "\n";
		String outname_Bipartcnt("Bipartition_Count");
		outname_Bipartcnt.make_stdname(paras);
		ofstream fout;
		fout.open((char *) outname_Bipartcnt);
		Bipart_32.print_Bipart(fout);
		fout.close();
		trees_32.release_tree();
	}
	else {
		Taxa.set_bitstr_size((u_64)1);
		Bipartition<u_64> Bipart_64(&Taxa, n_tree);
		std::cout << "Bipartition initialized...\n"; 
		TreeSet<u_64> trees_64(&Bipart_64, &Taxa);
		trees_64.ReadTree(fname, tree_pos, flag_label);
		end = std::clock();
		cout << "Read tree time(ms):\t\t" << (end - start) / 1000.0 << "\n";
		std::cout << "Computing bipartition...\n"; 
		start = std::clock();
		trees_64.compute_Bipart();
		end = std::clock();
		cout << "Compute Bipartition time(ms):\t" << (end - start) / 1000.0 << "\n";
		Bipart_64.print_summary(trees_64.get_size());
		std::cout << "Computing bipartition matrix...\n"; 
		start = std::clock();
		TreeObjects<u_64> treesobj_64(&(trees_64));
		treesobj_64.Compute_Bipart_Matrix();
		end = std::clock();
		cout << "Form Bipartition Matrix time(ms):\t" << (end - start) / 1000.0 << "\n";
		String outname_Bipartcnt("Bipartition_Count");
		outname_Bipartcnt.make_stdname(paras);
		ofstream fout;
		fout.open((char *) outname_Bipartcnt);
		Bipart_64.print_Bipart(fout);
		fout.close();
		trees_64.release_tree();
	}
	return true;

}
