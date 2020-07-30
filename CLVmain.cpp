
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

//CLVmain.h
//        command line version   March/11/2010
//                               by whuang
//        updated 2016-2-29      by whuang 
/*
*/

#include "randgen.h"
#include <iostream>
#include <fstream>
#include "warray.cpp"
#include "wstring.h"
#include "wmatrix.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include "wfile.h"
#include "wNLDR.h"
#include "wDimEst.h"
#include "Trees.h"
#include "TreeOPE.h"
#include <map>

using namespace std;

void Compute_Community(Trees *TreesData, map<String, String> &paras, String memorydata);
void Compute_Covariance(Trees *TreesData, map<String, String> &paras);
void Compute_Affinity(Trees *TreesData, map<String, String> &paras, String memorydata);
void Compute_Distance(Trees *TreesData, map<String, String> &paras);
void Compute_Consensus_Tree(Trees *TreesData, map<String, String> &paras);
void Compute_BipartMatrix(Trees *TreesData, map<String, String> &paras);
void trees_driver(map<String, String> &paras);
void dimest_driver(String fname, String Est, String Init, String para_fname);
void driver(String fname, String ftype, String dim, String cost, String algo, String init_md, String flag, long seed, String para_fname);
void nldr_driver(map<String, String> &paras);
void aff_driver(map<String, String> &paras);
String get_path(String fname);

int main(int argc, char* argv[])
{
    if(argc > 1 && (String) argv[1] == (String) "-dimest")
    {
        String dimest_default_paras[4] = {"ccatest.out", "NN_DIM", "DIS", "dimest_parameters.csv"};
        String dimest_options[4] = {"-f", "-e", "-i", "-p"};

        for(int i = 1; i < argc; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                if((String) argv[i] == dimest_options[j] && i + 1 < argc && argv[i + 1][0] != '-')
                {
                    i++;
                    dimest_default_paras[j] = argv[i];
                    break;
                }
            }
        }
        dimest_driver(dimest_default_paras[0], dimest_default_paras[1], dimest_default_paras[2], dimest_default_paras[3]);
        return 0;
    } else
    if(argc > 1 && (String) argv[1] == (String) "-nldr")
    {
        String default_paras[9] = {"trajectory1.out", "COR", "2", "CCA", "STOCHASTIC", "RAND", "", "1", "nldr_parameters.csv"};
        String options[9] = {"-f", "-t", "-d", "-c", "-a", "-i", "-o", "-s", "-p"};
    
        for(int i = 1; i < argc; i++)
        {
            for(int j = 0; j < 9; j++)
            {
                if((String) argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-')
                {
                    i++;
                    default_paras[j] = argv[i];
                    break;
                }
            }
        }
		map<String, String> paras;
		for (int i = 0; i < 9; i++)
		{
			paras[options[i]] = default_paras[i];
		}
		paras["-path"] = get_path(paras[String("-f")]);
        cout << "loading all parameters" << endl;
        //driver(default_paras[0], default_paras[1], default_paras[2], default_paras[3], default_paras[4], default_paras[5], default_paras[6], atoi(default_paras[7]), default_paras[8]);
		nldr_driver(paras);
		return 0;
    } else
    if(argc > 1 && (String) argv[1] == (String) "-trees")
    {
        String default_paras[24] = {"nuctrees.txt", "0", "0", "Community", "list", 
                                    "", "Majority", "Newick", "URF", "Exp",
                                    "Covariance", "CNM", "1", "0", "1", "0", "1", "0", "1", "0", "1", "0", "auto", "Trees"};
        String options[24] =       {"-f", "-w", "-r", "-o", "-bfm", 
                                    "-if", "-ct", "-cfm", "-dm", "-am",
                                    "-t", "-cm", "-lp", "-lps", "-lpe", "-lpiv", "-ln", "-lns", "-lne", "-lniv", "-hf", "-lf", "-lm", "-ft"};
        
        for(int i = 1; i < argc; i++)
        {
            for(int j = 0; j < 24; j++)
            {
                if((String) argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-')
                {
                    i++;
                    default_paras[j] = argv[i];
                    break;
                }
            }
        }
        map<String, String> paras;
        for(int i = 0; i < 24; i++)
        {
            paras[options[i]] = default_paras[i];
        }
		paras["-path"] = get_path(paras[String("-f")]);
        cout << "loading all parameters" << endl;
        
        trees_driver(paras);
    } 
	else if (argc > 1 && (String)argv[1] == (String) "-aff") {
		String default_paras[2] = { "", "Exp" };
		String options[2] = { "-f", "-am" };
		for (int i = 1; i < argc; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if ((String)argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-')
				{
					i++;
					default_paras[j] = argv[i];
					break;
				}
			}
		}
		map<String, String> paras;
		for (int i = 0; i < 2; i++)
		{
			paras[options[i]] = default_paras[i];
		}
		paras["-path"] = get_path(paras[String("-f")]);
		cout << "loading all parameters" << endl;

		aff_driver(paras);
		return 0;
	}else
    if(argc == 2 && ((String) argv[1] == (String) "-h" || (String) argv[1] == (String) "-help"))
    {
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "HELP OF COMMAND LINE VERSION\n";
        std::cout << "(Details also can be found in the user manual.)\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "There are three general run modesfor the command-line version of TreeScaper(“CLVTreeScaper”). Using one of these three flags as the first command-line argument sets the mode (e.g., CLVTreeScaper –trees).\n";
        std::cout << "\n";
        std::cout << "(1) -trees\n";
        std::cout << "\n";
        std::cout << "In this mode, users can compute a majority rule/strict consensus tree, distance matrix,bipartition matrix, covariance matrix, affinity matrix, or detect communitiesin an affinity or covariance network. Relevant arguments include:\n";
        std::cout << "\n";
        std::cout << "-f:Provide the name of the file that contains the data\n";
        std::cout << "-ft: The file type. Options are:\n";
        std::cout << "'Trees': the file contains trees. The tree format can be either Newick or Nexus.\n";
        std::cout << "'Dist': the file contains distance matrix which can be used to compute affinity matrix or communities.\n";
        std::cout << "'Cova': the file contains covariance matrix which can be used to compute communities.\n";
        std::cout << "-w: Indicate whether trees are weighted. Options are:\n";
        std::cout << "'1': weighted\n";
        std::cout << "'0': unweighted\n";
        std::cout << "-r: Indicate whether trees are rooted. Options are:\n";
        std::cout << "'1': rooted\n";
        std::cout << "'0': unrooted\n";
        std::cout << "-o: this option is used to indicate what output the user is interested in. Options are:\n";
        std::cout << "'BipartMatrix',\n";
        std::cout << "'Consensus',\n";
        std::cout << "'Dist',\n";
        std::cout << "'Affinity',\n";
        std::cout << "'Covariance',\n";
        std::cout << "'Community'.\n";
        std::cout << "\n";
        std::cout << "When outputting a bipartition matrix (-o BipartMatrix):\n";
        std::cout << "\n";
        std::cout << "-bfm: Bipartition matrix output type. Options are:\n";
        std::cout << "'list': output sparse matrix in the form (row, column, value)\n";
        std::cout << "'matrix': output as if it is a full matrix\n";
        std::cout << "\n";
        std::cout << "When computing a majority-rule or strict consensus tree (-o Consensus), use the -if, -ct, and/or-cfm flags:\n";
        std::cout << "\n";
        std::cout << "-if: The name of a list file. Consensus tree computations will only consider the trees indicated in the file.\n";
        std::cout << "-ct: The type of consensus tree to be computed. Options are:\n";
        std::cout << "'Majority': Majority consensus tree\n";
        std::cout << "'Strict': Strict consensus tree\n";
        std::cout << "-cfm: Format of the consensus tree file. Options are:\n";
        std::cout << "'Newick'\n";
        std::cout << "'Nexus'\n";
        std::cout << "\n";
        std::cout << "When computing a distance matrix (-o Dist):\n";
        std::cout << "\n";
        std::cout << "-dm: Indicates the distance metric. Options are:\n";
        std::cout << "	'URF': Unweighted Robinson-Foulds distance\n";
        std::cout << "	'RF': Weighted Robinson-Foulds distance\n";
        std::cout << "	'Mat': Matching distance\n";
        std::cout << "	‘SPR’: Subtree-Prune-Regraft\n";
        std::cout << "\n";
        std::cout << "When computing an affinity matrix (-o Affinity):\n";
        std::cout << "\n";
        std::cout << "-dm: Indicates the distance metric. Options are:\n";
        std::cout << "	'URF': Unweighted Robinson Foulds distance\n";
        std::cout << "	'RF': Weighted Robinson Foulds distance\n";
        std::cout << "	'Mat': Matching distance\n";
        std::cout << "	'SPR': Subtree-Prune and Regraft \n";
        std::cout << "-am: Indicates the distance to affinity transformation. Options are:\n";
        std::cout << "	'Rec': Reciprocal\n";
        std::cout << "	'Exp': Exponential\n";
        std::cout << "\n";
        std::cout << "When detecting communities (-o Community):\n";
        std::cout << "\n";
        std::cout << "-t: Target matrix used to compute communities. Options are:\n";
        std::cout << "'Affinity': affinity matrix\n";
        std::cout << "'Covariance': covariance matrix\n";
        std::cout << "-cm: Model used to compute communities. Options are:\n";
        std::cout << "'CNM': Configuration Null Model\n";
        std::cout << "'CPM': Constant Potts Model\n";
        std::cout << "'ERNM': Erdos-Rényi Null Model\n";
        std::cout << "'NNM': No Null Model\n";
        std::cout << "\n";
        std::cout << "-lm: Method of plateau detection. Options are:\n";
        std::cout << "'auto': automatically choose lambdas and find plateaus\n";
        std::cout << "'manu': specify intervals by users to find plateaus\n";
        std::cout << "\n";
        std::cout << "The following flags are used to specify values of lambda for manual searches:\n";
        std::cout << "\n";
        std::cout << "-lp: Specify a fixed value of λ+.Must be between 0 and 1.Used when -lpiv is zero (see below).\n";
        std::cout << "-lps, -lpe, -lpiv: Starting, ending, and sampling intervals for λ+. Used to explore a range of possible values for λ+.\n";
        std::cout << "\n";
        std::cout << "-ln: Specify a fixed value ofλ-. Must be between 0 and 1. Used when -lniv is zero (see below).\n";
        std::cout << "-lns, -lne, -lniv: Starting, ending, and sampling intervals for λ-. Used to explore a range of possible values for λ-.\n";
        std::cout << "\n";
        std::cout << "Note: Either λ+ or λ- must be fixed, because plateau detection is undefined when both vary.\n";
        std::cout << "\n";
        std::cout << "-hf: Frequency upper bound. A number between 0 and 1.Nodes with frequencies above this value are ignored.\n";
        std::cout << "-lf: Frequency lower bound. A number between 0 and 1. Nodes with frequencies below this value are ignored.\n";
        std::cout << "\n";
        std::cout << "Examples of command-line runs\n";
        std::cout << "Options specified by the are given inside braces. When specific alternatives are available, they are separated by commas (e.g., {option1,option2}). When numbers can be specified anywhere in a continuous range, the bounds of the range are separated by a dash (e.g., {0-1}).\n";
        std::cout << "\n";
        std::cout << "Compute a Bipartition Matrix:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o BipartMatrix -bfm \n";
        std::cout << "{list,matrix}\n";
        std::cout << "\n";
        std::cout << "Compute a Consensus Tree:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Consensus -if IndicesFileName -ct {Majority,Strict} -cfm {Newick,Nexus}\n";
        std::cout << "\n";
        std::cout << "Compute Distance Matrix:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Dist -dm {URF,RF,Mat,SPR}\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "Compute Affinity Matrix:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Affinity -dm {URF,RF,Mat,SPR}\n";
        std::cout << "-am {Exp,Rec}\n";
        std::cout << "\n";
        std::cout << "Compute Covariance Matrix:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Covariance\n";
        std::cout << "\n";
        std::cout << "Compute Communities with λ+ Fixed:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Community -t {Affinity,Covariance} -cm {CNM,CPM,ERNM,NNM} -lm manu -lp {AnyNumber} -lns {AnyNumber} -lne {AnyNumber} -lniv {AnyNumber} -hf {0-1} -lf {0-1}\n";
        std::cout << "\n";
        std::cout << "Compute Communities with λ- Fixed:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Community -t {Affinity,Covariance} -cm {CNM,CPM,ERNM,NNM} -lm manu -ln {AnyNumber} -lps {AnyNumber} -lpe {AnyNumber} -lpiv {AnyNumber} -hf {0-1} -lf {0-1}\n";
        std::cout << "\n";
        std::cout << "Compute Communities with Automatically Chosen Lambdas:\n";
        std::cout << "./CLVTreeScaper -trees -f {trees.txt} -ft Trees -w {1,0} -r {1,0} -o Community -t {Affinity/Covariance} -cm {CNM/CPM/ERNM/NNM} -lm auto -hf {0-1} -lf {0-1}\n";
        std::cout << "\n";
        std::cout << "Load Distances and Compute Affinity Matrix:\n";
        std::cout << "./CLVTreeScaper -trees -f {dist.txt} -ft Dist -o Affinity -am {Exp,Rec}\n";
        std::cout << "\n";
        std::cout << "Load Distances and Compute Affinity Communities Automatically:\n";
        std::cout << "./CLVTreeScaper -trees -f {dist.txt} -ft Dist -o Community -t Affinity -am {Exp,Rec} -cm {CNM,CPM,ERNM,NNM} -lm auto -hf {0-1} -lf {0-1}\n";
        std::cout << "\n";
        std::cout << "Load Covariances and Compute Communities using Automatic Search on Lambda\n";
        std::cout << "./CLVTreeScaper -trees -f {cova.txt} -ft Cova -o Community -cm {CNM,CPM,ERNM,NNM} -lm auto -hf {0-1} -lf {0-1}\n";
        std::cout << "\n";
        std::cout << "(2)–nldr\n";
        std::cout << "\n";
        std::cout << "In this mode, users can project trees into lower dimensional space using non-linear dimensionality reduction (NLDR). Relevant arguments include:\n";
        std::cout << "\n";
        std::cout << "-f: Name of the file containing distance data.\n";
        std::cout << "-t: The type of distances contained in the file. Options are:\n";
        std::cout << "'DIS':a lower triangle matrix of original distances\n";
        std::cout << "'COR':low-dimensional Euclidean coordinates (if already computed).\n";
        std::cout << "-d: The desired dimension of the Euclidean representation (usually 1, 2, or 3).\n";
        std::cout << "-c: The chosen cost function. Options are:\n";
        std::cout << "‘CLASSIC_MDS’\n";
        std::cout << "‘KRUSKAL1’\n";
        std::cout << "‘NORMALIZED’\n";
        std::cout << "‘SAMMON’\n";
        std::cout << "‘CCA’\n";
        std::cout << "-a: The chosen NLDRalgorithm. Options are:\n";
        std::cout << "‘LINEAR_ITERATION’\n";
        std::cout << "‘MAJORIZATION’\n";
        std::cout << "‘GAUSS_SEIDEL’\n";
        std::cout << "‘STOCHASTIC’\n";
        std::cout << "-i: The methodfor generating initial Euclidean coordinates. Options are:\n";
        std::cout << "'RAND': Randomly choose coordinates for each point.\n";
        std::cout << "'CLASSIC_MDS': Generate initial coordinates using classic multi-dimensional scaling (MDS).\n";
        std::cout << "-o: The suffix for output file names.\n";
        std::cout << "-s: A random seed, if initial coordinates are generated randomly.\n";
        std::cout << "\n";
        std::cout << "Example command:\n";
        std::cout << "./CLVTreeScaper -nldr -f {test.out} -t {DIS,COR} -d {1,2,3,…}\n";
        std::cout << "-c {CLASSIC_MDS,KRUSKAL1,NORMALIZED,SAMMON,CCA}\n";
        std::cout << "-a {LINEAR_ITERATION,MAJORIZATION,GAUSS_SEIDEL,STOCHASTIC} -i {RAND,CLASSIC_MDS} -o {run1} -s 1\n";
        std::cout << "\n";
        std::cout << "(3), -dimest\n";
        std::cout << "\n";
        std::cout << "In this mode, users can estimate the intrinsic dimensionality of their data. This estimate can help in deciding on an appropriate number of dimensions to use when performing NLDR projections.\n";
        std::cout << "\n";
        std::cout << "-f: Name of the file containing distance data.\n";
        std::cout << "-i: The type of distances contained in the file. Options are:\n";
        std::cout << "'DIS': a lower triangle matrix of original distances\n";
        std::cout << "'COR': low-dimensional Euclidean coordinates (if already computed).\n";
        std::cout << "-e: The chosen estimator. Options are:\n";
        std::cout << "'CORR_DIM': correlation dimension estimator\n";
        std::cout << "'NN_DIM': nearest neighbor estimator\n";
        std::cout << "'MLE_DIM': maximum likelihood estimator\n";
        std::cout << "\n";
        std::cout << "Example command:\n";
        std::cout << "./CLVTreeScaper -dimest -f {test.out} -i {DIS,COR} -e {CORR_DIM,NN_DIM,MLE_DIM}\n";
    } else
    {
        cout << "Error: input error!" << endl;
    }
    
    return 0;
}

void trees_driver(map<String, String> &paras)
{
    String fname = paras["-f"];
    string stdfname = (char *) fname;
    
//     std::cout << "\nTest1: fname = " << stdname << std::endl;
        
    File file(fname);
    if(! file.is_open())
    {
        cout << "Error: can not open the data file!" << endl;
        return;
    }

// 	ifstream file;
//     file.open(stdfname);
//     if(!file.is_open())
//     {
//         cout << "Error: can not open the data file!" << endl;
//         return;
//     }
    
    Trees *TreesData = new Trees;
    if(paras["-ft"] == (String) "Trees")
    {
        TreesData->initialTrees(stdfname);
        if(paras["-w"] == (String) "1")
        {
            TreesData->Settreeweighttype(true);
            std::cout << "weighted ";
        } else
        {
            TreesData->Settreeweighttype(false);
            std::cout << "unweighted ";
        }
        if(paras["-r"] == (String) "1")
        {
            TreesData->Settreeroottype(true);
            std::cout << "rooted ";
        } else
        {
            TreesData->Settreeroottype(false);
            std::cout << "unrooted ";
        }
        std::cout << "trees" << std::endl;
        TreesData->ReadTrees();
        TreesData->compute_numofbipart();
        cout << "Successfully read " << TreesData->Get_n_trees() << " trees from file: " << fname << "." << endl;
    }
    Compute_BipartMatrix(TreesData, paras);
    
    delete TreesData;
}

void Compute_BipartMatrix(Trees *TreesData, map<String, String> &paras)
{
    if(paras["-ft"] == (String) "Trees")
    {
        TreesData->Compute_Hash();
        TreesData->Compute_Bipart_Matrix();
    
        cout << "Successfully computed bipartitation matrix." << endl;
        
        String fname = paras["-f"];
        string stdfname = (char *) fname;
        string namebipartmatrix;
        ofstream outBipartMatrix;

		String** info = new String*[10];
		for (int i = 0; i < 10; i++)
			info[i] = new String[2];
		info[0][0] = "created";			info[0][1] = time_stamp();
		info[1][0] = "output_type";		info[1][1] = "Bipartition matrix";
		info[2][0] = "output_format";	info[2][1] = "";
		info[3][0] = "source";			info[3][1] = fname;


		if (paras["-bfm"] == (String) "list")
		{
			string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, (String) "List format");
			outBipartMatrix.open(namebipartmatrix.c_str());
			cout << "Output bipartition matrix in list format to " << namebipartmatrix << endl;
			TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
			outBipartMatrix.close();

			namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname);
			File file_Bipart(namebipartmatrix.c_str());
			file_Bipart.clean();
			info[2][1] = "List format";
			file_Bipart.insert_header(info, 4);
			file_Bipart.close();
			outBipartMatrix.open(namebipartmatrix.c_str());
			TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
		}
		else if (paras["-bfm"] == (String) "matrix")
		{
			string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, (String) "Matrix format");
			outBipartMatrix.open(namebipartmatrix.c_str());
			cout << "Output bipartition matrix in matrix format to " << namebipartmatrix << endl;
			TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
			outBipartMatrix.close();

			namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname);
			File file_Bipart(namebipartmatrix.c_str());
			file_Bipart.clean();
			info[2][1] = "Matrix format";
			file_Bipart.insert_header(info, 4);
			file_Bipart.close();
			outBipartMatrix.open(namebipartmatrix.c_str());
			TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
			}
		else
		{
			cout << "Error: setting of -bfm is not correct. Unable to output bipartition matrix" << endl;
			return;
		}
    }
      
    if(paras["-o"] == (String) "Consensus")
        Compute_Consensus_Tree(TreesData, paras);
    else
    if(paras["-o"] == (String) "Dist" || paras["-o"] == (String) "Affinity" || (paras["-o"] == (String) "Community" && paras["-t"] == (String) "Affinity"))
        Compute_Distance(TreesData, paras);
    else
    if(paras["-o"] == (String) "Covariance" || (paras["-o"] == (String) "Community" && paras["-t"] == (String) "Covariance"))
        Compute_Covariance(TreesData, paras);
}

void Compute_Covariance(Trees *TreesData, map<String, String> &paras)
{
	String** info = new String*[10];
	for (int i = 0; i < 10; i++)
		info[i] = new String[2];
	info[0][0] = "created";			info[0][1] = time_stamp();
	info[1][0] = "output_type";		info[1][1] = "Covariance matrix";
	info[2][0] = "label_type";		info[2][1] = "";
	info[3][0] = "source";			info[3][1] = paras["-f"];

    if(paras["-ft"] == (String) "Trees")
    {
        TreesData->Compute_Bipart_Covariance();
        cout << "Successfully computed covariance matrix of bipartition." << endl;     
        
        string outCovaName = TreesData->make_DISToutput_name("Covariance Matrix");
        TreesData->print_matrix("Covariance Matrix", outCovaName);
        cout << "Successfully printed Covariance Matrix matrix!" << endl;

		String outCovaName2 = paras["-path"];
		outCovaName2 += "Covariance.out";


		File file_Cova(outCovaName2);
		file_Cova.clean();
		info[2][1] = "Trees";
		file_Cova.insert_header(info, 4);
		file_Cova.close();

		TreesData->print_matrix2("Covariance Matrix", (char *)outCovaName2);

        if(paras["-o"] == (String) "Community")
        {
            Compute_Community(TreesData, paras, "Covariance Matrix");
        }
    }
    if(paras["-ft"] == (String) "Cova")
    {
        String fname = paras["-f"];
        string stdfname = (char *) fname;
        TreesData->load_covariancefile(stdfname);

		String outCovaName2 = paras["-path"];
		outCovaName2 += "Covariance.out";


		File file_Cova(outCovaName2);
		file_Cova.clean();
		info[2][1] = "Bipartition";
		file_Cova.insert_header(info, 4);
		file_Cova.close();

		TreesData->print_matrix2("Covariance Matrix", (char *)outCovaName2);

        if(paras["-o"] == (String) "Community")
        {
            Compute_Community(TreesData, paras, "File-covariance");
        }
    }
}

void Compute_Community(Trees *TreesData, map<String, String> &paras, String memorydata)
{
    int modelType = 0;
    if(paras["-cm"] == (String) "CNM")
        modelType = 3;
    else 
    if(paras["-cm"] == (String) "CPM")
        modelType = 4;
    else
    if(paras["-cm"] == (String) "ERNM")
        modelType = 2;
    else 
    if(paras["-cm"] == (String) "NNM")
        modelType = 1;

    string stdparam3 = (char *) paras["-hf"];
    string stdparam4 = (char *) paras["-lf"];

    if(paras["-lm"] == (String) "auto")
    {
        if(TreesData->compute_community_automatically(memorydata, modelType, stdparam3, stdparam4))
        {
            cout << "Successfully detected communities of " << memorydata << " by model: " << paras["-cm"] << " with high freq. bound:" << stdparam3 << ", low freq. bound:" << stdparam4 << "!" << endl;
            cout << "Lambdas are chosen automatically." << endl;
        }
        return;
    }

    Array<double> param1;
    double lpiv = atof((char *) paras["-lpiv"]);
    double lp = atof((char *) paras["-lp"]);
    double lps = atof((char *) paras["-lps"]);
    double lpe = atof((char *) paras["-lpe"]);
        
    if(0 == lpiv)
    {
        param1.resize(1);
        param1[0] = lp;
    } else
    {
        int size = (int) ((lpe - lps) / lpiv + 1);
        param1.resize(size);
        for(int i = 0; i < size; i++)
            param1[i] = lps + i * lpiv;
    }
    
    Array<double> param2;
    double lniv = atof((char *) paras["-lniv"]);
    double ln = atof((char *) paras["-ln"]);
    double lns = atof((char *) paras["-lns"]);
    double lne = atof((char *) paras["-lne"]);
        
    if(0 == lniv)
    {
        param2.resize(1);
        param2[0] = ln;
    } else
    {
        int size = (int) ((lne - lns) / lniv + 1);
        param2.resize(size);
        for(int i = 0; i < size; i++)
            param2[i] = lns + i * lniv;
    }

    if(TreesData->compute_community_manually(memorydata, modelType, param1, param2, stdparam3, stdparam4))
    {
        cout << "Successfully detected communities of " << memorydata << " by model: " << paras["-cm"] << " with high freq. bound:" << stdparam3 << ", low freq. bound:" << stdparam4 << "!" << endl;
        cout << "Lambda positive: " << param1 << endl;
        cout << "Lambda negative: " << param2 << endl;
    }
}

void Compute_Distance(Trees *TreesData, map<String, String> &paras)
{
    String memorydata;
	String** info = new String*[10];
	for (int i = 0; i < 10; i++)
		info[i] = new String[2];
	info[0][0] = "created";			info[0][1] = time_stamp();
	info[1][0] = "output_type";		info[1][1] = "Distance matrix of trees";
	info[2][0] = "rooted";			info[2][1] = paras["-r"];
	info[3][0] = "weighted";		info[3][1] = paras["-w"];
	info[4][0] = "distance_type";	info[4][1] = paras["-dm"];
	info[5][0] = "source";			info[5][1] = paras["-f"];

    if(paras["-ft"] == (String) "Trees")
    {
        bool dis;
        if(paras["-dm"] == (String) "URF")
        {
            memorydata = (String) "Unweighted RF-distance";
            dis = TreesData->Compute_RF_dist_by_hash(false);
        }
        else
        if(paras["-dm"] == (String) "RF")
        {
            memorydata = (String) "Weighted RF-distance";
            dis = TreesData->Compute_RF_dist_by_hash(true);
        }
        else
        if(paras["-dm"] == (String) "Mat")
        {
            memorydata = (String) "Matching-distance";
            dis = TreesData->Compute_Matching_dist();
        } 
        else
        if(paras["-dm"] == (String) "SPR")
        {
            memorydata = (String) "SPR-distance";
            dis = TreesData->Compute_SPR_dist();
        } 
        else
        {
              cout << "Error: Setting of -dm is not correct. Unable to compute distance matrix" << endl;
              return;
        }
    
        if(dis)
        {
            std::cout << "Successfully computed " << memorydata << " distance." << std::endl;
        } else
        {
            std::cout << "Error: Unable to compute " << memorydata << " distance." << std::endl;
            return;
        }
        
        string outDistName = TreesData->make_DISToutput_name(memorydata);
        TreesData->print_matrix(memorydata, outDistName);
        cout << "Successfully printed " << memorydata << " matrix!" << endl;



		String outDistName2 = paras["-path"];
		outDistName2 += "Distance.out";


		File file_Dist(outDistName2);
		file_Dist.clean();
		file_Dist.insert_header(info, 6);
		file_Dist.close();
		TreesData->print_matrix2(memorydata, (char *)outDistName2);
    }
    
    if(paras["-ft"] == (String) "Dist")
    {
        String fname = paras["-f"];
        string stdfname = (char *) fname;
        TreesData->load_distfile(stdfname);
        memorydata = "File-distance";
    }
    
    if(paras["-o"] == (String) "Affinity" || (paras["-o"] == (String) "Community" && paras["-t"] == (String) "Affinity"))
        Compute_Affinity(TreesData, paras, memorydata);
        
}

void Compute_Affinity(Trees *TreesData, map<String, String> &paras, String memorydata)
{
    if(paras["-am"] == (String) "Exp")
    {
        std::cout << "Applying exponential to distance matrix to obtain affinity matrix" << endl;
        TreesData->Compute_Affinity_dist(memorydata, 2);
    } else
    if(paras["-am"] == (String) "Rec")
    {
        std::cout << "Applying reciprocal to distance matrix to obtain affinity matrix" << endl;
        TreesData->Compute_Affinity_dist(memorydata, 1);
    } else
    {
        std::cout << "Error: setting of -am is not correct. Unable to compute Affinity matrix" << std::endl;
        return;
    }
    cout << "Successfully computed affinity matrix" << endl;
    
    if(paras["-ft"] == (String) "Trees")
    {
    	if(paras["-am"] == (String) "Rec")
    	{
    		if(paras["-dm"] == (String) "URF")
        	{
            	memorydata = (String) "Affinity-Reciprocal-URF";
        	}
        	else
        	if(paras["-dm"] == (String) "RF")
        	{
            	memorydata = (String) "Affinity-Reciprocal-RF";
        	}
        	else
        	if(paras["-dm"] == (String) "Mat")
        	{
            	memorydata = (String) "Affinity-Reciprocal-match";
        	}
        	else
        	if(paras["-dm"] == (String) "SPR")
        	{
            	memorydata = (String) "Affinity-Reciprocal-SPR";
        	}
        }
        else
        if(paras["-am"] == (String) "Exp")
        {
        	if(paras["-dm"] == (String) "URF")
        	{
            	memorydata = (String) "Affinity-Exponential-URF";
        	}
        	else
        	if(paras["-dm"] == (String) "RF")
        	{
            	memorydata = (String) "Affinity-Exponential-RF";
        	}
        	else
        	if(paras["-dm"] == (String) "Mat")
        	{
            	memorydata = (String) "Affinity-Exponential-match";
        	}
        	else
        	if(paras["-dm"] == (String) "SPR")
        	{
            	memorydata = (String) "Affinity-Exponential-SPR";
        	}
        }
    } else
    if(paras["-ft"] == (String) "Dist")
    {
        if(paras["-am"] == (String) "Rec")
            memorydata = (String) "Affinity-Reciprocal-filedist";
        else
        if(paras["-am"] == (String) "Exp")
            memorydata = (String) "Affinity-Exponential-filedist";
    } else
    {
        cout << "Error: Incorrect affinity matrix." << endl;
        return;
    }
    
    string outAffName = TreesData->make_DISToutput_name(memorydata);
    TreesData->print_matrix(memorydata, outAffName);
    cout << "Successfully printed " << memorydata << " matrix!" << endl;
    
    if(paras["-o"] == (String) "Community")
        Compute_Community(TreesData, paras, memorydata);;
}

void Compute_Consensus_Tree(Trees *TreesData, map<String, String> &paras)
{
	String memorydata;
	String** info = new String*[10];
	for (int i = 0; i < 10; i++)
		info[i] = new String[2];
	info[0][0] = "created";			info[0][1] = time_stamp();
	info[1][0] = "output_type";		info[1][1] = "Consensus trees";
	info[2][0] = "output_format";	info[2][1] = paras["-ctm"];
	info[3][0] = "tree_type";		info[3][1] = paras["-ct"];
	info[4][0] = "tree_index";		info[4][1] = paras["-if"];
	info[5][0] = "source";			info[5][1] = paras["-f"];
	Array<int> *treeidx = TreesData->getidxlist();
	if (paras["-if"] == (String) "")
	{
		treeidx->resize(TreesData->Get_n_trees());
		for (int i = 0; i < TreesData->Get_n_trees(); i++)
			(*treeidx)[i] = i;
		std::cout << "Consider all trees to compute the consensus tree." << std::endl;
	}
	else
	{
		String idxfname = paras["-if"];
		File file(idxfname);
		if (!file.is_open())
		{
			cout << "Error: Can not open the data file!" << endl;
			return;
		}
		file.seek(0);
		int num = file.lines();
		file.seek(0);
		treeidx->resize(num);

		for (int i = 0; i < num; i++)
		{
			file >> (*treeidx)[i];
		}
		std::cout << "Consider indices from file: " << idxfname << "; " << num << " trees." << std::endl;
	}

	if (paras["-ct"] == (String) "Majority")
	{
		if (TreesData->compute_consensus_tree(MAJORITYTREE, ""))
			std::cout << "Successfully computed the majority consensus tree!" << std::endl;
	}
	else
		if (paras["-ct"] == (String) "Strict")
		{
			if (TreesData->compute_consensus_tree(STRICTTREE, ""))
				std::cout << "Successfully computed the strict consensus tree!" << std::endl;
		}
		else
		{
			cout << "Error: Setting of -ct is not correct. Unable to compute consensus tree." << endl;
			return;
		}

	String confname = paras["-f"].before('.');
	confname += "_";
	if (paras["-if"] == (String) "")
	{
		confname += "All_";
	}
	else
	{
		confname += paras["-if"].before('.');
		confname += "_";
	}
	confname += paras["-ct"];
	confname += "_consensus_tree.out";
	string outName = (char *)confname;
	//string outName2 = (char *)paras["-f"];

	String outName2 = paras["-path"];
	outName2 += "Consensus.out";


	File file_Conf(outName2);
	

	if (paras["-cfm"] == (String) "Newick")
	{
		outConf.open(outName, ios::trunc);
		outConf.close();
		TreesData->WriteConsensusTree(outName, NEWICK);
		cout << "Successfully outputted Newick format trees to file: " << confname << endl;


		file_Conf.clean();
		file_Conf.insert_header(info, 6);
		file_Conf.close();
		TreesData->WriteConsensusTree(outName2, NEWICK);
	}
	else
		if (paras["-cfm"] == (String) "Nexus")
		{
			outConf.open(outName, ios::trunc);
			outConf.close();
			TreesData->WriteConsensusTree(outName, NEXUS);
			cout << "Successfully outputted Newick format trees to file: " << confname << endl;


			file_Conf.clean();
			file_Conf.insert_header(info, 6);
			file_Conf.close();
			TreesData->WriteConsensusTree(outName2, NEXUS);
		}
		else
		{
			cout << "warning: setting of -cfm is not correct. output consensus tree to Nexus format." << endl;
			outConf.open(outName, ios::trunc);
			outConf.close();
			TreesData->WriteConsensusTree(outName, NEXUS);
			cout << "Successfully outputted Newick format trees to file: " << confname << endl;

			info[2][1] = "Nexus";
			file_Conf.clean();
			file_Conf.insert_header(info, 6);
			file_Conf.close();
			TreesData->WriteConsensusTree(outName2, NEXUS);
			return;
		}
}

void dimest_driver(String fname, String Est, String Init, String para_fname)
{
    cout << "loading" << endl;
    DimEst dimest(fname, Est, Init, para_fname);
    cout << "computing" << endl;
    dimest.Compute_Dim();
    cout << "outputing" << endl;
    dimest.output_to_files();
};

void nldr_driver(map<String, String> &paras)
{
	//NLDR nLdr(fname, ftype, dim, cost, algo, init_md, flag, seed, para_fname);
	NLDR nLdr(paras);
	nLdr.Compute_NLDR();
	nLdr.result_analysis();
	nLdr.output_to_files();
}



int check_info(String** info, String content, int &lines) {
	int ind = lines;
	for (int i = 0; i < lines; i++) {
		if (info[i][0] == content) {
			ind = i;
			break;
		}
	}
	if (ind == lines)
		lines++;
	return ind;
}

void aff_driver(map<String, String> & paras) {
	Matrix<double> dist_mat;
	Matrix<double> aff_mat;
	int size;
	double eps = 1000000;
	double ratio = 0.1;

	String fname = paras["-f"];
	File D_file(fname);
	if (!D_file.is_open())
	{
		std::cout << "Error: File \"" << fname << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
		exit(0);
	}
	int pos = D_file.end_header();
	D_file.seek(pos);
	size = D_file.lines();
	D_file.seek(pos);

	//D.resize(size, size);
	//for (int i = 0; i < size; i++)
	//	for (int j = 0; j <= i; j++)
	//	{
	//		D_file >> D.matrix[i][j];
	//		D.matrix[j][i] = D.matrix[i][j];
	//	}
	String tree;
	double index;
	size--;
	dist_mat.resize(size, size);
	D_file >> tree;
	for (int i = 0; i < size; i++)
		D_file >> index;
	for (int i = 0; i < size; i++)
	{
		D_file >> index;
		for (int j = 0; j <= i; j++)
		{
			D_file >> dist_mat.matrix[i][j];
			dist_mat.matrix[j][i] = dist_mat.matrix[i][j];
		}
	}

	String** info = new String*[10];
	for (int i = 0; i < 10; i++)
		info[i] = new String[2];
	int lines = D_file.load_header(info);
	info[0][0] = "created";			info[0][1] = time_stamp();
	info[1][0] = "output_type";		info[1][1] = "Affinity matrix";
	info[2][0] = "source";			info[2][1] = paras["-f"];
	info[3][0] = "method";			info[3][1] = paras["-am"];


	aff_mat.resize(size, size);


	if (paras["-am"] == (String) "Exp")
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				if (dist_mat(i, j) > 0 && eps > dist_mat(i, j))
					eps = dist_mat(i, j);
		eps = eps * ratio;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				aff_mat.matrix[i][j] = 1.0 / (eps + dist_mat(i, j));
		}
	}
	else if ((paras["-am"] == (String) "Rec")) {
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				aff_mat.matrix[i][j] = exp(-dist_mat(i, j));
		}
	}

	String Aff_Filename = paras["-path"];
	Aff_Filename += "Affinity_Matrix.out";
	File file_Aff(Aff_Filename);
	file_Aff.clean();
	file_Aff.insert_header(info, 4);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j <= i; j++)
			file_Aff << aff_mat(i, j) << "\t";
		file_Aff << "" << std::endl;
	}
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