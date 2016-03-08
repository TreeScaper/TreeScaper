
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
        cout << "loading all parameters" << endl;
        driver(default_paras[0], default_paras[1], default_paras[2], default_paras[3], default_paras[4], default_paras[5], default_paras[6], atoi(default_paras[7]), default_paras[8]);
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
        cout << "loading all parameters" << endl;
        
        trees_driver(paras);
    } else
    {
        cout << "error: input error!" << endl;
    }
    
    return 0;
}

void trees_driver(map<String, String> &paras)
{
    String fname = paras["-f"];
    string stdfname = (char *) fname;
    File file(fname);
    if(! file.is_open())
    {
        cout << "Error: can not open the data file!" << endl;
        return;
    }
    
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
    
        cout << "successfully computed bipartitation matrix." << endl;
        
        String fname = paras["-f"];
        string stdfname = (char *) fname;
        string namebipartmatrix;
        ofstream outBipartMatrix;
        if(paras["-bfm"] == (String) "list")
        {
            string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, (String) "List format");
            outBipartMatrix.open(namebipartmatrix.c_str());
            cout << "Output bipartition matrix in list format to " << namebipartmatrix << endl;
            TreesData->OutputBipartitionMatrix(outBipartMatrix, RCVLIST);
        } else
        if(paras["-bfm"] == (String) "matrix")
        {
            string namebipartmatrix = TreesData->make_Bipart_Matrix_name(stdfname, (String) "Matrix format");
            outBipartMatrix.open(namebipartmatrix.c_str());
            cout << "Output bipartition matrix in matrix format to " << namebipartmatrix << endl;
            TreesData->OutputBipartitionMatrix(outBipartMatrix, FULLMATRIX);
        } else
        {
            cout << "error: setting of -bfm is not correct. Unable to output bipartition matrix" << endl;
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
    if(paras["-ft"] == (String) "Trees")
    {
        TreesData->Compute_Bipart_Covariance();
        cout << "Successfully computed covariance matrix of bipartition." << endl;
        TreesData->print_matrix("Covariance Matrix");
        cout << "Successfully printed Covariance Matrix matrix!" << endl;
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
              cout << "error: setting of -dm is not correct. Unable to compute distance matrix" << endl;
              return;
        }
    
        if(dis)
        {
            std::cout << "successfully computed " << memorydata << " distance." << std::endl;
        } else
        {
            std::cout << "error: Unable to compute " << memorydata << " distance." << std::endl;
            return;
        }
        
        
        TreesData->print_matrix(memorydata);
        cout << "Successfully printed " << memorydata << " matrix!" << endl;
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
        std::cout << "error: setting of -am is not correct. Unable to compute Affinity matrix" << std::endl;
        return;
    }
    cout << "Succeed in computing affinity matrix" << endl;
    if(paras["-ft"] == (String) "Trees")
    {
        if(paras["-dm"] == (String) "URF")
        {
            memorydata = (String) "Affinity-URF";
        }
        else
        if(paras["-dm"] == (String) "RF")
        {
            memorydata = (String) "Affinity-RF";
        }
        else
        if(paras["-dm"] == (String) "Mat")
        {
            memorydata = (String) "Affinity-match";
        }
        else
        if(paras["-dm"] == (String) "SPR")
        {
            memorydata = (String) "Affinity-SPR";
        } 
    } else
    if(paras["-ft"] == (String) "Dist")
    {
        memorydata = (String) "Affinity-filedist";
    } else
    {
        cout << "Error: incorrect affinity matrix." << endl;
        return;
    }
    
    TreesData->print_matrix(memorydata);
    cout << "Successfully printed " << memorydata << " matrix!" << endl;
    
    if(paras["-o"] == (String) "Community")
        Compute_Community(TreesData, paras, memorydata);;
}

void Compute_Consensus_Tree(Trees *TreesData, map<String, String> &paras)
{
    Array<int> *treeidx = TreesData->getidxlist();
    if(paras["-if"] == (String) "")
    {
        treeidx->resize(TreesData->Get_n_trees());
        for(int i = 0; i < TreesData->Get_n_trees(); i++)
            (*treeidx)[i] = i;
        std::cout << "Consider all trees to compute the consensus tree." << std::endl;
    } else
    {
        String idxfname = paras["-if"];
        File file(idxfname);
        if(! file.is_open())
        {
            cout << "Error: can not open the data file!" << endl;
            return;
        }
        file.seek(0);
        int num = file.lines();
        file.seek(0);
        treeidx->resize(num);

        for(int i = 0; i < num; i++)
        {
            file >> (*treeidx)[i];
        }
        std::cout << "Consider indices from file: " << idxfname << "; " << num << " trees." << std::endl;
    }
    
    if(paras["-ct"] == (String) "Majority")
    {
        if(TreesData->compute_consensus_tree(MAJORITYTREE, ""))
            std::cout << "successfully computed the majority consensus tree!" << std::endl;
    } else
    if(paras["-ct"] == (String) "Strict")
    {
        if(TreesData->compute_consensus_tree(STRICTTREE, ""))
            std::cout << "successfully computed the strict consensus tree!" << std::endl;
    } else
    {
          cout << "error: setting of -ct is not correct. Unable to compute consensus tree." << endl;
          return;
    }
    
    String confname = paras["-f"].before('.');
    confname += "_";
    if(paras["-if"] == (String) "")
    {
        confname += "All_";
    } else
    {
        confname += paras["-if"].before('.');
        confname += "_";
    }
    confname += paras["-ct"];
    confname += "_consensus_tree.out";
    string outName = (char *) confname;

    if(paras["-cfm"] == (String) "Newick")
    {
        TreesData->WriteConsensusTree(outName, NEWICK);
        cout << "successfully outputted Newick format trees to file: " << confname << endl;
    }
    else
    if(paras["-cfm"] == (String) "Nexus")
    {
        TreesData->WriteConsensusTree(outName, NEXUS);
        cout << "successfully outputted Nexus format trees to file: " << confname << endl;
    } else
    {
          cout << "warning: setting of -cfm is not correct. output consensus tree to Nexus format." << endl;
        TreesData->WriteConsensusTree(outName, NEXUS);
        cout << "successfully outputted Nexus format trees to file: " << confname << endl;
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

void driver(String fname, String ftype, String dim, String cost, String algo, String init_md, String flag, long seed, String para_fname)
{
    NLDR nLdr(fname, ftype, dim, cost, algo, init_md, flag, seed, para_fname);
    nLdr.Compute_NLDR();
    nLdr.result_analysis();
    nLdr.output_to_files();
}
