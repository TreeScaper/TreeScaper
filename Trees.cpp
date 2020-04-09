
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

// Trees.cpp
// Member function definitions for class Trees.h
//              by whuang gzhou


#ifndef TREES_CPP
#define TREES_CPP

#include "Trees.h"



//########################ZD comment########################################
//# This structure stores a set of trees and
//# associated data structure for 4 distances,
//# which is implemented by a void pointer.
//# Unless the particular routine for computing
//# distance is called, the cooresponding pointer
//# will be empty by default.
//########################ZD comment########################################

std::map<String, void ***> StrToDist;
//std::map<String, int ***> StrToIntDist;
//std::map<String, double ***> StrToDoubleDist;

Trees::Trees()
{
    String str_matrix("Unweighted RF-distance");
    StrToDist["Unweighted RF-distance"] = (void ***) &dist_URF;
    str_matrix = (String)("Matching-distance");
    StrToDist["Matching-distance"] = (void ***) &dist_match;
    str_matrix = (String)("SPR-distance");
    StrToDist["SPR-distance"] = (void ***) &dist_SPR;

    str_matrix = (String)("Weighted RF-distance");
    StrToDist["Weighted RF-distance"] = (void ***) &dist_RF;
    str_matrix = (String) "Geodesic-distance";
    StrToDist["Geodesic-distance"] = (void ***) &dist_geo;
    str_matrix = (String) "File-distance";
    StrToDist["File-distance"] = (void ***) &dist_file;
    str_matrix = (String) "File-coordinate";
    StrToDist["File-coordinate"] = (void ***) &coord_file;

    str_matrix = (String)("Covariance Matrix");
    StrToDist["Covariance Matrix"] = (void ***) &treecov;
    str_matrix = (String)("File-covariance");
    StrToDist["File-covariance"] = (void ***) &filecov;

    str_matrix = (String) "Affinity-Reciprocal-URF";
    StrToDist["Affinity-Reciprocal-URF"] = (void ***) &affi_Recip_URF;
    str_matrix = (String) "Affinity-Reciprocal-RF";
    StrToDist["Affinity-Reciprocal-RF"] = (void ***) &affi_Recip_RF;
    str_matrix = (String) "Affinity-Reciprocal-match";
    StrToDist["Affinity-Reciprocal-match"] = (void ***) &affi_Recip_match;
    str_matrix = (String) "Affinity-Reciprocal-SPR";
    StrToDist["Affinity-Reciprocal-SPR"] = (void ***) &affi_Recip_SPR;
    str_matrix = (String) "Affinity-Reciprocal-geodesic";
    StrToDist["Affinity-Reciprocal-geodesic"] = (void ***) &affi_Recip_geo;
    str_matrix = (String) "Affinity-Reciprocal-filedist";
    StrToDist["Affinity-Reciprocal-filedist"] = (void ***) &affi_Recip_file;

    str_matrix = (String) "Affinity-Exponential-URF";
    StrToDist["Affinity-Exponential-URF"] = (void ***) &affi_Exp_URF;
    str_matrix = (String) "Affinity-Exponential-RF";
    StrToDist["Affinity-Exponential-RF"] = (void ***) &affi_Exp_RF;
    str_matrix = (String) "Affinity-Exponential-match";
    StrToDist["Affinity-Exponential-match"] = (void ***) &affi_Exp_match;
    str_matrix = (String) "Affinity-Exponential-SPR";
    StrToDist["Affinity-Exponential-SPR"] = (void ***) &affi_Exp_SPR;
    str_matrix = (String) "Affinity-Exponential-geodesic";
    StrToDist["Affinity-Exponential-geodesic"] = (void ***) &affi_Exp_geo;
    str_matrix = (String) "Affinity-Exponential-filedist";
    StrToDist["Affinity-Exponential-filedist"] = (void ***) &affi_Exp_file;

    str_matrix = (String)("File-affinity");
    StrToDist["File-affinity"] = (void ***) &fileaffinity;

    isrooted = true;
    isweighted = false;
    treesfilename = "";
    n_trees = 0;
    treeset = NULL;

    // sparse representation of bipartition matrix
    sbipartmatrix = NULL;
    numberofbipartition = NULL;
    bipart_count = NULL;

//    bipartcovariance = NULL;
    treecov_size = 0;
    treecov = NULL;
    filecov_size = 0;
    filecov = NULL;
    dist_URF = NULL;
    dist_RF = NULL;
    dist_match = NULL;
    dist_SPR = NULL;
    dist_geo = NULL;

    affi_Recip_URF = NULL;
    affi_Recip_RF = NULL;
    affi_Recip_match = NULL;
    affi_Recip_SPR = NULL;
    affi_Recip_geo = NULL;

    affi_Exp_URF = NULL;
    affi_Exp_RF = NULL;
    affi_Exp_match = NULL;
    affi_Exp_SPR = NULL;
    affi_Exp_geo = NULL;

    affinityfile_size = 0;
    fileaffinity = NULL;

    file_distsize = 0;
    dist_file = NULL;
    affi_Recip_file = NULL;
    affi_Exp_file = NULL;
    file_coordinatesize = 0;
    file_coordinatedim = 0;
    coord_file = NULL;

    covariance_freeid_size = 0;
    covariance_freeid = NULL;
    covariance_nonfree_id_size = 0;
    covariance_nonfree_id = NULL;
    com_info = NULL;
    com_info_col = 0;

//    bipartFreq = NULL;
//    bipartFreqIdx = NULL;
}

void Trees::destructor()
{
    deletetrees();
    deleteConsensustree();
    delete [] numberofbipartition;
    numberofbipartition = NULL;
    delete [] bipart_count;
    bipart_count = NULL;
    vec_hashrf.hashrfmap_clear();

    delete_double_array(treecov, treecov_size);
    delete_double_array(filecov, filecov_size);
    delete_double_array(dist_URF, n_trees);
    delete_double_array(dist_RF, n_trees);
    delete_double_array(dist_match, n_trees);
    delete_double_array(dist_SPR, n_trees);
    delete_double_array(dist_geo, n_trees);

    delete_double_array(affi_Recip_URF, n_trees);
    delete_double_array(affi_Recip_RF, n_trees);
    delete_double_array(affi_Recip_match, n_trees);
    delete_double_array(affi_Recip_SPR, n_trees);
    delete_double_array(affi_Recip_geo, n_trees);

    delete_double_array(affi_Exp_URF, n_trees);
    delete_double_array(affi_Exp_RF, n_trees);
    delete_double_array(affi_Exp_match, n_trees);
    delete_double_array(affi_Exp_SPR, n_trees);
    delete_double_array(affi_Exp_geo, n_trees);

    delete_double_array(coord_file, file_coordinatesize);
    delete_double_array(dist_file, file_distsize);
    delete_double_array(affi_Recip_file, file_distsize);
    delete_double_array(affi_Exp_file, file_distsize);
    delete_double_array(fileaffinity, affinityfile_size);
    delete_double_array(com_info, covariance_nonfree_id_size + 4);

    delete sbipartmatrix;
    sbipartmatrix = NULL;
//    delete_double_array(bipartFreq, splits);
//    delete [] bipartFreqIdx;
//    bipartFreqIdx = NULL;

    for(std::map<unsigned long long, Array<char> *>::iterator itr = hash2bitstr.begin(); itr != hash2bitstr.end(); itr++)
    {
        delete itr->second;
        itr->second = NULL;
    }
    hash2bitstr.clear(); //


    // set to be default
    isrooted = true;
    isweighted = false;
    treesfilename = "";
    n_trees = 0;
    treecov_size = 0;
    filecov_size = 0;
    file_distsize = 0;
    file_coordinatesize = 0;
    file_coordinatedim = 0;
    affinityfile_size = 0;

    delete [] covariance_freeid;
    covariance_freeid = NULL;
    covariance_freeid_size = 0;

    delete [] covariance_nonfree_id;
    covariance_nonfree_id = NULL;
    covariance_nonfree_id_size = 0;
}

void Trees::deletetrees()
{
    //cout << "set address:" << (long) treeset << endl;//----
    if(treeset != NULL)
    {
        for(unsigned int i = 0; i < n_trees; i++)
        {
            //cout << "i-th tree:" << i << ", address:" << (long) treeset[i] << endl;//----
            if(treeset[i] != NULL)
            {
                TreeOPE::killnewicktree(treeset[i]);
                treeset[i] = NULL;
            }
        }
        delete [] treeset;
        treeset = NULL;
    }
//    if(consensustrees.is_empty() == 0)
//    {
//        for(unsigned int i = 0; i < consensustrees.get_length(); i++)
//        {
//            TreeOPE::killnewicktree(consensustrees[i]);
//            consensustrees[i] = NULL;
//        }
//        consensustrees.resize(0);
//    }
}

void Trees::deleteConsensustree()
{
    if(consensustrees.is_empty() == 0)
    {
        for(unsigned int i = 0; i < consensustrees.get_length(); i++)
        {
            TreeOPE::killnewicktree(consensustrees[i]);
            consensustrees[i] = NULL;
        }
        consensustrees.resize(0);
    }
}

void Trees::deleteBipartitionMatrix()
{
    if(treeset != NULL)
    {
        for(unsigned int i = 0; i < n_trees; i++)
        {
            if(treeset[i] != NULL)
            {
                TreeOPE::killnewicktreeBitstr(treeset[i]);
            }
        }
    }

    delete [] bipart_count;
    bipart_count = NULL;
    vec_hashrf.hashrfmap_clear();

    delete sbipartmatrix;
    sbipartmatrix = NULL;

    for(std::map<unsigned long long, Array<char> *>::iterator itr = hash2bitstr.begin(); itr != hash2bitstr.end(); itr++)
    {
        delete itr->second;
        itr->second = NULL;
    }
    hash2bitstr.clear(); //
}

Trees::~Trees()
{
    destructor();
}

template<class T>
void Trees::delete_double_array(T ** (&arr), int n)
{
    if(arr != NULL)
    {
        for(int i = 0; i < n; i++)
        {
            delete [] arr[i];
        }
        delete [] arr;
        arr = NULL;
    }
}

template<class T>
void Trees::delete_double_array(T *** arr, int n)
{
    if((*arr) != NULL)
    {
        for(int i = 0; i < n; i++)
        {
            delete [] (*arr)[i];
        }
        delete [] (*arr);
        (*arr) = NULL;
    }
}

void Trees::delete_matrix(String str_matrix)
{
    cout << "Deleting: " << str_matrix << "\n\n";
    if(str_matrix == (String)"Bipartition Matrix")
    {
        deleteBipartitionMatrix();
    } else if(str_matrix == (String)"Covariance Matrix")
        delete_double_array((double ***) StrToDist[str_matrix], treecov_size);
    else if(str_matrix == (String)"File-covariance")
        delete_double_array((double ***) StrToDist[str_matrix], filecov_size);
    else if(str_matrix == (String)"File-distance" || str_matrix == (String)"Affinity-filedist")
        delete_double_array((double ***) StrToDist[str_matrix], file_distsize);
    else if(str_matrix == (String)"File-coordinate")
        delete_double_array((double ***) StrToDist[str_matrix], file_coordinatesize);
    else if(str_matrix == (String)"File-affinity")
        delete_double_array((double ***) StrToDist[str_matrix], affinityfile_size);
    else if(str_matrix == (String)"Matching-distance" || str_matrix == (String)"SPR-distance")
        delete_double_array((int ***) StrToDist[str_matrix], n_trees);
    else
        delete_double_array((double ***) StrToDist[str_matrix], n_trees);
    cout << "Done!!\n\n";
}

string Trees::make_DISToutput_name(String str_matrix)
{
    String result = treesfilename.c_str();
    result = result.before('.');
    result += "_";

    if(treesfilename != "")
    {
        if(isrooted)
            result += "rooted_";
        else
            result += "unrooted_";

        if(isweighted)
            result += "weighted_";
        else
            result += "unweighted_";
    }

    if(str_matrix == (String) "Unweighted RF-distance")
        result += "RF-distance";
    else if(str_matrix == (String) "Weighted RF-distance")
        result += "RF-distance";
    else if(str_matrix == (String) "Covariance Matrix")
        result += "Covariance_Matrix";
    else
    {
        result += (char *) str_matrix;
    }

    result += ".out";
    return (char *) result;
}

template<class T>
void Trees::print_double_array(T *** arr, int n, string outfile)
{
    ofstream fout;
    if (outfile != "")
        fout.open(outfile.c_str());

    if((*arr) != NULL)
    {
        fout << setw(5) << "tree" << ' ';
        for (int i = 0, c = 1; i < n; i++, c++)
            fout << setw(5) << c << "\t" << ' ';
        fout << endl;

        for (int i = 0, c = 1; i < n; i++, c++)
        {
            fout << setw(5) << c << "\t" << ' ';
            for (int j = 0; j <= i; j++)
            {
                fout << setw(5) << (*arr)[i][j] << "\t" << ' ';
            }
            fout << endl;
        }
    }
    if (outfile != "")
        fout.close();
}

template<class T>
void Trees::print_coordinate_matrix(T *** arr, int n, int m, string outfile)
{
    ofstream fout;
    if (outfile != "")
        fout.open(outfile.c_str());

    if((*arr) != NULL)
    {
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < m; j++)
            {
                fout << setw(5) << (*arr)[i][j] << "\t";
            }
            fout << endl;
        }
    }
    if (outfile != "")
        fout.close();
}

void Trees::print_matrix(String str_matrix, string outfile)
{
    if (str_matrix == (String)"Covariance Matrix")
        print_double_array((double ***) StrToDist[str_matrix], treecov_size, outfile);
    else if(str_matrix == (String)"File-covariance")
        print_double_array((double ***) StrToDist[str_matrix], filecov_size, outfile);
    else if(str_matrix == (String)"Matching-distance" || str_matrix == (String)"SPR-distance")
        print_double_array((int ***) StrToDist[str_matrix], n_trees, outfile);
    else if(str_matrix == (String)"File-distance" || str_matrix == (String)"Affinity-filedist")
        print_double_array((double ***) StrToDist[str_matrix], file_distsize, outfile);
    else if(str_matrix == (String)"File-coordinate")
        print_coordinate_matrix((double ***) StrToDist[str_matrix], file_coordinatesize, file_coordinatedim, outfile);
    else if(str_matrix == (String)"File-affinity")
        print_double_array((double ***) StrToDist[str_matrix], affinityfile_size, outfile);
    else
    {
        print_double_array((double ***) StrToDist[str_matrix], n_trees, outfile);
    }
}

//########################ZD comment########################################
//# Initialize a set of trees by calling constructer
//# (parsetree) from TreeOPE for Newick tree, creating
//# leaveslabelsmaps for NEXUS tree. Complicated string
//#operations are involved here, which is unnecessary.
//########################ZD comment########################################

void Trees::initialTrees(string fname)
{
    treesfilename = fname;
    treeformat = CheckTreeFormat();
    n_trees = 0;
    if(treeformat == NEXUS)
    {
        leaveslabelsmaps.clear();
        ifstream treefile(treesfilename.c_str(), ios::binary);
        string line;
        int num_line,Line_taxa_begin,count;
        int a;
        string b;
        char namebuff[256];

        while (!treefile.eof())
        {
            for (num_line = 0; getline(treefile, line); num_line++)
            {
                if (line.find_first_not_of("\t\n ") == string::npos) continue;
                if (line.find("Translate", 0) != string::npos || line.find("TRANSLATE", 0) != string::npos || line.find("translate", 0) != string::npos)
                {
                    Line_taxa_begin = num_line;
                    for (count = Line_taxa_begin + 1; getline(treefile, line); count++)
                    {
                        if (line.find_first_not_of("\t\n ") == string::npos) continue;
                        if ((line.find(";", 0) == string::npos) || (line.find("(") == string::npos
                                                                    && find_if(line.begin(), line.end(), (int(*)(int))std::isdigit) != line.end()))
                        {
                            istringstream iss(line);
                            iss >> a >> b;
                            b.erase(remove(b.begin(),b.end(),','),b.end());
                            b.erase(remove(b.begin(),b.end(),';'),b.end());

                            if(b.length() > 255)
                            {
                                cout << "The longest name only can include 255 chars!\n\n";
                                exit(0);
                            }
                            strcpy(namebuff, b.c_str());
                            namebuff[b.length()] = NULL;

                            leaveslabelsmaps.push(string (namebuff));
                        }
                        else if(line.find("(",0) != string::npos)
                        {
                            n_trees++;
                        }
                    }
                }
            }
        }
        treefile.close();
    } else
    if(treeformat == NEWICK)
    {
//        ofstream outNexus;
        FILE *fp = NULL;
        string line, s;
        fp = fopen(treesfilename.c_str(), "r");
//        unsigned int Num_tree = 0;
        int err, count;
        NEWICKTREE *newickTree = NULL;

        ifstream treeFile(treesfilename.c_str(), ios::out);
        if (treeFile)
        {
            while (getline(treeFile, line) && (line.find_first_not_of("\t\n ") != string::npos))
            {
                n_trees ++;
            }
        }

        newickTree = TreeOPE::loadnewicktree2(fp, &err);
        if (!newickTree)
        {
            switch (err)
            {
                case -1:
                    printf("Out of memory \n");
                    break;
                case -2:
                    printf("Pase error \n");
                    break;
                case -3:
                    printf("Cannot load file \n");
                    break;
                default:
                    printf("Error %d\n", err);
            }
        }
        try {
            leaveslabelsmaps.clear();
            TreeOPE::GetTaxaLabels(newickTree->root, leaveslabelsmaps);
        }
        catch (LabelMap::AlreadyPushedEx ex) {
            cerr << "Error : The label ' " << ex.label << "' appeared twice!\n\n";
            exit(2);
        }
        TreeOPE::killnewicktree(newickTree);
        fclose(fp);
        treeFile.close();
    } else
    {
        cout << "Warning: Unregonized file format!\n\n";
    }
    selected_trees.resize(n_trees + 1);
    selected_trees[0] = -1;
}

void Trees::Settreeroottype(bool isrt)
{
    isrooted = isrt;
}

void Trees::Settreeweighttype(bool iswt)
{
    isweighted = iswt;
}

Treefileformat Trees::CheckTreeFormat()
{
    int num_line;
    string line,first_line;
    ifstream treefile(treesfilename.c_str(), ios::binary);
    Treefileformat format;

    if (!treefile)
    {
        cout << "Unable to open the file!\n\n";
        exit(0);
    }

    while (!treefile.eof())
    {
        for (num_line = 0; getline(treefile, line); num_line++)
            if (num_line == 0)
            {
                istringstream iss(line);
                iss >> first_line;
                if (first_line.compare("#NEXUS") == 0)
                    format = NEXUS;
                else
                    format = NEWICK;
            }
    }
    treefile.close();
    return format;
};

//########################ZD comment########################################
//# Another constructor that seems to be duplicated
//# version of initialTrees. This constructer manipulate
//# the string in a complicated way to add weighted 1
//# as default for trees that is unweighted. Also note
//# that for both NEXUS and Newic trees, they are now
//# stored in parsetree.
//########################ZD comment########################################

void Trees::ReadTrees() // newick nexus
{
    deletetrees();
    treeset = new NEWICKTREE *[n_trees];
    int error;

    // fixed a leaf to be one side of the root
    std::string leafroot = "1";//--leaveslabelsmaps.name(0); // normalized unrooted tree by WH
    NEWICKNODE *lrpt = NULL; // normalized unrooted tree by WH
    int indexchild = -1; // normalized unrooted tree by WH

    if(treeformat == NEXUS)
    {
        ifstream treefile(treesfilename.c_str(), ios::binary);
        string line;
        int num_line, Line_taxa_begin, count;
        int treesetidx = 0;
        while (!treefile.eof())
        {
            for (num_line = 0; getline(treefile, line); num_line++)
            {
                if (line.find_first_not_of("\t\n ") == string::npos) continue;
                if (line.find("Translate", 0) != string::npos || line.find("TRANSLATE", 0) != string::npos || line.find("translate", 0) != string::npos)
                {
                    Line_taxa_begin = num_line;
                    for (count = Line_taxa_begin + 1; getline(treefile, line); count++)
                    {
                        if (line.find_first_not_of("\t\n ") == string::npos) continue;
                        if(line.find("(",0) != string::npos)
                        {
                            size_t tstart = line.find_first_of("(", 0);
                            if(line.find(":", 0) != string::npos)
                            {
                                char *btree = new char[line.substr(tstart).length() + 1];
                                strcpy(btree, line.substr(tstart).c_str());
                                btree[line.substr(tstart).length()] = NULL;
                                treeset[treesetidx] = TreeOPE::parsetree(btree, &error, NULL);
                                delete btree;
                            }
                            else
                            {
                                //----------- make it weighted string ----------------
                                char *btree = new char[2 * line.substr(tstart).length() + 1];
                                int idx = 0;
//                                const char *stree = line.substr(tstart).c_str();      // For Some reason this did not work for some file types in CLV (MAC)

                                char *stree = new char[line.substr(tstart).length() + 1];
                                strcpy(stree, line.substr(tstart).c_str());
                                stree[line.substr(tstart).length()] = '\0';



                                for(int i = 0; i < line.substr(tstart).length(); i++)
                                {
                                    if((stree[i] == ')' || stree[i] == ',') && i < line.substr(tstart).length() - 1)
                                    {
                                        btree[idx++] = ':';
                                        btree[idx++] = '1';
                                    }
                                    btree[idx++] = stree[i];
                                }
                                btree[idx] = '\0';
                                //-------------work but not good way-------------
                                treeset[treesetidx] = TreeOPE::parsetree(btree, &error, NULL);
                                delete btree;
                            }

                            if(!isrooted)
                            {
                                lrpt = TreeOPE::findleaf(leafroot, treeset[treesetidx]->root, NULL, &indexchild); // normalized unrooted tree by WH
                                TreeOPE::normailzedTree(lrpt, treeset[treesetidx], indexchild); // normalized unrooted tree by WH
                            }
                            treesetidx++;
                        }
                    }
                }
            }
        }
        treefile.close();
    } else
    if(treeformat == NEWICK)
    {
        string s;
//        fp = fopen(fileNames.c_str(), "r");
        int count = 0;

        ifstream treeFile(treesfilename.c_str(), ios::out);

        int NN = 10000;
        char *gzbuff = (char*)malloc(NN);

        treeFile.clear();
        treeFile.seekg(0);
        while (getline(treeFile, s) && (s.find_first_not_of("\t\n ") != string::npos))
        {
//            std::cout << "count:" << count << std::endl;//-- WHtest
            gzbuff[0] = NULL;
            if(s.find(":", 0) != string::npos)
            {
                char *btree = new char[s.length() + 1];
                strcpy(btree, s.c_str());
                btree[s.length()] = NULL;
                treeset[count] = TreeOPE::parsetree(btree, &error, NULL);
                delete btree;
            }
            else
            {
                //----------- make it weighted string ----------------
                char *btree = new char[2 * s.length() + 1];
                int idx = 0;
                const char *stree = s.c_str();
                for(int i = 0; i < s.length(); i++)
                {
                    if((stree[i] == ')' || stree[i] == ',') && i < s.length() - 1)
                    {
                        btree[idx++] = ':';
                        btree[idx++] = '1';
                    }
                    btree[idx++] = stree[i];
                }
                btree[idx] = NULL;
                //-------------work but not good way-------------
                treeset[count] = TreeOPE::parsetree(btree, &error, NULL);
                delete btree;
            }
//            cout << count << ":!!" << btree << "!!" << endl;//---

            TreeOPE::Label_strint((treeset[count])->root, leaveslabelsmaps);
            if(!isrooted)
            {
                lrpt = TreeOPE::findleaf(leafroot, treeset[count]->root, NULL, &indexchild); // normalized unrooted tree by WH
                TreeOPE::normailzedTree(lrpt, treeset[count], indexchild); // normalized unrooted tree by WH
            }
            count++;
        }

        free(gzbuff);
        treeFile.close();
    }
}

//########################ZD comment########################################
//# This is associated to the sum of the degrees, s,
//# of all nodes and the number of leaves, n, :
//# s/2-n
//########################ZD comment########################################

void Trees::compute_numofbipart()
{
    if(numberofbipartition != NULL)
        delete [] numberofbipartition;
    numberofbipartition = new int [n_trees];
    for (int i = 0; i < n_trees; i++)
    {
        numberofbipartition[i] = TreeOPE::sumofdegree(treeset[i]->root, isrooted);
        numberofbipartition[i] = numberofbipartition[i]/2-leaveslabelsmaps.size();
    }
}

void Trees::WriteTrees(string &outfile, Treefileformat tf) // newick nexus
{
    if(tf == NEWICK)
    {
        ofstream outNewick;
        int NN = 10000;
        char *gzbuff = (char*) malloc(NN);
        outNewick.open(outfile.c_str());

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = '\0';
        }

        for(int i = 0; i < n_trees; i++)
        {
            gzbuff[0] = NULL;
            TreeOPE::printTree_new(treeset[i]->root, taxa_str, gzbuff, 0, NN);
            outNewick << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);
        outNewick.close();
    } else
    if(tf == NEXUS)
    {
        ofstream outNexus;
        outNexus.open(outfile.c_str());
        outNexus << "#NEXUS" << endl;
        outNexus << "BEGIN TAXA;" << endl;
        outNexus << "      dimensions ntax=" << leaveslabelsmaps.size() << ";" << endl;
        outNexus << "      taxlabels ";
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            outNexus<< leaveslabelsmaps.name(i) <<" ";
        }
        outNexus << ";" << endl;
        outNexus << "END;" << endl;
        outNexus << "BEGIN TREES;" << endl;
        outNexus << "      dimension ntree=" << n_trees  << ";" << endl;
        outNexus << "      " << "Translate" << endl;
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            if (i < leaveslabelsmaps.size() - 1)
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << "," << endl;
            else
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << endl;
        }
        outNexus << "                ;" << endl;

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = NULL;
        }

        int NN = 10000;
        char *gzbuff = (char*)malloc(NN);

        for(int count = 0; count < n_trees; count++)
        {
            gzbuff[0] = NULL;
            TreeOPE::printTree_nex(treeset[count]->root, leaveslabelsmaps.size(), gzbuff, 0, NN);
            outNexus << "TREE tree " << count << " = "  << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);

        outNexus << "End;" << endl;
        outNexus.close();
    }
}

void Trees::WriteConsensusTree(string &outfile, Treefileformat tf)
{
    if(tf == NEWICK)
    {
        ofstream outNewick;
        int NN = 10000;
        char *gzbuff = (char*) malloc(NN);
        outNewick.open(outfile.c_str());

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = '\0';
        }

        for(int i = 0; i < consensustrees.get_length(); i++)
        {
            gzbuff[0] = NULL;
            TreeOPE::printTree_new(consensustrees[i]->root, taxa_str, gzbuff, 0, NN);
            outNewick << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);
        outNewick.close();
    } else
    if(tf == NEXUS)
    {
        ofstream outNexus;
        outNexus.open(outfile.c_str());
        outNexus << "#NEXUS" << endl;
        outNexus << "BEGIN TAXA;" << endl;
        outNexus << "      dimensions ntax=" << leaveslabelsmaps.size() << ";" << endl;
        outNexus << "      taxlabels ";
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            outNexus<< leaveslabelsmaps.name(i) <<" ";
        }
        outNexus << ";" << endl;
        outNexus << "END;" << endl;
        outNexus << "BEGIN TREES;" << endl;
        outNexus << "      dimension ntree=" << n_trees  << ";" << endl;
        outNexus << "      " << "Translate" << endl;
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            if (i < leaveslabelsmaps.size() - 1)
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << "," << endl;
            else
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << endl;
        }
        outNexus << "                ;" << endl;

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = NULL;
        }

        int NN = 10000;
        char *gzbuff = (char*)malloc(NN);

        for(int count = 0; count < consensustrees.get_length(); count++)
        {
            gzbuff[0] = NULL;
            TreeOPE::printTree_nex(consensustrees[count]->root, leaveslabelsmaps.size(), gzbuff, 0, NN);
            outNexus << "TREE tree " << count << " = "  << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);

        outNexus << "End;" << endl;
        outNexus.close();
    }
}

string Trees::WriteTreesFilename(string fname, string type)
{
    string result;
    string rawname = fname.substr(0, fname.find_last_of("."));
    result = rawname;
    result += "_";

    if(isrooted)
        result += "rooted_";
    else
        result += "unrooted_";

    if (isweighted)
        result += "weighted_";
    else
        result += "unweighted_";

    type.erase(remove( type.begin(), type.end(), ' ' ), type.end() );
    result += type;

    result += ".out";
    return result;
}

string Trees::WriteConsensusTreeFilename(string fname, string type)
{
    string result;
    string rawname = fname.substr(0, fname.find_last_of("."));
    result = rawname;
    result += "_";

    if(isrooted)
        result += "rooted_";
    else
        result += "unrooted_";

    if (isweighted)
        result += "weighted_";
    else
        result += "unweighted_";

    result += "Consensus_tree_";

    type.erase(remove( type.begin(), type.end(), ' ' ), type.end() );
    result += type;

    result += ".out";
    return result;
}

void Trees::Printf(int *idx, int length)
{
    for(int i = 0; i < length; i++)
    {
        TreeOPE::printnewicktree(treeset[idx[i]]);
        cout << endl;
    }
}

//########################ZD comment########################################
//# Generate the hash table and compute the
//# hash value in tree.
//########################ZD comment########################################

void Trees::Compute_Hash()
{
    // Set a random number for m1 (= Initial size of hash table)
    // m1 is the closest value to (t*n)+(t*n*HASHTABLE_FACTOR)
    #define HASHTABLE_FACTOR                       0.2

    // the c value of m2 > c*t*n
    unsigned int C                          =1000;
    int32 NEWSEED                           = 0  ;

    if(!vec_hashrf._hashtab2.empty())
    {
        vec_hashrf.hashrfmap_clear();
    }

    unsigned long long M1 = 0;
    unsigned long long M2 = 0;
    if (NEWSEED != 1000)
    {
        vec_hashrf.uhashfunc_init(n_trees, leaveslabelsmaps.size(), C, NEWSEED);
    }
    else
        vec_hashrf.uhashfunc_init(n_trees, leaveslabelsmaps.size(), C);

    M1 = vec_hashrf._HF.getM1();
    M2 = vec_hashrf._HF.getM2();
    vec_hashrf._hashtab2.resize(M1);

    if(!hash2bitstr.empty())
    {
        hash2bitstr.clear();
    }

    for (unsigned int treeIdx = 0; treeIdx < n_trees; ++treeIdx)
    {
        unsigned int numBitstr = 0;
        TreeOPE::dfs_compute_hash(treeset[treeIdx]->root, leaveslabelsmaps, vec_hashrf, treeIdx,
                                  numBitstr, M1, M2, isweighted, leaveslabelsmaps.size(), hash2bitstr, numberofbipartition[treeIdx]);
    }

    /*map<unsigned long long, Array<char> *>::iterator it = hash2bitstr.begin();
    Array<char> *pt = NULL;
    int n_taxa = leaveslabelsmaps.size();
    for(;it != hash2bitstr.end(); it++)
    {
        cout << "key:" << it->first << endl;//---
        pt = it->second;
        pt->printbits(n_taxa);
        cout << endl;
    }*/
}

//########################ZD comment########################################
//# The arrays of indivial tree's hashvalue, tree index
//# and weight were combined and sorted. Since the hash
//# value represents the unique bipartition, the number
//# of unique bipartion can be counted via checking the
//# hash value. As a result, a sparse bipartition matrix
//# that stores weight of unique bipartition versus
//# trees is created. The "sort" which is different then
//# "Sort" is confusing here. Is it the default sort in c++?
//########################ZD comment########################################

void Trees::Compute_Bipart_Matrix()
{
    int n_taxa = leaveslabelsmaps.size();

    int idex = 0;
    unsigned long long *matrix_hv = new unsigned long long[(n_taxa - 1) * n_trees];
    unsigned int *matrix_treeIdx = new unsigned int[(n_taxa - 1) * n_trees];
    double *matrix_weight = new double[(n_taxa - 1) * n_trees];

    for (unsigned int treeIdx = 0; treeIdx < n_trees; ++treeIdx)
    {
        TreeOPE::bipart(treeset[treeIdx]->root, treeIdx, matrix_hv, matrix_treeIdx, matrix_weight, idex, 0, isrooted);
    }
    Sort(matrix_hv, matrix_treeIdx, matrix_weight, idex);

    unsigned long long *temp = new unsigned long long[idex];
    bipart_count = new unsigned int[idex]; // count numbers of each bipartition
    for (int i = 0; i < idex; i++)
    {
        temp[i] = matrix_hv[i];
    }
    sort(temp, temp+idex);
    unsigned long long *Unique_bipart = new unsigned long long[idex];

    int Unique_idx = 0;
    Unique_bipart[Unique_idx] = temp[0];
    bipart_count[Unique_idx] = 1;   // the first bipartition
    for (int i = 1; i < idex; i++)
    {
        if(temp[i] > temp[i-1])
        {
            Unique_idx++;
            Unique_bipart[Unique_idx] = temp[i];
            bipart_count[Unique_idx] = 1;
         }
         else
            bipart_count[Unique_idx] += 1;
    }
    delete [] temp;
    Unique_idx = Unique_idx + 1;
    cout << "Unique bipartition number: " << Unique_idx << endl;
    treecov_size = Unique_idx;

    //**************bipartition type****************************//
    cout << "Translation of taxa:" << endl;
    for (int i = 0; i < n_taxa; i++)
        cout << leaveslabelsmaps.name(i) << " , " << i+1 << endl;
    map<unsigned long long, Array<char> *>::iterator it = hash2bitstr.begin();
    Array<char> *pt = NULL;
    for(;it != hash2bitstr.end(); it++)
    {
        for (int j = 0; j < Unique_idx; j++)
        {
            if (it->first == Unique_bipart[j])
            {
                cout << "bipartition " << j + 1 << " : ";//---
                pt = it->second;
                pt->printbits(n_taxa);
                cout << ", appear times: " << bipart_count[j] << endl;
            }
            else
                continue;
        }
    }

    //==============Create Sparse Bipartition Matrix==================//
    double *Vals = new double [Unique_idx * n_trees];
    int *RowInds = new int [Unique_idx * n_trees];
    int *ColPtr = new int [n_trees+1];
    int idx = 0;
    int Col_ind = 0;
    ColPtr[Col_ind] = 0;
    if (matrix_treeIdx[0] != 1)
    {
        for (int i = 0; i < matrix_treeIdx[0]-1; i++)
        {
            ColPtr[++Col_ind] = idx;
        }
    }
    for (int j = 0; j < idex; j++)
    {
        for (int i = 0; i < Unique_idx; i++)
        {
            if (matrix_hv[j] == Unique_bipart[i])
            {
                if (isweighted)
                    Vals[idx] = matrix_weight[j];
                else
                {
                    Vals[idx] = 1;
                }
                RowInds[idx] = i;
                if (j > 1 && matrix_treeIdx[j] != matrix_treeIdx[j-1])
                {
                    for (int k = 0; k < (matrix_treeIdx[j] - matrix_treeIdx[j-1]); k++)
                    {
                        ColPtr[++Col_ind] = idx;
                    }
                }
                idx++;
            }
            else
                continue;
        }
    }
    ColPtr[++Col_ind] = idx;
    if (matrix_treeIdx[idex - 1] < n_trees)
    {
        for (int k = 0; k < n_trees - matrix_treeIdx[idex-1]; k++)
        {
            ColPtr[++Col_ind] = idx;
        }
    }
    sbipartmatrix = new SparseMatrix(Unique_idx, n_trees, Vals, RowInds, ColPtr);

//    for(int i = 0; i < Unique_idx; i++)
//    {
//        for(int j = 0; j < n_trees; j++)
//        {
//            cout << (*sbipartmatrix)(i,j) << " ";
//        }
//        cout << endl;
//    }

    delete [] Unique_bipart;
    delete [] matrix_hv;
    delete [] matrix_treeIdx;
    delete [] matrix_weight;
}


string Trees::make_Bipart_Matrix_name(string fname, String format)
{
    string result;
    string rawname = fname.substr(0, fname.find_last_of("."));
    result = rawname;
    result += "_";

    if(isrooted)
        result += "rooted_";
    else
        result += "unrooted_";

    if (isweighted)
        result += "weighted_";
    else
        result += "unweighted_";

    if(format == (String) "List format")
        result += "listFormat_";
    else if(format == (String) "Matrix format")
        result += "matrixFormat_";

    result += "bipartition.out";
    return result;
}

//########################ZD comment########################################
//# Sort weighted edges from the same tree by hash value.
//########################ZD comment########################################

void Trees::Sort(unsigned long long *matrix_hv,
                       unsigned int *matrix_treeIdx,
                       double *matrix_weight, int &idx)
{
    unsigned long long temp_hv;
    double temp_weight;
    for (int i = 0; i < idx; i++)
    {
        for (int j = i+1; j < idx; j++)
        {
            if (matrix_treeIdx[j] == matrix_treeIdx[i])
            {
                if (matrix_hv[i] > matrix_hv[j])
                {
                    temp_hv = matrix_hv[i];
                    matrix_hv[i] = matrix_hv[j];
                    matrix_hv[j] = temp_hv;

                    temp_weight = matrix_weight[i];
                    matrix_weight[i] = matrix_weight[j];
                    matrix_weight[j] = temp_weight;
                }
            }
        }
    }
}

//########################ZD comment########################################
//# This is rank-1 matrix given by two vectors.
//# It is confusing with the SparseMatrix::Multiply_vec
//# and should be integrated in Vector class.
//########################ZD comment########################################

double **Trees::Vec_multiply(const double* Vec1, const double* Vec2, int Unique_idx)
{
    double **result = new double* [Unique_idx];
    for (int i = 0; i < Unique_idx; i++)
        result[i] = new double[Unique_idx];

    for (int i = 0; i < Unique_idx; i++)
    {
        for (int j = 0; j < Unique_idx; j++)
        {
            result[i][j] = Vec1[i] * Vec2[j];
        }
    }
    return result;
}

void Trees::print_bipartitionofonetree(NEWICKNODE*currentnode, bool isrooted, int depth)
{
    depth++;
    int n_taxa = leaveslabelsmaps.size();
    if(depth > 1 && currentnode->Nchildren > 0 && isrooted)
    {
         currentnode->bitstr->printbits(n_taxa);
         cout << endl;
    }
    else if(depth > 2 && currentnode->Nchildren > 0)
    {
         currentnode->bitstr->printbits(n_taxa);
         cout << endl;
    }

    for (int i = 0; i < currentnode->Nchildren; i++)
    {
        if(currentnode->child[i]->Nchildren > 0)
        {
            print_bipartitionofonetree(currentnode->child[i], isrooted, depth);
        }
    }
}

void Trees::Get_bipartitionofonetree(NEWICKNODE*currentnode, bool isrooted, int depth, Array<Array<char> > &bitstrofatree, int &idx)
{
    depth++;
    //int n_taxa = leaveslabelsmaps.size();
    if(depth > 1 && currentnode->Nchildren > 0 && isrooted)
    {
        bitstrofatree[idx] = *(currentnode->bitstr);
        idx++;
    }
    else if(depth > 2 && currentnode->Nchildren > 0)
    {
        bitstrofatree[idx] = *(currentnode->bitstr);
        idx++;
    }

    for (int i = 0; i < currentnode->Nchildren; i++)
    {
        if(currentnode->child[i]->Nchildren > 0)
        {
            Get_bipartitionofonetree(currentnode->child[i], isrooted, depth, bitstrofatree, idx);
        }
    }
}

//########################ZD comment########################################
//# This routine compute the affinity distance, a,
//# from the given distance ,d. The formula is either
//#1/(rel_eps + d)
//# or e^(-d),
//# depending on the flag "type". It accepts unweighted/weighted
//# RF-distance, Matching-distance, SPR-distance or
//# distance given in file.
//########################ZD comment########################################

void Trees::Compute_Affinity_dist(String str_matrix, int type)
{
    double eps = 1000000;
    double ratio = 0.1;
    if (str_matrix == (String)"Unweighted RF-distance")
    {
        switch(type)
        {
        case 1:
            delete_double_array(affi_Recip_URF, n_trees);
            affi_Recip_URF = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Recip_URF[i] = new double [n_trees];

            for(int i = 0; i < n_trees; i++)
                for(int j = 0; j < n_trees; j++)
                    if(dist_URF[i][j] > 0 && eps > dist_URF[i][j])
                        eps = dist_URF[i][j];
            eps = eps * ratio;

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Recip_URF[i][j] = 1.0 / (eps + dist_URF[i][j]);
            }
            break;
        case 2:
            delete_double_array(affi_Exp_URF, n_trees);
            affi_Exp_URF = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Exp_URF[i] = new double [n_trees];

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Exp_URF[i][j] = exp(- dist_URF[i][j]);
            }
            break;
        }
    }
    else if (str_matrix == (String)"Weighted RF-distance")
    {
        switch(type)
        {
        case 1:
            delete_double_array(affi_Recip_RF, n_trees);
            affi_Recip_RF = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Recip_RF[i] = new double [n_trees];

            for(int i = 0; i < n_trees; i++)
                for(int j = 0; j < n_trees; j++)
                    if(dist_RF[i][j] > 0 && eps > dist_RF[i][j])
                        eps = dist_RF[i][j];
            eps = eps * ratio;

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Recip_RF[i][j] = 1.0 / (eps + dist_RF[i][j]);
            }
            break;
        case 2:
            delete_double_array(affi_Exp_RF, n_trees);
            affi_Exp_RF = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Exp_RF[i] = new double [n_trees];

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Exp_RF[i][j] = exp(- dist_RF[i][j]);
            }
            break;
        }
    }
    else if (str_matrix == (String)"Matching-distance")
    {
        switch(type)
        {
        case 1:
            delete_double_array(affi_Recip_match, n_trees);
            affi_Recip_match = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Recip_match[i] = new double [n_trees];

            for(int i = 0; i < n_trees; i++)
                for(int j = 0; j < n_trees; j++)
                    if(dist_match[i][j] > 0 && eps > dist_match[i][j])
                        eps = dist_match[i][j];
            eps = eps * ratio;

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Recip_match[i][j] = 1.0 / (eps + dist_match[i][j]);
            }
            break;
        case 2:
            delete_double_array(affi_Exp_match, n_trees);
            affi_Exp_match = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Exp_match[i] = new double [n_trees];

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Exp_match[i][j] = exp(- dist_match[i][j]);
            }
            break;
        }
    }
    else if (str_matrix == (String)"SPR-distance")
    {
        switch(type)
        {
        case 1:
            delete_double_array(affi_Recip_SPR, n_trees);
            affi_Recip_SPR = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Recip_SPR[i] = new double [n_trees];

            for(int i = 0; i < n_trees; i++)
                for(int j = 0; j < n_trees; j++)
                    if(dist_SPR[i][j] > 0 && eps > dist_SPR[i][j])
                        eps = dist_SPR[i][j];
            eps = eps * ratio;

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Recip_SPR[i][j] = 1.0 / (eps + dist_SPR[i][j]);
            }
            break;
        case 2:
            delete_double_array(affi_Exp_SPR, n_trees);
            affi_Exp_SPR = new double *[n_trees];
            for (int i = 0; i < n_trees; i++)
                affi_Exp_SPR[i] = new double [n_trees];

            for (int i = 0; i < n_trees; i++)
            {
                for (int j = 0; j < n_trees; j++)
                    affi_Exp_SPR[i][j] = exp(- dist_SPR[i][j]);
            }
            break;
        }
    }
    else if(str_matrix == (String)"File-distance")
    {
        switch(type)
        {
        case 1:
            delete_double_array(affi_Recip_file, file_distsize);
            affi_Recip_file = new double *[file_distsize];
            for (int i = 0; i < file_distsize; i++)
                affi_Recip_file[i] = new double [file_distsize];

            for(int i = 0; i < file_distsize; i++)
                for(int j = 0; j < file_distsize; j++)
                    if(dist_file[i][j] > 0 && eps > dist_file[i][j])
                        eps = dist_file[i][j];
            eps = eps * ratio;

            for (int i = 0; i < file_distsize; i++)
            {
                for (int j = 0; j < file_distsize; j++)
                    affi_Recip_file[i][j] = 1.0 / (eps + dist_file[i][j]);
            }
            break;
        case 2:
            delete_double_array(affi_Exp_file, file_distsize);
            affi_Exp_file = new double *[file_distsize];
            for (int i = 0; i < file_distsize; i++)
                affi_Exp_file[i] = new double [file_distsize];

            for (int i = 0; i < file_distsize; i++)
            {
                for (int j = 0; j < file_distsize; j++)
                    affi_Exp_file[i][j] = exp(- dist_file[i][j]);
            }
            break;
        }
    }
}

/*bool Trees::Compute_Matching_dist()
{
    delete_double_array(dist_match, n_trees);
    dist_match = new int *[n_trees];
    for (int i = 0; i < n_trees; i++)
        dist_match[i] = new int [n_trees];

    Ptree tree1, tree2;
    int num_leaves = leaveslabelsmaps.size();

    tree1.lchild = new int [2 * num_leaves];
    tree1.rchild = new int [2 * num_leaves];
    tree1.parent = new int [2 * num_leaves];
    tree1.edge = new int *[2 * num_leaves];
    for(int i = 0; i < num_leaves * 2; i++)
    {
        tree1.edge[i] = new int [2 * num_leaves];
    }
    tree2.lchild = new int [2 * num_leaves];
    tree2.rchild = new int [2 * num_leaves];
    tree2.parent = new int [2 * num_leaves];
    tree2.edge = new int *[2 * num_leaves];
    for(int i = 0; i < num_leaves * 2; i++)
    {
        tree2.edge[i] = new int [2 * num_leaves];
    }

    for (int j = 0; j < n_trees; j++)
    {
        dist_match[j][j] = 0;
        for (int i = 0; i < j; i++)
        {
            if (!TreeOPE::newick2lcbb(treeset[i], leaveslabelsmaps.size(), &tree1))
            {
                cout << "Tree " << i << "  is not a binary tree" << endl;
                for (int k = 0; k < n_trees; k++)
                    delete [] dist_match[k];
                delete [] dist_match;
                return false;
            }

            if (!TreeOPE::newick2lcbb(treeset[j], leaveslabelsmaps.size(), &tree2))
            {
                cout << "Tree " << j << "  is not a binary tree" << endl;
                for (int k = 0; k < n_trees; k++)
                    delete [] dist_match[k];
                delete [] dist_match;
                return false;
            }

            dist_match[i][j] = tree_mmdis(&tree1, &tree2, leaveslabelsmaps.size());
            dist_match[j][i] = dist_match[i][j];
        }

    }

    delete [] tree1.lchild;
    delete [] tree1.rchild;
    delete [] tree1.parent;
    for(int i = 0; i < num_leaves * 2; i++)
    {
        delete [] tree1.edge[i];
    }
    delete [] tree1.edge;

    delete [] tree2.lchild;
    delete [] tree2.rchild;
    delete [] tree2.parent;
    for(int i = 0; i < num_leaves * 2; i++)
    {
        delete [] tree2.edge[i];
    }
    delete [] tree2.edge;

    return true;
}*/

//########################ZD comment########################################
//# It constructs the edge matrix of a Ptree,
//# which should be implemented in TreeOPE.
//########################ZD comment########################################

void Trees:: pttree(struct Ptree *treeA, int node)
{
    if (((*treeA).lchild[node] == -1) &&  ((*treeA).rchild[node] == -1))
    {
        if (node != 0)
        {
            for(int i = 0; i < (*treeA).leaf_number; i++)
            {
                (*treeA).edge[node][i] = 0;
            }
            (*treeA).edge[node][node] = 1;
        }
        else
        {
            for(int i = 0; i < (*treeA).leaf_number; i++)
            {
                (*treeA).edge[node][i] = 1;
            }
            (*treeA).edge[node][node] = 0;
        }
    }
    else
    {
        if ((*treeA).lchild[node] != -1)
            pttree(treeA, (*treeA).lchild[node]);
        if ((*treeA).rchild[node] != -1)
            pttree(treeA, (*treeA).rchild[node]);
        for(int i = 0; i < (*treeA).leaf_number; i++)
        {
            (*treeA).edge[node][i] = ( (*treeA).edge[(*treeA).lchild[node]][i] + (*treeA).edge[(*treeA).rchild[node]][i] ) % 2;
        }
    }
}

//########################ZD comment########################################
//# This distance is given by the solution
//# of Hungarian algorithm of the cost matrix,
//# r, given by compute_matrix.
//########################ZD comment########################################

int Trees::tree_mmdis(struct Ptree *tree1, struct Ptree *tree2, int num_leaf)
{
    int* r;
    int** m;
    int mmdis;
    hungarian_problem_t p;

    int root = num_leaf;
    pttree(tree1,root);
    pttree(tree2,root);
    r = (int *)malloc(sizeof(int)*((root - 3)*(root - 3)));
    if (r == NULL){
        printf("malloc r failed!\n"); getchar();
    }
    compute_matrix(r,root-3, tree1, tree2);
    m = array_to_matrix(r,root-3,root-3);
    /* initialize the gungarian_problem using the cost matrix*/
    hungarian_init(&p, m , root-3, root-3, HUNGARIAN_MODE_MINIMIZE_COST) ;
    /* solve the assignement problem */
    mmdis = hungarian_solve(&p);
    /* free used memory */
    hungarian_free(&p);
    for(int i=0;i<root-3;i++)
    {
        free(m[i]);
    }
    free(m);
    free(r);
    return(mmdis);
}


//########################ZD comment########################################
//# It accumulates the number common edges from
//# two trees and store in a vectorized matrix, r.
//########################ZD comment########################################

void Trees::compute_matrix(int *r, int range, struct Ptree *tree1, struct Ptree *tree2)
{
    int e1,e2;
    int temp1, temp2;
    int row_r = 0, col_r = 0;
    int spe1, spe2;

    if ((*tree1).lchild[(*tree1).leaf_number] > (*tree1).leaf_number)
        spe1 = (*tree1).lchild[(*tree1).leaf_number];
    else
        spe1 = (*tree1).rchild[(*tree1).leaf_number];
    if ((*tree2).lchild[(*tree2).leaf_number] > (*tree2).leaf_number)
        spe2 = (*tree2).lchild[(*tree2).leaf_number];
    else
        spe2 = (*tree2).rchild[(*tree2).leaf_number];

    for(e1 = (*tree1).leaf_number + 1; e1 < 2 * (*tree1).leaf_number - 1; e1++)
    {
        if (e1 != spe1)
        {
            col_r = 0;
            for(e2 = (*tree2).leaf_number + 1; e2 < 2 * (*tree2).leaf_number - 1; e2++)
            {
                if (e2 != spe2)
                {
                    temp1 = 0;
                    temp2 = 0;
                    for(int i=0; i < (*tree2).leaf_number; i++)
                    {
                        if ( (*tree1).edge[e1][i] != (*tree2).edge[e2][i] )
                            temp1++;
                    }
                    for(int i=0; i < (*tree2).leaf_number; i++)
                    {
                        if ( (*tree1).edge[e1][i] == (*tree2).edge[e2][i] )
                            temp2++;
                    }

                    if (temp1 < temp2)
                        r[row_r*range+col_r] = temp1;
                    else
                        r[row_r*range+col_r] = temp2;
                    col_r++;
                }
            }
            if ( col_r !=  (*tree1).leaf_number - 3 )
            {
                cout << "Overflow of " << col_r << "!" << endl;
                getchar();
            }
            row_r++;
        }
    }

    if (row_r !=  (*tree1).leaf_number - 3)
    {
        cout << "Overflow of row_r " <<  row_r << "!" << endl;
        getchar();
    }
}

//########################ZD comment########################################
//# Convert a vectorizerd matrix to 2 dimensional
//# array, which should not be here.
//########################ZD comment########################################

int ** Trees:: array_to_matrix (int* m, int rows, int cols)
{
    int i,j;
    int** r;
    r = (int**)calloc(rows,sizeof(int*));
    if (r == NULL)
    {
        cout << "Alloc r failed!" << endl;
        getchar();
    }
    for(i=0;i<rows;i++)
    {
        r[i] = (int*)calloc(cols,sizeof(int));
        if (r[i] == NULL)
        {
            cout << "Alloc r failed!" << endl;
            getchar();
        }
        for(j=0;j<cols;j++)
        {
            r[i][j] = m[i*cols+j];
        }
    }
    return r;
}

void Trees::OutputBipartitionMatrix(std::ostream &output, SparseMatrixOutputType smtype)
{
    sbipartmatrix->OutputSparseMatrix(output, smtype);
}

void Trees::load_distfile(string fname)
{
    ifstream distfile(fname.c_str(), ios::binary);

    if (!distfile)
    {
        cout << "Unable to open the file!\n\n";
        exit(0);
    }

    file_distsize = 0;
    string line;
    while (!distfile.eof())
    {
        getline(distfile, line);
        if (!line.empty())
        {
            file_distsize++;
        }
    }
    file_distsize--;

    delete_double_array(dist_file, file_distsize);
    dist_file = new double *[file_distsize];
    double index;
    distfile.clear();
    distfile.seekg(0, ios::beg);
    distfile >> line;
    for (int i = 0; i < file_distsize; i++)
        distfile >> index;
    for (int i = 0; i < file_distsize; i++)
    {
        distfile >> index;
        dist_file[i] = new double [file_distsize];
        for (int j = 0; j <= i; j++)
        {
            distfile >> dist_file[i][j];
            dist_file[j][i] = dist_file[i][j];
        }
    }
};

void Trees::load_coordinatefile(string fname)
{
    ifstream coordfile(fname.c_str(), ios::binary);

    if (!coordfile)
    {
        cout << "Unable to open the file!\n\n";
        exit(0);
    }

    file_coordinatesize = 0;
    string line;
    while (!coordfile.eof())
    {
        getline(coordfile, line);
        if (!line.empty())
        {
            file_coordinatesize++;
        }
    }

    coordfile.clear();
    coordfile.seekg(0, ios::beg);
    file_coordinatedim = 0;
    int pos = 0;
    bool pre_is_table = true;
    getline(coordfile, line);
    for(int i = 0; i < line.size(); i++)
    {
        pos = line.find('\t', pos);
        if(i == pos)
        {
            if(!pre_is_table)
            {
                file_coordinatedim++;
                pos += 1;
            }
            pre_is_table = true;
        }
        else
            pre_is_table = false;
    }
    if(!pre_is_table)
        file_coordinatedim++;

    delete_double_array(coord_file, file_coordinatesize);
    coord_file = new double *[file_coordinatesize];
    for(int i = 0; i < file_coordinatesize; i++)
        coord_file[i] = new double [file_coordinatedim];

    coordfile.clear();
    coordfile.seekg(0, ios::beg);
    for (int i = 0; i < file_coordinatesize; i++)
    {
        for (int j = 0; j < file_coordinatedim; j++)
        {
            coordfile >> coord_file[i][j];
        }
    }
};

void Trees::load_affinityfile(string fname)
{
    ifstream afile(fname.c_str(), ios::binary);

    if (!afile)
    {
        cout << "Unable to open the file!\n\n";
        exit(0);
    }

    affinityfile_size = 0;
    string line;
    while (!afile.eof())
    {
        getline(afile, line);
        if (!line.empty())
        {
            affinityfile_size++;
        }
    }
    affinityfile_size--;

    delete_double_array(fileaffinity, affinityfile_size);
    fileaffinity = new double *[affinityfile_size];
    double index;
    afile.clear();
    afile.seekg(0, ios::beg);
    afile >> line;
    for (int i = 0; i < affinityfile_size; i++)
        afile >> index;
    for (int i = 0; i < affinityfile_size; i++)
    {
        afile >> index;
        fileaffinity[i] = new double [affinityfile_size];
        for (int j = 0; j <= i; j++)
        {
            afile >> fileaffinity[i][j];
            fileaffinity[j][i] = fileaffinity[i][j];
        }
    }
};

void Trees::load_covariancefile(string fname)
{
    ifstream covfile(fname.c_str(), ios::binary);

    if (!covfile)
    {
        cout << "Unable to open the file!\n\n";
        exit(0);
    }

    filecov_size = 0;
    string line;
    while (!covfile.eof())
    {
        getline(covfile, line);
        if (!line.empty())
        {
            filecov_size++;
        }
    }
    filecov_size--;

    delete_double_array(filecov, filecov_size);
    filecov = new double *[filecov_size];
    double index;
    covfile.clear();
    covfile.seekg(0, ios::beg);
    covfile >> line;
    for (int i = 0; i < filecov_size; i++)
        covfile >> index;
    for (int i = 0; i < filecov_size; i++)
    {
        covfile >> index;
        filecov[i] = new double [filecov_size];
        for (int j = 0; j <= i; j++)
        {
            covfile >> filecov[i][j];
            filecov[j][i] = filecov[i][j];
        }
    }
};



bool Trees::bipartmatrixIsexisting()
{
    if(sbipartmatrix == NULL)
        return false;

    return true;
}


bool Trees::covarianceMatrixIsexisting()
{
    if(treecov == NULL)
        return false;

    return true;
}


bool Trees::treesAreexisting()
{
    if(treeset == NULL)
        return false;

    return true;
}

bool Trees::consensusTreeIsexisting()
{
    if(consensustrees.is_empty() == 1)
        return false;

    return true;
}

string Trees::create_temp_name(String str_matrix)
 {
     char* rawname = (char *)str_matrix;
     string result(rawname);
     result += "_temp.txt";
     return result;
 }

 template<class T>
 void Trees::print_comm_array(T *** arr, int n, string outfile, bool arr_is_covariance, double highfreq, double lowfreq)
 {
     ofstream fout;
     if (outfile != "")
     fout.open(outfile.c_str());

//    if(arr_is_covariance)
//    {
        covariance_freeid = new int[n];
        covariance_freeid_size = 0;
        covariance_nonfree_id = new int[n];
        covariance_nonfree_id_size = 0;
//    }

   if((*arr) != NULL)
   {
        if(arr_is_covariance)
        {
            for (int i = 0; i < n; i++)
            {
                bool crit;
                if(n_trees == 0)
                    crit = ((*arr)[i][i] <= lowfreq * (1 - lowfreq) || (*arr)[i][i] <= highfreq * (1 - highfreq));
                else
                    crit = (((double) bipart_count[i] / n_trees) <= lowfreq || ((double) bipart_count[i] / n_trees) >= highfreq);

                if(crit)   //if((*arr)[i][i] != 0)
                {
                    covariance_freeid[covariance_freeid_size] = i;
                    covariance_freeid_size++;
                }
                else
                {
                    covariance_nonfree_id[covariance_nonfree_id_size] = i;
                    covariance_nonfree_id_size++;
                }
            }
            fout << ">" << endl;
            for (int i = 0; i < covariance_nonfree_id_size; i++)
            {
                fout << i << " " << i << endl;
//                fout << covariance_nonfree_id[i] << " " << covariance_nonfree_id[i] << endl;
            }
            fout << ">" << endl;
            fout << "1 1" << endl;
            fout << ">" << endl;
            for (int i = 0; i < covariance_nonfree_id_size; i++)
            {
                for (int j = 0; j < covariance_nonfree_id_size; j++)
                {
                    fout << j  << " " << i << " " << (*arr)[covariance_nonfree_id[j]][covariance_nonfree_id[i]] << " " << "1" << endl;
//                    fout << covariance_nonfree_id[j]  << " " << covariance_nonfree_id[i] << " " << (*arr)[covariance_nonfree_id[j]][covariance_nonfree_id[i]] << " " << "1" << endl;
                }
            }
        }
        else
        {
            fout << ">" << endl;
            for (int i = 0; i < n; i++)
            {
                fout << i << " " << i << endl;
            }
            fout << ">" << endl;
            fout << "1 1" << endl;
            fout << ">" << endl;

            for(int i = 0; i < n; i++)
                covariance_nonfree_id[i] = i;
            covariance_nonfree_id_size = n;

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    fout << j << " " << i << " " << (*arr)[i][j] << " " << "1" << endl;
        }
    }
    if (outfile != "")
        fout.close();
 }

 void Trees::print_community_file(String str_matrix, string outfile, double highfreq, double lowfreq)
 {
     if (str_matrix == (String)"Covariance Matrix")
        print_comm_array((double ***) StrToDist[str_matrix], treecov_size, outfile, true, highfreq, lowfreq);
     else if (str_matrix == (String)"File-covariance")
        print_comm_array((double ***) StrToDist[str_matrix], filecov_size, outfile, true, highfreq, lowfreq);
     else if(str_matrix == (String)"Affinity-Reciprocal-filedist")
        print_comm_array((double ***) StrToDist[str_matrix], file_distsize, outfile, false, highfreq, lowfreq);
     else if(str_matrix == (String)"Affinity-Exponential-filedist")
        print_comm_array((double ***) StrToDist[str_matrix], file_distsize, outfile, false, highfreq, lowfreq);
     else if(str_matrix == (String)"File-affinity")
        print_comm_array((double ***) StrToDist[str_matrix], affinityfile_size, outfile, false, highfreq, lowfreq);
     else
        print_comm_array((double ***) StrToDist[str_matrix], n_trees, outfile, false, highfreq, lowfreq);
 }

 string Trees::create_out_name(String str_matrix)
 {
     char* rawname = (char *)str_matrix;
     string result(rawname);
     result += ".bin";
     return result;
 }

 string Trees::create_node_name(String str_matrix)
 {
     char* rawname = (char *)str_matrix;
     string result(rawname);
     result += "_node_map.txt";
     return result;
 }

 string Trees::create_conf_name(String str_matrix)
 {
     char* rawname = (char *)str_matrix;
     string result(rawname);
     result += ".conf";
     return result;
 }

 int Trees::read_conf(char* filename, int* &conf, int* &sign)
 {
     ifstream finput;
     finput.open(filename, fstream::in | fstream::binary);
     if (finput.fail())
     {
         cout << "Could not read config file: " << filename << "\n\n";
         exit(-1);
     }

     int nb_layers;
     //read number of layers
     finput.read((char *)&nb_layers, sizeof(int));
     conf = (int*)malloc((long)nb_layers*sizeof(int));
     sign = (int*)malloc((long)nb_layers*sizeof(int));

     // read conf
     finput.read((char *)conf, nb_layers*sizeof(int));
     if (finput.fail() || finput.eof())
     {
         cout << "Error!\n\n";
         exit(-1);
         //throw new exception();
     }
     finput.read((char *)sign, nb_layers*sizeof(int));

     return nb_layers;
 }

 void Trees::create_resolution(double lp, double ln, int nb_layers, int* sign, double* &lambda)
 {
     for (int i = 0; i < nb_layers; i++)
     {
         if (sign[i] == POSITIVE)
             lambda[i] = lp;
         else if (sign[i] == NEGATIVE)
             lambda[i] = ln;
         else
             cout << "Error: Layer not recognized!\n\n";
     }
 }

 string Trees::create_comm_name(String str_matrix)
 {
    char* rawname = (char *)str_matrix;
    string result(rawname);
    result += "_community.out";
    delete rawname;

    return result;
 }


 bool Trees::compute_community_automatically(String str_matrix, int modelType, string highfreq, string lowfreq)
  {
     if(com_info != NULL)
         delete_double_array(com_info, covariance_nonfree_id_size + 4);
     srand(time(NULL));
     covariance_freeid_size = 0;
     if(covariance_freeid != NULL)
         delete [] covariance_freeid;
     covariance_freeid = NULL;
     covariance_nonfree_id_size = 0;
     if(covariance_nonfree_id != NULL)
         delete [] covariance_nonfree_id;
     covariance_nonfree_id = NULL;
     string temp_file = create_temp_name(str_matrix);
     double highfrequence = atof(highfreq.c_str());
     double lowfrequence = atof(lowfreq.c_str());
     if ((str_matrix == (String)"Covariance Matrix" || str_matrix == (String) "File-covariance") && (highfrequence > 1.0 || highfrequence < 0.0
         || lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
     {
         cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
         return false;
     }
     print_community_file(str_matrix, temp_file, highfrequence, lowfrequence);
     char *infile = strdup(temp_file.c_str());
     string outf = create_out_name(str_matrix);
     char *outfile = strdup(outf.c_str());
     string nodef = create_node_name(str_matrix);
     char *node_map_file = strdup(nodef.c_str());
     string conff = create_conf_name(str_matrix);
     char *conf_file = strdup(conff.c_str());
     int is_weighted = 1;
     int is_directed = 1;
     int is_single_slice = 0;
     double interslice_weight = 1.0;

     int* conf = NULL;
     int* sign = NULL;
     double* lambda = NULL;

     Slicer s(infile, (double)interslice_weight, is_directed, is_weighted, is_single_slice);
     Graph* g = s.get_graph();
     g->display_binary(outfile);
     delete g;
     s.display_node_mapping(node_map_file);
     s.display_conf(conf_file, modelType);

     int layers = read_conf(conf_file, conf, sign);
     lambda = new double[layers];


     //----------------fixed lambda neg, find two lambda_+ ----------
     //----------------such that the numbers of community  ----------
     //----------------are minimum or maximum.             ----------

     double lambda_neg = 0;
     double lambda_pos_min = -1, lambda_pos_max = 1;
     Community * community;
     map<double, Community *> LamCommunities;
     map<double, Community *>::iterator it, it2;
     map<double, double> mods;
     int times = 0;
     int numNodes = 0;

     while(times <= 20)
     {
         create_resolution(lambda_pos_min, lambda_neg, layers, sign, lambda);
         community = new Community(outfile, conf, sign, lambda);
         int stochastic = 0;
         GreedyLouvain::iterate_randomly = stochastic;
         GreedyLouvain::detect_communities(community);
         if(community->nb_comm == 1)
         {
             numNodes = community->g->nb_nodes;
             mods[lambda_pos_min] = community->modularity();
             LamCommunities[lambda_pos_min] = community;
             break;
         } else
         {
             delete community;
             community = NULL;
             lambda_pos_min *= 2;
         }
         times++;
     }
     if(times == 11)
     {
         cout << "Error: Cannot find lambda pos min!\n\n";
         return false;
     }

     times = 0;
     while(times <= 20)
     {
         create_resolution(lambda_pos_max, lambda_neg, layers, sign, lambda);
         community = new Community(outfile, conf, sign, lambda);
         int stochastic = 0;
         GreedyLouvain::iterate_randomly = stochastic;
         GreedyLouvain::detect_communities(community);

         if(str_matrix == (String) "Affinity-Reciprocal-URF" || str_matrix == (String) "Affinity-Exponential-URF"
                 || str_matrix == (String) "Affinity-Reciprocal-RF" || str_matrix == (String) "Affinity-Exponential-RF"
                 || str_matrix == (String) "Affinity-Reciprocal-match" || str_matrix == (String) "Affinity-Exponential-match"
                 || str_matrix == (String) "Affinity-Reciprocal-SPR" || str_matrix == (String) "Affinity-Exponential-SPR"
                 || str_matrix == (String) "Affinity-Reciprocal-geodesic" || str_matrix == (String) "Affinity-Exponential-geodesic"
                 || str_matrix == (String) "Affinity-Reciprocal-filedist" || str_matrix == (String) "Affinity-Exponential-filedist"
                  || str_matrix == (String) "File-affinity")
         {
             bool allsametopo = true;
             double **affimatrix = *((double ***) StrToDist[(String) "Unweighted RF-distance"]);
             if(affimatrix != NULL)
             {
                 double samevalue = affimatrix[0][0];
                 for(int i = 0; i < community->nb_comm; i++)
                 {
                     int repidx = -1;
                     for(int j = 0; j < covariance_nonfree_id_size; j++)
                     {
                         if(community->n2c[j] == i)
                         {
                             if(repidx == -1)
                             {
                                 repidx = j;
                             } else
                             {
                                 if(fabs(affimatrix[repidx][j] - samevalue) >= 1e-10)
                                 {
                                     allsametopo = false;
                                 }
                             }
                         }
                     }
                 }
                 if(allsametopo)
                 {
                     mods[lambda_pos_max] = community->modularity();
                     LamCommunities[lambda_pos_max] = community;
                     break;
                 }
             }
         }
         if(community->nb_comm == covariance_nonfree_id_size)
         {
             mods[lambda_pos_max] = community->modularity();
             LamCommunities[lambda_pos_max] = community;
             break;
         } else
         {
             delete community;
             community = NULL;
             lambda_pos_max *= 2;
         }
         times++;
     }
     if(times == 11)
     {
         cout << "Error: Cannot find lambda pos max!\n\n";
         return false;
     }

     cout << "lambda neg: " << lambda_neg << endl;
     cout << "lambda pos min: " << lambda_pos_min << ", lambda pos max: " << lambda_pos_max << endl;

     queue<pair<double, double> > qq;
     pair<double, double> p, pp;
     vector<pair<double, double> > plateausLb(0);
     vector<pair<double, double> > plateausUb(0);
     bool SameAsFirst, SameAsSecond;
     double lambda_pos, plateaubound = 0;
     int atleastsearchnum = 0;
     p = make_pair(lambda_pos_min, lambda_pos_max);
     qq.push(p);

     // width first search
     cout << "\nThe testing values of lambda pos are:" << endl;
     cout << lambda_pos_min << ", " << lambda_pos_max << ", ";
     while(!qq.empty())
     {
         pp = qq.front();
         qq.pop();
         lambda_pos = (pp.first + pp.second) / 2;
         cout << lambda_pos << ", ";
         std::cout.flush();
         create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
         community = new Community(outfile, conf, sign, lambda);
         int stochastic = 0;
         GreedyLouvain::iterate_randomly = stochastic;
         GreedyLouvain::detect_communities(community);

         SameAsFirst = false;
         if(lambda_pos_min == pp.first)
         {
             if(community->nb_comm == 1)
             {
                 SameAsFirst = true;
                 lambda_pos_min = lambda_pos;
             }
         } else
         {
             it = LamCommunities.find(pp.first);
             if(it != LamCommunities.end())
             {
                 if(Community::IsSameCommunity(community, it->second))
                     SameAsFirst = true;
             }
         }
         SameAsSecond = false;
         if(lambda_pos_max == pp.second)
         {
             if(community->nb_comm == covariance_nonfree_id_size)
             {
                 SameAsSecond = true;
                 lambda_pos_max = lambda_pos;
             }
         } else
         {
             it = LamCommunities.find(pp.second);
             if(it != LamCommunities.end())
             {
                 if(Community::IsSameCommunity(community, it->second))
                     SameAsSecond = true;
             }
         }
         mods[lambda_pos] = community->modularity();
         LamCommunities[lambda_pos] = community;

         // added new interval for testing
         if(!SameAsFirst && ((lambda_pos - pp.first) > plateaubound || plateaubound == 0 || LamCommunities.size() + qq.size() < atleastsearchnum))
         {
             qq.push(make_pair(pp.first, lambda_pos));
         }
         if(!SameAsSecond && ((pp.second - lambda_pos) > plateaubound || plateaubound == 0 || LamCommunities.size() + qq.size() < atleastsearchnum))
         {
             qq.push(make_pair(lambda_pos, pp.second));
         }

         // find plateaus
         if(SameAsFirst && community->nb_comm != 1)
         {
             if(plateaubound == 0)
                 plateaubound = (lambda_pos - pp.first) / 2;
             int i;
             for(i = 0; i < plateausLb.size(); i++)
             {
                 if(plateausLb[i].second == pp.first)
                 {
                     plateausLb[i].second = lambda_pos;
                     break;
                 }
             }
             if(i == plateausLb.size())
             {
                 plateausLb.resize(plateausLb.size() + 1);
                 plateausLb[plateausLb.size() - 1] = make_pair(pp.first, lambda_pos);
             }
         }
         if(SameAsSecond && community->nb_comm != covariance_nonfree_id_size)
         {
             if(plateaubound == 0)
                 plateaubound = (pp.second - lambda_pos) / 2;
             int i;
             for(i = 0; i < plateausLb.size(); i++)
             {
                 if(plateausLb[i].first == pp.second)
                 {
                     plateausLb[i].first = lambda_pos;
                     break;
                 }
             }
             if(i == plateausLb.size())
             {
                 plateausLb.resize(plateausLb.size() + 1);
                 plateausLb[plateausLb.size() - 1] = make_pair(lambda_pos, pp.second);
             }
         }
     }
     cout << endl;

     cout << "\nThe found plateaus are:" << endl;
     for(int i = 0; i < plateausLb.size(); i++)
         cout << "[" << plateausLb[i].first << ", " << plateausLb[i].second
              << "], length: " << plateausLb[i].second - plateausLb[i].first << ", number of communities: " <<
                 LamCommunities[plateausLb[i].first]->nb_comm << endl;

     cout << "\nDetailed check lambda pos are (if necessary):" << endl;

     plateausUb.resize(plateausLb.size());
     for(int i = 0; i < plateausLb.size(); i++)
         plateausUb[i] = make_pair(plateausLb[i].first, plateausLb[i].second);

     double extstart, extend, max_length = 0;
     vector<int> totestidix(0);

     do{
         totestidix.resize(0);

         for(int i = 0; i < plateausUb.size(); i++)
         {
             extstart = plateausLb[i].first - (plateausLb[i].second - plateausLb[i].first);
             extend = plateausLb[i].second + (plateausLb[i].second - plateausLb[i].first);
             for(it = LamCommunities.begin(); it != LamCommunities.end(); it++)
             {
                 if(it->first > extstart && it->first < plateausLb[i].first)
                     extstart = it->first;
                 if(it->first < extend && it->first > plateausLb[i].second)
                     extend = it->first;
             }
             plateausUb[i].first = extstart;
             plateausUb[i].second = extend;

             if(plateausLb[i].second - plateausLb[i].first > max_length)
             {
                 max_length = plateausLb[i].second - plateausLb[i].first;
             }
         }

         for(int i = 0; i < plateausLb.size(); i++)
         {
             if(((plateausUb[i].second - plateausUb[i].first) - max_length) >= 0.000001)
             {
                 totestidix.resize(totestidix.size() + 1);
                 totestidix[totestidix.size() - 1] = i;
             }
         }
         if(totestidix.size() <= 1)
             break;
         else
         {
             for(int i = 0; i < totestidix.size(); i++)
             {
                 int idx = totestidix[i];

                 lambda_pos = (plateausUb[idx].first + plateausLb[idx].first) / 2;
                 if(LamCommunities.find(lambda_pos) == LamCommunities.end())
                 {
                     cout << lambda_pos << ", ";
                     std::cout.flush();
                     create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
                     community = new Community(outfile, conf, sign, lambda);
                     int stochastic = 0;
                     GreedyLouvain::iterate_randomly = stochastic;
                     GreedyLouvain::detect_communities(community);
                     mods[lambda_pos] = community->modularity();
                     LamCommunities[lambda_pos] = community;
                 }
                 it = LamCommunities.find(lambda_pos);
                 it2 = LamCommunities.find(plateausLb[idx].first);
                 if(Community::IsSameCommunity(it2->second, it->second))
                     plateausLb[idx].first = lambda_pos;
                 else
                     plateausUb[idx].first = lambda_pos;
                 lambda_pos = (plateausLb[idx].second + plateausUb[idx].second) / 2;
                 if(LamCommunities.find(lambda_pos) == LamCommunities.end())
                 {
                     cout << lambda_pos << ", ";
                     std::cout.flush();
                     create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
                     community = new Community(outfile, conf, sign, lambda);
                     int stochastic = 0;
                     GreedyLouvain::iterate_randomly = stochastic;
                     GreedyLouvain::detect_communities(community);
                     mods[lambda_pos] = community->modularity();
                     LamCommunities[lambda_pos] = community;
                 }
                 it = LamCommunities.find(lambda_pos);
                 it2 = LamCommunities.find(plateausLb[idx].second);
                 if(Community::IsSameCommunity(it2->second, it->second))
                     plateausLb[idx].second = lambda_pos;
                 else
                     plateausUb[idx].second = lambda_pos;
             }
             cout << endl;
         }
     } while(true);

     cout << "\nThe found plateaus are:" << endl;
     for(int i = 0; i < plateausLb.size(); i++)
     {
//         cout << "[" << plateausUb[i].first << ":" << plateausLb[i].first << ", " << plateausLb[i].second
//              << ":" << plateausUb[i].second << "], length lower bound:" <<
//                 plateausLb[i].second - plateausLb[i].first << ", length upper bound:" <<
//                 plateausUb[i].second - plateausUb[i].first << ", number of communities:" <<
//                 LamCommunities[plateausLb[i].first]->nb_comm << endl;

         cout << "Transition lower bound: [" << plateausUb[i].first << ", " << plateausLb[i].first << "]" << endl;
         cout << "Transition upper bound: [" << plateausLb[i].second << ", " << plateausUb[i].second << "]" << endl;
         cout << "Updated plateau: [" << plateausLb[i].first << ", " << plateausLb[i].second << "], length: " <<
                 plateausLb[i].second - plateausLb[i].first << ", number of communities: " <<
                 LamCommunities[plateausLb[i].first]->nb_comm << endl << endl;
     }

     com_info_col = LamCommunities.size() + 1;
     com_info = new double *[covariance_nonfree_id_size + 4];
     for(int i = 0; i < covariance_nonfree_id_size + 4; i++)
         com_info[i] = new double [com_info_col];

     com_info[0][0] = numNodes;

     com_info[1][0] = lambda_neg;
     com_info[2][0] = 0;
     com_info[3][0] = 0;
     for(int i = 0; i < covariance_nonfree_id_size; i++)
         com_info[i + 4][0] = covariance_nonfree_id[i];

     // compute labels of communities
     int k = 0;
     com_info[0][1] = 0;
     it = LamCommunities.begin();
     community = it->second;

     com_info[1][1] = it->first;
     com_info[2][1] = it->second->nb_comm;
     com_info[3][1] = mods[it->first];//---it->second->modularity();
     for(int i = 4; i < covariance_nonfree_id_size + 4; i++)
         com_info[i][1] = it->second->n2c[i - 4];

     it++;
     for(; it != LamCommunities.end(); it++)
     {
         k++;
         if(Community::IsSameCommunity(community, it->second)) // problem
             com_info[0][k + 1] = com_info[0][k];
         else
             com_info[0][k + 1] = com_info[0][k] + 1;
         community = it->second;
         com_info[1][k + 1] = it->first;
         com_info[2][k + 1] = it->second->nb_comm;
         com_info[3][k + 1] = mods[it->first];
         for(int i = 4; i < covariance_nonfree_id_size + 4; i++)
             com_info[i][k + 1] = it->second->n2c[i - 4];
     }
 /*
     int maxcommlength = 0;
     int currentcommlength = 1;
     int plateauidx = 1;
     for(int i = 2; i < lambdasize + 1; i++)
     {
         if(com_info[0][i] == com_info[0][i - 1])
             currentcommlength++;
         else
             currentcommlength = 1;
         if(currentcommlength > maxcommlength)
         {
             maxcommlength = currentcommlength;
             plateauidx = i;
         }
     }
     for(int i = 0; i < covariance_nonfree_id_size + 4; i++)
     {
         com_info[i][lambdasize + 1] = com_info[i][plateauidx];
     }
 */

     // output plateaus
#ifdef COMMAND_LINE_VERSION
     string stdnamepla = treesfilename.c_str();
     string stdoutfname = treesfilename.c_str();
#else
     string stdnamepla = commfilename.c_str();
     string stdoutfname = commfilename.c_str();
#endif
     stdnamepla = stdnamepla.substr(0, stdnamepla.find_last_of("."));
     stdnamepla += "_";

     if (str_matrix != (String)"File-covariance" && str_matrix != (String) "Affinity-Reciprocal-filedist"
             && str_matrix != (String) "Affinity-Exponential-filedist" && str_matrix != (String) "File-affinity")
     {
        if(isrooted)
            stdnamepla += "rooted_";
        else
            stdnamepla += "unrooted_";

        if(isweighted)
            stdnamepla += "weighted_";
        else
            stdnamepla += "unweighted_";
     }

     if(str_matrix == (String) "Unweighted RF-distance")
         stdnamepla += "RF-distance";
     else if(str_matrix == (String) "Weighted RF-distance")
         stdnamepla += "RF-distance";
     else if(str_matrix == (String) "Covariance Matrix")
         stdnamepla += "Covariance_Matrix";
     else
         stdnamepla += str_matrix;

     if(modelType == 3)
         stdnamepla += "_Configuration_Null_Model";
     else if(modelType == 4)
         stdnamepla += "_Constant_Potts_Model";
     else if(modelType == 2)
         stdnamepla += "_Erdos-Renyi_Null_Model";
     else if(modelType == 1)
         stdnamepla += "_No_Null_Model";

     stdnamepla += "_communities_auto_plateaus.out";

     String fnamepla(stdnamepla.c_str());
     File filepla(fnamepla);
     filepla.clean();
     if(! filepla.is_open())
     {
         cout << "Unable to open the file: " << fnamepla << "\n\n";
         return false;
     }
     filepla << "lengthLB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausLb[i].second - plateausLb[i].first << "\t";
     filepla << "\n";
     filepla << "lengthUB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausUb[i].second - plateausUb[i].first << "\t";
     filepla << "\n";
     filepla << "startUB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausUb[i].first << "\t";
     filepla << "\n";
     filepla << "startLB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausLb[i].first << "\t";
     filepla << "\n";
     filepla << "endLB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausLb[i].second << "\t";
     filepla << "\n";
     filepla << "endUB:\t";
     for(int i = 0; i < plateausLb.size(); i++)
         filepla << plateausUb[i].second << "\t";
     filepla << "\n";
     if (str_matrix == (String) "Covariance Matrix")
        filepla << "Bipartition index" << "\t" << "Community Index" << "\n";
     else
         filepla << "Tree index" << "\t" << "Community Index" << "\n";
     for(int i = 0; i < covariance_nonfree_id_size; i++)
     {
         filepla << covariance_nonfree_id[i] << "\t";
         for(int j = 0; j < plateausLb.size(); j++)
         {
             filepla << LamCommunities[plateausLb[j].first]->n2c[i] << "\t";
         }
         filepla << "\n";
     }

     // output communities information
     stdoutfname = stdoutfname.substr(0, stdoutfname.find_last_of("."));
     stdoutfname += "_";

     if (str_matrix != (String)"File-covariance" && str_matrix != (String) "Affinity-filedist"
             && str_matrix != (String) "File-affinity")
     {
         if(isrooted)
             stdoutfname += "rooted_";
         else
             stdoutfname += "unrooted_";

         if(isweighted)
             stdoutfname += "weighted_";
         else
             stdoutfname += "unweighted_";
     }

     if(str_matrix == (String) "Unweighted RF-distance")
         stdoutfname += "RF-distance";
     else if(str_matrix == (String) "Weighted RF-distance")
         stdoutfname += "RF-distance";
     else if(str_matrix == (String) "Covariance Matrix")
         stdoutfname += "Covariance_Matrix";
     else
        stdoutfname += str_matrix;

     if(modelType == 3)
         stdoutfname += "_Configuration_Null_Model";
     else if(modelType == 4)
         stdoutfname += "_Constant_Potts_Model";
     else if(modelType == 2)
         stdoutfname += "_Erdos-Renyi_Null_Model";
     else if(modelType == 1)
         stdoutfname += "_No_Null_Model";

     stdoutfname += "_community_auto_results.out";

     String outfname(stdoutfname.c_str());
     File file(outfname);
     file.clean();
     if(! file.is_open())
     {
         cout << "Unable to open the file: " << outfname << "\n\n";
         return false;
     }

     for(int i = 0; i < covariance_nonfree_id_size + 4; i++)
     {
         if (i == 0)
         {
             if (str_matrix == (String) "Covariance Matrix")
                file << "Same community as previous or not (first number is number of bipartitions)" << "\n";
             else
                file << "Same community as previous or not (first number is number of trees)" << "\n";
         }
         else if (i == 1)
             file << "Value of lambda: "<< "\n";
         else if (i == 2)
             file << "Number of communities: " << "\n";
         else if (i == 3)
             file << "Value of modularity: " << "\n";
         else if (i == 4)
         {
             if (str_matrix == (String) "Covariance Matrix")
                file << "Community index (first column is bipartition index): " << "\n";
             else
                 file << "Community index (first column is tree index): " << "\n";
         }
         for(int j = 0; j < com_info_col; j++)
             file << com_info[i][j] << "\t";
         file << "\n";
     }

     //delete temporary file
     free(infile);
     free(outfile);
     free(node_map_file);
     free(conf_file);

     free(conf);
     free(sign);
     if (lambda != NULL)
         delete [] lambda;

     for(it = LamCommunities.begin(); it != LamCommunities.end(); it++)
         delete it->second;

     const char *tempfile0 = temp_file.c_str();
     remove(tempfile0);
     const char *tempfile1 = outf.c_str();
     remove(tempfile1);
     const char *tempfile2 = nodef.c_str();
     remove(tempfile2);
     const char *tempfile3 = conff.c_str();
     remove(tempfile3);

     cout << "Output community results to file: " << outfname << endl;
     cout << "and " << fnamepla << "\n\n";
     return true;
  }


bool Trees::compute_community_fixedlambda(String str_matrix, int modelType, double lambdapos, double lambdaneg, string highfreq, string lowfreq) // lambda negative = 0
{
    cout << "Compute community for new lambda pair: lambdapos = " << lambdapos << ", lambdaneg = " << lambdaneg << endl;

    srand(time(NULL));
    string temp_file = create_temp_name(str_matrix);
    double highfrequence = atof(highfreq.c_str());
    double lowfrequence = atof(lowfreq.c_str());
    if ((str_matrix == (String)"Covariance Matrix" || str_matrix == (String) "File-covariance") && (highfrequence > 1.0 || highfrequence < 0.0
        || lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
    {
        cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
        return false;
    }
    print_community_file(str_matrix, temp_file, highfrequence, lowfrequence);
    char *infile = strdup(temp_file.c_str());
    string outf = create_out_name(str_matrix);
    char *outfile = strdup(outf.c_str());
    string nodef = create_node_name(str_matrix);
    char *node_map_file = strdup(nodef.c_str());
    string conff = create_conf_name(str_matrix);
    char *conf_file = strdup(conff.c_str());
    int is_weighted = 1;
    int is_directed = 1;
    int is_single_slice = 0;
    double interslice_weight = 1.0;

    int* conf = NULL;
    int* sign = NULL;
    double* lambda = NULL;

    Slicer s(infile, (double)interslice_weight, is_directed, is_weighted, is_single_slice);
    Graph* g = s.get_graph();
    g->display_binary(outfile);
    delete g;
    s.display_node_mapping(node_map_file);
    s.display_conf(conf_file, modelType);

    int layers = read_conf(conf_file, conf, sign);
    lambda = new double[layers];

    Community *community = NULL;
    create_resolution(lambdapos, lambdaneg, layers, sign, lambda);
    community = new Community(outfile, conf, sign, lambda);
    int stochastic = 0;
    GreedyLouvain::iterate_randomly = stochastic;
    GreedyLouvain::detect_communities(community);

    cout << "Number of communities is " << community->nb_comm
         << "and modularity is " << community->modularity() << endl;

    int com_info_row = covariance_nonfree_id_size + 4;
    int com_info_col_old = com_info_col;
    double **com_info_old = com_info;
    double lambdavar;

    com_info_col++;
    com_info = new double*[com_info_row];
    for(int i = 0; i < com_info_row; i++)
        com_info[i] = new double [com_info_col];

    lambdavar = (com_info_old[1][0] == lambdaneg) ? lambdapos : lambdaneg;
    for(int j = 0; j < com_info_row; j++)
        com_info[j][0] = com_info_old[j][0];
    int newidx = com_info_col_old;
    for(int i = 1; i < com_info_col_old; i++)
    {
        if(com_info_old[1][i] < lambdavar)
            for(int j = 0; j < com_info_row; j++)
                com_info[j][i] = com_info_old[j][i];
        else
        {
            newidx = i;
            break;
        }
    }
    bool sameaspre = true;
    if(newidx > 1)
    {
        for(int i = 0; i < covariance_nonfree_id_size; i++)
        {
            if(community->n2c[i] != com_info_old[i + 4][newidx - 1])
            {
                sameaspre = false;
                break;
            }
        }
    }
    bool sameaslatter = true;
    if(newidx != com_info_col_old)
    {
        for(int i = 0; i < covariance_nonfree_id_size; i++)
        {
            if(community->n2c[i] != com_info_old[i + 4][newidx])
            {
                sameaslatter = false;
                break;
            }
        }
    }

    if(newidx == 1)
        com_info[0][newidx] = 0;
    else if(sameaspre)
        com_info[0][newidx] = com_info[0][newidx - 1];
    else
        com_info[0][newidx] = com_info[0][newidx - 1] + 1;
    com_info[1][newidx] = lambdavar;
    com_info[2][newidx] = community->nb_comm;
    com_info[3][newidx] = community->modularity();
    for(int i = 4; i < com_info_row; i++)
        com_info[i][newidx] = community->n2c[i-4];
    if(newidx != com_info_col_old)
    {
        if(sameaslatter)
            com_info[0][newidx + 1] = com_info[0][newidx];
        else
            com_info[0][newidx + 1] = com_info[0][newidx] + 1;
        for(int i = 1; i < com_info_row; i++)
            com_info[i][newidx + 1] = com_info_old[i][newidx];
    }
    for(int i = newidx + 2; i < com_info_col; i++)
    {
        com_info[0][i] = com_info[0][i - 1] + (com_info_old[0][i - 1] - com_info_old[0][i - 2]);
        for(int j = 0; j < com_info_row; j++)
            com_info[j][i] = com_info_old[j][i - 1];
    }

//    for(int i = 0; i < com_info_row; i++)
//    {
//        for(int j = 0; j < com_info_col; j++)
//            cout << com_info[i][j] << "\t";
//        cout << endl;
//    }

    free(infile);
    free(outfile);
    free(node_map_file);
    free(conf_file);

    if(com_info_old != NULL)
        delete_double_array(com_info_old, com_info_row);

    delete community;
    free(conf);
    free(sign);
    if (lambda != NULL)
        delete [] lambda;

    const char *tempfile0 = temp_file.c_str();
    remove(tempfile0);
    const char *tempfile1 = outf.c_str();
    remove(tempfile1);
    const char *tempfile2 = nodef.c_str();
    remove(tempfile2);
    const char *tempfile3 = conff.c_str();
    remove(tempfile3);
}

bool Trees::compute_community_manually(String str_matrix, int modelType, Array<double> param1, Array<double> param2, string highfreq, string lowfreq)
{
    if(com_info != NULL)
        delete_double_array(com_info, covariance_nonfree_id_size + 4);
    srand(time(NULL));
    covariance_freeid_size = 0;
    covariance_freeid = NULL;
    covariance_nonfree_id_size = 0;
    covariance_nonfree_id = NULL;
    string temp_file = create_temp_name(str_matrix);
    double highfrequence = atof(highfreq.c_str());
    double lowfrequence = atof(lowfreq.c_str());
   if ((str_matrix == (String)"Covariance Matrix" || str_matrix == (String) "File-covariance") && (highfrequence > 1.0 || highfrequence < 0.0
        || lowfrequence > 1.0 || lowfrequence < 0.0 || (highfrequence - lowfrequence) <= 0.0))
    {
        cout << "Warning: The high and low frequencies must be between 0 and 1!\n\n";
        return false;
    }

    print_community_file(str_matrix, temp_file, highfrequence, lowfrequence);
    char *infile = strdup(temp_file.c_str());
    string outf = create_out_name(str_matrix);
    char *outfile = strdup(outf.c_str());
    string nodef = create_node_name(str_matrix);
    char *node_map_file = strdup(nodef.c_str());
    string conff = create_conf_name(str_matrix);
    char *conf_file = strdup(conff.c_str());
    int is_weighted = 1;
    int is_directed = 1;
    int is_single_slice = 0;
    double interslice_weight = 1.0;

    int* conf = NULL;
    int* sign = NULL;
    double* lambda = NULL;

    Slicer s(infile, (double)interslice_weight, is_directed, is_weighted, is_single_slice);
    Graph* g = s.get_graph();
    g->display_binary(outfile);
    delete g;
    s.display_node_mapping(node_map_file);
    s.display_conf(conf_file, modelType);

    int layers = read_conf(conf_file, conf, sign);
    lambda = new double[layers];
    double lambda_pos;// = atof(param1.c_str());
    double lambda_neg;// = atof(param2.c_str());

#ifdef COMMAND_LINE_VERSION
     string stdoutfname = treesfilename.c_str();
#else
     string stdoutfname = stdoutfname.c_str();
#endif
    stdoutfname = stdoutfname.substr(0, stdoutfname.find_last_of("."));
    stdoutfname += "_";

    if (str_matrix != (String)"File-covariance" && str_matrix != (String) "Affinity-Reciprocal-filedist"
            && str_matrix != (String) "Affinity-Exponential-filedist" && str_matrix != (String) "File-affinity")
    {
        if(isrooted)
            stdoutfname += "rooted_";
        else
            stdoutfname += "unrooted_";

        if(isweighted)
            stdoutfname += "weighted_";
        else
            stdoutfname += "unweighted_";
    }

    if(str_matrix == (String) "Unweighted RF-distance")
        stdoutfname += "RF-distance";
    else if(str_matrix == (String) "Weighted RF-distance")
        stdoutfname += "RF-distance";
    else if(str_matrix == (String) "Covariance Matrix")
        stdoutfname += "Covariance_Matrix";
    else
       stdoutfname += str_matrix;

    if(modelType == 3)
        stdoutfname += "_Configuration_Null_Model";
    else if(modelType == 4)
        stdoutfname += "_Constant_Potts_Model";
    else if(modelType == 2)
        stdoutfname += "_Erdos-Renyi_Null_Model";
    else if(modelType == 1)
        stdoutfname += "_No_Null_Model";

    stdoutfname += "_community_manual_results.out";

    String outfname(stdoutfname.c_str());
    File file(outfname);
    file.clean();
    if(! file.is_open())
    {
        cout << "Unable to open the file: " << outfname << "\n\n";
        return false;
    }
    int lambdasize = 0;

    if(param1.get_length() == 1)
        lambdasize = param2.get_length();
    else if(param2.get_length() == 1)
        lambdasize = param1.get_length();
    else
    {
        cout << "Error: The length of the array of lambda negative or the length of the array lambda positive should be one!\n\n";
        return false;
    };
    com_info_col = lambdasize + 1;

    com_info = new double *[covariance_nonfree_id_size + 4];
    for(int i = 0; i < covariance_nonfree_id_size + 4; i++)
        com_info[i] = new double [com_info_col];

    Community** communities = new Community *[lambdasize];

    for(int k = 0; k < lambdasize; k++)
    {
        if(param1.get_length() == 1)
        {
            lambda_pos = param1[0];
            lambda_neg = param2[k];
        } else
        {
            lambda_pos = param1[k];
            lambda_neg = param2[0];
        }
        if (modelType != 1)
        {
            cout << "Using lambda+ " << lambda_pos << ", lambda- " << lambda_neg << endl;
        }
        create_resolution(lambda_pos, lambda_neg, layers, sign, lambda);
        communities[k] = new Community(outfile, conf, sign, lambda);
        int stochastic        = 0;
        GreedyLouvain::iterate_randomly = stochastic;
        if(stochastic)
            cout << "Using random node order" << endl;
        cout << "Network has "
        << communities[k]->g->nb_nodes << " nodes, "
        << communities[k]->g->nb_links << " links, " << endl;

        GreedyLouvain::detect_communities(communities[k]);
        double mod = communities[k]->modularity();

        cout << "Value of modularity:" << mod << endl;
        cout << "Number of communities: " << communities[k]->nb_comm << endl;
        //co->display_comm2node();

        if(covariance_freeid != NULL && covariance_freeid_size != 0)
        {
            for (int i = 0; i < communities[k]->nb_comm; i++)
            {
                std::cout << "Community " << i + 1 << " includes nodes: ";
                for (int j = 0; j < communities[k]->size; j++)
                {
                    if(communities[k]->n2c[j] == i)
                        std::cout << covariance_nonfree_id[j] << ",";
                }
                std::cout << "\n";
            }
            std::cout << "Free node index:" << endl;
            for (int i = 0; i < covariance_freeid_size; i++)
                std::cout << covariance_freeid[i] << ", ";
            std::cout << "\n";
        }
        else
        {
            for (int i = 0; i < communities[k]->nb_comm; i++)
            {
                std::cout << "Community " << i + 1 << " includes nodes: ";
                for (int j = 0; j < communities[k]->size; j++)
                {
                    if(communities[k]->n2c[j] == i)
                        std::cout << j << ",";
                }
                std::cout << "\n";
            }
        }

        if(k == 0)
        {
            com_info[0][0] = communities[k]->g->nb_nodes;

            if(param1.get_length() == 1)
                com_info[1][0] = param1[0];
            else
                com_info[1][0] = param2[0];
            com_info[2][0] = 0;
            com_info[3][0] = 0;
            for(int i = 0; i < covariance_nonfree_id_size; i++)
                com_info[i + 4][0] = covariance_nonfree_id[i];
        }

        if(k == 0)
            com_info[0][1] = 0;
        else
        {
            if(Community::IsSameCommunity(communities[k - 1], communities[k]))
                com_info[0][k + 1] = com_info[0][k];
            else
                com_info[0][k + 1] = com_info[0][k] + 1;
        }
        if(param1.get_length() == 1)
            com_info[1][k + 1] = param2[k];
        else
            com_info[1][k + 1] = param1[k];

        com_info[2][k + 1] = communities[k]->nb_comm;
        com_info[3][k + 1] = mod;

        for(int i = 4; i < covariance_nonfree_id_size + 4; i++)
            com_info[i][k + 1] = communities[k]->n2c[i - 4];
    }

    for(int i = 0; i < covariance_nonfree_id_size + 4; i++)
    {
        if (i == 0)
        {
            if (str_matrix == (String) "Covariance Matrix")
               file << "Same community as previous or not (first number is number of bipartitions)" << "\n";
            else
               file << "Same community as previous or not (first number is number of trees)" << "\n";
        }
        else if (i == 1)
            file << "Value of lambda: "<< "\n";
        else if (i == 2)
            file << "Number of communities: " << "\n";
        else if (i == 3)
            file << "Value of modularity: " << "\n";
        else if (i == 4)
        {
            if (str_matrix == (String) "Covariance Matrix")
               file << "Community index (first column is bipartition index): " << "\n";
            else
                file << "Community index (first column is tree index): " << "\n";
        }
        for(int j = 0; j < lambdasize + 1; j++)
            file << com_info[i][j] << "\t";
        file << "\n";
    }

    //output into a file
    /*string out_comm = create_comm_name(str_matrix);
    char *out_filename;
    out_filename = strdup(out_comm.c_str());
    ostream* comm_output = &cout;
    ofstream* foutput = new ofstream();
    if (out_filename != NULL)
    {
        foutput->open(out_filename,fstream::out);
        comm_output = foutput;
    }
    co->display_partition(*comm_output);

    if (out_filename != NULL)
    {
        foutput->close();
    }
    delete foutput;*/

    //delete temporary file
    free(infile);
    free(outfile);
    free(node_map_file);
    free(conf_file);

    free(conf);
    free(sign);
    if (lambda != NULL)
        delete [] lambda;
    for(int i = 0; i < lambdasize; i++)
    {
        delete communities[i];
    }
    delete communities;
    const char *tempfile0 = temp_file.c_str();
    remove(tempfile0);
    const char *tempfile1 = outf.c_str();
    remove(tempfile1);
    const char *tempfile2 = nodef.c_str();
    remove(tempfile2);
    const char *tempfile3 = conff.c_str();
    remove(tempfile3);

    cout << "Output community results to file: " << outfname << "\n\n";
    return true;
 }

const NEWICKTREE *Trees::get_tree(int idx)
{
    if(idx < 0 || idx >= n_trees)
    {
        std::cout << idx << "-th tree does not exist in the tree set!" << std::endl << endl;
        return NULL;
    }

    return treeset[idx];
}

bool Trees::compute_consensus_tree(ConsensusTree type, const char *listname)
{
    deleteConsensustree();

    int num = idxlist.get_length();

    map<unsigned long long, unsigned long long> bipcount;
    for(int i = 0; i < idxlist.get_length(); i++)
    {
        int idx = idxlist[i];
        if(i >= n_trees || i < 0)
        {
            std::cout << "Warning:" << n_trees << " trees in data set! " << i << " is out of range of trees!" << std::endl << endl;
            continue;
        }
        TreeOPE::obtainbipartcount(treeset[idx], isrooted, bipcount);
    }

    map<unsigned long long, unsigned long long>::iterator iter;
    map<unsigned long long, Array<char> *>::iterator h2biter;
    Array<char> *contreebitstr = new Array<char> [leaveslabelsmaps.size()];
    Array<double> confreq(leaveslabelsmaps.size());
    Array<unsigned long long> conhash(leaveslabelsmaps.size());
    unsigned int conbipnum = 0;

    for(iter = bipcount.begin(); iter != bipcount.end(); iter++)
    {
        if( (type == MAJORITYTREE && iter->second > (double) 0.5 * num) ||
            (type == STRICTTREE && iter->second == num))
        {
            h2biter = hash2bitstr.find(iter->first);
            if(h2biter != hash2bitstr.end())
            {
                contreebitstr[conbipnum] = *(h2biter->second);
                confreq[conbipnum] = (double) (iter->second) / num;
                conhash[conbipnum] = iter->first;
                conbipnum++;
            }
        }
    }
    std::cout << "Obtain " << conbipnum << " bipartitions in the consensus tree!" << std::endl << endl;
    int contreenum = 0;//---consensustrees.get_length();
    consensustrees.resize(contreenum + 1);
    consensuslist.resize(contreenum + 1);
    consensuslist[contreenum] = listname;

    if(!TreeOPE::buildconsensustree(consensustrees[contreenum],
                                    confreq,
                                    conhash,
                                    contreebitstr,
                                    conbipnum,
                                    leaveslabelsmaps.size()))
    {
        std::cout << "Error: Incompatible bipartitions! Unable to compute consensus tree!\n\n";
        delete [] contreebitstr;
        return false;
    }

    delete [] contreebitstr;
    return true;
}

//void Trees::Compute_Cumulative(double** bipartFreq, int* bipartFreqIdx, int splits, int increments)
//{
////    if(bipartFreq != NULL)
////        delete_double_array(bipartFreq, splits);

//    int Unique_idx = treecov_size;
//    int burnIn = round(n_trees * 0.1);       // 10% burn-in
//    int subIntvls = ceil((n_trees - burnIn) / increments);

//    // save bipart_count into a temporary array
//    unsigned int* tmp = NULL;
//    tmp = new unsigned int[Unique_idx];
//    for(int i = 0; i < Unique_idx; i++)
//        tmp[i] = bipart_count[i];

//    // determine bipartitions with max number of frequencies
//    int* maxVal = NULL;
//    maxVal = new int[splits];
//    for(int i = 0; i < splits; i++)
//        maxVal[i] = 0;

//    for(int i = 0; i < splits; i++)
//    {
//        maxVal[i] = tmp[0];
//        for(int j = 0; j < Unique_idx; j++)
//        {
//            if(maxVal[i] < tmp[j])
//            {
//                maxVal[i] = tmp[j];
//                bipartFreqIdx[i] = j;
//            }
//        }
//        tmp[bipartFreqIdx[i]] = 0;
//    }

//    // compute relative cumulative frequency
//    int endSubIntvl;
//    int k = 0;
//    double total = subIntvls;
//    int** tmpBipartFreq = NULL;
//    tmpBipartFreq = new int*[splits];
//    for(int i = 0; i < splits; i++)
//        tmpBipartFreq[i] = new int[increments];

//    for(int i = 0; i < splits; i++)
//        for(int j = 0; j < increments; j++)
//            tmpBipartFreq[i][j] = 0;

//    cout << "Burning in 10% of the trees." << endl;
//    for(int j = burnIn; j < n_trees; j += subIntvls)
//    {
//        for(int i = 0; i < splits; i++)
//        {
//            if(k > 0)
//                tmpBipartFreq[i][k] += tmpBipartFreq[i][k-1];

//            if((j + subIntvls) > n_trees || (j + subIntvls) == (n_trees))
//                endSubIntvl = n_trees;
//            else
//                endSubIntvl = j + subIntvls;

//            for(int l = j; l < endSubIntvl; l++)
//            {
//                if((*sbipartmatrix)(bipartFreqIdx[i],l) != 0)
//                {
//                    tmpBipartFreq[i][k] += 1;
//                }
//            }
//            bipartFreq[i][k] = (double) ((tmpBipartFreq[i][k]/total)*100);
//        }
//        if((j + subIntvls) > n_trees || (j + subIntvls) == (n_trees))
//            total += (n_trees - burnIn);
//        else
//            total += subIntvls;
//        k += 1;
//        if(k >= increments)
//            break;
//    }

//    // output
////    cout << "Cumuative sum for max " << splits << " bipart freq: " << endl;
////    for(int i = 0; i < splits; i++)
////    {
////        for(int j = 0; j < increments; j++)
////            cout << bipartFreq[i][j] << " ";
////        cout << endl;
////    }
////    cout << endl;


//    // free allocated memory
//    for(int i = 0; i < splits; i++)
//        delete [] tmpBipartFreq[i];
//    delete [] tmpBipartFreq;
//    delete [] tmp;
//    delete [] maxVal;
//}

//void Trees::Compute_Slide(double** bipartFreq, int* bipartFreqIdx, int splits, int increments)
//{
//    int Unique_idx = treecov_size;
//    int burnIn = round(n_trees*0.1);       // 10% burn-in
//    int subIntvls = ceil((n_trees - burnIn) / increments);

////    cout << "Bipartition Freq: ";
////    for(int i = 0; i < Unique_idx; i++)
////        cout << bipart_count[i] << " ";
////    cout << endl;

//    // save bipart_count into a temporary array
//    unsigned int* tmp = NULL;
//    tmp = new unsigned int[Unique_idx];
//    for(int i = 0; i < Unique_idx; i++)
//        tmp[i] = bipart_count[i];

//    // determine bipartitions with max number of frequencies
//    int* maxVal = NULL;
//    maxVal = new int[splits];
//    for(int i = 0; i < splits; i++)
//        maxVal[i] = 0;

////    int* bipartFreqIdx = NULL;
////    bipartFreqIdx = new int[splits];
////    for(int i = 0; i < splits; i++)
////        bipartFreqIdx[i] = 0;

//    for(int i = 0; i < splits; i++)
//    {
//        maxVal[i] = tmp[0];
//        for(int j = 0; j < Unique_idx; j++)
//        {
//            if(maxVal[i] < tmp[j])
//            {
//                maxVal[i] = tmp[j];
//                bipartFreqIdx[i] = j;
//            }
//        }
//        tmp[bipartFreqIdx[i]] = 0;
//    }

//    cout << "Max " << splits << " Bipart Vals: ";
//    for(int i = 0; i < splits; i++)
//        cout << maxVal[i] <<" ";
//    cout << endl;

////    cout << "Max " << splits << " Bipart Vals Indices: ";
////    for(int i = 0; i < splits; i++)
////        cout << bipartFreqIdx[i] <<" ";
////    cout << endl;

//    // compute frequency for interval
//    int endSubIntvl;
//    int k = 0;
//    cout << "Burning in 10% of the trees." << endl;
//    for(int j = burnIn; j < n_trees; j += subIntvls)
//    {
//        for(int i = 0; i < splits; i++)
//        {
//            if((j + subIntvls) > n_trees || (j + subIntvls) == (n_trees - 1))
//                endSubIntvl = n_trees;
//            else
//                endSubIntvl = j + subIntvls;

//            for(int l = j; l < endSubIntvl; l++)
//            {
////                cout << "sbipartmatrix( " << bipartFreqIdx[i] << " , " << l << " ) = " << (*sbipartmatrix)(bipartFreqIdx[i],l) << endl;
//                if((*sbipartmatrix)(bipartFreqIdx[i],l) != 0)
//                    bipartFreq[i][k] += 1;
//            }
//        }
//        k += 1;
//        if(k >= increments)
//            break;
//    }

//    // output
////    cout << "Cumuative sum for max " << splits << " bipart freq: " << endl;
////    for(int i = 0; i < splits; i++)
////    {
////        for(int j = 0; j < increments; j++)
////            cout << bipartFreq[i][j] << " ";
////        cout << endl;
////    }
////    cout << endl;


//    // free allocated memory
//    delete [] tmp;
//    delete [] maxVal;
//}

#ifndef COMMAND_LINE_VERSION
QString Trees::get_treefilename_without_path()
{
    QString qtname = treesfilename.c_str();
    QStringList list = qtname.split("/");
    return list.value(list.length() - 1);
}
#endif

string Trees::Print_selected_indices()
{
#ifndef COMMAND_LINE_VERSION
    QString qtname = get_treefilename_without_path();
    string stdname = qtname.toStdString();
#else
    string stdname = treesfilename;
#endif

    String outputname(stdname.c_str());

    int i = 0;
    while(outputname[i] != '.' && outputname[i] != '\0')
        i++;

    outputname = outputname(0, i);
    outputname += (String) "_selected_trees_indices.out";

    stdname = (char *) outputname;

    File fout(outputname);
    fout.clean();

    for(int i = 0; i < selected_trees.get_length(); i++)
    {
        if(selected_trees[i] == -1)
            break;
        fout << selected_trees[i] << endl;
    }

    return stdname;
}


string Trees::Print_selected_trees(Treefileformat tf)
{
#ifndef COMMAND_LINE_VERSION
    QString qtname = get_treefilename_without_path();
    string stdname = qtname.toStdString();
#else
    string stdname = treesfilename;
#endif

    String outputname(stdname.c_str());

    int i = 0;
    while(outputname[i] != '.' && outputname[i] != '\0')
        i++;

    outputname = outputname(0, i);

    if(tf == NEWICK)
    {
        outputname += (String) "_selected_trees_NEWICK.out";
        stdname = (char *) outputname;

        WriteSelectedTrees(stdname, NEWICK);
    }
    else
    if(tf == NEXUS)
    {
        outputname += (String) "_selected_trees_NEXUS.out";
        stdname = (char *) outputname;

        WriteSelectedTrees(stdname, NEXUS);
    }

    return stdname;
}

string Trees::WriteSelectedTreesFilename(string type)
{
#ifndef COMMAND_LINE_VERSION
    QString qtname = get_treefilename_without_path();
    string outputname = qtname.toStdString();
#else
    string outputname = treesfilename;
#endif

    outputname = outputname.substr(0, outputname.find_last_of("."));

//    int i = 0;
//    while(outputname[i] != '.' && outputname[i] != '\0')
//        i++;

//    outputname = outputname(0, i);
    outputname += "_selected_trees_";
    outputname += type.c_str();
    outputname += ".out";
//    stdname = (char *) outputname;

    return outputname;
}

void Trees::WriteSelectedTrees(string &outfile, Treefileformat tf) // newick nexus
{
    if(tf == NEWICK)
    {
        ofstream outNewick;
        int NN = 10000;
        char *gzbuff = (char*) malloc(NN);
        outNewick.open(outfile.c_str());

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = '\0';
        }

        for(int i = 0; i < selected_trees.get_length(); i++)
        {
            if(selected_trees[i] == -1)
                break;

            gzbuff[0] = NULL;
            TreeOPE::printTree_new(treeset[selected_trees[i]]->root, taxa_str, gzbuff, 0, NN);
            outNewick << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);
        outNewick.close();
    } else
    if(tf == NEXUS)
    {
        int n_subtrees;
        for(int i = 0; i < selected_trees.get_length(); i++)
        {
            if(selected_trees[i] == -1)
                break;
            n_subtrees++;
        }

        ofstream outNexus;
        outNexus.open(outfile.c_str());
        outNexus << "#NEXUS" << endl;
        outNexus << "BEGIN TAXA;" << endl;
        outNexus << "      dimensions ntax=" << leaveslabelsmaps.size() << ";" << endl;
        outNexus << "      taxlabels ";
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            outNexus<< leaveslabelsmaps.name(i) <<" ";
        }
        outNexus << ";" << endl;
        outNexus << "END;" << endl;
        outNexus << "BEGIN TREES;" << endl;
        outNexus << "      dimension ntree=" << n_subtrees  << ";" << endl;
        outNexus << "      " << "Translate" << endl;
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            if (i < leaveslabelsmaps.size() - 1)
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << "," << endl;
            else
                outNexus << "                " << i+1 << " " << leaveslabelsmaps.name(i) << endl;
        }
        outNexus << "                ;" << endl;

        char **taxa_str = new char*[leaveslabelsmaps.size()];
        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            string temp = leaveslabelsmaps.name(i);
            taxa_str[i] = new char[temp.length()+1];
            for (int j = 0; j < temp.length(); j++)
                taxa_str[i][j] = temp[j];
            taxa_str[i][temp.length()] = NULL;
        }

        int NN = 10000;
        char *gzbuff = (char*)malloc(NN);

        for(int count = 0; count < selected_trees.get_length(); count++)
        {
            if(selected_trees[count] == -1)
                break;

            gzbuff[0] = NULL;
            TreeOPE::printTree_nex(treeset[selected_trees[count]]->root, leaveslabelsmaps.size(), gzbuff, 0, NN);
            outNexus << "TREE tree " << selected_trees[count] << " = "  << gzbuff << ";" << endl;
        }

        for (int i = 0; i < leaveslabelsmaps.size(); i++)
        {
            delete[] taxa_str[i];
            taxa_str[i] = NULL;
        }
        delete[] taxa_str;
        taxa_str = NULL;

        free(gzbuff);

        outNexus << "End;" << endl;
        outNexus.close();
    }
}



#ifdef COMMAND_LINE_VERSION

//########################ZD comment########################################
//# Compute the bipartition covariance matrix from
//# the matrix, C, created by Compute_Bipart_Matrix,
//# M: M1 = MM^T, v1 = mean(M), v2 = sum(M),
//# M2=v2v1^T, M3 = v1v1^T, C = (M1-M2-M2^T+n*M3)/(n-1).
//# Note that it is implemented via sparse matrix-vector
//# multiplication.
//########################ZD comment########################################

void Trees::Compute_Bipart_Covariance()
{
    SparseMatrix* trans = sbipartmatrix->transpose();
    SparseMatrix* result = sbipartmatrix->Multiply(*trans);
    double* M = sbipartmatrix->Mean(n_trees);
    double* Ones = new double[n_trees];
    for (int i = 0; i < n_trees; i++)
        Ones[i] = 1;
    double* tmp = sbipartmatrix->Multiply_vec(Ones);
    int Unique_idx = treecov_size;
    double ** temp1 = Vec_multiply(tmp, M, Unique_idx);
    double ** temp2 = Vec_multiply(M, M, Unique_idx);
    delete_double_array(treecov, treecov_size);
    treecov = new double*[Unique_idx];
    for (int i = 0; i < Unique_idx; i++)
        treecov[i] = new double[Unique_idx];

    for (int i = 0; i < Unique_idx; i++)
        for (int j = 0; j < Unique_idx; j++)
            treecov[i][j] = 0;

   for (int i = 0; i < Unique_idx; i++)
    {
        for (int j = 0; j < Unique_idx; j++)
        {
            treecov[i][j] = treecov[i][j] + (*result)(i,j) - temp1[i][j] - temp1[j][i] + n_trees * temp2[i][j];
            treecov[i][j] = treecov[i][j] / (n_trees - 1);
        }
    }

    delete [] M;
    M = NULL;

    delete [] Ones;
    Ones = NULL;

    delete [] tmp;
    tmp = NULL;

    for (int i = 0; i < Unique_idx; i++)
    {
            delete [] temp1[i];
            delete [] temp2[i];
    }
    delete [] temp1;
    temp1 = NULL;
    delete [] temp2;
    temp2 = NULL;

    if(trans != NULL)
    {
        delete trans;
        trans = NULL;
    }
    if(result != NULL)
    {
        delete result;
        result = NULL;
    }
}

//########################ZD comment########################################
//# Compute the unweighted/weighted RF distance.
//# For the unweighted distance, accumulate the
//# number of each unique bipartition's presence
//# in each tree, f_ij, and the # of bipartitions,
//# n_i,then d_ij = (n_i+n_j-2f_ij)/2. For weighted
//# case, it is more complicated.
//########################ZD comment########################################


bool Trees::Compute_RF_dist_by_hash(bool ISWEIGHTED)
{
    if(!isweighted && ISWEIGHTED)
    {
        cout << "Warning: The trees in memory are unweighted! Unable to compute weighted RF distances!\n\n";
        return false;
    }

    bool bUbid = false; // for counting the number pf unique bipartitions
    unsigned long uBID = 0;

    if (ISWEIGHTED == false)
    {
        delete_double_array(dist_URF, n_trees);
        dist_URF = new double *[n_trees];
        for (int i = 0; i < n_trees; i++)
            dist_URF[i] = new double [n_trees];
        for(int i = 0; i < n_trees; i++)
        {
            for(int j = 0; j < n_trees; j++)
                dist_URF[i][j] = 0;
        }
        for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti)
        {
            unsigned int sizeVec = vec_hashrf._hashtab2[hti].size();
            if (sizeVec)
            {
                uBID += sizeVec;
                if (!bUbid)
                {
                    for (unsigned int i = 0; i < sizeVec; ++i)
                    {
                        unsigned int sizeTreeIdx = vec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
                        if (sizeTreeIdx > 1)
                        {
                            for (unsigned int j = 0; j < sizeTreeIdx; ++j)
                            {
                                for (unsigned int k = 0; k < sizeTreeIdx; ++k)
                                {
                                    if (j == k)
                                        continue;
                                    else
                                        dist_URF[vec_hashrf._hashtab2[hti][i]._vec_treeidx[j]][vec_hashrf._hashtab2[hti][i]._vec_treeidx[k]] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < n_trees; i++)
        {
            for (int j = 0; j < i; j++)
            {
                dist_URF[i][j] = (double) ((numberofbipartition[i] + numberofbipartition[j] - 2 * dist_URF[i][j]) / 2);
                dist_URF[j][i] = dist_URF[i][j];
            }
        }
    }
    else if (ISWEIGHTED == true)
    {
        delete_double_array(dist_RF, n_trees);
        dist_RF = new double *[n_trees];
        for (int i = 0; i < n_trees; i++)
            dist_RF[i] = new double [n_trees];
        for(int i = 0; i < n_trees; i++)
            for(int j = 0; j < n_trees; j++)
                dist_RF[i][j] = 0;
        vec_hashrf._hashtab.resize(vec_hashrf._hashtab2.size());

        for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti)
        {
           unsigned int sizeLinkedList = vec_hashrf._hashtab2[hti].size();
            if (sizeLinkedList > 0)
            {
                for (unsigned int i1 = 0; i1 < sizeLinkedList; ++i1)
                {
                    unsigned int bidi = vec_hashrf._hashtab2[hti][i1]._vec_treeidx.size();
                    for (unsigned int i2 = 0; i2 < bidi; ++i2)
                    {
                        BUCKET_STRUCT_T bk;
                        bk._hv2 = vec_hashrf._hashtab2[hti][i1]._hv2;
                        bk._t_i = vec_hashrf._hashtab2[hti][i1]._vec_treeidx[i2];
                        bk._dist = vec_hashrf._hashtab2[hti][i1]._vec_dist[i2];
                        vec_hashrf._hashtab[hti].push_back(bk);
                    }
                }
            }
        }
        vec_hashrf._hashtab2.clear();
        for (unsigned int hti = 0; hti < vec_hashrf._hashtab.size(); ++hti)
        {
           unsigned int sizeLinkedList = vec_hashrf._hashtab[hti].size();
            if (sizeLinkedList > 1)
            {
                vector<unsigned long> vec_hv2;
                vector<unsigned long>::iterator itr_vec_hv2;
                // Collect unique hv2 values in the linked list
                for (unsigned int i = 0; i < sizeLinkedList; ++i)
                {
                    unsigned long hv2 = vec_hashrf._hashtab[hti][i]._hv2;
                    if (vec_hv2.empty())
                        vec_hv2.push_back(hv2);
                    else
                    {
                        itr_vec_hv2 = find(vec_hv2.begin(), vec_hv2.end(), hv2);
                        if (itr_vec_hv2 == vec_hv2.end())
                            vec_hv2.push_back(hv2);
                    }
                }
                // distance
                vector<vector<float> > vec_dist(vec_hv2.size(), vector<float>(n_trees, 0));

                // Set the distance array with distance at proper tree index
                for (unsigned int i = 0; i < sizeLinkedList; ++i)
                {
                    for (unsigned int j = 0; j < vec_hv2.size(); ++j)
                    {
                        if (vec_hashrf._hashtab[hti][i]._hv2 == vec_hv2[j])
                            vec_dist[j][vec_hashrf._hashtab[hti][i]._t_i] = vec_hashrf._hashtab[hti][i]._dist;
                    }
                }
                // Update floatsim matrix using vec_dist
                for (unsigned int i = 0; i < vec_dist.size(); ++i)
                {
                    for (unsigned int j = 0; j < vec_dist[i].size(); ++j)
                    {
                        for (unsigned int k = 0; k < vec_dist[i].size(); ++k)
                        {
                            if (j == k)
                                continue;
                            else
                                dist_RF[j][k] += (vec_dist[i][j] - vec_dist[i][k] > 0) ? (vec_dist[i][j] - vec_dist[i][k]) : (vec_dist[i][k] - vec_dist[i][j]);
                        }
                    }
                }
                vec_hv2.clear();
                vec_dist.clear();
            }
            else if (sizeLinkedList == 1)
            {
                // propagate the dist value to other tree's distance
                for (unsigned int i = 0; i < n_trees; ++i)
                {
                    if (i == vec_hashrf._hashtab[hti][0]._t_i)
                        continue;
                    else
                    {
                        dist_RF[i][vec_hashrf._hashtab[hti][0]._t_i] += vec_hashrf._hashtab[hti][0]._dist;
                        dist_RF[vec_hashrf._hashtab[hti][0]._t_i][i] += vec_hashrf._hashtab[hti][0]._dist;
                    }
                }
            }
        }
        for (int i = 0; i < n_trees; i++)
        {
            for (int j = 0; j < i; j++)
            {
                dist_RF[i][j] = (dist_RF[i][j] + dist_RF[j][i]) / 4;
                dist_RF[j][i] = dist_RF[i][j];
            }
        }
    }

    return true;
}

//########################ZD comment########################################
//# The matching distance is given by the solution
//# to Hungarian algorithm on the table with entries
//# of number of XOR element in "bitstrofatree",
//# which is all possible bipartition of one tree.
//# Line 1415 may have a bug.
//########################ZD comment########################################

bool Trees::Compute_Matching_dist()
{
    delete_double_array(dist_match, n_trees);
    dist_match = new int *[n_trees];
    for (int i = 0; i < n_trees; i++)
        dist_match[i] = new int [n_trees];

    if (!bipartmatrixIsexisting())
    {
        cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
        return false;
    }
    int n_taxa = leaveslabelsmaps.size();
    int max_numberofbipartition = 0;
    int idx = 0;

    for (int i = 0; i < n_trees; i++)
    {
        dist_match[i][i] = 0;
        Array<Array<char> > bitstrofatree(numberofbipartition[i], Array<char> ());
        idx = 0;
        Get_bipartitionofonetree(treeset[i]->root, isrooted, 0, bitstrofatree, idx);

        for (int j = 0; j < i; j++)
        {
            Array<Array<char> > bitstrofatreej(numberofbipartition[j], Array<char> ());
            idx = 0;
            Get_bipartitionofonetree(treeset[j]->root, isrooted, 0, bitstrofatreej, idx);

            if (numberofbipartition[i] <= numberofbipartition[j])
                max_numberofbipartition = numberofbipartition[j];
            else
                max_numberofbipartition = numberofbipartition[i];
            int ** strdist = new int *[max_numberofbipartition];
            for (int k = 0; k < max_numberofbipartition; k++)
                strdist[k] = new int[max_numberofbipartition];
            hungarian_problem_t p;

            for (int k = 0; k < max_numberofbipartition; k++)
            {
                for (int l = 0; l < max_numberofbipartition; l++)
                {
                    int result = 0;
                    if (k < numberofbipartition[i] && l < numberofbipartition[j])
                    {
                        if(Array<char>::bitstrXOR(bitstrofatree[k], bitstrofatreej[l], result))
                            strdist[k][l] = result;
                    }
                    else if (k >= numberofbipartition[i] && l < numberofbipartition[j])
                    {
                        if(Array<char>::onebitstrXOR(bitstrofatreej[l], result))
                            strdist[k][l] = result;
                    }
                    else
                    {
                        if(Array<char>::onebitstrXOR(bitstrofatree[k], result))
                            strdist[k][l] = result;
                    }

                    if (((n_taxa % 2 == 0 && strdist[k][l] > ((int) n_taxa / 2))
                         || (n_taxa % 2 != 0 && strdist[k][l] >= ((int) n_taxa / 2) + 1))
                            && !isrooted)
                        strdist[k][l] = n_taxa - strdist[k][l];
                }
            }

            hungarian_init(&p, strdist , max_numberofbipartition, max_numberofbipartition, HUNGARIAN_MODE_MINIMIZE_COST) ;
            int mdist = hungarian_solve(&p);
            dist_match[i][j] = mdist;
            dist_match[j][i] = dist_match[i][j];
            hungarian_free(&p);
            for(int k = 0; k < max_numberofbipartition; k++)
            {
                delete [] strdist[k];
            }
            delete [] strdist;
        }
    }

    return true;
}

bool Trees::Compute_SPR_dist()
{
    delete_double_array(dist_SPR, n_trees);
    dist_SPR = new int *[n_trees];
    for (int i = 0; i < n_trees; i++)
        dist_SPR[i] = new int [n_trees];

    int NN = 10000;
    char *gzbuff = (char*) malloc(NN);

    char **taxa_str = new char*[leaveslabelsmaps.size()];
    for (int i = 0; i < leaveslabelsmaps.size(); i++)
    {
        string temp = leaveslabelsmaps.name(i);
        taxa_str[i] = new char[temp.length()+1];
        for (int j = 0; j < temp.length(); j++)
            taxa_str[i][j] = temp[j];
        taxa_str[i][temp.length()] = '\0';
    }

    std::map<std::string, int> label_map = std::map<std::string, int>();
    std::map<int, std::string> reverse_label_map = std::map<int, std::string>();

    //set random seed
    srand(unsigned(time(0)));
    std::vector<SPRNode *> sprtrees = std::vector<SPRNode *>();
    std::vector<std::string> names = std::vector<std::string>();
    for(int i = 0; i < n_trees; i++)
    {
        gzbuff[0] = NULL;
        TreeOPE::printTree_new(treeset[i]->root, taxa_str, gzbuff, 0, NN);
        std::string s = std::string(gzbuff);
        size_t loc = s.find_first_of("(");
        if (loc != string::npos)
        {
            std::string name = "";
            if(loc != 0)
            {
                name = s.substr(0, loc);
                s.erase(0, loc);
            }

            if (!isrooted)
                s = SPR_root(s);

            SPRNode *T2 = spr_building_tree(s);
            T2->labels_to_numbers(&label_map, &reverse_label_map);
            names.push_back(name);
            sprtrees.push_back(T2);
        }
    }

    for (int i = 0; i < leaveslabelsmaps.size(); i++)
    {
        delete[] taxa_str[i];
        taxa_str[i] = NULL;
    }
    delete[] taxa_str;
    taxa_str = NULL;
    free(gzbuff);

    for(int i = 0; i < n_trees; i++)
    {
        dist_SPR[i][i] = 0;
        for(int j = 0; j < i; j++)
        {
            dist_SPR[i][j] = rSPR_branch_and_bound_simple_clustering(sprtrees[i], sprtrees[j]);
            dist_SPR[j][i] = dist_SPR[i][j];
        }
    }

    // clean
    for (vector<SPRNode *>::iterator T2 = sprtrees.begin(); T2 != sprtrees.end(); T2++)
        (*T2)->delete_tree();

    return true;
}
#else

void Trees::Compute_Bipart_Covariance()
{
    SparseMatrix* trans = sbipartmatrix->transpose();
    SparseMatrix* result = sbipartmatrix->Multiply(*trans);
    double* M = sbipartmatrix->Mean(n_trees);
    double* Ones = new double[n_trees];
    for (int i = 0; i < n_trees; i++)
        Ones[i] = 1;
    double* tmp = sbipartmatrix->Multiply_vec(Ones);
    int Unique_idx = treecov_size;
    double ** temp1 = Vec_multiply(tmp, M, Unique_idx);
    double ** temp2 = Vec_multiply(M, M, Unique_idx);
    delete_double_array(treecov, treecov_size);
    treecov = new double*[Unique_idx];
    for (int i = 0; i < Unique_idx; i++)
        treecov[i] = new double[Unique_idx];

    for (int i = 0; i < Unique_idx; i++)
        for (int j = 0; j < Unique_idx; j++)
            treecov[i][j] = 0;

    /* Set up progress bar */
    int value = 0, maxProgressSize;
    maxProgressSize = (int) Unique_idx * Unique_idx;
    QProgressDialog dlg;
    dlg.setLabelText("Computing covariance matrix...");
    dlg.setWindowModality(Qt::WindowModal);
    dlg.setMinimum(0);
    dlg.setMaximum(maxProgressSize);
    dlg.setMinimumDuration(500);
    /*---------------------*/

    for (int i = 0; i < Unique_idx; i++)
    {
        for (int j = 0; j < Unique_idx; j++)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            treecov[i][j] = treecov[i][j] + (*result)(i,j) - temp1[i][j] - temp1[j][i] + n_trees * temp2[i][j];
            treecov[i][j] = treecov[i][j] / (n_trees - 1);

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButton;
            else
                value++;
            /*---------------------*/
        }
    }

    goto clearMemory;

cancelButton:
    /* Progress was canceled */

    // Delete any values that may have been saved in treecov
    delete_double_array(treecov, treecov_size);

    goto clearMemory;

clearMemory:
    delete [] M;
    M = NULL;

    delete [] Ones;
    Ones = NULL;

    delete [] tmp;
    tmp = NULL;

    for (int i = 0; i < Unique_idx; i++)
    {
            delete [] temp1[i];
            delete [] temp2[i];
    }
    delete [] temp1;
    temp1 = NULL;
    delete [] temp2;
    temp2 = NULL;

    if(trans != NULL)
    {
        delete trans;
        trans = NULL;
    }
    if(result != NULL)
    {
        delete result;
        result = NULL;
    }

    /* Close Progress bar */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/
}


bool Trees::Compute_RF_dist_by_hash(bool ISWEIGHTED)
{
    if(!isweighted && ISWEIGHTED)
    {
        cout << "Warning: The trees in memory are unweighted! Unable to compute weighted RF distances!\n\n";
        return false;
    }

    bool bUbid = false; // for counting the number pf unique bipartitions
    unsigned long uBID = 0;

    QProgressDialog dlg;    // Declare progress bar

    if (ISWEIGHTED == false)
    {
        /* Set up progress bar */
        int value = 0, maxProgressSize;
        maxProgressSize = (int) vec_hashrf._hashtab2.size() + n_trees*((n_trees + 1) / 2);
        dlg.setLabelText("Computing Unweighted RF distance...");
        dlg.setWindowModality(Qt::WindowModal);
        dlg.setMinimum(0);
        dlg.setMaximum(maxProgressSize);
        dlg.setMinimumDuration(500);
        /*---------------------*/


        delete_double_array(dist_URF, n_trees);
        dist_URF = new double *[n_trees];
        for (int i = 0; i < n_trees; i++)
            dist_URF[i] = new double [n_trees];
        for(int i = 0; i < n_trees; i++)
        {
            for(int j = 0; j < n_trees; j++)
                dist_URF[i][j] = 0;
        }
        for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            unsigned int sizeVec = vec_hashrf._hashtab2[hti].size();
            if (sizeVec)
            {
                uBID += sizeVec;
                if (!bUbid)
                {
                    for (unsigned int i = 0; i < sizeVec; ++i)
                    {
                        unsigned int sizeTreeIdx = vec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
                        if (sizeTreeIdx > 1)
                        {
                            for (unsigned int j = 0; j < sizeTreeIdx; ++j)
                            {
                                for (unsigned int k = 0; k < sizeTreeIdx; ++k)
                                {
                                    if (j == k)
                                        continue;
                                    else
                                        dist_URF[vec_hashrf._hashtab2[hti][i]._vec_treeidx[j]][vec_hashrf._hashtab2[hti][i]._vec_treeidx[k]] += 1;
                                }
                            }
                        }
                    }
                }
            }

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButtonURFsect1;
            else
                value++;
            /*---------------------*/
        }
        for (int i = 0; i < n_trees; i++)
        {
            for (int j = 0; j < i; j++)
            {
                /* Update progress bar */
                dlg.setValue(value);
                /*---------------------*/

                dist_URF[i][j] = (double) ((numberofbipartition[i] + numberofbipartition[j] - 2 * dist_URF[i][j]) / 2);
                dist_URF[j][i] = dist_URF[i][j];

                /* Update progress bar */
                value++;
                /*---------------------*/
            }
        }

        /* Close Progress bar */
        dlg.setValue(maxProgressSize);
        dlg.reset();
        dlg.close();
        /*--------------------*/
    }
    else if (ISWEIGHTED == true)
    {
        delete_double_array(dist_RF, n_trees);
        dist_RF = new double *[n_trees];
        for (int i = 0; i < n_trees; i++)
            dist_RF[i] = new double [n_trees];
        for(int i = 0; i < n_trees; i++)
            for(int j = 0; j < n_trees; j++)
                dist_RF[i][j] = 0;
        vec_hashrf._hashtab.resize(vec_hashrf._hashtab2.size());

        /* Set up progress bar */
        int value = 0, maxProgressSize;
        maxProgressSize = (int) vec_hashrf._hashtab2.size() + vec_hashrf._hashtab.size() + n_trees*((n_trees + 1) / 2);
        dlg.setLabelText("Computing Weighted RF distance...");
        dlg.setWindowModality(Qt::WindowModal);
        dlg.setMinimum(0);
        dlg.setMaximum(maxProgressSize);
        dlg.setMinimumDuration(500);
        /*---------------------*/

        /* Section 1 */
        for (unsigned int hti = 0; hti < vec_hashrf._hashtab2.size(); ++hti)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            unsigned int sizeLinkedList = vec_hashrf._hashtab2[hti].size();
            if (sizeLinkedList > 0)
            {
                for (unsigned int i1 = 0; i1 < sizeLinkedList; ++i1)
                {
                    unsigned int bidi = vec_hashrf._hashtab2[hti][i1]._vec_treeidx.size();
                    for (unsigned int i2 = 0; i2 < bidi; ++i2)
                    {
                        BUCKET_STRUCT_T bk;
                        bk._hv2 = vec_hashrf._hashtab2[hti][i1]._hv2;
                        bk._t_i = vec_hashrf._hashtab2[hti][i1]._vec_treeidx[i2];
                        bk._dist = vec_hashrf._hashtab2[hti][i1]._vec_dist[i2];
                        vec_hashrf._hashtab[hti].push_back(bk);
                    }
                }
            }

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButtonRFsect1;
            else
                value++;
            /*---------------------*/
        }

        /* Section 2 */
        vec_hashrf._hashtab2.clear();
        for (unsigned int hti = 0; hti < vec_hashrf._hashtab.size(); ++hti)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            unsigned int sizeLinkedList = vec_hashrf._hashtab[hti].size();
            if (sizeLinkedList > 1)
            {
                vector<unsigned long> vec_hv2;
                vector<unsigned long>::iterator itr_vec_hv2;
                // Collect unique hv2 values in the linked list
                for (unsigned int i = 0; i < sizeLinkedList; ++i)
                {
                    unsigned long hv2 = vec_hashrf._hashtab[hti][i]._hv2;
                    if (vec_hv2.empty())
                        vec_hv2.push_back(hv2);
                    else
                    {
                        itr_vec_hv2 = find(vec_hv2.begin(), vec_hv2.end(), hv2);
                        if (itr_vec_hv2 == vec_hv2.end())
                            vec_hv2.push_back(hv2);
                    }
                }
                // distance
                vector<vector<float> > vec_dist(vec_hv2.size(), vector<float>(n_trees, 0));

                // Set the distance array with distance at proper tree index
                for (unsigned int i = 0; i < sizeLinkedList; ++i)
                {
                    for (unsigned int j = 0; j < vec_hv2.size(); ++j)
                    {
                        if (vec_hashrf._hashtab[hti][i]._hv2 == vec_hv2[j])
                            vec_dist[j][vec_hashrf._hashtab[hti][i]._t_i] = vec_hashrf._hashtab[hti][i]._dist;
                    }
                }
                // Update floatsim matrix using vec_dist
                for (unsigned int i = 0; i < vec_dist.size(); ++i)
                {
                    for (unsigned int j = 0; j < vec_dist[i].size(); ++j)
                    {
                        for (unsigned int k = 0; k < vec_dist[i].size(); ++k)
                        {
                            if (j == k)
                                continue;
                            else
                                dist_RF[j][k] += (vec_dist[i][j] - vec_dist[i][k] > 0) ? (vec_dist[i][j] - vec_dist[i][k]) : (vec_dist[i][k] - vec_dist[i][j]);
                        }
                    }
                }
                vec_hv2.clear();
                vec_dist.clear();
            }
            else if (sizeLinkedList == 1)
            {
                // propagate the dist value to other tree's distance
                for (unsigned int i = 0; i < n_trees; ++i)
                {
                    if (i == vec_hashrf._hashtab[hti][0]._t_i)
                        continue;
                    else
                    {
                        dist_RF[i][vec_hashrf._hashtab[hti][0]._t_i] += vec_hashrf._hashtab[hti][0]._dist;
                        dist_RF[vec_hashrf._hashtab[hti][0]._t_i][i] += vec_hashrf._hashtab[hti][0]._dist;
                    }
                }
            }

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButtonRFsect2;
            else
                value++;
            /*---------------------*/
        }

        /* Section 3 */
        for (int i = 0; i < n_trees; i++)
        {
            for (int j = 0; j < i; j++)
            {
                /* Update progress bar */
                dlg.setValue(value);
                /*---------------------*/

                dist_RF[i][j] = (dist_RF[i][j] + dist_RF[j][i]) / 4;
                dist_RF[j][i] = dist_RF[i][j];

                /* Update progress bar */
                if(dlg.wasCanceled())
                    goto cancelButtonRFsect3;
                else
                    value++;
                /*---------------------*/
            }
        }

        /* Close Progress bar */
        dlg.setValue(maxProgressSize);
        dlg.reset();
        dlg.close();
        /*--------------------*/
    }

    return true;
cancelButtonURFsect1:
    /* Progress was canceled for Unweighted RF*/

    // Delete any values in dist_URF
    delete_double_array(dist_URF, n_trees);

    /* Close Progress bar */
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;

cancelButtonRFsect1:
    /* Progress was canceled in section 1 for Weighted RF*/

    vec_hashrf._hashtab2.clear();

    // Delete any values in dist_URF
    delete_double_array(dist_RF, n_trees);
    vec_hashrf._hashtab.clear();

    /* Close Progress bar */
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;

cancelButtonRFsect2:
    /* Progress was canceled in section 2 for Weighted RF*/

    // Delete any values in dist_URF
    delete_double_array(dist_RF, n_trees);
    vec_hashrf._hashtab.clear();

    /* Close Progress bar */
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;

cancelButtonRFsect3:
    /* Progress was canceled in section 3 for Weighted RF*/

    // Delete any values in dist_URF
    delete_double_array(dist_RF, n_trees);
    vec_hashrf._hashtab.clear();

    /* Close Progress bar */
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;

}

bool Trees::Compute_Matching_dist()
{
    delete_double_array(dist_match, n_trees);
    dist_match = new int *[n_trees];
    for (int i = 0; i < n_trees; i++)
        dist_match[i] = new int [n_trees];

    if (!bipartmatrixIsexisting())
    {
        cout << "Warning: There is no bipartition matrix in the memory! Please compute it first.\n\n";
        return false;
    }
    int n_taxa = leaveslabelsmaps.size();
    int max_numberofbipartition = 0;
    int idx = 0;

    /* Set up progress bar */
    int value = 0, maxProgressSize;
    maxProgressSize = (int) n_trees + n_trees*((n_trees + 1) / 2);
    QProgressDialog dlg;
    dlg.setLabelText("Computing Matching distance...");
    dlg.setWindowModality(Qt::WindowModal);
    dlg.setMinimum(0);
    dlg.setMaximum(maxProgressSize);
    dlg.setMinimumDuration(500);
    /*---------------------*/

    for (int i = 0; i < n_trees; i++)
    {
        dist_match[i][i] = 0;
        Array<Array<char> > bitstrofatree(numberofbipartition[i], Array<char> ());
        idx = 0;
        Get_bipartitionofonetree(treeset[i]->root, isrooted, 0, bitstrofatree, idx);

        for (int j = 0; j < i; j++)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            Array<Array<char> > bitstrofatreej(numberofbipartition[j], Array<char> ());
            idx = 0;
            Get_bipartitionofonetree(treeset[j]->root, isrooted, 0, bitstrofatreej, idx);

            if (numberofbipartition[i] <= numberofbipartition[j])
                max_numberofbipartition = numberofbipartition[j];
            else
                max_numberofbipartition = numberofbipartition[i];
            int ** strdist = new int *[max_numberofbipartition];
            for (int k = 0; k < max_numberofbipartition; k++)
                strdist[k] = new int[max_numberofbipartition];
            hungarian_problem_t p;

            for (int k = 0; k < max_numberofbipartition; k++)
            {
                for (int l = 0; l < max_numberofbipartition; l++)
                {
                    int result = 0;
                    if (k < numberofbipartition[i] && l < numberofbipartition[j])
                    {
                        if(Array<char>::bitstrXOR(bitstrofatree[k], bitstrofatreej[l], result))
                            strdist[k][l] = result;
                    }
                    else if (k >= numberofbipartition[i] && l < numberofbipartition[j])
                    {
                        if(Array<char>::onebitstrXOR(bitstrofatreej[l], result))
                            strdist[k][l] = result;
                    }
                    else
                    {
                        if(Array<char>::onebitstrXOR(bitstrofatree[k], result))
                            strdist[k][l] = result;
                    }

                    if (((n_taxa % 2 == 0 && strdist[k][l] > ((int) n_taxa / 2))
                         || (n_taxa % 2 != 0 && strdist[k][l] >= ((int) n_taxa / 2) + 1))
                            && !isrooted)
                        strdist[k][l] = n_taxa - strdist[k][l];
                }
            }

            hungarian_init(&p, strdist , max_numberofbipartition, max_numberofbipartition, HUNGARIAN_MODE_MINIMIZE_COST) ;
            int mdist = hungarian_solve(&p);
            dist_match[i][j] = mdist;
            dist_match[j][i] = dist_match[i][j];
            hungarian_free(&p);
            for(int k = 0; k < max_numberofbipartition; k++)
            {
                delete [] strdist[k];
            }
            delete [] strdist;

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButton;
            else
                value++;
            /*---------------------*/
        }
    }

    /* Close Progress bar */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return true;

cancelButton:
    /* Progress was canceled */

    // Delete any values that may have been saved in dist_match
    delete_double_array(dist_match, n_trees);

    /* Close Progress bar */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;
}

bool Trees::Compute_SPR_dist()
{
    delete_double_array(dist_SPR, n_trees);
    dist_SPR = new int *[n_trees];
    for (int i = 0; i < n_trees; i++)
        dist_SPR[i] = new int [n_trees];

    int NN = 10000;
    char *gzbuff = (char*) malloc(NN);

    /* Set up progress bar */
    int value = 0, maxProgressSize;
    maxProgressSize = (int) n_trees + n_trees*((n_trees + 1) / 2);
    QProgressDialog dlg;
    dlg.setLabelText("Computing SPR distance...");
    dlg.setWindowModality(Qt::WindowModal);
    dlg.setMinimum(0);
    dlg.setMaximum(maxProgressSize);
    dlg.setMinimumDuration(500);
    /*---------------------*/

    char **taxa_str = new char*[leaveslabelsmaps.size()];
    for (int i = 0; i < leaveslabelsmaps.size(); i++)
    {
        string temp = leaveslabelsmaps.name(i);
        taxa_str[i] = new char[temp.length()+1];
        for (int j = 0; j < temp.length(); j++)
            taxa_str[i][j] = temp[j];
        taxa_str[i][temp.length()] = '\0';
    }

    std::map<std::string, int> label_map = std::map<std::string, int>();
    std::map<int, std::string> reverse_label_map = std::map<int, std::string>();

    /* Section 1----------------------------------------------------------*/
    //set random seed
    srand(unsigned(time(0)));
    std::vector<SPRNode *> sprtrees = std::vector<SPRNode *>();
    std::vector<std::string> names = std::vector<std::string>();
    for(int i = 0; i < n_trees; i++)
    {
        /* Update progress bar */
        dlg.setValue(value);
        /*---------------------*/

        gzbuff[0] = NULL;
        TreeOPE::printTree_new(treeset[i]->root, taxa_str, gzbuff, 0, NN);
        std::string s = std::string(gzbuff);
        size_t loc = s.find_first_of("(");
        if (loc != string::npos)
        {
            std::string name = "";
            if(loc != 0)
            {
                name = s.substr(0, loc);
                s.erase(0, loc);
            }

            if (!isrooted)
                s = SPR_root(s);

            SPRNode *T2 = spr_building_tree(s);
            T2->labels_to_numbers(&label_map, &reverse_label_map);
            names.push_back(name);
            sprtrees.push_back(T2);
        }

        /* Update progress bar */
        if(dlg.wasCanceled())
            goto cancelButtonSect1;
        else
            value++;
        /*---------------------*/
    }

    for (int i = 0; i < leaveslabelsmaps.size(); i++)
    {
        delete[] taxa_str[i];
        taxa_str[i] = NULL;
    }
    delete[] taxa_str;
    taxa_str = NULL;
    free(gzbuff);

    /* Section 2: Compute SPR----------------------------------------------------*/
    for(int i = 0; i < n_trees; i++)
    {
        dist_SPR[i][i] = 0;
        for(int j = 0; j < i; j++)
        {
            /* Update progress bar */
            dlg.setValue(value);
            /*---------------------*/

            dist_SPR[i][j] = rSPR_branch_and_bound_simple_clustering(sprtrees[i], sprtrees[j]);
            dist_SPR[j][i] = dist_SPR[i][j];

            /* Update progress bar */
            if(dlg.wasCanceled())
                goto cancelButtonSect2;
            else
                value++;
            /*---------------------*/
        }
    }


    // clean
    for (vector<SPRNode *>::iterator T2 = sprtrees.begin(); T2 != sprtrees.end(); T2++)
        (*T2)->delete_tree();


    /* Close Progress Box */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return true;

cancelButtonSect1:
     /* If cancelled in section 1, release data and return false */
    for (int i = 0; i < leaveslabelsmaps.size(); i++)
    {
        delete[] taxa_str[i];
        taxa_str[i] = NULL;
    }
    delete[] taxa_str;
    taxa_str = NULL;
    free(gzbuff);


    // clean
    for (vector<SPRNode *>::iterator T2 = sprtrees.begin(); T2 != sprtrees.end(); T2++)
        (*T2)->delete_tree();


    /* Close Progress Box */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;

cancelButtonSect2:
     /* If cancelled in section 2 (when computing SPR), release data and return false */

    // Delete any values that may have been saved in dist_SPR
    delete_double_array(dist_SPR, n_trees);

    // clean
    for (vector<SPRNode *>::iterator T2 = sprtrees.begin(); T2 != sprtrees.end(); T2++)
        (*T2)->delete_tree();


    /* Close Progress Box */
    dlg.setValue(maxProgressSize);
    dlg.reset();
    dlg.close();
    /*--------------------*/

    return false;
}

#endif

#endif
