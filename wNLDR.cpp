
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

// wNLDR.cpp
//  cost function: Linear MDS, Metric MDS, CCA, SAMMON
//  algorithm    : Gauss Seidel, Majorization, stachatic gradient.
//          March/11/2010
//                    by whuang

#ifndef WNLDR_CPP
#define WNLDR_CPP

#include "wNLDR.h"

void NLDR::init_NLDR(String fname, double **dist, int sizeinput, String dim, String cost, String algo, String init_md, String flag, long seed, String para_fname)
{
	// load initial info.
//	NLDR_init_parameters(para_fname);  this is for command line version only.
	time_cost = -1;
	dim_str = dim;
	cost_function = cost;
	algorithm = algo;
	file_flag = flag;
	STRESS = -1;

    File D_file(fname);

    D_prefname = fname;//D_file.prefix_name_lastof();
    D_postfname = fname;//D_file.postfix_name_lastof();

    int dim_int = atoi(dim_str);
    if(dim_int < 1)
    {
        std::cout << "Error: The dimension must be equal to or greater than 1! Please use command -h to see help." << std::endl;
        exit(0);
    }

    if(dist == NULL && sizeinput == 0)
    {
        if(!D_file.is_open())
        {
            std::cout << "Error: File \"" << fname << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
            exit(0);
        }
        size = D_file.lines();
        this->NLDR_load_D();
    }
    else if(dist != NULL && sizeinput > 0)
    {
/*        if(D.get_row() != 0)
        {
            if(D.get_col() != 0)
                for(int i = 0; i < D.get_row(); i++)
                    delete [] D.matrix[i];
            delete [] D.matrix;
        }
        D.matrix = dist;*/
        size = sizeinput;
        D.resize(size, size);
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                D.matrix[i][j] = (double) dist[i][j];
    } else
    {
        std::cout << "Error: Incorrect input data parameters" << std::endl;
        exit(0);
    }

	if(seed < -1)
	{
        std::cout << "Error: The seed of random generator must be equal to or greater than -1 and also be a integer! Please use command -h to see help." << std::endl;
		exit(0);
    }
    this->NLDR_init_X(init_md, seed);
}

#ifdef COMMAND_LINE_VERSION

void NLDR::init_NLDR(String fname, String ftype, String dim, String cost, String algo, String init_md, String flag, long seed, String para_fname)
{
	// load initial info.
    NLDR_init_parameters(para_fname);
	time_cost = -1;
	dim_str = dim;
	cost_function = cost;
	algorithm = algo;
	file_flag = flag;
    STRESS = -1;
    File D_file(fname);

	if(!D_file.is_open())
	{
        std::cout << "Error: File \"" << fname << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
		exit(0);
    }

    D_prefname = D_file.prefix_name_lastof();
    D_postfname = D_file.postfix_name_lastof();
	size = D_file.lines();
	D_file.seek(0);
	int dim_int = atoi(dim_str);

	if(dim_int < 1)
	{
        std::cout << "Error: The dimension must be equal to or greater than 1! Please use command -h to see help." << std::endl;
		exit(0);
	}

	if(ftype == (String) "DIS")
	{
		this->NLDR_load_D();
		this->NLDR_init_X(init_md, seed);
	} else
	if(ftype == (String) "COR")
	{
		this->NLDR_load_MX();
		this->NLDR_init_X(init_md, seed);
	} else
        std::cout << "Warning: Undefined matrix type" << std::endl;

	if(seed < -1)
	{
        std::cout << "Error: The seed of random generator must be equal to or greater than -1 and also be a integer! Please use command -h to see help." << std::endl;
		exit(0);
	}
}

void NLDR::NLDR_init_parameters(String para_filename)
{
	parameters.distance_file_type = 0;
	parameters.length_tru = 11;
    parameters.interval_tru = 5;
    parameters.length_con = 11;
	parameters.iterval_con = 5;
	parameters.random_start = -100;
	parameters.random_end = 100;

	parameters.KRU_LIN_e = 0.00001;
	parameters.KRU_LIN_max_iter = 1000;
	parameters.KRU_LIN_run_time = 0;
	parameters.KRU_MAJ_e = 0.00001;
	parameters.KRU_MAJ_max_iter = 1000;
	parameters.KRU_MAJ_run_time = 0;
	parameters.KRU_GAU_e = 0.0001;
	parameters.KRU_GAU_max_iter = 1000;
	parameters.KRU_GAU_run_time = 0;
	parameters.KRU_STO_epochs = 100;
	parameters.KRU_STO_alpha0 = 0.5;
	parameters.KRU_STO_alphan = 0.01;

	parameters.NOR_MAJ_e = 0.00001;
	parameters.NOR_MAJ_max_iter = 1000;
	parameters.NOR_MAJ_run_time = 0;
	parameters.NOR_GAU_e = 0.0001;
	parameters.NOR_GAU_max_iter = 1000;
	parameters.NOR_GAU_run_time = 0;
	parameters.NOR_STO_epochs = 100;
	parameters.NOR_STO_alpha0 = 0.25;
	parameters.NOR_STO_alphan = 0.01;

	parameters.NLM_MAJ_e = 0.0001;
	parameters.NLM_MAJ_max_iter = 1000;
	parameters.NLM_MAJ_run_time = 0;
	parameters.NLM_GAU_e = 0.0001;
	parameters.NLM_GAU_max_iter = 1000;
	parameters.NLM_GAU_run_time = 0;
	parameters.NLM_STO_epochs = 100;
	parameters.NLM_STO_alpha0 = 0.25;
	parameters.NLM_STO_alphan = 0.01;

	parameters.CCA_MAJ_e = 0.0001;
	parameters.CCA_MAJ_max_iter = 1000;
	parameters.CCA_MAJ_run_time = 0;
	parameters.CCA_MAJ_lambda0 = 12000;
	parameters.CCA_MAJ_lambdan = 1;
	parameters.CCA_GAU_e = 0.0001;
	parameters.CCA_GAU_max_iter = 1000;
	parameters.CCA_GAU_run_time = 0;
	parameters.CCA_GAU_lambda0 = 12000;
	parameters.CCA_GAU_lambdan = 1;
	parameters.CCA_STO_epochs = 100;
	parameters.CCA_STO_lambda0 = 12000;
	parameters.CCA_STO_lambdan = 1;
	parameters.CCA_STO_alpha0 = 0.5;
	parameters.CCA_STO_alphan = 0.01;

	if(para_filename == (String) "")
        return;
	Array<Mapping<Mix, Mix> > para_arr;
	Mapping<Mix, Mix> para;

	File fpara(para_filename);
	if(!fpara.is_open())
	{
        std::cout << "Error: File \"" << para_filename << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
		exit(0);
	}

    para_arr = LPIMPORT_FORM->Read(para_filename);
	int num = para_arr.get_length();

	Mix Key, Value;
	for(int i = 0; i < num; i++)
    {
        Key = para_arr[i][(String) "name"];
        Value = para_arr[i][(String) "value"];
        para[Key] = Value;
	}

	if(para.Include_Key((String) "length_tru"))
		parameters.length_tru = para[(String) "length_tru"];
	if(para.Include_Key((String) "interval_tru"))
		parameters.interval_tru = para[(String) "interval_tru"];
	if(para.Include_Key((String) "length_con"))
		parameters.length_con = para[(String) "length_con"];
	if(para.Include_Key((String) "iterval_con"))
		parameters.iterval_con = para[(String) "iterval_con"];
	if(para.Include_Key((String) "random_start"))
		parameters.random_start = para[(String) "random_start"];
	if(para.Include_Key((String) "random_end"))
		parameters.random_end = para[(String) "random_end"];

	if(para.Include_Key((String) "KRU_LIN_e"))
		parameters.KRU_LIN_e = para[(String) "KRU_LIN_e"];
	if(para.Include_Key((String) "KRU_LIN_max_iter"))
		parameters.KRU_LIN_max_iter = para[(String) "KRU_LIN_max_iter"];
	if(para.Include_Key((String) "KRU_LIN_run_time"))
		parameters.KRU_LIN_run_time = (long) ((int) para[(String) "KRU_LIN_run_time"]);
	if(para.Include_Key((String) "KRU_MAJ_e"))
		parameters.KRU_MAJ_e = para[(String) "KRU_MAJ_e"];
	if(para.Include_Key((String) "KRU_MAJ_max_iter"))
		parameters.KRU_MAJ_max_iter = para[(String) "KRU_MAJ_max_iter"];
	if(para.Include_Key((String) "KRU_MAJ_run_time"))
		parameters.KRU_MAJ_run_time = (long) ((int) para[(String) "KRU_MAJ_run_time"]);
	if(para.Include_Key((String) "KRU_GAU_e"))
		parameters.KRU_GAU_e = para[(String) "KRU_GAU_e"];
	if(para.Include_Key((String) "KRU_GAU_max_iter"))
		parameters.KRU_GAU_max_iter = para[(String) "KRU_GAU_max_iter"];
	if(para.Include_Key((String) "KRU_GAU_run_time"))
		parameters.KRU_GAU_run_time = (long) ((int)para[(String) "KRU_GAU_run_time"]);
	if(para.Include_Key((String) "KRU_STO_epochs"))
		parameters.KRU_STO_epochs = para[(String) "KRU_STO_epochs"];
	if(para.Include_Key((String) "KRU_STO_alpha0"))
		parameters.KRU_STO_alpha0 = para[(String) "KRU_STO_alpha0"];
	if(para.Include_Key((String) "KRU_STO_alphan"))
		parameters.KRU_STO_alphan = para[(String) "KRU_STO_alphan"];

	if(para.Include_Key((String) "NOR_MAJ_e"))
		parameters.NOR_MAJ_e = para[(String) "NOR_MAJ_e"];
	if(para.Include_Key((String) "NOR_MAJ_max_iter"))
		parameters.NOR_MAJ_max_iter = para[(String) "NOR_MAJ_max_iter"];
	if(para.Include_Key((String) "NOR_MAJ_run_time"))
		parameters.NOR_MAJ_run_time = (long) ((int) para[(String) "NOR_MAJ_run_time"]);
	if(para.Include_Key((String) "NOR_GAU_e"))
		parameters.NOR_GAU_e = para[(String) "NOR_GAU_e"];
	if(para.Include_Key((String) "NOR_GAU_max_iter"))
		parameters.NOR_GAU_max_iter = para[(String) "NOR_GAU_max_iter"];
	if(para.Include_Key((String) "NOR_GAU_run_time"))
		parameters.NOR_GAU_run_time = (long) ((int) para[(String) "NOR_GAU_run_time"]);
	if(para.Include_Key((String) "NOR_STO_epochs"))
		parameters.NOR_STO_epochs = para[(String) "NOR_STO_epochs"];
	if(para.Include_Key((String) "NOR_STO_alpha0"))
		parameters.NOR_STO_alpha0 = para[(String) "NOR_STO_alpha0"];
	if(para.Include_Key((String) "NOR_STO_alphan"))
		parameters.NOR_STO_alphan = para[(String) "NOR_STO_alphan"];

	if(para.Include_Key((String) "NLM_MAJ_e"))
		parameters.NLM_MAJ_e = para[(String) "NLM_MAJ_e"];
	if(para.Include_Key((String) "NLM_MAJ_max_iter"))
		parameters.NLM_MAJ_max_iter = para[(String) "NLM_MAJ_max_iter"];
	if(para.Include_Key((String) "NLM_MAJ_run_time"))
		parameters.NLM_MAJ_run_time = (long) ((int) para[(String) "NLM_MAJ_run_time"]);
	if(para.Include_Key((String) "NLM_GAU_e"))
		parameters.NLM_GAU_e = para[(String) "NLM_GAU_e"];
	if(para.Include_Key((String) "NLM_GAU_max_iter"))
		parameters.NLM_GAU_max_iter = para[(String) "NLM_GAU_max_iter"];
	if(para.Include_Key((String) "NLM_GAU_run_time"))
		parameters.NLM_GAU_run_time = (long) ((int) para[(String) "NLM_GAU_run_time"]);
	if(para.Include_Key((String) "NLM_STO_epochs"))
		parameters.NLM_STO_epochs = para[(String) "NLM_STO_epochs"];
	if(para.Include_Key((String) "NLM_STO_alpha0"))
		parameters.NLM_STO_alpha0 = para[(String) "NLM_STO_alpha0"];
	if(para.Include_Key((String) "NLM_STO_alphan"))
		parameters.NLM_STO_alphan = para[(String) "NLM_STO_alphan"];

	if(para.Include_Key((String) "CCA_MAJ_e"))
		parameters.CCA_MAJ_e = para[(String) "CCA_MAJ_e"];
	if(para.Include_Key((String) "CCA_MAJ_max_iter"))
		parameters.CCA_MAJ_max_iter = para[(String) "CCA_MAJ_max_iter"];
	if(para.Include_Key((String) "CCA_MAJ_run_time"))
		parameters.CCA_MAJ_run_time = (long) ((int) para[(String) "CCA_MAJ_run_time"]);
	if(para.Include_Key((String) "CCA_MAJ_lambda0"))
		parameters.CCA_MAJ_lambda0 = para[(String) "CCA_MAJ_lambda0"];
	if(para.Include_Key((String) "CCA_MAJ_lambdan"))
		parameters.CCA_MAJ_lambdan = para[(String) "CCA_MAJ_lambdan"];
	if(para.Include_Key((String) "CCA_GAU_e"))
		parameters.CCA_GAU_e = para[(String) "CCA_GAU_e"];
	if(para.Include_Key((String) "CCA_GAU_max_iter"))
		parameters.CCA_GAU_max_iter = para[(String) "CCA_GAU_max_iter"];
	if(para.Include_Key((String) "CCA_GAU_run_time"))
		parameters.CCA_GAU_run_time = (long) ((int) para[(String) "CCA_GAU_run_time"]);
	if(para.Include_Key((String) "CCA_GAU_lambda0"))
		parameters.CCA_GAU_lambda0 = para[(String) "CCA_GAU_lambda0"];
	if(para.Include_Key((String) "CCA_GAU_lambdan"))
		parameters.CCA_GAU_lambdan = para[(String) "CCA_GAU_lambdan"];
	if(para.Include_Key((String) "CCA_STO_epochs"))
		parameters.CCA_STO_epochs = para[(String) "CCA_STO_epochs"];
	if(para.Include_Key((String) "CCA_STO_lambda0"))
		parameters.CCA_STO_lambda0 = para[(String) "CCA_STO_lambda0"];
	if(para.Include_Key((String) "CCA_STO_lambdan"))
		parameters.CCA_STO_lambdan = para[(String) "CCA_STO_lambdan"];
	if(para.Include_Key((String) "CCA_STO_alpha0"))
		parameters.CCA_STO_alpha0 = para[(String) "CCA_STO_alpha0"];
	if(para.Include_Key((String) "CCA_STO_alphan"))
        parameters.CCA_STO_alphan = para[(String) "CCA_STO_alphan"];
}

void NLDR::NLDR_load_MX()
{
	String fname = D_prefname;
	fname += ".";
	fname += D_postfname;
	File D_file(fname);
	if(!D_file.is_open())
	{
        std::cout << "Error: File \"" << fname << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
        exit(0);
    }
	int dim = D_file.cols();
	D_file.seek(0);

    std::cout << "size : " << size << ", dim:" << dim << std::endl;//-------
	MX.resize(size, dim);
	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
			D_file >> MX.matrix[i][j];

	D.resize(size, size);
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			for(int k = 0; k < dim; k++)
				D.matrix[i][j] += (MX.matrix[i][k] - MX.matrix[j][k]) * (MX.matrix[i][k] - MX.matrix[j][k]);
			D.matrix[i][j] = sqrt(D.matrix[i][j]);
			D.matrix[j][i] = D.matrix[i][j];
		}

	File out_D("test_D.out");///-----
	out_D.clean();
	out_D.seek(0);
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			out_D << D.matrix[i][j] << "\t";
		}
		out_D << "\n";//---
	}
	out_D.close();///------------
};
#endif

void NLDR::Compute_NLDR()
{
    std::cout << "Start program: " << D_prefname << "_" << dim_str << "D_" << cost_function << "_" << algorithm << std::endl;
	long start = clock();
	long end;
    std::cout << "Compute start time:" << start << std::endl;
	if(cost_function == (String) "CLASSIC_MDS")
		this->CLASSIC_MDS();
	else
	if(cost_function == (String) "KRUSKAL1")
		this->KRUSKAL1();
	else
	if(cost_function == (String) "NORMALIZED")
		this->NORMALIZED();
	else
	if(cost_function == (String) "SAMMON")
		this->SAMMON();
	else
	if(cost_function == (String) "CCA")
		this->CCA();
	end = clock();
	time_cost = end - start;
    std::cout << "Compute end time:" << end << std::endl;
    std::cout << "Compute time cost:" << time_cost << std::endl;
}

void NLDR::result_analysis()
{
	long start = clock();
	long end;
    std::cout << "Analysis start time:" << start << std::endl;
	Trustworthiness.resize(parameters.length_tru, 2);
	Continuity.resize(parameters.length_con, 2);
	for(int i = 0; i < parameters.length_tru; i++)
		Trustworthiness(i, 0) = (i * parameters.interval_tru > 0) ? i * parameters.interval_tru : 1;

	for(int i = 0; i < parameters.length_con; i++)
		Continuity(i, 0) = (i * parameters.iterval_con > 0) ? i * parameters.iterval_con : 1;

	Trustworthiness_analysis();

	Continuity_analysis();

	oneNN_analysis();
	end = clock();
    std::cout << "Analysis end time:" << end << std::endl;
    std::cout << "Analysis time cost:" << end - start << std::endl;
};

void NLDR::output_to_files()
{
	String filename_COR, filename_DIS, filename_STR, filename_TIM, filename_TRU, filename_CON, filename_1NN;
	this->make_output_file_names(filename_COR, filename_DIS, filename_STR, filename_TIM, filename_TRU, filename_CON, filename_1NN);
	File file_COR(filename_COR);
	File file_DIS(filename_DIS);
	File file_STR(filename_STR);
	File file_TIM(filename_TIM);
	File file_TRU(filename_TRU);
	File file_CON(filename_CON);
	File file_1NN(filename_1NN);
	if(COR.get_row() > 0 && COR.get_col() > 0)
	{
		file_COR.clean();
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < atoi(dim_str); j++)
				file_COR << COR(i, j) << "\t";
            file_COR << "" << std::endl;
		}
	} else
        std::cout << "Warning: Have not computed COR yet" << std::endl;

	if(DIS.get_row() > 0 && DIS.get_col() > 0)
	{
		file_DIS.clean();
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j <= i; j++)
				file_DIS << DIS(i, j) << "\t";
            file_DIS << "" << std::endl;
		}
	} else
        std::cout << "Warning: Have not computed DIS yet" << std::endl;

	if(STRESS != -1)
	{
		file_STR.clean();
        file_STR << STRESS << std::endl;
	} else
        std::cout << "Warning: Have not computed STRESS yet" << std::endl;

	if(time_cost != -1)
	{
		file_TIM.clean();
        file_TIM << time_cost << std::endl;
	} else
        std::cout << "Warning: Have not recorded time cost yet" << std::endl;

	if(Trustworthiness.get_row() > 0 && Trustworthiness.get_col() > 0)
	{
		file_TRU.clean();
		for(int i = 0; i < Trustworthiness.get_row(); i++)
		{
			for(int j = 0; j < Trustworthiness.get_col(); j++)
				file_TRU << Trustworthiness(i, j) << "\t";
            file_TRU << "" << std::endl;
		}
	}

	if(Continuity.get_row() > 0 && Continuity.get_col() > 0)
	{
		file_CON.clean();
		for(int i = 0; i < Continuity.get_row(); i++)
		{
			for(int j = 0; j < Continuity.get_col(); j++)
				file_CON << Continuity(i, j) << "\t";
            file_CON << "" << std::endl;
		}
	}

	if(oneNN.get_row() > 0 && oneNN.get_col() > 0)
	{
		file_1NN.clean();
        file_1NN << oneNN.matrix[0][0] << std::endl;
        file_1NN << oneNN.matrix[1][0] << std::endl;
	} else
        std::cout << "Warning: Have not computed 1NN yet" << std::endl;

    std::cout << "Please see the help for details of the introduction of output files." << std::endl;
}

void NLDR::make_output_file_names(String &filename_COR, String &filename_DIS, String &filename_STR, String &filename_TIM, String &filename_TRU, String &filename_CON, String &filename_1NN)
{
	if(cost_function == (String) "CLASSIC_MDS")
	{
		filename_COR = make_filename(D_prefname, dim_str, cost_function, "COR", "SVD", "");
		filename_DIS = make_filename(D_prefname, dim_str, cost_function, "DIS", "SVD", "");
		filename_STR = make_filename(D_prefname, dim_str, cost_function, "STR", "SVD", "");
		filename_TIM = make_filename(D_prefname, dim_str, cost_function, "TIM", "SVD", "");
		filename_TRU = make_filename(D_prefname, dim_str, cost_function, "TRU", "SVD", "");
		filename_CON = make_filename(D_prefname, dim_str, cost_function, "CON", "SVD", "");
		filename_1NN = make_filename(D_prefname, dim_str, cost_function, "1NN", "SVD", "");
	}else
	{
		filename_COR = make_filename(D_prefname, dim_str, cost_function, "COR", algorithm, file_flag);
		filename_DIS = make_filename(D_prefname, dim_str, cost_function, "DIS", algorithm, file_flag);
		filename_STR = make_filename(D_prefname, dim_str, cost_function, "STR", algorithm, file_flag);
		filename_TIM = make_filename(D_prefname, dim_str, cost_function, "TIM", algorithm, file_flag);
		filename_TRU = make_filename(D_prefname, dim_str, cost_function, "TRU", algorithm, file_flag);
		filename_CON = make_filename(D_prefname, dim_str, cost_function, "CON", algorithm, file_flag);
		filename_1NN = make_filename(D_prefname, dim_str, cost_function, "1NN", algorithm, file_flag);
	}
};

void NLDR::NLDR_load_D()
{
	String fname = D_prefname;
	fname += ".";
    fname += D_postfname;
	File D_file(fname);
	if(!D_file.is_open())
	{
        std::cout << "Error: File \"" << fname << "\" cannot be opened! Please check if this file exists or is readable." << std::endl;
        exit(0);
    }
	D_file.seek(0);

    if(parameters.distance_file_type == 0)
    {
        D.resize(size, size);
        for(int i = 0; i < size; i++)
            for(int j = 0; j <= i; j++)
            {
                D_file >> D.matrix[i][j];
                D.matrix[j][i] = D.matrix[i][j];
            }
    } else
    if(parameters.distance_file_type == 1)
    {
        String tree;
        double index;
        size--;
        D.resize(size, size);
        D_file >> tree;
        for(int i = 0; i < size; i++)
            D_file >> index;
        for(int i = 0; i < size; i++)
        {
            D_file >> index;
            for(int j = 0; j <= i; j++)
            {
                D_file >> D.matrix[i][j];
                D.matrix[j][i] = D.matrix[i][j];
            }
        }
    } else
    {
        std::cout << "Warning: This distance file format does not exist in this program!" << std::endl;
    }
};

void NLDR::NLDR_init_X(String init_md, long seed)
{
	File mds_result;
	String init_fname;
    int dim_int = atoi(dim_str);
    X.resize(size, dim_int);

	if(init_md == (String) "RAND")
    {
		if(seed == -1)
			init_genrand((unsigned) time(NULL));
		else
            init_genrand(seed);
		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim_int; j++)
				X.matrix[i][j] = genrand_real2() * (parameters.random_end - parameters.random_start) + parameters.random_start;
	} else
	if(init_md == (String) "CLASSIC_MDS")
	{
		init_fname = make_filename(D_prefname, dim_str, "CLASSIC_MDS", "COR", "SVD", "");
		mds_result.open(init_fname);
		if(mds_result.is_open() && cost_function != (String) "CLASSIC_MDS")
		{
			for(int i = 0; i < size; i++)
				for(int j = 0; j < dim_int; j++)
					mds_result >> X.matrix[i][j];
		}
		if(!mds_result.is_open() && cost_function != (String) "CLASSIC_MDS")
		{
			this->CLASSIC_MDS();
			X = COR;
		}
	}
};

String NLDR::make_filename(String name_D, String dimension, String cost_f, String output, String algorithm, String flag)
{
	String result;
	result = name_D;
	result += "_";
	result += dimension;
	result += "D_";
	result += cost_f;
	result += "_";
	result += output;
	result += "_";
	result += algorithm;
	result += flag;
	result += ".out";
	return result;
}
/*
void NLDR::NLDR_init_parameters(String para_filename)
{
	parameters.length_tru = 11;
    parameters.interval_tru = 5;
    parameters.length_con = 11;
	parameters.iterval_con = 5;
	parameters.random_start = -100;
	parameters.random_end = 100;

	parameters.KRU_LIN_e = 0.00001;
	parameters.KRU_LIN_max_iter = 1000;
	parameters.KRU_LIN_run_time = 0;
	parameters.KRU_MAJ_e = 0.00001;
	parameters.KRU_MAJ_max_iter = 1000;
	parameters.KRU_MAJ_run_time = 0;
	parameters.KRU_GAU_e = 0.0001;
	parameters.KRU_GAU_max_iter = 1000;
	parameters.KRU_GAU_run_time = 0;
	parameters.KRU_STO_epochs = 100;
	parameters.KRU_STO_alpha0 = 0.5;
	parameters.KRU_STO_alphan = 0.01;

	parameters.NOR_MAJ_e = 0.00001;
	parameters.NOR_MAJ_max_iter = 1000;
	parameters.NOR_MAJ_run_time = 0;
	parameters.NOR_GAU_e = 0.0001;
	parameters.NOR_GAU_max_iter = 1000;
	parameters.NOR_GAU_run_time = 0;
	parameters.NOR_STO_epochs = 100;
	parameters.NOR_STO_alpha0 = 0.25;
	parameters.NOR_STO_alphan = 0.01;

	parameters.NLM_MAJ_e = 0.0001;
	parameters.NLM_MAJ_max_iter = 1000;
	parameters.NLM_MAJ_run_time = 0;
	parameters.NLM_GAU_e = 0.0001;
	parameters.NLM_GAU_max_iter = 1000;
	parameters.NLM_GAU_run_time = 0;
	parameters.NLM_STO_epochs = 100;
	parameters.NLM_STO_alpha0 = 0.25;
	parameters.NLM_STO_alphan = 0.01;

	parameters.CCA_MAJ_e = 0.0001;
	parameters.CCA_MAJ_max_iter = 1000;
	parameters.CCA_MAJ_run_time = 0;
	parameters.CCA_MAJ_lambda0 = 12000;
	parameters.CCA_MAJ_lambdan = 1;
	parameters.CCA_GAU_e = 0.0001;
	parameters.CCA_GAU_max_iter = 1000;
	parameters.CCA_GAU_run_time = 0;
	parameters.CCA_GAU_lambda0 = 12000;
	parameters.CCA_GAU_lambdan = 1;
	parameters.CCA_STO_epochs = 100;
	parameters.CCA_STO_lambda0 = 12000;
	parameters.CCA_STO_lambdan = 1;
	parameters.CCA_STO_alpha0 = 0.5;
	parameters.CCA_STO_alphan = 0.01;

	if(para_filename == (String) "")
		return;

	Array<Mapping<Mix, Mix> > para_arr;
	Mapping<Mix, Mix> para;

	File fpara(para_filename);
	if(!fpara.is_open())
	{
        std::cout << "error: file \"" << para_filename << "\" can not be open, please check if this file exists or is readable." << std::endl;
		exit(0);
	}

	para_arr = LPIMPORT_FORM->Read(para_filename);
	int num = para_arr.get_length();
	Mapping<Mix, Mix> m;
	Mix Key, Value;
	for(int i = 0; i < num; i++)
	{
		m = para_arr[i];
		Key = m[(String) "name"];
		Value = m[(String) "value"];
		para[Key] = Value;
	}

	if(para.Include_Key((String) "length_tru"))
		parameters.length_tru = para[(String) "length_tru"];
	if(para.Include_Key((String) "interval_tru"))
		parameters.interval_tru = para[(String) "interval_tru"];
	if(para.Include_Key((String) "length_con"))
		parameters.length_con = para[(String) "length_con"];
	if(para.Include_Key((String) "iterval_con"))
		parameters.iterval_con = para[(String) "iterval_con"];
	if(para.Include_Key((String) "random_start"))
		parameters.random_start = para[(String) "random_start"];
	if(para.Include_Key((String) "random_end"))
		parameters.random_end = para[(String) "random_end"];

	if(para.Include_Key((String) "KRU_LIN_e"))
		parameters.KRU_LIN_e = para[(String) "KRU_LIN_e"];
	if(para.Include_Key((String) "KRU_LIN_max_iter"))
		parameters.KRU_LIN_max_iter = para[(String) "KRU_LIN_max_iter"];
	if(para.Include_Key((String) "KRU_LIN_run_time"))
		parameters.KRU_LIN_run_time = (long) ((int) para[(String) "KRU_LIN_run_time"]);
	if(para.Include_Key((String) "KRU_MAJ_e"))
		parameters.KRU_MAJ_e = para[(String) "KRU_MAJ_e"];
	if(para.Include_Key((String) "KRU_MAJ_max_iter"))
		parameters.KRU_MAJ_max_iter = para[(String) "KRU_MAJ_max_iter"];
	if(para.Include_Key((String) "KRU_MAJ_run_time"))
		parameters.KRU_MAJ_run_time = (long) ((int) para[(String) "KRU_MAJ_run_time"]);
	if(para.Include_Key((String) "KRU_GAU_e"))
		parameters.KRU_GAU_e = para[(String) "KRU_GAU_e"];
	if(para.Include_Key((String) "KRU_GAU_max_iter"))
		parameters.KRU_GAU_max_iter = para[(String) "KRU_GAU_max_iter"];
	if(para.Include_Key((String) "KRU_GAU_run_time"))
		parameters.KRU_GAU_run_time = (long) ((int)para[(String) "KRU_GAU_run_time"]);
	if(para.Include_Key((String) "KRU_STO_epochs"))
		parameters.KRU_STO_epochs = para[(String) "KRU_STO_epochs"];
	if(para.Include_Key((String) "KRU_STO_alpha0"))
		parameters.KRU_STO_alpha0 = para[(String) "KRU_STO_alpha0"];
	if(para.Include_Key((String) "KRU_STO_alphan"))
		parameters.KRU_STO_alphan = para[(String) "KRU_STO_alphan"];

	if(para.Include_Key((String) "NOR_MAJ_e"))
		parameters.NOR_MAJ_e = para[(String) "NOR_MAJ_e"];
	if(para.Include_Key((String) "NOR_MAJ_max_iter"))
		parameters.NOR_MAJ_max_iter = para[(String) "NOR_MAJ_max_iter"];
	if(para.Include_Key((String) "NOR_MAJ_run_time"))
		parameters.NOR_MAJ_run_time = (long) ((int) para[(String) "NOR_MAJ_run_time"]);
	if(para.Include_Key((String) "NOR_GAU_e"))
		parameters.NOR_GAU_e = para[(String) "NOR_GAU_e"];
	if(para.Include_Key((String) "NOR_GAU_max_iter"))
		parameters.NOR_GAU_max_iter = para[(String) "NOR_GAU_max_iter"];
	if(para.Include_Key((String) "NOR_GAU_run_time"))
		parameters.NOR_GAU_run_time = (long) ((int) para[(String) "NOR_GAU_run_time"]);
	if(para.Include_Key((String) "NOR_STO_epochs"))
		parameters.NOR_STO_epochs = para[(String) "NOR_STO_epochs"];
	if(para.Include_Key((String) "NOR_STO_alpha0"))
		parameters.NOR_STO_alpha0 = para[(String) "NOR_STO_alpha0"];
	if(para.Include_Key((String) "NOR_STO_alphan"))
		parameters.NOR_STO_alphan = para[(String) "NOR_STO_alphan"];

	if(para.Include_Key((String) "NLM_MAJ_e"))
		parameters.NLM_MAJ_e = para[(String) "NLM_MAJ_e"];
	if(para.Include_Key((String) "NLM_MAJ_max_iter"))
		parameters.NLM_MAJ_max_iter = para[(String) "NLM_MAJ_max_iter"];
	if(para.Include_Key((String) "NLM_MAJ_run_time"))
		parameters.NLM_MAJ_run_time = (long) ((int) para[(String) "NLM_MAJ_run_time"]);
	if(para.Include_Key((String) "NLM_GAU_e"))
		parameters.NLM_GAU_e = para[(String) "NLM_GAU_e"];
	if(para.Include_Key((String) "NLM_GAU_max_iter"))
		parameters.NLM_GAU_max_iter = para[(String) "NLM_GAU_max_iter"];
	if(para.Include_Key((String) "NLM_GAU_run_time"))
		parameters.NLM_GAU_run_time = (long) ((int) para[(String) "NLM_GAU_run_time"]);
	if(para.Include_Key((String) "NLM_STO_epochs"))
		parameters.NLM_STO_epochs = para[(String) "NLM_STO_epochs"];
	if(para.Include_Key((String) "NLM_STO_alpha0"))
		parameters.NLM_STO_alpha0 = para[(String) "NLM_STO_alpha0"];
	if(para.Include_Key((String) "NLM_STO_alphan"))
		parameters.NLM_STO_alphan = para[(String) "NLM_STO_alphan"];

	if(para.Include_Key((String) "CCA_MAJ_e"))
		parameters.CCA_MAJ_e = para[(String) "CCA_MAJ_e"];
	if(para.Include_Key((String) "CCA_MAJ_max_iter"))
		parameters.CCA_MAJ_max_iter = para[(String) "CCA_MAJ_max_iter"];
	if(para.Include_Key((String) "CCA_MAJ_run_time"))
		parameters.CCA_MAJ_run_time = (long) ((int) para[(String) "CCA_MAJ_run_time"]);
	if(para.Include_Key((String) "CCA_MAJ_lambda0"))
		parameters.CCA_MAJ_lambda0 = para[(String) "CCA_MAJ_lambda0"];
	if(para.Include_Key((String) "CCA_MAJ_lambdan"))
		parameters.CCA_MAJ_lambdan = para[(String) "CCA_MAJ_lambdan"];
	if(para.Include_Key((String) "CCA_GAU_e"))
		parameters.CCA_GAU_e = para[(String) "CCA_GAU_e"];
	if(para.Include_Key((String) "CCA_GAU_max_iter"))
		parameters.CCA_GAU_max_iter = para[(String) "CCA_GAU_max_iter"];
	if(para.Include_Key((String) "CCA_GAU_run_time"))
		parameters.CCA_GAU_run_time = (long) ((int) para[(String) "CCA_GAU_run_time"]);
	if(para.Include_Key((String) "CCA_GAU_lambda0"))
		parameters.CCA_GAU_lambda0 = para[(String) "CCA_GAU_lambda0"];
	if(para.Include_Key((String) "CCA_GAU_lambdan"))
		parameters.CCA_GAU_lambdan = para[(String) "CCA_GAU_lambdan"];
	if(para.Include_Key((String) "CCA_STO_epochs"))
		parameters.CCA_STO_epochs = para[(String) "CCA_STO_epochs"];
	if(para.Include_Key((String) "CCA_STO_lambda0"))
		parameters.CCA_STO_lambda0 = para[(String) "CCA_STO_lambda0"];
	if(para.Include_Key((String) "CCA_STO_lambdan"))
		parameters.CCA_STO_lambdan = para[(String) "CCA_STO_lambdan"];
	if(para.Include_Key((String) "CCA_STO_alpha0"))
		parameters.CCA_STO_alpha0 = para[(String) "CCA_STO_alpha0"];
	if(para.Include_Key((String) "CCA_STO_alphan"))
		parameters.CCA_STO_alphan = para[(String) "CCA_STO_alphan"];
}
*/

// zd_coded: the code for extracting COR from SVD is causing memory leak.
// It has been fixed now. zd --2/21/20
void NLDR::CLASSIC_MDS()
{
	String filename_U = make_filename(D_prefname, "ALL",cost_function, "U", "SVD", "");
	String filename_S = make_filename(D_prefname, "ALL", cost_function, "S", "SVD", "");
	String filename_Vt = make_filename(D_prefname, "ALL", cost_function, "Vt", "SVD", "");

	File file_U(filename_U);
	File file_S(filename_S);
	File file_Vt(filename_Vt);
	Matrix<double> Scalar_Product;
	Scalar_Product = D.compute_scalar_product_matrix();
	U.resize(size, size), S.resize(size, size), Vt.resize(size, size);

	if(!file_U.is_open() || !file_S.is_open() || !file_Vt.is_open())
        {
                if(!Scalar_Product.SVD_LIB(U, S, Vt))
                {
                        std::cout << "Error: Singular Value Decomposition failed." << std::endl;
                        exit(0);
                }

                file_U.clean();
                file_S.clean();
                file_Vt.clean();

		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < size; j++)
			{
				file_U << U(i, j) << "\t";
				file_Vt << Vt(i, j) << "\t";
			}
            file_S << S(i, i) << std::endl;
            file_U << "" << std::endl;
            file_Vt << "" << std::endl;
		}
	}
	else
	{
		for(int i = 0; i < size; i++)
			file_S >> S(i, i);
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
			{
				file_Vt >> Vt(i, j);
				file_U >> U(i, j);
			}
	}

	COR.resize(size, atoi(dim_str));

	int flag = -1;
	for(int i = 0; i < atoi(dim_str); i++)
	{
		flag++;
		while(flag < size && Vt(flag, 0) * U(0, flag) < 0) // Skip eigenvectors associated to negative eigenvalues.
			flag++;
		if(flag >= size)
			break;
		for(int j = 0; j < size; j++)
		{
			COR(j, i) = sqrt(S(flag, flag)) * Vt(flag, j);
		}
	}

	DIS = COR.compute_Distance_Matrix();

	STRESS = CLASSIC_MDS_stress_function(Scalar_Product, DIS.compute_scalar_product_matrix());
}

double NLDR::CLASSIC_MDS_stress_function(const Matrix<double> &S1, const Matrix<double> &S2)
{
	double stress_value = 0;

	double dif = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
		{
			dif = S1.matrix[i][j] - S2.matrix[i][j];
			stress_value += dif * dif;
		}

	return stress_value;
};

void NLDR::KRUSKAL1()
{
	if(algorithm == (String) "LINEAR_ITERATION")
		this->KRUSKAL1_LINEAR_ITERATION();
	else
	if(algorithm == (String) "MAJORIZATION")
		this->KRUSKAL1_MAJORIZATION();
	else
	if(algorithm == (String) "GAUSS_SEIDEL")
		this->KRUSKAL1_GAUSS_SEIDEL();
    else
    if(algorithm == (String) "STOCHASTIC")
        this->KRUSKAL1_STOCHASTIC();
    else
    if(algorithm == (String) "METROPOLIS")
        this->KRUSKAL1_METROPOLIS();
	else
	{
        std::cout << "Error: Invalid algorithm for KRUSKAL1!" << std::endl;
        std::cout << "For KRUSKAL1, the available algorithms are LINEAR_ITERATION, MAJORIZATION, GAUSS_SEIDEL and STOCHASTIC." << std::endl;
		exit(0);
	}

};

void NLDR::KRUSKAL1_LINEAR_ITERATION()
{
	int dim = atoi(dim_str);
	Matrix<double> COR_change(size, dim);
	Matrix<double> d(dim, 1);
	COR = X;
	double d_length = 0;
	double resid = 0;
	double stress_diff_accum = 0;
	double stress_norm_accum = 0;
	double err = 1000;
	double times = 0;
	double stress1 = 0;
	double stress2 = 0;
	double step_size = (double) 1 / size;

	long start = clock();
	long cur_time = start;
	long run_time = parameters.KRU_LIN_run_time;

	while(err > parameters.KRU_LIN_e && (times < parameters.KRU_LIN_max_iter || !parameters.KRU_LIN_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		stress_diff_accum = 0;
		stress_norm_accum = 0;
		COR_change.generate(ZERO_GEN);
		for(int i = 0; i < size; i++)
			for(int j = 0; j < i; j++)
			{
				for(int k = 0; k < dim; k++)
					d.matrix[k][0] = COR.matrix[j][k] - COR.matrix[i][k];

				d_length = 0;
				for(int k = 0; k < dim; k++)
					d_length += d.matrix[k][0] * d.matrix[k][0];
				d_length = sqrt(d_length);
				if(d_length > 0.00000001)
				{
					resid = d_length - D.matrix[i][j];
					for(int k = 0; k < dim; k++)
						d.matrix[k][0] *= resid / d_length;

					for(int k = 0; k < dim; k++)
					{
						COR_change.matrix[i][k] += d.matrix[k][0];
						COR_change.matrix[j][k] -= d.matrix[k][0];
					}
					stress_diff_accum += resid * resid;
					stress_norm_accum += d_length * d_length;
				}
			}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim; j++)
				COR.matrix[i][j] += step_size * COR_change.matrix[i][j];

		stress2 = stress_diff_accum / stress_norm_accum;
        std::cout << "stress:" << stress2 << std::endl;
		err = fabs((stress1 - stress2) / stress2);
        std::cout << "rel error : " << err << std::endl;
		stress1 = stress2;
		times++;
        std::cout << "times : " << times << std::endl;
	}
	DIS = COR.compute_Distance_Matrix();
	STRESS = stress2;
	return;
};

void NLDR::KRUSKAL1_MAJORIZATION()
{
	double stress1 = 0;
	double stress2 = 0;
	double err = 1000;
	int times = 0;
	int dim_int = atoi(dim_str);
	COR = X;

	long start = clock();
	long cur_time = start;
	long run_time = parameters.KRU_MAJ_run_time;

	Matrix<double> DX(size, size);
	Matrix<double> BX(size, size);
	stress1 = this->KRUSKAL1_stress_function(DX);
    std::cout << "stress:" << stress1 << std::endl;

	this->KRUSKAL1_compute_BX(DX, BX);

	while(err > parameters.KRU_MAJ_e && (times < parameters.KRU_MAJ_max_iter || !parameters.KRU_MAJ_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		COR = this->KRUSKAL1_compute_Z(BX);
		stress2 = this->KRUSKAL1_rescale(DX);
        std::cout << "stress:" << stress2 << std::endl;
		err = fabs(((double) stress1 - stress2) / stress2);
        std::cout << "rel error : " << err << std::endl;
		stress1 = stress2;
		this->KRUSKAL1_compute_BX(DX, BX);
		times++;
        std::cout << "times : " << times << std::endl;
	}

	DIS = COR.compute_Distance_Matrix();

	STRESS = this->KRUSKAL1_stress_function(DIS);

	return;
};

Matrix<double> NLDR::KRUSKAL1_compute_Z(const Matrix<double> &BX)
{
	int dim = atoi(dim_str);
	Matrix<double> Z(size, dim);

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			Z.matrix[i][j] = 0;
			for(int k = 0; k < size; k++)
			{
				Z.matrix[i][j] += BX.matrix[i][k] * COR.matrix[k][j];
			}
			Z.matrix[i][j] /= size;
		}

	return Z;
};

double NLDR::KRUSKAL1_stress_function(Matrix<double> &DX)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double stress_value = 0;
	double eta = 0;
	double rho = 0;
	double eta_d = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			DX.matrix[i][j] = 0;
			DX.matrix[j][i] = 0;
			for(int k = 0; k < dim; k++)
			{
				dif = COR.matrix[i][k] - COR.matrix[j][k];
				DX.matrix[i][j] += dif * dif;
			}
			DX.matrix[i][j] = sqrt(DX.matrix[i][j]);
			DX.matrix[j][i] = DX.matrix[i][j];
			eta += DX.matrix[j][i] * DX.matrix[j][i];
			rho += DX.matrix[j][i] * D.matrix[j][i];
			eta_d += D.matrix[j][i] * D.matrix[j][i];
		}
	para1 = eta;
	para2 = eta_d;
	para3 = rho;

	stress_value = (eta_d + eta - 2 * rho) / (eta);

	return stress_value;
};

double NLDR::KRUSKAL1_rescale(Matrix<double> &DX)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double stress_value = 0;
	double eta = 0;
	double rho = 0;
	double eta_d = 0;
	double beta = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			DX.matrix[i][j] = 0;
			DX.matrix[j][i] = 0;
			for(int k = 0; k < dim; k++)
			{
				dif = COR.matrix[i][k] - COR.matrix[j][k];
				DX.matrix[i][j] += dif * dif;
			}
			DX.matrix[i][j] = sqrt(DX.matrix[i][j]);
			DX.matrix[j][i] = DX.matrix[i][j];
			eta += DX.matrix[j][i] * DX.matrix[j][i];
			rho += DX.matrix[j][i] * D.matrix[j][i];
			eta_d += D.matrix[j][i] * D.matrix[j][i];
		}
	beta = eta_d / rho;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			DX.matrix[i][j] *= beta;
			DX.matrix[j][i] = DX.matrix[i][j];
		}
		for(int j = 0; j < dim; j++)
			COR.matrix[i][j] *= beta;
	}

	stress_value = (eta_d + eta * beta * beta - 2 * rho * beta) / (eta * beta * beta);

	return stress_value;
};

double NLDR::KRUSKAL1_stress_function_by_matrix(const Matrix<double> &DX)
{
	return 0;
};

void NLDR::KRUSKAL1_compute_BX(const Matrix<double> &DX, Matrix<double> &BX)
{
	double *row_sum = new double[size];
	for(int i = 0; i < size; i++)
		row_sum[i] = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			if(DX.matrix[i][j] != 0)
			{
				BX.matrix[i][j] = -D.matrix[i][j] / DX.matrix[i][j];
				BX.matrix[j][i] = BX.matrix[i][j];
				row_sum[i] += BX.matrix[i][j];
				row_sum[j] += BX.matrix[j][i];
			}
		}

	for(int i = 0; i < size; i++)
		BX.matrix[i][i] = - row_sum[i];

	delete [] row_sum;
};

void NLDR::KRUSKAL1_GAUSS_SEIDEL()
{
	double alpha = 0.2;
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	int times = 0;
	int dim = atoi(dim_str);
	double eta = 0;
	double rho = 0;
	double eta_d = 0;
	double sum_dtx_x;
	double sum_x;
	double n1, n2, n3, n4;
	Matrix<double> DX(size, size);
	Matrix<double> XD(size, dim);
	Matrix<double> XD2(size, dim);
	Matrix<double> DTX(size, 1);
	Matrix<double> DTX3(size, 1);
	Matrix<double> e1(dim, 1);
	Matrix<double> e2(dim, 1);
	Matrix<double> c(dim, 1);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.KRU_GAU_run_time;
	COR = X;

	stress2 = KRUSKAL1_stress_function(DX);
    std::cout << "\nstress value : " << stress2 << std::endl;

	eta = para1;
	eta_d = para2;
	rho = para3;

	while(err > parameters.KRU_GAU_e && (times < parameters.KRU_GAU_max_iter || !parameters.KRU_GAU_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		stress1 = stress2;

		for(int j = 0; j < size; j++)
		{
			for(int h = 0; h < size; h++)
				for(int k = 0; k < dim; k++)
				{
					XD.matrix[h][k] = COR.matrix[j][k] - COR.matrix[h][k];
					XD2.matrix[h][k] = XD.matrix[h][k] * XD.matrix[h][k];
				}

			for(int h = 0; h < size; h++)
			{
				if(fabs(DX.matrix[j][h]) > 0.00000001)
				{
					DTX.matrix[h][0] = D.matrix[j][h] / DX.matrix[j][h];
					DTX3.matrix[h][0] = D.matrix[j][h] / DX.matrix[j][h] / DX.matrix[j][h] / DX.matrix[j][h];
				} else
				{
					DTX.matrix[h][0] = 0;
					DTX3.matrix[h][0] = 0;
				}
			}

			for(int h = 0; h < dim; h++)
			{
				sum_dtx_x = 0;
				for(int k = 0; k < size; k++)
					sum_dtx_x += DTX.matrix[k][0] * XD.matrix[k][h];
				n1 = sum_dtx_x * eta;
				sum_x = 0;
				for(int k = 0; k < size; k++)
					sum_x += XD.matrix[k][h];
				n2 = sum_x * eta_d;
				n3 = rho * sum_x;
				e1.matrix[h][0] = -2 * (n1 + n2 - 2 * n3) / eta / eta;
				n4 = 0;
				for(int k = 0; k < size; k++)
					n4 += - DTX3.matrix[k][0] * XD2.matrix[k][h] + DTX.matrix[k][0];
				n4 *= eta;
				e2.matrix[h][0] = -2 * ((n4 + (size - 1) * eta_d - 2 * (size - 1) * rho) * eta * eta - (n1 + n2 - 2 * n3) * 2 * eta * sum_x * 2) / eta / eta / eta / eta;
			}

			for(int h = 0; h < dim; h++)
			{
				if(! this->search_first_condition(DX, stress2, - e1.matrix[h][0] / fabs(e2.matrix[h][0]), j, h, - e1.matrix[h][0] * e1.matrix[h][0] / fabs(e2.matrix[h][0]), 0, alpha, 1, 1, stress2))
                    ;//---std::cout << "!";
				COR.matrix[j][h] = COR.matrix[j][h] - alpha * e1.matrix[h][0] / fabs(e2.matrix[h][0]);
				this->update_distance_matrix(DX, j, h, COR.matrix[j][h], 0);
				eta = para1;
				eta_d = para2;
				rho = para3;
			}
		}

		for(int i = 0; i < dim; i++)
		{
			c.matrix[i][0] = 0;
			for(int j = 0; j < size; j++)
			{
				c.matrix[i][0] = c.matrix[i][0] + COR.matrix[j][i];
			}
			c.matrix[i][0] = c.matrix[i][0] / size;
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim; j++)
			{
				COR.matrix[i][j] = COR.matrix[i][j] - c.matrix[j][0];
			}

		stress2 = KRUSKAL1_stress_function(DX);
		err = fabs((stress1 - stress2) / stress2);
        std::cout << "\nstress value : " << stress2 << std::endl;
        std::cout << "rel error : " << err << std::endl;

		times = times + 1;
        std::cout << "times : " << times << std::endl;
	}

    std::cout << "\nstress value : " << stress2 << std::endl;

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::KRUSKAL1_STOCHASTIC()
{
	double eta = 0;
	double eta_d = 0;
	double rho = 0;
	double n1, n2, n3;
	int epochs = parameters.KRU_STO_epochs;
	int train_len = epochs * size;
	int ind = 0;
	int dim = atoi(dim_str);
	double alpha0 = parameters.KRU_STO_alpha0;
	double stress_value;
	Matrix<double> alpha(train_len, 1);
	Matrix<int> sample_inds(train_len, 1);
	Matrix<double> DX(size, size);
	Matrix<double> DE(size, dim);
	COR = X;

	stress_value = KRUSKAL1_stress_function(DX);
    std::cout << "\nstress value : " << stress_value << std::endl;
	eta = para1;
	eta_d = para2;
	rho = para3;

	for(int i = 0; i < train_len; i++)
	{
		alpha.matrix[i][0] = alpha0 * eta_d * pow((double) parameters.KRU_STO_alphan, (double) i / (train_len - 1));
	}

	for(int i = 0; i < train_len; i++)
		sample_inds.matrix[i][0] = genrand_int32() % size;

	for(int i = 0; i < train_len; i++)
	{
		ind = sample_inds.matrix[i][0];

		n2 = 0;
		for(int h = 0; h < size; h++)
			n2 += (D.matrix[ind][h] - DX.matrix[ind][h]) * (D.matrix[ind][h] - DX.matrix[ind][h]);

		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
			{
				if(DX.matrix[ind][j] > 0.0000001)
					n1 = - 2 * (D.matrix[ind][j] - DX.matrix[ind][j]) / DX.matrix[ind][j] * (COR.matrix[j][k] - COR.matrix[ind][k]) * eta;
				else
					n1 = 0;
				n3 = 0;
				for(int h = 0; h < size; h++)
					n3 += COR.matrix[j][k] - COR.matrix[h][k];
				DE.matrix[j][k] = (n1 - n2 * n3) / eta / eta;
			}

		//--- update ---
		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
			{
				COR.matrix[j][k] = COR.matrix[j][k] - alpha.matrix[i][0] * DE.matrix[j][k];
			}

		stress_value = KRUSKAL1_stress_function(DX);
		eta = para1;
		eta_d = para2;
		rho = para3;

		if(i % size == 0)
		{
            std::cout << "epochs : " << (int) i / size << ", stress : " << stress_value << std::endl;
		}
	}

	STRESS = stress_value;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::KRUSKAL1_METROPOLIS()
{
    int dim = atoi(dim_str);
    double T = 10000;
    double delta = 0;
    Matrix<double> gradient(size, dim);
    Matrix<double> candidate(size, dim);
    Matrix<double> DX(size, size);
    Matrix<double> ADX(size, size);
    Matrix<double> CORX(size, dim);
    double str1 = 0;
    double str2 = 0;
    CORX = X;
    COR = CORX;
    str1 = KRUSKAL1_stress_function(DX);
    std::cout << "start : accept stress : " << str1 << std::endl;
    ADX = DX;
    KRUSKAL1_GRADIENT(gradient, CORX, ADX);
    delta = gradient.norm(M_FRO_NORM);
    std::cout << "gradient norm:" << delta << std::endl;
    delta = (double) 1 / (delta * delta);
    std::cout << "delta:" << delta << std::endl;
    for(int i = 0; i < 2000; i++)
    {
        for(int j = 0; j < size; j++)
            for(int k = 0; k < dim; k++)
            {
                candidate.matrix[j][k] = CORX.matrix[j][k] - delta * gradient.matrix[j][k] + sqrt(T * delta / size / dim) * Normal_random_generator();
            }
        COR = candidate;
        str2 = KRUSKAL1_stress_function(DX);
        if(genrand_real2() < exp((str1 - str2) / (T)))
        {
            CORX = candidate;
            ADX = DX;
            str1 = str2;
            KRUSKAL1_GRADIENT(gradient, CORX, ADX);
            std::cout << i << " : accept stress : " << str1 << std::endl;
        }

        T = T * 0.99;
    }
    COR = CORX;
    STRESS = str1;
    DIS = COR.compute_Distance_Matrix();
};

void NLDR::KRUSKAL1_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX)
{
    int dim = atoi(dim_str);
    double n1, n2;
    double eta = para1;
    double eta_d = para2;
    double rho = para3;

    for(int i = 0; i < size; i++)
        for(int j = 0; j < dim; j++)
        {
            n1 = 0;
            n2 = 0;
            for(int k = 0; k < size; k++)
            {
                n1 += CORX.matrix[i][j] - CORX.matrix[k][j];
                if(ADX.matrix[i][k] > 0.00000001)
                    n2 += D.matrix[i][k] / ADX.matrix[i][k] * (CORX.matrix[i][j] - CORX.matrix[k][j]);
            }
            gradient.matrix[i][j] = (double) -2 / eta / eta * (n1 * eta + eta_d * n2 - 2 * rho * n2);
        }
};

void NLDR::NORMALIZED()
{
	if(algorithm == (String) "MAJORIZATION")
		this->NORMALIZED_MAJORIZATION();
	else
	if(algorithm == (String) "GAUSS_SEIDEL")
		this->NORMALIZED_GAUSS_SEIDEL();
	else
	if(algorithm == (String) "STOCHASTIC")
		this->NORMALIZED_STOCHASTIC();
    else
    if(algorithm == (String) "METROPOLIS")
        this->NORMALIZED_METROPOLIS();
	else
	{
        std::cout << "Error: Invalid algorithm for NORMALIZED STRESS!" << std::endl;
        std::cout << "For NORMALIZED STRESS, the available algorithms are MAJORIZATION, GAUSS_SEIDEL and STOCHASTIC." << std::endl;
		exit(0);
	}
};

void NLDR::NORMALIZED_MAJORIZATION()
{
	double stress1 = 0;
	double stress2 = 0;
	double err = 1000;
	int times = 0;
	int dim_int = atoi(dim_str);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.NOR_MAJ_run_time;
	COR = X;

	Matrix<double> DX(size, size);
	Matrix<double> BX(size, size);
	stress1 = this->NORMALIZED_stress_function(DX);
    std::cout << "stress:" << stress1 << std::endl;

	this->NORMALIZED_compute_BX(DX, BX);

	while(err > parameters.NOR_MAJ_e && (times < parameters.NOR_MAJ_max_iter || !parameters.NOR_MAJ_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		COR = this->NORMALIZED_compute_Z(BX);
		stress2 = this->NORMALIZED_stress_function(DX);
        std::cout << "stress:" << stress2 << std::endl;
		err = fabs(((double) stress1 - stress2) / stress2);
        std::cout << "rel error : " << err << std::endl;
		stress1 = stress2;
		this->NORMALIZED_compute_BX(DX, BX);
		times++;
        std::cout << "times : " << times << std::endl;
	}

	DIS = COR.compute_Distance_Matrix();

	STRESS = this->NORMALIZED_stress_function(DIS);

	return;
};

Matrix<double> NLDR::NORMALIZED_compute_Z(const Matrix<double> &BX)
{
	int dim = atoi(dim_str);
	Matrix<double> Z(size, dim);

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			Z.matrix[i][j] = 0;
			for(int k = 0; k < size; k++)
			{
				Z.matrix[i][j] += BX.matrix[i][k] * COR.matrix[k][j];
			}
			Z.matrix[i][j] /= size;
		}

	return Z;
}

double NLDR::NORMALIZED_stress_function(Matrix<double> &DX)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double stress_value = 0;
	double normalize_value = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			DX.matrix[i][j] = 0;
			DX.matrix[j][i] = 0;
			for(int k = 0; k < dim; k++)
			{
				dif = COR.matrix[i][k] - COR.matrix[j][k];
				DX.matrix[i][j] += dif * dif;
			}
			DX.matrix[i][j] = sqrt(DX.matrix[i][j]);
			DX.matrix[j][i] = DX.matrix[i][j];
			dif = DX.matrix[j][i] - D.matrix[j][i];
			stress_value += dif * dif;
			normalize_value += D.matrix[j][i] * D.matrix[j][i];
		}
	stress_value /= normalize_value;

	return stress_value;
}

double NLDR::NORMALIZED_stress_function_by_matrix(const Matrix<double> &DX)
{
	double dif = 0;
	double stress_value = 0;

	stress_value = pow((DX - D).norm(M_FRO_NORM), 2) / 2;

	return stress_value;
}

void NLDR::NORMALIZED_compute_BX(const Matrix<double> &DX, Matrix<double> &BX)
{
	double *row_sum = new double[size];
	for(int i = 0; i < size; i++)
		row_sum[i] = 0;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			if(DX.matrix[i][j] != 0)
			{
				BX.matrix[i][j] = -D.matrix[i][j] / DX.matrix[i][j];
				BX.matrix[j][i] = BX.matrix[i][j];
				row_sum[i] += BX.matrix[i][j];
				row_sum[j] += BX.matrix[j][i];
			}
		}

	for(int i = 0; i < size; i++)
		BX.matrix[i][i] = - row_sum[i];

	delete [] row_sum;
}

void NLDR::NORMALIZED_GAUSS_SEIDEL()
{
	double alpha = 0.2;
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	int times = 0;
	int dim = atoi(dim_str);
	double sum_dtx;
	double sum_dtx_x;
	double sum_dtx3_xx;
	double weight = 0;
	Matrix<double> DX(size, size);
	Matrix<double> XD(size, dim);
	Matrix<double> XD2(size, dim);
	Matrix<double> DTX(size, 1);
	Matrix<double> DTX3(size, 1);
	Matrix<double> e1(dim, 1);
	Matrix<double> e2(dim, 1);
	Matrix<double> c(dim, 1);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.NOR_GAU_run_time;
	COR = X;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
			weight += D.matrix[i][j] * D.matrix[i][j];

	stress2 = NORMALIZED_stress_function(DX);
    std::cout << "\nstress value : " << stress2 << std::endl;

	while(err > parameters.NOR_GAU_e && (times < parameters.NOR_GAU_max_iter || !parameters.NOR_GAU_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		stress1 = stress2;

		for(int j = 0; j < size; j++)
		{
			for(int h = 0; h < size; h++)
				for(int k = 0; k < dim; k++)
				{
					XD.matrix[h][k] = COR.matrix[j][k] - COR.matrix[h][k];
					XD2.matrix[h][k] = XD.matrix[h][k] * XD.matrix[h][k];
				}

			for(int h = 0; h < size; h++)
			{
				if(fabs(DX.matrix[j][h]) > 0.00000001)
				{
					DTX.matrix[h][0] = D.matrix[j][h] / DX.matrix[j][h];
					DTX3.matrix[h][0] = D.matrix[j][h] / DX.matrix[j][h] / DX.matrix[j][h] / DX.matrix[j][h];
				} else
				{
					DTX.matrix[h][0] = 0;
					DTX3.matrix[h][0] = 0;
				}
			}

			for(int h = 0; h < dim; h++)
			{
				sum_dtx = 0;
				for(int k = 0; k < size; k++)
					sum_dtx += DTX.matrix[k][0] - 1;
				sum_dtx_x = 0;
				for(int k = 0; k < size; k++)
					sum_dtx_x += (DTX.matrix[k][0] - 1) * XD.matrix[k][h];
				e1.matrix[h][0] = - 2 * sum_dtx_x / weight;

				sum_dtx3_xx = 0;
				for(int k = 0; k < size; k++)
					sum_dtx3_xx += DTX3.matrix[k][0] * XD2.matrix[k][h];

				e2.matrix[h][0] = - 2 * (sum_dtx - sum_dtx3_xx) / weight;
			}

			for(int h = 0; h < dim; h++)
			{
                if(! this->search_first_condition(DX, stress2, - e1.matrix[h][0] / fabs(e2.matrix[h][0]), j, h, - e1.matrix[h][0] * e1.matrix[h][0] / fabs(e2.matrix[h][0]), weight, alpha, 1, 1, stress2))
                    ;//---std::cout << "!";
				COR.matrix[j][h] = COR.matrix[j][h] - alpha * e1.matrix[h][0] / fabs(e2.matrix[h][0]);
				this->update_distance_matrix(DX, j, h, COR.matrix[j][h], 0);
			}
		}

		for(int i = 0; i < dim; i++)
		{
			c.matrix[i][0] = 0;
			for(int j = 0; j < size; j++)
			{
				c.matrix[i][0] = c.matrix[i][0] + COR.matrix[j][i];
			}
			c.matrix[i][0] = c.matrix[i][0] / size;
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim; j++)
			{
				COR.matrix[i][j] = COR.matrix[i][j] - c.matrix[j][0];
			}

		stress2 = NORMALIZED_stress_function(DX);
		err = fabs((stress1 - stress2) / stress2);
        std::cout << "\nstress value : " << stress2 << std::endl;
        std::cout << "rel error : " << err << std::endl;

		times = times + 1;
        std::cout << "times : " << times << std::endl;
	}

    std::cout << "\nstress value : " << stress2 << std::endl;

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::NORMALIZED_STOCHASTIC()
{
	int epochs = parameters.NOR_STO_epochs;
	int train_len = epochs * size;
	int ind = 0;
	int dim = atoi(dim_str);
	double alpha0 = 0;
	double stress_value;
	double weight = 0;
	Matrix<double> alpha(train_len, 1);
	Matrix<int> sample_inds(train_len, 1);
	Matrix<double> XD(size, dim);
	Matrix<double> DX(size, 1);
	Matrix<double> FY(size, 1);
	Matrix<double> result_dis(size, size);
	COR = X;
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
			weight += D.matrix[i][j] * D.matrix[i][j];
	alpha0 = weight * parameters.NOR_STO_alpha0;

	for(int i = 0; i < train_len; i++)
	{
		alpha.matrix[i][0] = alpha0 * pow((double) parameters.NOR_STO_alphan, (double) i / (train_len - 1));
	}

	for(int i = 0; i < train_len; i++)
		sample_inds.matrix[i][0] = genrand_int32() % size;

	for(int i = 0; i < train_len; i++)
	{
		ind = sample_inds.matrix[i][0];

		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
				XD.matrix[j][k] = COR.matrix[j][k] - COR.matrix[ind][k];

		for(int j = 0; j < size; j++)
		{
			DX.matrix[j][0] = 0;
			for(int k = 0; k < dim; k++)
			{
				DX.matrix[j][0] += XD.matrix[j][k] * XD.matrix[j][k];
			}
			DX.matrix[j][0] = sqrt((double) DX.matrix[j][0]);
			if(fabs(DX.matrix[j][0]) < 0.0000000001)
				DX.matrix[j][0] = 1;
			FY.matrix[j][0] = (D.matrix[j][ind] / DX.matrix[j][0] - 1);
		}

		//--- update ---
		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
				COR.matrix[j][k] = COR.matrix[j][k] + 2 / weight * alpha.matrix[i][0] * FY.matrix[j][0] * XD.matrix[j][k];

		if(i % size == 0)
		{
			stress_value = NORMALIZED_stress_function(result_dis);
            std::cout << "epochs : " << (int) i / size << ", stress : " << stress_value << std::endl;
		}
	}

	STRESS = stress_value;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::NORMALIZED_METROPOLIS()
{
    int dim = atoi(dim_str);
    double T = 10000;
    double delta = 0;
    Matrix<double> gradient(size, dim);
    Matrix<double> candidate(size, dim);
    Matrix<double> DX(size, size);
    Matrix<double> ADX(size, size);
    Matrix<double> CORX(size, dim);
    double str1 = 0;
    double str2 = 0;
    CORX = X;
    COR = CORX;
    str1 = NORMALIZED_stress_function(DX);
    std::cout << "start : accept stress : " << str1 << std::endl;
    ADX = DX;
    NORMALIZED_GRADIENT(gradient, CORX, ADX);
    delta = gradient.norm(M_FRO_NORM);
    std::cout << "gradient norm:" << delta << std::endl;
    delta = (double) 1 / (delta * delta);
    std::cout << "delta:" << delta << std::endl;
    for(int i = 0; i < 2000; i++)
    {
        for(int j = 0; j < size; j++)
            for(int k = 0; k < dim; k++)
            {
                candidate.matrix[j][k] = CORX.matrix[j][k] - delta * gradient.matrix[j][k] + sqrt(T * delta / size / dim) * Normal_random_generator();
            }
        COR = candidate;
        str2 = NORMALIZED_stress_function(DX);
        if(genrand_real2() < exp((str1 - str2) / (T)))
        {
            CORX = candidate;
            ADX = DX;
            str1 = str2;
            NORMALIZED_GRADIENT(gradient, CORX, ADX);
            std::cout << i << " : accept stress : " << str1 << std::endl;
        }

        T = T * 0.99;
    }
    COR = CORX;
    STRESS = str1;
    DIS = COR.compute_Distance_Matrix();
};

void NLDR::NORMALIZED_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX)
{
    int dim = atoi(dim_str);
    double c = 0;
    for(int i = 0; i < size; i++)
        for(int j = 0; j < i; j++)
            c += D.matrix[i][j] * D.matrix[i][j];

    for(int i = 0; i < size; i++)
        for(int j = 0; j < dim; j++)
        {
            gradient.matrix[i][j] = 0;
            for(int k = 0; k < size; k++)
            {
                if(ADX.matrix[i][k] > 0.00000001)
                    gradient.matrix[i][j] += (ADX.matrix[i][k] - D.matrix[i][k]) / ADX.matrix[i][k] * (CORX.matrix[i][j] - CORX.matrix[k][j]);
            }
            gradient.matrix[i][j] /= c;
        }
};

double NLDR::Normal_random_generator()
{
    static int status = 0;
    static double result1;
    static double result2;
    double u1;
    double u2;

    status = status % 2;
    if(status == 0)
    {
        u1 = genrand_real2();
        u2 = genrand_real2();
        result1 = sqrt(-2 * log(u1)) * cos(2 * 3.141592653589793 * u2);
        result2 = sqrt(-2 * log(u2)) * sin(2 * 3.141592653589793 * u2);
        status++;
        return result1;
    }
    status++;
    return result2;
};

void NLDR::SAMMON()
{
	if(algorithm == (String) "MAJORIZATION")
		this->SAMMON_MAJORIZATION();
	else
	if(algorithm == (String) "GAUSS_SEIDEL")
		this->SAMMON_GAUSS_SEIDEL();
	else
	if(algorithm == (String) "STOCHASTIC")
		this->SAMMON_STOCHASTIC();
    else
    if(algorithm == (String) "METROPOLIS")
        this->SAMMON_METROPOLIS();
	else
	{
        std::cout << "Error: Invalid algorithm for SAMMON's STRESS(NLM)!" << std::endl;
        std::cout << "For SAMMON's STRESS, the available algorithms are MAJORIZATION, GAUSS_SEIDEL and STOCHASTIC." << std::endl;
		exit(0);
	}
};

void NLDR::SAMMON_MAJORIZATION()
{
	Matrix<double> V(size, size);
	Matrix<double> DX(size, size);
	Matrix<double> BX(size, size);
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	int times = 0;
	long start = clock();
	long cur_time = start;
	long run_time = parameters.NLM_MAJ_run_time;
	COR = X;

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < i; j++)
		{
			if(D.matrix[i][j] < 0.00001)
				V.matrix[i][j] = - 100000;
			else
				V.matrix[i][j] = - (double) 1 / D.matrix[i][j];
			V.matrix[j][i] = V.matrix[i][j];
			V.matrix[i][i] = V.matrix[i][i] - V.matrix[i][j];
			V.matrix[j][j] = V.matrix[j][j] - V.matrix[j][i];
		}
	}

	V = V + 1;

	if(! V.compute_inverse_matrix(V))
	{
        std::cout << "Inverse computation failed!" << std::endl;
		return;
	}

	V = V - (double) 1 / size / size;

	stress1 = SAMMON_stress_function(DX);

    std::cout << "stress : " << stress1 << std::endl;

	SAMMON_compute_BX(DX, BX);

	while(err > parameters.NLM_MAJ_e && (times < parameters.NLM_MAJ_max_iter || !parameters.NLM_MAJ_max_iter) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		COR = SAMMON_compute_Z(BX, V);
		stress2 = SAMMON_stress_function(DX);
        std::cout << "stress : " << stress2 << std::endl;
		err = fabs(((double) stress1 - stress2) / stress2);
        std::cout << "rel error : " << err << std::endl;
		stress1 = stress2;
		SAMMON_compute_BX(DX, BX);
		times++;
        std::cout << "times : " << times << std::endl;
	}

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

Matrix<double> NLDR::SAMMON_compute_Z(const Matrix<double> &BX, const Matrix<double> &V)
{
	int dim = atoi(dim_str);
	double *A = new double[size * dim];
	Matrix<double> Z(size, dim);
	for(int i = 0; i < size * dim; i++)
		A[i] = 0;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			A[i * dim + j] = 0;
			for(int k = 0; k < size; k++)
			{
				A[i * dim + j] = A[i * dim + j] + BX.matrix[i][k] * COR.matrix[k][j];
			}
		}

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			Z.matrix[i][j] = 0;
			for(int k = 0; k < size; k++)
			{
				Z.matrix[i][j] = Z.matrix[i][j] + V.matrix[i][k] * A[k * dim + j];
			}
		}

	delete [] A;
	return Z;
};

double NLDR::SAMMON_stress_function(Matrix<double> &DX)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double stress_value = 0;
	double total = 0;
	for(int i = 0; i < size; i++)
	{
		DX.matrix[i][i] = 0;
		for(int j = 0; j < i; j++)
		{
			DX.matrix[i][j] = 0;
			DX.matrix[j][i] = 0;
			for(int k = 0; k < dim; k++)
			{
				dif = COR.matrix[i][k] - COR.matrix[j][k];
				DX.matrix[i][j] = DX.matrix[i][j] + dif * dif;
			}
			DX.matrix[i][j] = sqrt((double) DX.matrix[i][j]);
			DX.matrix[j][i] = DX.matrix[i][j];
			dif = DX.matrix[j][i] - D.matrix[j][i];
			if(D.matrix[j][i] > 0.00000001)
				stress_value = stress_value + dif * dif / D.matrix[j][i];
			total += D.matrix[i][j];
		}
	}

	return stress_value / total;
};

double NLDR::SAMMON_stress_function_by_matrix(const Matrix<double> &DX)
{
	double dif = 0;
	double stress_value = 0;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < i; j++)
		{
			dif = DX(j, i) - D(j, i);
			if(D(j, i) > 0.00000001)
				stress_value = stress_value + dif * dif / D(j, i);
		}
	}

	return stress_value;
};

void NLDR::SAMMON_compute_BX(const Matrix<double> &DX, Matrix<double> &BX)
{
	double *row_sum = new double[size];
	for(int i = 0; i < size; i++)
		row_sum[i] = 0;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			if(DX.matrix[i][j] != 0)
			{
				BX.matrix[i][j] = - (double) 1 / DX(i, j);
				BX.matrix[j][i] = BX.matrix[i][j];
				row_sum[i] = row_sum[i] + BX.matrix[i][j];
				row_sum[j] = row_sum[j] + BX.matrix[j][i];
			}
		}

	for(int i = 0; i < size; i++)
		BX.matrix[i][i] = - row_sum[i];

	delete [] row_sum;
};

void NLDR::SAMMON_GAUSS_SEIDEL()
{
	double alpha = 0.2;
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	int times = 0;
	int dim = atoi(dim_str);
	double sum_dtxt;
	double sum_dtx_xt;
	double sum_dx3_xx;
	double weight = 0;
	Matrix<double> DX(size, size);
	Matrix<double> XD(size, dim);
	Matrix<double> XD2(size, dim);
	Matrix<double> DTX(size, 1);
	Matrix<double> DTX3(size, 1);
	Matrix<double> e1(dim, 1);
	Matrix<double> e2(dim, 1);
	Matrix<double> c(dim, 1);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.NLM_GAU_run_time;
	COR = X;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
			weight += D.matrix[i][j];

	stress2 = SAMMON_stress_function(DX);
    std::cout << "\nstress value : " << stress2 << std::endl;

	while(err > parameters.NLM_GAU_e && (times < parameters.NLM_GAU_max_iter || !parameters.NLM_GAU_max_iter) && (cur_time - start < run_time || !run_time))
    {
		cur_time = clock();
		stress1 = stress2;

		for(int j = 0; j < size; j++)
		{
			for(int h = 0; h < size; h++)
				for(int k = 0; k < dim; k++)
				{
					XD.matrix[h][k] = COR.matrix[j][k] - COR.matrix[h][k];
					XD2.matrix[h][k] = XD.matrix[h][k] * XD.matrix[h][k];
                }

			for(int h = 0; h < size; h++)
			{
				if(fabs(DX.matrix[j][h]) > 0.00000001 && fabs(D.matrix[j][h]) > 0.00000001)
				{
					DTX.matrix[h][0] = (D.matrix[j][h] - DX.matrix[j][h]) / DX.matrix[j][h] / D.matrix[j][h];
                    DTX3.matrix[h][0] = (double) 1 / DX.matrix[j][h] / DX.matrix[j][h] / DX.matrix[j][h];
				} else
                {
					DTX.matrix[h][0] = 0;
					DTX3.matrix[h][0] = 0;
				}
            }

			for(int h = 0; h < dim; h++)
			{
				sum_dtxt = 0;
				for(int k = 0; k < size; k++)
                    sum_dtxt += DTX.matrix[k][0];
				sum_dtx_xt = 0;
				for(int k = 0; k < size; k++)
                    sum_dtx_xt += DTX.matrix[k][0] * XD.matrix[k][h];
                e1.matrix[h][0] = (double) - 2 * sum_dtx_xt / weight;

				sum_dx3_xx = 0;
				for(int k = 0; k < size; k++)
                    sum_dx3_xx += DTX3.matrix[k][0] * XD2.matrix[k][h];

                e2.matrix[h][0] = (double) - 2 * (sum_dtxt - sum_dx3_xx) / weight;
            }

			for(int h = 0; h < dim; h++)
            {
                if(! this->search_first_condition(DX, stress2, - e1.matrix[h][0] / fabs(e2.matrix[h][0]), j, h, - e1.matrix[h][0] * e1.matrix[h][0] / fabs(e2.matrix[h][0]), weight, alpha, 1, 1, stress2))
                    ;//---std::cout << "!";
				COR.matrix[j][h] = COR.matrix[j][h] - alpha * e1.matrix[h][0] / fabs(e2.matrix[h][0]);
				this->update_distance_matrix(DX, j, h, COR.matrix[j][h], 0);
			}
        }

		for(int i = 0; i < dim; i++)
		{
			c.matrix[i][0] = 0;
			for(int j = 0; j < size; j++)
			{
				c.matrix[i][0] = c.matrix[i][0] + COR.matrix[j][i];
			}
			c.matrix[i][0] = c.matrix[i][0] / size;
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim; j++)
			{
				COR.matrix[i][j] = COR.matrix[i][j] - c.matrix[j][0];
			}

        stress2 = SAMMON_stress_function(DX);
		err = fabs((stress1 - stress2) / stress2);
        std::cout << "\nstress value : " << stress2 << std::endl;
        std::cout << "rel error : " << err << std::endl;

		times = times + 1;
        std::cout << "times : " << times << std::endl;
	}

    std::cout << "\nstress value : " << stress2 << std::endl;

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::SAMMON_STOCHASTIC()
{
	int epochs = parameters.NLM_STO_epochs;
	int train_len = epochs * size;
	int ind = 0;
	int dim = atoi(dim_str);
	double alpha0 = 0;
	double stress_value = 0;
	double weight = 0;
	Matrix<double> alpha(train_len, 1);
	Matrix<int> sample_inds(train_len, 1);
	Matrix<double> XD(size, dim);
	Matrix<double> DX(size, 1);
	Matrix<double> FY(size, 1);
	Matrix<double> result_dis(size, size);
	COR = X;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
			weight += D.matrix[i][j];
	alpha0 = weight * parameters.NLM_STO_alpha0;

	for(int i = 0; i < train_len; i++)
		alpha.matrix[i][0] = alpha0 * pow((double) parameters.NLM_STO_alphan, (double) i / (train_len - 1));

	for(int i = 0; i < train_len; i++)
		sample_inds.matrix[i][0] = genrand_int32() % size;

	for(int i = 0; i < train_len; i++)
	{
		ind = sample_inds.matrix[i][0];

		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
				XD.matrix[j][k] = COR.matrix[j][k] - COR.matrix[ind][k];

		for(int j = 0; j < size; j++)
		{
			DX.matrix[j][0] = 0;
			for(int k = 0; k < dim; k++)
			{
				DX.matrix[j][0] += XD.matrix[j][k] * XD.matrix[j][k];
			}
			DX.matrix[j][0] = sqrt((double) DX.matrix[j][0]);
			if(fabs(DX.matrix[j][0]) < ZERO)
				DX.matrix[j][0] = 1;
		}

		for(int j = 0; j < size; j++)
		{
			if(fabs(D.matrix[j][ind]) > ZERO)
				FY.matrix[j][0] = 1 / DX.matrix[j][0] - 1 / D.matrix[j][ind];
			else
				FY.matrix[j][0] = 1 / DX.matrix[j][0] - 1;
		}
		//--- update ---
		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
			{
				COR.matrix[j][k] = COR.matrix[j][k] + alpha.matrix[i][0] * 2 / weight * FY.matrix[j][0] * XD.matrix[j][k];
			}
		if(i % size == 0)
		{
			stress_value = SAMMON_stress_function(result_dis);
            std::cout << "epochs : " << (int) i / size  << ", stress : " << stress_value << std::endl;
		}
	}

	STRESS = stress_value;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::SAMMON_METROPOLIS()
{
    int dim = atoi(dim_str);
    double T = 10000;
    double delta = 0;
    Matrix<double> gradient(size, dim);
    Matrix<double> candidate(size, dim);
    Matrix<double> DX(size, size);
    Matrix<double> ADX(size, size);
    Matrix<double> CORX(size, dim);
    double str1 = 0;
    double str2 = 0;
    CORX = X;
    COR = CORX;
    str1 = SAMMON_stress_function(DX);
    std::cout << "start : accept stress : " << str1 << std::endl;
    ADX = DX;
    SAMMON_GRADIENT(gradient, CORX, ADX);
    delta = gradient.norm(M_FRO_NORM);
    std::cout << "gradient norm:" << delta << std::endl;
    delta = (double) 1 / (delta * delta);
    std::cout << "delta:" << delta << std::endl;
    for(int i = 0; i < 2000; i++)
    {
        for(int j = 0; j < size; j++)
            for(int k = 0; k < dim; k++)
            {
                candidate.matrix[j][k] = CORX.matrix[j][k] - delta * gradient.matrix[j][k] + sqrt(T * delta / size / dim) * Normal_random_generator();
            }
        COR = candidate;
        str2 = SAMMON_stress_function(DX);
        if(genrand_real2() < exp((str1 - str2) / (T)))
        {
            CORX = candidate;
            ADX = DX;
            str1 = str2;
            SAMMON_GRADIENT(gradient, CORX, ADX);
            std::cout << i << " : accept stress : " << str1 << std::endl;
        }

        T = T * 0.99;
    }
    COR = CORX;
    STRESS = str1;
    DIS = COR.compute_Distance_Matrix();
};

void NLDR::SAMMON_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX)
{
    int dim = atoi(dim_str);
    double c = 0;
    for(int i = 0; i < size; i++)
        for(int j = 0; j < i; j++)
            c += D.matrix[i][j];

    for(int i = 0; i < size; i++)
        for(int j = 0; j < dim; j++)
        {
            gradient.matrix[i][j] = 0;
            for(int k = 0; k < size; k++)
            {
                if(ADX.matrix[i][k] > 0.00000001 && D.matrix[i][k] > 0.00000001)
                    gradient.matrix[i][j] += (ADX.matrix[i][k] - D.matrix[i][k]) / ADX.matrix[i][k] / D.matrix[i][k] * (CORX.matrix[i][j] - CORX.matrix[k][j]);
            }
            gradient.matrix[i][j] /= c;
        }
};

void NLDR::CCA()
{
	if(algorithm == (String) "MAJORIZATION")
		this->CCA_MAJORIZATION();
	else
	if(algorithm == (String) "GAUSS_SEIDEL")
		this->CCA_GAUSS_SEIDEL();
	else
	if(algorithm == (String) "STOCHASTIC")
		this->CCA_STOCHASTIC();
    else
    if(algorithm == (String) "METROPOLIS")
        this->CCA_METROPOLIS();
	else
	{
        std::cout << "Error: Invalid algorithm for CCA's STRESS!" << std::endl;
        std::cout << "For CCA's STRESS, the available algorithms are MAJORIZATION, GAUSS_SEIDEL and STOCHASTIC." << std::endl;
		exit(0);
	}
};

void NLDR::CCA_MAJORIZATION()
{
	int train_len = parameters.CCA_MAJ_max_iter;
	double lambda0 = parameters.CCA_MAJ_lambda0;
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	int times = 0;
	Matrix<double> lambda(train_len, 1);
	Matrix<double> DX(size, size);
	Matrix<double> BX(size, size);
	Matrix<double> V(size, size);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.CCA_MAJ_run_time;
	COR = X;

	for(int i = 0; i < train_len; i++)
	{
		lambda.matrix[i][0] = lambda0 * pow((double) parameters.CCA_MAJ_lambdan / (lambda0), (double) i / (train_len - 1));
	}

	stress1 = CCA_stress_function(DX, lambda.matrix[0][0]);

	CCA_compute_BX(DX, lambda.matrix[0][0], BX);

	while(err > parameters.CCA_MAJ_e && (times < train_len || !train_len) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		lambda0 = 2 * compute_lambda(DX);
		if(lambda0 > lambda.matrix[times][0])
			lambda.matrix[times][0] = lambda0;

		for(int i = 0; i < size; i++)
			V.matrix[i][i] = 0;
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < i; j++)
			{
				V.matrix[i][j] = - exp(- DX.matrix[i][j] / lambda.matrix[times][0]);
				V.matrix[j][i] = V.matrix[i][j];
				V.matrix[i][i] = V.matrix[i][i] - V.matrix[i][j];
				V.matrix[j][j] = V.matrix[j][j] - V.matrix[j][i];
			}
		}

		V = V + 1;

		if(!V.compute_inverse_matrix(V))
		{
            std::cout << "Inverse computation failed!" << std::endl;
			return;
		}

		V = V - (double) 1 / size / size;

		COR = CCA_compute_Z(BX, V);
		stress2 = CCA_stress_function(DX, lambda.matrix[times][0]);
        std::cout << "stress : " << stress2 << std::endl;
		err = fabs(((double) stress1 - stress2) / stress2);
        std::cout << "rel error : " << err << std::endl;

		CCA_compute_BX(DX, lambda.matrix[times][0], BX);

		if(stress2 > stress1)
			break;

		stress1 = stress2;
		times = times + 1;
        std::cout << "times : " << times << std::endl;
	}

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::CCA_compute_BX(const Matrix<double> &DX, double lambda, Matrix<double> &BX)
{
	double *row_sum = new double[size];
	for(int i = 0; i < size; i++)
		row_sum[i] = 0;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			if(DX.matrix[i][j] != 0)
			{
				BX.matrix[i][j] = - (double) D.matrix[i][j] / DX.matrix[i][j] * exp(- DX.matrix[i][j] / lambda);
				BX.matrix[j][i] = BX.matrix[i][j];
				row_sum[i] = row_sum[i] + BX.matrix[i][j];
				row_sum[j] = row_sum[j] + BX.matrix[j][i];
			}
		}

	for(int i = 0; i < size; i++)
		BX.matrix[i][i] = - row_sum[i];

	delete [] row_sum;
}

Matrix<double> NLDR::CCA_compute_Z(const Matrix<double> &BX, const Matrix<double> &V)
{
	int dim = atoi(dim_str);
	Matrix<double> Z(size, dim);
	double *A = new double[size * dim];
	for(int i = 0; i < size * dim; i++)
		A[i] = 0;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			A[i * dim + j] = 0;
			for(int k = 0; k < size; k++)
			{
				A[i * dim + j] = A[i * dim + j] + BX.matrix[i][k] * COR.matrix[k][j];
			}
		}

	for(int i = 0; i < size; i++)
		for(int j = 0; j < dim; j++)
		{
			Z.matrix[i][j] = 0;
			for(int k = 0; k < size; k++)
			{
				Z.matrix[i][j] = Z.matrix[i][j] + V.matrix[i][k] * A[k * dim + j];
			}
		}

	delete [] A;
	return Z;
}

double NLDR::compute_lambda(const Matrix<double> &DX)
{
	double result = 0;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			if(result < (DX.matrix[i][j] - D.matrix[i][j]))
				result = (DX.matrix[i][j] - D.matrix[i][j]) / 2;
		}
	}

	return result;
}

double NLDR::CCA_stress_function(Matrix<double> &DX, double lambda)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double stress_value = 0;
	for(int i = 0; i < size; i++)
	{
		DX.matrix[i][i] = 0;
		for(int j = 0; j < i; j++)
		{
			DX.matrix[i][j] = 0;
			DX.matrix[j][i] = 0;
			for(int k = 0; k < dim; k++)
			{
				dif = COR.matrix[i][k] - COR.matrix[j][k];
				DX.matrix[i][j] = DX.matrix[i][j] + dif * dif;
			}
			DX.matrix[i][j] = sqrt((double) DX.matrix[i][j]);
			DX.matrix[j][i] = DX.matrix[i][j];
			dif = DX.matrix[j][i] - D.matrix[j][i];
			if(D.matrix[j][i] > 0.00000001)
				stress_value = stress_value + dif * dif  * exp(- DX.matrix[i][j] / lambda);
		}
	}

	return stress_value;
}

double NLDR::CCA_stress_function_by_matrix(const Matrix<double> &DX, double lambda)
{
	double dif = 0;
	double stress_value = 0;

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			dif = D(i, j) - DX(i, j);
			stress_value = stress_value + dif * dif * exp(- DX(i, j) / lambda);
		}

	return stress_value;
}

void NLDR::CCA_GAUSS_SEIDEL()
{
	int dim = atoi(dim_str);
	int train_len = parameters.CCA_GAU_max_iter;
	int times = 0;
	double lambda0 = parameters.CCA_GAU_lambda0;
	double stress1 = 0, stress2 = 0;
	double err = 1000;
	double alpha = 0.2;
	Matrix<double> lambda(train_len * size, 1);
	Matrix<double> DX(size, size);
	Matrix<double> XD(size, dim);
	Matrix<double> XD2(size, dim);
	Matrix<double> DPJ(size, 1);
	Matrix<double> DQ(size, 1);
	Matrix<double> DR(size, 1);
	Matrix<double> term(size, 1);
	Matrix<double> term2(size, 1);
	Matrix<double> e1(dim, 1);
	Matrix<double> e2(dim, 1);
	Matrix<double> c(dim, 1);
	long start = clock();
	long cur_time = start;
	long run_time = parameters.CCA_GAU_run_time;
	COR = X;

	for(int i = 0; i < train_len; i++)
		for(int j = 0; j < size; j++)
			lambda.matrix[i * size + j][0] = lambda0 * pow((double) parameters.CCA_GAU_lambdan / (lambda0), (double) (i * size + j) / (train_len * size - 1));

	stress2 = CCA_stress_function(DX, lambda.matrix[0][0]);
    std::cout << "stress value : " << stress2 << std::endl;

	while(err > parameters.CCA_GAU_e && (times < train_len || !train_len) && (cur_time - start < run_time || !run_time))
	{
		cur_time = clock();
		stress1 = stress2;
		for(int j = 0; j < size; j++)
		{
			lambda0 = lambda.matrix[times * size + j][0];
			for(int h = 0; h < size; h++)
				for(int k = 0; k < dim; k++)
				{
					XD.matrix[h][k] = COR.matrix[j][k] - COR.matrix[h][k];
					XD2.matrix[h][k] = XD.matrix[h][k] * XD.matrix[h][k];
				}

			for(int h = 0; h < size; h++)
			{
				DPJ.matrix[h][0] = 0;
				for(int k = 0; k < dim; k++)
				{
					DPJ.matrix[h][0] = DPJ.matrix[h][0] + XD2.matrix[h][k];
				}
				DPJ.matrix[h][0] = sqrt((double) DPJ.matrix[h][0]);
			}

			for(int h = 0; h < size; h++)
			{
				DQ.matrix[h][0] = D.matrix[h][j] - DPJ.matrix[h][0];
				if(lambda0 < - DQ.matrix[h][0])
					lambda0 = - DQ.matrix[h][0];
				lambda.matrix[times * size + j][0] = lambda0;
			}

			for(int h = 0; h < size; h++)
			{
				DR.matrix[h][0] = (2 + DQ.matrix[h][0] / lambda.matrix[times * size + j][0]) * exp(-D.matrix[h][j] / lambda.matrix[times * size + j][0]);
			}

			for(int h = 0;h < size; h++)
			{
				if(fabs(DPJ.matrix[h][0]) > 0.0000001)
					term.matrix[h][0] = DQ.matrix[h][0] * DR.matrix[h][0] / DPJ.matrix[h][0];
				else
					term.matrix[h][0] = 0;
			}

			for(int h = 0; h < dim; h++)
			{
				e1.matrix[h][0] = 0;
				for(int k = 0; k < size; k++)
				{
					e1.matrix[h][0] = e1.matrix[h][0] - XD.matrix[k][h] * term.matrix[k][0];
				}
			}

			for(int h = 0; h < size; h++)
			{
				if(fabs(DPJ.matrix[h][0]) > 0.00000000001)
					term2.matrix[h][0] = ((DPJ.matrix[h][0] * DPJ.matrix[h][0] *(2 * D.matrix[h][j] + 3 * lambda.matrix[times * size + j][0] - DPJ.matrix[h][0]) - (DPJ.matrix[h][0] + lambda.matrix[times * size + j][0]) * D.matrix[h][j] * (D.matrix[h][j] + 2 * lambda.matrix[times * size + j][0])) / DPJ.matrix[h][0] / DPJ.matrix[h][0] / DPJ.matrix[h][0] / lambda.matrix[times * size + j][0] / lambda.matrix[times * size + j][0]) * exp(-D.matrix[h][j] / lambda.matrix[times * size + j][0]);
				else
					term2.matrix[h][0] = 0;
			}

			for(int h = 0; h < dim; h++)
			{
				e2.matrix[h][0] = 0;
				for(int k = 0; k < size; k++)
				{
					e2.matrix[h][0] = e2.matrix[h][0] - term.matrix[k][0] - XD2.matrix[k][h] * term2.matrix[k][0];
				}
			}

			for(int h = 0; h < dim; h++)
			{

				if(! search_first_condition(DX, stress2, - e1.matrix[h][0] / fabs(e2.matrix[h][0]), j, h, - e1.matrix[h][0] * e1.matrix[h][0] / fabs(e2.matrix[h][0]), lambda.matrix[times * size + j][0], alpha, 1, 1, stress2))
                    ;//---std::cout << "!";
				COR.matrix[j][h] = COR.matrix[j][h] - alpha * e1.matrix[h][0] / fabs(e2.matrix[h][0]);
				update_distance_matrix(DX, j, h, COR.matrix[j][h], 0);
			}
		}

		for(int i = 0; i < dim; i++)
		{
			c.matrix[i][0] = 0;
			for(int j = 0; j < size; j++)
			{
				c.matrix[i][0] = c.matrix[i][0] + COR.matrix[j][i];
			}
			c.matrix[i][0] = c.matrix[i][0] / size;
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < dim; j++)
			{
				COR.matrix[i][j] = COR.matrix[i][j] - c.matrix[j][0];
			}

		stress2 = CCA_stress_function_by_matrix(DX, lambda.matrix[times * size + size - 1][0]);
		err = fabs((stress1 - stress2) / stress2);
        std::cout << "stress value : " << stress2 << std::endl;
        std::cout << "rel error : " << err << std::endl;

		times = times + 1;
        std::cout << "times : " << times << std::endl;
	}

	STRESS = stress2;

	DIS = COR.compute_Distance_Matrix();

	return;
};

void NLDR::CCA_STOCHASTIC()
{
	int dim = atoi(dim_str);
	int epochs = parameters.CCA_STO_epochs;
	int train_len = epochs * size;
	int ind = 0;
	double alpha0 = parameters.CCA_STO_alpha0;
	double lambda0 = parameters.CCA_STO_lambda0;
	Matrix<double> alpha(train_len, 1);
	Matrix<double> lambda(train_len, 1);
	Matrix<int> sample_inds(train_len, 1);
	Matrix<double> DX(size, 1);
	Matrix<double> Y(dim, 1);
	Matrix<double> DY(size, dim);
	Matrix<double> dy(size ,1);
	Matrix<double> FY(size, 1);
    COR = X;
	for(int i = 0; i < train_len; i++)
	{
		alpha.matrix[i][0] = alpha0 * pow((double) parameters.CCA_STO_alphan, (double) i / (train_len - 1));
		lambda.matrix[i][0] = lambda0 * pow((double) parameters.CCA_STO_lambdan / (lambda0), (double) i / (train_len - 1));
	}

	for(int i = 0; i < train_len; i++)
		sample_inds.matrix[i][0] = genrand_int32() % size;

	for(int i = 0; i < train_len; i++)
	{
        ind = sample_inds.matrix[i][0];
		for(int j = 0; j < size; j++)
		{
			DX.matrix[j][0] = D.matrix[j][ind];
		}

		for(int j = 0; j < dim; j++)
		{
			Y.matrix[j][0] = COR.matrix[ind][j];
		}

		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
			{
				DY.matrix[j][k] = COR.matrix[j][k] - Y.matrix[k][0];
			}

		for(int j = 0; j < size; j++)
		{
			dy.matrix[j][0] = 0;
			for(int k = 0; k < dim; k++)
			{
				dy.matrix[j][0] = dy.matrix[j][0] + DY.matrix[j][k] * DY.matrix[j][k];
			}
			dy.matrix[j][0] = sqrt((double) dy.matrix[j][0]);
			if(fabs(dy.matrix[j][0]) < 0.0000000001)
				dy.matrix[j][0] = 1;
		}

		for(int j = 0; j < size; j++)
		{
			FY.matrix[j][0] =  exp( - dy.matrix[j][0] / lambda.matrix[i][0]) * (DX.matrix[j][0] / dy.matrix[j][0] - 1);
		}

		//--- update ---
		for(int j = 0; j < size; j++)
			for(int k = 0; k < dim; k++)
			{
				COR.matrix[j][k] = COR.matrix[j][k] + alpha.matrix[i][0] * FY.matrix[j][0] * DY.matrix[j][k];
			}
		if(i % size == 0)
            std::cout << "epochs : " << (int) i / size << std::endl;
	}

	DIS = COR.compute_Distance_Matrix();

	STRESS = CCA_stress_function_by_matrix(DIS, lambda.matrix[train_len - 1][0]);

    std::cout << "stress value : " << STRESS << std::endl;

	return;
};

void NLDR::CCA_METROPOLIS()
{
    int dim = atoi(dim_str);
    double T = 10000;
    double delta = 0;
    int train_len = 2000;
    Matrix<double> lambda(train_len, 1);
    double lambda0 = 12000;
    double lambdan = 1;
    Matrix<double> gradient(size, dim);
    Matrix<double> candidate(size, dim);
    Matrix<double> DX(size, size);
    Matrix<double> ADX(size, size);
    Matrix<double> CORX(size, dim);
    double str1 = 0;
    double str2 = 0;
    CORX = X;
    COR = CORX;

    double weight = 0;
    for(int i = 0; i < size; i++)
        for(int j = 0; j < i; j++)
            weight += D.matrix[i][j] * D.matrix[i][j];

    for(int i = 0; i < train_len; i++)
        lambda.matrix[i][0] = lambda0 * pow(lambdan / lambda0, (double) i / (train_len - 1));

    str1 = CCA_stress_function(DX, lambda.matrix[0][0]) / weight;
    std::cout << "start : accept stress : " << str1 << std::endl;
    ADX = DX;
    CCA_GRADIENT(gradient, CORX, ADX, lambda.matrix[0][0]);
    delta = gradient.norm(M_FRO_NORM) / weight;
    std::cout << "gradient norm:" << delta << std::endl;
    delta = (double) 1 / (delta * delta);
    std::cout << "delta:" << delta << std::endl;

    for(int i = 0; i < train_len; i++)
    {
        lambda0 = compute_lambda(ADX) * 2;
        if(lambda0 > lambda.matrix[i][0])
            lambda.matrix[i][0] = lambda0;
        COR = CORX;
        str1 = CCA_stress_function_by_matrix(ADX, lambda.matrix[i][0]) / weight;
        CCA_GRADIENT(gradient, CORX, ADX, lambda.matrix[i][0]);
        for(int j = 0; j < size; j++)
            for(int k = 0; k < dim; k++)
            {
                candidate.matrix[j][k] = CORX.matrix[j][k] - delta * gradient.matrix[j][k] / weight + sqrt(T * delta / size / dim) * Normal_random_generator();
            }
        COR = candidate;
        str2 = CCA_stress_function(DX, lambda.matrix[i][0]) / weight;
        if(genrand_real2() < exp((str1 - str2) / (T)))
        {
            CORX = candidate;
            ADX = DX;
            std::cout << i << " : accept stress : " << str2 << std::endl;
        }

        T = T * 0.99;
    }
    COR = CORX;
    STRESS = str2;
    DIS = COR.compute_Distance_Matrix();
};

void NLDR::CCA_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX, double lambda)
{
    int dim = atoi(dim_str);

    for(int i = 0; i < size; i++)
        for(int j = 0; j < dim; j++)
        {
            gradient.matrix[i][j] = 0;
            for(int k = 0; k < size; k++)
            {
                if(ADX.matrix[i][k] > 0.00000001 && D.matrix[i][k] > 0.00000001)
                    gradient.matrix[i][j] += (ADX.matrix[i][k] - D.matrix[i][k]) / ADX.matrix[i][k] * (CORX.matrix[i][j] - CORX.matrix[k][j]) * (2 + (D.matrix[i][k] - ADX.matrix[i][k]) / lambda) * exp(- ADX.matrix[i][k] / lambda);
            }
        }
};

bool NLDR::search_first_condition(const Matrix<double> &DX, double fxc, const Matrix<double> &dir, int i, int j, double initslope, double para, double &step_size, double lambda, double pre_step_size, double &fxnew)
{
	int dim = atoi(dim_str);
	double minstepsize = 0.01;
	double alpha = 0.1;
	double sizetemp = 0;
	double a = 0, b = 0;
	double disc = 0;
	double fxnew_pre = fxnew;
	step_size = lambda;
	Matrix<double> xc(dim, 1);
	Matrix<double> xcnew(dim, 1);

	for(int k = 0; k < dim; k++)
		xc(k) = COR(j, k);

	while(true)
	{
		for(int k = 0; k < dim; k++)
			xcnew(k) = xc(k) + step_size * dir(k);

		fxnew = update_Ei_cost_function(fxc, DX, i, j, para, xcnew);

		if(fxnew <= fxc + alpha * step_size * initslope)
		{
			return true;
		} else
		if(step_size < minstepsize)
		{
			step_size = minstepsize;
			return false;
		} else
		{
			if(step_size == 1)
			{
				sizetemp = - initslope / (2 * (fxnew - fxc - initslope));
			} else
			{
				a = ((fxnew - fxc - step_size * initslope) / step_size / step_size - (fxnew_pre - fxc - pre_step_size * initslope) / pre_step_size / pre_step_size) / (step_size - pre_step_size);
				b = (- pre_step_size * (fxnew - fxc - step_size * initslope) / step_size / step_size + step_size * (fxnew_pre - fxc - pre_step_size * initslope) / pre_step_size / pre_step_size) / (step_size - pre_step_size);

				disc = b * b - 3 * a * initslope;
				if(fabs(a) < 0.0000001)
					sizetemp = - initslope / (2 * b);
				else
					sizetemp = (-b + sqrt(disc)) / (3 * a);

				if(sizetemp > 0.5 * step_size)
					sizetemp = 0.5 * step_size;
			}

			pre_step_size = step_size;
			fxnew_pre = fxnew;
			if(sizetemp <= 0.1 * step_size)
				step_size = 0.1 * step_size;
			else
				step_size = sizetemp;
		}
	}

	return true;
}

bool NLDR::search_first_condition(const Matrix<double> &DX, double fxc, double dir, int i, int j, double initslope, double para, double &step_size, double lambda, double pre_step_size, double &fxnew)
{
	double minstepsize = 0.01;
	double xc = COR(i, j);
	double xcnew = 0;
	double alpha = 0.00001;
	double sizetemp = 0;
	double a = 0;
	double b = 0;
	double disc = 0;
	double fxnew_pre = fxnew;
	step_size = lambda;
	while(true)
	{
		xcnew = xc + step_size * dir;

		fxnew = this->update_cost_function(DX, fxc, i, j, xcnew, para);

		if(fxnew <= fxc + alpha * step_size * initslope)
		{
			return true;
		} else
		if(step_size < minstepsize)
		{
			step_size = minstepsize;
			return false;
		} else
		{
			if(step_size == 1)
			{
				sizetemp = - initslope / (2 * (fxnew - fxc - initslope));
			} else
			{
				a = ((fxnew - fxc - step_size * initslope) / step_size / step_size - (fxnew_pre - fxc - pre_step_size * initslope) / pre_step_size / pre_step_size) / (step_size - pre_step_size);
				b = (- pre_step_size * (fxnew - fxc - step_size * initslope) / step_size / step_size + step_size * (fxnew_pre - fxc - pre_step_size * initslope) / pre_step_size / pre_step_size) / (step_size - pre_step_size);

				disc = b * b - 3 * a * initslope;
				if(fabs(a) < 0.0000001)
					sizetemp = - initslope / (2 * b);
				else
					sizetemp = (-b + sqrt(disc)) / (3 * a);

				if(sizetemp > 0.5 * step_size)
					sizetemp = 0.5 * step_size;
			}

			pre_step_size = step_size;
			fxnew_pre = fxnew;
			if(sizetemp <= 0.1 * step_size)
				step_size = 0.1 * step_size;
			else
				step_size = sizetemp;
		}
	}

	return true;
}

void NLDR::update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew, double para)
{
	if(cost_function == (String) "KRUSKAL1")
		return KRUSKAL1_update_distance_matrix(DX, i, j, xcnew);

	if(cost_function == (String) "NORMALIZED")
		return NORMALIZED_update_distance_matrix(DX, i, j, xcnew);

	if(cost_function == (String) "SAMMON")
		return SAMMON_update_distance_matrix(DX, i, j, xcnew);

	if(cost_function == (String) "CCA")
		return CCA_update_distance_matrix(DX, i, j, xcnew, para);

    std::cout << "Warning: This cost function does not exist in this program!" << std::endl;
};

void NLDR::KRUSKAL1_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew)
{
	double eta = para1;
	double rho = para3;
	int dim = atoi(dim_str);
	double dif = 0;
	for(int h = 0; h < size; h++)
	{
		eta -= DX.matrix[i][h] * DX.matrix[i][h];
		rho -= DX.matrix[i][h] * D.matrix[i][h];
		DX.matrix[i][h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			DX.matrix[i][h] = DX.matrix[i][h] + dif * dif;
		}
		DX.matrix[i][h] = sqrt(DX.matrix[i][h]);
		DX.matrix[h][i] = DX.matrix[i][h];
		eta += DX.matrix[i][h] * DX.matrix[i][h];
		rho += DX.matrix[i][h] * D.matrix[i][h];
	}
	DX.matrix[i][i] = 0;
	para1 = eta;
	para3 = rho;
}

void NLDR::NORMALIZED_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew)
{
	int dim = atoi(dim_str);
	double dif = 0;
	for(int h = 0; h < size; h++)
	{
		DX.matrix[i][h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			DX.matrix[i][h] = DX.matrix[i][h] + dif * dif;
		}
		DX.matrix[i][h] = sqrt(DX.matrix[i][h]);
		DX.matrix[h][i] = DX.matrix[i][h];
	}
	DX.matrix[i][i] = 0;
};

void NLDR::SAMMON_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew)
{
	int dim = atoi(dim_str);
	double dif = 0;
	for(int h = 0; h < size; h++)
	{
		DX.matrix[i][h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			DX.matrix[i][h] = DX.matrix[i][h] + dif * dif;
		}
		DX.matrix[i][h] = sqrt(DX.matrix[i][h]);
		DX.matrix[h][i] = DX.matrix[i][h];
	}
	DX.matrix[i][i] = 0;
};

void NLDR::CCA_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew, double para)
{
	int dim = atoi(dim_str);
	double dif = 0;
	for(int h = 0; h < size; h++)
	{
		DX.matrix[i][h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			DX.matrix[i][h] = DX.matrix[i][h] + dif * dif;
		}
		DX.matrix[i][h] = sqrt(DX.matrix[i][h]);
		DX.matrix[h][i] = DX.matrix[i][h];
	}
	DX.matrix[i][i] = 0;
};

double NLDR::update_cost_function(const Matrix<double> &DX, double fxc, int i, int j, double xcnew, double para)
{
	if(cost_function == (String) "KRUSKAL1")
		return KRUSKAL1_update_stress_function(DX, i, j, xcnew, fxc);

	if(cost_function == (String) "NORMALIZED")
        return NORMALIZED_update_stress_function(DX, i, j, xcnew, fxc, para);

	if(cost_function == (String) "SAMMON")
        return SAMMON_update_stress_function(DX, i, j, xcnew, fxc, para);

	if(cost_function == (String) "CCA")
		return CCA_update_stress_function(DX, i, j, xcnew, fxc, para);

    std::cout << "Warning: This cost function does not exist in this program!" << std::endl;
	return 0;
};

double NLDR::KRUSKAL1_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc)
{
	double eta = para1;
	double eta_d = para2;
	double rho = para3;
	int dim = atoi(dim_str);
	double *di_ = new double[size];
	double dif = 0;

	for(int h = 0; h < size; h++)
	{
		eta -= DX.matrix[i][h] * DX.matrix[i][h];
		rho -= DX.matrix[i][h] * D.matrix[i][h];
		di_[h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			di_[h] = di_[h] + dif * dif;
		}
		di_[h] = sqrt(di_[h]);
		eta += di_[h] * di_[h];
		rho += di_[h] * D.matrix[i][h];
	}

	fxc = (eta + eta_d - 2 * rho) / eta;

	delete [] di_;

	return fxc;
};

double NLDR::NORMALIZED_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double weight)
{
	int dim = atoi(dim_str);
	double *di_ = new double[size];
	double dif = 0;

	for(int h = 0; h < size; h++)
	{
		di_[h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			di_[h] = di_[h] + dif * dif;
		}
		di_[h] = sqrt(di_[h]);
	}

	for(int h = 0; h < size; h++)
	{
		if(h != i)
            fxc = fxc + (2 * D.matrix[i][h] - di_[h] - DX.matrix[i][h]) * (DX.matrix[i][h] - di_[h]) / weight;
	}

	delete [] di_;

	return fxc;
};

double NLDR::SAMMON_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double weight)
{
	int dim = atoi(dim_str);
	double *di_ = new double[size];
	double dif = 0;

	for(int h = 0; h < size; h++)
	{
		di_[h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			di_[h] = di_[h] + dif * dif;
		}
		di_[h] = sqrt(di_[h]);
	}

	for(int h = 0; h < size; h++)
	{
		if(h != i && D.matrix[i][h] >= 0.0000000001)
            fxc = fxc + (2 * D.matrix[i][h] - di_[h] - DX.matrix[i][h]) * (DX.matrix[i][h] - di_[h]) / D.matrix[i][h] / weight;
	}

	delete [] di_;

	return fxc;
};

double NLDR::CCA_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double para)
{
	int dim = atoi(dim_str);
	double *di_ = new double[size];
	double dif = 0;

	for(int h = 0; h < size; h++)
	{
		di_[h] = 0;
		for(int k = 0; k < dim; k++)
		{
			if(k != j)
				dif = COR.matrix[i][k] - COR.matrix[h][k];
			else
				dif = xcnew - COR.matrix[h][k];

			di_[h] = di_[h] + dif * dif;
		}
		di_[h] = sqrt(di_[h]);
	}

	for(int h = 0; h < size; h++)
	{
		if(h != i)
		{
			dif = D.matrix[i][h] - di_[h];
			fxc += (dif * dif) * exp(- di_[h] / para);
			dif = D.matrix[i][h] - DX.matrix[i][h];
			fxc -= dif * dif * exp(- DX.matrix[i][h] / para);
		}
	}

	delete [] di_;

	return fxc;
};

double NLDR::update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, double para, const Matrix<double> &xcnew)
{
	if(cost_function == (String) "NORMALIZED_Ei")
		return NORMALIZED_update_Ei_cost_function(fxc, DX, ind, j, xcnew);

	if(cost_function == (String) "SAMMON_Ei")
		return SAMMON_update_Ei_cost_function(fxc, DX, ind, j, xcnew);

	if(cost_function == (String) "CCA_Ei")
		return CCA_update_Ei_cost_function(fxc, DX, ind, j, para, xcnew);

    std::cout << "Warning: No such cost function" << std::endl;
	return 0;
}

double NLDR::NORMALIZED_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, const Matrix<double> &xcnew)
{
	int dim = atoi(dim_str);
	double Ei = fxc;
	double dijnew = 0;
	double dif = 0;
	for(int k = 0; k < dim; k++)
	{
		dif = COR(ind, k) - xcnew(k);
		dijnew += dif * dif;
	}
	dijnew = sqrt(dijnew);

	Ei += (2 * D(ind, j) - dijnew - DX(ind, j)) * (DX(ind, j) - dijnew) / 2;
	return Ei;
}

double NLDR::SAMMON_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, const Matrix<double> &xcnew)
{
	int dim = atoi(dim_str);
	double Ei = fxc;
	double dijnew = 0;
	double dif = 0;
	for(int k = 0; k < dim; k++)
	{
		dif = COR(ind, k) - xcnew(k);
		dijnew += dif * dif;
	}
	dijnew = sqrt(dijnew);

	Ei += (2 * D(ind, j) - dijnew - DX(ind, j)) * (DX(ind, j) - dijnew) / 2 / D(ind, j);
	return Ei;
}

double NLDR::CCA_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, double para, const Matrix<double> &xcnew)
{
	int dim = atoi(dim_str);
	double dif = 0;
	double Ei = fxc;
	double dijnew = 0;
	for(int k = 0; k < dim; k++)
	{
		dif = COR(ind, k) - xcnew(k);
		dijnew += dif * dif;
	}
	dijnew = sqrt(dijnew);

	dif = D(ind, j) - dijnew;
	Ei += (dif * dif) * exp(- dijnew / para);
	dif = D(ind, j) - DX(ind, j);
	Ei -= dif * dif * exp(- DX(ind, j) / para);

	return Ei;
}

void NLDR::Trustworthiness_analysis()
{
	int row = Trustworthiness.get_row();
	int n = (int) Trustworthiness(row - 1, 0) + 1;
	double min_j = 0;
	int min_jk = 0;
	int order = 0;
	Matrix<double> MININ(n, 2);
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
			if(DIS.matrix[i][j] < ZERO)
			{
				MININ.matrix[0][0] = j;
				break;
			}

		// find j-th nearest point for point i for imbeding D.
		for(int j = 1; j < n; j++)
		{
			min_j = 10000;
			for(int k = 0; k < size; k++)
			{
				if(MININ.matrix[j - 1][1] == DIS.matrix[i][k] && k > MININ.matrix[j - 1][0] && DIS.matrix[i][k] < min_j)
				{
					min_j = DIS.matrix[i][k];
					min_jk = k;
					break;
				}
				if(MININ.matrix[j - 1][1] < DIS.matrix[i][k] && DIS.matrix[i][k] < min_j)
				{
					min_j = DIS.matrix[i][k];
					min_jk = k;
				}
			}
			MININ.matrix[j][0] = min_jk;
			MININ.matrix[j][1] = min_j;
		}

		// find each imbeding neighbor's order in original D.
		for(int j = 1; j < n; j++)
		{
			order = 0;
			for(int k = 0; k < size; k++)
			{
				if(D.matrix[i][k] < D.matrix[i][(int) MININ.matrix[j][0]])
					order++;
			}

			for(int k = 0; k < Trustworthiness.get_row(); k++)
			{
				if(j <= Trustworthiness.matrix[k][0] && order > Trustworthiness.matrix[k][0])
					Trustworthiness.matrix[k][1] += order - Trustworthiness.matrix[k][0];
			}
		}
	}
	for(int i = 0; i < row; i++)
		Trustworthiness.matrix[i][1] = 1 - (double) 2 * Trustworthiness.matrix[i][1] / (size * Trustworthiness.matrix[i][0] * (2 * size - 3 * Trustworthiness.matrix[i][0] - 1));
};

void NLDR::Continuity_analysis()
{
	int row = Continuity.get_row();
	int n = (int) Continuity(row - 1, 0) + 1;
	double min_j = 0;
	int min_jk = 0;
	int order = 0;
	Matrix<double> MININ(n, 2);
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
			if(DIS.matrix[i][j] < ZERO)
			{
				MININ.matrix[0][0] = j;
				break;
			}

		// find j-th nearest point for point i for imbeding D.
		for(int j = 1; j < n; j++)
		{
			min_j = 10000;
			for(int k = 0; k < size; k++)
			{
				if(MININ.matrix[j - 1][1] == D.matrix[i][k] && k > MININ.matrix[j - 1][0] && D.matrix[i][k] < min_j)
				{
					min_j = D.matrix[i][k];
					min_jk = k;
					break;
				}
				if(MININ.matrix[j - 1][1] < D.matrix[i][k] && D.matrix[i][k] < min_j)
				{
					min_j = D.matrix[i][k];
					min_jk = k;
				}
			}
			MININ.matrix[j][0] = min_jk;
			MININ.matrix[j][1] = min_j;
		}

		// find each imbeding neighbor's order in original D.
		for(int j = 1; j < n; j++)
		{
			order = 0;
			for(int k = 0; k < size; k++)
			{
				if(DIS.matrix[i][k] < DIS.matrix[i][(int) MININ.matrix[j][0]])
					order++;
			}
			for(int k = 0; k < Continuity.get_row(); k++)
			{
				if(j <= Continuity.matrix[k][0] && order > Continuity.matrix[k][0])
					Continuity.matrix[k][1] += order - Continuity.matrix[k][0];
			}
		}
	}

	for(int i = 0; i < row; i++)
		Continuity.matrix[i][1] = 1 - (double) 2 * Continuity.matrix[i][1] / (size * Continuity.matrix[i][0] * (2 * size - 3 * Continuity.matrix[i][0] - 1));
};

void NLDR::oneNN_analysis()
{
	oneNN.resize(2, 1);
	int sal_index[15] = {119, 225, 381, 1164, 1270, 1466, 1615, 1737, 1848, 1959, 2314, 2422, 2800, 2903, 3011};
	int Mam_index[15] = {219, 365, 905, 1267, 1495, 1928, 2482, 2677, 2847, 2976, 4535, 4685, 5741, 5855, 6001};
	int Ather_index[15] = {256, 461, 876, 1815, 2201, 2645, 3288, 3523, 4030, 4401, 5091, 5310, 6672, 6860, 7022};
	int salshu_index[15] = {107, 210, 321, 1089, 1191, 1312, 1426, 1533, 1640, 1745, 1981, 2085, 2356, 2458, 2559};
	int *index = NULL;
	double min = 0;
	int min_j = 0;
	int i_index = 0;
	int j_index = 0;
	int correct_num = 0;

	if(D_prefname == (String) "Sal_allBStrees")
		index = sal_index;
	else
	if(D_prefname == (String) "Mam_allBStrees")
		index = Mam_index;
	else
	if(D_prefname == (String) "Ather_allBStrees")
		index = Ather_index;
	else
	if(D_prefname == (String) "Sal_allBStrees_shuffled")
		index = salshu_index;
	else
		return;

	i_index = 0;
	correct_num = 0;
	for(int i = 0; i < size; i++)
	{
		if(i >= index[i_index])
			i_index++;
		min_j = (i == 0) ? 1 : 0;
		min = DIS.matrix[i][min_j];
		for(int j = 0; j < size; j++)
		{
			if(DIS.matrix[i][j] < min && i != j)
			{
				min = DIS.matrix[i][j];
				min_j = j;
			}
		}

		for(int j = 0; j < 15; j++)
		{
			if(min_j < index[j])
			{
				j_index = j;
				break;
			}
		}
		if(j_index == i_index)
			correct_num++;
	}
	oneNN.matrix[1][0] = (double) correct_num / size;

	i_index = 0;
	correct_num = 0;
	for(int i = 0; i < size; i++)
	{
		if(i >= index[i_index])
			i_index++;
		min_j = (i == 0) ? 1 : 0;
		min = D.matrix[i][min_j];
		for(int j = 0; j < size; j++)
		{
			if(D.matrix[i][j] < min && i != j)
			{
				min = D.matrix[i][j];
				min_j = j;
			}
		}

		for(int j = 0; j < 15; j++)
		{
			if(min_j < index[j])
			{
				j_index = j;
				break;
			}
		}
		if(j_index == i_index)
			correct_num++;
	}

	oneNN.matrix[0][0] = (double) correct_num / size;
};

#endif
