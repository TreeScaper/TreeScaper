
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

// wNLDR.h
//  cost function: Linear MDS, Metric MDS, CCA, SAMMON
//  algorithm    : Gauss Seidel, Majorization, stachatic gradient.
//          March/11/2010
//                    by whuang

#ifndef WNLDR_H
#define WNLDR_H

#include "wdef.h"
#include "warray.cpp"
#include "wstring.h"
#include "wmix.h"
#include "wimport_form.h"
#include "wmapping.cpp"
#include "randgen.h"
#include "wfile.h"
#include "wmatrix.cpp"

struct nldr_parameters{
    int distance_file_type; //
	int length_tru; // 11
	int interval_tru; // 11
	int length_con; // 5
	int iterval_con; // 5
	double random_start; // -100
	double random_end; // 100

	double KRU_LIN_e; // 0.00001
	int KRU_LIN_max_iter; // 1000
	long KRU_LIN_run_time; // 0
	double KRU_MAJ_e; // 0.00001
	int KRU_MAJ_max_iter; // 1000
	long KRU_MAJ_run_time; // 0
	double KRU_GAU_e; // 0.0001
	int KRU_GAU_max_iter; // 1000
	long KRU_GAU_run_time; // 0
	int KRU_STO_epochs; // 100
	double KRU_STO_alpha0; // 0.5
	double KRU_STO_alphan; // 0.01

	double NOR_MAJ_e; // 0.00001
	int NOR_MAJ_max_iter; // 1000
	long NOR_MAJ_run_time; // 0
	double NOR_GAU_e; // 0.0001
	int NOR_GAU_max_iter; // 1000
	long NOR_GAU_run_time; // 0
	int NOR_STO_epochs; // 100
	double NOR_STO_alpha0; // 0.25
	double NOR_STO_alphan; // 0.01

	double NLM_MAJ_e; // 0.0001
	int NLM_MAJ_max_iter; // 1000
	long NLM_MAJ_run_time; // 0
	double NLM_GAU_e; // 0.0001
	int NLM_GAU_max_iter; // 1000
	long NLM_GAU_run_time; // 0
	int NLM_STO_epochs; // 100
	double NLM_STO_alpha0; // 0.25
	double NLM_STO_alphan; // 0.01

	double CCA_MAJ_e; // 0.0001
	int CCA_MAJ_max_iter; // 1000
	long CCA_MAJ_run_time; // 0
	double CCA_MAJ_lambda0; // 12000
	double CCA_MAJ_lambdan; // 1
	double CCA_GAU_e; // 0.0001
	int CCA_GAU_max_iter; // 1000
	long CCA_GAU_run_time; // 0
	double CCA_GAU_lambda0; // 12000
	double CCA_GAU_lambdan; // 1
	int CCA_STO_epochs; // 100
	double CCA_STO_lambda0; // 12000
	double CCA_STO_lambdan; // 1
	double CCA_STO_alpha0; // 0.5
	double CCA_STO_alphan; // 0.01
};

class NLDR{
public:
	NLDR(){};
	~NLDR(){};

	NLDR(String fname, String dim, String cost, String algo = (String) "", String init_md = (String) "CLASSIC_MDS", String flag = (String) "", long seed = -1, String para_fname = "")
	{
        double **dist = NULL;
        init_NLDR(fname, dist, 0, dim, cost, algo, init_md, flag, seed, para_fname);
	}

  void init_NLDR(String fname, double **dist, int size, String dim, String cost, String algo = (String) "", String init_md = (String) "CLASSIC_MDS", String flag = (String) "", long seed = -1, String para_fname = "");

#ifdef COMMAND_LINE_VERSION
	NLDR(String fname, String ftype, String dim, String cost, String algo = (String) "", String init_md = (String) "CLASSIC_MDS", String flag = (String) "", long seed = -1, String para_fname = "")
	{
		init_NLDR(fname, ftype, dim, cost, algo, init_md, flag, seed, para_fname);
	}

	void init_NLDR(String fname, String ftype, String dim, String cost, String algo = (String) "", String init_md = (String) "CLASSIC_MDS", String flag = (String) "", long seed = -1, String para_fname = "");
#endif

	void Compute_NLDR();

	void result_analysis();

	void output_to_files();

    nldr_parameters parameters;

private:
//functions:

#ifdef COMMAND_LINE_VERSION
	void NLDR_init_parameters(String para_filename);

	void NLDR_load_MX();
#endif

	void make_output_file_names(String &filename_COR, String &filename_DIS, String &filename_STR, String &filename_TIM, String &filename_TRU, String &filename_CON, String &filename_1NN);

	void NLDR_load_D();

	void NLDR_init_X(String init_md, long seed);

	String make_filename(String name_D, String dimension, String cost_f, String output, String algorithm, String flag);

//	void NLDR_init_parameters(String para_filename);

	void CLASSIC_MDS();

	double CLASSIC_MDS_stress_function(const Matrix<double> &S1, const Matrix<double> &S2);

	void KRUSKAL1();
	
	void KRUSKAL1_LINEAR_ITERATION();
	
	void KRUSKAL1_MAJORIZATION();

	Matrix<double> KRUSKAL1_compute_Z(const Matrix<double> &BX);

	double KRUSKAL1_stress_function(Matrix<double> &DX);

	double KRUSKAL1_rescale(Matrix<double> &DX);

	double KRUSKAL1_stress_function_by_matrix(const Matrix<double> &DX);

	void KRUSKAL1_compute_BX(const Matrix<double> &DX, Matrix<double> &BX);

	void KRUSKAL1_GAUSS_SEIDEL();

	void KRUSKAL1_STOCHASTIC();

    void KRUSKAL1_METROPOLIS();

    void KRUSKAL1_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX);

	void NORMALIZED();

	void NORMALIZED_MAJORIZATION();

	Matrix<double> NORMALIZED_compute_Z(const Matrix<double> &BX);

	double NORMALIZED_stress_function(Matrix<double> &DX);

	double NORMALIZED_stress_function_by_matrix(const Matrix<double> &DX);

	void NORMALIZED_compute_BX(const Matrix<double> &DX, Matrix<double> &BX);

	void NORMALIZED_GAUSS_SEIDEL();

	void NORMALIZED_STOCHASTIC();

    void NORMALIZED_METROPOLIS();

    void NORMALIZED_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX);

    double Normal_random_generator();

	void SAMMON();

	void SAMMON_MAJORIZATION();

	Matrix<double> SAMMON_compute_Z(const Matrix<double> &BX, const Matrix<double> &V);
	
	double SAMMON_stress_function(Matrix<double> &DX);
	
	double SAMMON_stress_function_by_matrix(const Matrix<double> &DX);
	
	void SAMMON_compute_BX(const Matrix<double> &DX, Matrix<double> &BX);

	void SAMMON_GAUSS_SEIDEL();

	void SAMMON_STOCHASTIC();

    void SAMMON_METROPOLIS();

    void SAMMON_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX);

	void CCA();

	void CCA_MAJORIZATION();

	void CCA_compute_BX(const Matrix<double> &DX, double lambda, Matrix<double> &BX);

	Matrix<double> CCA_compute_Z(const Matrix<double> &BX, const Matrix<double> &V);

	double compute_lambda(const Matrix<double> &DX);

	double CCA_stress_function(Matrix<double> &DX, double lambda);

	double CCA_stress_function_by_matrix(const Matrix<double> &DX, double lambda);

	void CCA_GAUSS_SEIDEL();

	void CCA_STOCHASTIC();

    void CCA_METROPOLIS();

    void CCA_GRADIENT(Matrix<double> &gradient, Matrix<double> &CORX, Matrix<double> &ADX, double lambda);

	bool search_first_condition(const Matrix<double> &DX, double fxc, const Matrix<double> &dir, int i, int j, double initslope, double para, double &step_size, double lambda, double pre_step_size, double &fxnew);

	bool search_first_condition(const Matrix<double> &DX, double fxc, double dir, int i, int j, double initslope, double para, double &step_size, double lambda, double pre_step_size, double &fxnew);

	void update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew, double para);

	void KRUSKAL1_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew);

	void NORMALIZED_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew);

	void SAMMON_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew);

	void CCA_update_distance_matrix(Matrix<double> &DX, int i, int j, double xcnew, double para);

	double update_cost_function(const Matrix<double> &DX, double fxc, int i, int j, double xcnew, double para);

	double KRUSKAL1_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc);

    double NORMALIZED_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double weight);

    double SAMMON_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double weight);

	double CCA_update_stress_function(const Matrix<double> &DX, int i, int j, double xcnew, double fxc, double para);

	double update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, double para, const Matrix<double> &xcnew);

	double NORMALIZED_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, const Matrix<double> &xcnew);

	double SAMMON_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, const Matrix<double> &xcnew);

	double CCA_update_Ei_cost_function(double fxc, const Matrix<double> &DX, int ind, int j, double para, const Matrix<double> &xcnew);

	void Trustworthiness_analysis();

	void Continuity_analysis();

	void oneNN_analysis();

//storage:

//initial information:
	Matrix<double> D;
	Matrix<double> MX;
	Matrix<double> X;
	int size;
	String D_prefname;
	String D_postfname;
	String dim_str;
	String cost_function;
	String algorithm;
	String file_flag;

//output information:
	Matrix<double> U;
	Matrix<double> S;
	Matrix<double> Vt;
	Matrix<double> COR;
	Matrix<double> DIS;
	long time_cost;
	double STRESS;
	// result analysis
	Matrix<double> Trustworthiness;
	Matrix<double> Continuity;
	Matrix<double> oneNN;
	double para1;
	double para2;
    double para3;
};

#endif
