
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

// dimension estimator file
//            April/20/2010
//                    whuang

#ifndef WDIMEST_CPP
#define WDIMEST_CPP

#include "wDimEst.h"

void DimEst::initial_DimEst(String filename, double **dist, int sizeinput, int diminput, String meth, String dataform, String para_fname)
{
    method = meth;
    num_r = parameters.cor_n;// only for this interface version
    File D_file(filename);
    D_prefname = filename;//D_file.prefix_name();
    D_postfname = filename;//D_file.postfix_name();

    // get distance matrix
    if(dataform == (String) "DIS")
    {
        if(dist == NULL && sizeinput == 0)
        {
            if(!D_file.is_open())
            {
                cout << "Error: File \"" << filename << "\" cannot be open! Please check if this file exists or is readable." << endl;
                exit(0);
            }

            size = D_file.lines();
            D.resize(size, size);
            X.resize(num_r, 1);
            D_file.seek(0);

            if(parameters.distance_file_type == 0)
            {
                for(int i = 0; i < size; i++)
                    for(int j = 0; j <= i; j++)
                    {
                        D_file >> D.matrix[i][j];
                        D.matrix[j][i] = D.matrix[i][j];
                    }
            }
            else if(parameters.distance_file_type == 1)
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
                cout << "Warning: This is not a type of distance file format!" << endl;
            }
        }
        else if(dist != NULL && sizeinput > 0)
        {
            size = sizeinput;
            D.resize(size, size);
            X.resize(num_r, 1);

            for(int i = 0; i < size; i++)
                for(int j = 0; j < size; j++)
                    D.matrix[i][j] = (double) dist[i][j];
        }
        else
        {
            cout << "Error: Incorrect input data parameters" << endl;
            exit(0);
        }
    }
    else if(dataform == (String) "COR")
    {
        if(dist == NULL && sizeinput == 0 && diminput == 0)
        {
            if(!D_file.is_open())
            {
                cout << "Error: File \"" << filename << "\" cannot be open! Please check if this file exists or is readable." << endl;
                exit(0);
            }

            size = D_file.lines();
            D.resize(size, size);
            X.resize(num_r, 1);
            D_file.seek(0);

            int dim = D_file.cols();
            Matrix<double> X(size, dim);
            for(int i = 0; i < size; i++)
                for(int j = 0; j < dim; j++)
                    D_file >> X.matrix[i][j];
            D = X.compute_Distance_Matrix();

        }
        else if(dist != NULL && sizeinput > 0 && diminput > 0)
        {
            size = sizeinput;
            int dim = diminput;
            D.resize(size, size);
            X.resize(num_r, 1);

            Matrix<double> X(size, dim);
            for(int i = 0; i < size; i++)
                for(int j = 0; j < dim; j++)
                    X.matrix[i][j] = (double) dist[i][j];
            D = X.compute_Distance_Matrix();
        }
        else
        {
            cout << "Error: Incorrect input data parameters" << endl;
            exit(0);
        }
    }
    else
    {
        cout << "Warning: Undefined matrix type!\n\n";
        exit(0);
    }
};

#ifdef COMMAND_LINE_VERSION
void DimEst::initial_DimEst(String filename, String meth, String dataform, String para_fname)
{
    method = meth;
    init_parameters(para_fname); //for command version
    File D_file(filename);
    if(!D_file.is_open())
    {
        cout << "Error: File \"" << filename << "\" cannot be opened! Please check if this file exists or is readable." << endl;
        exit(0);
    }
    size = D_file.lines();
    D_prefname = D_file.prefix_name_lastof();
    D_postfname = D_file.postfix_name_lastof();
    D.resize(size, size);
    X.resize(num_r, 1);
    D_file.seek(0);

    // get distance matrix
    if(dataform == (String) "DIS")
    {
        if(parameters.distance_file_type == 0)
        {
            for(int i = 0; i < size; i++)
                for(int j = 0; j <= i; j++)
                {
                    D_file >> D.matrix[i][j];
                    D.matrix[j][i] = D.matrix[i][j];
                }
        } else if(parameters.distance_file_type == 1)
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
            cout << "Warning: This is not a type of distance file format!" << endl;
        }
    } else if(dataform == (String) "COR")
    {
        int dim = D_file.cols();
        Matrix<double> X(size, dim);
        for(int i = 0; i < size; i++)
            for(int j = 0; j < dim; j++)
                D_file >> X.matrix[i][j];
        D = X.compute_Distance_Matrix();
    }
    else
    {
        cout << "Warning: Undefined matrix type!\n\n";
        exit(0);
    }
};

void DimEst::init_parameters(String para_filename)
{
	parameters.cor_n = 200;
	parameters.mle_n = 200;
	parameters.nn_n = 200;
    parameters.distance_file_type = 0;

	if(para_filename == (String) "")
	{
		num_r = 200;
		return;
	}

	Array<Mapping<Mix, Mix> > para_arr;
    Mapping<Mix, Mix> para;
    para_arr = LPIMPORT_FORM->Read(para_filename);

	int num = para_arr.get_length();
    Mix Key, Value;
	for(int i = 0; i < num; i++)
    {
        Key = para_arr[i][(String) "name"];
        Value = para_arr[i][(String) "value"];
        para[Key] = Value;;
    }

	if(para.Include_Key((String) "cor_n"))
		parameters.cor_n = para[(String) "cor_n"];
	if(para.Include_Key((String) "mle_n"))
		parameters.mle_n = para[(String) "mle_n"];
	if(para.Include_Key((String) "nn_n"))
		parameters.nn_n = para[(String) "nn_n"];

	if(method == (String) "CORR_DIM")
		num_r = parameters.cor_n;
	else
	if(method == (String) "NN_DIM")
		num_r = parameters.nn_n;
	else
	if(method == (String) "MLE_DIM")
		num_r = parameters.mle_n;
	else
	{
        cout << "Error: This estimator does not exist in this program!!" << endl;
		exit(0);
	}
};
#endif

void DimEst::Compute_Dim()
{
    std::cout << "Start program: " << D_prefname << "_" << method << std::endl;
    if(method == (String) "CORR_DIM")
        return Compute_Cor_Dim();

    if(method == (String) "NN_DIM")
        return Compute_NN_Dim();

	if(method == (String) "EIG_DIM")
		return Compute_Eig_Dim();

    if(method == (String) "MLE_DIM")
        return Compute_MLE_Dim();

    cout << "Error: This estimator does not exist in this program!!" << endl;
	exit(0);
};

void DimEst::Compute_Cor_Dim()
{
	Dim_C.resize(num_r, 1);
    Dim_D.resize(num_r - 2, 1);
	double D_min = D.matrix[1][0];
	double D_max = D.matrix[0][0];

	for(int i = 0; i < size; i++)
		for(int j = 0; j < i; j++)
		{
			if(D_min > D.matrix[i][j])
				D_min = D.matrix[i][j];
			if(D_max < D.matrix[i][j])
				D_max = D.matrix[i][j];
		}

	D_min = log(D_min + 0.001);
	D_max = log(D_max + 1);

	// initialize X
	for(int i = 0; i < num_r; i++)
		X.matrix[i][0] = (double) i / (num_r - 1) * (D_max - D_min) + D_min;

	double r = 0;
	for(int i = 0; i < num_r; i++)
	{
		r = exp(X.matrix[i][0]);
		Dim_C.matrix[i][0] = Cor_Dim_C(r);
		if(fabs(Dim_C.matrix[i][0]) > ZERO)
			Dim_C.matrix[i][0] = log(Dim_C.matrix[i][0]);
	}

	for(int i = 0; i < num_r - 2; i++)
	{
		Dim_D.matrix[i][0] = (Dim_C.matrix[i + 2][0] - Dim_C.matrix[i][0]) / (X.matrix[i + 2][0] - X.matrix[i][0]);
	}
};

void DimEst::Compute_NN_Dim()
{
	Dim_C.resize(num_r, 1);
	Dim_D.resize(num_r - 2, 1);
	Matrix<double> NND(num_r, 1);
	// initialize X
	for(int i = 0; i < num_r; i++)
		X.matrix[i][0] = i + 1;

	Compute_NN_Dist(NND);

	for(int i = 0; i < num_r; i++)
		X.matrix[i][0] = log(X.matrix[i][0]);

	for(int i = 0; i < num_r; i++)
		Dim_C.matrix[i][0] = log(NND.matrix[i][0]);

	for(int i = 0; i < num_r - 2; i++)
	{
		Dim_D.matrix[i][0] = (Dim_C.matrix[i + 2][0] - Dim_C.matrix[i][0]) / (log(X.matrix[i + 2][0]) - log(X.matrix[i][0]));
	}
};

void DimEst::Compute_NN_Dist(Matrix<double> &NND)
{
	int n = (int) (X.matrix[num_r - 1][0] + 1);
	double min_j = 0;
	int min_jk = 0;
	Matrix<double> MININ(n, 2);

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
			if(D.matrix[i][j] < ZERO)
			{
				MININ.matrix[0][0] = j;
				break;
			}

		for(int j = 1; j < n; j++)
		{
			min_j = 100000000;
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
			MININ.matrix[j][1] = min_j;
			MININ.matrix[j][0] = min_jk;
		}

		for(int j = 0; j < num_r; j++)
		{
			NND.matrix[j][0] += MININ.matrix[(int) X.matrix[j][0]][1];
		}
	}

	for(int i = 0; i < num_r; i++)
		NND.matrix[i][0] /= size;
};

void DimEst::Compute_Eig_Dim()
{
};

void DimEst::Compute_MLE_Dim()
{
	Dim_C.resize(num_r, 1);
	Dim_D.resize(1, 1);
	Matrix<double> Mk(num_r, 1);

	// initialize X
	for(int i = 0; i < num_r; i++)
		X.matrix[i][0] = i + 1;

	Compute_MLE_Mk(Mk);

	for(int i = 0; i < num_r; i++)
	{
		Dim_C.matrix[i][0] = Mk.matrix[i][0];
		Dim_D.matrix[0][0] += Mk.matrix[i][0];
	}
	Dim_D.matrix[0][0] /= num_r;
};

void DimEst::Compute_MLE_Mk(Matrix<double> &Mk)
{
	int n = (int) (X.matrix[num_r - 1][0] + 2);
	double min_j = 0;
	double mk = 0;
	int min_jk = 0;
	int end_k = 0;
	Matrix<double> MININ(n, 2);

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
			if(D.matrix[i][j] < ZERO)
			{
				MININ.matrix[0][0] = j;
				break;
			}

		for(int j = 1; j < n; j++)
		{
			min_j = 100000000;
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
			MININ.matrix[j][1] = min_j;
			MININ.matrix[j][0] = min_jk;
		}

		for(int j = 0; j < num_r; j++)
		{
			mk = 0;
			end_k = (int) (X.matrix[j][0]);
			for(int k = 1; k <= end_k; k++)
			{
				if(MININ.matrix[k][1] > ZERO)
					mk += log(MININ.matrix[end_k + 1][1] / MININ.matrix[k][1]);
			}
			if(mk > ZERO)
				mk = (end_k - 1) / mk;
			Mk.matrix[j][0] += mk;
		}
	}

	for(int i = 0; i < num_r; i++)
		Mk.matrix[i][0] /= size;
};

void DimEst::output_to_files()
{
	if(method == (String) "CORR_DIM" || method == (String) "NN_DIM")
	{
		String filename_Dim_C, filename_Dim_D;
		filename_Dim_C = makefilename("logvslog");
		filename_Dim_D = makefilename("deri_logvslog");
		File result1(filename_Dim_C);
		File result2(filename_Dim_D);
		result1.clean();
		result2.clean();

		if(method == (String) "CORR_DIM")
		{
			for(int i = 0; i < num_r; i++)
				result1 << X.matrix[i][0] << "\t" << Dim_C.matrix[i][0] << endl;

			for(int i = 0; i < num_r - 2; i++)
				result2 << X.matrix[i + 1][0] << "\t" << Dim_D.matrix[i][0] << endl;
		} else
		{
			for(int i = 0; i < num_r; i++)
				result1 << Dim_C.matrix[i][0] << "\t" << X.matrix[i][0] << endl;

			for(int i = 0; i < num_r - 2; i++)
				result2 << X.matrix[i + 1][0] << "\t" << (double) 1 / Dim_D.matrix[i][0] << endl;
		}
	} else
	if(method == (String) "MLE_DIM")
	{
		String filename_Dim_C, filename_Dim_D;
		filename_Dim_C = makefilename("k");
		filename_Dim_D = makefilename("dim");
		File result1(filename_Dim_C);
		File result2(filename_Dim_D);
		result1.clean();
		result2.clean();

		for(int i = 0; i < num_r; i++)
			result1 << X.matrix[i][0] << "\t" << Dim_C.matrix[i][0] << endl;

		result2 << Dim_D.matrix[0][0] << endl;
	}
	cout << "Please see the help for details of the introduction of output files." << endl;
};

String DimEst::makefilename(String output)
{
	String result;
	result = D_prefname;
	result += "_";
	result += method;
	result += "_";
	result += output;
	result += ".out";
	return result;
};

double DimEst::Cor_Dim_C(double r)
{
	int n = 0;
	for(int i = 0; i < size; i++)
		for(int j = i + 1; j < size; j++)
		{
			if(D.matrix[i][j] <= r)
				n++;
		}

	return (double) 2 * n / (size * (size - 1));
};

#endif
