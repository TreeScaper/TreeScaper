
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

// wmatrix.cpp
//         Member function definitions for class wvector.h
//                          by whuang Nov/15/2008

#ifndef WMATRIX_CPP
#define WMATRIX_CPP

#undef max
#undef min

#define mymin(a,b) ((a)<(b)?(a):(b))
#define mymax(a,b) ((a)>(b)?(a):(b))

#include <cstdlib>
#include <ctime>
#include <cmath>
#include "wmatrix.h"

template<class T>
const Matrix<T> &Matrix<T>::operator=(const Matrix<T> &right)
{
	//---matrix = right.matrix;
	(*this).destruct();
	row = right.row;
	col = right.col;
	matrix = new T *[row];
	for(int i = 0; i < row; i++)
		matrix[i] = new T [col];
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			matrix[i][j] = right(i, j);

	return *this;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T> &right) const
{
	if(row != right.row || col != right.col)
		return false;

	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			if(matrix[i][j] != right.matrix[i][j])
				return false;

	return true;
}

template<class T>
void Matrix<T>::resize(int r, int c)
{
	int min_row = (row > r) ? r : row;
	int min_col = (col > c) ? c : col;
	T *arr = new T[row * col];
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			arr[i * col + j] = this->matrix[i][j];
	this->destruct();
	matrix = new T *[r];
	for(int i = 0; i < r; i++)
		matrix[i] = new T [c];

	row = r;
	col = c;

	for(int i = 0; i < min_row; i++)
		for(int j = 0; j < min_col; j++)
			matrix[i][j] = arr[i * col + j];

	for(int i = min_row; i < r; i++)
		for(int j = 0; j < c; j++)
			matrix[i][j] = (T) NULL;
	
	for(int i = 0; i < min_row; i++)
		for(int j = min_col; j < c; j++)
			matrix[i][j] = (T) NULL;
	delete [] arr;
}

template<class T>
Matrix<T> Matrix<T>::get_vector(const int r, const int type) const
{
	Matrix<T> result;
	if(type == ROW)
	{
		result.resize(1, col);
		for(int i = 0; i < col; i++)
			result(0, i) = (*this)(r, i);
	} else
	if(type == COL)
	{
		result.resize(row, 1);
		for(int i = 0; i < row; i++)
			result(i, 0) = (*this)(i, r);
	}

	return result;
};

template<class T>
long double Matrix<T>::norm(const long double p) const           //  compute column vector p-norm.
{
	assert(col == 1);

	long double result = 0;

	if(p != INF_NORM)
	{
		for(int i = 0; i < row; i++)
			result += pow(fabs((long double) (matrix[i][0])), p);

		return  pow(result, (double) 1 / p);
	}

	for(int i = 0; i < row; i++)
        result = mymax(result, fabs((long double) matrix[i][0]));

	return result;
}

template<class T>
long double Matrix<T>::norm(const String &str) const                              // only compute p-norm p = 1, p = infinity and F
{
	Matrix<long double> mat;
	long double result = 0;

	if(str == (String) M_ONE_NORM)
	{
		mat.resize(col, 1);
		for(int i = 0; i < col; i++)
			for(int j = 0; j < row; j++)
				mat(i, 0) += fabs((long double) (matrix[j][i]));
		return mat.norm(INF_NORM);
	}
	else
	if(str == (String) M_INF_NORM)
	{
		mat.resize(row, 1);
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				mat(i, 0) += fabs((long double) (matrix[i][j]));
		return mat.norm(INF_NORM);
	}
	else
	if(str == (String) M_FRO_NORM)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				result += matrix[i][j] * matrix[i][j];

		return sqrt(result);
	}
	cout << "warning: wrong norm type\n";
	return 0;
}

template<class T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> result(col, row);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			result(j, i) = matrix[i][j];

	return result;
};

template<class T>
Matrix<T> Matrix<T>::generate(const int type, const int begin_row)
{
	if(type == ZERO_GEN)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = 0;

		return (*this);
	}else
	if(type == RAND_GEN)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = ((T) (rand() % 2001) - 1000) / 10;

		return (*this);
    } else
	if(type == E_GEN)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				if(i == j + begin_row)
					matrix[i][j] = 1;
				else
					matrix[i][j] = 0;

		return (*this);
	}else
	if(type == INPUT_GEN)
	{
		for(int i = 0; i < row; i++)
		{
			for(int j = 0; j < col; j++)
				cin >> matrix[i][j];
		}
	}else
	if(type == POS_DEF_DIA_GEN)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
			{
				if(i == j)
					matrix[i][j] = (((T) (rand() % 1000)) + 1) / 10;
				else
					matrix[i][j] = 0;
			}
	}else
	if(type == ORTH_GEN)
	{
		if(col < row)
			(*this).transpose();

		Matrix<T> r(1, col);

		(*this).generate(RAND_GEN);
		for(int i = 0; i < row; i++)
		{
			r = (*this).get_vector(i, ROW);

			for(int j = 0; j < i; j++)
				r = r - (((*this).get_vector(j, ROW))*(r.transpose()))(0,0) / 
				(((*this).get_vector(j, ROW))*(((*this).get_vector(j, ROW)).transpose()))(0, 0) *
				((*this).get_vector(j, ROW));

			r = r / (r.transpose().norm(2));

			for(int j = 0; j < col; j++)
				matrix[i][j] = r(0, j);
		}

		if(col < row)
			(*this).transpose();
	}else
	if(type == ORTH_Q_QAQ_GEN)
	{
		assert(row == col);
		Matrix<T> orth(row, col);
		orth.generate(ORTH_GEN);
		(*this) =  (orth) * (*this) * (orth.transpose());
	}else
	if(type == ONES_GEN)
	{
		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = 1;
	}

//POS_DEF_HH_GEN 7

	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::sub_matrix(const int row_begin, const int row_end, const int col_begin, const int col_end)
{
	Matrix<T> result(row_end - row_begin + 1, col_end - col_begin + 1);

	for(int i = 0; i < result.row; i++)
		for(int j = 0; j < result.col; j++)
			result(i, j) = matrix[row_begin + i][col_begin + j];

	return result;
}

template<class T>
Matrix<T> Matrix<T>::LU_factorization(Matrix<int> &P, Matrix<int> &Q, const int type)
{
	assert(col == row);
	P.resize(row, 1);
	Q.resize(row, 1);
	for(int i = 0; i < row; i++)
	{
		P(i, 0) = i;
		Q(i, 0) = i;
	}

	if(type == WITHOUT_PIVOTING)
		(*this).LU_without_pivoting();
	else
	if(type == PARTIAL_PIVOTING)
		(*this).LU_with_partial_pivoting(P);
	else
	if(type == COMPLETE_PIVOTING)
		(*this).LU_with_complete_pivoting(P, Q);
	else
		cout << "warning: there isn't that type.\n";
	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::LU_without_pivoting()
{
	for(int i = 0; i < row - 1; i++)
		(*this).i_step_of_gauss_elimination(i);
	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::LU_with_partial_pivoting(Matrix<int> &P)
{
	int max_row;

	for(int i = 0; i < row; i++)
	{
		max_row = i;
		for(int j = i + 1; j < row; j++)                  //select the abs maximum of the i column and i row to dim row
			if(fabs(matrix[j][i]) > fabs(matrix[max_row][i]))
				max_row = j;

		for(int j = 0; j < row; j++)                          //change i row with max_row row
			(*this).exchange(matrix[max_row][j], matrix[i][j]);

		(*this).exchange(P(i, 0), P(max_row, 0));                       //record the rows which changed

		(*this).i_step_of_gauss_elimination(i);
	}

	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::LU_with_complete_pivoting(Matrix<int> &P, Matrix<int> &Q)
{
	int max_row, max_column;

	for(int i = 0; i < row - 1; i++)
	{
		max_row = i;
		max_column = i;
		for(int j = i; j < row; j++)                  //select the abs maximum of the i--dim * i--dim
			for(int k = i; k < row; k++)
			{
				if(fabs(matrix[j][k]) > fabs(matrix[max_row][max_column]))
				{
					max_row = j;
					max_column = k;
				}
			}

		for(int j = 0; j < row; j++)                          //change i row with max_row
			(*this).exchange(matrix[max_row][j], matrix[i][j]);

		for(int j = 0; j < row; j++)
			(*this).exchange(matrix[j][max_column], matrix[j][i]);

		(*this).exchange(P(i, 0), P(max_row, 0));                       //record the rows which changed

		(*this).exchange(Q(i, 0), Q(max_column, 0));

		(*this).i_step_of_gauss_elimination(i);
	}

	return (*this);
}

template<class T>
void Matrix<T>::i_step_of_gauss_elimination(int i)
{
	int dim = row;
	if(fabs((*this)(i, i)) > ZERO)
		for(int j = i + 1; j < dim; j++)
			for(int k = i; k < dim; k++)
			{
				if(k > i)
					(*this)(j, k) = (*this)(j, k) - (*this)(i, k) * (*this)(j, i);
				else
					(*this)(j, k) = (*this)(j, k) / (*this)(i, i);   // k = i here
			}
	else
		cout << "warning:" << i << "-th step is divide by zero\n";
}

template<class T>
bool Matrix<T>::solve_tridiagonal_system(Matrix<T> diagonal, Matrix<T> upper_diagonal, Matrix<T> low_diagonal, Matrix<T> vec)
{
	int n = diagonal.get_row();
	double coef = 0;
	if(n - 1 != upper_diagonal.get_row() || n - 1 != low_diagonal.get_row() || n != vec.get_row())
	{
		cout << "warning: the dimensions doesn't match." << endl;
		return false;
	}
	this->resize(n, 1);

	for(int i = 0; i < n - 1; i++)
	{
		coef = - low_diagonal.matrix[i][0] / diagonal.matrix[i][0];
		low_diagonal.matrix[i][0] += diagonal.matrix[i][0] * coef;
		diagonal.matrix[i + 1][0] += upper_diagonal.matrix[i][0] * coef;
		vec.matrix[i + 1][0] += vec.matrix[i][0] * coef;
	}
	cout << "low:" << low_diagonal << endl;//--
	cout << "diagonal:" << diagonal << endl;//--
	cout << "upper:" << upper_diagonal << endl;//--
	if(fabs(diagonal.matrix[n - 1][0]) < ZERO)
	{
		cout << "warning: this function doesn't deal with the case when diagonal is too small" << endl;
		cout << n - 1 << " diagonal entry: " << diagonal.matrix[n - 1][0] << endl;
		return false;
	}
	this->matrix[n - 1][0] = vec.matrix[n - 1][0] / diagonal.matrix[n - 1][0];
	for(int i = n - 2; i >= 0; i--)
	{
		if(fabs(diagonal.matrix[i][0]) < ZERO)
		{
			cout << "warning: this function doesn't deal with the case when diagonal is too small" << endl;
			cout << i << " diagonal entry: " << diagonal.matrix[i][0] << endl;
			return false;
		}
		this->matrix[i][0] = (vec.matrix[i][0] - this->matrix[i + 1][0] * upper_diagonal.matrix[i][0]) / diagonal.matrix[i][0];
	}

	return true;
}

template<class T>
Matrix<T> Matrix<T>::solve_LU_linear_system(Matrix<T> &b, Matrix<int> &P, Matrix<int> &Q)
{
	b.permutation_row_vector(P);
	(*this).solve_low_triangular(b, LU);
	(*this).solve_upper_triangular(b);
	b.permutation_column_vector(Q);

	return b;
}

template<class T>
Matrix<T> Matrix<T>::permutation_row_vector(Matrix<int> &P)
{
	Matrix<T> result(row, 1);
	for(int i = 0; i < row; i++)
		result(i, 0) = matrix[P(i, 0)][0];

	return (*this) = result;
}

template<class T>
Matrix<T> Matrix<T>::solve_low_triangular(Matrix<T> &b, int type)
{
    (*this).resize(mymin(row, col), mymin(row, col));
	b.resize(col, 1);

    for(int i = 0; i < mymin(row, col); i++)
	{
		for(int j = 0; j < i; j++)
			b(i, 0) = b(i, 0) - b(j, 0) * (*this)(i, j);
		if(type != LU)
			b(i, 0) = b(i, 0) / (*this)(i, i);
	}
	return b;
}

template<class T>
Matrix<T> Matrix<T>::solve_upper_triangular(Matrix<T> &b)
{
    (*this).resize(mymin(row, col), mymin(row, col));
	b.resize(col, 1);

    for(int i = mymin(row, col) - 1; i >= 0; i--)
	{
        for(int j = mymin(row, col) - 1; j > i; j--)
			b(i, 0) = b(i, 0) - b(j, 0) * (*this)(i, j);
		b(i, 0) = b(i, 0) / (*this)(i, i);
	}

	return b;
}

template<class T>
Matrix<T> Matrix<T>::permutation_column_vector(Matrix<int> &Q)
{
	Matrix<T> result(row, 1);
	for(int i = 0; i < row; i++)
		result(Q(i, 0), 0) = matrix[i][0];

	return (*this) = result;
}

template<class T>
bool Matrix<T>::householder_reflector(Matrix<T> &b, Matrix<T> &norm, Matrix<T> &ext)
{  // norm store the norm of ||x||, ext store the first element of x of H. b is the constant of the problem,
	Matrix<T> x, v;
	norm.resize(col, 1);
	ext.resize(col, 1);

    for(int i = 0; i < mymin(row, col); i++)
	{
		x = ((*this).sub_matrix(i, row - 1, i, i));
		if(x.norm(2) <= ZERO)
		{
			cout << "warning: " << i << "-th norm almost equal zero";
			return false;
		}
		x(0, 0) = x(0, 0) + (*this)(i, i) * x.norm(2) / fabs((*this)(i, i));

		for(int j = 0; j < col - i; j++)
		{
			v = (*this).sub_matrix(i, row - 1, i + j, i + j) - 
				(2 * ((x.transpose()) * (*this).sub_matrix(i, row - 1, i + j, i + j))(0, 0)) / (pow(x.norm(2), 2)) * x;

			(*this).sub_matrix(i, row - 1, i + j, i + j) = v;
			for(int k = 0; k < row - i; k++)
			{
				(*this)(k + i, j + i) = v(k, 0);
			}
		}

		v = b.sub_matrix(i, row - 1, 0, 0) - 
			(2 * ((x.transpose()) * b.sub_matrix(i, row - 1, 0, 0))(0, 0)) / (pow(x.norm(2), 2)) * x;
		for(int k = 0; k < row - i; k++)
			b(k + i, 0) = v(k, 0);

		ext(i, 0) = x(0, 0);
		for(int k = i + 1; k < row; k++)
			(*this)(k, i) = x(k - i, 0);

		norm(i, 0) = pow(x.norm(2), 2);
	}
	return true;
}

template<class T>
bool Matrix<T>::get_bidiagonal(Matrix<T> &U, Matrix<T> &V, int end_col)
{
	Matrix<T> v, x;
	U.resize((*this).get_row(), (*this).get_row());
	V.resize((*this).get_col(), (*this).get_col());
	U.generate(E_GEN);
	V.generate(E_GEN);

	int row = (*this).get_row(), col = (*this).get_col();

	// get bidiagonal matrix
    for(int i = 0; i < mymin(row, col) - end_col; i++)
	{
		// operation for columns
		if(i < row - 1)
		{
			x = ((*this).sub_matrix(i, row - 1, i, i));
			if(x.norm(2) <= ZERO)
			{
				cout << "warning: " << i << "-th column norm almost equal zero";
				return false;
			}
			x(0, 0) = x(0, 0) + (*this)(i, i) * x.norm(2) / fabs((*this)(i, i));

			for(int j = 0; j < col - i; j++)
			{
				v = (*this).sub_matrix(i, row - 1, i + j, i + j) - 
					(2 * ((x.transpose()) * (*this).sub_matrix(i, row - 1, i + j, i + j))(0, 0)) / (pow(x.norm(2), 2)) * x;

				(*this).sub_matrix(i, row - 1, i + j, i + j) = v;
				for(int k = 0; k < row - i; k++)
				{
					(*this)(k + i, j + i) = v(k, 0);
				}
			}

			for(int j = 0; j < row; j++)
			{
				v = U.sub_matrix(i, row - 1, j, j) - 
					(2 * ((x.transpose()) * U.sub_matrix(i, row - 1, j, j))(0, 0)) / (pow(x.norm(2), 2)) * x;
				U.sub_matrix(i, row - 1, j, j) = v;
				for(int k = 0; k < row - i; k++)
				{
					U(k + i, j) = v(k, 0);
				}
			}
		}

		// then operation for rows
        if(i < col - 2 - mymax(0, end_col - 1))
		{
			x = ((*this).sub_matrix(i, i, i + 1, col - 1));
			if((x.transpose()).norm(2) <= ZERO)
			{
				cout << "warning: " << i << "-th row norm almost equal zero";
				return false;
			}

			x(0, 0) = x(0, 0) + (*this)(i, i) * ((x.transpose()).norm(2)) / fabs((*this)(i, i));

			for(int j = 0; j < row - i; j++)
			{
				v = (*this).sub_matrix(i + j, i + j, i + 1, col - 1) - 
					(2 * ((*this).sub_matrix(i + j, i + j, i + 1, col - 1) * (x.transpose()))(0, 0)) / (pow((x.transpose()).norm(2), 2)) * x;

				(*this).sub_matrix(i + j, i + j, i + 1, col - 1) = v;
				for(int k = 0; k < col - i - 1; k++)
				{
					(*this)(j + i, k + i + 1) = v(0, k);
				}
			}

			for(int j = 0; j < col; j++)
			{
				v = V.sub_matrix(j, j, i + 1, col - 1) - 
					(2 * (V.sub_matrix(j, j, i + 1, col - 1) * (x.transpose()))(0, 0)) / (pow((x.transpose()).norm(2), 2)) * x;
				V.sub_matrix(j, j, i + 1, col - 1) = v;
				for(int k = 0; k < col - i - 1; k++)
				{
					V(j, k + i + 1) = v(0, k);
				}
			}
		}
	}

	return true;
}

template<class T>
bool Matrix<T>::SVD(Matrix<T> &U, Matrix<T> &V)
{
	double u = 0;
	double p = 0, r = 0, o = 0, b1 = 0, b2 = 0;
	double update_data1 = 0, update_data2 = 0;
	double new_last_element = 0, old_last_element = 0;
	int row = (*this).get_row();
	int col = (*this).get_col();
	Matrix<T> Ui, Vi;

	(*this).get_bidiagonal(U, V);

	// for row >= column
	// start iteration
	for(int k = 1; k < col; k++)
	{
		do
		{
			old_last_element = (*this)(col - k - 1, col - k);
			// determine shift  (Wilkinson shift)
			u = pow((*this)(col - k, col - k), 2.0) + pow((*this)(col - k - 1, col - k), 2.0);

			// compute G1 = {{r, o}, {-o, r}}
			b1 = pow((*this)(0, 0), 2.0) - u;
			b2 = (*this)(0, 0) * (*this)(0, 1);
			p = sqrt(pow(b1, 2.0) + pow(b2, 2.0));
			r = b1 / p;
			o = -b2 / p;

			// compute B * G1
			update_data1 = (*this)(0, 0) * r - (*this)(0, 1) * o;
			update_data2 = (*this)(0, 0) * o + (*this)(0, 1) * r;
			(*this)(0, 0) = update_data1;
			(*this)(0, 1) = update_data2;

			update_data1 = (*this)(1, 0) * r - (*this)(1, 1) * o;
			update_data2 = (*this)(1, 0) * o + (*this)(1, 1) * r;
			(*this)(1, 0) = update_data1;
			(*this)(1, 1) = update_data2;

			// record G1 on V
			for(int i = 0; i < V.get_row(); i++)
			{
				update_data1 = V(i, 0) * r - V(i, 1) * o;
				update_data2 = V(i, 0) * o + V(i, 1) * r;
				V(i, 0) = update_data1;
				V(i, 1) = update_data2;
			}
			
			// get bidiagonal again....................can be optimized
			(*this).get_bidiagonal(Ui, Vi, k);

			// record change on V and U
			U = Ui * U;
			V = V * Vi;
			
			new_last_element = (*this)(col - k - 1, col - k);

		}while(fabs(new_last_element - old_last_element) > ZERO);
	}

	for(int i = 0; i < (*this).get_col(); i++)
	{
		if((*this)(i, i) <= 0)
		{
			for(int j = 0; j < (*this).get_row(); j++)
				(*this)(j, i) = - (*this)(j, i);
			for(int k = 0; k < (*this).get_col(); k++)
				V(k, i) = - V(k, i);
		}
	}

	return true;
}

template<class T>
bool Matrix<T>::SVD_LIB(Matrix<T> &U, Matrix<T> &S, Matrix<T> &Vt)
{
	integer M = row;
	integer N = col;
	integer LDA = M;
	integer LDU = M;
	integer LDVT = N;
	integer LWORK;
	integer INFO;
    integer mn = mymin( M, N );
    integer MN = mymax( M, N );
	int LWORKSIZE = 5 * MN;
	char JOBU = 'A';
	char JOBVT = 'A';
	double *wk = new double[LWORKSIZE];
	LWORK = LWORKSIZE;

	double *A = new double[row * col];
	double *ut = new double[row * row];
	double *v = new double[col * col];
    double *s = new double[mymin(row, col)];
	U.resize(row, row);
	S.resize(row, col);
	Vt.resize(col, col);
	
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			A[j * col + i] = matrix[i][j];

	// compute SVD
	cout << "computing SVD" << endl;
	dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, s, ut, &LDU, v, &LDVT, wk, &LWORK, &INFO);
	
	delete [] wk;
	
	for(int i = 0; i < row; i++)
		for(int j = 0; j < row; j++)
			U(i, j) = ut[j * row + i];

	for(int i = 0; i < col; i++)
		for(int j = 0; j < col; j++)
			Vt(i, j) = v[j * col + i];

        for(int i = 0; i < mn; i++)
		S(i, i) = s[i];

	delete [] A;
	delete [] ut;
	delete [] v;
	delete [] s;

	if(INFO == 0)
	{
		cout << "SVD successfully" << endl;
		return true;
	}
	else
	{
		cout << "warning: SVD fail!" << endl;
		return false;
	}
}

template<class T>
bool Matrix<T>::compute_inverse_matrix(Matrix<T> &inverse_matrix)
{
	assert(row == col);
	inverse_matrix.resize(row, col);
	real *V = new real[row * col];
	real *Vn = new real[row * col];
	integer *ipiv = new integer[row * col];
	integer M = row;
	integer lda = row;
	integer lwork = row;
	integer INFO;

	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			V[i * col + j] = matrix[i][j];

	sgetrf_(&M, &M, V, &M, ipiv, &INFO);
	if(INFO != 0)
	{
		cout << "warning: step1 fail" << endl;
		delete [] V;
		delete [] Vn;
		delete [] ipiv;
		return false;
	}

	sgetri_(&M, V, &lda, ipiv, Vn, &lwork, &INFO);

	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			inverse_matrix.matrix[i][j] = V[i * col + j];

	delete [] V;
	delete [] Vn;
	delete [] ipiv;
	if(INFO == 0)
	{
		return true;
	}
	else
	{
		cout << "warning: step2 fail" << endl;
		return false;
	}
}

template<class T>
int Matrix<T>::speedest_descent(Matrix<T> &b, Matrix<T> &x, double zero, int type, int type2, Matrix<T> &solution)
{//x is first point    b is problem vector    zero is the stop value
	int times = 0;

	if(type == NONE_PRECON)
		times = (*this).speedest_descent_none_precon(b, x, zero, type2, solution);
	else
	if(type == JACOBI_PRECON)
		times = (*this).speedest_descent_jacobi_precon(b, x, zero, type2, solution);
	else
	if(type == SYM_GS_PRECON)
		times = (*this).speedest_descent_sym_gs_precon(b, x, zero, type2, solution);

	return times;
}

template<class T>
int Matrix<T>::speedest_descent_none_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r;                       //r is residual
	Matrix<T> p;                       //p is direction
	double a;                          //a is distant
	int times = 0;

	r = b - (*this) * x;
	while(r.norm(2) > zero)
	{
		p = r;
		a = ((p.transpose() * r)(0, 0)) / ((p.transpose() * (*this) * p)(0, 0));

		x = x + a * p;
		r = b - (*this) * x;

		times++;
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}

	return times;
}

template<class T>
int Matrix<T>::speedest_descent_jacobi_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r;                       //r is residual
	Matrix<T> p;                       //p is direction
	Matrix<T> D(row, col);             //diagonal matrix
	double a;                          //a is distant
	int times = 0;

    for(int i = 0; i < mymin(row, col); i++)
		D(i, i) = matrix[i][i];

	r = b - (*this) * x;

	while(r.norm(2) > zero)
	{
		p = r;
		D.solve_upper_triangular(p);
		a = ((p.transpose() * r)(0, 0)) / ((p.transpose() * (*this) * p)(0, 0));
		x = x + a * p;
		r = b - (*this) * x;
		times++;
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}

	return times;
}

template<class T>
int Matrix<T>::speedest_descent_sym_gs_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r;                                               //r is residual
	Matrix<T> p;                                               //p is direction
	Matrix<T> D(row, col), D_L(row, col), D_U(row, col);       //diagonal matrix
	double a;                                                  //a is distant
	int times = 0;

	for(int i = 0; i < row; i++)                               //initial D, D - L and D - U
		for(int j = 0; j < col; j++)
		{
			if(i == j)
				D(i, i) = matrix[i][i];
			if(i >= j)
				D_L(i, j) = matrix[i][j];
			if(i <= j)
				D_U(i, j) = matrix[i][j];
		}

	r = b - (*this) * x;

	while(r.norm(2) > zero)
	{
		p = r;
		D_L.solve_low_triangular(p);                                  //this part is used
		D.solve_upper_triangular(p);                                  //to solve
		D_U.solve_upper_triangular(p);                                //(D - L)'D'(D_U)'r, here ' mean inverse
		a = ((p.transpose() * r)(0, 0)) / ((p.transpose() * (*this) * p)(0, 0));
		x = x + a * p;
		r = b - (*this) * x;
		times++;                                                      //record how many steps
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}

	return times;
}

template<class T>
Matrix<T> Matrix<T>::conjugate_gradient(Matrix<T> &b, Matrix<T> &x, double zero, int type, int type2, Matrix<T> &solution)
{
	if(type == NONE_PRECON)
		(*this).conjugate_gradient_none_precon(b, x, zero, type2, solution);
	else
	if(type == JACOBI_PRECON)
		(*this).conjugate_gradient_jacobi_precon(b, x, zero, type2, solution);
	else
	if(type == SYM_GS_PRECON)
		(*this).conjugate_gradient_sym_gs_precon(b, x, zero, type2, solution);

	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::conjugate_gradient_none_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r, r1;
	Matrix<T> p;
	Matrix<T> v;
	double alpha, beta;

	r = b - (*this) * x;
	p = r;

	for(int i = 0; i < row && x.norm(2) > zero; i++)
	{
		v = (*this) * p;
		alpha = (r.transpose() * r)(0, 0) / ((p.transpose() * v)(0, 0));
		x = x + alpha * p;
		r1 = r - alpha * v;
		beta = (r1.transpose() * r1)(0, 0) / (r.transpose() * r)(0, 0);
		p = r1 + beta * p;
		r = r1;
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}
	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::conjugate_gradient_jacobi_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r, r1;
	Matrix<T> z, z1;
	Matrix<T> D(row, col);
	Matrix<T> p;
	Matrix<T> v;
	double alpha, beta;

    for(int i = 0; i < mymin(row, col); i++)
		D(i, i) = matrix[i][i];

	r = b - (*this) * x;
	z = r;
	D.solve_upper_triangular(z);
	p = z;

	for(int i = 0; i < row && x.norm(2) > zero; i++)
	{
		v = (*this) * p;
		alpha = (r.transpose() * z)(0, 0) / ((p.transpose() * v)(0, 0));
		x = x + alpha * p;
		r1 = r - alpha * v;

		z1 = r1;
		D.solve_upper_triangular(z1);

		beta = (r1.transpose() * z1)(0, 0) / (r.transpose() * z)(0, 0);
		r = r1;
		z = z1;
		p = z + beta * p;
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}
	return (*this);
}

template<class T>
Matrix<T> Matrix<T>::conjugate_gradient_sym_gs_precon(Matrix<T> &b, Matrix<T> &x, double zero, int type, Matrix<T> &solution)
{
	Matrix<T> r, r1;
	Matrix<T> z, z1;
	Matrix<T> D(row, col), D_L(row, col), D_U(row, col);
	Matrix<T> p;
	Matrix<T> v;
	double alpha, beta;

	for(int i = 0; i < row; i++)                               //initial D, D - L and D - U
		for(int j = 0; j < col; j++)
		{
			if(i == j)
				D(i, i) = matrix[i][i];
			if(i >= j)
				D_L(i, j) = matrix[i][j];
			if(i <= j)
				D_U(i, j) = matrix[i][j];
		}

	r = b - (*this) * x;
	z = r;
	D_L.solve_low_triangular(z);                                  //this part is used
	D.solve_upper_triangular(z);                                  //to solve
	D_U.solve_upper_triangular(z);                                //(D - L)'D'(D_U)'r, here ' mean inverse
	p = z;

	for(int i = 0; i < row && x.norm(2) > zero; i++)
	{
		v = (*this) * p;
		alpha = (r.transpose() * z)(0, 0) / ((p.transpose() * v)(0, 0));
		x = x + alpha * p;
		r1 = r - alpha * v;
		z1 = r1;
		D_L.solve_low_triangular(z1);                                  //this part is used
		D.solve_upper_triangular(z1);                                  //to solve
		D_U.solve_upper_triangular(z1);                                //(D - L)'D'(D_U)'r, here ' mean inverse
		beta = (r1.transpose() * z1)(0, 0) / (r.transpose() * z)(0, 0);
		r = r1;
		z = z1;
		p = z + beta * p;
		if(type == TEST && solution.norm(2) != 0)
			cout << "||x - sol|| / ||sol||  : " << (x - solution).norm(2) / solution.norm(2) << "\n";
		else
		if(type == TEST)
			cout << "||x||  : " << (x).norm(2) << "\n";
	}

	return (*this);
}

#endif
