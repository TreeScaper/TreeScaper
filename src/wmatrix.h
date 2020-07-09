
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

//wmatrix.h
//        defination of some operation in mathecation
//                      by whuang

#ifndef WMATRIX_H
#define WMATRIX_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>
#include "wstring.h"

extern "C"
{
#include "f2c.h"
//#ifndef _WIN32
#include "clapack.h"
//#else
//#include <dgesvd.h>
//#include <sgetrf.h>
//#include <sgetri.h>
//#endif
}
#undef abs

#define ZERO pow(10.0, -8)

#define ROW 1
#define COL 2

#define INF_NORM 0

#define M_ONE_NORM "1"
#define M_INF_NORM "infinity"
#define M_FRO_NORM "frobenius"

#define ZERO_GEN 0
#define RAND_GEN 1
#define E_GEN 2
#define INPUT_GEN 3
#define POS_DEF_DIA_GEN 4
#define ORTH_GEN 5
#define ORTH_Q_QAQ_GEN 6
#define POS_DEF_HH_GEN 7
#define ONES_GEN 8

#define WITHOUT_PIVOTING 1
#define PARTIAL_PIVOTING 2
#define COMPLETE_PIVOTING 3

#define NONE_PRECON 1
#define JACOBI_PRECON 2
#define SYM_GS_PRECON 3

#define LU 1

#define SOLVE 1
#define TEST 2

using namespace std;

template<class T>
class Matrix{

	friend istream &operator>>(istream &input, Matrix<T> &mat)
	{
		return input;
	};

	friend ostream &operator<<(ostream &output, const Matrix<T> &mat)
	{
		cout << "{(" << mat.row << ", " << mat.col << ")\n";
		for(int i = 0; i < mat.row; i++)
		{
			for(int j = 0; j < mat.col; j++)
				if(j != mat.col - 1)
					cout << mat(i, j) << ", ";
				else
					cout << mat(i, j);
			if(i != mat.row - 1)
				cout << "\n";
		}
		cout << "\n}\n";

		return output;
	};

	friend Matrix operator+(Matrix<T> left, Matrix<T> right)
	{
		int max_row = max(left.get_row(), right.get_row());
		int max_col = max(left.get_col(), right.get_col());

		Matrix<T> result(max_row, max_col);
		left.resize(max_row, max_col);
		right.resize(max_row, max_col);

		for(int i = 0; i < max_row; i++)
			for(int j = 0; j < max_col; j++)
			{
				result(i, j) = left(i, j) + right(i, j);
			}

		return result;
	};
	
	template<class S>
	friend Matrix operator+(Matrix<T> left, S right)
	{
		Matrix<T> result(left.get_row(), left.get_col());

		for(int i = 0; i < left.get_row(); i++)
			for(int j = 0; j < left.get_col(); j++)
			{
				result.matrix[i][j] = left.matrix[i][j] + right;
			}

		return result;
	};

	template<class S>
	friend Matrix operator+(S left, Matrix<T> right)
	{
		Matrix<T> result(right.get_row(), right.get_col());

		for(int i = 0; i < right.get_row(); i++)
			for(int j = 0; j < right.get_col(); j++)
			{
				result(i, j) = left + right(i, j);
			}

		return result;
	};

	friend Matrix<T> operator-(Matrix<T> left, Matrix<T> right)
	{
		int max_row = max(left.get_row(), right.get_row());
		int max_col = max(left.get_col(), right.get_col());

		Matrix<T> result(max_row, max_col);
		left.resize(max_row, max_col);
		right.resize(max_row, max_col);

		for(int i = 0; i < max_row; i++)
			for(int j = 0; j < max_col; j++)
			{
				result(i, j) = left(i, j) - right(i, j);
			}

		return result;
	};

	template<class S>
	friend Matrix operator-(Matrix<T> left, S right)
	{
		Matrix<T> result(left.get_row(), left.get_col());

		for(int i = 0; i < left.get_row(); i++)
			for(int j = 0; j < left.get_col(); j++)
			{
				result.matrix[i][j] = left.matrix[i][j] - right;
			}

		return result;
	};

	template<class S>
	friend Matrix operator-(S left, Matrix<T> right)
	{
		Matrix<T> result(right.get_row(), right.get_col());

		for(int i = 0; i < right.get_row(); i++)
			for(int j = 0; j < right.get_col(); j++)
			{
				result(i, j) = left - right(i, j);
			}

		return result;
	};

	friend Matrix operator*(const Matrix<T> &left, const Matrix<T> &right)
	{
		assert(left.col == right.row);
		Matrix<T> result(left.row, right.col);
		for(int i = 0; i < right.col; i++)
			for(int j = 0; j < left.row; j++)
				for(int k = 0; k < left.col; k++)
					result.matrix[j][i] += left.matrix[j][k] * right.matrix[k][i];

		return result;
	};

	template<class S>
	friend Matrix operator*(S value, Matrix<T> mat)
	{
		Matrix<T> result(mat.row, mat.col);
		for(int i = 0; i < mat.row; i++)
			for(int j = 0; j < mat.col; j++)
				result(i, j) = value * mat(i, j);

		return result;
	};

	template<class S>
	friend Matrix operator*(Matrix<T> mat, S value)
	{
		Matrix<T> result = value * mat;
		return result;
	};

	template<class S>
	friend Matrix operator/(Matrix<T> mat, S value)
	{
		Matrix<T> result(mat.row, mat.col);
		for(int i = 0; i < mat.row; i++)
			for(int j = 0; j < mat.col; j++)
				result(i, j) = mat(i, j) / value;

		return result;
	};

public:

	template <class S, int N>
	Matrix(const S (&arr)[N])                                                   // intialize and copy a nomal array to a Array
	{
		int r = SIZE <N> ::cnt;
		int c = length_array(*arr);
		row = 0;
		col = 0;
		(*this).resize(r, c);

		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = (T) arr[i][j];
	};

	Matrix(Matrix &M)
	{
		row = M.row;
		col = M.col;
		matrix = new T *[row];
		for(int i = 0; i < row; i++)
			matrix[i] = new T [col];

		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = M.matrix[i][j];
	};

	Matrix(const Matrix &M)
	{
		row = M.row;
		col = M.col;
		matrix = new T *[row];
		for(int i = 0; i < row; i++)
			matrix[i] = new T [col];

		for(int i = 0; i < row; i++)
			for(int j = 0; j < col; j++)
				matrix[i][j] = M.matrix[i][j];
	}

	Matrix(const int r, const int c)
	{
		row = 0;
		col = 0;

		(*this).resize(r, c);
	};

	Matrix()
	{
		row = 0;
		col = 0;
	};

	~Matrix()
	{
		this->destruct();
	};

	void destruct()
	{
		if(row == 0 || col == 0)
			return;
		for(int i = 0; i < row; i++)
			delete [] matrix[i];
		delete [] matrix;
	};

	Matrix<T> conjugate_gradient(Matrix<T> &, Matrix<T> &, double, int = NONE_PRECON, int = SOLVE, Matrix<T> & = Matrix<T> ());

	Matrix<T> conjugate_gradient_none_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	Matrix<T> conjugate_gradient_jacobi_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	Matrix<T> conjugate_gradient_sym_gs_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	int speedest_descent(Matrix<T> &, Matrix<T> &, double, int = NONE_PRECON, int = SOLVE, Matrix<T> & = Matrix<T> ());

	int speedest_descent_none_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	int speedest_descent_jacobi_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	int speedest_descent_sym_gs_precon(Matrix<T> &, Matrix<T> &, double, int = SOLVE, Matrix<T> & = Matrix<T> ());

	bool householder_reflector(Matrix<T> &, Matrix<T> &, Matrix<T> &);

	bool get_bidiagonal(Matrix<T> &, Matrix<T> &, int = 0);

	bool SVD(Matrix<T> &, Matrix<T> &);

	bool SVD_LIB(Matrix<T> &, Matrix<T> &, Matrix<T> &);

	bool compute_inverse_matrix(Matrix<T> &);

	bool solve_tridiagonal_system(Matrix<T> diagonal, Matrix<T> upper_diagonal, Matrix<T> low_diagonal, Matrix<T> vec);

	Matrix<T> solve_LU_linear_system(Matrix<T> &, Matrix<int> &, Matrix<int> &);

	Matrix<T> permutation_row_vector(Matrix<int> &);

	Matrix<T> solve_low_triangular(Matrix<T> &, int = 0);

	Matrix<T> solve_upper_triangular(Matrix<T> &);

	Matrix<T> permutation_column_vector(Matrix<int> &);

	Matrix<T> LU_factorization(Matrix<int> &, Matrix<int> &, const int = WITHOUT_PIVOTING);

	Matrix<T> LU_without_pivoting();

	Matrix<T> LU_with_partial_pivoting(Matrix<int> &);

	Matrix<T> LU_with_complete_pivoting(Matrix<int> &, Matrix<int> &);

	void i_step_of_gauss_elimination(int);

	Matrix<T> sub_matrix(const int, const int, const int, const int);

	Matrix<T> generate(const int, const int = 0);

	Matrix<T> transpose();

	Matrix<double> compute_scalar_product_matrix()
	{
		assert(row == col);
		int size = row;
		Matrix<double> Scalar_Product(size, size);

		double *row_sum = new double[size];
		double *col_sum = new double[size];
		double total_sum = 0;

		for(int i = 0; i < size; i++)
		{
			row_sum[i] = 0;
			col_sum[i] = 0;
		}

		// compute S
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				Scalar_Product.matrix[i][j] = matrix[i][j] * matrix[i][j];//--D[i * size + j] * D[i * size + j];
				Scalar_Product.matrix[j][i] = Scalar_Product.matrix[i][j];
				row_sum[i] += Scalar_Product.matrix[i][j];
				col_sum[j] += Scalar_Product.matrix[i][j];
				row_sum[j] += Scalar_Product.matrix[j][i];
				col_sum[i] += Scalar_Product.matrix[j][i];
				total_sum += 2 * Scalar_Product.matrix[i][j];
			}
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
				Scalar_Product.matrix[i][j] = -0.5 * (Scalar_Product.matrix[i][j] - row_sum[i] / size - col_sum[j] / size + total_sum / size / size);

		delete [] row_sum;
		delete [] col_sum;
		return Scalar_Product;
	};

	Matrix<double> compute_Distance_Matrix()
	{
		Matrix<double> COR = (*this);
		int size = COR.get_row();
		int dim = COR.get_col();
		Matrix<double> D(size, size);
		double ddd = 0;
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < i; j++)
			{
				D.matrix[i][j] = 0;
				D.matrix[i][j] = 0;
				for(int k = 0; k < dim; k++)
				{
					ddd = COR.matrix[i][k] - COR.matrix[j][k];
					D(i, j) += ddd * ddd;
				}
				D.matrix[i][j] = sqrt((double) D.matrix[i][j]);
				D.matrix[j][i] = D.matrix[i][j];
			}
		}
		return D;
	};

	long double norm(const long double) const;                              // only compute vector, that is the number of matrix' column is 1

	long double norm(const String &) const;                              // only compute p-norm p = 1, p = infinity and F for a matrix

	Matrix<T> get_vector(const int , const int) const;

	const Matrix &operator=(const Matrix<T> &);

	bool operator==(const Matrix<T> &) const;

	bool operator!=(const Matrix<T> &right) const{return ! ((*this) == right);};

	T &operator()(const int r, const int c = 0){return matrix[r][c];};

	const T &operator()(const int r, const int c = 0) const{return matrix[r][c];};

	bool is_empty(){return (row == 0 && col == 0);};

	int get_row(){return row;};

	int get_row() const{return row;};

	int get_col(){return col;};

	int get_col() const{return col;};

	void resize(int, int);

	T **matrix;
private:

	template<int N>
	struct SIZE{                                                               // as a tool to store the length of a array
		static const int cnt = N;
	};

	template<class S, int N>
	int length_array(const S (&arr)[N])
	{
		return SIZE<N> :: cnt;
	};

	template<class S>
	S exchange(S &left, S &right)
	{
		S medium;
		medium = left;
		left = right;
		right = medium;
		return left;
	};

	int row, col;
};

#endif
