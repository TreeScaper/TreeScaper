#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>

enum SparseMatrixOutputType{RCVLIST, FULLMATRIX};

class SparseMatrixWH
{
public:
    SparseMatrixWH(){rows = 0; cols = 0; colptr = NULL; rowind = NULL; vals = NULL;}
    SparseMatrixWH(int rows, int columns, double *Vals, int *RowInds, int *ColPtr);
	~SparseMatrixWH();

    void destructor();
	
	// Returns the element at this row and column.
	double operator ()(int row, int column) const;
	
	// Return the transpose of the sparse matrix.
	SparseMatrixWH* transpose() const;
	
	// Return the product of two sparse matrix.
	SparseMatrixWH* Multiply(const SparseMatrixWH& mat) const;
	
	// Return a mean vector
	double* Mean(const int& num_tree);
	
	// Multiply a sparse matrix with a column vector
	double* Multiply_vec(const double* vec) const;

    // added by WH
    void OutputSparseMatrix(std::ostream &output, SparseMatrixOutputType smtype);

//    // add by MM
//    void OutputSparseMatrix2(std::ostream &output, SparseMatrixOutputType smtype);

private:

	int rows;
	int cols;

	int* colptr;	// length cols + 1
    int* rowind;	// length = number of entries = colptr[cols]
    double* vals;	// length = number of entries = colptr[cols]

	SparseMatrixWH(int rows, int columns, int nnz)
	{
		this->rows = rows;
		cols = columns;
		
		vals = new double[nnz];
		colptr = new int[cols + 1];
		rowind = new int[nnz];
	}
	
	SparseMatrixWH(int rows, int columns)
	{
		this->rows = rows;
		cols = columns;
		
		vals = NULL;
		rowind = NULL;
		colptr = new int[cols + 1];
	}
	
	inline void setNNZ(int nnz)
	{
		vals = new double[nnz];
		rowind = new int[nnz];
	}
	
	template <typename T>
	class Vector
	{
	public:
		inline Vector(int initialSize)
		{
			allocated = initialSize;
			count = 0;
			values = new T[allocated];
		}
		
		inline ~Vector()
		{
			delete[] values;
		}
		
		void add(const T& elem)
		{
			if (count == allocated)
			{
				T *temp = new T[allocated];
				for (int i = 0; i < allocated; i++)
				{
					temp[i] = values[i];
				}
				delete[] values;
				values = new T[allocated * 2];
				for (int i = 0; i < allocated; i++)
				{
					values[i] = temp[i];
				}
				delete[] temp;
				allocated *= 2;
			}
			values[count++] = elem;
		}
		
		inline const T& operator[] (int index) const
		{
			return values[index];
		}
		
		inline T& operator[] (int index)
		{
			return values[index];
		}
		
		inline int size() const
		{
			return count;
		}
		
	private:
		T *values;
		int count;
		int allocated;
	};
	
		
	
};

#endif // SPARSE_MATRIX_H
