
#include "Sparse_matrix.h"

SparseMatrix::SparseMatrix(int rows, int columns, double* Vals, int* RowInds, int*ColPtr)
{
	this->rows = rows;
    cols = columns;
	this->vals = Vals;
	colptr = ColPtr;
	rowind = RowInds;
}

SparseMatrix::~SparseMatrix()
{
    destructor();
}

void SparseMatrix::destructor()
{
    if (colptr != NULL)
    {
        if (colptr == rowind)
            rowind = NULL;
        delete[] colptr;
        colptr = NULL;
    }

    if (vals != NULL)
    {
        delete[] vals;
        vals = NULL;
    }

    if (rowind != NULL)
    {
        delete[] rowind;
        rowind = NULL;
    }
}

double SparseMatrix::operator()(int row, int column) const
{
	if (column <= cols / 2)
	{
		int i = colptr[column];
		while (i < colptr[column + 1] && rowind[i] < row)
		{
			i++;
		}
		if (i < colptr[column + 1] && rowind[i] == row)
		{
			return vals[i];
		}
	}
	else
	{
		int i = colptr[column + 1] - 1;
		while (i >= colptr[column] && rowind[i] > row)
		{
			i--;
		}
		if (i >= colptr[column] && rowind[i] == row)
		{
			return vals[i];
		}
	}
	return 0;
}

SparseMatrix* SparseMatrix::transpose() const
{
	int *temp = new int[rows];
    double* values = new double[colptr[cols]];
	int* newColptr = new int[rows + 1];
    int* newRowind = new int[colptr[cols]];
	for (int i = 0; i < rows; i++)
        temp[i] = 0;

    for (int i = 0; i < colptr[cols]; i++)
    {
        temp[rowind[i]]++;
    }

	int total = 0;
	for (int i = 0; i < rows; i++)
	{
		newColptr[i] = total;
		total += temp[i];
		temp[i] = newColptr[i];
	}
	newColptr[rows] = total;
    for (int i = 0; i < cols; i++)
	{
        for (int j = colptr[i]; j < colptr[i+1]; j++)
        {
            int k = temp[rowind[j]]++;
            newRowind[k] = i;
            values[k] = vals[j];
        }
    }
	delete[] temp;
    return new SparseMatrix(cols, rows, values, newRowind, newColptr);
}

SparseMatrix* SparseMatrix::Multiply(const SparseMatrix& mat) const
{
	double *resultcol = new double[rows];
	Vector<double> results(rows * mat.cols / 2);
	Vector<int> rowInds(rows * mat.cols / 2);
	SparseMatrix *result = new SparseMatrix(rows, mat.cols);
	result->colptr[0] = 0;
	for (int i = 0; i < rows; i++)
		resultcol[i] = 0;
	for (int i = 0; i < mat.cols; i++)
	{
		for (int k = mat.colptr[i]; k < mat.colptr[i+1]; k++)
		{
			for (int l = colptr[mat.rowind[k]]; l < colptr[mat.rowind[k]+1]; l++)
				resultcol[rowind[l]] += mat.vals[k] * vals[l];
		}
		for (int k = 0; k < rows; k++)
		{
			if (resultcol[k] != 0)
			{
				results.add(resultcol[k]);
				rowInds.add(k);
				resultcol[k] = 0;
			}
		}
		result->colptr[i+1] = results.size();
	}
	delete[] resultcol;
	result->setNNZ(results.size());
	for (int i = 0; i < results.size(); i++)
	{
		result->vals[i] = results[i];
		result->rowind[i] = rowInds[i];
	}
	return result;
}

double* SparseMatrix::Mean(const unsigned int& num_tree)
{
	double* result = new double[rows];
	for (int i = 0; i < rows; i++)
		result[i] = 0;
    for (int i = 0; i < cols; i++)
	{
       for (int j = colptr[i]; j < colptr[i+1]; j++)
            result[rowind[j]] += vals[j];

	}
	for (int i = 0; i < rows; i++)
		result[i] /= num_tree;
	return result;
}

double* SparseMatrix::Multiply_vec(const double* vec) const
{
	double* result = new double[rows];
	
	for (int i = 0; i < rows; i++)
		result[i] = 0;
	
    for (int i = 0; i < cols; i++)
	{
        if (vec[i] != 0 && colptr[i] < colptr[cols])
		{
			for (int j = colptr[i]; j < colptr[i+1]; j++)
				result[rowind[j]] += vec[i] * vals[j];
		}
	}
	return result;
}

// added by WH
void SparseMatrix::OutputSparseMatrix(std::ostream &output, SparseMatrixOutputType smtype)
{
    if(smtype == RCVLIST)
    {
        for (int i = 0; i < cols; i++)
        {
            for (int j = colptr[i]; j < colptr[i+1]; j++)
                output << rowind[j] << "\t" << i << "\t" << vals[j] << std::endl;
        }
    } else
    if(smtype == FULLMATRIX)
    {
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
                output << (*this)(i, j) << "\t";
            output << std::endl;
        }
    }
}

//// added by MM: Adds name of the columns when bipartition matrix is outputted in list format
//void SparseMatrix::OutputSparseMatrix2(std::ostream &output, SparseMatrixOutputType smtype)
//{
//    if(smtype == RCVLIST)
//    {
//        output << "Bipart \t Tree \t Weight" << std::endl;
//        for (int i = 0; i < cols; i++)
//        {
//            for (int j = colptr[i]; j < colptr[i+1]; j++)
//                output << rowind[j] + 1 << "\t" << i + 1 << "\t" << vals[j] << std::endl;
//        }
//    } else
//    if(smtype == FULLMATRIX)
//    {
//        for(int i = 0; i < rows; i++)
//        {
//            for(int j = 0; j < cols; j++)
//                output << (*this)(i, j) << "\t";
//            output << std::endl;
//        }
//    }
//}
