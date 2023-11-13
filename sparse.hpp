#pragma once
#include "array.hpp"
#include "SpecMat.hpp"
#include <cassert>

enum SparseMatrixOutputType
{
	RCVLIST,
	FULLMATRIX
};

//class dSparseMatrix {
//	// This is a extented column compressed storage(CCS) of sparse matrix with
//	// additional reduced row compressed storage: ColRef.
//	//
//	// The ColPtr in classic CCS is implemented in the data structure of RowInd,
//	// which is a linked-dynamic-block-list, darray_dc. Each block is a dynamic array
//	// storing row indices of nonzeros in a column. The pointer to the column
//	// is then the pointer to each block.
//public:
//	dSparseMatrix() : row(0), col(0), Val(nullptr), RowInd(nullptr), RowCCSInd(nullptr), ColInd(nullptr), flag_RCS(false) {};
//	dSparseMatrix(INTEGER_TYPE c) : dSparseMatrix() {
//		col = c;
//	}
//	dSparseMatrix(INTEGER_TYPE r, INTEGER_TYPE c, bool RCS = false) : dSparseMatrix() {
//		row = r;
//		col = c;
//		flag_RCS = RCS;
//	}
//	dSparseMatrix(const dSparseMatrix& src) : dSparseMatrix(src.row, src.col, src.flag_RCS) {
//		Val = new darray_dc<PRECISION>(*(src.Val));
//		RowInd = new darray_dc<INTEGER_TYPE>(*(src.RowInd));
//		//ColRef = (src.ColRef == nullptr) ? nullptr : new darray_dc<PRECISION*>(*src.ColRef);
//		ColInd = (src.ColInd == nullptr) ? nullptr : new darray_dc<INTEGER_TYPE>(*src.ColInd);
//	}
//
//	friend void swap(dSparseMatrix& lhs, dSparseMatrix& rhs) {
//		using std::swap;
//		swap(lhs.row, rhs.row);
//		swap(lhs.col, rhs.col);
//		swap(lhs.Val, rhs.Val);
//		swap(lhs.flag_RCS, rhs.flag_RCS);
//		swap(lhs.RowInd, rhs.RowInd);
//		//swap(lhs.ColRef, rhs.ColRef);
//		swap(lhs.ColInd, rhs.ColInd);
//	}
//
//	dSparseMatrix& operator=(const dSparseMatrix& rhs) {
//		dSparseMatrix temp = dSparseMatrix(rhs);
//		swap(*this, temp);
//		return *this;
//	}
//
//	class ColIterator {
//	public:
//		typedef typename dSparseMatrix::ColIterator self_type;
//		typedef typename darray_dc<PRECISION>::iterator val_it_type;
//		typedef typename darray_dc<INTEGER_TYPE>::iterator ind_it_type;
//		typedef PRECISION& reference;
//		typedef PRECISION* pointer;
//
//		dSparseMatrix::ColIterator() : val_it(val_it_type()), row_ind_it(ind_it_type()) {};
//		dSparseMatrix::ColIterator(val_it_type v_it, ind_it_type ri_it) : val_it(v_it), row_ind_it(ri_it) {};
//		dSparseMatrix::ColIterator(const self_type& src) : val_it(val_it_type(src.val_it)), row_ind_it(ind_it_type(src.row_ind_it)) {};
//
//		friend void swap(self_type& lhs, self_type& rhs) {
//			using std::swap;
//			swap(lhs.val_it, rhs.val_it);
//			swap(lhs.row_ind_it, rhs.row_ind_it);
//		}
//
//		self_type& operator=(const self_type& rhs) {
//			self_type temp(rhs);
//			swap(*this, temp);
//			return *this;
//		}
//
//		self_type operator++(int) { val_it++; row_ind_it++; return *this; };
//		self_type operator++() { self_type i = *this; val_it++; row_ind_it++; return i; };
//		// Currently, we do not support reverse order iterator for simplicity reason.
//
//		reference operator*() { return *val_it; };
//		pointer operator->() { return val_it.operator->(); };
//		INTEGER_TYPE get_row_ind() { return *row_ind_it; };
//		INTEGER_TYPE get_col_ind() { return val_it.get_block_ind(); };
//		INTEGER_TYPE get_row_CCS_ind() { return val_it.get_entry_ind(); };
//
//		bool operator==(const self_type& rhs) { return ((val_it == rhs.val_it) && (row_ind_it == rhs.row_ind_it)); }
//		bool operator!=(const self_type& rhs) { return ((val_it != rhs.val_it) || (row_ind_it != rhs.row_ind_it)); }
//
//	private:
//		val_it_type val_it;
//		ind_it_type row_ind_it;
//	};
//
//	class RowIterator {
//	public:
//		typedef typename dSparseMatrix::RowIterator self_type;
//		typedef typename darray_dc<PRECISION*>::iterator ref_it_type;
//		typedef typename darray_dc<INTEGER_TYPE>::iterator ind_it_type;
//		typedef darray_dc<PRECISION> val_arr_type;
//		typedef PRECISION& reference;
//		typedef PRECISION* pointer;
//
//		dSparseMatrix::RowIterator() : col_ind_it(ind_it_type()) , row_CCS_ind_it(ind_it_type()), Val(nullptr) {};
//		dSparseMatrix::RowIterator(ind_it_type ci_it, ind_it_type rci_it, val_arr_type* v) : col_ind_it(ci_it), row_CCS_ind_it(rci_it), Val(v) {};
//		dSparseMatrix::RowIterator(const self_type& src) : col_ind_it(ind_it_type(src.col_ind_it)), row_CCS_ind_it(src.row_CCS_ind_it), Val(src.Val) {};
//
//		friend void swap(self_type& lhs, self_type& rhs) {
//			using std::swap;
//			//swap(lhs.val_ref_it, rhs.val_ref_it);
//			swap(lhs.col_ind_it, rhs.col_ind_it);
//			swap(lhs.row_CCS_ind_it, rhs.row_CCS_ind_it);
//			swap(lhs.Val, rhs.Val);
//		}
//
//		self_type& operator=(const self_type& rhs) {
//			self_type temp(rhs);
//			swap(*this, temp);
//			return *this;
//		}
//
//		self_type operator++(int) { col_ind_it++; row_CCS_ind_it++; return *this; };
//		self_type operator++() { self_type i = *this; col_ind_it++; row_CCS_ind_it++; return i; };
//		// Currently, we do not support reverse order iterator for simplicity reason.
//
//
//
//		INTEGER_TYPE get_row_ind() { return col_ind_it.get_block_ind(); };
//		INTEGER_TYPE get_row_CCS_ind() { return *row_CCS_ind_it; };
//		INTEGER_TYPE get_col_ind() { return *col_ind_it; };
//
//		pointer operator->() { return Val->ele_ptr(get_col_ind(), get_row_CCS_ind()); };
//		reference operator*() { return *operator->(); };
//
//
//		bool operator==(const self_type& rhs) { return (col_ind_it == rhs.col_ind_it); }
//		bool operator!=(const self_type& rhs) { return (col_ind_it != rhs.col_ind_it); }
//
//	private:
//		//ref_it_type val_ref_it;
//		ind_it_type col_ind_it;
//		ind_it_type row_CCS_ind_it;
//		darray_dc<PRECISION>* Val;
//	};
//
//	typedef typename dSparseMatrix::RowIterator row_it_type;
//	typedef typename dSparseMatrix::ColIterator col_it_type;
//
//	col_it_type CCS_begin() { return col_it_type(Val->begin(), RowInd->begin()); };
//	col_it_type CCS_end() { return col_it_type(Val->end(), RowInd->end()); };
//	row_it_type RCS_begin() { return row_it_type(ColInd->begin(), RowCCSInd->begin(), Val); };
//	row_it_type RCS_end() { return row_it_type(ColInd->end(), RowCCSInd->begin(), Val); };
//
//	PRECISION* ele_ptr_CCS(INTEGER_TYPE i, INTEGER_TYPE j) {
//		// get the j-th nonzero in the i-th column
//		return Val->ele_ptr(i, j);
//	}
//
//	PRECISION& ele_CCS(INTEGER_TYPE i, INTEGER_TYPE j) {
//		// get the j-th nonzero in the i-th column
//		return *ele_ptr_CCS(i, j);
//	}
//
//	void initialize_CCS() {
//		assert(col > 0);
//		if (Val != nullptr)
//			delete Val;
//		Val = new darray_dc<PRECISION>(col);
//		if (RowInd != nullptr)
//			delete RowInd;
//		RowInd = new darray_dc<INTEGER_TYPE>(col);
//	}
//
//	void initialize_RCS() {
//		assert(row > 0);
//		//if (ColRef != nullptr)
//		//	ColRef->~darray_dc();
//		//ColRef = new darray_dc<PRECISION*>(row);
//		if (ColInd != nullptr)
//			delete ColInd;
//		ColInd = new darray_dc<INTEGER_TYPE>(row);
//		if (RowCCSInd != nullptr)
//			delete RowCCSInd;
//		RowCCSInd = new darray_dc<INTEGER_TYPE>(row);
//	}
//
//	void set_CCS_col(INTEGER_TYPE j, darray<PRECISION>* val_j, darray<INTEGER_TYPE>* rowind_j);
//
//	void push_CCS_col(darray<PRECISION>* val_j, darray<INTEGER_TYPE>* rowind_j);
//
//	void set_row(INTEGER_TYPE r) { row = r; };
//
//	void set_RCS();
//
//	void set_RCS_row(INTEGER_TYPE i, darray<INTEGER_TYPE>* val_ref_i, darray<INTEGER_TYPE>* colind_i);
//
//	void print(std::ostream& output, SparseMatrixOutputType smtype);
//
//private:
//	INTEGER_TYPE row;
//	INTEGER_TYPE col;
//	bool flag_RCS;
//	darray_dc<PRECISION>* Val;
//	darray_dc<INTEGER_TYPE>* RowInd;
//	//darray_dc<PRECISION&>* ColRef;
//	//darray_dc<PRECISION*>* ColRef;
//	darray_dc<INTEGER_TYPE>* RowCCSInd;
//	darray_dc<INTEGER_TYPE>* ColInd;
//};

class SparseMatrix
{
	typedef SparseMatrix self_type;
	typedef Array2D_tuple<INTEGER_TYPE, PRECISION> CCS_arr_type;
	typedef Array2D_tuple<INTEGER_TYPE, INTEGER_TYPE> RCS_arr_type;
	typedef int ind_type;
	typedef PRECISION val_type;
	typedef int add_type;

private:
	INTEGER_TYPE row;
	INTEGER_TYPE col;
	bool flag_RCS;
	CCS_arr_type *CCS;
	RCS_arr_type *RCS;
	// This is a extented column compressed storage(CCS) of sparse matrix with
	// additional reduced row compressed storage.
	//
	// The ColPtr in classic CCS is implemented in the data structure of RowInd,
	// which is a linked-dynamic-block-list, darray_dc. Each block is a dynamic array
	// storing row indices of nonzeros in a column. The pointer to the column
	// is then the pointer to each block.
public:
	SparseMatrix(INTEGER_TYPE r, INTEGER_TYPE c, CCS_arr_type *ptr_CCS, RCS_arr_type *ptr_RCS, bool RCS) : row(r),
																							   col(c),
																							   CCS(ptr_CCS),
																							   RCS(ptr_RCS),
																							   flag_RCS(RCS)
	{
		row = r;
		col = c;
		flag_RCS = RCS;
	}
	SparseMatrix() : SparseMatrix(0, 0, nullptr, nullptr, false){};

	SparseMatrix(INTEGER_TYPE c) : SparseMatrix(0, c, nullptr, nullptr, false) {}

	SparseMatrix(INTEGER_TYPE r, INTEGER_TYPE c) : SparseMatrix(r, c, nullptr, nullptr, false) {}

	SparseMatrix(const SparseMatrix &src) : SparseMatrix(src.row, src.col)
	{
		CCS = new CCS_arr_type(*(src.CCS));
		if (src.flag_RCS)
		{
			RCS = new RCS_arr_type(*(src.RCS));
			flag_RCS = true;
		}
	}

	friend void swap(SparseMatrix &lhs, SparseMatrix &rhs)
	{
		using std::swap;
		swap(lhs.row, rhs.row);
		swap(lhs.col, rhs.col);
		swap(lhs.flag_RCS, rhs.flag_RCS);
		swap(lhs.CCS, rhs.CCS);
		swap(lhs.RCS, rhs.RCS);
	}

	self_type &operator=(const self_type &rhs)
	{
		self_type temp = self_type(rhs);
		swap(*this, temp);
		return *this;
	}

	Array<val_type> &operator[](INTEGER_TYPE i)
	{
		return *(CCS->get_val_ptr(i));
	}

	Array<ind_type> &get_CCS_ind_c(INTEGER_TYPE i)
	{
		return *(CCS->get_key_ptr(i));
	}

	Array<val_type> &get_CCS_val_c(INTEGER_TYPE i)
	{
		return *(CCS->get_val_ptr(i));
	}

	Array<ind_type> *get_CCS_ind_c_ptr(INTEGER_TYPE i)
	{
		return CCS->get_key_ptr(i);
	}

	Array<val_type> *get_CCS_val_c_ptr(INTEGER_TYPE i)
	{
		return CCS->get_val_ptr(i);
	}

	Array<ind_type> &get_RCS_ind_c(INTEGER_TYPE i)
	{
		return *(RCS->get_key_ptr(i));
	}

	Array<add_type> &get_RCS_val_c(INTEGER_TYPE i)
	{
		return *(RCS->get_val_ptr(i));
	}

	Array<ind_type> *get_RCS_ind_c_ptr(INTEGER_TYPE i)
	{
		return RCS->get_key_ptr(i);
	}

	Array<add_type> *get_RCS_val_c_ptr(INTEGER_TYPE i)
	{
		return RCS->get_val_ptr(i);
	}

	INTEGER_TYPE get_row() { return row; };
	INTEGER_TYPE get_col() { return col; };

	INTEGER_TYPE get_nonzero_size()
	{
		return CCS->get_size();
	}

	val_type &operator()(INTEGER_TYPE j, INTEGER_TYPE e)
	{
		return CCS->get_val_ptr(j)->ele(e);
	}

	val_type &CCS_ele(INTEGER_TYPE j, INTEGER_TYPE e)
	{
		return CCS->get_val_ptr(j)->ele(e);
	}

	val_type &RCS_ele(INTEGER_TYPE i, INTEGER_TYPE e)
	{
		INTEGER_TYPE col_ind = RCS->get_key_ptr(i)->ele(e);
		INTEGER_TYPE row_CCS_ind = RCS->get_val_ptr(i)->ele(e);
		return CCS_ele(col_ind, row_CCS_ind);
	}

	val_type dense_ele(INTEGER_TYPE i, INTEGER_TYPE j)
	{
		assert(i < row && j < col);
		auto row_ind_ptr = CCS->get_key_ptr(j);
		INTEGER_TYPE col_size = row_ind_ptr->get_size();
		INTEGER_TYPE row_ind = 0;
		// Bisection search should be done here.
		for (INTEGER_TYPE it = 0; it < col_size; it++)
		{
			row_ind = (*row_ind_ptr)[it];
			if (row_ind == i)
				return this->RCS_ele(j, it);
			else if (row_ind > i)
				return 0;
		}
		return 0;
	}

	void transpose_shallow() { this->set_RCS(); };

	void transpose_inplace();

	SparseMatrix transpose();

	void initialize_CCS()
	{
		assert(col > 0);
		if (CCS != nullptr)
			delete CCS;
		CCS = new CCS_arr_type(col);
		//CCS->set_c_size(0);
	}

	void initialize_RCS()
	{
		assert(row > 0);
		if (RCS != nullptr)
			delete RCS;
		RCS = new RCS_arr_type(row);
		RCS->resize(row);
		for (int i = 0; i < row; i++)
			RCS->init_c(i, col);
	}

	void push_CCS_col(Array<INTEGER_TYPE> *rowind_ptr, Array<PRECISION> *val_ptr)
	{
		assert(CCS->get_size() < col);
		assert(rowind_ptr->get_size() == val_ptr->get_size());
		CCS->push_c(rowind_ptr, val_ptr);
	};

	void set_CCS_col(INTEGER_TYPE j, Array<INTEGER_TYPE> *rowind_ptr, Array<PRECISION> *val_ptr) { CCS->set_c(j, rowind_ptr, val_ptr); };

	void set_row(INTEGER_TYPE r) { row = r; };

	void set_RCS()
	{
		assert(row > 0);
		if (flag_RCS)
			return;
		if (RCS == nullptr)
			initialize_RCS();
		else
			RCS->earse();

		Array<PRECISION> *val_c_ptr;
		Array<INTEGER_TYPE> *row_ind_c_ptr;
		PRECISION *val_ptr;
		INTEGER_TYPE *row_ind_ptr;
		INTEGER_TYPE container_size;

		for (INTEGER_TYPE i = 0; i < col; i++)
		{
			row_ind_c_ptr = CCS->get_key_ptr(i);
			val_c_ptr = CCS->get_val_ptr(i);

			row_ind_ptr = row_ind_c_ptr->get_vec();
			val_ptr = val_c_ptr->get_vec();

			container_size = row_ind_c_ptr->get_size();
			for (INTEGER_TYPE j = 0; j < container_size; j++)
			{
				INTEGER_TYPE row_ind = row_ind_ptr[j];
				RCS->key()[row_ind].push(i);
				RCS->val()[row_ind].push(j);
			}
		}
		flag_RCS = true;
	};

	void set_RCS_row(INTEGER_TYPE i, Array<INTEGER_TYPE> *colind_ptr, Array<add_type> *valadd_ptr) { RCS->set_c(i, colind_ptr, valadd_ptr); };

	void print(std::ostream &output, SparseMatrixOutputType smtype)
	{
		if (smtype == RCVLIST)
		{
			Array<PRECISION> *val_c_ptr;
			Array<INTEGER_TYPE> *row_ind_c_ptr;
			PRECISION *val_ptr;
			INTEGER_TYPE *row_ind_ptr;
			INTEGER_TYPE container_size;

			for (INTEGER_TYPE i = 0; i < col; i++)
			{
				row_ind_c_ptr = CCS->get_key_ptr(i);
				val_c_ptr = CCS->get_val_ptr(i);

				row_ind_ptr = row_ind_c_ptr->get_vec();
				val_ptr = val_c_ptr->get_vec();

				container_size = row_ind_c_ptr->get_size();
				for (INTEGER_TYPE j = 0; j < container_size; j++)
					output << row_ind_ptr[j] << ' ' << i << ' ' << val_ptr[j] << '\n';
			}
		}
		else
		{
			this->set_RCS();
			for (int i = 0; i < row; i++)
			{
				int RCS_row_it = 0;
				Array<ind_type> *RCS_key_c_ptr = RCS->get_key_ptr(i);
				Array<add_type> *RCS_val_c_ptr = RCS->get_val_ptr(i);
				int RCS_row_size = (RCS_key_c_ptr == nullptr ? 0 : RCS_key_c_ptr->get_size());
				int cur_col = (RCS_row_size == 0) ? -1 : RCS_key_c_ptr->ele(RCS_row_it);
				add_type cur_CCS_row_ind = RCS_val_c_ptr->ele(RCS_row_it);

				for (int j = 0; j < col; j++)
				{
					if (j == cur_col)
					{
						output << CCS->ele_val(j, cur_CCS_row_ind) << '\t';
						if (RCS_row_it < RCS_row_size - 1)
						{
							RCS_row_it++;
							cur_col = RCS_key_c_ptr->ele(RCS_row_it);
							cur_CCS_row_ind = RCS_val_c_ptr->ele(RCS_row_it);
						}
					}
					else
						output << 0 << '\t';
				}
				output << '\n';
			}
		}
	}

	void sort_by_row_ind()
	{
		for (int i = 0; i < col; i++)
			CCS->sort_by_key(i);
	}

	Array<PRECISION> col_sum()
	{
		assert(flag_RCS);
		Array<PRECISION> Ans(row, row, 0.0);
		for (int i = 0; i < row; i++)
		{
			auto val_col_ind = RCS->get_key_ptr(i)->get_vec();
			auto val_CCS_ind = RCS->get_val_ptr(i)->get_vec();
			int row_size = RCS->get_val_ptr(i)->get_size();
			for (int j = 0; j < row_size; j++)
				Ans[i] += CCS->ele_val(val_col_ind[j], val_CCS_ind[j]);
		}
		return Ans;
	}

	Array<PRECISION> col_mean()
	{
		assert(flag_RCS);
		Array<PRECISION> Ans(row, row, 0.0);
		for (int i = 0; i < row; i++)
		{
			auto val_col_ind = RCS->get_key_ptr(i)->get_vec();
			auto val_CCS_ind = RCS->get_val_ptr(i)->get_vec();
			int row_size = RCS->get_val_ptr(i)->get_size();
			for (int j = 0; j < row_size; j++)
				Ans[i] += CCS->ele_val(val_col_ind[j], val_CCS_ind[j]);
		}
		return Ans;
	}

	Matrix<PRECISION> outer_prod_XXT()
	{
		Matrix<PRECISION> Ans(row, row, (PRECISION)0);

		PRECISION v_t = 0;

		for (int cur_col = 0; cur_col < col; cur_col++)
		{
			Array<ind_type> *CCS_key_c_ptr = CCS->get_key_ptr(cur_col);
			Array<val_type> *CCS_val_c_ptr = CCS->get_val_ptr(cur_col);
			int CCS_col_size = (CCS_key_c_ptr == nullptr ? 0 : CCS_key_c_ptr->get_size());
			for (int i = 0; i < CCS_col_size; i++)
			{
				int cur_row = (*CCS_key_c_ptr)[i];
				PRECISION cur_val = (*CCS_val_c_ptr)[i];
				Ans(cur_row, cur_row) += cur_val * cur_val;
				for (int j = i + 1; j < CCS_col_size; j++)
				{
					int next_row = (*CCS_key_c_ptr)[j];
					PRECISION next_val = (*CCS_val_c_ptr)[j];
					v_t = cur_val * next_val;
					Ans(cur_row, next_row) += v_t;
					Ans(next_row, cur_row) += v_t;
				}
			}
		}
		return Ans;
	}

	void outer_prod_XXT_inplace(Matrix<PRECISION> &M)
	{
		assert((M.get_row() == row) && (M.get_col() == row));
		PRECISION v_t = 0;

		for (int cur_col = 0; cur_col < col; cur_col++)
		{
			Array<ind_type> *CCS_key_c_ptr = CCS->get_key_ptr(cur_col);
			Array<val_type> *CCS_val_c_ptr = CCS->get_val_ptr(cur_col);
			int CCS_col_size = (CCS_key_c_ptr == nullptr ? 0 : CCS_key_c_ptr->get_size());
			for (int i = 0; i < CCS_col_size; i++)
			{
				int cur_row = (*CCS_key_c_ptr)[i];
				PRECISION cur_val = (*CCS_val_c_ptr)[i];
				M(cur_row, cur_row) += cur_val * cur_val;
				for (int j = i + 1; j < CCS_col_size; j++)
				{
					int next_row = (*CCS_key_c_ptr)[j];
					PRECISION next_val = (*CCS_val_c_ptr)[j];
					v_t = cur_val * next_val;
					M(cur_row, next_row) += v_t;
					M(next_row, cur_row) += v_t;
				}
			}
		}
	}

	void outer_prod_XXT_inplace(SpecMat::LowerTri<PRECISION> &M)
	{
		//assert((M.get_row() == row) && (M.get_col() == row));
		PRECISION v_t = 0;

		for (int cur_col = 0; cur_col < col; cur_col++)
		{
			Array<ind_type> *CCS_key_c_ptr = CCS->get_key_ptr(cur_col);
			Array<val_type> *CCS_val_c_ptr = CCS->get_val_ptr(cur_col);
			int CCS_col_size = (CCS_key_c_ptr == nullptr ? 0 : CCS_key_c_ptr->get_size());
			for (int i = 0; i < CCS_col_size; i++)
			{
				int cur_row = (*CCS_key_c_ptr)[i];
				PRECISION cur_val = (*CCS_val_c_ptr)[i];
				M(cur_row, cur_row) += cur_val * cur_val;
				for (int j = i + 1; j < CCS_col_size; j++)
				{
					int next_row = (*CCS_key_c_ptr)[j];
					PRECISION next_val = (*CCS_val_c_ptr)[j];
					v_t = cur_val * next_val;
					M(next_row, cur_row) += v_t;
				}
			}
		}
	}

	SparseMatrix *subMat_col(const Array<ind_type> &col_subset, Array<ind_type> &ind_mapping);

	SparseMatrix *subMat_row(const Array<ind_type> &row_subset, Array<ind_type> &ind_mapping);

	void compress_indices(Array<ind_type> **ind, INTEGER_TYPE n, Array<ind_type> &ind_mapping);
};
