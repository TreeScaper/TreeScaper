#include "sparse.hpp"



void SparseMatrix::transpose_inplace() {
	this->set_RCS();

	Array<ind_type>** col_ind_ptr = new Array<ind_type> * [this->row];
	Array<val_type>** val_ptr = new Array<val_type> * [this->row];
	for (int i = 0; i < this->row; i++) {
		col_ind_ptr[i] = nullptr;
		val_ptr[i] = nullptr;
	}

	for (int i = 0; i < this->row; i++) {
		auto cur_col_ind_ptr = this->RCS->pop_key_c(i);
		auto cur_CCS_ind_ptr = this->RCS->pop_val_c(i);
		if (cur_col_ind_ptr != nullptr) {
			size_t row_size = cur_col_ind_ptr->get_size();
			if (row_size) {
				auto val_in_row = new Array<val_type>(0, row_size);
				for (int j = 0; j < row_size; j++)
					(*val_in_row).push((*this)((*cur_col_ind_ptr)[j], (*cur_CCS_ind_ptr)[j]));
				col_ind_ptr[i] = cur_col_ind_ptr;
				val_ptr[i] = val_in_row;
			}
			else {
				col_ind_ptr[i] = new Array<ind_type>();
				val_ptr[i] = new Array<val_type>();
			}
		}
		else {
			col_ind_ptr[i] = new Array<ind_type>();
			val_ptr[i] = new Array<val_type>();
		}
 	}
	// Pop all RCS records, record number of nonempty rows in nonempty_row, col_ind in 

	SparseMatrix* trans = new SparseMatrix(this->col, this->row);
	(*trans).initialize_CCS();
	for (int i = 0; i < this->row; i++) 
		(*trans).push_CCS_col(col_ind_ptr[i], val_ptr[i]);
	(*trans).set_RCS();
	// Push the RCS records in the CCS record of temp

	swap(*this, *trans);
	
	delete trans;

	delete[] col_ind_ptr;
	delete[] val_ptr;

}

SparseMatrix SparseMatrix::transpose() {
	Array<ind_type>** col_ind_ptr = new Array<ind_type> * [this->row];
	Array<val_type>** val_ptr = new Array<val_type> * [this->row];
	for (int i = 0; i < this->row; i++) {
		col_ind_ptr = nullptr;
		val_ptr = nullptr;
	}

	for (int i = 0; i < this->row; i++) {
		auto temp_col_ind_ptr = this->RCS->get_key_ptr(i);
		auto temp_row_CCS_ptr = this->RCS->get_val_ptr(i);
		if (temp_col_ind_ptr != nullptr) {
			size_t row_size = temp_col_ind_ptr->get_size();
			if (row_size) {
				auto col_ind_in_row = new Array<ind_type>(0, row_size);
				auto val_in_row = new Array<val_type>(0, row_size);
				
				for (int j = 0; j < row_size; j++){
					(*col_ind_in_row).push((*temp_col_ind_ptr)[j]);
					(*val_in_row).push((*this)((*temp_col_ind_ptr)[j], (*temp_row_CCS_ptr)[j]));
				}
				col_ind_ptr[i] = col_ind_in_row;
				val_ptr[i] = val_in_row;
			}
			else {
				col_ind_ptr[i] = new Array<ind_type>();
				val_ptr[i] = new Array<val_type>();
			}
		}
		else {
			col_ind_ptr[i] = new Array<ind_type>();
			val_ptr[i] = new Array<val_type>();
		}
	}
	// Retreat all RCS records

	SparseMatrix trans = SparseMatrix(this->col, this->row);
	trans.initialize_CCS();
	for (int i = 0; i < this->row; i++)
		trans.push_CCS_col(col_ind_ptr[i], val_ptr[i]);
	trans.set_RCS();
	// Push the RCS records in the CCS record of temp

	delete[] col_ind_ptr;
	delete[] val_ptr;

	return trans;

}

SparseMatrix* SparseMatrix::subMat_col(const Array<ind_type>& col_subset, Array<ind_type>& ind_mapping) {
	size_t sub_col = col_subset.get_size();
	SparseMatrix* Ans = new SparseMatrix(sub_col);
	
	
	Ans->initialize_CCS();

	Array<ind_type>** sub_row_ind_ptr = new Array<ind_type> * [sub_col];
	for (int i = 0; i < sub_col; i++)
		sub_row_ind_ptr[i] = new Array<ind_type>(*(this->CCS->get_key_ptr(col_subset(i))));

	compress_indices(sub_row_ind_ptr, sub_col, ind_mapping);

	Ans->set_row(ind_mapping.get_size());


	for (int i = 0; i < sub_col; i++)
		Ans->push_CCS_col(sub_row_ind_ptr[i], new Array<val_type>(*(this->CCS->get_val_ptr(col_subset(i)))));
	
	// A lot of empty row expected, but we want t

	if (this->flag_RCS) {
		Ans->set_RCS();
	}
	
	delete[] sub_row_ind_ptr;

	return Ans;
}

SparseMatrix* SparseMatrix::subMat_row(const Array<ind_type>& row_subset, Array<ind_type>& ind_mapping) {
	size_t sub_row = row_subset.get_size();
	SparseMatrix* Ans = new SparseMatrix(sub_row);

	Ans->initialize_CCS();

	Array<ind_type>** sub_col_ind_ptr = new Array<ind_type> * [sub_row];
	Array<val_type>** sub_val_ptr = new Array<val_type> * [sub_row];
	for (int i = 0; i < sub_row; i++){
		sub_col_ind_ptr[i] = new Array<ind_type>(*(this->RCS->get_key_ptr(row_subset(i))));
		sub_val_ptr[i] = new Array<val_type>(0, sub_col_ind_ptr[i]->get_size());
		auto CCS_ind_ptr = this->RCS->get_val_ptr(row_subset(i));
		auto col_ind_ptr = sub_col_ind_ptr[i];
		for (int j = 0; j < col_ind_ptr->get_size(); j++)
			(*sub_val_ptr[i]).push(this->CCS_ele((*col_ind_ptr)[j], (*CCS_ind_ptr)[j]));
	}

	compress_indices(sub_col_ind_ptr, sub_row, ind_mapping);

	Ans->set_row(ind_mapping.get_size());

	for (int i = 0; i < sub_row; i++)
		Ans->push_CCS_col(sub_col_ind_ptr[i], sub_val_ptr[i]);

	Ans->transpose_inplace();
	
	delete[] sub_col_ind_ptr;
	delete[] sub_val_ptr;

	return Ans;
}

void SparseMatrix::compress_indices(Array<ind_type>** ind, size_t n, Array<ind_type>& ind_mapping) {
	size_t* row_ind = new size_t[n];
	size_t* row_val = new size_t[n];
	memset(row_ind, 0, n * sizeof(size_t));
	memset(row_val, 0, n * sizeof(size_t));

	Array<ind_type> active(0, n);
	Array<ind_type> indices(0, n);
	ind_type pos = 0;
	ind_type* CCS_ind = new ind_type[n];
	memset(CCS_ind, 0, n * sizeof(ind_type));
	for (int i = 0; i < n; i++) {
		pos = active.order_push((*ind[i])[CCS_ind[i]], false);
		indices.push(i, pos);
	}

	size_t r = 0;

	while (!active.is_empty()) {
		int temp = active.back();
		int col_ind = -1;
		if (temp == r) {
			while (active.back() == temp && !active.is_empty()) {
				active.pop();
				col_ind = indices.pop();
				if (CCS_ind[col_ind] < (ind[col_ind]->get_size() - 1)) {
					CCS_ind[col_ind]++;
					pos = active.order_push((*ind[col_ind])[CCS_ind[col_ind]], false);
					indices.push(col_ind, pos);
				}

			}
			ind_mapping.push(temp);
			r++;
		}
		else {
			while (active.back() == temp && !active.is_empty()) {
				active.pop();
				col_ind = indices.pop();
				auto cur_CCS_ind = CCS_ind[col_ind];
				(*ind[col_ind])[cur_CCS_ind] = r;
				if (CCS_ind[col_ind] < (ind[col_ind]->get_size() - 1)) {
					CCS_ind[col_ind]++;
					cur_CCS_ind = CCS_ind[col_ind];
					pos = active.order_push((*ind[col_ind])[cur_CCS_ind], false);
					indices.push(col_ind, pos);
				}
			}
			ind_mapping.push(temp);
			r++;
		}
	}

	delete[] CCS_ind;

	ind_mapping.align();

}