#pragma once

#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>

#include "array.hpp"
#include "zdtree.hpp"
#include "sparse.hpp"


template <class T>
class TreeSetObjects {
	// This class includes essential objects of trees that require further compuatations.
private:
	TreeSet<T>* Trees;
	Bipartition<T>* Bipart;
	TaxonList* Taxa;

	SparseMatrix* sb2t_mat;
	Matrix<precision>* cov_mat;
	Matrix<precision>* dis_mat;

	int n_trees;
	int n_leaves;
	int all_bipart;
	int	unique_bipart;

	bool ISWEIGHTED;
	bool ISROOTED;
	bool ISSORTED_Bipart;
public:
	TreeSetObjects() : Trees(nullptr), sb2t_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
		n_trees(0), all_bipart(0), unique_bipart(0), ISWEIGHTED(false), ISROOTED(false), ISSORTED_Bipart(false){};

	TreeSetObjects(TreeSet<T>* trees, bool rooted, bool weighted) : Trees(trees), sb2t_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
		ISWEIGHTED(weighted), ISROOTED(rooted) {
		set_basic_para();
	};

	void set_basic_para() {
		Bipart = Trees->get_bipart_ptr();
		Taxa = Trees->get_taxa_ptr();

		n_trees = Trees->get_size();
		n_leaves = Taxa->size; // This must be initialized after Taxa.

		all_bipart = Trees->get_all_bipart_size();
		unique_bipart = Trees->get_unique_bipart_size();

		ISSORTED_Bipart = false;
	}

	void set_weighted(bool weighted) { ISWEIGHTED = weighted; };

	//void sort_BipartId2Tree_by_treeid() {
	//	// Bubble sort implemented here.
	//	int temp;

	//	if (!ISSORTED_Bipart) {
	//		for (int i = 0; i < unique_bipart; i++) {
	//			size_t n_trees_Bipart_i = BipartId2Tree->get_c(i).get_size();
	//			if (BipartId2Tree->ele(i, 0) == -1)
	//				continue;
	//			else {
	//				while (n_trees_Bipart_i > 1) {
	//					size_t new_n = 0;
	//					for (int j = 1; j < n_trees_Bipart_i; j++) {
	//						if (BipartId2Tree->ele(i, j - 1) > BipartId2Tree->ele(i, j)) {
	//							// Swapped
	//							temp = BipartId2Tree->ele(i, j - 1);
	//							BipartId2Tree->ele(i, j - 1) = BipartId2Tree->ele(i, j);
	//							BipartId2Tree->ele(i, j) = temp;
	//							new_n = j;
	//						}
	//					}
	//					n_trees_Bipart_i = new_n;
	//				}
	//			}
	//		}
	//		ISSORTED_Bipart = true;
	//	}

	//}

	void Compute_Bipart() {
		if (Bipart->is_empty()) {
			for (int i = 0; i < Taxa->size; i++) {
				BitString<T> temp = BitString<T>(Taxa->bitstr_size, Taxa->size, i);
				Bipart->push(temp);
			}
		}
		
		for (int i = 0; i < n_trees; i++)
			Trees->ele(i).compute_bitstring(Bipart);

		// Update unique bipartition number
		unique_bipart = Bipart->get_size();
	};

	void Compute_Bipart_Matrix() {
		if (sb2t_mat != nullptr) {
			delete sb2t_mat;
			sb2t_mat = nullptr;
		}

		sb2t_mat = new SparseMatrix(unique_bipart, n_trees);

		sb2t_mat->initialize_CCS();

		for (int i = 0; i < n_trees; i++)
			sb2t_mat->push_CCS_col(Trees->ele(i).pop_t2b(), Trees->ele(i).pop_weight());

		sb2t_mat->set_RCS();
	}

	void Compute_Covariance_Matrix() {
		assert(sb2t_mat != nullptr);
		assert(unique_bipart > 0);
		if (cov_mat != nullptr)
			delete cov_mat;
		cov_mat = new Matrix<precision>(unique_bipart, unique_bipart, (precision)0);

		Array<precision> c_mean = sb2t_mat->col_mean();
		Array<precision> c_sum = sb2t_mat->col_sum();


		for (int i = 0; i < unique_bipart; i++)
		{
			(*cov_mat)(i, i) += (n_trees * c_mean[i] * c_mean[i] - 2 * c_sum[i] * c_mean[i]);
			for (int j = i + 1; j < unique_bipart; j++)
				(*cov_mat)(i, j) += (n_trees * c_mean[i] * c_mean[j] - c_sum[i] * c_mean[j] - c_sum[j] * c_mean[i]);
		}

		sb2t_mat->outer_prod_XXT_inplace((*cov_mat));

		for (int i = 0; i < unique_bipart; i++)
			for (int j = i + 1; j < unique_bipart; j++)
				(*cov_mat)(i, j) /= (n_trees - 1);
	};

	void Compute_Distance_Matrix() {};

	void Compute_RF_Distance_Matrix() {
		// Update formula: 
		// 1) D[j, k] += abs(M[i, j] - M[i, k]), for i in 1:row.
		// 2) D[j, k] = (D[j, k] + D[k, j]) / 4.
		// For 2), the computation of 1) guarantee that D is symmetric, which simply makes 2) become D = D / 2.
		assert(sb2t_mat != nullptr);
		assert(n_trees > 0);
		if (dis_mat != nullptr)
			delete dis_mat;
		dis_mat = new Matrix<precision>(n_trees, n_trees, (precision)0);
		sb2t_mat->set_RCS();

		for (int i = 0; i < unique_bipart; i++) {
			// The i-th row in B2T matrix.
			auto RCS_col_ind_ptr = sb2t_mat->get_RCS_ind_c_ptr(i);
			auto RCS_val_add_ptr = sb2t_mat->get_RCS_val_c_ptr(i);

			size_t RCS_row_size = (RCS_col_ind_ptr == nullptr) ? 0 : RCS_col_ind_ptr->get_size();

			// Handle 0 <= j < cur_col  with M[i, j] = 0 first.
			for (int j = 0; j < (*RCS_col_ind_ptr)[0]; j++) {
				size_t next_col;
				precision next_val;
				for (int RCS_next_ind = 0; RCS_next_ind < RCS_row_size; RCS_next_ind++) {
					// Found all nonzero behind j.
					next_col = (*RCS_col_ind_ptr)[RCS_next_ind];
					next_val = *((*RCS_val_add_ptr)[RCS_next_ind]) / 2;
					(*dis_mat)(j, next_col) += next_val;
					(*dis_mat)(next_col, j) += next_val;
				}
			}
			
			// Handle j = cur_col with M[i, j] != 0.
			for (int RCS_ind = 0; RCS_ind < RCS_row_size; RCS_ind++) {
				// RCS_ind pointing at nonzero entries on the i-th row.
				size_t cur_col = (*RCS_col_ind_ptr)[RCS_ind];
				precision cur_val = *((*RCS_val_add_ptr)[RCS_ind]) / 2;

				size_t RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				size_t next_col = (*RCS_col_ind_ptr)[RCS_next_ind];
				precision next_val = *((*RCS_val_add_ptr)[RCS_next_ind]) / 2;
				// Record the next nonzero. If this is the last nonzero in the row, record itself
				// which will not be triggered in later computation.

				for (int j = cur_col + 1; j < next_col; j++) {
					// Handle cur_col < j < next_col  with M[i, j] = 0 first.
					// For here, we only need M[i, k] != 0 with j < k.
					for (RCS_next_ind; RCS_next_ind < RCS_row_size; RCS_next_ind++) {
						// Found all nonzero behind j.
						next_col = (*RCS_col_ind_ptr)[RCS_next_ind];
						next_val = *((*RCS_val_add_ptr)[RCS_next_ind]) / 2;
						(*dis_mat)(next_col, j) += next_val;
						(*dis_mat)(j, next_col) += next_val;
					}
				}

				RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				next_col = (*RCS_col_ind_ptr)[RCS_next_ind];
				next_val = *((*RCS_val_add_ptr)[RCS_next_ind]) / 2;
				// Restore the record of the next nonzero.

				// Handle j = cur_col
				// For here, we need to consider all j < k.
				for (int k = cur_col + 1; k < n_trees; k++) {
					// Update D[k, cur_col] with k > cur_col, i.e. the lower triangle part.
					// Do the same with the upper triangle part, i.e., D[cur_col, k]
					// D[k, cur_col] += abs(M[i, k] - M[i, cur_col]);

					if (k != next_col) {
						// M[i, k] = 0;
						(*dis_mat)(k, cur_col) += cur_val;
						(*dis_mat)(cur_col, k) += cur_val;
					}
					else {
						// M[i, k] = next_val
						precision temp = (cur_val > next_val) ? (cur_val - next_val) : (next_val - cur_val);
						(*dis_mat)(k, cur_col) += temp;
						(*dis_mat)(cur_col, k) += temp;

						// update next nonzero
						RCS_next_ind = RCS_next_ind + (RCS_next_ind < RCS_row_size - 1);
						next_col = (*RCS_col_ind_ptr)[RCS_next_ind];
						next_val = *((*RCS_val_add_ptr)[RCS_next_ind]) / 2;
					}
				}
			}
		}
	}

	int get_all_bipart_num() { return all_bipart; };
	int get_unique_bipart_num() { return unique_bipart; };
	int get_tree_set_size() { return n_trees; };
	int get_leaf_set_size() { return n_leaves; };
	


	void print_Bipart_List(ostream& fout) {
		std::cout << "Printing Bipartition matrix!\n";
		std::cout << n_leaves << '\t' << unique_bipart << '\n';
		//fout << n_leaves << '\t' << unique_bipart << '\n';

		std::cout << 0 << " ";
		(*Bipart)[0].print_BitString(std::cout);
		std::cout << ' ' <<  (*sb2t_mat).get_RCS_ind_c_ptr(0)->get_size() << '\n';
		for (int i = 0; i < unique_bipart; i++){
			std::cout << i << " ";
			(*Bipart)[i].print_BitString(std::cout);
			std::cout << ' ' <<  (*sb2t_mat).get_RCS_ind_c_ptr(i)->get_size() << '\n';
		}
	};

	void print_Bipart2Tree_Matrix(ostream& fout, SparseMatrixOutputType smtype) { (*sb2t_mat).print(fout, smtype); };

	void print_Covariance_Matrix(ostream& fout) { (*cov_mat).print(fout, LOWERTRI); };

	void print_Distance_Matrix(ostream& fout) { (*dis_mat).print(fout, LOWERTRI); };

	void print_summary() {
		// int max_collision_size = 0, collision_cnt_5 = 0, collision_cnt_10 = 0, empty_cnt = 0;
		// int cur_size;
		// for (int i = 0; i < hash_bound; i++) {
		// 	if (!Hash2Id->get_c_ptr(i)->is_empty()) {
		// 		cur_size = Hash2Id->get_c_ptr(i)->get_size();
		// 		if (cur_size > max_collision_size)
		// 			max_collision_size = cur_size;
		// 		if (cur_size >= 5)
		// 			collision_cnt_5++;
		// 		if (cur_size >= 10)
		// 			collision_cnt_10++;
		// 		if (cur_size == 0)
		// 			empty_cnt++;
		// 		//std::cout << "Collusion beam # " << i << ":\t";
		// 		//for (int j = 0; j < cur_size; j++)
		// 		//	std::cout << Hash2Id[i]->operator[](j) <<"  ";
		// 		//std::cout << '\n';
		// 	}
		// 	else {
		// 		empty_cnt++;
		// 	}
		// }
		std::cout << "-------------Bipartition info summary------------------\n";
		std::cout << "\tTree set size: " << n_trees << ",\n";
		std::cout << "\tLeaf set size: " << n_leaves << ",\n";
		std::cout << "\tAll bipartitions encounted: " << all_bipart << ",\n";
		std::cout << "\tDistinct bipartition number: " << unique_bipart << ".\n";
		// std::cout << "\tHash container size: " << hash_bound << ",\n";
		// std::cout << "\tEmpty container: " << empty_cnt << ",\n";
		// std::cout << "\tMaximum Collision beam size: " << max_collision_size << ",\n";
		// std::cout << "\tCount of Collision beams over 5: " << collision_cnt_5 << ",\n";
		// std::cout << "\tCount of Collision beams over 10: " << collision_cnt_10 << ".\n";
		std::cout << "-------------Bipartition info summary------------------\n";
	}

};

