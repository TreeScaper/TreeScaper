#pragma once

#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>

#include "array.hpp"
#include "sparse.hpp"
#include "SpecMat.hpp"
#include "zdtree.hpp"
#include "zdfileio.hpp"

template <class T>
class TreeSetObjects
{
	// This class includes essential objects of trees that require further compuatations.
private:
	TreeSet<T> *Trees;
	Bipartition<T> *Bipart;
	TaxonList *Taxa;

	SparseMatrix *sb2t_mat;
	// Matrix<PRECISION>* cov_mat;
	// Matrix<PRECISION>* dis_mat;

	SpecMat::LowerTri<PRECISION> *cov_mat;
	SpecMat::LowerTri<PRECISION> *dis_mat;

	int n_trees;
	int n_leaves;
	int all_bipart;
	int unique_bipart;

	int n_sub_bipart;
	int n_sub_trees;

	bool ISWEIGHTED;
	bool ISROOTED;
	bool ISSORTED_Bipart;

public:
	TreeSetObjects() : Trees(nullptr), sb2t_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
					   n_trees(0), all_bipart(0), unique_bipart(0), n_sub_bipart(0), n_sub_trees(0), ISWEIGHTED(false), ISROOTED(false), ISSORTED_Bipart(false){};

	TreeSetObjects(TreeSet<T> *trees, bool rooted, bool weighted) : Trees(trees), sb2t_mat(nullptr), cov_mat(nullptr), dis_mat(nullptr),
																	ISWEIGHTED(weighted), ISROOTED(rooted)
	{
		set_basic_para();
	};

	~TreeSetObjects()
	{
		if (Trees != nullptr)
			delete Trees;
		if (Bipart != nullptr)
			delete Bipart;
		if (sb2t_mat != nullptr)
			delete sb2t_mat;
		if (cov_mat != nullptr)
			delete cov_mat;
		if (dis_mat != nullptr)
			delete dis_mat;
	}

	void set_basic_para()
	{
		Bipart = Trees->get_bipart_ptr();
		Taxa = Trees->get_taxa_ptr();

		n_trees = Trees->get_size();
		n_leaves = Taxa->size; // This must be initialized after Taxa.

		all_bipart = Trees->get_all_bipart_size();
		unique_bipart = Trees->get_unique_bipart_size();

		n_sub_bipart = unique_bipart;
		n_sub_trees = n_trees;

		ISSORTED_Bipart = false;
	}

	void set_weighted(bool weighted) { ISWEIGHTED = weighted; };

	// void sort_BipartId2Tree_by_treeid() {
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

	int Compute_Strict_Consensus_Tree(int *bipart_id, PRECISION *weights)
	{
		sb2t_mat->set_RCS();
		int bipart_num = sb2t_mat->get_row();
		int tree_num = sb2t_mat->get_col();

		// Array<int> *consensus_bipart_id = new Array<int>(0, row);
		// Array<int> *consensus_bipart_id = new Array<int>(0, bipart_num);
		int num_con_edge = 0;

		for (int i = 0; i < bipart_num; i++)
			if (sb2t_mat->get_RCS_ind_c_ptr(i)->get_size() == tree_num)
			{
				bipart_id[num_con_edge] = i;
				weights[num_con_edge] = 1.0;
				num_con_edge += 1;
			}
		return num_con_edge;
	}

	int Compute_Major_Consensus_Tree(int *bipart_id, PRECISION *weights)
	{
		sb2t_mat->set_RCS();
		int bipart_num = sb2t_mat->get_row();
		int tree_num = sb2t_mat->get_col();

		// Array<int> *consensus_bipart_id = new Array<int>(0, row);
		// Array<int> *consensus_bipart_id = new Array<int>(0, bipart_num);
		int num_con_edge = 0;
		PRECISION w = 0.0;

		for (int i = 0; i < bipart_num; i++)
		{
			w = sb2t_mat->get_RCS_ind_c_ptr(i)->get_size() / tree_num;
			if (w > 0.5)
			{
				bipart_id[num_con_edge] = i;
				weights[num_con_edge] = w;
				num_con_edge += 1;
			}
		}
		return num_con_edge;
	}

	// void print_Consensus_Tree(std::ostream &out, int *bipart_id, PRECISION *weights, int num_edge)
	// {
	// 	assert(Trees->issameleaf());
	// 	auto tree_ll = TreeLinkedList<T>(*Trees->get_uniform_leafset(), bipart_id, weights, num_edge, Bipart);
	// 	tree_ll.print_NewickStr(out, *Taxa);
	// }

	void Compute_Bipart()
	{
		std::cout << "-----Enter Compute_Bipart-----\n\n";


		if (Bipart->is_empty())
		{
			std::cout << "Initializing bipartition set...";
			for (int i = 0; i < Taxa->size; i++)
			{
				BitString<T> temp = BitString<T>(Taxa->bitstr_size, Taxa->size, i);
				Bipart->push(temp);
			}
		}

		if(Trees->issameleaf())
		{
			std::cout << "Computing bitstring of the uniform leaf set...";
			BitString<T>* uniform_leafset = new BitString<T>(Taxa->bitstr_size, Taxa->size);
			auto tree = Trees->ele(0);
			auto bitstrs = *Trees->ele(0).get_t2b();
			// for (int i = 0; i < tree.get_size(); i++)
			// 	(*uniform_leafset) |= (*Bipart)(bitstrs(i));
			
			Trees->set_uniform_leafset(uniform_leafset);

			for (int i = 0; i < n_trees; i++)
				Trees->ele(i).compute_bitstring(Bipart);
		}
		else
		{
			Trees->init_leafsets();
			for (int i = 0; i < n_trees; i++)
				Trees->set_leafsets(i, Trees->ele(i).compute_bitstring(Bipart));
		}

		// Update unique bipartition number
		unique_bipart = Bipart->get_size();

		std::cout << "-----Leave Compute_Bipart-----\n\n";

	};

	void Compute_Bipart_Matrix()
	{
		if (sb2t_mat != nullptr)
		{
			delete sb2t_mat;
			sb2t_mat = nullptr;
		}

		if (Trees->issameleaf())
			sb2t_mat = new SparseMatrix(unique_bipart, n_trees);
		else
			sb2t_mat = new SparseMatrix(2 * unique_bipart, n_trees);

		sb2t_mat->initialize_CCS();

		for (int i = 0; i < n_trees; i++)
			if (Trees->issameleaf())
				sb2t_mat->push_CCS_col(Trees->ele(i).pop_t2b(), Trees->ele(i).pop_weight());
			else
			{
				auto bipart_lower = Trees->ele(i).pop_t2b();
				auto bipart_upper = Trees->ele(i).pop_t2b_upper();

				auto weight_lower = Trees->ele(i).pop_weight();

				// pointers above will be released while existing.

				auto full_treevec = new Array<size_t>(join(*bipart_lower, *bipart_upper));
				auto full_weights = new Array<PRECISION>(duplicate(2, *weight_lower));

				// pointers above will be sent to the sparse matrix and held there.

				sb2t_mat->push_CCS_col(full_treevec, full_weights);

				delete bipart_lower;
				delete bipart_upper;
				delete weight_lower;
			}

		sb2t_mat->set_RCS();
	}

	void Compute_Covariance_Matrix()
	{
		assert(sb2t_mat != nullptr);
		assert(unique_bipart > 0);
		if (cov_mat != nullptr)
			delete cov_mat;
		// cov_mat = new Matrix<PRECISION>(unique_bipart, unique_bipart, (PRECISION)0);
		cov_mat = new SpecMat::LowerTri<PRECISION>(unique_bipart);
		cov_mat->form_row_ptr();

		Array<PRECISION> c_mean = sb2t_mat->col_mean();
		Array<PRECISION> c_sum = sb2t_mat->col_sum();

		for (int i = 0; i < unique_bipart; i++)
		{
			(*cov_mat)(i, i) += (n_trees * c_mean[i] * c_mean[i] - 2 * c_sum[i] * c_mean[i]);
			for (int j = 0; j < i; j++)
				(*cov_mat)(j, i) += (n_trees * c_mean[i] * c_mean[j] - c_sum[i] * c_mean[j] - c_sum[j] * c_mean[i]);
		}

		sb2t_mat->outer_prod_XXT_inplace((*cov_mat));

		for (int i = 0; i < unique_bipart; i++)
			for (int j = 0; j <= i; j++)
				(*cov_mat)(j, i) /= (n_trees - 1);
	};

	void Compute_Covariance_Matrix(const Array<size_t> &sub_bipart_id, Array<size_t> &bipart_id_mapping)
	{
		assert(sb2t_mat != nullptr);

		std::cout << "\tComputing sparse sub-matrix...\n";

		SparseMatrix *sb2t_sub_mat = sb2t_mat->subMat_row(sub_bipart_id, bipart_id_mapping);

		std::cout << "\tComputing covariance matrix of the sparse sub-matrix...\n";

		size_t sub_row = sb2t_sub_mat->get_row(), sub_col = sb2t_sub_mat->get_col();

		this->n_sub_bipart = sub_row;
		this->n_sub_trees = sub_col;

		if (cov_mat != nullptr)
			delete cov_mat;
		// cov_mat = new Matrix<PRECISION>(unique_bipart, unique_bipart, (PRECISION)0);
		cov_mat = new SpecMat::LowerTri<PRECISION>(sub_row);
		cov_mat->form_row_ptr();

		// Array<PRECISION> c_mean = sb2t_sub_mat->col_mean();
		// auto sum = sb2t_mat->col_sum();
		Array<PRECISION> c_sum = sb2t_sub_mat->col_sum();

		for (int i = 0; i < sub_row; i++)
		{
			// (*cov_mat)(i, i) += (sub_col * c_mean[i] * c_mean[i] - 2 * c_sum[i] * c_mean[i]);
			// for (int j = 0; j < i; j++)
			// 	(*cov_mat)(i, j) += (sub_col * c_mean[i] * c_mean[j] - c_sum[i] * c_mean[j] - c_sum[j] * c_mean[i]);
			for (int j = 0; j <= i; j++)
				(*cov_mat)(i, j) -= c_sum[i] * c_sum[j] / this->n_trees;
		}

		sb2t_sub_mat->outer_prod_XXT_inplace((*cov_mat));

		for (int i = 0; i < sub_row; i++)
			for (int j = 0; j <= i; j++)
				(*cov_mat)(i, j) /= (this->n_trees - 1);
		delete sb2t_sub_mat;
	};

	void Compute_RF_Distance_Matrix()
	{
		// Update formula:
		// 1) D[j, k] += abs(M[i, j] - M[i, k]), for i in 1:row.
		// 2) D[j, k] = (D[j, k] + D[k, j]) / 4.
		// For 2), the computation of 1) guarantee that D is symmetric, which simply makes 2) become D = D / 2.
		assert(sb2t_mat != nullptr);
		assert(n_trees > 0);
		if (dis_mat != nullptr)
			delete dis_mat;
		// dis_mat = new Matrix<PRECISION>(n_trees, n_trees, (PRECISION)0);
		dis_mat = new SpecMat::LowerTri<PRECISION>(n_trees);
		dis_mat->form_row_ptr();
		sb2t_mat->set_RCS();

		for (int i = 0; i < unique_bipart; i++)
		{
			// The i-th row in B2T matrix.
			auto col_ind_ptr = sb2t_mat->get_RCS_ind_c_ptr(i);
			auto CCS_ind_ptr = sb2t_mat->get_RCS_val_c_ptr(i);

			size_t RCS_row_size = (col_ind_ptr == nullptr) ? 0 : col_ind_ptr->get_size();

			// Handle 0 <= j < cur_col  with M[i, j] = 0 first.
			for (int j = 0; j < (*col_ind_ptr)[0]; j++)
			{
				for (int RCS_next_ind = 0; RCS_next_ind < RCS_row_size; RCS_next_ind++)
				{
					// Found all nonzero behind j.
					size_t next_col = (*col_ind_ptr)[RCS_next_ind];
					PRECISION next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
					(*dis_mat)(next_col, j) += next_val;
					//(*dis_mat)(j, next_col) += next_val;
				}
			}

			// Handle j = cur_col with M[i, j] != 0.
			for (int RCS_ind = 0; RCS_ind < RCS_row_size; RCS_ind++)
			{
				// RCS_ind pointing at nonzero entries on the i-th row.
				size_t cur_col = (*col_ind_ptr)[RCS_ind];
				PRECISION cur_val = sb2t_mat->CCS_ele(cur_col, (*CCS_ind_ptr)[RCS_ind]);

				size_t RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				size_t next_col = (*col_ind_ptr)[RCS_next_ind];
				PRECISION next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
				// Record the next nonzero. If this is the last nonzero in the row, record itself
				// which will not be triggered in later computation.

				for (int j = cur_col + 1; j < next_col; j++)
				{
					// Handle cur_col < j < next_col  with M[i, j] = 0 first.
					// For here, we only need M[i, k] != 0 with j < k.
					for (RCS_next_ind; RCS_next_ind < RCS_row_size; RCS_next_ind++)
					{
						// Found all nonzero behind j.
						next_col = (*col_ind_ptr)[RCS_next_ind];
						next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
						(*dis_mat)(next_col, j) += next_val;
						//(*dis_mat)(j, next_col) += next_val;
					}
				}

				RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				next_col = (*col_ind_ptr)[RCS_next_ind];
				next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
				// Restore the record of the next nonzero.

				// Handle j = cur_col
				// For here, we need to consider all j < k.
				for (int k = cur_col + 1; k < n_trees; k++)
				{
					// Update D[k, cur_col] with k > cur_col, i.e. the lower triangle part.
					// Do the same with the upper triangle part, i.e., D[cur_col, k]
					// D[k, cur_col] += abs(M[i, k] - M[i, cur_col]);

					if (k != next_col)
					{
						// M[i, k] = 0;
						(*dis_mat)(k, cur_col) += cur_val;
						//(*dis_mat)(cur_col, k) += cur_val;
					}
					else
					{
						// M[i, k] = next_val
						PRECISION temp = (cur_val > next_val) ? (cur_val - next_val) : (next_val - cur_val);
						(*dis_mat)(k, cur_col) += temp;
						//(*dis_mat)(cur_col, k) += temp;

						// update next nonzero
						RCS_next_ind = RCS_next_ind + (RCS_next_ind < RCS_row_size - 1);
						next_col = (*col_ind_ptr)[RCS_next_ind];
						next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
					}
				}
			}
		}
		if (!Trees->issameleaf())
		{
			for (auto j = 0; j < n_trees; j++)
				for (auto i = j + 1; i < n_trees; i++)
					(*dis_mat)(i, j) = (*dis_mat)(i, j) / 2.0;
		}
	}

	void Compute_RF_Distance_Matrix(const Array<size_t> &sub_tree_id, Array<size_t> &bipart_id_mapping)
	{
		// Update formula:
		// 1) D[j, k] += abs(M[i, j] - M[i, k]), for i in 1:row.
		// 2) D[j, k] = (D[j, k] + D[k, j]) / 4.
		// For 2), the computation of 1) guarantee that D is symmetric, which simply makes 2) become D = D / 2.
		assert(sb2t_mat != nullptr);
		assert(sub_tree_id.get_size() > 0);

		std::cout << "\tComputing sparse sub-matrix...\n";

		SparseMatrix *sb2t_sub_mat = sb2t_mat->subMat_col(sub_tree_id, bipart_id_mapping);

		std::cout << "\tComputing distance matrix of the sparse sub-matrix...\n";

		size_t sub_row = sb2t_sub_mat->get_row(), sub_col = sb2t_sub_mat->get_col();

		this->n_sub_bipart = sub_row;
		this->n_sub_trees = sub_col;

		if (dis_mat != nullptr)
			delete dis_mat;
		// dis_mat = new Matrix<PRECISION>(n_trees, n_trees, (PRECISION)0);
		dis_mat = new SpecMat::LowerTri<PRECISION>(sub_col);
		dis_mat->form_row_ptr();

		sb2t_sub_mat->set_RCS();

		for (int i = 0; i < sub_row; i++)
		{
			// The i-th row in B2T matrix.
			auto col_ind_ptr = sb2t_sub_mat->get_RCS_ind_c_ptr(i);
			auto CCS_ind_ptr = sb2t_sub_mat->get_RCS_val_c_ptr(i);

			size_t RCS_row_size = (col_ind_ptr == nullptr) ? 0 : col_ind_ptr->get_size();

			// Handle 0 <= j < cur_col  with M[i, j] = 0 first.
			for (int j = 0; j < (*col_ind_ptr)[0]; j++)
			{
				for (int RCS_next_ind = 0; RCS_next_ind < RCS_row_size; RCS_next_ind++)
				{
					// Found all nonzero behind j.
					size_t next_col = (*col_ind_ptr)[RCS_next_ind];
					PRECISION next_val = sb2t_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
					(*dis_mat)(next_col, j) += next_val;
					//(*dis_mat)(j, next_col) += next_val;
				}
			}

			// Handle j = cur_col with M[i, j] != 0.
			for (int RCS_ind = 0; RCS_ind < RCS_row_size; RCS_ind++)
			{
				// RCS_ind pointing at nonzero entries on the i-th row.
				size_t cur_col = (*col_ind_ptr)[RCS_ind];
				PRECISION cur_val = sb2t_sub_mat->CCS_ele(cur_col, (*CCS_ind_ptr)[RCS_ind]);

				size_t RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				size_t next_col = (*col_ind_ptr)[RCS_next_ind];
				PRECISION next_val = sb2t_sub_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
				// Record the next nonzero. If this is the last nonzero in the row, record itself
				// which will not be triggered in later computation.

				for (int j = cur_col + 1; j < next_col; j++)
				{
					// Handle cur_col < j < next_col  with M[i, j] = 0 first.
					// For here, we only need M[i, k] != 0 with j < k.
					for (RCS_next_ind; RCS_next_ind < RCS_row_size; RCS_next_ind++)
					{
						// Found all nonzero behind j.
						next_col = (*col_ind_ptr)[RCS_next_ind];
						next_val = sb2t_sub_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
						(*dis_mat)(next_col, j) += next_val;
						//(*dis_mat)(j, next_col) += next_val;
					}
				}

				RCS_next_ind = (RCS_ind < RCS_row_size - 1) ? RCS_ind + 1 : RCS_ind;
				next_col = (*col_ind_ptr)[RCS_next_ind];
				next_val = sb2t_sub_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
				// Restore the record of the next nonzero.

				// Handle j = cur_col
				// For here, we need to consider all j < k.
				for (int k = cur_col + 1; k < sub_col; k++)
				{
					// Update D[k, cur_col] with k > cur_col, i.e. the lower triangle part.
					// Do the same with the upper triangle part, i.e., D[cur_col, k]
					// D[k, cur_col] += abs(M[i, k] - M[i, cur_col]);

					if (k != next_col)
					{
						// M[i, k] = 0;
						(*dis_mat)(k, cur_col) += cur_val;
						//(*dis_mat)(cur_col, k) += cur_val;
					}
					else
					{
						// M[i, k] = next_val
						PRECISION temp = (cur_val > next_val) ? (cur_val - next_val) : (next_val - cur_val);
						(*dis_mat)(k, cur_col) += temp;
						//(*dis_mat)(cur_col, k) += temp;

						// update next nonzero
						RCS_next_ind = RCS_next_ind + (RCS_next_ind < RCS_row_size - 1);
						next_col = (*col_ind_ptr)[RCS_next_ind];
						next_val = sb2t_sub_mat->CCS_ele(next_col, (*CCS_ind_ptr)[RCS_next_ind]);
					}
				}
			}
		}

		if (!Trees->issameleaf())
		{
			for (auto j = 0; j < sub_col; j++)
				for (auto i = j + 1; i < sub_col; i++)
					(*dis_mat)(i, j) = (*dis_mat)(i, j) / 2.0;
		}
	};

	int get_all_bipart_num() { return all_bipart; };
	int get_unique_bipart_num() { return unique_bipart; };
	int get_tree_set_size() { return n_trees; };
	int get_leaf_set_size() { return n_leaves; };
	auto get_sb2t() const {return this->sb2t_mat;}
	SpecMat::LowerTri<PRECISION> *get_dis_mat() { return dis_mat; };
	SpecMat::LowerTri<PRECISION> *get_cov_mat() { return cov_mat; };
	TreeSet<T> *pop_trees()
	{
		auto temp = Trees;
		Trees = nullptr;
		return temp;
	}

	void print_Bipart_List(ostream &fout)
	{

		fout << n_leaves << '\t' << unique_bipart << '\n';
		for (int i = 0; i < unique_bipart; i++)
		{
			fout << i << " ";
			(*Bipart)[i].print_BitString(fout);
			fout << ' ' << (*sb2t_mat).get_RCS_ind_c_ptr(i)->get_size() << '\n';
		}
	};

	void print_Bipart_List(ostream &fout, const Array<int> &indices)
	{
		int tree_num = sb2t_mat->get_col();
		int ind = -1;
		for (int i = 0; i < indices.get_size(); i++)
		{
			ind = indices(i);
			(*Bipart)[ind].print_BitString(fout);
			fout << ' ' << ((PRECISION)sb2t_mat->get_RCS_ind_c_ptr(ind)->get_size()) / tree_num << '\n';
		}
	};

	void print_Bipart_List(ostream &fout, int *indices, int num_bipart)
	{
		int tree_num = sb2t_mat->get_col();
		int ind = -1;
		for (int i = 0; i < num_bipart; i++)
		{
			ind = indices[i];
			(*Bipart)[ind].print_BitString(fout);
			fout << ' ' << ((PRECISION)sb2t_mat->get_RCS_ind_c_ptr(ind)->get_size()) / tree_num << '\n';
		}
	};

	void print_Bipart2Tree_Matrix(ostream &fout, SparseMatrixOutputType smtype) { (*sb2t_mat).print(fout, smtype); };

	void print_Covariance_Matrix(ostream &fout) { (*cov_mat).print(fout); };

	void print_Distance_Matrix(ostream &fout) { (*dis_mat).print(fout); };

	void print_summary()
	{
		int max_collision_size = 0, collision_cnt_5 = 0, collision_cnt_10 = 0, empty_cnt = 0;
		int cur_size;
		auto Hash2Id = Bipart->get_hash_table();
		auto hash_bound = Bipart->get_hash_bound();
		for (int i = 0; i < hash_bound; i++)
		{
			if (Hash2Id->get_c_ptr(i) != nullptr)
			{
				cur_size = Hash2Id->get_c_ptr(i)->get_size();
				if (cur_size > max_collision_size)
					max_collision_size = cur_size;
				if (cur_size >= 5)
					collision_cnt_5++;
				if (cur_size >= 10)
					collision_cnt_10++;
				if (cur_size == 0)
					empty_cnt++;
			}
			else
			{
				empty_cnt++;
			}
		}
		std::cout << "-------------Infomation summary------------------\n";
		std::cout << "\tTree set size:\t\t\t\t" << n_trees << ",\n";
		std::cout << "\tLeaf set size:\t\t\t\t" << n_leaves << ",\n";
		std::cout << "\tAll bipartitions encounted:\t\t" << all_bipart << ",\n";
		std::cout << "\tNumber of Distinct bipartition:\t\t" << unique_bipart << ",\n";
		std::cout << "\tNumber of computed bipartitions:\t" << n_sub_bipart << ",\n";
		std::cout << "\tNumber of computed trees:\t\t" << n_sub_trees << ".\n\n";
		std::cout << "\tHash container size:\t\t\t" << hash_bound << ",\n";
		std::cout << "\tEmpty container:\t\t\t" << empty_cnt << ",\n";
		std::cout << "\tMaximum Collision beam size:\t\t" << max_collision_size << ",\n";
		std::cout << "\tCount of Collision beams over 5:\t" << collision_cnt_5 << ",\n";
		std::cout << "\tCount of Collision beams over 10:\t" << collision_cnt_10 << ".\n";
		std::cout << "-------------Infomation summary------------------\n";
	}

	bool save_treeobj(String fname, bool B2T, bool DIS, bool COV)
	{

		// This function output following informations in given format to a binary file. (* is the bool indicator, 0 for nothing and 1 followed with the contents.)
		// n_tree n_leaves all_bipart unique_bipart ISWEIGHTED ISROOTED
		// * sb2t_mat;
		// * cov_mat;
		// * dis_mat;

		std::ofstream fout;
		fout.open(fname, std::ios::out | std::ios::binary);

		fout.write((char *)&n_trees, sizeof(int));
		fout.write((char *)&n_leaves, sizeof(int));
		fout.write((char *)&all_bipart, sizeof(int));
		fout.write((char *)&unique_bipart, sizeof(int));
		fout.write((char *)&ISWEIGHTED, sizeof(bool));
		fout.write((char *)&ISROOTED, sizeof(bool));

		fout.write((char *)&B2T, sizeof(bool));
		if (B2T)
			binary_print_smat(fout, sb2t_mat);

		fout.write((char *)&DIS, sizeof(bool));
		if (DIS && dis_mat != nullptr)
			binary_print_lowertri(fout, dis_mat);

		fout.write((char *)&COV, sizeof(bool));
		if (COV && cov_mat != nullptr)
			binary_print_lowertri(fout, cov_mat);

		fout.close();

		return true;
	}

	bool read_treeobj(String fname)
	{
		std::ifstream fin;
		fin.open(fname, std::ios::in | std::ios::binary);
		if (!fin.is_open())
		{
			std::cout << "Error! Cannot reload the binary file\n";
			throw(1);
		}

		bool B2T, DIS, COV;

		fin.read((char *)&n_trees, sizeof(int));
		fin.read((char *)&n_leaves, sizeof(int));
		fin.read((char *)&all_bipart, sizeof(int));
		fin.read((char *)&unique_bipart, sizeof(int));
		fin.read((char *)&ISWEIGHTED, sizeof(bool));
		fin.read((char *)&ISROOTED, sizeof(bool));

		fin.read((char *)&B2T, sizeof(bool));
		if (B2T)
		{
			if (sb2t_mat != nullptr)
				delete sb2t_mat;
			sb2t_mat = binary_read_smat(fin);
		}

		// fin.read((char *)&DIS, sizeof(bool));
		// if (DIS)
		// {
		// 	if (dis_mat != nullptr)
		// 		delete dis_mat;
		// 	dis_mat = binary_read_lowertri(fin);
		// }

		// fin.read((char *)&COV, sizeof(bool));
		// if (COV)
		// {
		// 	if (cov_mat != nullptr)
		// 		delete cov_mat;
		// 	cov_mat = binary_read_lowertri(fin);
		// }

		fin.close();
		return true;
	}

	bool save_treeset(String fname)
	{
		std::ofstream fout;
		fout.open(fname, std::ios::out | std::ios::binary);

		binary_print_taxon(fout, Taxa);
		binary_print_bipart(fout, Bipart);
		fout.close();
		return true;
	}

	bool read_treeset(String fname, T dummy)
	{
		std::ifstream fin;
		fin.open(fname, std::ios::in | std::ios::binary);

		if (Trees != nullptr)
			delete Trees;
		if (Bipart != nullptr)
			delete Bipart;
		if (Taxa != nullptr)
			delete Taxa;

		Taxa = binary_read_taxon(fin);
		Bipart = binary_read_bipart(fin, dummy, Taxa);

		return true;
	}

	void bipart_frequency_check(double lf, double hf, Array<size_t> &bipart_id)
	{
		assert(sb2t_mat != nullptr);

		sb2t_mat->set_RCS();

		double check = 0;

		for (int i = 0; i < unique_bipart; i++)
		{
			check = sb2t_mat->get_RCS_ind_c(i).get_size() / (double)n_trees;
			if (lf < check && check < hf)
				bipart_id.push(i);
		}

		bipart_id.align();
	}
};
