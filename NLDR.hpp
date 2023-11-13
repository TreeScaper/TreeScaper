#pragma once

#include "SpecMat.hpp"
#include "optim.hpp"
#include <random>

#ifndef INTEGER_TYPE
#define INTEGER_TYPE int
#endif

using SpecMat::LowerTri;

enum NLDR_STRESS_TYPE
{
	KRUSKAL1,
	NORMALIZED,
	SAMMON,
	CCA
};

namespace NLDR
{
	class NLDR_Var
	{
	public:
		Matrix<PRECISION> *Cor;
		LowerTri<PRECISION> *Dis;
		PRECISION *eta;
		PRECISION *eta_d;
		PRECISION *rho;

		NLDR_Var() : Cor(nullptr), Dis(nullptr), eta(nullptr), eta_d(nullptr), rho(nullptr){};
		NLDR_Var(Matrix<PRECISION> *C, LowerTri<PRECISION> *D, PRECISION *e, PRECISION *e_d, PRECISION *r) : Cor(C), Dis(D), eta(e), eta_d(e_d), rho(r){};
	};

	class NLDR_Var_Diff
	{
	public:
		INTEGER_TYPE activated_row;
		INTEGER_TYPE activated_coo;
		PRECISION eta;
		PRECISION rho;
		PRECISION *newD_row;

		NLDR_Var_Diff(INTEGER_TYPE i, INTEGER_TYPE k, PRECISION e, PRECISION r, PRECISION *nD_row)
			: activated_row(i), activated_coo(k), eta(e), rho(r), newD_row(nD_row){};
	};

	class NLDR_Dir
	{
	public:
		Matrix<PRECISION> *Delta;
		// For general direction, Delta has the same dimension with Coordiante matrix, n_pts \times dim.
		// For Gauss-Seidel method, Delta will be of size 1 \times dim, (a row).
		INTEGER_TYPE activated_row;
		INTEGER_TYPE activated_coo;

		NLDR_Dir(Matrix<PRECISION> *D, INTEGER_TYPE a_r, INTEGER_TYPE a_c) : Delta(D), activated_row(a_r), activated_coo(a_c){};
	};

	class NLDR_Paras
	{
		// weight = 2 eta_d
	public:
		LowerTri<PRECISION> *Dis;
		LowerTri<PRECISION> *BX;
		Matrix<PRECISION> *V;
		INTEGER_TYPE dim;
		INTEGER_TYPE pts;
		INTEGER_TYPE nk;
		NLDR_Var_Diff *VarDiff;
		PRECISION *lambda;
		PRECISION cur_lambda;
		PRECISION pre_stress;
		PRECISION total_d;
		PRECISION eta_d;
		PRECISION T;
		NLDR_STRESS_TYPE nst;
		INTEGER_TYPE *iter;
		INTEGER_TYPE *s_ind;

		Array<PRECISION> *row_sum;
		Array<PRECISION> *A;
		Matrix<PRECISION> *CX_temp;

		Matrix<PRECISION> *perturb;
		std::default_random_engine *eng;

		NLDR_Paras() : Dis(nullptr), dim(0){};
		NLDR_Paras(LowerTri<PRECISION> *D, LowerTri<PRECISION> *matBX, Matrix<PRECISION> *matV, INTEGER_TYPE d, NLDR_Var_Diff *VD, PRECISION *l, PRECISION ps)
			: Dis(D), BX(matBX), V(matV), dim(d), VarDiff(VD), lambda(l), pre_stress(ps) { this->compute_total_eta(); };

		void compute_total_eta()
		{
			total_d = 0;
			eta_d = 0;
			auto n = Dis->dimension();
			for (auto i = 0; i < n; i++)
			{
				auto D_i = (*Dis)[i];
				for (auto j = 0; j < i; j++)
				{
					total_d += D_i[j];
					eta_d += D_i[j] * D_i[j];
				}
			}
		};

		void set_stress_type(NLDR_STRESS_TYPE st) { nst = st; };
	};

	class NLDR_Paras_GAU
	{
	public:
		INTEGER_TYPE pts; //n
		INTEGER_TYPE dim; //d
		LowerTri<PRECISION> *Dis;
		NLDR_Var_Diff *VarDiff;
		PRECISION total_d;
		PRECISION eta_d;
		Matrix<PRECISION> *XD;
		Matrix<PRECISION> *XD2;
		Array<PRECISION> *DTX;
		Array<PRECISION> *DTX3;
		Array<PRECISION> *DR;
		PRECISION *lambda;
		PRECISION *cur_lambda;
		OptimLib::Wolfe_1st_StepSize *wolfe;

		NLDR_Paras_GAU(LowerTri<PRECISION> *dis, INTEGER_TYPE d)
		{
			pts = dis->dimension();
			dim = d;
			Dis = dis;
			VarDiff = nullptr;
			lambda = nullptr;
			total_d = 0;
			eta_d = 0;
			XD = nullptr;
			XD2 = nullptr;
			DTX = nullptr;
			DTX3 = nullptr;

			for (auto i = 0; i < pts; i++)
			{
				auto D_i = (*Dis)[i];
				for (auto j = 0; j < i; j++)
				{
					total_d += D_i[j];
					eta_d += D_i[j] * D_i[j];
				}
			}
		}

	};

	void NLDR_Var_Copy(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	PRECISION stress_raw(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void stress(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void Cor2Dis(Matrix<PRECISION> &C, LowerTri<PRECISION> &D);

	void Cor2Dis_with_paras(NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void Cor2Dis_no_paras(NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void update_by_index(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_distance_row_from_VarDiff(NLDR_Var *X, NLDR_Var *Y, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras);

	void compute_BX(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void compute_V(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void update_cor_by_BX(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	void update_cor_by_row(NLDR_Var *X, NLDR_Dir *Dir);

	void update_coordinates_inplace(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_coordinates_notinplace(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_Var_inplace(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void Init_X(NLDR::NLDR_Var *X, PRECISION *stressX, NLDR::NLDR_Paras *paras,
				PRECISION (*Cost)(NLDR::NLDR_Var *, NLDR::NLDR_Paras *), LowerTri<PRECISION> *Dis_0);

	void Init_X(NLDR::NLDR_Var *X, LowerTri<PRECISION> *DX0, LowerTri<PRECISION> *D0);

	void update_MAJORIZATION(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	bool accept_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	bool loop_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws);

	void direction_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void slope_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_GAUSS_SEIDEL_iteration(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void copy_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void MAJORIZATION(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
					  NLDR_STRESS_TYPE stress_type, INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol,
					  PRECISION *lambda);

	void STOCHASTIC(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
					NLDR_STRESS_TYPE stress_type, INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol,
					PRECISION *alpha, PRECISION *lambda, INTEGER_TYPE *s_ind);

	void METROPOLIS(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
					NLDR_STRESS_TYPE stress_type, INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol,
					PRECISION *lambda, std::default_random_engine *eng, PRECISION T);

	void GAUSS_SEIDEL(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
					  NLDR_STRESS_TYPE stress_type, INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol,
					  PRECISION *lambda);

	/****************************************KRUSKAL1****************************************/

	void KRUSKAL1_distance_stress_function(NLDR::NLDR_Var *X, PRECISION *f, NLDR::NLDR_Paras *Paras);

	// Compute stress from distance matrix
	PRECISION KRUSKAL1_stress_by_distance(const LowerTri<PRECISION> &DX, const LowerTri<PRECISION> &D);

	// Compute stress from parameters eta, eta_d and rho
	PRECISION KRUSKAL1_stress_by_parameter(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *paras);

	PRECISION KRUSKAL1_stress_by_parameter(NLDR_Var *X);

	void KRUSKAL1_rescale(NLDR::NLDR_Var *X);

	void compute_BX_KRUSKAL1(NLDR_Var *X, NLDR_Paras *Paras);

	void update_cor_by_BX_KRUSKAL1(NLDR_Var *X, NLDR_Paras *Paras);

	void update_KRUSKAL1_MAJORIZATION(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void KRUSKAL1_MAJORIZATION(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
							   INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol);

	void direction_KRUSKAL1_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_KRUSKAL1_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);



	void direction_KRUSKAL1_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_KRUSKAL1_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_KRUSKAL1_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void pre_update_by_index_KRUSKAL1_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras);

	void direction_KRUSKAL1_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, INTEGER_TYPE a_row,
										 PRECISION **XD_, PRECISION **XD2_, PRECISION *DTX_, PRECISION *DTX3_);

	void update_KRUSKAL1_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	/****************************************KRUSKAL1****************************************/

	/****************************************NORMALIZED**************************************/

	void NORMALIZED_distance_stress_function(NLDR::NLDR_Var *X, PRECISION *f, NLDR::NLDR_Paras *Paras);

	PRECISION NORMALIZED_stress_by_distance(const LowerTri<PRECISION> &DX, const LowerTri<PRECISION> &D);

	PRECISION NORMALIZED_stress_by_parameter(NLDR::NLDR_Var *X, NLDR_Paras *paras);

	PRECISION NORMALIZED_stress_by_parameter(NLDR_Var *X);

	void compute_BX_NORMALIZED(NLDR_Var *X, NLDR_Paras *Paras);

	void update_cor_by_BX_NORMALIZED(NLDR_Var *X, NLDR_Paras *Paras);

	void update_NORMALIZED_MAJORIZATION(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void NORMALIZED_MAJORIZATION(LowerTri<PRECISION> *Dis, INTEGER_TYPE n, INTEGER_TYPE d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
								 INTEGER_TYPE MaxIter, INTEGER_TYPE MaxTime, PRECISION ErrTol, PRECISION DifTol);

	void direction_NORMALIZED_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_NORMALIZED_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_NORMALIZED_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_NORMALIZED_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void pre_update_by_index_NORMALIZED_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras);

	void update_NORMALIZED_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_NORMALIZED_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_NORMALIZED_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, INTEGER_TYPE a_row,
        								   PRECISION** XD_, PRECISION** XD2_, PRECISION* DTX_, PRECISION* DTX3_);

	void update_NORMALIZED_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	/****************************************NORMALIZED**************************************/
	/****************************************SAMMON******************************************/

	PRECISION SAMMON_stress_by_distance(NLDR_Var *X, NLDR::NLDR_Paras *Paras);

	PRECISION SAMMON_stress_by_distance(NLDR_Var *X, LowerTri<PRECISION> *D, PRECISION total_d);

	void SAMMON_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras *Paras);

	void SAMMON_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras_GAU *Paras);

	void compute_BX_SAMMON(NLDR_Var *X, NLDR_Paras *Paras);

	void compute_V_SAMMON(NLDR_Var *X, NLDR_Paras *Paras);

	void update_cor_by_BX_SAMMON(NLDR_Var *X, NLDR_Paras *Paras);

	void update_SAMMON_MAJORIZATION(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_SAMMON_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_SAMMON_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_SAMMON_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_SAMMON_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void pre_update_by_index_SAMMON_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras);

	void update_SAMMON_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_SAMMON_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_SAMMON_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, INTEGER_TYPE a_row,
        							   PRECISION** XD_, PRECISION** XD2_, PRECISION* DTX_, PRECISION* DTX3_);
	
	void update_SAMMON_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	/****************************************SAMMON******************************************/
	/****************************************CCA*********************************************/

	PRECISION compute_lambda(NLDR_Var *X, NLDR_Paras *Paras);

	PRECISION CCA_stress_by_distance(NLDR_Var *X, NLDR_Paras *Paras);

	PRECISION CCA_stress_by_distance_lambda(LowerTri<PRECISION> *DX, LowerTri<PRECISION> *D, PRECISION lambda);

	void CCA_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras *Paras);

	void CCA_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras_GAU *Paras);

	void compute_BX_CCA(NLDR_Var *X, NLDR_Paras *Paras);

	void compute_V_CCA(NLDR_Var *X, NLDR_Paras *Paras);

	void update_cor_by_BX_CCA(NLDR_Var *X, NLDR_Paras *Paras);

	void update_CCA_MAJORIZATION(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_CCA_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void update_CCA_STOCHASTIC(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_CCA_METROPOLIS(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_CCA_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void pre_update_by_index_CCA_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras);

	void update_CCA_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void cost_CCA_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	void direction_CCA_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, INTEGER_TYPE a_row,
        							PRECISION** XD_, PRECISION** XD2_, PRECISION* DPJ_, PRECISION* DQ_, PRECISION* DR_);

	void update_CCA_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, INTEGER_TYPE *ind);

	/****************************************CCA*********************************************/
}