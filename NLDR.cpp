#include "NLDR.hpp"

void NLDR::NLDR_Var_Copy(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *src = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *f_src = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));
    NLDR_Var *des = reinterpret_cast<NLDR_Var *>(kws.arg(ind[2]));
    PRECISION *f_des = reinterpret_cast<PRECISION *>(kws.arg(ind[3]));

    *(des->eta) = *(src->eta);
    *(des->rho) = *(src->rho);

    unsigned n = des->Cor->get_row();
    unsigned d = des->Cor->get_col();
    PRECISION **des_Cor_ptr = des->Cor->get_vec();
    PRECISION **src_Cor_ptr = src->Cor->get_vec();
    PRECISION *des_Dis_ptr = des->Dis->get_vec();
    PRECISION *src_Dis_ptr = src->Dis->get_vec();
    unsigned n_dis = n * (n + 1) / 2;
    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < d; k++)
            des_Cor_ptr[i][k] = src_Cor_ptr[i][k];

    for (auto i = 0; i < n_dis; i++)
        des_Dis_ptr[i] = src_Dis_ptr[i];

    *f_des = *f_src;
};

PRECISION NLDR::stress_raw(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras)
{
    switch (Paras->nst)
    {
    case KRUSKAL1:
        return KRUSKAL1_stress_by_parameter(X, Paras);
    case NORMALIZED:
        return NORMALIZED_stress_by_parameter(X, Paras);
    case SAMMON:
        return SAMMON_stress_by_distance(X, Paras);
    case CCA:
        return CCA_stress_by_distance(X, Paras);
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }
}

void NLDR::stress(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *f = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    *f = stress_raw(X, paras);
}

void NLDR::Cor2Dis(Matrix<PRECISION> &C, LowerTri<PRECISION> &D)
{
    unsigned n = C.get_row();
    unsigned dim = C.get_col();
    PRECISION temp;

    auto C_ = C.get_vec();

    for (auto i = 0; i < n; i++)
    {
        auto D_i = D[i];
        for (auto j = 0; j < i; j++)
        {
            D_i[j] = 0;
            for (auto k = 0; k < dim; k++)
            {
                temp = C_[i][k] - C_[j][k];
                D_i[j] += temp * temp;
            }
            D_i[j] = sqrt(D_i[j]);
        }
        D_i[i] = 0;
    }
}

void NLDR::Cor2Dis_with_paras(NLDR_Var *X, NLDR::NLDR_Paras *Paras)
{
    Matrix<PRECISION> &CX = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    unsigned n = CX.get_row();
    unsigned dim = CX.get_col();
    PRECISION temp = 0;
    PRECISION &eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION &rho = *(X->rho);
    eta = 0;
    rho = 0;

    auto CX_ = CX.get_vec();

    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        for (auto j = 0; j < i; j++)
        {
            DX_i[j] = 0;
            for (auto k = 0; k < dim; k++)
            {
                temp = CX_[i][k] - CX_[j][k];
                DX_i[j] += temp * temp;
            }
            DX_i[j] = sqrt(DX_i[j]);
            eta += DX_i[j] * DX_i[j];
            rho += DX_i[j] * D_i[j];
        }
        DX_i[i] = 0;
    }
}

void NLDR::Cor2Dis_no_paras(NLDR_Var *X, NLDR::NLDR_Paras *Paras)
{
    Matrix<PRECISION> &CX = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    unsigned n = CX.get_row();
    unsigned dim = CX.get_col();
    PRECISION temp = 0;

    auto CX_ = CX.get_vec();

    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        for (auto j = 0; j < i; j++)
        {
            DX_i[j] = 0;
            for (auto k = 0; k < dim; k++)
            {
                temp = CX_[i][k] - CX_[j][k];
                DX_i[j] += temp * temp;
            }
            DX_i[j] = sqrt(DX_i[j]);
        }
        DX_i[i] = 0;
    }
}

// void NLDR::pre_update_by_index(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
// {
//     Matrix<PRECISION> &C = *(X->Cor);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     Matrix<PRECISION> &Delta = *(Dir->Delta);

//     unsigned dim = C.get_col(), n = C.get_row();
//     unsigned i = Dir->activated_row;
//     unsigned k = Dir->activated_coo;
//     PRECISION e1 = Delta(0, k);
//     PRECISION e2 = Delta(1, k);
//     PRECISION xcnew = C(i, k) - stepsize * e1 / e2;

//     PRECISION *newD_i = Paras->VarDiff->newD_row;
//     Paras->VarDiff->activated_row = i;
//     Paras->VarDiff->activated_coo = k;
//     PRECISION eta = *(X->eta);
//     PRECISION eta_d = *(X->eta_d);
//     PRECISION rho = *(X->rho);

//     PRECISION DX_ij, D_ij;
//     auto C_ = C.get_vec();

//     PRECISION temp = 0;
//     if (st == KRUSKAL1 || st == NORMALIZED)
//         for (auto j = 0; j < n; j++)
//         {
//             DX_ij = DX(i, j);
//             D_ij = D(i, j);
//             eta -= DX_ij * DX_ij;
//             rho -= DX_ij * D_ij;
//             newD_i[j] = 0;
//             for (auto l = 0; l < dim; l++)
//             {
//                 if (l != k)
//                     temp = C_[i][l]- C_[j][l];
//                 else
//                 {
//                     if (i != j)
//                         temp = xcnew - C_[j][k];
//                     else
//                         temp = 0;
//                 }
//                 newD_i[j] += temp * temp;
//             }
//             newD_i[j] = sqrt(newD_i[j]);
//             eta += newD_i[j] * newD_i[j];
//             rho += newD_i[j] * D_ij;
//         }
//     else
//         for (auto j = 0; j < n; j++)
//         {
//             newD_i[j] = 0;
//             for (auto l = 0; l < dim; l++)
//             {
//                 if (l != k)
//                     temp = C_[i][k]- C_[j][k];
//                 else
//                     temp = xcnew - C_[j][k];
//                 newD_i[j] += temp * temp;
//             }
//             newD_i[j] = sqrt(newD_i[j]);
//         }
//     newD_i[i] = 0;
//     Paras->VarDiff->eta = eta;
//     Paras->VarDiff->rho = rho;
// }

void NLDR::update_distance_row_from_VarDiff(NLDR_Var *X, NLDR_Var *Y, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
{
    unsigned i = Paras->VarDiff->activated_row;
    unsigned k = Paras->VarDiff->activated_coo;
    unsigned n = X->Cor->get_row();
    unsigned d = X->Cor->get_col();

    Matrix<PRECISION> &CY_ = *Y->Cor;
    LowerTri<PRECISION> &DY_ = *Y->Dis;
    Matrix<PRECISION> &Delta_ = *Dir->Delta;

    CY_(i, k) = (*X->Cor)(i, k) - stepsize * Delta_(0, k) / Delta_(1, k);
    auto newD = Paras->VarDiff->newD_row;
    for (auto j = 0; j < i; j++)
        DY_(i, j) = newD[j];
    DY_(i, i) = 0;
    for (auto j = i + 1; j < n; j++)
        DY_(j, i) = newD[j];

    *Y->eta = Paras->VarDiff->eta;
    *Y->eta_d = Paras->eta_d;
    *Y->rho = Paras->VarDiff->rho;
}

// void NLDR::compute_V(NLDR_Var *X, NLDR_Paras *Paras)
// {
//     unsigned n;
//     LowerTri<PRECISION> *D;
//     Matrix<PRECISION> *V;
//     LowerTri<PRECISION> *DX;

//     PRECISION ** V_;
//     switch (Paras->nst)
//     {
//     case SAMMON:
//         D = Paras->Dis;
//         V = Paras->V;
//         V_ = V->get_vec();
//         n = D->dimension();
//         for (auto i = 0; i < n; i++)
//         {
//             auto D_i = (*D)[i];
//             for (auto j = 0; j < i; j++)
//             {
//                 if (D_i[j] < 1e-5)
//                     V_[i][j] = -1e5;
//                 else
//                     V_[i][j] = -1.0 / D_i[j];
//                 V_[j][i] = V_[i][j];
//                 V_[i][i] -= V_[i][j];
//                 V_[j][j] -= V_[j][i];
//             }
//         }
//         shift(*V, 1.0);
//         inverse(*V);
//         shift(*V, -1.0 / n / n);
//         break;
//     case CCA:
//         if (!X)
//             break;
//         D = Paras->Dis;
//         DX = X->Dis;
//         V = Paras->V;
//         V_ = V->get_vec();
//         n = D->dimension();
//         for (auto i = 0; i < n; i++)
//         {
//             auto D_i = (*D)[i];
//             auto DX_i = (*DX)[i];
//             PRECISION c_l = Paras->lambda[*Paras->iter];
//             for (auto j = 0; j < i; j++)
//             {
//                 V_[i][j] = -exp(-DX_i[j] / c_l);
//                 V_[j][i] = V_[i][j];
//                 V_[i][i] -= V_[i][j];
//                 V_[j][j] -= V_[j][i];
//             }
//         }
//         shift(*V, 1.0);
//         inverse(*V);
//         shift(*V, -1.0 / n / n);
//         break;
//     default:
//         break;
//     }
// };

// void NLDR::compute_BX(NLDR_Var *X, NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     LowerTri<PRECISION> &BX = *(Paras->BX);
//     PRECISION *row_sum = new PRECISION[n];
//     memset(row_sum, 0, n * sizeof(PRECISION));
//     for (auto i = 0; i < n; i++)
//     {
//         auto DX_i = DX[i];
//         auto D_i = D[i];
//         auto BX_i = BX[i];
//         for (auto j = 0; j < i; j++)
//             if (DX_i[j] != 0)
//             {
//                 switch (Paras->nst)
//                 {
//                 case KRUSKAL1:
//                     BX_i[j] = -D_i[j] / DX_i[j];
//                     break;
//                 case NORMALIZED:
//                     BX_i[j] = -D_i[j] / DX_i[j];
//                     break;
//                 case SAMMON:
//                     BX_i[j] = -1.0 / DX_i[j];
//                     break;
//                 case CCA:
//                     BX_i[j] = -D_i[j] / DX_i[j] * exp(-DX_i[j] / Paras->lambda[*Paras->iter]);
//                     break;
//                 default:
//                     std::cout << "Cost function type in NLDR_Paras not found!\n";
//                     throw(1);
//                 }

//                 row_sum[i] += BX_i[j];
//                 row_sum[j] += BX_i[j];
//             }
//             else
//                 BX_i[j] = 0;
//     }

//     for (auto i = 0; i < n; i++)
//         BX(i, i) = -row_sum[i];
// }

// void NLDR::update_cor_by_BX(NLDR_Var *X, NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     LowerTri<PRECISION> &BX = *(Paras->BX);
//     PRECISION **BX_ = new PRECISION*[n];
//     for(auto i = 0; i < n; i++)
//         BX_[i] = BX[i];
//     auto XC_ = Cor.get_vec();

//     if (Paras->nst == CCA || Paras->nst == SAMMON)
//     {
//         Matrix<PRECISION> &V = *(Paras->V);
//         auto V_ = V.get_vec();
//         PRECISION *A = new PRECISION[n * dim];
//         memset(A, 0, n * dim * sizeof(PRECISION));
//         for (auto i = 0; i < n; i++)
//             for (auto k = 0; k < dim; k++){
//                 for (auto j = 0; j <= i; j++)
//                     A[i * dim + k] += BX_[i][j] * XC_[j][k];
//                 for (auto j = i + 1; j < n; j++)
//                     A[i * dim + k] += BX_[j][i] * XC_[j][k];
//                 // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
//             }

//         for (auto i = 0; i < n; i++)
//             for (auto k = 0; k < dim; k++)
//             {
//                 XC_[i][k] = 0;
//                 for (auto j = 0; j < n; j++)
//                     XC_[i][k] += V_[i][j] * A[j * dim + k];
//             }
//         delete[] A;
//     }
//     else
//     {
//         Matrix<PRECISION> C_temp(Cor);
//         auto C_temp_ = C_temp.get_vec();
//         for (auto i = 0; i < n; i++)
//             for (auto k = 0; k < dim; k++)
//             {
//                 XC_[i][k] = 0;
//                 for (auto j = 0; j <= i; j++)
//                     XC_[i][k] += BX_[i][j] * C_temp_[j][k];
//                 for (auto j = i + 1; j < n; j++)
//                     XC_[i][k] += BX_[j][i] * C_temp_[j][k];
//                 // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
//                 XC_[i][k] /= n;
//             }
//     }
//     delete[] BX_;
// }

void NLDR::update_coordinates_inplace(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    PRECISION t = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));

    Matrix<PRECISION> &CX = *(X->Cor);
    Matrix<PRECISION> &DeltaX = *(Dir->Delta);

    shift(CX, DeltaX, t);
}

void NLDR::update_coordinates_notinplace(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    PRECISION t = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    NLDR_Var *Y = reinterpret_cast<NLDR_Var *>(kws.arg(ind[3]));

    unsigned forward_cp_ind[4] = {ind[0], ind[4], ind[3], ind[5]};
    NLDR_Var_Copy(kws, forward_cp_ind);

    Matrix<PRECISION> &CY = *(Y->Cor);
    Matrix<PRECISION> &DeltaY = *(Dir->Delta);

    shift(CY, DeltaY, t);
}

void NLDR::update_Var_inplace(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    update_coordinates_inplace(kws, ind);
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    Cor2Dis_with_paras(X, paras);
}

void NLDR::Init_X(NLDR::NLDR_Var *X, PRECISION *stressX, NLDR::NLDR_Paras *paras,
                  PRECISION (*Cost)(NLDR::NLDR_Var *, NLDR::NLDR_Paras *), LowerTri<PRECISION> *Dis_0)
{

    std::cout << "\t\tCompute components in X...\n";
    PRECISION &eta = *(X->eta);
    PRECISION &eta_d = *(X->eta_d);
    PRECISION &rho = *(X->rho);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);
    eta = 0;
    eta_d = 0;
    rho = 0;
    unsigned n = D.dimension();
    PRECISION total = 0;

    if (Dis_0 && Dis_0->dimension() != 0)
    {
        memcpy(DX.get_vec(), Dis_0->get_vec(), (n * (n + 1) / 2) * sizeof(PRECISION));
        for (auto i = 0; i < n; i++)
        {
            auto DX_i = DX[i];
            auto D_i = D[i];
            for (auto j = 0; j < i; j++)
            {
                eta += DX_i[j] * DX_i[j];
                eta_d += D_i[j] * D_i[j];
                rho += D_i[j] * DX_i[j];
            }
        }
    }
    else
        Cor2Dis_with_paras(X, paras);
    for (auto i = 0; i < n; i++)
    {
        auto D_i = D[i];
        for (auto j = 0; j < i; j++)
            total += D_i[j];
    }
    paras->total_d = total;
    paras->eta_d = eta_d;
    *stressX = Cost(X, paras);
}

void NLDR::Init_X(NLDR::NLDR_Var *X, LowerTri<PRECISION> *DX0, LowerTri<PRECISION> *D0)
{

    std::cout << "\t\tCompute components in X...\n";
    PRECISION &eta = *(X->eta);
    PRECISION &eta_d = *(X->eta_d);
    PRECISION &rho = *(X->rho);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *D0;
    eta = 0;
    eta_d = 0;
    rho = 0;
    unsigned n = D.dimension();
    PRECISION total = 0;

    if (DX0 && DX0->dimension() != 0)
        memcpy(DX.get_vec(), DX0->get_vec(), (n * (n + 1) / 2) * sizeof(PRECISION));
    else
        Cor2Dis(*X->Cor, *X->Dis);

    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        for (auto j = 0; j < i; j++)
        {
            eta += DX_i[j] * DX_i[j];
            eta_d += D_i[j] * D_i[j];
            rho += D_i[j] * DX_i[j];
        }
    }
}

void NLDR::update_MAJORIZATION(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    switch (paras->nst)
    {
    case KRUSKAL1:
        update_KRUSKAL1_MAJORIZATION(kws, ind);
        break;
    case NORMALIZED:
        update_NORMALIZED_MAJORIZATION(kws, ind);
        break;
    case SAMMON:
        update_SAMMON_MAJORIZATION(kws, ind);
        break;
    case CCA:
        update_CCA_MAJORIZATION(kws, ind);
        break;
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }
};

void NLDR::update_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{

    // in-place update
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    PRECISION stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    update_coordinates_inplace(kws, ind);
    if (paras->nst != KRUSKAL1)
        return;
    else
    {
        Cor2Dis_with_paras(X, paras);
        paras->pre_stress = KRUSKAL1_stress_by_parameter(X, paras);
    }
};

void NLDR::direction_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    switch (paras->nst)
    {
    case KRUSKAL1:
        direction_KRUSKAL1_STOCHASTIC(kws, ind);
        break;
    case NORMALIZED:
        direction_NORMALIZED_STOCHASTIC(kws, ind);
        break;
    case SAMMON:
        direction_SAMMON_STOCHASTIC(kws, ind);
        break;
    case CCA:
        direction_CCA_STOCHASTIC(kws, ind);
        break;
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }
};

void NLDR::update_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // not-in-place update
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    NLDR_Var *Y = reinterpret_cast<NLDR_Var *>(kws.arg(ind[3]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    PRECISION T = paras->T;
    PRECISION r_scalar = 0;
    unsigned n = paras->Dis->dimension();
    unsigned d = paras->dim;
    r_scalar = sqrt(T * stepsize / n / d);

    auto eng = paras->eng;
    std::normal_distribution<PRECISION> rand(0.0, 1.0);
    // Matrix<PRECISION> Perturb(n, d);

    auto &Perturb = *paras->perturb;
    auto Perturb_ = Perturb.get_vec();
    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < d; k++)
            Perturb_[i][k] = rand(*eng);

    update_coordinates_notinplace(kws, ind);
    shift(*(Y->Cor), Perturb, r_scalar);
    switch (paras->nst)
    {
    case KRUSKAL1:
        Cor2Dis_with_paras(Y, paras);
        break;
    case NORMALIZED:
        Cor2Dis_with_paras(Y, paras);
        break;
    case SAMMON:
        Cor2Dis_no_paras(Y, paras);
        break;
    case CCA:
        Cor2Dis_no_paras(Y, paras);
        break;
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }
};

void NLDR::direction_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    bool accepted = *reinterpret_cast<bool *>(kws.arg(12));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    switch (paras->nst)
    {
    case KRUSKAL1:
        if (accepted)
            direction_KRUSKAL1_METROPOLIS(kws, ind);
        break;
    case NORMALIZED:
        if (accepted)
            direction_NORMALIZED_METROPOLIS(kws, ind);
        break;
    case SAMMON:
        if (accepted)
            direction_SAMMON_METROPOLIS(kws, ind);
        break;
    case CCA:
        direction_CCA_METROPOLIS(kws, ind); // Always recompute direction (with different lambda)
        break;
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }
};

bool NLDR::accept_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    PRECISION f_cur = *reinterpret_cast<PRECISION *>(kws.arg(3));
    PRECISION f_new = *reinterpret_cast<PRECISION *>(kws.arg(5));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    auto eng = paras->eng;
    std::uniform_real_distribution<PRECISION> rand(0.0, 1.0);

    bool accepted = false;
    switch (paras->nst)
    {
    case KRUSKAL1:
        accepted = (rand(*eng) < exp((f_cur - f_new) / paras->T));
        break;
    case NORMALIZED:
        accepted = (rand(*eng) < exp((f_cur - f_new) / paras->T));
        break;
    case SAMMON:
        accepted = (rand(*eng) < exp((f_cur - f_new) / paras->T));
        break;
    case CCA:
        accepted = (rand(*eng) < exp((f_cur - f_new) / (paras->T * paras->eta_d)));
        break;
    default:
        std::cout << "Cost function not found!\n";
        throw(1);
    }

    paras->T = paras->T * 0.99;
}

// void NLDR::direction_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind){
//     NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));

//     auto d = Dir->Delta->get_col();
//     if (Dir->activated_coo != d - 1)
//     {
//         Dir->activated_coo = Dir->activated_coo + 1;
//         return;
//     }

//     auto paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
//     auto n = paras->Dis->dimension();

//     Dir->activated_coo = 0;
//     Dir->activated_row = (Dir->activated_row + 1) % n;

//     //std::cout << Dir->activated_row << "\t" << Dir->activated_coo << "\n";

//     auto nst = paras->nst;
//     switch (nst)
//     {
//     case KRUSKAL1:
//         direction_KRUSKAL1_GAUSS_SEIDEL(kws, ind);
//         break;
//     case NORMALIZED:
//         direction_NORMALIZED_GAUSS_SEIDEL(kws, ind);
//         break;
//     case SAMMON:
//         direction_SAMMON_GAUSS_SEIDEL(kws, ind);
//         break;
//     case CCA:
//         direction_CCA_GAUSS_SEIDEL(kws, ind);
//         break;
//     default:
//         std::cout << "Cost function not found!\n";
//         throw(1);
//     }
// }

void NLDR::slope_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Dir *Delta = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    PRECISION &s = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));

    auto k = Delta->activated_coo;
    auto e1 = (*Delta->Delta)(0, k);
    auto e2 = (*Delta->Delta)(1, k);
    s = -e1 * e1 / e2;
}

// void NLDR::update_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
// {
//     // This update routine is called in Line_search wolfe.
//     // Perform not inplace update from VarDiff, no update on cost.
//     auto X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
//     auto Dir = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
//     auto stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
//     auto Y = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[3]));
//     auto Paras = reinterpret_cast<NLDR::NLDR_Paras_GAU *>(kws.carg(0));

//     pre_update_by_index_KRUSKAL1_GAUSS_SEIDEL(X, Dir, stepsize, Paras);

//     unsigned i = Paras->VarDiff->activated_row;
//     unsigned k = Paras->VarDiff->activated_coo;
//     unsigned n = X->Cor->get_row();
//     unsigned d = X->Cor->get_col();

//     Matrix<PRECISION> &CY_ = *Y->Cor;
//     LowerTri<PRECISION> &DY_ = *Y->Dis;
//     Matrix<PRECISION> &Delta_ = *Dir->Delta;

//     CY_(i, k) = (*X->Cor)(i, k) - stepsize * Delta_(0, k) / Delta_(1, k);
//     auto newD = Paras->VarDiff->newD_row;
//     for (auto j = 0; j < i; j++)
//         DY_(i, j) = newD[j];
//     DY_(i, i) = newD[i];
//     for (auto j = i + 1; j < n; j++)
//         DY_(j, i) = newD[j];

//     *Y->eta = Paras->VarDiff->eta;
//     *Y->eta_d = Paras->eta_d;
//     *Y->rho = Paras->VarDiff->rho;
// }

void NLDR::cost_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    auto *X = reinterpret_cast<NLDR_Var *>(kws.arg(1));
    auto *fx = reinterpret_cast<PRECISION *>(kws.arg(2));
    // This is dangerous. The updated cost need the old point for fast computation
    // but on the scope of OptimLib_Search_StepSize::Wolfe_1st_Condition::Cost,
    // only the current point is visible, referenced in ind[0] = 6. The old point,
    // on the global scope of OptimLib_Search_StepSize is referenced in 1.
    NLDR_Var *Y = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *f = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));
    NLDR::NLDR_Paras *Paras = reinterpret_cast<NLDR::NLDR_Paras *>(kws.carg(0));

    if (Paras->nst == KRUSKAL1 || Paras->nst == NORMALIZED)
        *f = stress_raw(Y, Paras); // Use new point
    else if (Paras->nst == SAMMON)
        SAMMON_stress_by_index(X, *f, Paras); // Use old point with update info stored in Paras->Var_Diff.
    else if (Paras->nst == CCA)
        CCA_stress_by_index(X, *f, Paras); // Use old point with update info stored in Paras->Var_Diff.
    else
    {
        std::cout << "Cost function type in NLDR_Paras not found!\n";
        throw(1);
    }
}

void NLDR::cost_GAUSS_SEIDEL_iteration(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // Every k-iterations(a row get updated), CCA adopt new lambda and need to recompute cost function,
    // Every nk-iterations(entire coordinates get updated), compute the error.
    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto fx = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));
    auto Y = reinterpret_cast<NLDR_Var *>(kws.arg(10));

    auto Paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));
    auto a_c = Paras->VarDiff->activated_coo;
    auto a_r = Paras->VarDiff->activated_row;
    if (a_c == (Paras->dim - 1) && Paras->nst == CCA)
        *fx = CCA_stress_by_distance_lambda(X->Dis, Paras->Dis, Paras->lambda[*Paras->iter / Paras->dim + 1]);

    if (a_r == (Paras->Dis->dimension() - 1) && a_c == (Paras->dim - 1))
    {
        centralize(*X->Cor);
        auto CX_ = X->Cor->get_vec();
        auto CY_ = Y->Cor->get_vec();
        for (auto i = 0; i < X->Cor->get_row(); i++)
            for (auto j = 0; j < X->Cor->get_col(); j++)
                CY_[i][j] = CX_[i][j];
        // Shift corrdinates so that they center at origin.
        auto err = reinterpret_cast<PRECISION *>(kws.arg(2));
        *err = (Paras->pre_stress - *fx) / *fx;
        Paras->pre_stress = *fx;
        std::cout << "Current stress:\t" << *fx << ", Relative error:\t" << *err << '\n';
    }
}

void NLDR::copy_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    auto src = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto des = reinterpret_cast<NLDR_Var *>(kws.arg(ind[1]));
    auto Paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));

    auto i = Paras->VarDiff->activated_row;
    auto k = Paras->VarDiff->activated_coo;
    auto &Dis_des = *des->Dis;
    auto &Dis_src = *src->Dis;
    (*des->Cor)(i, k) = (*src->Cor)(i, k);
    auto n = Paras->Dis->dimension();
    for (auto j = 0; j < i; j++)
        Dis_des(i, j) = Dis_src(i, j);
    for (auto j = i + 1; j < n; j++)
        Dis_des(j, i) = Dis_src(j, i);
    *(des->eta) = *(src->eta);
    *(des->rho) = *(src->rho);
}

void NLDR::MAJORIZATION(LowerTri<PRECISION> *Dis, unsigned n, unsigned d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
                        NLDR_STRESS_TYPE stress_type, unsigned MaxIter, unsigned MaxTime, PRECISION ErrTol, PRECISION DifTol,
                        PRECISION *lambda)
{
    using namespace OptimLib;

    std::cout << "\tInitial parameters of solver...\n";

    unsigned it = 0, ti = 0;
    PRECISION err = 100, diff = 100, nD = 0;

    Optim_Paras o_paras(MaxIter, MaxTime, ErrTol, DifTol);

    LowerTri<PRECISION> BX(n, true);
    Matrix<PRECISION> *V = nullptr;
    Matrix<PRECISION> *CX_temp = nullptr;
    Array<PRECISION> *A = nullptr;

    Array<PRECISION> row_sum(n);

    NLDR_Paras nldr_paras(Dis, &BX, nullptr, d, nullptr, lambda, 100);
    nldr_paras.set_stress_type(stress_type);
    nldr_paras.iter = &it;

    Matrix<PRECISION> CX(*Cor_0);
    LowerTri<PRECISION> DX(n, true);
    PRECISION etaX = 0, eta_d = 0, rhoX = 0, stressX;
    NLDR_Var X(&CX, &DX, &etaX, &eta_d, &rhoX);
    Init_X(&X, &stressX, &nldr_paras, stress_raw, Dis_0);
    nldr_paras.pre_stress = stressX;
    nldr_paras.row_sum = &row_sum;

    std::cout << "\tCompute V matrix...\n";

    if (stress_type == KRUSKAL1)
    {
        CX_temp = new Matrix<PRECISION>(n, d, 0.0);
        nldr_paras.CX_temp = CX_temp;
    }
    else if (stress_type == NORMALIZED)
    {
        CX_temp = new Matrix<PRECISION>(n, d, 0.0);
        nldr_paras.CX_temp = CX_temp;
    }
    else if (stress_type == SAMMON)
    {
        V = new Matrix<PRECISION>(n, n, 0.0);
        nldr_paras.V = V;
        A = new Array<PRECISION>(n * d);
        nldr_paras.A = A;
        compute_V_SAMMON(nullptr, &nldr_paras);
    }
    else if (stress_type == CCA)
    {
        V = new Matrix<PRECISION>(n, n, 0.0);
        nldr_paras.V = V;
        A = new Array<PRECISION>(n * d);
        nldr_paras.A = A;
    }

    std::cout << "\tInitial components in X...\n";

    // Get inital info.

    std::cout << "\tInitial stress:\t" << stressX << "\n";

    std::cout << "\tInitial Optim_Iter object...\n";

    Optim_KwArg kws_iter(OPTIM_ITERATION, 1);
    kws_iter.init_iter_routine(nullOptimFunc, nullOptimFunc, update_MAJORIZATION, nullOptimFunc, nullOptimFunc);
    kws_iter.init_iter_variable(&it, &ti, &err, &diff, (void *)&X, &stressX, nullptr, nullptr, &nD, nullptr);
    kws_iter.set_carg(&nldr_paras, 0);
    Optim_Iter o_iter(&kws_iter);

    std::cout << "\tInitial Optim_Solver object...\n";

    Optim_Update_Solver o_solver(&o_paras, &o_iter);

    std::cout << "\tStart Solving...\n";

    o_solver.solve();
    it -= 1;
    stressX = stress_raw(&X, &nldr_paras);
    std::cout << "Final stress:\t" << stressX << "\t computed time:\t" << (double)ti / CLOCKS_PER_SEC << ".\n";

    (*Cor_0) = Matrix<PRECISION>(*(X.Cor));
    // if (Dis_0)
    //     (*Dis_0) = LowerTri<PRECISION>(*(X.Dis));

    if (V)
        delete[] V;
    if (A)
        delete[] A;
    if (CX_temp)
        delete[] CX_temp;
}

void NLDR::STOCHASTIC(LowerTri<PRECISION> *Dis, unsigned n, unsigned d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
                      NLDR_STRESS_TYPE stress_type, unsigned MaxIter, unsigned MaxTime, PRECISION ErrTol, PRECISION DifTol,
                      PRECISION *alpha, PRECISION *lambda, unsigned *s_ind)
{
    using namespace OptimLib;

    unsigned it = 0, ti = 0;
    PRECISION err = 100, diff = 100, nD = 0;
    Optim_Paras o_paras(MaxIter, MaxTime, ErrTol, DifTol);

    NLDR_Paras nldr_paras(Dis, nullptr, nullptr, d, nullptr, lambda, 100);
    nldr_paras.set_stress_type(stress_type);
    nldr_paras.iter = &it;
    nldr_paras.s_ind = s_ind;

    Matrix<PRECISION> CX(*Cor_0);
    LowerTri<PRECISION> DX(n, true);
    PRECISION etaX = 0, eta_d = 0, rhoX = 0, stressX;
    NLDR_Var X(&CX, &DX, &etaX, &eta_d, &rhoX);
    Init_X(&X, &stressX, &nldr_paras, stress_raw, Dis_0);
    nldr_paras.pre_stress = stressX;

    // Get inital info.

    Matrix<PRECISION> DeltaX(n, d, 0.0);
    NLDR_Dir DirX(&DeltaX, -1, -1);

    PRECISION stepsize = 1.0;

    Optim_Computed_StepSize o_stepsize(alpha);

    Optim_KwArg kws_iter(OPTIM_ITERATION, 1);
    kws_iter.init_iter_routine(nullOptimFunc, direction_STOCHASTIC, update_STOCHASTIC, nullOptimFunc, nullOptimFunc);
    kws_iter.init_iter_variable(&it, &ti, &err, &diff, (void *)&X, &stressX, (void *)&DirX, nullptr, &nD, &stepsize);
    kws_iter.set_carg(&nldr_paras, 0);
    Optim_Iter o_iter(&kws_iter);

    Optim_Stepping_Solver o_solver(&o_paras, &o_iter, &o_stepsize);

    o_solver.solve();

    Cor2Dis_with_paras(&X, &nldr_paras);
    it -= 1;
    stressX = stress_raw(&X, &nldr_paras);
    std::cout << "Final stress:\t" << stressX << "\t computed time:\t" << (double)ti / CLOCKS_PER_SEC << ".\n";

    (*Cor_0) = Matrix<PRECISION>(*(X.Cor));
    // if (Dis_0)
    //     (*Dis_0) = LowerTri<PRECISION>(*(X.Dis));
}

void NLDR::METROPOLIS(LowerTri<PRECISION> *Dis, unsigned n, unsigned d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
                      NLDR_STRESS_TYPE stress_type, unsigned MaxIter, unsigned MaxTime, PRECISION ErrTol, PRECISION DifTol,
                      PRECISION *lambda, std::default_random_engine *eng, PRECISION T)
{
    using namespace OptimLib;

    unsigned it = 0, ti = 0;
    PRECISION err = 100, diff = 100, nD = 0;

    Optim_Paras o_paras(MaxIter, MaxTime, ErrTol, DifTol);

    NLDR_Paras nldr_paras(Dis, nullptr, nullptr, d, nullptr, lambda, 100);
    // No need to store V matrix, Var_Diff for update by index, lambda and weight.
    nldr_paras.set_stress_type(stress_type);
    nldr_paras.iter = &it;
    nldr_paras.T = T;
    nldr_paras.eng = eng;

    Matrix<PRECISION> CX(*Cor_0);
    LowerTri<PRECISION> DX(n, true);
    PRECISION etaX = 0, eta_d = 0, rhoX = 0, stressX;
    NLDR_Var X(&CX, &DX, &etaX, &eta_d, &rhoX);
    Init_X(&X, &stressX, &nldr_paras, stress_raw, Dis_0);

    // Get inital info.

    Matrix<PRECISION> CY(*Cor_0);
    LowerTri<PRECISION> DY(n, true);
    PRECISION etaY = 0, rhoY = 0, stressY;
    NLDR_Var Y(&CY, &DY, &etaY, &eta_d, &rhoY);
    // Initialize trial point Y.

    Matrix<PRECISION> DeltaX(n, d, 0.0);
    NLDR_Dir DirX(&DeltaX, -1, -1);

    Matrix<PRECISION> perturb(n, d);
    nldr_paras.perturb = &perturb;

    PRECISION stepsize;
    bool accepted = true;

    Optim_KwArg kws_iter(OPTIM_ITERATION, 1);
    kws_iter.init_iter_routine(stress, direction_METROPOLIS, update_METROPOLIS, nullOptimFunc, NLDR_Var_Copy);
    kws_iter.init_iter_variable(&it, &ti, &err, &diff, (void *)&X, &stressX, (void *)&DirX,
                                nullptr, &nD, &stepsize, (void *)&Y, &stressY, &accepted);
    kws_iter.set_carg(&nldr_paras, 0);
    Optim_Iter o_iter(&kws_iter);

    unsigned first_gradient[2] = {4, 6};

    direction_METROPOLIS(kws_iter, first_gradient); // Compute the first gradient.

    switch (stress_type)
    {
    case KRUSKAL1:
        stepsize = 1.0 / DeltaX.norm2(M_FOR_NORM);
        break;
    case NORMALIZED:
        stepsize = 1.0 / DeltaX.norm2(M_FOR_NORM);
        break;
    case SAMMON:
        stepsize = 1.0 / DeltaX.norm2(M_FOR_NORM);
        break;
    case CCA:
        stepsize = eta_d / DeltaX.norm2(M_FOR_NORM);
        break;
    default:
        std::cout << "Stress function not found!\n";
        throw(1);
    }
    // compute const stepsize

    Optim_Const_StepSize o_stepsize(stepsize);

    Optim_Trial_Stepping_Solver o_solver(&o_paras, &o_iter, &o_stepsize, accept_METROPOLIS, &accepted);

    o_solver.solve();

    Cor2Dis_with_paras(&X, &nldr_paras);
    it -= 1;
    stressX = stress_raw(&X, &nldr_paras);
    std::cout << "Final stress:\t" << stressX << "\t computed time:\t" << (double)ti / CLOCKS_PER_SEC << ".\n";

    (*Cor_0) = Matrix<PRECISION>(*(X.Cor));
    // if (Dis_0)
    //     (*Dis_0) = LowerTri<PRECISION>(*(X.Dis));
}

void NLDR::GAUSS_SEIDEL(LowerTri<PRECISION> *Dis, unsigned n, unsigned d, Matrix<PRECISION> *Cor_0, LowerTri<PRECISION> *Dis_0,
                        NLDR_STRESS_TYPE stress_type, unsigned MaxIter, unsigned MaxTime, PRECISION ErrTol, PRECISION DifTol,
                        PRECISION *lambda)
{
    using namespace OptimLib;

    unsigned it = 0, ti = 0, oit = 0;
    PRECISION err = 100, diff = 100, nD = 0;

    Optim_Paras o_paras(MaxIter, MaxTime, ErrTol, DifTol);

    NLDR_Paras_GAU nldr_paras(Dis, d);

    PRECISION *new_distance_row = new PRECISION[n];
    memset(new_distance_row, 0, n * sizeof(PRECISION));
    NLDR_Var_Diff NVD(0, 0, 0, 0, new_distance_row);
    nldr_paras.VarDiff = &NVD;

    Matrix<PRECISION> XD(n, d);
    Matrix<PRECISION> XD2(n, d);
    Array<PRECISION> DTX(n);
    Array<PRECISION> DTX3(n);
    Array<PRECISION> DR(n);
    nldr_paras.XD = &XD;
    nldr_paras.XD2 = &XD2;
    nldr_paras.DTX = &DTX;
    nldr_paras.DTX3 = &DTX3;
    nldr_paras.DR = &DR;

    Matrix<PRECISION> CX(*Cor_0);
    LowerTri<PRECISION> DX(n);
    DX.form_row_ptr();
    PRECISION etaX = 0, eta_d = nldr_paras.eta_d, rhoX = 0, stressX;
    NLDR_Var X(&CX, &DX, &etaX, &eta_d, &rhoX);
    Init_X(&X, Dis_0, Dis);

    Matrix<PRECISION> CY(*Cor_0);
    LowerTri<PRECISION> DY(DX);
    PRECISION etaY = 0, rhoY = 0, stressY;
    NLDR_Var Y(&CY, &DY, &etaY, &eta_d, &rhoY);
    DY.form_row_ptr();
    // Initialize trial point Y for wolfe_1st search.

    Matrix<PRECISION> DeltaX(2, d, 0.0);
    NLDR_Dir DirX(&DeltaX, -1, d - 1);

    PRECISION stepsize, slope;

    Optim_KwArg kws_iter(OPTIM_ITERATION, 1);
    kws_iter.init_iter_variable(&it, &ti, &err, &diff, (void *)&X, &stressX, (void *)&DirX,
                                &slope, &nD, &stepsize, (void *)&Y, &stressY, nullptr);
    kws_iter.set_carg(&nldr_paras, 0);

    Optim_KwArg kws_wolfe(WOLFE_1ST_STEPSIZE_SEARCH, 1);
    kws_wolfe.init_wolfe1st_variable(&err, (void *)&X, &stressX, (void *)&DirX, &slope, &stepsize, (void *)&Y, &stressY);
    kws_wolfe.set_carg(&nldr_paras, 0);

    if (stress_type == KRUSKAL1)
    {
        kws_iter.init_iter_routine(nullOptimFunc, nullOptimFunc, update_KRUSKAL1_GAUSS_SEIDEL, nullOptimFunc, nullOptimFunc);
        kws_wolfe.init_wolfe1st_routine(cost_KRUSKAL1_GAUSS_SEIDEL_line_search, update_KRUSKAL1_GAUSS_SEIDEL_line_search, copy_GAUSS_SEIDEL_line_search);
        stressX = KRUSKAL1_stress_by_parameter(&X);
        stressY = stressX;
    }
    else if (stress_type == NORMALIZED)
    {
        kws_iter.init_iter_routine(nullOptimFunc, nullOptimFunc, update_NORMALIZED_GAUSS_SEIDEL, nullOptimFunc, nullOptimFunc);
        kws_wolfe.init_wolfe1st_routine(cost_NORMALIZED_GAUSS_SEIDEL_line_search, update_NORMALIZED_GAUSS_SEIDEL_line_search, copy_GAUSS_SEIDEL_line_search);
        stressX = NORMALIZED_stress_by_parameter(&X);
        stressY = stressX;
    }
    else if (stress_type == SAMMON)
    {
        kws_iter.init_iter_routine(nullOptimFunc, nullOptimFunc, update_SAMMON_GAUSS_SEIDEL, nullOptimFunc, nullOptimFunc);
        kws_wolfe.init_wolfe1st_routine(cost_SAMMON_GAUSS_SEIDEL_line_search, update_SAMMON_GAUSS_SEIDEL_line_search, copy_GAUSS_SEIDEL_line_search);
        stressX = SAMMON_stress_by_distance(&X, Dis, nldr_paras.total_d);
        stressY = stressX;
    }
    else if (stress_type == CCA)
    {
        kws_iter.init_iter_routine(nullOptimFunc, nullOptimFunc, update_CCA_GAUSS_SEIDEL, nullOptimFunc, nullOptimFunc);
        kws_wolfe.init_wolfe1st_routine(cost_CCA_GAUSS_SEIDEL_line_search, update_CCA_GAUSS_SEIDEL_line_search, copy_GAUSS_SEIDEL_line_search);
        nldr_paras.lambda = lambda;
        nldr_paras.cur_lambda = lambda;
        stressX = CCA_stress_by_distance_lambda(&DX, Dis, *lambda);
        stressY = stressX;
    }
    else
    {
        std::cout << "Stress type not found!\n";
        throw(1);
    }

    Optim_Iter o_iter(&kws_iter);

    Wolfe_1st_StepSize o_stepsize(&kws_wolfe, 1e-5, 1.0, 0.01);
    nldr_paras.wolfe = &o_stepsize;

    Optim_Update_Solver o_solver(&o_paras, &o_iter);

    // Optim_Stepping_Solver o_solver(&o_paras, &o_iter, &o_stepsize);
    // Optim_Stepping_Solver o_solver(&o_paras, &o_iter, &o_cstepsize);

    o_solver.solve();

    std::cout << "Final stress:\t" << stressX << "\t computed time:\t" << (double)ti / CLOCKS_PER_SEC << ".\n";

    (*Cor_0) = Matrix<PRECISION>(*(X.Cor));
    // if (Dis_0)
    //     (*Dis_0) = LowerTri<PRECISION>(*(X.Dis));

    delete[] new_distance_row;
}

/****************************************KRUSKAL1****************************************/

void NLDR::KRUSKAL1_distance_stress_function(NLDR::NLDR_Var *X, PRECISION *f, NLDR::NLDR_Paras *Paras)
{
    PRECISION &stress = *f;
    Cor2Dis_with_paras(X, Paras);
    stress = KRUSKAL1_stress_by_parameter(X, Paras);
}

PRECISION NLDR::KRUSKAL1_stress_by_parameter(NLDR_Var *X, NLDR_Paras *Paras)
{
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);
    return (eta_d + eta - 2 * rho) / eta;
}

PRECISION NLDR::KRUSKAL1_stress_by_parameter(NLDR_Var *X)
{
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);
    return (eta_d + eta - 2 * rho) / eta;
}

void NLDR::KRUSKAL1_rescale(NLDR::NLDR_Var *X)
{
    Matrix<PRECISION> &C = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    PRECISION &eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION &rho = *(X->rho);
    unsigned dim = C.get_col();
    unsigned n = C.get_row();

    PRECISION beta = eta_d / rho;
    for (auto i = 0; i < n; i++)
    {
        PRECISION *DX_i = DX[i];
        for (auto j = 0; j < i; j++)
            DX_i[j] *= beta;
        for (auto k = 0; k < dim; k++)
            C(i, k) *= beta;
    }
    eta *= (beta * beta);
    rho *= beta;
}

void NLDR::compute_BX_KRUSKAL1(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &BX = *(Paras->BX);
    PRECISION *row_sum = Paras->row_sum->get_vec();
    memset(row_sum, 0, n * sizeof(PRECISION));

    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        auto BX_i = BX[i];
        for (auto j = 0; j < i; j++)
            if (DX_i[j] != 0)
            {
                BX_i[j] = -D_i[j] / DX_i[j];
                row_sum[i] += BX_i[j];
                row_sum[j] += BX_i[j];
            }
            else
                BX_i[j] = 0;
    }

    for (auto i = 0; i < n; i++)
        BX(i, i) = -row_sum[i];
}

void NLDR::update_cor_by_BX_KRUSKAL1(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    auto &Cor = *(X->Cor);
    auto &BX = *(Paras->BX);
    auto XC_ = Cor.get_vec();

    auto &CX_temp = *(Paras->CX_temp);
    auto CX_temp_ = CX_temp.get_vec();

    for (auto i = 0; i < n; i++)
        memcpy(CX_temp_[i], XC_[i], dim * sizeof(PRECISION));

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            XC_[i][k] = 0;
            for (auto j = 0; j <= i; j++)
                XC_[i][k] += BX(i, j) * CX_temp_[j][k];
            for (auto j = i + 1; j < n; j++)
                XC_[i][k] += BX(j, i) * CX_temp_[j][k];
            // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
            XC_[i][k] /= n;
        }
}

void NLDR::update_KRUSKAL1_MAJORIZATION(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // ip_Update(X, ..., paras)
    // Taking X, stress, error | paras
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *stress = reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    PRECISION *error = reinterpret_cast<PRECISION *>(kws.arg(ind[4]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    compute_BX_KRUSKAL1(X, paras);
    update_cor_by_BX_KRUSKAL1(X, paras);
    Cor2Dis_with_paras(X, paras);
    KRUSKAL1_rescale(X);
    *stress = KRUSKAL1_stress_by_parameter(X, paras);

    *error = (paras->pre_stress - *stress) / (*stress);
    paras->pre_stress = *(stress);
}

void NLDR::direction_KRUSKAL1_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned i = paras->s_ind[*(paras->iter)];
    LowerTri<PRECISION> &D = *(paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    Matrix<PRECISION> &CX = *(X->Cor);
    Matrix<PRECISION> &Delta = *(Dir->Delta);
    unsigned n = D.dimension();
    unsigned dim = CX.get_col();

    PRECISION D_ij, DX_ij;

    PRECISION **CX_ = CX.get_vec();
    PRECISION **Delta_ = Delta.get_vec();

    PRECISION a = 0, b = 0, c = 0, temp = 0, eta = *(X->eta);
    // time_t start = clock();
    //  for (auto D_i = D.row_begin(a_row), DX_i = DX.row_begin(a_row), end = D.row_end(a_row) ; D_i != end; D_i++, DX_i++)
    //  {
    //      //temp = (D(a_row, j) - DX(a_row, j));
    //      temp = (*D_i) - (*DX_i);
    //      b += temp * temp;
    //  }

    for (int j = 0; j < i; j++)
    {
        temp = D(i, j) - DX(i, j);
        b += temp * temp;
    }

    for (int j = i + 1; j < n; j++)
    {
        temp = D(j, i) - DX(j, i);
        b += temp * temp;
    }

    // for (auto j = 0; j < n; j++)
    // for (auto D_i = D.row_begin(a_row), DX_i = DX.row_begin(a_row), end = D.row_end(a_row);
    //      D_i != end;
    //      D_i++, DX_i++)
    //     for (auto k = 0; k < dim; k++)
    //     {
    //         unsigned j = D_i.get_col();
    //         D_ij = *D_i;
    //         DX_ij = *DX_i;
    //         if (DX_ij > 1e-7)
    //             a = -2.0 * (D_ij - DX_ij) / DX_ij * (CX_[j][k] - CX_[a_row][k]) * eta;
    //         else
    //             a = 0;
    //         c = 0;
    //         for (auto i = 0; i < n; i++)
    //             c += (CX_[j][k] - CX_[a_row][k]);
    //         Delta_[j][k] = -(a - b * c) / eta / eta;
    //     }

    for (int j = 0; j < i; j++)
        for (auto k = 0; k < dim; k++)
        {
            D_ij = D(i, j);
            DX_ij = DX(i, j);
            if (DX_ij > 1e-7)
                a = -2.0 * (D_ij - DX_ij) / DX_ij * (CX_[j][k] - CX_[i][k]) * eta;
            else
                a = 0;
            c = 0;
            for (auto i = 0; i < n; i++)
                c += (CX_[j][k] - CX_[i][k]);
            Delta_[j][k] = -(a - b * c) / eta / eta;
        }

    for (auto k = 0; k < dim; k++)
        Delta_[i][k] = 0;

    for (int j = i + 1; j < n; j++)
        for (auto k = 0; k < dim; k++)
        {
            D_ij = D(j, i);
            DX_ij = DX(j, i);
            if (DX_ij > 1e-7)
                a = -2.0 * (D_ij - DX_ij) / DX_ij * (CX_[j][k] - CX_[i][k]) * eta;
            else
                a = 0;
            c = 0;
            for (auto i = 0; i < n; i++)
                c += (CX_[j][k] - CX_[i][k]);
            Delta_[j][k] = -(a - b * c) / eta / eta;
        }
    // std::cout << "\t Computation time: " << (clock() - start) / (double) CLOCKS_PER_SEC << "\n";
}

void NLDR::direction_KRUSKAL1_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    PRECISION n1, n2;
    double eta = *(X->eta);
    double eta_d = *(X->eta_d);
    double rho = *(X->rho);
    Matrix<PRECISION> &Cor = *(X->Cor);
    Matrix<PRECISION> &DeltaX = *(Dir->Delta);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    auto CX_ = Cor.get_vec();
    auto DeltaX_ = DeltaX.get_vec();
    PRECISION D_ij, DX_ij;

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            n1 = 0;
            n2 = 0;
            for (auto j = 0; j < i; j++)
            {
                n1 += CX_[i][k] - CX_[j][k];
                D_ij = D(i, j);
                DX_ij = DX(i, j);
                if (DX_ij > 1e-8)
                    n2 += D_ij / DX_ij * (CX_[i][k] - CX_[j][k]);
            }

            for (auto j = i + 1; j < n; j++)
            {
                n1 += CX_[i][k] - CX_[j][k];
                D_ij = D(j, i);
                DX_ij = DX(j, i);
                if (DX_ij > 1e-8)
                    n2 += D_ij / DX_ij * (CX_[i][k] - CX_[j][k]);
            }
            DeltaX_[i][k] = 2.0 / eta / eta * (n1 * eta + eta_d * n2 - 2 * rho * n2);
        }
}

// KRUSKAL1_GAUSS_SEIDEL

void NLDR::pre_update_by_index_KRUSKAL1_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
{
    Matrix<PRECISION> &C = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    Matrix<PRECISION> &Delta = *(Dir->Delta);

    unsigned dim = C.get_col(), n = C.get_row();
    unsigned i = Dir->activated_row;
    unsigned k = Dir->activated_coo;
    PRECISION e1 = Delta(0, k);
    PRECISION e2 = Delta(1, k);
    PRECISION xcnew = C(i, k) - stepsize * e1 / e2;

    PRECISION *newD_i = Paras->VarDiff->newD_row;
    Paras->VarDiff->activated_row = i;
    Paras->VarDiff->activated_coo = k;
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);

    PRECISION DX_ij, D_ij;
    auto C_ = C.get_vec();

    PRECISION temp = 0;
    for (auto j = 0; j < i; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        eta -= DX_ij * DX_ij;
        rho -= DX_ij * D_ij;
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
        eta += newD_i[j] * newD_i[j];
        rho += newD_i[j] * D_ij;
    }
    newD_i[i] = 0;
    for (auto j = i + 1; j < n; j++)
    {
        DX_ij = DX(j, i);
        D_ij = D(j, i);
        eta -= DX_ij * DX_ij;
        rho -= DX_ij * D_ij;
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
        eta += newD_i[j] * newD_i[j];
        rho += newD_i[j] * D_ij;
    }
    Paras->VarDiff->eta = eta;
    Paras->VarDiff->rho = rho;
}

void NLDR::update_KRUSKAL1_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This update routine is called in Line_search wolfe.
    // Perform not inplace update from VarDiff, no update on cost.
    auto X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
    auto stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    auto Y = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[3]));
    auto Paras = reinterpret_cast<NLDR::NLDR_Paras_GAU *>(kws.carg(0));

    pre_update_by_index_KRUSKAL1_GAUSS_SEIDEL(X, Dir, stepsize, Paras);

    update_distance_row_from_VarDiff(X, Y, Dir, stepsize, Paras);
}

void NLDR::cost_KRUSKAL1_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *Y = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *f = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));

    *f = KRUSKAL1_stress_by_parameter(Y);
}

void NLDR::direction_KRUSKAL1_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, unsigned a_row,
                                           PRECISION **XD_, PRECISION **XD2_, PRECISION *DTX_, PRECISION *DTX3_)
{

    auto n = paras->pts;
    auto d = paras->dim;

    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    Dir->activated_coo = 0;
    Dir->activated_row = a_row;
    auto i = a_row;

    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, n2, n3, n4, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho);

    /**************************Form XD, XD2**************************/
    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }
    /**************************Form XD, XD2**************************/

    /**************************Form DTX, DTX3**************************/
    for (auto j = 0; j < i; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        if (DX_ij > 1e-8)
        {
            DTX_[j] = D_ij / DX_ij;
            DTX3_[j] = D_ij / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }

    DTX_[i] = 0;
    DTX3_[i] = 0;

    for (auto j = i + 1; j < n; j++)
    {
        DX_ij = DX(j, i);
        D_ij = D(j, i);
        if (DX_ij > 1e-8)
        {
            DTX_[j] = D_ij / DX_ij;
            DTX3_[j] = D_ij / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }
    /**************************Form DTX, DTX32**************************/

    for (auto k = 0; k < d; k++)
    {
        a = 0;
        for (auto j = 0; j < n; j++)
            a += DTX_[j] * XD_[j][k];
        a *= eta;
        b = 0;
        for (auto j = 0; j < n; j++)
            b += XD_[j][k];
        n2 = b * eta_d;
        n3 = b * rho;
        Delta_[0][k] = -2 * (a + n2 - 2 * n3) / eta / eta;
        n4 = 0;
        for (auto j = 0; j < n; j++)
            n4 += -DTX3_[j] * XD2_[j][k] + DTX_[j];
        n4 *= eta;
        Delta_[1][k] = -2 * ((n4 + (n - 1) * eta_d - 2 * (n - 1) * rho) * eta * eta - (a + n2 - 2 * n3) * 4 * eta * b) / eta / eta / eta / eta;
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

void NLDR::update_KRUSKAL1_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This is called by the Optim_Update_Solver.Optim::Iter, ind = x, Dir, stepsize, stress, error

    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    auto &f = *reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    auto &e = *reinterpret_cast<PRECISION *>(kws.arg(ind[4]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));
    auto &wolfe = *paras->wolfe;
    auto &wolfe_kws = *wolfe.get_status();
    auto &slope = *reinterpret_cast<PRECISION *>(wolfe_kws.arg(4));

    auto f_pre = f;
    auto n = X->Cor->get_row(), d = X->Cor->get_col();

    auto XD_ = paras->XD->get_vec(), XD2_ = paras->XD2->get_vec();
    auto DTX_ = paras->DTX->get_vec(), DTX3_ = paras->DTX3->get_vec();

    auto Delta_ = Dir->Delta->get_vec();

    for (auto i = 0; i < n; i++)
    {
        direction_KRUSKAL1_GAUSS_SEIDEL(X, Dir, paras, i, XD_, XD2_, DTX_, DTX3_);
        for (auto k = 0; k < d; k++)
        {
            Dir->activated_coo = k;
            slope = -Delta_[0][k] * Delta_[0][k] / Delta_[1][k];
            wolfe.stepsize();
        }
    }
    centralize(*X->Cor);
    e = (f_pre - f) / f;
    std::cout << "stress value : " << f << '\n';
    std::cout << "rel error : " << e << '\n'
              << '\n';
}

/****************************************KRUSKAL1****************************************/

/****************************************NORMALIZED**************************************/
void NLDR::NORMALIZED_distance_stress_function(NLDR_Var *X, PRECISION *f, NLDR_Paras *Paras)
{
    PRECISION &stress = *f;
    Cor2Dis_with_paras(X, Paras);
    stress = NORMALIZED_stress_by_parameter(X, Paras);
}

PRECISION NLDR::NORMALIZED_stress_by_parameter(NLDR_Var *X, NLDR_Paras *paras)
{
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);
    return (eta_d + eta - 2 * rho) / eta_d;
}

PRECISION NLDR::NORMALIZED_stress_by_parameter(NLDR_Var *X)
{
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);
    return (eta_d + eta - 2 * rho) / eta_d;
}

void NLDR::compute_BX_NORMALIZED(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &BX = *(Paras->BX);
    PRECISION *row_sum = Paras->row_sum->get_vec();
    memset(row_sum, 0, n * sizeof(PRECISION));
    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        auto BX_i = BX[i];
        for (auto j = 0; j < i; j++)
            if (DX_i[j] != 0)
            {
                BX_i[j] = -D_i[j] / DX_i[j];

                row_sum[i] += BX_i[j];
                row_sum[j] += BX_i[j];
            }
            else
                BX_i[j] = 0;
    }

    for (auto i = 0; i < n; i++)
        BX(i, i) = -row_sum[i];
}

void NLDR::update_cor_by_BX_NORMALIZED(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    auto &Cor = *(X->Cor);
    auto &BX = *(Paras->BX);
    auto XC_ = Cor.get_vec();

    auto &CX_temp = *(Paras->CX_temp);
    auto CX_temp_ = CX_temp.get_vec();

    for (auto i = 0; i < n; i++)
        memcpy(CX_temp_[i], XC_[i], dim * sizeof(PRECISION));

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            XC_[i][k] = 0;
            for (auto j = 0; j <= i; j++)
                XC_[i][k] += BX(i, j) * CX_temp_[j][k];
            for (auto j = i + 1; j < n; j++)
                XC_[i][k] += BX(j, i) * CX_temp_[j][k];
            // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
            XC_[i][k] /= n;
        }
}

void NLDR::update_NORMALIZED_MAJORIZATION(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // ip_Update(X, ..., paras)
    // Taking X, stress, error | paras
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *stress = reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    PRECISION *error = reinterpret_cast<PRECISION *>(kws.arg(ind[4]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    compute_BX_NORMALIZED(X, paras);
    update_cor_by_BX_NORMALIZED(X, paras);
    Cor2Dis_with_paras(X, paras);
    KRUSKAL1_rescale(X);
    *stress = NORMALIZED_stress_by_parameter(X);

    *error = (paras->pre_stress - *stress) / (*stress);
    paras->pre_stress = *(stress);
}

void NLDR::direction_NORMALIZED_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned i = paras->s_ind[*(paras->iter)];
    LowerTri<PRECISION> &D = *(paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    Matrix<PRECISION> &CX = *(X->Cor);
    Matrix<PRECISION> &Delta = *(Dir->Delta);
    unsigned n = D.dimension();
    unsigned dim = CX.get_col();

    auto CX_ = CX.get_vec();
    auto Delta_ = Delta.get_vec();

    PRECISION scalar = 0;

    for (auto j = 0; j < n; j++)
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] = CX_[j][k] - CX_[i][k];

    for (auto j = 0; j < i; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-10)
            scalar = 1.0;
        scalar = (D(i, j) / scalar) - 1.0;
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] *= 2 * scalar;
    }

    scalar = -1.0;
    for (auto k = 0; k < dim; k++)
        Delta_[i][k] *= 2 * scalar;

    for (auto j = i + 1; j < n; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-10)
            scalar = 1.0;
        scalar = (D(j, i) / scalar) - 1.0;
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] *= 2 * scalar;
    }
}

void NLDR::direction_NORMALIZED_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    Matrix<PRECISION> &Cor = *(X->Cor);
    Matrix<PRECISION> &DeltaX = *(Dir->Delta);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    auto CX_ = Cor.get_vec();
    auto DeltaX_ = DeltaX.get_vec();
    PRECISION D_ij, DX_ij;

    PRECISION c = paras->eta_d;

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            DeltaX_[i][k] = 0;
            for (auto j = 0; j < i; j++)
            {
                D_ij = D(i, j);
                DX_ij = DX(i, j);
                if (DX_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij * (CX_[i][k] - CX_[j][k]);
            }
            for (auto j = i + 1; j < n; j++)
            {
                D_ij = D(j, i);
                DX_ij = DX(j, i);
                if (DX_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij * (CX_[i][k] - CX_[j][k]);
            }
            DeltaX_[i][k] /= c;
        }
}

void NLDR::direction_NORMALIZED_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    auto n = X->Cor->get_row();
    auto d = X->Cor->get_col();
    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    auto i = Dir->activated_row;

    Matrix<PRECISION> XD(n, d);
    Matrix<PRECISION> XD2(n, d);
    Array<PRECISION> DTX(n);
    Array<PRECISION> DTX3(n);

    auto XD_ = XD.get_vec();
    auto XD2_ = XD2.get_vec();
    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho);

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < n; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        if (DX_ij > 1e-8)
        {
            DTX[j] = D_ij / DX_ij;
            DTX3[j] = D_ij / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX[j] = 0;
            DTX3[j] = 0;
        }
    }

    for (auto k = 0; k < d; k++)
    {
        a = 0;
        for (auto j = 0; j < n; j++)
            a += DTX[j] - 1;
        b = 0;
        for (auto j = 0; j < n; j++)
            b += (DTX[j] - 1) * XD_[j][k];
        Delta_[0][k] = -2 * b / eta_d;
        c = 0;
        for (auto j = 0; j < n; j++)
            c += DTX3[j] * XD2_[j][k];
        Delta_[1][k] = -2 * (a - c) / eta_d;
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

// NORMALIZED_GAUSS_SEIDEL

void NLDR::pre_update_by_index_NORMALIZED_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
{
    Matrix<PRECISION> &C = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    Matrix<PRECISION> &Delta = *(Dir->Delta);

    unsigned dim = C.get_col(), n = C.get_row();
    unsigned i = Dir->activated_row;
    unsigned k = Dir->activated_coo;
    PRECISION e1 = Delta(0, k);
    PRECISION e2 = Delta(1, k);
    PRECISION xcnew = C(i, k) - stepsize * e1 / e2;

    PRECISION *newD_i = Paras->VarDiff->newD_row;
    Paras->VarDiff->activated_row = i;
    Paras->VarDiff->activated_coo = k;
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);

    PRECISION DX_ij, D_ij;
    auto C_ = C.get_vec();

    PRECISION temp = 0;
    for (auto j = 0; j < i; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        eta -= DX_ij * DX_ij;
        rho -= DX_ij * D_ij;
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
        eta += newD_i[j] * newD_i[j];
        rho += newD_i[j] * D_ij;
    }
    newD_i[i] = 0;
    for (auto j = i + 1; j < n; j++)
    {
        DX_ij = DX(j, i);
        D_ij = D(j, i);
        eta -= DX_ij * DX_ij;
        rho -= DX_ij * D_ij;
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
        eta += newD_i[j] * newD_i[j];
        rho += newD_i[j] * D_ij;
    }
    Paras->VarDiff->eta = eta;
    Paras->VarDiff->rho = rho;
}

void NLDR::update_NORMALIZED_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This update routine is called in Line_search wolfe.
    // Perform not inplace update from VarDiff, no update on cost.
    auto X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
    auto stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    auto Y = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[3]));
    auto Paras = reinterpret_cast<NLDR::NLDR_Paras_GAU *>(kws.carg(0));

    pre_update_by_index_NORMALIZED_GAUSS_SEIDEL(X, Dir, stepsize, Paras);

    update_distance_row_from_VarDiff(X, Y, Dir, stepsize, Paras);
}

void NLDR::cost_NORMALIZED_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *Y = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *f = reinterpret_cast<PRECISION *>(kws.arg(ind[1]));

    *f = NORMALIZED_stress_by_parameter(Y);
}

void NLDR::direction_NORMALIZED_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, unsigned a_row,
                                             PRECISION **XD_, PRECISION **XD2_, PRECISION *DTX_, PRECISION *DTX3_)
{

    auto n = paras->pts;
    auto d = paras->dim;

    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    Dir->activated_coo = 0;
    Dir->activated_row = a_row;
    auto i = a_row;

    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho);

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < i; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        if (DX_ij > 1e-8)
        {
            DTX_[j] = D_ij / DX_ij;
            DTX3_[j] = D_ij / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }
    DTX_[i] = 0;
    DTX3_[i] = 0;

    for (auto j = i + 1; j < n; j++)
    {
        DX_ij = DX(j, i);
        D_ij = D(j, i);
        if (DX_ij > 1e-8)
        {
            DTX_[j] = D_ij / DX_ij;
            DTX3_[j] = D_ij / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }

    for (auto k = 0; k < d; k++)
    {
        a = 0;
        for (auto j = 0; j < n; j++)
            a += DTX_[j] - 1;
        b = 0;
        for (auto j = 0; j < n; j++)
            b += (DTX_[j] - 1) * XD_[j][k];
        Delta_[0][k] = -2 * b / eta_d;
        c = 0;
        for (auto j = 0; j < n; j++)
            c += DTX3_[j] * XD2_[j][k];
        Delta_[1][k] = -2 * (a - c) / eta_d;
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

void NLDR::update_NORMALIZED_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This is called by the Optim_Update_Solver.Optim::Iter, ind = x, Dir, stepsize, stress, error

    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    auto &f = *reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    auto &e = *reinterpret_cast<PRECISION *>(kws.arg(ind[4]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));
    auto &wolfe = *paras->wolfe;
    auto &wolfe_kws = *wolfe.get_status();
    auto &slope = *reinterpret_cast<PRECISION *>(wolfe_kws.arg(4));

    auto f_pre = f;
    auto n = X->Cor->get_row(), d = X->Cor->get_col();

    auto XD_ = paras->XD->get_vec(), XD2_ = paras->XD2->get_vec();
    auto DTX_ = paras->DTX->get_vec(), DTX3_ = paras->DTX3->get_vec();

    auto Delta_ = Dir->Delta->get_vec();

    for (auto i = 0; i < n; i++)
    {
        direction_NORMALIZED_GAUSS_SEIDEL(X, Dir, paras, i, XD_, XD2_, DTX_, DTX3_);
        for (auto k = 0; k < d; k++)
        {
            Dir->activated_coo = k;
            slope = -Delta_[0][k] * Delta_[0][k] / Delta_[1][k];
            wolfe.stepsize();
        }
    }
    centralize(*X->Cor);
    e = (f_pre - f) / f;
    std::cout << "stress value : " << f << '\n';
    std::cout << "rel error : " << e << '\n'
              << '\n';
}

/****************************************NORMALIZED**************************************/

/****************************************SAMMON******************************************/

PRECISION NLDR::SAMMON_stress_by_distance(NLDR_Var *X, NLDR_Paras *Paras)
{
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);

    unsigned dim = Paras->dim;
    unsigned n = DX.dimension();

    PRECISION temp = 0;
    PRECISION stress = 0;

    for (auto i = 0; i < n; i++)
    {
        PRECISION *DX_i = DX.ele(i);
        PRECISION *D_i = D[i];
        for (auto j = 0; j < i; j++)
        {
            temp = DX_i[j] - D_i[j];
            if (D_i[j] > 1e-8)
                stress += (temp * temp) / D_i[j];
        }
    }

    // std::cout << "Total_d:\t" << Paras->total_d << "\n";
    stress /= Paras->total_d;

    return stress;
}

PRECISION NLDR::SAMMON_stress_by_distance(NLDR_Var *X, LowerTri<PRECISION> *D, PRECISION total_d)
{
    LowerTri<PRECISION> &DX = *(X->Dis);

    unsigned n = DX.dimension();

    PRECISION temp = 0;
    PRECISION stress = 0;

    for (auto i = 0; i < n; i++)
    {
        PRECISION *DX_i = DX.ele(i);
        PRECISION *D_i = D->ele(i);
        for (auto j = 0; j < i; j++)
        {
            temp = DX_i[j] - D_i[j];
            if (D_i[j] > 1e-8)
                stress += (temp * temp) / D_i[j];
        }
    }

    stress /= total_d;

    return stress;
}

void NLDR::SAMMON_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras *Paras)
{
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    PRECISION *newD_i = Paras->VarDiff->newD_row;
    PRECISION total_d = Paras->total_d;
    unsigned i = Paras->VarDiff->activated_row;
    unsigned n = X->Cor->get_row();

    PRECISION temp = 0;
    for (auto j = 0; j < i; j++)
    {
        auto D_ij = D(i, j);
        auto DX_ij = DX(i, j);
        if (D_ij >= 1e-10)
            temp += (2 * D_ij - newD_i[j] - DX_ij) * (DX_ij - newD_i[j]) / D_ij;
    }
    for (auto j = i + 1; j < n; j++)
    {
        auto D_ij = D(j, i);
        auto DX_ij = DX(j, i);
        if (D_ij >= 1e-10)
            temp += (2 * D_ij - newD_i[j] - DX_ij) * (DX_ij - newD_i[j]) / D_ij;
    }
    stress += temp / total_d;
}

void NLDR::SAMMON_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras_GAU *Paras)
{
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    PRECISION *newD_i = Paras->VarDiff->newD_row;
    PRECISION total_d = Paras->total_d;
    unsigned i = Paras->VarDiff->activated_row;
    unsigned n = X->Cor->get_row();

    PRECISION temp = 0;
    for (auto j = 0; j < i; j++)
    {
        auto D_ij = D(i, j);
        auto DX_ij = DX(i, j);
        if (D_ij >= 1e-10)
            temp += (2 * D_ij - newD_i[j] - DX_ij) * (DX_ij - newD_i[j]) / D_ij;
    }
    for (auto j = i + 1; j < n; j++)
    {
        auto D_ij = D(j, i);
        auto DX_ij = DX(j, i);
        if (D_ij >= 1e-10)
            temp += (2 * D_ij - newD_i[j] - DX_ij) * (DX_ij - newD_i[j]) / D_ij;
    }
    stress += temp / total_d;
}

void NLDR::compute_BX_SAMMON(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &BX = *(Paras->BX);
    PRECISION *row_sum = Paras->row_sum->get_vec();
    memset(row_sum, 0, n * sizeof(PRECISION));
    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        auto BX_i = BX[i];
        for (auto j = 0; j < i; j++)
            if (DX_i[j] != 0)
            {
                BX_i[j] = -1.0 / DX_i[j];
                row_sum[i] += BX_i[j];
                row_sum[j] += BX_i[j];
            }
            else
                BX_i[j] = 0;
    }

    for (auto i = 0; i < n; i++)
        BX(i, i) = -row_sum[i];
}

void NLDR::compute_V_SAMMON(NLDR_Var *X, NLDR_Paras *Paras)
{
    auto D = Paras->Dis;
    auto V = Paras->V;
    auto V_ = V->get_vec();
    auto n = D->dimension();

    for (auto i = 0; i < n; i++)
    {
        auto D_i = (*D)[i];
        for (auto j = 0; j < i; j++)
        {
            if (D_i[j] < 1e-5)
                V_[i][j] = -1e5;
            else
                V_[i][j] = -1.0 / D_i[j];
            V_[j][i] = V_[i][j];
            V_[i][i] -= V_[i][j];
            V_[j][j] -= V_[j][i];
        }
    }
    shift(*V, 1.0);
    inverse(*V);
    shift(*V, -1.0 / n / n);
};

void NLDR::update_cor_by_BX_SAMMON(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    auto &Cor = *(X->Cor);
    auto &BX = *(Paras->BX);

    auto XC_ = Cor.get_vec();
    auto &V = *(Paras->V);
    auto V_ = V.get_vec();

    auto A = Paras->A->get_vec();
    memset(A, 0, n * dim * sizeof(PRECISION));
    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            for (auto j = 0; j <= i; j++)
                A[i * dim + k] += BX(i, j) * XC_[j][k];
            for (auto j = i + 1; j < n; j++)
                A[i * dim + k] += BX(j, i) * XC_[j][k];
            // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
        }

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            XC_[i][k] = 0;
            for (auto j = 0; j < n; j++)
                XC_[i][k] += V_[i][j] * A[j * dim + k];
        }
}

void NLDR::update_SAMMON_MAJORIZATION(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // ip_Update(X, ..., paras)
    // Taking X, stress, error ; paras

    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *stress = reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    PRECISION *error = reinterpret_cast<PRECISION *>(kws.arg(ind[4]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    compute_BX_SAMMON(X, paras);
    update_cor_by_BX_SAMMON(X, paras);
    Cor2Dis_no_paras(X, paras);
    *stress = SAMMON_stress_by_distance(X, paras);

    *error = (paras->pre_stress - *stress) / (*stress);
    paras->pre_stress = *(stress);
}

void NLDR::direction_SAMMON_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned i = paras->s_ind[*(paras->iter)];
    LowerTri<PRECISION> &D = *(paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    Matrix<PRECISION> &CX = *(X->Cor);
    Matrix<PRECISION> &Delta = *(Dir->Delta);
    unsigned n = D.dimension();
    unsigned dim = CX.get_col();

    auto CX_ = CX.get_vec();
    auto Delta_ = Delta.get_vec();

    PRECISION scalar = 0;

    for (auto j = 0; j < i; j++)
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] = CX_[j][k] - CX_[i][k];

    for (auto j = 0; j < i; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-8)
            scalar = 1.0;
        if (D(i, j) > 1e-8)
            scalar = (1.0 / scalar) - (1.0 / D(i, j));
        else
            scalar = (1.0 / scalar) - 1.0;
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] *= 2 * scalar;
    }

    scalar = 0;
    for (auto k = 0; k < dim; k++)
        scalar += Delta_[i][k] * Delta_[i][k];
    scalar = sqrt(scalar);
    if (scalar < 1e-8)
        scalar = 1.0;
    scalar = (1.0 / scalar) - 1.0;
    for (auto k = 0; k < dim; k++)
        Delta_[i][k] *= 2 * scalar;

    for (auto j = i + 1; j < n; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-8)
            scalar = 1.0;
        if (D(j, i) > 1e-8)
            scalar = (1.0 / scalar) - (1.0 / D(j, i));
        else
            scalar = (1.0 / scalar) - 1.0;

        for (auto k = 0; k < dim; k++)
            Delta_[j][k] *= 2 * scalar;
    }
}

void NLDR::direction_SAMMON_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    Matrix<PRECISION> &Cor = *(X->Cor);
    Matrix<PRECISION> &DeltaX = *(Dir->Delta);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    auto CX_ = Cor.get_vec();
    auto DeltaX_ = DeltaX.get_vec();
    PRECISION D_ij, DX_ij;

    PRECISION c = paras->total_d;

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            DeltaX_[i][k] = 0;
            for (auto j = 0; j < i; j++)
            {
                D_ij = D(i, j);
                DX_ij = DX(i, j);
                if (DX_ij > 1e-8 && D_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij / D_ij * (CX_[i][k] - CX_[j][k]);
            }

            for (auto j = i + 1; j < n; j++)
            {
                D_ij = D(j, i);
                DX_ij = DX(j, i);
                if (DX_ij > 1e-8 && D_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij / D_ij * (CX_[i][k] - CX_[j][k]);
            }
            DeltaX_[i][k] /= c;
        }
}

void NLDR::direction_SAMMON_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    if (Dir->activated_coo != 0)
        return;

    auto n = X->Cor->get_row();
    auto d = X->Cor->get_col();
    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    auto i = Dir->activated_row;

    Matrix<PRECISION> XD(n, d);
    Matrix<PRECISION> XD2(n, d);
    Array<PRECISION> DTX(n);
    Array<PRECISION> DTX3(n);

    auto XD_ = XD.get_vec();
    auto XD2_ = XD2.get_vec();
    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho), total_d = paras->total_d;

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < n; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        if (DX_ij > 1e-8 && D_ij > 1e-8)
        {
            DTX[j] = (D_ij - DX_ij) / DX_ij / D_ij;
            DTX3[j] = 1.0 / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX[j] = 0;
            DTX3[j] = 0;
        }
    }

    for (auto k = 0; k < d; k++)
    {
        a = 0;
        for (auto j = 0; j < n; j++)
            a += DTX[j];
        b = 0;
        for (auto j = 0; j < n; j++)
            b += DTX[j] * XD_[j][k];
        Delta_[0][k] = -2.0 * b / total_d;
        c = 0;
        for (auto j = 0; j < n; j++)
            c += DTX3[j] * XD2_[j][k];
        Delta_[1][k] = -2.0 * (a - c) / total_d;
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

// SAMMON_GAUSS_SEIDEL

void NLDR::pre_update_by_index_SAMMON_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
{
    Matrix<PRECISION> &C = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    Matrix<PRECISION> &Delta = *(Dir->Delta);

    unsigned dim = C.get_col(), n = C.get_row();
    unsigned i = Dir->activated_row;
    unsigned k = Dir->activated_coo;
    PRECISION e1 = Delta(0, k);
    PRECISION e2 = Delta(1, k);
    PRECISION xcnew = C(i, k) - stepsize * e1 / e2;

    PRECISION *newD_i = Paras->VarDiff->newD_row;
    Paras->VarDiff->activated_row = i;
    Paras->VarDiff->activated_coo = k;
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);

    PRECISION DX_ij, D_ij;
    auto C_ = C.get_vec();

    PRECISION temp = 0;

    for (auto j = 0; j < n; j++)
    {
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            // if (l != k)
            //     temp = C_[i][k]- C_[j][k];
            // else
            //     temp = xcnew - C_[j][k];

            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
    }
}

void NLDR::update_SAMMON_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This update routine is called in Line_search wolfe.
    // Perform not inplace update from VarDiff, no update on cost.
    auto X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
    auto stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    auto Y = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[3]));
    auto Paras = reinterpret_cast<NLDR::NLDR_Paras_GAU *>(kws.carg(0));

    pre_update_by_index_SAMMON_GAUSS_SEIDEL(X, Dir, stepsize, Paras);

    update_distance_row_from_VarDiff(X, Y, Dir, stepsize, Paras);
}

void NLDR::cost_SAMMON_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(1));
    auto &fY = *reinterpret_cast<PRECISION *>(kws.arg(ind[1]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));

    SAMMON_stress_by_index(X, fY, paras);
}

void NLDR::direction_SAMMON_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, unsigned a_row,
                                         PRECISION **XD_, PRECISION **XD2_, PRECISION *DTX_, PRECISION *DTX3_)
{

    auto n = paras->pts;
    auto d = paras->dim;

    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    Dir->activated_coo = 0;
    Dir->activated_row = a_row;
    auto i = a_row;

    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho), total_d = paras->total_d;

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < i; j++)
    {
        DX_ij = DX(i, j);
        D_ij = D(i, j);
        if (DX_ij > 1e-8 && D_ij > 1e-8)
        {
            DTX_[j] = (D_ij - DX_ij) / DX_ij / D_ij;
            DTX3_[j] = 1.0 / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }

    DTX_[i] = 0;
    DTX3_[i] = 0;

    for (auto j = i + 1; j < n; j++)
    {
        DX_ij = DX(j, i);
        D_ij = D(j, i);
        if (DX_ij > 1e-8 && D_ij > 1e-8)
        {
            DTX_[j] = (D_ij - DX_ij) / DX_ij / D_ij;
            DTX3_[j] = 1.0 / DX_ij / DX_ij / DX_ij;
        }
        else
        {
            DTX_[j] = 0;
            DTX3_[j] = 0;
        }
    }

    for (auto k = 0; k < d; k++)
    {
        a = 0;
        for (auto j = 0; j < n; j++)
            a += DTX_[j];
        b = 0;
        for (auto j = 0; j < n; j++)
            b += DTX_[j] * XD_[j][k];
        Delta_[0][k] = -2.0 * b / total_d;
        c = 0;
        for (auto j = 0; j < n; j++)
            c += DTX3_[j] * XD2_[j][k];
        Delta_[1][k] = -2.0 * (a - c) / total_d;
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

void NLDR::update_SAMMON_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This is called by the Optim_Update_Solver.Optim::Iter, ind = x, Dir, stepsize, stress, error

    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    auto &f = *reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    auto &e = *reinterpret_cast<PRECISION *>(kws.arg(ind[4]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));
    auto &wolfe = *paras->wolfe;
    auto &wolfe_kws = *wolfe.get_status();
    auto &slope = *reinterpret_cast<PRECISION *>(wolfe_kws.arg(4));

    auto f_pre = f;
    auto n = X->Cor->get_row(), d = X->Cor->get_col();

    auto XD_ = paras->XD->get_vec(), XD2_ = paras->XD2->get_vec();
    auto DTX_ = paras->DTX->get_vec(), DTX3_ = paras->DTX3->get_vec();

    auto Delta_ = Dir->Delta->get_vec();

    for (auto i = 0; i < n; i++)
    {
        direction_SAMMON_GAUSS_SEIDEL(X, Dir, paras, i, XD_, XD2_, DTX_, DTX3_);
        for (auto k = 0; k < d; k++)
        {
            Dir->activated_coo = k;
            slope = -Delta_[0][k] * Delta_[0][k] / Delta_[1][k];
            wolfe.stepsize();
        }
    }
    centralize(*X->Cor);
    e = (f_pre - f) / f;
    std::cout << "stress value : " << f << '\n';
    std::cout << "rel error : " << e << '\n'
              << '\n';
}

/****************************************SAMMON******************************************/
/****************************************CCA*********************************************/

PRECISION NLDR::compute_lambda(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    PRECISION lambda = 0;
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);

    PRECISION **DX_ = new PRECISION *[n];
    PRECISION **D_ = new PRECISION *[n];
    for (auto i = 0; i < n; i++)
    {
        DX_[i] = DX[i];
        D_[i] = D[i];
    }

    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j <= i; j++)
            if (lambda < (DX_[i][j] - D_[i][j]))
                lambda = (DX_[i][j] - D_[i][j]) / 2.0;
        for (auto j = i + 1; j < n; j++)
            if (lambda < (DX_[j][i] - D_[j][i]))
                lambda = (DX_[j][i] - D_[j][i]) / 2.0;
        // Efficient row scan
    }

    delete[] DX_;
    delete[] D_;

    return lambda;
}

PRECISION NLDR::CCA_stress_by_distance(NLDR_Var *X, NLDR_Paras *Paras)
{
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);

    unsigned dim = Paras->dim;
    unsigned n = DX.dimension();

    PRECISION temp = 0;
    PRECISION stress = 0;
    PRECISION c_l = Paras->lambda[(*Paras->iter)];
    // PRECISION c_l = Paras->cur_lambda;

    for (auto i = 0; i < n; i++)
    {
        PRECISION *DX_i = DX.ele(i);
        PRECISION *D_i = D.ele(i);
        for (auto j = 0; j < i; j++)
        {
            temp = DX_i[j] - D_i[j];
            stress += temp * temp * exp(-DX_i[j] / c_l);
        }

        // for (auto DX_ij = DX.row_begin(i), D_ij = D.row_begin(i), DX_end = DX.row_end(i);
        //     DX_ij != DX_end; DX_ij++, D_ij++)
        //     {
        //         temp = (*DX_ij - *D_ij);
        //         if (*D_ij > 1e-8)
        //         stress += temp * temp * exp(- (*DX_ij) / c_l);
        //     }
    }

    return stress;
}

PRECISION NLDR::CCA_stress_by_distance_lambda(LowerTri<PRECISION> *DX, LowerTri<PRECISION> *D, PRECISION lambda)
{
    PRECISION temp = 0;
    PRECISION stress = 0;
    unsigned n = DX->dimension();
    for (auto i = 0; i < n; i++)
    {
        PRECISION *DX_i = DX->ele(i);
        PRECISION *D_i = D->ele(i);
        for (auto j = 0; j < i; j++)
        {
            temp = DX_i[j] - D_i[j];
            stress += temp * temp * exp(-DX_i[j] / lambda);
        }
    }
}

void NLDR::CCA_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras *Paras)
{
    // This can only be used in CCA_GUASS_SEIDEL
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    PRECISION *newD_i = Paras->VarDiff->newD_row;
    PRECISION temp;
    unsigned i = Paras->VarDiff->activated_row;
    unsigned n = X->Cor->get_row();
    unsigned dim = Paras->dim;

    PRECISION c_l = Paras->lambda[(*Paras->iter) / dim];
    // PRECISION c_l = Paras->cur_lambda;

    for (auto j = 0; j < i; j++)
    {
        auto D_ij = D(i, j);
        auto DX_ij = DX(i, j);
        temp = D_ij - newD_i[j];
        stress += temp * temp * exp(-newD_i[j] / c_l);
        temp = D_ij - DX_ij;
        stress -= temp * temp * exp(-DX_ij / c_l);
    }

    for (auto j = i + 1; j < n; j++)
    {
        auto D_ij = D(j, i);
        auto DX_ij = DX(j, i);
        temp = D_ij - newD_i[j];
        stress += temp * temp * exp(-newD_i[j] / c_l);
        temp = D_ij - DX_ij;
        stress -= temp * temp * exp(-DX_ij / c_l);
    }
}

void NLDR::CCA_stress_by_index(NLDR_Var *X, PRECISION &stress, NLDR_Paras_GAU *Paras)
{
    // This can only be used in CCA_GUASS_SEIDEL
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    PRECISION *newD_i = Paras->VarDiff->newD_row;
    PRECISION temp;
    unsigned i = Paras->VarDiff->activated_row;
    unsigned n = X->Cor->get_row();

    PRECISION c_l = *Paras->cur_lambda;

    for (auto j = 0; j < i; j++)
    {
        auto D_ij = D(i, j);
        auto DX_ij = DX(i, j);
        temp = D_ij - newD_i[j];
        stress += temp * temp * exp(-newD_i[j] / c_l);
        temp = D_ij - DX_ij;
        stress -= temp * temp * exp(-DX_ij / c_l);
    }

    for (auto j = i + 1; j < n; j++)
    {
        auto D_ij = D(j, i);
        auto DX_ij = DX(j, i);
        temp = D_ij - newD_i[j];
        stress += temp * temp * exp(-newD_i[j] / c_l);
        temp = D_ij - DX_ij;
        stress -= temp * temp * exp(-DX_ij / c_l);
    }
}

void NLDR::compute_BX_CCA(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    LowerTri<PRECISION> &BX = *(Paras->BX);
    PRECISION *row_sum = Paras->row_sum->get_vec();
    memset(row_sum, 0, n * sizeof(PRECISION));
    for (auto i = 0; i < n; i++)
    {
        auto DX_i = DX[i];
        auto D_i = D[i];
        auto BX_i = BX[i];
        for (auto j = 0; j < i; j++)
            if (DX_i[j] != 0)
            {

                BX_i[j] = -D_i[j] / DX_i[j] * exp(-DX_i[j] / Paras->lambda[*Paras->iter]);

                row_sum[i] += BX_i[j];
                row_sum[j] += BX_i[j];
            }
            else
                BX_i[j] = 0;
    }

    for (auto i = 0; i < n; i++)
        BX(i, i) = -row_sum[i];
}

void NLDR::compute_V_CCA(NLDR_Var *X, NLDR_Paras *Paras)
{
    auto D = Paras->Dis;
    auto V = Paras->V;
    auto DX = X->Dis;

    auto V_ = V->get_vec();
    auto n = D->dimension();
    for (auto i = 0; i < n; i++)
    {
        auto D_i = (*D)[i];
        auto DX_i = (*DX)[i];
        PRECISION c_l = Paras->lambda[*Paras->iter];
        for (auto j = 0; j < i; j++)
        {
            V_[i][j] = -exp(-DX_i[j] / c_l);
            V_[j][i] = V_[i][j];
            V_[i][i] -= V_[i][j];
            V_[j][j] -= V_[j][i];
        }
    }
    shift(*V, 1.0);
    inverse(*V);
    shift(*V, -1.0 / n / n);
};

void NLDR::update_cor_by_BX_CCA(NLDR_Var *X, NLDR_Paras *Paras)
{
    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    auto &Cor = *(X->Cor);
    auto &BX = *(Paras->BX);

    auto XC_ = Cor.get_vec();
    auto &V = *(Paras->V);
    auto V_ = V.get_vec();

    auto A = Paras->A->get_vec();
    memset(A, 0, n * dim * sizeof(PRECISION));
    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            for (auto j = 0; j <= i; j++)
                A[i * dim + k] += BX(i, j) * XC_[j][k];
            for (auto j = i + 1; j < n; j++)
                A[i * dim + k] += BX(j, i) * XC_[j][k];
            // Efficient scan. (using operator () that converts the indices will bring in multiple conversions and if-check's)
        }

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            XC_[i][k] = 0;
            for (auto j = 0; j < n; j++)
                XC_[i][k] += V_[i][j] * A[j * dim + k];
        }
}

void NLDR::update_CCA_MAJORIZATION(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // ip_Update(X, ..., paras)
    // Taking X, stress, error ; paras

    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    PRECISION *stress = reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    PRECISION *error = reinterpret_cast<PRECISION *>(kws.arg(ind[4]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    PRECISION &cur_lambda = paras->lambda[*(paras->iter)];
    PRECISION computed_lambda = 2 * compute_lambda(X, paras);
    if (cur_lambda < computed_lambda)
        cur_lambda = computed_lambda;

    compute_BX_CCA(X, paras);
    compute_V_CCA(X, paras);
    update_cor_by_BX_CCA(X, paras);
    Cor2Dis_no_paras(X, paras);
    *stress = CCA_stress_by_distance(X, paras);

    *error = (paras->pre_stress - *stress) / (*stress);
    if (*error < 0)
    {
        std::cout << "Stress increased in CCA_MAJORIZATION. Terminate immediately.\n";
        *error = 0;
    }
    paras->pre_stress = *(stress);
}

void NLDR::direction_CCA_STOCHASTIC(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned i = paras->s_ind[*(paras->iter)];
    LowerTri<PRECISION> &D = *(paras->Dis);
    LowerTri<PRECISION> &DX = *(X->Dis);
    Matrix<PRECISION> &CX = *(X->Cor);
    Matrix<PRECISION> &Delta = *(Dir->Delta);
    unsigned n = D.dimension();
    unsigned dim = CX.get_col();

    auto CX_ = CX.get_vec();
    auto Delta_ = Delta.get_vec();

    PRECISION scalar = 0;
    PRECISION lambda = paras->lambda[*(paras->iter)];

    for (auto j = 0; j < n; j++)
        for (auto k = 0; k < dim; k++)
            Delta_[j][k] = CX_[j][k] - CX_[i][k];

    for (auto j = 0; j < i; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-8)
            scalar = 1.0;
        scalar = exp(-scalar / lambda) * (D(i, j) / scalar - 1.0);
        for (auto k = 0; k < dim; k++)
            Delta(j, k) *= scalar;
    }

    scalar = 0;
    for (auto k = 0; k < dim; k++)
        scalar += Delta_[i][k] * Delta_[i][k];
    scalar = sqrt(scalar);
    if (scalar < 1e-8)
        scalar = 1.0;
    scalar = -exp(-scalar / lambda);
    for (auto k = 0; k < dim; k++)
        Delta(i, k) *= scalar;

    for (auto j = i + 1; j < n; j++)
    {
        scalar = 0;
        for (auto k = 0; k < dim; k++)
            scalar += Delta_[j][k] * Delta_[j][k];
        scalar = sqrt(scalar);
        if (scalar < 1e-8)
            scalar = 1.0;
        scalar = exp(-scalar / lambda) * (D(j, i) / scalar - 1.0);
        for (auto k = 0; k < dim; k++)
            Delta(j, k) *= scalar;
    }
}

void NLDR::direction_CCA_METROPOLIS(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    unsigned n = X->Cor->get_row();
    unsigned dim = X->Cor->get_col();
    Matrix<PRECISION> &Cor = *(X->Cor);
    Matrix<PRECISION> &DeltaX = *(Dir->Delta);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);
    PRECISION lambda = paras->lambda[*(paras->iter)];

    auto CX_ = Cor.get_vec();
    auto DeltaX_ = DeltaX.get_vec();
    PRECISION D_ij, DX_ij;

    for (auto i = 0; i < n; i++)
        for (auto k = 0; k < dim; k++)
        {
            DeltaX_[i][k] = 0;
            for (auto j = 0; j < i; j++)
            {
                D_ij = D(i, j);
                DX_ij = DX(i, j);
                if (DX_ij > 1e-8 && D_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij * (CX_[i][k] - CX_[j][k]) * (2.0 + (D_ij - DX_ij) / lambda) * exp(-DX_ij / lambda);
            }

            for (auto j = i + 1; j < n; j++)
            {
                D_ij = D(j, i);
                DX_ij = DX(j, i);
                if (DX_ij > 1e-8 && D_ij > 1e-8)
                    DeltaX_[i][k] += (D_ij - DX_ij) / DX_ij * (CX_[i][k] - CX_[j][k]) * (2.0 + (D_ij - DX_ij) / lambda) * exp(-DX_ij / lambda);
            }
        }
}

void NLDR::direction_CCA_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    NLDR_Var *X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    NLDR_Dir *Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    NLDR_Paras *paras = reinterpret_cast<NLDR_Paras *>(kws.carg(0));

    if (Dir->activated_coo != 0)
        return;

    auto n = X->Cor->get_row();
    auto d = X->Cor->get_col();
    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);
    PRECISION &c_l = paras->lambda[(*paras->iter) / d];
    // auto lambda = paras->lambda;
    // auto cur_iter = *(paras->iter);

    auto i = Dir->activated_row;

    Matrix<PRECISION> XD(n, d);
    Matrix<PRECISION> XD2(n, d);
    Array<PRECISION> DPJ(n);
    Array<PRECISION> DQ(n);
    Array<PRECISION> DR(n);

    auto XD_ = XD.get_vec();
    auto XD2_ = XD2.get_vec();
    auto Delta_ = Dir->Delta->get_vec();

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho), total_d = paras->total_d;

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < n; j++)
    {
        DPJ[j] = 0;
        for (auto k = 0; k < d; k++)
            DPJ[j] += XD2_[j][k];
        DPJ[j] = sqrt(DPJ[j]);
    }

    for (auto j = 0; j < n; j++)
    {
        DQ[j] = D(j, i) - DPJ[j];
        if (c_l < -DQ[j])
            c_l = -DQ[j];
    }

    for (auto j = 0; j < n; j++)
        DR[j] = (2 + DQ[j] / c_l) * exp(-D(j, i) / c_l);

    for (auto j = 0; j < n; j++)
        if (fabs(DPJ[j]) > 1e-7)
            DQ[j] *= (DR[j] / DPJ[j]);
        else
            DQ[j] = 0;

    for (auto k = 0; k < d; k++)
    {
        Delta_[0][k] = 0;
        for (auto j = 0; j < n; j++)
            Delta_[0][k] -= XD_[j][k] * DQ[j];
    }

    for (auto j = 0; j < n; j++)
    {
        D_ij = D(j, i);
        if (fabs(DPJ[j]) > 1e-11)
            DR[j] = (DPJ[j] * DPJ[j] * (2 * D_ij + 3 * c_l - DPJ[j]) - (DPJ[j] + c_l) * D_ij * (D_ij + 2 * c_l)) / DPJ[j] / DPJ[j] / DPJ[j] / c_l / c_l * exp(-D_ij / c_l);
        else
            DR[j] = 0;
    }
    for (auto k = 0; k < d; k++)
    {
        Delta_[1][k] = 0;
        for (auto j = 0; j < n; j++)
            Delta_[1][k] -= (DQ[j] + XD2_[j][k] * DR[j]);
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

// CCA_GAUSS_SEIDEL

void NLDR::pre_update_by_index_CCA_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, PRECISION stepsize, NLDR_Paras_GAU *Paras)
{
    Matrix<PRECISION> &C = *(X->Cor);
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(Paras->Dis);
    Matrix<PRECISION> &Delta = *(Dir->Delta);

    unsigned dim = C.get_col(), n = C.get_row();
    unsigned i = Dir->activated_row;
    unsigned k = Dir->activated_coo;
    PRECISION e1 = Delta(0, k);
    PRECISION e2 = Delta(1, k);
    PRECISION xcnew = C(i, k) - stepsize * e1 / e2;

    PRECISION *newD_i = Paras->VarDiff->newD_row;
    Paras->VarDiff->activated_row = i;
    Paras->VarDiff->activated_coo = k;
    PRECISION eta = *(X->eta);
    PRECISION eta_d = *(X->eta_d);
    PRECISION rho = *(X->rho);

    PRECISION DX_ij, D_ij;
    auto C_ = C.get_vec();

    PRECISION temp = 0;

    for (auto j = 0; j < n; j++)
    {
        newD_i[j] = 0;
        for (auto l = 0; l < dim; l++)
        {
            // if (l != k)
            //     temp = C_[i][k]- C_[j][k];
            // else
            //     temp = xcnew - C_[j][k];

            if (l != k)
                temp = C_[i][l] - C_[j][l];
            else
            {
                if (i != j)
                    temp = xcnew - C_[j][k];
                else
                    temp = 0;
            }
            newD_i[j] += temp * temp;
        }
        newD_i[j] = sqrt(newD_i[j]);
    }
}

void NLDR::update_CCA_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This update routine is called in Line_search wolfe.
    // Perform not inplace update from VarDiff, no update on cost.
    auto X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
    auto stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
    auto Y = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[3]));
    auto Paras = reinterpret_cast<NLDR::NLDR_Paras_GAU *>(kws.carg(0));

    pre_update_by_index_CCA_GAUSS_SEIDEL(X, Dir, stepsize, Paras);

    update_distance_row_from_VarDiff(X, Y, Dir, stepsize, Paras);
}

void NLDR::cost_CCA_GAUSS_SEIDEL_line_search(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(1));
    auto &fY = *reinterpret_cast<PRECISION *>(kws.arg(ind[1]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));

    CCA_stress_by_index(X, fY, paras);
}

void NLDR::direction_CCA_GAUSS_SEIDEL(NLDR_Var *X, NLDR_Dir *Dir, NLDR_Paras_GAU *paras, unsigned a_row,
                                      PRECISION **XD_, PRECISION **XD2_, PRECISION *DPJ_, PRECISION *DQ_, PRECISION *DR_)
{

    auto n = paras->pts;
    auto d = paras->dim;

    auto C_ = X->Cor->get_vec();
    LowerTri<PRECISION> &DX = *(X->Dis);
    LowerTri<PRECISION> &D = *(paras->Dis);

    Dir->activated_coo = 0;
    Dir->activated_row = a_row;
    auto i = a_row;

    auto Delta_ = Dir->Delta->get_vec();
    PRECISION &c_l = *paras->cur_lambda;

    PRECISION a, b, c, DX_ij, D_ij, eta = *(X->eta), eta_d = paras->eta_d, rho = *(X->rho), total_d = paras->total_d;

    for (auto j = 0; j < n; j++)
    {
        for (auto k = 0; k < d; k++)
        {
            XD_[j][k] = C_[i][k] - C_[j][k];
            XD2_[j][k] = XD_[j][k] * XD_[j][k];
        }
    }

    for (auto j = 0; j < n; j++)
    {
        DPJ_[j] = 0;
        for (auto k = 0; k < d; k++)
            DPJ_[j] += XD2_[j][k];
        DPJ_[j] = sqrt(DPJ_[j]);
    }

    for (auto j = 0; j < i; j++)
    {
        DQ_[j] = D(i, j) - DPJ_[j];
        if (c_l < -DQ_[j])
            c_l = -DQ_[j];
    }

    DQ_[i] = -DPJ_[i];
    if (c_l < -DQ_[i])
        c_l = -DQ_[i];

    for (auto j = i + 1; j < n; j++)
    {
        DQ_[j] = D(j, i) - DPJ_[j];
        if (c_l < -DQ_[j])
            c_l = -DQ_[j];
    }

    for (auto j = 0; j < i; j++)
        DR_[j] = (2 + DQ_[j] / c_l) * exp(-D(i, j) / c_l);

    DR_[i] = (2 + DQ_[i] / c_l);

    for (auto j = i + 1; j < n; j++)
        DR_[j] = (2 + DQ_[j] / c_l) * exp(-D(j, i) / c_l);

    for (auto j = 0; j < n; j++)
        if (fabs(DPJ_[j]) > 1e-7)
            DQ_[j] *= (DR_[j] / DPJ_[j]);
        else
            DQ_[j] = 0;

    for (auto k = 0; k < d; k++)
    {
        Delta_[0][k] = 0;
        for (auto j = 0; j < n; j++)
            Delta_[0][k] -= XD_[j][k] * DQ_[j];
    }

    for (auto j = 0; j < i; j++)
    {
        D_ij = D(i, j);
        if (fabs(DPJ_[j]) > 1e-11)
            DR_[j] = (DPJ_[j] * DPJ_[j] * (2 * D_ij + 3 * c_l - DPJ_[j]) - (DPJ_[j] + c_l) * D_ij * (D_ij + 2 * c_l)) / DPJ_[j] / DPJ_[j] / DPJ_[j] / c_l / c_l * exp(-D_ij / c_l);
        else
            DR_[j] = 0;
    }

    // DPJ[i] = 0;
    DR_[i] = 0;

    for (auto j = i + 1; j < n; j++)
    {
        D_ij = D(j, i);
        if (fabs(DPJ_[j]) > 1e-11)
            DR_[j] = (DPJ_[j] * DPJ_[j] * (2 * D_ij + 3 * c_l - DPJ_[j]) - (DPJ_[j] + c_l) * D_ij * (D_ij + 2 * c_l)) / DPJ_[j] / DPJ_[j] / DPJ_[j] / c_l / c_l * exp(-D_ij / c_l);
        else
            DR_[j] = 0;
    }

    for (auto k = 0; k < d; k++)
    {
        Delta_[1][k] = 0;
        for (auto j = 0; j < n; j++)
            Delta_[1][k] -= (DQ_[j] + XD2_[j][k] * DR_[j]);
        Delta_[1][k] = fabs(Delta_[1][k]);
    }
}

void NLDR::update_CCA_GAUSS_SEIDEL(OptimLib::Optim_KwArg &kws, unsigned *ind)
{
    // This is called by the Optim_Update_Solver.Optim::Iter, ind = x, Dir, stepsize, stress, error

    auto X = reinterpret_cast<NLDR_Var *>(kws.arg(ind[0]));
    auto Dir = reinterpret_cast<NLDR_Dir *>(kws.arg(ind[1]));
    auto &f = *reinterpret_cast<PRECISION *>(kws.arg(ind[3]));
    auto &e = *reinterpret_cast<PRECISION *>(kws.arg(ind[4]));

    auto paras = reinterpret_cast<NLDR_Paras_GAU *>(kws.carg(0));
    auto &wolfe = *paras->wolfe;
    auto &wolfe_kws = *wolfe.get_status();
    auto &slope = *reinterpret_cast<PRECISION *>(wolfe_kws.arg(4));

    auto f_pre = f;
    auto n = X->Cor->get_row(), d = X->Cor->get_col();

    auto XD_ = paras->XD->get_vec(), XD2_ = paras->XD2->get_vec();
    auto DPJ_ = paras->DTX->get_vec(), DQ_ = paras->DTX3->get_vec(), DR_ = paras->DR->get_vec();

    auto Delta_ = Dir->Delta->get_vec();

    for (auto i = 0; i < n; i++)
    {
        direction_CCA_GAUSS_SEIDEL(X, Dir, paras, i, XD_, XD2_, DPJ_, DQ_, DR_);
        for (auto k = 0; k < d; k++)
        {
            Dir->activated_coo = k;
            slope = -Delta_[0][k] * Delta_[0][k] / Delta_[1][k];
            wolfe.stepsize();
        }
        paras->cur_lambda++;
    }
    centralize(*X->Cor);
    CCA_stress_by_distance_lambda(X->Dis, paras->Dis, *paras->cur_lambda);
    e = (f_pre - f) / f;
    std::cout << "stress value : " << f << '\n';
    std::cout << "rel error : " << e << '\n'
              << '\n';
}

/****************************************CCA*********************************************/

// void NLDR::CCA_distance_stress_function(NLDR::NLDR_Var *X, PRECISION *f, NLDR::NLDR_Paras *Paras, PRECISION lambda)
// {
//     Matrix<PRECISION> &C = *(X->Cor);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     unsigned dim = Paras->dim;
//     unsigned n = C.get_row();

//     PRECISION temp = 0;
//     PRECISION &stress = *f;
//     stress = 0;

//     for (auto i = 0; i < n; i++)
//     {
//         PRECISION *DX_i = DX[i];
//         PRECISION *D_i = D[i];
//         for (auto j = 0; j < i; j++)
//         {
//             DX_i[j] = 0;
//             for (auto k = 0; k < dim; k++)
//             {
//                 temp = C(i, k) - C(j, k);
//                 DX_i[j] += temp * temp;
//             }
//             DX_i[j] = sqrt(DX_i[j]);
//             temp = DX_i[j] - D_i[j];
//             if (D_i[j] > 1e-8)
//                 stress += temp * temp * exp(-DX_i[j] / lambda);
//         }
//     }
// }

// void NLDR::inplace_update_distance_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, NLDR::NLDR_Paras *Paras, PRECISION stepsize)
// {
//     Matrix<PRECISION> &C = *(X->Cor);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     Matrix<PRECISION> &Delta = *(Dir->Delta);
//     unsigned dim = Delta.get_col(), n = Delta.get_row();
//     unsigned i = Dir->activated_row;
//     unsigned k = Dir->activated_coo;
//     PRECISION xcnew = C(i, k) + stepsize * Delta(i, k);

//     PRECISION temp = 0;

//     for (auto j = 0; j < n; j++)
//     {
//         DX(i, j) = 0;
//         for (auto l = 0; l < dim; l++)
//         {
//             if (l != k)
//                 temp = C(i, k) - C(j, k);
//             else
//                 temp = xcnew - C(j, k);
//             DX(i, j) += temp * temp;
//         }
//         DX(i, j) = sqrt(DX(i, j));
//     }
//     DX(i, i) = 0;
// }

// void NLDR::pre_update_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, PRECISION stepsize, NLDR::NLDR_Paras *Paras)
// {
//     Matrix<PRECISION> &C = *(X->Cor);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     Matrix<PRECISION> &Delta = *(Dir->Delta);

//     unsigned dim = Delta.get_col(), n = Delta.get_row();
//     unsigned i = Dir->activated_row;
//     unsigned k = Dir->activated_coo;
//     PRECISION xcnew = C(i, k) + stepsize * Delta(i, k);

//     PRECISION *newD_i = Paras->VarDiff->newD_row;
//     Paras->VarDiff->activated_row = i;
//     Paras->VarDiff->activated_coo = k;
//     PRECISION eta = *(X->eta);
//     PRECISION eta_d = *(X->eta_d);
//     PRECISION rho = *(X->rho);

//     PRECISION temp = 0;

//     for (auto j = 0; j < n; j++)
//     {
//         eta -= DX(i, j) * DX(i, j);
//         rho -= DX(i, j) * D(i, j);
//         newD_i[j] = 0;
//         for (auto l = 0; l < dim; l++)
//         {
//             if (l != k)
//                 temp = C(i, k) - C(j, k);
//             else
//                 temp = xcnew - C(j, k);
//             newD_i[j] += temp * temp;
//         }
//         eta += newD_i[j];
//         newD_i[j] = sqrt(newD_i[j]);
//         rho += newD_i[j] * D(i, j);
//     }
//     newD_i[i] = 0;
//     Paras->VarDiff->eta = eta;
//     Paras->VarDiff->eta_d = eta_d;
//     Paras->VarDiff->rho = rho;
// }

// void NLDR::inplace_update_by_index(OptimLib::Optim_KwArg &kws, unsigned *ind);

// void NLDR::KRUSKAL1_update_stress_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, PRECISION stepsize, PRECISION &stress, NLDR::NLDR_Paras *Paras)
// {
//     NLDR::pre_update_by_index(X, Dir, stepsize, Paras);
//     // Update info is computed and stored at Paras->VarDiff.

//     PRECISION eta, eta_d, rho;

//     if (Dir->activated_coo == 0)
//         Paras->pre_stress = stress;

//     eta = Paras->VarDiff->eta;
//     eta_d = Paras->VarDiff->eta_d;
//     rho = Paras->VarDiff->rho;
//     stress = (eta + eta_d - 2 * rho) / eta;
// }

// void NLDR::NORMALIZED_update_stress_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, PRECISION stepsize, PRECISION &stress, NLDR::NLDR_Paras *Paras)
// {
//     NLDR::pre_update_by_index(X, Dir, stepsize, Paras);

//     if (Dir->activated_coo == 0)
//         Paras->pre_stress = stress;

//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     PRECISION *newD_i = Paras->VarDiff->newD_row;
//     PRECISION weight = Paras->weight;
//     unsigned i = Dir->activated_row;
//     unsigned n = Dir->Delta->get_row();

//     for (auto j = 0; j < n; j++)
//         if (j != i)
//             stress += (2 * D(i, j) - newD_i[j] - DX(i, j)) * (DX(i, j) - newD_i[j]) / weight;
// }

// void NLDR::SAMMON_update_stress_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, PRECISION stepsize, PRECISION &stress, NLDR::NLDR_Paras *Paras)
// {
//     NLDR::pre_update_by_index(X, Dir, stepsize, Paras);

//     if (Dir->activated_coo == 0)
//         Paras->pre_stress = stress;

//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     PRECISION *newD_i = Paras->VarDiff->newD_row;
//     PRECISION weight = Paras->weight;
//     unsigned i = Dir->activated_row;
//     unsigned n = Dir->Delta->get_row();

//     for (auto j = 0; j < n; j++)
//         if (j != i && D(i, j) >= 1e-10)
//             stress += (2 * D(i, j) - newD_i[j] - DX(i, j)) * (D(i, j) - newD_i[j]) / D(i, j) / weight;
// }

// void NLDR::CCA_update_stress_by_index(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Dir, PRECISION stepsize, PRECISION &stress, NLDR::NLDR_Paras *Paras)
// {
//     NLDR::pre_update_by_index(X, Dir, stepsize, Paras);

//     if (Dir->activated_coo == 0)
//         Paras->pre_stress = stress;

//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     PRECISION *newD_i = Paras->VarDiff->newD_row;
//     PRECISION weight = Paras->weight;
//     PRECISION temp;
//     unsigned i = Dir->activated_row;
//     unsigned n = Dir->Delta->get_row();

//     for (auto j = 0; j < n; j++)
//     {
//         if (j != i)
//         {
//             temp = D(i, j) - newD_i[j];
//             stress += temp * temp * exp(-newD_i[j] / weight);
//             temp = D(i, j) - D(i, j);
//             stress -= temp * temp * exp(-D(i, j) / weight);
//         }
//     }
// }

// void NLDR::inplace_update_with_pre_update_by_index(OptimLib::Optim_KwArg &kws, unsigned *ind)
// {
//     // Perform update from VarDiff, no update on cost.
//     NLDR::NLDR_Var *X = reinterpret_cast<NLDR::NLDR_Var *>(kws.arg(ind[0]));
//     NLDR::NLDR_Dir *Delta = reinterpret_cast<NLDR::NLDR_Dir *>(kws.arg(ind[1]));
//     PRECISION stepsize = *reinterpret_cast<PRECISION *>(kws.arg(ind[2]));
//     PRECISION stress = *reinterpret_cast<PRECISION *>(kws.arg(ind[3])); // Already updated
//     NLDR::NLDR_Paras *Paras = reinterpret_cast<NLDR::NLDR_Paras *>(kws.carg(0));

//     unsigned i = Paras->VarDiff->activated_row;
//     unsigned k = Paras->VarDiff->activated_coo;
//     unsigned n = Paras->dim;

//     (*X->Cor)(i, k) += stepsize * (*Delta->Delta)(0, k);
//     for (auto j = 0; j < n; j++)
//         (*X->Dis)(i, j) = Paras->VarDiff->newD_row[j];
//     *X->eta = Paras->VarDiff->eta;
//     *X->eta_d = Paras->VarDiff->eta_d;
//     *X->rho = Paras->VarDiff->rho;
// }

// void NLDR::update_cor_dis_by_BX(NLDR::NLDR_Var *X, const LowerTri<PRECISION> &BX)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     Matrix<PRECISION> C_temp(Cor);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             Cor(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//                 Cor(i, k) += BX.ele(i, j) * C_temp(j, k);
//             Cor(i, k) /= n;
//         }
//     NLDR::Cor2Dis(X->Cor, X->Dis);
// }

// void NLDR::update_cor_dis_by_BX_V(NLDR::NLDR_Var *X, const LowerTri<PRECISION> &BX, const Matrix<PRECISION> &V)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     PRECISION *A = new PRECISION[n * dim];
//     memset(A, 0, n * dim * sizeof(PRECISION));
//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//             for (auto j = 0; j < m; j++)
//                 A[i * dim + k] += BX.ele(i, j) * Cor.ele(j, k);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             Cor(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//                 Cor(i, k) += V.ele(i, j) * A[j * dim + k];
//         }
//     NLDR::Cor2Dis(X->Cor, X->Dis);
//     delete[] A;
// }

// void NLDR::update_cor_by_BX_V(NLDR::NLDR_Var *X, const LowerTri<PRECISION> &BX, const Matrix<PRECISION> &V)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     PRECISION *A = new PRECISION[n * dim];
//     memset(A, 0, n * dim * sizeof(PRECISION));
//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//             for (auto j = 0; j < m; j++)
//                 A[i * dim + k] += BX.ele(i, j) * Cor.ele(j, k);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             Cor(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//                 Cor(i, k) += V.ele(i, j) * A[j * dim + k];
//         }
//     delete[] A;
// }

// void NLDR::SAMMOM_compute_BX(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras, LowerTri<PRECISION> &BX)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     PRECISION *row_sum = new PRECISION[n];
//     memset(row_sum, 0, n * sizeof(PRECISION));

//     for (auto i = 0; i < n; i++)
//         for (auto j = 0; j < i; j++)
//             if (DX(i, j) != 0)
//             {
//                 BX(i, j) = -1.0 / DX(i, j);
//                 row_sum[i] += BX(i, j);
//                 row_sum[j] += BX(i, j);
//             }

//     for (auto i = 0; i < n; i++)
//         BX(i, i) = -row_sum[i];
// }

// void NLDR::SAMMOM_compute_BX(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras, LowerTri<PRECISION> &BX)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     PRECISION *row_sum = new PRECISION[n];
//     memset(row_sum, 0, n * sizeof(PRECISION));

//     for (auto i = 0; i < n; i++)
//         for (auto j = 0; j < i; j++)
//             if (DX(i, j) != 0)
//             {
//                 BX(i, j) = -1.0 / DX(i, j);
//                 row_sum[i] += BX(i, j);
//                 row_sum[j] += BX(i, j);
//             }

//     for (auto i = 0; i < n; i++)
//         BX(i, i) = -row_sum[i];
// }

// void NLDR::CCA_compute_BX(NLDR::NLDR_Var *X, NLDR::NLDR_Paras *Paras, LowerTri<PRECISION> &BX, PRECISION lambda)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     PRECISION *row_sum = new PRECISION[n];
//     memset(row_sum, 0, n * sizeof(PRECISION));

//     for (auto i = 0; i < n; i++)
//         for (auto j = 0; j < i; j++)
//             if (DX(i, j) != 0)
//             {
//                 BX(i, j) = -D(i, j) / DX(i, j) * exp(-DX(i, j) / lambda);
//                 row_sum[i] += BX(i, j);
//                 row_sum[j] += BX(i, j);
//             }

//     for (auto i = 0; i < n; i++)
//         BX(i, i) = -row_sum[i];
// }

// void NLDR::KRUSKAL1_Gradient(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Delta, NLDR::NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     PRECISION n1, n2;
//     double eta = *(X->eta);
//     double eta_d = *(X->eta_d);
//     double rho = *(X->rho);
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     Matrix<PRECISION> &DeltaX = *(Delta->Delta);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             n1 = 0;
//             n2 = 0;
//             for (auto j = 0; j < n; j++)
//             {
//                 n1 += Cor(i, k) - Cor(j, k);
//                 if (DX(i, j) > 1e-8)
//                     n2 += D(i, j) / DX(i, j) * (Cor(i, k) - Cor(j, k));
//             }
//             DeltaX(i, k) = -2.0 / eta / eta * (n1 * eta + eta_d * n2 - 2 * rho * n2);
//         }
// }

// void NLDR::NORMALIZED_Gradient(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Delta, NLDR::NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     Matrix<PRECISION> &DeltaX = *(Delta->Delta);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);

//     PRECISION c = 0;
//     for (auto i = 0; i < n; i++)
//         for (auto j = 0; j < i; j++)
//             c += D(i, j) * D(i, j);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             DeltaX(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//             {

//                 if (DX(i, j) > 1e-8)

//                     DeltaX(i, k) += (DX(i, j) - D(i, j)) / DX(i, j) * (Cor(i, k) - Cor(j, k));
//             }
//             DeltaX(i, k) /= c;
//         }
// }

// void NLDR::SAMMON_Gradient(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Delta, NLDR::NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     Matrix<PRECISION> &DeltaX = *(Delta->Delta);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);

//     PRECISION c = 0;
//     for (auto i = 0; i < n; i++)
//         for (auto j = 0; j < i; j++)
//             c += D(i, j);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             DeltaX(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//             {

//                 if (DX(i, j) > 1e-8 && D(i, j) > 1e-8)
//                     DeltaX(i, k) += (DX(i, j) - D(i, j)) / DX(i, j) / D(i, j) * (Cor(i, k) - Cor(j, k));
//             }
//             DeltaX(i, k) /= c;
//         }
// }

// void NLDR::CCA_Gradient(NLDR::NLDR_Var *X, NLDR::NLDR_Dir *Delta, NLDR::NLDR_Paras *Paras)
// {
//     unsigned n = X->Cor->get_row();
//     unsigned dim = X->Cor->get_col();
//     Matrix<PRECISION> &Cor = *(X->Cor);
//     Matrix<PRECISION> &DeltaX = *(Delta->Delta);
//     LowerTri<PRECISION> &DX = *(X->Dis);
//     LowerTri<PRECISION> &D = *(Paras->Dis);
//     PRECISION lambda = *(Paras->lambda);

//     for (auto i = 0; i < n; i++)
//         for (auto k = 0; k < dim; k++)
//         {
//             DeltaX(i, k) = 0;
//             for (auto j = 0; j < n; j++)
//                 if (DX(i, j) > 1e-8 && D(i, j) > 1e-8)
//                     DeltaX(i, k) += (DX(i, j) - D(i, j)) / DX(i, j) * (Cor(i, k) - Cor(j, k)) * (2.0 + (D(i, j) - DX(i, j)) / lambda) * exp(-DX(i, j) / lambda);
//         }
// }
