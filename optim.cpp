#include "optim.hpp"

PRECISION OptimLib::Wolfe_1st_StepSize::stepsize()
{
    PRECISION &f = *reinterpret_cast<PRECISION *>(kws->arg(2));
    PRECISION slope = *reinterpret_cast<PRECISION *>(kws->arg(4));

    PRECISION a = 0, b = 0, disc = 0;
    PRECISION &fn = *reinterpret_cast<PRECISION *>(kws->arg(7));
    PRECISION &l = *reinterpret_cast<PRECISION *>(kws->arg(5));

    fn = f;
    l = init_l; //  Initialize stepsize
    PRECISION fnp = f, l_temp = 0, lp = init_l;

    while (true)
    {
        fn = f;

        (*U)(*kws, u_ind);     // Update point.
        (*C)(*kws, new_c_ind); // Update cost function. This may be done in update process. In such cases, C should be nullfunction.

        PRECISION t = l;

        //std::cout << "Stepsize:\t" << t << ", Slope:\t" << slope << ", Attempted cost:\t" << *fn << ", Difference:\t" << (*f - *fn) << '\n';

        if (fn <= f + this->alpha * t * slope)
        {
            // Found
            (*CP)(*kws, cp_rev_ind); // Write the updated result back to X for next step.
            f = fn;
            return l;
        }
        else if (t < min_l)
        {
            // Not found. Using minimal step
            //std::cout << "Warning! OptimLib::Wolfe_1st_StepSize fails to find qualified stepsize, using default: " << min_l << ".\n";
            fn = f;
            l = min_l;              // Use minimal stepsize
            (*U)(*kws, u_ind);       // Update point.
            (*C)(*kws, new_c_ind);   // Update cost function.
            (*CP)(*kws, cp_rev_ind); // Write the updated result back to X for next step.
            f = fn;
            return l;
        }
        else
        {

            if (t == lp)
                l_temp = -slope / (2 * (fn - f - slope));
            else
            {
                a = ((fn - f - t * slope) / t / t - (fnp - f - lp * slope) / lp / lp) / (t - lp);
                b = (-lp * (fn - f - t * slope) / t / t + t * (fnp - f - lp * slope) / lp / lp) / (t - lp);

                disc = b * b - 3 * a * slope;
                if (fabs(a) < 1e-7)
                    l_temp = -slope / (2 * b);
                else
                    l_temp = (sqrt(disc) - b) / (3 * a);

                if (l_temp > 0.5 * t)
                    l_temp = 0.5 * t;
            }
            lp = t;
            fnp = fn;
            if (l_temp <= 0.1 * t)
                l = 0.1 * t;
            else
                l = l_temp;
        }
    }

    return -1;
}

void OptimLib::Optim_Iter::next_step(Optim_StepSize *S)
        {
            D(*kws, d_ind);                 // Get update direction;
            SL(*kws, sl_ind);               // Get slope;
            *cur_step = S->stepsize();      // Get stepsize;
            U(*kws, u_ip_ind);              // Update solution;
            *diff = (*cur_step) * (*nD);    // Compute update difference
            C(*kws, c_x_ind);               // Update cost function
            // For solver using inplace update, the error update should be done within inplace update or cost evaluation.
        }

void OptimLib::Optim_Iter::next_step_no_stepsize()
        {
            D(*kws, d_ind);                 // Get update direction;
            U(*kws, u_ip_ind);              // Update solution;
            C(*kws, c_x_ind);               // Update cost function
            // For solver using inplace update, the error update should be done within inplace update or cost evaluation.
        }

void OptimLib::Optim_Iter::trial_step(Optim_StepSize *S)
{
    D(*kws, d_ind);
    *accepted = false;
    SL(*kws, sl_ind);
    if (S)
        *cur_step = S->stepsize();
    U(*kws, u_nip_ind);
    if (nD)
        *diff = (*cur_step) * (*nD);
    C(*kws, c_y_ind);
};

void OptimLib::Optim_Iter::accept_trial_step()
{
    PRECISION fx = *(reinterpret_cast<PRECISION *>(kws->arg(5)));
    PRECISION fy = *(reinterpret_cast<PRECISION *>(kws->arg(11)));
    if (fabs(fx) > 1e-8)
        *err = (fx - fy) / fx;
    CP(*kws, cp_ind_rev);
};

bool OptimLib::Optim_Solver::terminate()
{
    // std::cout << ((*cur_iter) >= paras->MaxIter && paras->MaxIter) << "\t" << paras->MaxIter << "\t" << *cur_iter << "\n";
    // std::cout << ((*cur_time) / (double)CLOCKS_PER_SEC > paras->MaxTime && paras->MaxTime) << "\t" << paras->MaxTime << "\t" << *cur_time << "\n";
    // std::cout << ((*cur_err) < paras->ErrTol && paras->ErrTol > 0) << "\t" << paras->ErrTol << "\t" << *cur_err << "\n";
    // std::cout << ((*cur_diff) < paras->DifTol && paras->DifTol > 0) << "\t" << paras->DifTol << "\t" << *cur_diff << "\n";
    return ((*cur_iter) >= paras->MaxIter && paras->MaxIter) ||
           ((*cur_time) / (double)CLOCKS_PER_SEC > paras->MaxTime && paras->MaxTime) ||
           (fabs(*cur_err) < paras->ErrTol && paras->ErrTol > 0) ||
           (fabs(*cur_diff) < paras->DifTol && paras->DifTol > 0);
}

void OptimLib::Optim_Update_Solver::solve()
{

    time_t start, end;

    while (!terminate())
    {
        start = clock();
        iter->next_step_no_stepsize();
        end = clock();
        *cur_time += (end - start);
        *cur_iter += 1;
    }
}

void OptimLib::Optim_Stepping_Solver::solve()
{

    time_t start, end;

    // PRECISION *f = reinterpret_cast<double *>(cur_states->arg(4));

    while (!terminate())
    {
        start = clock();
        iter->next_step(ss);
        end = clock();
        *cur_time += (end - start);
        *cur_iter += 1;

        // std::cout << *cur_iter << "\t" << *f << "\n";
    }
}

void OptimLib::Optim_Gauss_Seidel_Solver::solve()
{

    time_t start, end;

    // PRECISION *f = reinterpret_cast<double *>(cur_states->arg(4));

    while (!terminate())
    {
        start = clock();
        do 
        {
            iter->next_step(ss);
            *inner_iter += 1;
        }
        while (Loop_Check(*this->cur_states));
        end = clock();
        *cur_time += (end - start);
        *cur_iter += 1;


        // std::cout << *cur_iter << "\t" << *f << "\n";
    }
}


void OptimLib::Optim_Stepping_Solver::one_step()
{

    time_t start, end;

    // PRECISION *f = reinterpret_cast<double *>(cur_states->arg(4));

    start = clock();
    iter->next_step(ss);
    end = clock();
    *cur_time += (end - start);
    *cur_iter += 1;

    // std::cout << *cur_iter << "\t" << *f << "\n";
}


void OptimLib::Optim_Trial_Stepping_Solver::solve()
{

    time_t start, end;

    while (!terminate())
    {
        start = clock();
        iter->trial_step(ss);
        *accepted = A(*cur_states, ind_A);
        if (*accepted)
            iter->accept_trial_step();
        end = clock();
        *cur_time += (end - start);
        *cur_iter += 1;
    }
}

void OptimLib::nullOptimFunc(Optim_KwArg &, int *){};
