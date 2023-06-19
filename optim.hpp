#pragma once

#undef min
#undef max
#undef abs

#include <ctime>
#include <random>
#include <iostream>
#include <fstream>

#include "array.hpp"
#include "kwarg.hpp"

#define DEFAULT_OPTIM_MAXITER 1000
#define DEFAULT_OPTIM_MAXTIME 100
#define DEFAULT_OPTIM_ERRTOL 1e-6
#define DEFAULT_OPTIM_DIFTOL 1e-6

// This libiary require a strict rules in ordering pointers. In general, type KwArg kept null pointers to all parameters needed for specific objects.
// For objects in optimization library, the following partitions on pointers holds:
//      1) The first part keeps all methods/function needed for this object. They all take KwArg type and several unsigned indices as argument.
//      2) The second part keeps variable/points and all associated mutable data structure

#define PROBLEM_ROUT_NUM 6
#define PROBLEM_VARI_NUM 0
#define GENERAL_SEARCH_ROUT_NUM 1
#define GENERAL_SEARCH_VARI_NUM 4
#define WOLFE_SEARCH_ROUT_NUM 3
#define WOLFE_SEARCH_VARI_NUM 8
#define ITERATION_ROUT_NUM 5
#define ITERATION_VARI_NUM 13

enum OPTIM_TYPE
{
    OPTIM_PROB,
    NULL_STEPSIZE,
    CONST_STEPSIZE,
    COMPUTED_STEPSIZE,
    GENERAL_STEPSIZE_SEARCH,
    WOLFE_1ST_STEPSIZE_SEARCH,
    OPTIM_ITERATION
};

namespace OptimLib
{
    class Optim_KwArg
    {

    public:
        KwArg routine;
        KwArg arg;
        unsigned int CUS_IND;

        Optim_KwArg() : routine(KwArg()), arg(KwArg()), CUS_IND(0){};

        Optim_KwArg(unsigned int m_num, unsigned int a_num, unsigned int ca_num)
        {
            routine = KwArg(m_num);
            arg = KwArg(a_num + ca_num);
            CUS_IND = a_num;
        }

        Optim_KwArg(OPTIM_TYPE ot, unsigned int cus_var_num)
        {
            unsigned int m_n, a_n, ca_n;
            ca_n = cus_var_num;
            switch (ot)
            {
            case OPTIM_PROB:
                m_n = PROBLEM_ROUT_NUM;
                a_n = PROBLEM_VARI_NUM;
                break;
            case GENERAL_STEPSIZE_SEARCH:
                m_n = GENERAL_SEARCH_ROUT_NUM;
                a_n = GENERAL_SEARCH_VARI_NUM;
                break;
            case WOLFE_1ST_STEPSIZE_SEARCH:
                m_n = WOLFE_SEARCH_ROUT_NUM;
                a_n = WOLFE_SEARCH_VARI_NUM;
                break;
            case OPTIM_ITERATION:
                m_n = ITERATION_ROUT_NUM;
                a_n = ITERATION_VARI_NUM;
                break;
            default:
                std::cout << "Error! Fail to initialize super type Optim_KwArg for stepsize object.\n";
                throw(1);
            }

            routine = KwArg(m_n);
            arg = KwArg(a_n + ca_n);
            CUS_IND = a_n;
        };

        void init_iter_routine(void *cost, void *dir, void *ip_update, void *slope, void *copy)
        {
            this->set_routine(cost, 0);
            this->set_routine(dir, 1);
            this->set_routine(ip_update, 2);
            this->set_routine(slope, 3);
            this->set_routine(copy, 4);
        };

        void init_iter_variable(unsigned int *c_iter, unsigned int *c_time, PRECISION *c_err, PRECISION *c_diff,
                                void *x, PRECISION *f, void *D, PRECISION *slope, PRECISION *D_norm, PRECISION *stepsize)
        {
            this->set_arg(c_iter, 0);
            this->set_arg(c_time, 1);
            this->set_arg(c_err, 2);
            this->set_arg(c_diff, 3);
            this->set_arg(x, 4);
            this->set_arg(f, 5);
            this->set_arg(D, 6);
            this->set_arg(slope, 7);
            this->set_arg(D_norm, 8);
            this->set_arg(stepsize, 9);
        };

        void init_iter_variable(unsigned int *c_iter, unsigned int *c_time, PRECISION *c_err, PRECISION *c_diff,
                                void *x, PRECISION *fx, void *D, PRECISION *slope, PRECISION *D_norm, PRECISION *stepsize,
                                void *y, PRECISION *fy, bool *accepted)
        {
            this->set_arg(c_iter, 0);
            this->set_arg(c_time, 1);
            this->set_arg(c_err, 2);
            this->set_arg(c_diff, 3);
            this->set_arg(x, 4);
            this->set_arg(fx, 5);
            this->set_arg(D, 6);
            this->set_arg(slope, 7);
            this->set_arg(D_norm, 8);
            this->set_arg(stepsize, 9);
            this->set_arg(y, 10);
            this->set_arg(fy, 11);
            this->set_arg(accepted, 12);
        };

        void init_wolfe1st_routine(void *cost, void *g_update, void *copy)
        {
            this->set_routine(cost, 0);
            this->set_routine(g_update, 1);
            this->set_routine(copy, 2);
        };

        void init_wolfe1st_variable(PRECISION *c_err, void *x, PRECISION *f, void *D,
                                    PRECISION *slope, PRECISION *step, void *y, PRECISION *f_y)
        {
            this->set_arg(c_err, 0);
            this->set_arg(x, 1);
            this->set_arg(f, 2);
            this->set_arg(D, 3);
            this->set_arg(slope, 4);
            this->set_arg(step, 5);
            this->set_arg(y, 6);
            this->set_arg(f_y, 7);
        }

        void init_general_search_routine(void *search) { this->set_routine(search, 0); };

        void init_general_search_variable(void *x, void *D, PRECISION *stepsize)
        {
            this->set_arg(x, 0);
            this->set_arg(D, 1);
            this->set_arg(stepsize, 2);
        };

        void init_prob_routine(void *cost, void *dir, void *ip_update, void *g_update, void *slope, void *copy)
        {
            this->set_routine(cost, 0);
            this->set_routine(dir, 1);
            this->set_routine(ip_update, 2);
            this->set_routine(g_update, 3);
            this->set_routine(slope, 4);
            this->set_routine(copy, 5);
        };

        void *carg(unsigned i) { return this->arg[CUS_IND + i]; };

        template <class T>
        void set_routine(T *arg, unsigned i) { routine.set(arg, i); };
        template <class T>
        void set_arg(T *arg_ptr, unsigned i) { arg.set(arg_ptr, i); };
        template <class T>
        void set_carg(T *arg_ptr, unsigned i) { arg.set(arg_ptr, i + CUS_IND); };
    };

    class Optim_StepSize
    {
    protected:
        OPTIM_TYPE s_type;

    public:
        Optim_StepSize() : s_type(NULL_STEPSIZE){};
        Optim_StepSize(OPTIM_TYPE st) : s_type(st){};

        virtual PRECISION stepsize() = 0;
    };
    class Optim_Const_StepSize : public Optim_StepSize
    {
    protected:
        PRECISION a;

    public:
        Optim_Const_StepSize() : Optim_StepSize(CONST_STEPSIZE), a(1.0){};
        Optim_Const_StepSize(PRECISION x) : Optim_StepSize(CONST_STEPSIZE), a(x){};

        PRECISION stepsize()
        {
            return a;
        };
    };
    class Optim_Computed_StepSize : public Optim_StepSize
    {
    protected:
        PRECISION *a;
        unsigned ind;

    public:
        Optim_Computed_StepSize(PRECISION *x) : Optim_StepSize(COMPUTED_STEPSIZE), a(x), ind(0){};

        PRECISION stepsize()
        {
            ind++;
            return a[ind - 1];
        };
    };
    class Optim_Search_Stepsize : public Optim_StepSize
    {
        typedef PRECISION (*Search)(Optim_KwArg &, unsigned int *ind);

    protected:
        Optim_KwArg *kws;
        // kws must maintain the following order of pointers:
        // method :     Search
        // arg :        X, Delta, step, others.
        // Note that Cost and Error may be different from the general Cost and Error update formula,
        // if there is an efficient implemention restricted to given direction, provided with current status.
    private:
        Search search;
        unsigned int ind_stepsize[3] = {0, 1, 2};
        // For general Optim_Search_Stepsize, the search implementation must be provided from outside.
        // Note that this member wont be inherited to derived Stepsize structure, as they usually have
        // a more structual searching model implemented themselves.

    public:
        Optim_Search_Stepsize(Optim_KwArg *kw, OPTIM_TYPE ost, Search s) : Optim_StepSize(ost), kws(kw), search(s){};
        Optim_Search_Stepsize(Optim_KwArg *kw, Search s) : Optim_Search_Stepsize(kw, GENERAL_STEPSIZE_SEARCH, s){};
        PRECISION stepsize() { return search((*kws), ind_stepsize); };
        Optim_KwArg *get_status() { return kws; };
    };
    class Wolfe_1st_StepSize : public Optim_Search_Stepsize
    {
        typedef void (*Cost)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Update)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Copy)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Slope)(Optim_KwArg &, unsigned int *ind);

    private:
        unsigned int new_c_ind[2] = {6, 7};            // Cost value of y: y, f_y
        unsigned int cp_ind[2] = {1, 6};               // Copy x to y: x, y
        unsigned int cp_rev_ind[2] = {6, 1};           // Copy y to x: y, x
        unsigned int u_ind[7] = {1, 3, 5, 6, 2, 7, 0}; // general update: x, Delta, t, y, f_y, e

    protected:
        PRECISION alpha;
        PRECISION init_l;
        PRECISION min_l;
        Cost C;
        Update U;
        Copy CP;

        // kws must maintain the following order of pointers:
        // method :     Cost, Update, Var_Copy
        // arg :        error, diff, X, f, Delta, initial_slope, step, X_new, f_new , others.
        // Note that the U must not be inplace update, as the search process perform different updates.

    public:
        Wolfe_1st_StepSize(Optim_KwArg *kw, PRECISION a, PRECISION i_l, PRECISION m_l) : Optim_Search_Stepsize(kw, WOLFE_1ST_STEPSIZE_SEARCH, nullptr),
                                                                                         alpha(a), init_l(i_l), min_l(m_l)
        {
            C = reinterpret_cast<Cost>(kw->routine[0]);
            U = reinterpret_cast<Update>(kw->routine[1]);
            CP = reinterpret_cast<Copy>(kw->routine[2]);
        };

        PRECISION stepsize();
    };

    class Optim_Iter
    {
        typedef void (*Cost)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Direction)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Update)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Slope)(Optim_KwArg &, unsigned int *ind);
        typedef void (*Copy)(Optim_KwArg &, unsigned int *ind);

    private:
        unsigned c_x_ind[2] = {4, 5};                    // Cost of x: x, fx
        unsigned c_y_ind[2] = {10, 11};                  // Cost of y: y fy
        unsigned d_ind[3] = {4, 6, 8};                   // Get update direction: x, Delta, n
        unsigned u_ip_ind[5] = {4, 6, 9, 5, 2};          // Inplace update: x, Delta, t, f, e
        unsigned u_nip_ind[7] = {4, 6, 9, 10, 5, 11, 2}; // Inplace update: x, Delta, t, y, fx, fy, e
        unsigned sl_ind[3] = {4, 6, 7};                  // Compute slope: x, Delta, s
        unsigned cp_ind[4] = {4, 5, 10, 11};             // Copy x to y : x, fx, y, fy
        unsigned cp_ind_rev[4] = {10, 11, 4, 5};         // Copy y to x : y, fy, x, fx

    protected:
        Optim_KwArg *kws;

        Cost C;
        Direction D;
        Update U;
        Slope SL;
        Copy CP;

        unsigned int *iter;
        unsigned int *acc_time;
        PRECISION *err;
        PRECISION *diff;
        PRECISION *cur_step;
        PRECISION *nD;
        bool *accepted;

        // kws must maintain the following order of pointers:
        // method :     Cost, Direction, Update, Slope, Copy
        // arg :        cur_iter, cur_time, cur_err, cur_diff,  X, f, Delta, slope, |Delta|, stepsize, Y, fY, others.
        // Note that the update function here should be inplace update on X, if needed.
        // Note that update function should also update the cost value and relative error.

    public:
        Optim_Iter(Optim_KwArg *kw) : kws(kw)
        {
            C = reinterpret_cast<Cost>(kws->routine(0));
            D = reinterpret_cast<Direction>(kws->routine(1));
            U = reinterpret_cast<Update>(kws->routine(2));
            SL = reinterpret_cast<Slope>(kws->routine(3));
            CP = reinterpret_cast<Copy>(kws->routine(4));

            iter = reinterpret_cast<unsigned int *>(kws->arg(0));
            acc_time = reinterpret_cast<unsigned int *>(kws->arg(1));
            err = reinterpret_cast<PRECISION *>(kws->arg(2));
            diff = reinterpret_cast<PRECISION *>(kws->arg(3));
            nD = reinterpret_cast<PRECISION *>(kws->arg(8));
            cur_step = reinterpret_cast<PRECISION *>(kws->arg(9));
            accepted = reinterpret_cast<bool *>(kws->arg(12));
        }
        ~Optim_Iter(){};

        void next_step(Optim_StepSize *S);

        void next_step_no_stepsize();

        void trial_step(Optim_StepSize *S);

        void trial_step_no_stepsize(Optim_StepSize *S);

        void accept_trial_step();

        Optim_KwArg *get_current_states() { return kws; };
        unsigned int *get_cur_iter() { return iter; };
        unsigned int *get_cur_time() { return acc_time; };
        PRECISION *get_cur_err() { return err; };
        PRECISION *get_cur_diff() { return diff; };

        void set_cur_iter(unsigned int *arg)
        {
            kws->set_arg(arg, 0);
            iter = arg;
        };
        void set_cur_time(unsigned int *arg)
        {
            kws->set_arg(arg, 1);
            acc_time = arg;
        };
        void set_cur_err(PRECISION *arg)
        {
            kws->set_arg(arg, 2);
            err = arg;
        };
        void set_cur_diff(PRECISION *arg)
        {
            kws->set_arg(arg, 3);
            diff = arg;
        };
    };
    class Optim_Paras
    {
    public:
        int MaxIter;
        int MaxTime;
        PRECISION ErrTol;
        PRECISION DifTol;

        Optim_Paras() : MaxIter(DEFAULT_OPTIM_MAXITER), MaxTime(DEFAULT_OPTIM_MAXTIME), ErrTol(DEFAULT_OPTIM_ERRTOL),
                        DifTol(DEFAULT_OPTIM_DIFTOL){};
        Optim_Paras(int i, int t, PRECISION e, PRECISION d) : ErrTol(e == 0 ? DEFAULT_OPTIM_ERRTOL : e), DifTol(d == 0 ? DEFAULT_OPTIM_DIFTOL : d)
        {
            MaxIter = (i < 0) ? 0 : ((i == 0) ? DEFAULT_OPTIM_MAXITER : i);
            MaxTime = (t < 0) ? 0 : ((t == 0) ? DEFAULT_OPTIM_MAXTIME : t);
        };
    };
    class Optim_Solver
    {
    protected:
        Optim_Paras *paras;
        Optim_Iter *iter;
        Optim_StepSize *ss;

        Optim_KwArg *cur_states;
        unsigned int *cur_iter;
        unsigned int *cur_time;
        PRECISION *cur_err;
        PRECISION *cur_diff;

    public:
        Optim_Solver(Optim_Paras *pa, Optim_Iter *it, Optim_StepSize *s) : paras(pa), iter(it), ss(s)
        {
            if (it)
            {
                cur_states = iter->get_current_states();
                cur_iter = iter->get_cur_iter();
                cur_time = iter->get_cur_time();
                cur_err = iter->get_cur_err();
                cur_diff = iter->get_cur_diff();
            }
        };

        bool terminate();

        virtual void solve() = 0;
    };
    class Optim_Update_Solver : public Optim_Solver
    {

    public:
        Optim_Update_Solver(Optim_Paras *pa, Optim_Iter *it) : Optim_Solver(pa, it, nullptr){};

        void solve();
    };
    class Optim_Stepping_Solver : public Optim_Solver
    {

    public:
        Optim_Stepping_Solver(Optim_Paras *pa, Optim_Iter *it, Optim_StepSize *s) : Optim_Solver(pa, it, s){};

        void solve();
        // Basis solve with Direction function D and constant stepsize. If D is not provided, use gradient function.

        void one_step();
    };
    class Optim_Gauss_Seidel_Solver : public Optim_Solver
    {
    private:
        typedef bool (*Loop)(Optim_KwArg &);
        Loop Loop_Check;
        unsigned *inner_iter;

    public:
        Optim_Gauss_Seidel_Solver(Optim_Paras *pa, Optim_Iter *it, Optim_StepSize *s, Loop loop_check, unsigned *loop_iter)
            : Optim_Solver(pa, it, s)
        {
            // Optim_Iter keeps the inner iterator
            // Solver keeps the inner iterator in inner_iter and the outer iterator in cur_iter.
            inner_iter = cur_iter;
            cur_iter = loop_iter;
            Loop_Check = loop_check;
        }

        void solve();
    };

    class Optim_Trial_Stepping_Solver : public Optim_Stepping_Solver
    {
        typedef bool (*Accept)(Optim_KwArg &kws, unsigned *ind);

    protected:
        Accept A;
        bool *accepted;
        unsigned ind_A[6] = {0, 2, 4, 5, 10, 11}; // cur_iter, cur_err, x, fx, y, fy

    public:
        Optim_Trial_Stepping_Solver(Optim_Paras *pa, Optim_Iter *it, Optim_StepSize *s,
                                    Accept acc, bool *acc_bool) : Optim_Stepping_Solver(pa, it, s), A(acc), accepted(acc_bool){};

        void solve();
        // Basis solve with Direction function D and constant stepsize. If D is not provided, use gradient function.
    };

    // class Optim_Trial_Solver : public Optim_Solver
    // {
    //     typedef bool (*Accept)(KwArg &);
    //     typedef void (*Direction)(KwArg &);

    // protected:
    //     Accept *A;

    // public:
    //     Optim_Trial_Solver(Optim_Paras *pa, Optim_Prob *pr, Optim_Iter *it, Accept *ac) : Optim_Solver(pa, pr, it), A(ac){};

    //     bool solve();
    // };

    void nullOptimFunc(Optim_KwArg &, unsigned *);
}