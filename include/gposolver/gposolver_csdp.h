/*
    GpoSolver, Library for Global Polynomial Optimization
    Copyright (C) 2014-2015 Jan Heller, <hellej1@cmp.felk.cvut.cz>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GPOSOLVER_CSDP_H
#define GPOSOLVER_CSDP_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#ifndef NOSHORTS
// CSDP library must be compiled with -DNOSHORTS
#define NOSHORTS
#endif

extern "C"
{
#include <declarations.h>
}

#include "gposolver_base.h"

namespace GpoSolver
{
    using namespace std;
    using namespace Eigen;

    template <class T>
    class GpoSolverCsdp : public GpoSolverBase<T>
    {
    private:
        using GpoSolverBase<T>::evaluateSdpVariables;
        using GpoSolverBase<T>::problem;
        using GpoSolverBase<T>::cvars_ptr;
        using GpoSolverBase<T>::cvmul;
        using GpoSolverBase<T>::par_verbose;
        using typename GpoSolverBase<T>::SdpStatus;

        struct constraintmatrix *cons_csdp;
        struct blockmatrix C_csdp, X_csdp, Z_csdp, T_csdp;
        double *y_csdp;
        double pobj_csdp, dobj_csdp;
        double csdp_flag;

        double par_axtol;
        double par_atytol;
        double par_objtol;
        double par_pinftol;
        double par_dinftol;
        int par_maxiter;
        double par_minstepfrac;
        double par_maxstepfrac;
        double par_minstepp;
        double par_minstepd;
        int par_usexzgap;
        int par_tweakgap;
        int par_affine;
        int par_printlevel;
        double par_perturbobj;
        int par_fastmode;

    protected:
        void getZMatrix(struct blockmatrix &M, vector<MatrixXd> &Z)
        {
            int num_blocks = problem._num_blocks;

            Z.resize(num_blocks);
            for (int i = 0; i < num_blocks; i++)
            {
                int bsize = problem._block_sizes[i];
                Z[i] = Map<MatrixXd>(M.blocks[i + 1].data.mat, bsize, bsize);
            }
        }

        void freeBlockMatrix(struct blockmatrix &M)
        {
            for (int i = 1; i <= M.nblocks; i++)
            {
                if (M.blocks[i].blockcategory == DIAG)
                    free(M.blocks[i].data.vec);
                else if (M.blocks[i].blockcategory == MATRIX)
                    free(M.blocks[i].data.mat);
            }

            free(M.blocks);
        }

        void freeConstraints()
        {
            struct sparseblock *ptr;
            struct sparseblock *oldptr;

            if (cons_csdp != NULL)
            {
                for (int i = 1; i <= problem._num_cons; i++)
                {
                    ptr = cons_csdp[i].blocks;
                    while (ptr != NULL)
                    {
                        free(ptr->entries);
                        free(ptr->iindices);
                        free(ptr->jindices);
                        oldptr = ptr;
                        ptr = ptr->next;
                        free(oldptr);
                    }
                }

                free(cons_csdp);
            }
        }

        /*
         * Then free space for the block records.
         */

        SdpStatus solveSdp(double &obj, VectorXd &y, vector<MatrixXd> &Z)
        {
            int prob_size = problem._problem_size;
            int num_cons = problem._num_cons;

            VectorXd objective;
            evaluateSdpVariables(objective);

            ofstream pfile("param.csdp");

            if (!pfile.is_open())
                LOG_FATAL << "Cannot open file 'param.csdp' in the current directory for writing";

            pfile << "axtol=" << par_axtol << endl;
            pfile << "atytol=" << par_atytol << endl;
            pfile << "objtol=" << par_objtol << endl;
            pfile << "pinftol=" << par_pinftol << endl;
            pfile << "dinftol=" << par_dinftol << endl;
            pfile << "maxiter=" << par_maxiter << endl;
            pfile << "minstepfrac=" << par_minstepfrac << endl;
            pfile << "maxstepfrac=" << par_maxstepfrac << endl;
            pfile << "minstepp=" << par_minstepp << endl;
            pfile << "minstepd=" << par_minstepd << endl;
            pfile << "usexzgap=" << par_usexzgap << endl;
            pfile << "tweakgap=" << par_tweakgap << endl;
            pfile << "affine=" << par_affine << endl;
            pfile << "printlevel=" << par_printlevel << endl;
            pfile << "perturbobj=" << par_perturbobj << endl;
            pfile << "fastmode=" << par_fastmode << endl;

            pfile.close();

            // Solve SDP
            initsoln(prob_size, num_cons, C_csdp, objective.data(), cons_csdp, &X_csdp, &y_csdp, &Z_csdp);
            int ret = easy_sdp(prob_size, num_cons, C_csdp, objective.data(), cons_csdp, 0.0,
                               &X_csdp, &y_csdp, &Z_csdp, &pobj_csdp, &dobj_csdp);
            csdp_flag = true;
            remove("param.csdp");

            if (ret > 3)
                return GpoSolverBase<T>::SDP_FAILURE;

            VectorXd b = Map<VectorXd>(objective.data() + 1, num_cons);
            y = Map<VectorXd>(y_csdp + 1, num_cons);
            obj = (b.dot(y) + objective(0));
            getZMatrix(Z_csdp, Z);

            free_prob(problem._problem_size, problem._num_cons, T_csdp, NULL, NULL, X_csdp, y_csdp, Z_csdp);

            return GpoSolverBase<T>::SDP_SUCCESS;
        }

    public:
        GpoSolverCsdp()
        {
            GpoSolverBase<T>::init_function();
            cons_csdp = NULL;
            y_csdp = NULL;
            csdp_flag = false;

            pobj_csdp = 0;
            dobj_csdp = 0;

            // Parameters
            par_axtol = 1.0e-8;
            par_atytol = 1.0e-8;
            par_objtol = 1.0e-8;
            par_pinftol = 1.0e8;
            par_dinftol = 1.0e8;
            par_maxiter = 100;
            par_minstepfrac = 0.90;
            par_maxstepfrac = 0.97;
            par_minstepp = 1.0e-8;
            par_minstepd = 1.0e-8;
            par_usexzgap = 1;
            par_tweakgap = 0;
            par_affine = 0;
            par_printlevel = 0;
            par_perturbobj = 1;
            par_fastmode = 0;

            // matrices X,Z,T
            X_csdp.nblocks = 0;
            Z_csdp.nblocks = 0;
            T_csdp.nblocks = 0;
            T_csdp.blocks = NULL;

            // matrix C
            C_csdp.nblocks = problem._num_blocks;
            C_csdp.blocks = (struct blockrec *)calloc((problem._num_blocks + 1), sizeof(struct blockrec));
            if (C_csdp.blocks == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 1; i <= problem._num_blocks; i++)
            {
                int bsize = problem._block_sizes[i - 1];
                C_csdp.blocks[i].blocksize = bsize;
                C_csdp.blocks[i].blockcategory = MATRIX;
                C_csdp.blocks[i].data.mat = (double *)calloc(bsize * bsize, sizeof(double));
                if (C_csdp.blocks[i].data.mat == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                int ndata;
                const GpoProblem::conData *data;
                problem.getBlock(i - 1, 0, data, ndata);

                for (int j = 0; j < ndata; j++)
                {
                    int m = data[j].i + 1;
                    int n = data[j].j + 1;
                    C_csdp.blocks[i].data.mat[ijtok(m, n, bsize)] = -data[j].val;
                }

                const list<pair<int, int>> *cvars = problem.getVarsByBlock(i - 1, 0);

                if (cvars == NULL)
                    continue;

                for (list<pair<int, int>>::const_iterator it = cvars->begin(); it != cvars->end(); it++)
                {
                    int j = it->second;
                    int k = it->first;
                    int m = data[j].i + 1;
                    int n = data[j].j + 1;
                    cvars_ptr[k] = &(C_csdp.blocks[i].data.mat[ijtok(m, n, bsize)]);
                    cvmul(k) = -1.0;
                }
            }

            // Constraints
            struct sparseblock **cons_blocks = (struct sparseblock **)calloc(problem._num_cons + 1, sizeof(struct sparseblock *));
            if (cons_blocks == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cons_csdp = (struct constraintmatrix *)malloc((problem._num_cons + 1) * sizeof(struct constraintmatrix));
            if (cons_csdp == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 1; i <= problem._num_cons; i++)
            {
                const set<int> *cons_idx = problem.getBlocksByConstraint(i);

                if (!cons_idx)
                    continue;

                struct sparseblock **sb = &(cons_csdp[i].blocks);
                for (set<int>::const_iterator it = cons_idx->begin(); it != cons_idx->end(); it++)
                {
                    int ndata;
                    const GpoProblem::conData *data;
                    int j = *it;
                    problem.getBlock(j, i, data, ndata);

                    *sb = (struct sparseblock *)calloc(1, sizeof(struct sparseblock));
                    if (*sb == NULL)
                        LOG_FATAL << "Cannot allocate memory";

                    if (cons_blocks[i] == NULL)
                        cons_blocks[i] = *sb;
                    else
                    {
                        cons_blocks[i]->nextbyblock = *sb;
                        cons_blocks[i] = *sb;
                    }

                    (*sb)->constraintnum = i;
                    (*sb)->blocknum = j + 1;
                    (*sb)->numentries = ndata;
                    (*sb)->blocksize = problem._block_sizes[j];
                    (*sb)->issparse = ((ndata > 0.25 * (*sb)->blocksize) && ndata > 15) ? 0 : 1;
                    (*sb)->iindices = (int *)malloc((ndata + 1) * sizeof(int));
                    (*sb)->jindices = (int *)malloc((ndata + 1) * sizeof(int));
                    (*sb)->entries = (double *)malloc((ndata + 1) * sizeof(double));
                    if (((*sb)->iindices == NULL) || ((*sb)->jindices == NULL) || ((*sb)->entries == NULL))
                        LOG_FATAL << "Cannot allocate memory";

                    for (int k = 1; k <= ndata; k++)
                    {
                        (*sb)->iindices[k] = (int)data[k - 1].i + 1;
                        (*sb)->jindices[k] = (int)data[k - 1].j + 1;
                        (*sb)->entries[k] = (int)data[k - 1].val;
                    }

                    const list<pair<int, int>> *cvars = problem.getVarsByBlock(j, i);

                    if (!cvars)
                    {
                        sb = &((*sb)->next);
                        continue;
                    }

                    for (list<pair<int, int>>::const_iterator cit = cvars->begin(); cit != cvars->end(); cit++)
                        cvars_ptr[cit->first] = &((*sb)->entries[cit->second + 1]);

                    sb = &((*sb)->next);
                }
            }

            free(cons_blocks);
        }

        GpoSolverCsdp(const int n)
        {
            GpoSolverBase<T>::init_function(n);
            cons_csdp = NULL;
            y_csdp = NULL;
            csdp_flag = false;

            pobj_csdp = 0;
            dobj_csdp = 0;

            // Parameters
            par_axtol = 1.0e-8;
            par_atytol = 1.0e-8;
            par_objtol = 1.0e-8;
            par_pinftol = 1.0e8;
            par_dinftol = 1.0e8;
            par_maxiter = 100;
            par_minstepfrac = 0.90;
            par_maxstepfrac = 0.97;
            par_minstepp = 1.0e-8;
            par_minstepd = 1.0e-8;
            par_usexzgap = 1;
            par_tweakgap = 0;
            par_affine = 0;
            par_printlevel = 0;
            par_perturbobj = 1;
            par_fastmode = 0;

            // matrices X,Z,T
            X_csdp.nblocks = 0;
            Z_csdp.nblocks = 0;
            T_csdp.nblocks = 0;
            T_csdp.blocks = NULL;

            // matrix C
            C_csdp.nblocks = problem._num_blocks;
            C_csdp.blocks = (struct blockrec *)calloc((problem._num_blocks + 1), sizeof(struct blockrec));
            if (C_csdp.blocks == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 1; i <= problem._num_blocks; i++)
            {
                int bsize = problem._block_sizes[i - 1];
                C_csdp.blocks[i].blocksize = bsize;
                C_csdp.blocks[i].blockcategory = MATRIX;
                C_csdp.blocks[i].data.mat = (double *)calloc(bsize * bsize, sizeof(double));
                if (C_csdp.blocks[i].data.mat == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                int ndata;
                const GpoProblem::conData *data;
                problem.getBlock(i - 1, 0, data, ndata);

                for (int j = 0; j < ndata; j++)
                {
                    int m = data[j].i + 1;
                    int n = data[j].j + 1;
                    C_csdp.blocks[i].data.mat[ijtok(m, n, bsize)] = -data[j].val;
                }

                const list<pair<int, int>> *cvars = problem.getVarsByBlock(i - 1, 0);

                if (cvars == NULL)
                    continue;

                for (list<pair<int, int>>::const_iterator it = cvars->begin(); it != cvars->end(); it++)
                {
                    int j = it->second;
                    int k = it->first;
                    int m = data[j].i + 1;
                    int n = data[j].j + 1;
                    cvars_ptr[k] = &(C_csdp.blocks[i].data.mat[ijtok(m, n, bsize)]);
                    cvmul(k) = -1.0;
                }
            }

            // Constraints
            struct sparseblock **cons_blocks = (struct sparseblock **)calloc(problem._num_cons + 1, sizeof(struct sparseblock *));
            if (cons_blocks == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cons_csdp = (struct constraintmatrix *)malloc((problem._num_cons + 1) * sizeof(struct constraintmatrix));
            if (cons_csdp == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 1; i <= problem._num_cons; i++)
            {
                const set<int> *cons_idx = problem.getBlocksByConstraint(i);

                if (!cons_idx)
                    continue;

                struct sparseblock **sb = &(cons_csdp[i].blocks);
                for (set<int>::const_iterator it = cons_idx->begin(); it != cons_idx->end(); it++)
                {
                    int ndata;
                    const GpoProblem::conData *data;
                    int j = *it;
                    problem.getBlock(j, i, data, ndata);

                    *sb = (struct sparseblock *)calloc(1, sizeof(struct sparseblock));
                    if (*sb == NULL)
                        LOG_FATAL << "Cannot allocate memory";

                    if (cons_blocks[i] == NULL)
                        cons_blocks[i] = *sb;
                    else
                    {
                        cons_blocks[i]->nextbyblock = *sb;
                        cons_blocks[i] = *sb;
                    }

                    (*sb)->constraintnum = i;
                    (*sb)->blocknum = j + 1;
                    (*sb)->numentries = ndata;
                    (*sb)->blocksize = problem._block_sizes[j];
                    (*sb)->issparse = ((ndata > 0.25 * (*sb)->blocksize) && ndata > 15) ? 0 : 1;
                    (*sb)->iindices = (int *)malloc((ndata + 1) * sizeof(int));
                    (*sb)->jindices = (int *)malloc((ndata + 1) * sizeof(int));
                    (*sb)->entries = (double *)malloc((ndata + 1) * sizeof(double));
                    if (((*sb)->iindices == NULL) || ((*sb)->jindices == NULL) || ((*sb)->entries == NULL))
                        LOG_FATAL << "Cannot allocate memory";

                    for (int k = 1; k <= ndata; k++)
                    {
                        (*sb)->iindices[k] = (int)data[k - 1].i + 1;
                        (*sb)->jindices[k] = (int)data[k - 1].j + 1;
                        (*sb)->entries[k] = (int)data[k - 1].val;
                    }

                    const list<pair<int, int>> *cvars = problem.getVarsByBlock(j, i);

                    if (!cvars)
                    {
                        sb = &((*sb)->next);
                        continue;
                    }

                    for (list<pair<int, int>>::const_iterator cit = cvars->begin(); cit != cvars->end(); cit++)
                        cvars_ptr[cit->first] = &((*sb)->entries[cit->second + 1]);

                    sb = &((*sb)->next);
                }
            }

            free(cons_blocks);
        }

        virtual ~GpoSolverCsdp()
        {
            freeBlockMatrix(C_csdp);
            freeConstraints();
        }

        double getParameter(const string &param_name)
        {
            if (param_name == "axtol")
                return par_axtol;
            else if (param_name == "atytol")
                return par_atytol;
            else if (param_name == "objtol")
                return par_objtol;
            else if (param_name == "pinftol")
                return par_pinftol;
            else if (param_name == "dinftol")
                return par_dinftol;
            else if (param_name == "maxiter")
                return par_maxiter;
            else if (param_name == "minstepfrac")
                return par_minstepfrac;
            else if (param_name == "maxstepfrac")
                return par_maxstepfrac;
            else if (param_name == "minstepp")
                return par_minstepp;
            else if (param_name == "minstepd")
                return par_minstepd;
            else if (param_name == "usexzgap")
                return par_usexzgap;
            else if (param_name == "tweakgap")
                return par_tweakgap;
            else if (param_name == "affine")
                return par_affine;
            else if (param_name == "printlevel")
                return par_printlevel;
            else if (param_name == "perturbobj")
                return par_perturbobj;
            else if (param_name == "fastmode")
                return par_fastmode;

            return GpoSolverBase<T>::getParameter(param_name);
        }

        void setParameter(const string &param_name, const double val)
        {
            if (param_name == "verbose")
            {
                par_printlevel = (val) ? 1 : 0;
                par_verbose = val;
            }
            else if (param_name == "axtol")
                par_axtol = val;
            else if (param_name == "atytol")
                par_atytol = val;
            else if (param_name == "objtol")
                par_objtol = val;
            else if (param_name == "pinftol")
                par_pinftol = val;
            else if (param_name == "dinftol")
                par_dinftol = val;
            else if (param_name == "maxiter")
                par_maxiter = (int)val;
            else if (param_name == "minstepfrac")
                par_minstepfrac = val;
            else if (param_name == "maxstepfrac")
                par_maxstepfrac = val;
            else if (param_name == "minstepp")
                par_minstepp = val;
            else if (param_name == "minstepd")
                par_minstepd = val;
            else if (param_name == "usexzgap")
                par_usexzgap = (val == 0.0) ? 0 : 1;
            else if (param_name == "tweakgap")
                par_tweakgap = (val == 0.0) ? 0 : 1;
            else if (param_name == "affine")
                par_affine = (val == 0.0) ? 0 : 1;
            else if (param_name == "printlevel")
                par_printlevel = (int)val;
            else if (param_name == "perturbobj")
                par_perturbobj = val;
            else if (param_name == "fastmode")
                par_fastmode = (val == 0.0) ? 0 : 1;
            else
                GpoSolverBase<T>::setParameter(param_name, val);
        }

    }; // GpoSolverCsdp

} // GpoSolver namespace

#endif // GPOSOLVER_CSDP_H
