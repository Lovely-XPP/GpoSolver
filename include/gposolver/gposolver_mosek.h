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

#ifndef GPOSOLVER_MOSEK_H
#define GPOSOLVER_MOSEK_H

#include <stdlib.h>
#include <string>

#include "gposolver_base.h"
#include <mosek.h>

#include <Eigen/SparseCore>

namespace GpoSolver
{
    using namespace std;
    using namespace Eigen;

    template <class T>
    class GpoSolverMosek : public GpoSolverBase<T>
    {
    private:
        using GpoSolverBase<T>::evaluateSdpVariables;
        using GpoSolverBase<T>::problem;
        using GpoSolverBase<T>::cvars_ptr;
        using GpoSolverBase<T>::cvmul;
        using GpoSolverBase<T>::par_verbose;
        using typename GpoSolverBase<T>::SdpStatus;

        typedef SparseMatrix<double> SparseMatrixd;

        class MosekMatrix
        {
        private:
            MSKint32t *is;
            MSKint32t *js;
            MSKrealt *entries;
            int it;
            int reserve;
            int dim;

            void deleteAll(void)
            {
                if (is)
                {
                    delete[] is;
                    is = NULL;
                }
                if (js)
                {
                    delete[] js;
                    js = NULL;
                }
                if (entries)
                {
                    delete[] entries;
                    entries = NULL;
                }
            }

        public:
            MosekMatrix() : is(NULL),
                            js(NULL),
                            entries(NULL),
                            it(0),
                            reserve(0),
                            dim(0) {}

            ~MosekMatrix()
            {
                deleteAll();
            }

            void resize(int s)
            {
                deleteAll();

                is = new MSKint32t[s];
                if (is == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                js = new MSKint32t[s];
                if (js == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                entries = new MSKrealt[s];
                if (entries == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                it = 0;
                reserve = 0;
            }

            void resize(void)
            {
                resize(reserve);
            }

            void setDimension(int d)
            {
                dim = d;
            }

            void addReserve(int s)
            {
                reserve += s;
            }

            void pushEntry(int m, int n, double v)
            {
                if (m < n)
                {
                    int t = m;
                    m = n;
                    n = t;
                }

                is[it] = MSKint32t(m);
                js[it] = MSKint32t(n);
                entries[it] = MSKrealt(v);
                it++;
            }

            void pushBlock(const GpoProblem::conData *data, int ndata, int off, bool mflag = false)
            {
                if (mflag)
                {
                    for (int c = 0; c < ndata; c++)
                        pushEntry(data[c].i + off, data[c].j + off, -data[c].val);
                }
                else
                {
                    for (int c = 0; c < ndata; c++)
                        pushEntry(data[c].i + off, data[c].j + off, data[c].val);
                }
            }

            MSKrealt *getEntries(void)
            {
                return entries;
            }

            int getNumEntries(void)
            {
                return it;
            }

            MSKint32t *getIndicesI(void)
            {
                return is;
            }

            MSKint32t *getIndicesJ(void)
            {
                return js;
            }

            void getSparseMatrix(SparseMatrixd &SpMat)
            {
                vector<Triplet<double>> tripletList;
                tripletList.reserve(it * 2);

                for (int i = 0; i < it; i++)
                {
                    tripletList.push_back(Triplet<double>(is[i], js[i], entries[i]));
                    if (is[i] != js[i])
                        tripletList.push_back(Triplet<double>(js[i], is[i], entries[i]));
                }

                SpMat.resize(dim, dim);
                SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
            }
        };

        MSKenv_t env_msk;
        MSKtask_t ptask_msk;
        VectorXi blk_off;
        vector<MosekMatrix> Abar;
        MSKboundkeye *bkc;
        MSKstreamfunc streamfunc;

        int subToInd(int i, int j, int stride) const
        {
            if (i >= j)
                return j * stride - j * (j - 1) / 2 + (i - j);
            else
                return i * stride - i * (i - 1) / 2 + (j - i);
        }

        void testMSKStatus(const MSKtask_t task_msk, MSKrescodee r) const
        {
            char lastmsg[1000];
            MSKint32t lastmsglen;

            if (r != MSK_RES_OK)
            {
                MSK_getlasterror(task_msk, &r, 1000, &lastmsglen, lastmsg);
                LOG_FATAL << "Mosek Solver Error: " << lastmsg;
            }
        }

        void getZMatrix(const SparseMatrixd &X, const VectorXd &y, vector<MatrixXd> &Z)
        {
            SparseMatrixd Y, W;

            Abar[0].getSparseMatrix(Y);
            for (int i = 1; i <= problem._num_cons; i++)
            {
                Abar[i].getSparseMatrix(W);
                Y -= y(i - 1) * W;
            }

            Z.resize(problem._num_blocks);
            for (int b = 0; b < problem._num_blocks; b++)
                Z[b] = Y.block(blk_off(b), blk_off(b), problem._block_sizes[b], problem._block_sizes[b]);
        }

        void getXMatrix(double *barx, SparseMatrixd &X) const
        {
            vector<Triplet<double>> tripletList;
            int num_entries = 0;

            for (int b = 0; b < problem._num_blocks; b++)
                num_entries += problem._block_sizes[b] * problem._block_sizes[b];
            tripletList.reserve(num_entries);

            for (int b = 0; b < problem._num_blocks; b++)
            {
                for (int i = 0; i < problem._block_sizes[b]; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        double val = barx[subToInd(i + blk_off(b), j + blk_off(b), problem._problem_size)];
                        if (val != 0.0)
                        {
                            tripletList.push_back(Triplet<double>(i + blk_off(b), j + blk_off(b), val));
                            if (i != j)
                                tripletList.push_back(Triplet<double>(j + blk_off(b), i + blk_off(b), val));
                        }
                    }
                }
            }

            X.resize(problem._problem_size, problem._problem_size);
            X.setFromTriplets(tripletList.begin(), tripletList.end());
        }

        SdpStatus solveSdp(double &obj, VectorXd &y, vector<MatrixXd> &Z)
        {
            int prob_size = problem._problem_size;
            int num_cons = problem._num_cons;
            int num_sdp_vars = (prob_size * (prob_size + 1)) / 2;

            VectorXd objective;
            evaluateSdpVariables(objective);

            // Construct Mosek problem
            MSKtask_t task_msk;
            MSKint32t dim[1];
            MSKint64t idx;
            MSKrescodee trmcode;
            MSKsolstae solsta;
            double *barx = NULL;
            double falpha = 1.0;
            dim[0] = prob_size;

            testMSKStatus(ptask_msk, MSK_clonetask(ptask_msk, &task_msk));

            if (streamfunc)
                testMSKStatus(task_msk, MSK_linkfunctotaskstream(task_msk, MSK_STREAM_LOG, NULL, streamfunc));

            testMSKStatus(task_msk, MSK_appendcons(task_msk, num_cons));
            testMSKStatus(task_msk, MSK_appendbarvars(task_msk, 1, dim));

            for (int i = 1; i <= num_cons; i++)
                testMSKStatus(task_msk, MSK_putconbound(task_msk, i - 1, bkc[i - 1], double(objective(i)), double(objective(i))));

            testMSKStatus(task_msk, MSK_appendsparsesymmat(task_msk, dim[0], MSKint64t(Abar[0].getNumEntries()),
                                                           Abar[0].getIndicesI(), Abar[0].getIndicesJ(), Abar[0].getEntries(), &idx));
            testMSKStatus(task_msk, MSK_putbarcj(task_msk, 0, 1, &idx, &falpha));

            for (int i = 1; i <= problem._num_cons; i++)
            {
                testMSKStatus(task_msk, MSK_appendsparsesymmat(task_msk, dim[0], MSKint64t(Abar[i].getNumEntries()),
                                                               Abar[i].getIndicesI(), Abar[i].getIndicesJ(), Abar[i].getEntries(), &idx));
                testMSKStatus(task_msk, MSK_putbaraij(task_msk, i - 1, 0, 1, &idx, &falpha));
            }

            testMSKStatus(task_msk, MSK_optimizetrm(task_msk, &trmcode));

            testMSKStatus(task_msk, MSK_solutionsummary(task_msk, MSK_STREAM_MSG));

            testMSKStatus(task_msk, MSK_getsolsta(task_msk, MSK_SOL_ITR, &solsta));

            if ((solsta == MSK_SOL_STA_OPTIMAL) || (solsta == MSK_SOL_STA_NEAR_OPTIMAL) || (solsta == MSK_SOL_STA_UNKNOWN))
            {
                barx = (double *)MSK_calloctask(task_msk, num_sdp_vars, sizeof(MSKrealt));
                if (barx == NULL)
                    LOG_FATAL << "Cannot allocate memory";

                testMSKStatus(task_msk, MSK_getbarxj(task_msk, MSK_SOL_ITR, 0, barx));

                y.resize(num_cons);
                testMSKStatus(task_msk, MSK_gety(task_msk, MSK_SOL_ITR, y.data()));

                SparseMatrixd X;
                getXMatrix(barx, X);
                getZMatrix(X, y, Z);

                VectorXd b = Map<VectorXd>(objective.data() + 1, num_cons);
                y *= -1.0;
                obj = b.dot(y) + objective(0);

                MSK_freetask(task_msk, barx);
                MSK_deletetask(&task_msk);

                return GpoSolverBase<T>::SDP_SUCCESS;
            }
            else
            {
                MSK_deletetask(&task_msk);
                return GpoSolverBase<T>::SDP_FAILURE;
            }
        }

        static void MSKAPI printToStdout(void *handle, MSKCONST char str[])
        {
            cout << str;
        }

    public:
        GpoSolverMosek() : env_msk(NULL),
                           ptask_msk(NULL),
                           bkc(NULL),
                           streamfunc(NULL)
        {
            GpoSolverBase<T>::init_function();
            // SDP objective (+ objective constant at the 0-th position)
            bkc = new MSKboundkeye[problem._num_cons];
            if (bkc == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 0; i < problem._num_cons; i++)
                bkc[i] = MSK_BK_FX;

            // Abars (+ Cbar at the 0-th position)
            blk_off.resize(problem._num_blocks);
            Abar.resize(problem._num_cons + 1);

            blk_off(0) = 0;
            for (int j = 1; j < problem._num_blocks; j++)
                blk_off(j) = blk_off(j - 1) + problem._block_sizes[j - 1];

            for (int i = 0; i <= problem._num_cons; i++)
            {
                for (int j = 0; j < problem._num_blocks; j++)
                    Abar[i].addReserve(problem.getNumBlockEntries(j, i));
                Abar[i].resize();
                Abar[i].setDimension(problem._problem_size);
            }

            for (int i = 0; i <= problem._num_cons; i++)
            {
                for (int j = 0; j < problem._num_blocks; j++)
                {
                    int ndata;
                    const GpoProblem::conData *data;
                    problem.getBlock(j, i, data, ndata);
                    Abar[i].pushBlock(data, ndata, int(blk_off(j)), false);

                    const list<pair<int, int>> *cvars = problem.getVarsByBlock(j, i);

                    if (cvars == NULL)
                        continue;

                    double *entries = Abar[i].getEntries();
                    int it = Abar[i].getNumEntries();

                    for (list<pair<int, int>>::const_iterator cit = cvars->begin(); cit != cvars->end(); cit++)
                    {
                        cvars_ptr[cit->first] = entries + it - ndata + cit->second;
                        cvmul(cit->first) = 1.0;
                    }
                }
            }

            MSKrescodee r;

            r = MSK_makeenv(&env_msk, NULL);
            if (r != MSK_RES_OK)
                LOG_FATAL << "Mosek Solver Error";

            r = MSK_makeemptytask(env_msk, &ptask_msk);
            if (r != MSK_RES_OK)
                LOG_FATAL << "Mosek Solver Error";
        }

        GpoSolverMosek(const in n) : env_msk(NULL),
                           ptask_msk(NULL),
                           bkc(NULL),
                           streamfunc(NULL)
        {
            GpoSolverBase<T>::init_function(n);
            // SDP objective (+ objective constant at the 0-th position)
            bkc = new MSKboundkeye[problem._num_cons];
            if (bkc == NULL)
                LOG_FATAL << "Cannot allocate memory";

            for (int i = 0; i < problem._num_cons; i++)
                bkc[i] = MSK_BK_FX;

            // Abars (+ Cbar at the 0-th position)
            blk_off.resize(problem._num_blocks);
            Abar.resize(problem._num_cons + 1);

            blk_off(0) = 0;
            for (int j = 1; j < problem._num_blocks; j++)
                blk_off(j) = blk_off(j - 1) + problem._block_sizes[j - 1];

            for (int i = 0; i <= problem._num_cons; i++)
            {
                for (int j = 0; j < problem._num_blocks; j++)
                    Abar[i].addReserve(problem.getNumBlockEntries(j, i));
                Abar[i].resize();
                Abar[i].setDimension(problem._problem_size);
            }

            for (int i = 0; i <= problem._num_cons; i++)
            {
                for (int j = 0; j < problem._num_blocks; j++)
                {
                    int ndata;
                    const GpoProblem::conData *data;
                    problem.getBlock(j, i, data, ndata);
                    Abar[i].pushBlock(data, ndata, int(blk_off(j)), false);

                    const list<pair<int, int>> *cvars = problem.getVarsByBlock(j, i);

                    if (cvars == NULL)
                        continue;

                    double *entries = Abar[i].getEntries();
                    int it = Abar[i].getNumEntries();

                    for (list<pair<int, int>>::const_iterator cit = cvars->begin(); cit != cvars->end(); cit++)
                    {
                        cvars_ptr[cit->first] = entries + it - ndata + cit->second;
                        cvmul(cit->first) = 1.0;
                    }
                }
            }

            MSKrescodee r;

            r = MSK_makeenv(&env_msk, NULL);
            if (r != MSK_RES_OK)
                LOG_FATAL << "Mosek Solver Error";

            r = MSK_makeemptytask(env_msk, &ptask_msk);
            if (r != MSK_RES_OK)
                LOG_FATAL << "Mosek Solver Error";
        }

        virtual ~GpoSolverMosek()
        {
            if (bkc)
            {
                delete[] bkc;
                bkc = NULL;
            }

            MSK_deletetask(&ptask_msk);
            MSK_deleteenv(&env_msk);
        }

        double getParameter(const string &param_name)
        {
            if (param_name.find("MSK_IPAR", 0) == 0)
            {
                MSKint32t ipar;
                testMSKStatus(ptask_msk, MSK_getnaintparam(ptask_msk, param_name.c_str(), &ipar));
                return double(ipar);
            }
            else if (param_name.find("MSK_DPAR", 0) == 0)
            {
                MSKrealt rpar;
                testMSKStatus(ptask_msk, MSK_getnadouparam(ptask_msk, param_name.c_str(), &rpar));
                return double(rpar);
            }
            else if (param_name.find("MSK_SPAR", 0) == 0)
                LOG_FATAL << "Cannot get Mosek string parameters: " << param_name;

            return GpoSolverBase<T>::getParameter(param_name);
        }

        void setParameter(const string &param_name, const double val)
        {
            if (param_name == "verbose")
            {
                if (val == 1)
                    streamfunc = printToStdout;
                else
                    streamfunc = NULL;
                par_verbose = val;
            }
            else if (param_name.find("MSK_IPAR", 0) == 0)
                testMSKStatus(ptask_msk, MSK_putnaintparam(ptask_msk, param_name.c_str(), MSKint32t(val)));
            else if (param_name.find("MSK_DPAR", 0) == 0)
                testMSKStatus(ptask_msk, MSK_putnadouparam(ptask_msk, param_name.c_str(), MSKrealt(val)));
            else if (param_name.find("MSK_SPAR", 0) == 0)
                LOG_FATAL << "Cannot set Mosek string parameters: " << param_name;
            else
                GpoSolverBase<T>::setParameter(param_name, val);
        }

        void setMSKstreamfunc(MSKstreamfunc func)
        {
            streamfunc = func;
        }

    }; // GPOSolverMosek

} // GPOSolver namespace

#endif // GPOSOLVER_MOSEK_H
