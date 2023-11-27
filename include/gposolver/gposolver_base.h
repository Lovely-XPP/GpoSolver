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

#ifndef GPOSOLVER_BASE_H
#define GPOSOLVER_BASE_H

#include <iostream>
#include <string>
#include <fstream>
#include <cfloat>

#include <set>
#include <map>
#include <list>
#include <vector>

#include <Eigen/Core>

#include "logging.h"
#include "gpoproblem.h"

namespace GpoSolver
{
    using namespace std;
    using namespace Eigen;

    // Problem status
    typedef enum
    {
        SUCCESS_GLOPT = 0,
        SUCCESS = 1,
        UNBOUNDED_SOLUTION = 2,
        INFEASIBLE_SDP = 3,
        INCOMPLETE_BASIS = 4,
        FAILURE = 5,
    } GpoSolverStatus;

    // Problem solutions
    typedef list<vector<double>> GpoSolverSolutions;

    // GpoSolverBase
    template <class T>
    class GpoSolverBase
    {
    public:
    protected:
        typedef enum
        {
            SDP_SUCCESS = 0,
            SDP_FAILURE = 1,
        } SdpStatus;

        // Solver definition and data =============================================

        T problem_definition;
        GpoProblem &problem;
        int num_res;
        const double *rvals;
        const double *pvals;
        double **rvars_ptr;
        double **cvars_ptr;
        VectorXd cvmul;
        VectorXd obj_res;

        // Parameters =============================================================

        double par_sdptol;
        double par_maxnorm;
        double par_restol;
        double par_ranktol;
        double par_pivtol;
        int par_verbose;

        //  SDP problem ===========================================================

        virtual SdpStatus solveSdp(double &, VectorXd &, vector<MatrixXd> &) = 0;

        void evaluateSdpVariables(VectorXd &objective)
        {
            int num_rvars = problem._num_rpars;
            int num_cons = problem._num_cons;

            // Evaluate Objective
            objective.resize(num_cons + 1);
            objective.setZero();

            if (num_res == 0)
            {
                objective = obj_res;
            }
            else
            {
                for (int i = 0; i < num_res; i++)
                {
                    problem.evalSdpObjectiveVariables(rvals + i * num_rvars, pvals, rvars_ptr);
                    objective += obj_res;
                }
            }

            // Evaluate Constraints
            problem.evalSdpConstraintsVariables(pvals, cvmul.data(), cvars_ptr);
        }

        // Solution validation and extraction =====================================

        int nChoosek(int n, int k) const
        {
            if (k > n)
                return 0;

            if (k * 2 > n)
                k = n - k;

            if (k == 0)
                return 1;

            int r = n;

            for (int i = 2; i <= k; i++)
            {
                r = r * (n - i + 1);
                r = r / i;
            }

            return r;
        }

        double matrixCond(MatrixXd &M) const
        {
            JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
            const VectorXd &s = svd.singularValues();
            int s_min = s.minCoeff();

            if (s_min == 0.0)
                return DBL_MAX;
            else
                return s.maxCoeff() / s_min;
        }

        void givensRotation(const double a, const double b, Matrix2d &R) const
        {
            double t, c, s;

            if (b == 0)
            {
                c = 1;
                s = 0;
            }
            else
            {
                if (std::abs(b) > std::abs(a))
                {
                    t = -a / b;
                    s = 1 / std::sqrt(1 + t * t);
                    c = s * t;
                }
                else
                {
                    t = -b / a;
                    c = 1 / std::sqrt(1 + t * t);
                    s = c * t;
                }
            }

            R << c, s, -s, c;
        }

        bool orderedSchur(const MatrixXd &M, MatrixXd &U, MatrixXd &K) const
        {
            ComplexSchur<MatrixXd> schur(M.cols());

            schur.compute(M, true);

            if (schur.matrixT().imag().diagonal().norm() > 1e-8)
                return false;

            U = schur.matrixU().real();
            K = schur.matrixT().real();

            int n = K.rows();
            bool ordered = false;
            while (!ordered)
            {
                ordered = true;
                for (int k = 0; k < n - 1; k++)
                {
                    if (K(k, k) > K(k + 1, k + 1))
                    {
                        ordered = false;
                        Matrix2d R;
                        MatrixXd Q;
                        givensRotation(double(K(k, k + 1)), double(K(k + 1, k + 1) - K(k, k)), R);

                        Q = R.transpose() * K.block(k, k, 2, n - k);
                        K.block(k, k, 2, n - k) = Q;

                        Q = K.block(0, k, k + 2, 2) * R;
                        K.block(0, k, k + 2, 2) = Q;

                        Q = U.block(0, k, n, 2) * R;
                        U.block(0, k, n, 2) = Q;
                    }
                }
            }

            return true;
        }

        void colEchelonForm(MatrixXd &M, list<int> &b, double pivtol) const
        {
            int n = M.rows();
            int m = M.cols();
            int i = 0, j = 0, k = 0;
            int col = 0;
            double p, tp;

            b.clear();

            while ((i < m) && (j < n))
            {
                p = DBL_MIN;
                col = i;
                for (k = i; k < m; k++)
                {
                    tp = std::abs(M(j, k));
                    if (tp > p)
                    {
                        p = tp;
                        col = k;
                    }
                }

                if (p < pivtol)
                {
                    M.block(j, i, 1, m - i).setZero();
                    j++;
                }
                else
                {
                    b.push_back(j);

                    if (col != i)
                        M.block(j, i, n - j, 1).swap(M.block(j, col, n - j, 1));

                    M.block(j + 1, i, n - j - 1, 1) /= M(j, i);
                    M(j, i) = 1.0;

                    for (k = 0; k < m; k++)
                    {
                        if (k == i)
                            continue;

                        M.block(j, k, n - j, 1) -= M(j, k) * M.block(j, i, n - j, 1);
                    }

                    i++;
                    j++;
                }
            }
        }

        int getPowers(MatrixXi &M, int row, int col, int deg) const
        {
            int roff = 0;

            if (col + 1 == M.cols())
            {
                M(row, col) = deg;
                return 1;
            }
            else
            {
                if (deg == 0)
                    return 1;

                for (int i = deg; i >= 0; i--)
                {
                    int r = getPowers(M, row + roff, col + 1, deg - i);
                    M.block(row + roff, col, r, 1).fill(i);
                    roff += r;
                }
            }

            return roff;
        }

        void getPolyBasis(MatrixXi &B, int num_vars, int deg) const
        {
            int num_rows = nChoosek(num_vars + deg, deg);

            B.resize(num_rows, num_vars);
            B.setZero();

            int row = 0;
            for (int i = 0; i <= deg; i++)
                row += getPowers(B, row, 0, i);
        }

        void getBinomialCoeffs(MatrixXi &C, int num_vars, int deg) const
        {
            C.resize(deg + 2, num_vars);

            for (int i = 0; i < num_vars; i++)
            {
                int sum = 0;
                for (int j = 0; j < deg + 2; j++)
                {
                    if (i == 0)
                    {
                        C(j, i) = j;
                    }
                    else
                    {
                        sum += C(j, i - 1);
                        C(j, i) = sum;
                    }
                }
            }
        }

        double getSdpResidual(const vector<MatrixXd> &Z) const
        {
            double res = DBL_MAX;

            for (unsigned int i = 0; i < Z.size(); i++)
                res = std::min(res, Z[i].eigenvalues().real().minCoeff());

            return res;
        }

        int getRank(const MatrixXd &M, MatrixXd *V = NULL) const
        {
            JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
            const VectorXd &s = svd.singularValues();

            int rank = s.size();
            for (int i = 1; i < s.size(); i++)
            {
                if ((s(i) / (DBL_MIN + s(i - 1))) < par_ranktol)
                {
                    rank = i;
                    break;
                }
            }

            if (V != NULL)
            {
                const MatrixXd &U = svd.matrixU();
                *V = U.block(0, 0, U.rows(), rank) * s.block(0, 0, rank, 1).cwiseSqrt().asDiagonal();
            }

            return rank;
        }

        bool isSolutionFeasible(const VectorXd &sol, const double sdp_obj) const
        {
            int num_rvars = problem._num_rpars;
            double poly_obj = 0;

            if (num_res == 0)
                poly_obj = problem.evalPolyObjective(NULL, pvals, sol.data());
            else
                for (int i = 0; i < num_res; i++)
                    poly_obj += problem.evalPolyObjective(rvals + i * num_rvars, pvals, sol.data());

            if (std::abs(poly_obj - sdp_obj) > par_restol)
                return false;

            VectorXd poly_cons(problem._num_pcons);
            problem.evalPolyConstraints(sol.data(), pvals, poly_cons.data());

            for (int i = 0; i < problem._num_pcons; i++)
            {
                if (poly_cons(i) < -par_restol)
                    return false;
            }

            return true;
        }

        bool doRankTest(const MatrixXd &M, MatrixXd *V) const
        {
            int rank, prev_rank = 1, rank_diff = 0;
            int rank_shift = problem._rank_shift;
            int relax_order = problem._relax_order;
            int num_vars = problem._num_pvars;
            bool glopt = false;

            for (int rorder = 1; rorder <= relax_order; rorder++)
            {
                int dim = nChoosek(num_vars + rorder, num_vars);
                MatrixXd N = M.block(0, 0, dim, dim);

                if (rorder < relax_order)
                    rank = getRank(N);
                else
                    rank = getRank(N, V);

                if (prev_rank >= rank)
                    rank_diff++;
                else
                    rank_diff = 0;

                if (rank_diff >= rank_shift)
                {
                    glopt = true;
                    break;
                }

                prev_rank = rank;
            }

            return glopt;
        }

        bool extractSolutions(MatrixXd &V, const double sdp_obj, GpoSolverSolutions &sols) const
        {
            list<int> vbasis;
            MatrixXi pbasis, vpbasis, bcoeffs;
            int num_vars = problem._num_pvars;
            int degree = problem._relax_order * 2;
            vector<MatrixXd> N;
            N.resize(num_vars);

            getPolyBasis(pbasis, num_vars, degree);
            getBinomialCoeffs(bcoeffs, num_vars, degree);

            // Stabilizing hack as seen in YALMIP
            if (par_pivtol > 0)
                colEchelonForm(V, vbasis, par_pivtol);
            else
            {
                MatrixXd W = V;
                double tol = 1e-10;
                colEchelonForm(W, vbasis, tol);
                double cW = matrixCond(W);
                double cV = matrixCond(V);

                while ((tol < 1) && ((cW / cV) > 1e4))
                {
                    tol *= 5;
                    W = V;
                    colEchelonForm(W, vbasis, tol);
                    cW = matrixCond(W);
                }
                V = W;
            }

            int c = 0;
            bool fail_flag = false;
            int num_sols = vbasis.size();
            vpbasis.resize(num_sols, num_vars);
            for (list<int>::iterator it = vbasis.begin(); it != vbasis.end(); it++)
                vpbasis.row(c++) = pbasis.row(*it);

            MatrixXd M(num_sols, V.cols());
            M.setZero();
            VectorXd cfs(num_vars);
            cfs.setRandom();
            cfs = (cfs + cfs.Ones(num_vars)) / 2.0;
            double scfs = cfs.sum();
            cfs /= scfs;

            for (int i = 0; (i < num_vars) && (!fail_flag); i++)
            {
                MatrixXi vpb = vpbasis;
                VectorXi midx = vpb.Ones(num_sols, 1);
                vpb.col(i) += midx;

                for (int j = 0; j < num_vars; j++)
                {
                    VectorXi idx = vpb.rightCols(num_vars - j).rowwise().sum();
                    for (int k = 0; k < num_sols; k++)
                        midx(k) += bcoeffs(idx(k), num_vars - j - 1);
                }

                N[i].resize(num_sols, V.cols());
                for (int j = 0; j < num_sols; j++)
                {
                    if (int(midx(j)) <= V.rows())
                        N[i].row(j) = V.row(int(midx(j)) - 1);
                    else
                        fail_flag = true;
                }

                M += cfs(i) * N[i];
            }

            if (fail_flag)
                return false;

            MatrixXd U, K;
            if (!orderedSchur(M, U, K))
                return false;

            for (int i = 0; i < num_sols; i++)
            {
                VectorXd sol(num_vars);

                for (int j = 0; j < num_vars; j++)
                    sol(j) = U.col(i).transpose() * N[j] * U.col(i);

                if (isSolutionFeasible(sol, sdp_obj))
                    addSolution(sol, sols);
            }

            return true;
        }

        void addSolution(const VectorXd &sol, GpoSolverSolutions &sols) const
        {
            sols.push_back(vector<double>());
            vector<double> &v = sols.back();
            v.resize(sol.rows());

            for (int i = 0; i < sol.rows(); i++)
                v[i] = double(sol(i));
        }

        // Public methods =========================================================

    public:
        GpoSolverBase()
            : problem(problem_definition), num_res(0),
              rvals(NULL),
              pvals(NULL),
              rvars_ptr(NULL),
              cvars_ptr(NULL),
              par_sdptol(1e-03),
              par_maxnorm(1e+06),
              par_restol(1e-02),
              par_ranktol(1e-03),
              par_pivtol(-1),
              par_verbose(0)
        {}

        void init_function(const int n)
        {
            problem.init(n);

            rvars_ptr = new double *[problem._num_obj_vars];
            if (rvars_ptr == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cvars_ptr = new double *[problem._num_con_vars];
            if (cvars_ptr == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cvmul.resize(problem._num_con_vars);
            cvmul.setOnes();

            obj_res.resize(problem._num_cons + 1);
            obj_res.setZero();

            const GpoProblem::objData *odata = problem._obj_data;

            for (int i = 0; i < problem._num_obj_data; i++)
                obj_res(odata[i].mon) = odata[i].val;

            for (int i = 0; i < problem._num_obj_vars; i++)
                rvars_ptr[i] = obj_res.data() + problem._obj_vars[i].mon;
        }

        void init_function(void)
        {
            problem.init();

            rvars_ptr = new double *[problem._num_obj_vars];
            if (rvars_ptr == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cvars_ptr = new double *[problem._num_con_vars];
            if (cvars_ptr == NULL)
                LOG_FATAL << "Cannot allocate memory";

            cvmul.resize(problem._num_con_vars);
            cvmul.setOnes();

            obj_res.resize(problem._num_cons + 1);
            obj_res.setZero();

            const GpoProblem::objData *odata = problem._obj_data;

            for (int i = 0; i < problem._num_obj_data; i++)
                obj_res(odata[i].mon) = odata[i].val;

            for (int i = 0; i < problem._num_obj_vars; i++)
                rvars_ptr[i] = obj_res.data() + problem._obj_vars[i].mon;
        }

        virtual ~GpoSolverBase()
        {
            if (rvars_ptr)
            {
                delete[] rvars_ptr;
                rvars_ptr = NULL;
            }
            if (cvars_ptr)
            {
                delete[] cvars_ptr;
                cvars_ptr = NULL;
            }
        }

        GpoSolverStatus solve(const int num_res, const double *rvals, const double *pvals, GpoSolverSolutions &sols)
        {
            VectorXd y, sol;
            vector<MatrixXd> Z;
            double sdp_obj;

            // Setup problem instance
            this->num_res = num_res;
            this->rvals = rvals;
            this->pvals = pvals;

            if ((num_res != 0) && (rvals == NULL))
                LOG_FATAL << "(num_res != 0) && (ovals == NULL)";

            if ((rvals == NULL) && (problem._num_rpars != 0))
                LOG_FATAL << "Objective function parameters expected: ovals == NULL";

            if ((pvals == NULL) && (problem._num_ppars != 0))
                LOG_FATAL << "Constraints parameters expected: pvals == NULL";

            sols.clear();

            // Solve SDP
            SdpStatus ret = solveSdp(sdp_obj, y, Z);

            // Test solution

            if (ret != SDP_SUCCESS)
                return INFEASIBLE_SDP;

            double sdp_res = getSdpResidual(Z);

            if (sdp_res < -par_sdptol)
                return INFEASIBLE_SDP;

            if (y.norm() > par_maxnorm)
                return UNBOUNDED_SOLUTION;

            // Check global optimality of one solution
            sol = y.block(0, 0, problem._num_pvars, 1);
            bool feas = isSolutionFeasible(sol, sdp_obj);

            if (feas)
            {
                addSolution(sol, sols);
                return SUCCESS_GLOPT;
            }

            // Check global optimality for multiple solutions
            MatrixXd V(0, 0);
            bool glopt_flag = doRankTest(Z[0], &V);

            // Extract solutions
            if (V.cols() == 0)
                getRank(Z[0], &V);

            bool sol_flag = extractSolutions(V, sdp_obj, sols);

            if (!sol_flag)
                return INCOMPLETE_BASIS;

            if (glopt_flag)
                return SUCCESS_GLOPT;
            else
                return SUCCESS;
        }

        void readProblemDataTxt(const char *fname, MatrixXd &ovals, VectorXd &pvals) const
        {
            int num_ovals;
            int num_pvals;
            int num_res;

            ifstream bFile(fname);
            if (!bFile)
                LOG_FATAL << "Cannot open TXT parameter file: " << fname;

            bFile >> num_res;
            bFile >> num_ovals;
            bFile >> num_pvals;

            if (num_ovals != problem._num_rpars)
            {
                bFile.close();
                LOG_FATAL << "Number of residual parameters does not match";
            }

            if (num_pvals != problem._num_ppars)
            {
                bFile.close();
                LOG_FATAL << "Number of problem parameters does not match";
            }

            ovals.resize(num_ovals, num_res);
            pvals.resize(num_pvals);

            for (int j = 0; j < num_res; j++)
                for (int i = 0; i < num_ovals; i++)
                    bFile >> ovals(i, j);

            for (int i = 0; i < num_pvals; i++)
                bFile >> pvals(i);

            bFile.close();
        }

        // Parameters =============================================================

        virtual double getParameter(const string &param_name)
        {
            if (param_name == "sdptol")
                return par_sdptol;
            else if (param_name == "maxnorm")
                return par_maxnorm;
            else if (param_name == "restol")
                return par_restol;
            else if (param_name == "ranktol")
                return par_ranktol;
            else if (param_name == "pivtol")
                return par_pivtol;
            else if (param_name == "verbose")
                return par_verbose;

            LOG_FATAL << "Unrecognized parameter name: " << param_name;
            return 0; // should not reach
        }

        virtual void setParameter(const string &param_name, const double val)
        {
            if (param_name == "sdptol")
                par_sdptol = val;
            else if (param_name == "maxnorm")
                par_maxnorm = val;
            else if (param_name == "restol")
                par_restol = val;
            else if (param_name == "ranktol")
                par_ranktol = val;
            else if (param_name == "pivtol")
                par_pivtol = val;
            else
                LOG_FATAL << "Unrecognized parameter name: " << param_name;
        }

    }; // GPOSolverBase

} // GPOSolver namespace

#endif // GPOSOLVER_BASE_H
