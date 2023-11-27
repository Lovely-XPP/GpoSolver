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

#ifndef GPOSOLVER_SDPA_H
#define GPOSOLVER_SDPA_H

#include <stdlib.h>

#include <sdpa_call.h>

#include "gposolver_base.h"

namespace GpoSolver
{
    using namespace std;
    using namespace Eigen;

    template <class T>
    class GpoSolverSdpa : public GpoSolverBase<T>
    {
    private:
        using GpoSolverBase<T>::evaluateSdpVariables;
        using GpoSolverBase<T>::problem;
        using GpoSolverBase<T>::cvars_ptr;
        using GpoSolverBase<T>::cvmul;
        using GpoSolverBase<T>::par_verbose;
        using typename GpoSolverBase<T>::SdpStatus;

        SDPA sdpaParams;
        VectorXd cdata_sdpa;

        void getZMatrix(SDPA &sdpaProblem, vector<MatrixXd> &Z)
        {
            int num_blocks = problem._num_blocks;

            Z.resize(num_blocks);
            for (int b = 0; b < num_blocks; b++)
            {
                int bsize = problem._block_sizes[b];
                Z[b] = Map<MatrixXd>(sdpaProblem.getResultXMat(b + 1), bsize, bsize);
            }
        }

    protected:
        SdpStatus solveSdp(double &obj, VectorXd &y, vector<MatrixXd> &Z)
        {
            int num_blocks = problem._num_blocks;
            int num_cons = problem._num_cons;
            int num_cdata = problem._num_con_data;

            // Objective + Constraints variables
            VectorXd objective;
            evaluateSdpVariables(objective);

            SDPA sdpaProblem = sdpaParams;

            sdpaProblem.inputConstraintNumber(num_cons);
            sdpaProblem.inputBlockNumber(num_blocks);

            for (int i = 0; i < num_blocks; i++)
            {
                sdpaProblem.inputBlockSize(i + 1, problem._block_sizes[i]);
                sdpaProblem.inputBlockType(i + 1, SDPA::SDP);
            }

            sdpaProblem.initializeUpperTriangleSpace();

            // Objective
            for (int i = 1; i <= num_cons; i++)
                sdpaProblem.inputCVec(i, double(objective(i)));

            // Constraints
            const GpoProblem::conData *cdata = problem._con_data;
            for (int i = 0; i < num_cdata; i++)
            {
                int ip, jp;
                if (cdata[i].i <= cdata[i].j)
                {
                    ip = cdata[i].i;
                    jp = cdata[i].j;
                }
                else
                {
                    ip = cdata[i].j;
                    jp = cdata[i].i;
                }

                sdpaProblem.inputElement(int(cdata[i].mon), int(cdata[i].con) + 1, ip + 1, jp + 1, double(cdata_sdpa(i)));
            }

            sdpaProblem.initializeUpperTriangle();
            sdpaProblem.initializeSolve();
            sdpaProblem.solve();

            SDPA::PhaseType sol = sdpaProblem.getPhaseValue();

            if ((sol == SDPA::pdOPT) || (sol == SDPA::pdFEAS) || (sol == SDPA::pFEAS) || (sol == SDPA::dFEAS))
            {
                y = Map<VectorXd>(sdpaProblem.getResultXVec(), num_cons);
                obj = sdpaProblem.getDualObj() + objective(0);
                getZMatrix(sdpaProblem, Z);

                sdpaProblem.terminate();
                return GpoSolverBase<T>::SDP_SUCCESS;
            }
            else
            {
                sdpaProblem.terminate();
                return GpoSolverBase<T>::SDP_FAILURE;
            }
        }

    public:
        GpoSolverSdpa()
        {
            GpoSolverBase<T>::init_function();
            const GpoProblem::conData *cdata = problem._con_data;
            const GpoProblem::conVar *cvars = problem._con_vars;
            int num_cdata = problem._num_con_data;
            int num_cvars = problem._num_con_vars;

            cdata_sdpa.resize(num_cdata);

            for (int i = 0; i < num_cdata; i++)
                cdata_sdpa(i) = (cdata[i].mon == 0) ? -cdata[i].val : cdata[i].val;

            for (int i = 0; i < num_cvars; i++)
            {
                const GpoProblem::conData *cdatum;
                int ndata;

                problem.getBlock(cvars[i].con, cvars[i].mon, cdatum, ndata);

                cvmul(i) = (cvars[i].mon == 0) ? -1.0 : 1.0;
                cvars_ptr[i] = cdata_sdpa.data() + (cdatum - cdata) + cvars[i].i;
            }

            sdpaParams.setParameterType(SDPA::PARAMETER_DEFAULT);
            setParameter("verbose", 0);
        }

        GpoSolverSdpa(const int n)
        {
            GpoSolverBase<T>::init_function(n);
            const GpoProblem::conData *cdata = problem._con_data;
            const GpoProblem::conVar *cvars = problem._con_vars;
            int num_cdata = problem._num_con_data;
            int num_cvars = problem._num_con_vars;

            cdata_sdpa.resize(num_cdata);

            for (int i = 0; i < num_cdata; i++)
                cdata_sdpa(i) = (cdata[i].mon == 0) ? -cdata[i].val : cdata[i].val;

            for (int i = 0; i < num_cvars; i++)
            {
                const GpoProblem::conData *cdatum;
                int ndata;

                problem.getBlock(cvars[i].con, cvars[i].mon, cdatum, ndata);

                cvmul(i) = (cvars[i].mon == 0) ? -1.0 : 1.0;
                cvars_ptr[i] = cdata_sdpa.data() + (cdatum - cdata) + cvars[i].i;
            }

            sdpaParams.setParameterType(SDPA::PARAMETER_DEFAULT);
            setParameter("verbose", 0);
        }

        virtual ~GpoSolverSdpa() {}

        double getParameter(const string &param_name)
        {
            if (param_name == "MaxIteration")
                return (int)sdpaParams.getParameterMaxIteration();
            else if (param_name == "EpsilonStar")
                return sdpaParams.getParameterEpsilonStar();
            else if (param_name == "LambdaStar")
                return sdpaParams.getParameterLambdaStar();
            else if (param_name == "OmegaStar")
                return sdpaParams.getParameterOmegaStar();
            else if (param_name == "LowerBound")
                return sdpaParams.getParameterLowerBound();
            else if (param_name == "UpperBound")
                return sdpaParams.getParameterUpperBound();
            else if (param_name == "BetaStar")
                return sdpaParams.getParameterBetaStar();
            else if (param_name == "BetaBar")
                return sdpaParams.getParameterBetaBar();
            else if (param_name == "GammaStar")
                return sdpaParams.getParameterGammaStar();
            else if (param_name == "EpsilonDash")
                return sdpaParams.getParameterEpsilonDash();

            return GpoSolverBase<T>::getParameter(param_name);
        }

        void setParameter(const string &param_name, const double val)
        {
            if (param_name == "verbose")
            {
                (val) ? sdpaParams.setDisplay(stdout) : sdpaParams.setDisplay(NULL);
                par_verbose = val;
            }
            else if (param_name == "MaxIteration")
                sdpaParams.setParameterMaxIteration(int(val));
            else if (param_name == "EpsilonStar")
                sdpaParams.setParameterEpsilonStar(val);
            else if (param_name == "LambdaStar")
                sdpaParams.setParameterLambdaStar(val);
            else if (param_name == "OmegaStar")
                sdpaParams.setParameterOmegaStar(val);
            else if (param_name == "LowerBound")
                sdpaParams.setParameterLowerBound(val);
            else if (param_name == "UpperBound")
                sdpaParams.setParameterUpperBound(val);
            else if (param_name == "BetaStar")
                sdpaParams.setParameterBetaStar(val);
            else if (param_name == "BetaBar")
                sdpaParams.setParameterBetaBar(val);
            else if (param_name == "GammaStar")
                sdpaParams.setParameterGammaStar(val);
            else if (param_name == "EpsilonDash")
                sdpaParams.setParameterEpsilonDash(val);
            else
                GpoSolverBase<T>::setParameter(param_name, val);
        }

        void setDisplay(FILE *display)
        {
            sdpaParams.setDisplay(display);
        }

    }; // GpoSolverSdpa

} // GpoSolver namespace

#endif // GPOSOLVER_SDPA_H
