/*
   This file was automatically generated on November 25, 2023 11:44:53
   by GpoSolver, (c) 2014--2015, Jan Heller, hellej1@cmp.felk.cvut.cz,
   Czech Technical University, Prague, Czech Republic

   This solver takes:
       0 residual parameters:
       0 problem parameters :
*/

#include "stress_matrix_problem.h"

using namespace GpoSolver;
using namespace Eigen;

void StressMatrixProblem::updateAuxiliaries(const double *pars, double *aux) const
{
    int n = mat_size;
    for (size_t i = 0; i < n * n; i++)
    {
        aux[2 * i] = pars[i];
        aux[2 * i + 1] = pars[i + n * n];
    }
}

double StressMatrixProblem::evalPolyObjective(const double *rpars, const double *ppars, const double *vars) const
{
    double t0 = vars[2];

    return t0;
}

void StressMatrixProblem::evalPolyConstraints(const double *vars, const double *pars, double *cons) const
{
    int n = mat_size;
    double aux[1 + 2 * n * n];

    updateAuxiliaries(pars, aux);

    double t0 = -vars[0] * vars[0] - vars[1] * vars[1] + 1.0;
    cons[0] = t0;
    MatrixXd M;
    int count = 0;
    M.resize(n, n);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (j == i)
            {
                *(M.data() + count) = vars[0] * aux[2 * count] + vars[1] * aux[2 * count + 1] + vars[2];
                //printf("vars[%d] * aux[%d] + vars[%d] * aux[%d] + vars[2]\n", 0, 2 * count, 1, 2 * count + 1);
            }
            else
            {
                *(M.data() + count) = vars[0] * aux[2 * count] + vars[1] * aux[2 * count + 1];
                //printf("vars[%d] * aux[%d] + vars[%d] * aux[%d]\n", 0, 2 * count, 1, 2 * count + 1);
            }
            count++;
        }
    }
    cons[1] = M.eigenvalues().real().minCoeff();
}

void StressMatrixProblem::evalSdpObjectiveVariables(const double *rpars, const double *ppars, double **vars) const
{
}

void StressMatrixProblem::evalSdpConstraintsVariables(const double *pars, const double *vmul, double **vars) const
{

    int n = mat_size;
    double aux[1 + 2 * n * n];

    updateAuxiliaries(pars, aux);
    int count = 0;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i; j < n; j++)
        {
            *(vars[count]) = vmul[count] * (aux[2 * i + 2 * n * j]);
            //printf("*(vars[%d]) = vmul[%d] * (aux[%d])\n", count, count, 2 * i + 2 * n * j);
            count++;
        }
    }
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i; j < n; j++)
        {
            *(vars[count]) = vmul[count] * (aux[2 * i + 1 + 2 * n * j]);
            //printf("*(vars[%d]) = vmul[%d] * (aux[%d])\n", count, count, 2 * i + 2 * n * j + 1);
            count++;
        }
    }
}

const int StressMatrixProblem::num_blocks = 3;
const int StressMatrixProblem::num_cons = 9;
const int StressMatrixProblem::num_pcons = 2;
const int StressMatrixProblem::num_pvars = 3;
const int StressMatrixProblem::num_rpars = 0;
const int StressMatrixProblem::num_ppars = 0;
const int StressMatrixProblem::problem_size = 0;
const int StressMatrixProblem::relax_order = 1;
const int StressMatrixProblem::max_degree = 2;
const int StressMatrixProblem::rank_shift = 1;
const int StressMatrixProblem::block_sizes[] = {4, 1, 0};

const int StressMatrixProblem::num_obj_data = 1;
const int StressMatrixProblem::num_obj_vars = 0;
const StressMatrixProblem::objData StressMatrixProblem::obj_data[] = {
    {3, 1}};
const StressMatrixProblem::objVar StressMatrixProblem::obj_vars[] = {};
const int StressMatrixProblem::num_con_data = 0;
const int StressMatrixProblem::num_con_vars = 0;
const StressMatrixProblem::conData StressMatrixProblem::con_data[] = {
    {0, 0, 0, 0, 1.0}, {0, 1, 0, 1, 1.0}, {0, 2, 0, 2, 1.0}, {0, 3, 0, 3, 1.0}, {0, 4, 1, 1, 1.0}, {0, 5, 1, 2, 1.0}, {0, 6, 1, 3, 1.0}, {0, 7, 2, 2, 1.0}, {0, 8, 2, 3, 1.0}, {0, 9, 3, 3, 1.0}, {1, 0, 0, 0, 1.0}, {1, 4, 0, 0, -1.0}, {1, 7, 0, 0, -1.0}};
const StressMatrixProblem::conVar StressMatrixProblem::con_vars[] = {};
