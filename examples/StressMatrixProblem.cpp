/*
   This file was automatically generated on November 25, 2023 11:44:53
   by GpoSolver, (c) 2014--2015, Jan Heller, hellej1@cmp.felk.cvut.cz,
   Czech Technical University, Prague, Czech Republic

   Edited by 2023 Peng Yi, <yipeng3@mail2.sysu.edu.cn>
*/

#include <StressMatrixProblem.hpp>

using namespace GpoSolver;
using namespace Eigen;

void StressMatrixProblem::initUserDefined(std::vector<Eigen::MatrixXd> M_mat)
{
    // int vars
    _M_mat = M_mat;
    unsigned int mat_size = _M_mat[0].cols();
    unsigned int n = _M_mat.size();
    unsigned int var_count = n;

    // define problem
    num_blocks = 3;
    num_cons = (var_count + 2 + 1) * (var_count + 2) / 2 - 1;
    num_pcons = 2;
    num_pvars = var_count + 1;
    num_rpars = 0;
    num_ppars = 0;
    problem_size = n + 2 + mat_size + 1;
    relax_order = 1;
    max_degree = 2;
    rank_shift = 1;

    // block sizes
    block_sizes = new int[3];
    block_sizes[0] = var_count + 2;
    block_sizes[1] = 1;
    block_sizes[2] = mat_size;

    // obj data
    num_obj_data = 1;
    num_obj_vars = 0;
    obj_data = new objData[1];
    obj_vars = new objVar[1];
    obj_data[0] = {var_count + 1, 1};
    obj_vars[0] = {0};

    // generate con_data and con_vars
    // data size: var_count * mat_size (1, x1, x2, ..., t) * mat_size
    unsigned int all_var_count = var_count + 2; // add var 1
    std::vector<conData> vec_con_data = {};

    // genertare matrix: (1, x1, x2, ..., t)
    unsigned int con_data_idx = 0; // idx
    for (unsigned int i = 0; i < all_var_count; i++)
    {
        for (unsigned int j = i; j < all_var_count; j++)
        {
            conData tmp_con = {0, con_data_idx, i, j, 1.0};
            vec_con_data.emplace_back(tmp_con);
            con_data_idx++;
        }
    }
    // generate matrix: condition 2 ||x|| < 1
    conData tmp_con = {1, 0, 0, 0, 1.0};
    vec_con_data.emplace_back(tmp_con);
    for (unsigned int i = all_var_count; i > 2; i--)
    {
        if (i == all_var_count)
        {
            con_data_idx = i;
            conData tmp_con = {1, con_data_idx, 0, 0, -1.0};
            vec_con_data.emplace_back(tmp_con);
            continue;
        }
        con_data_idx += i;
        conData tmp_con = {1, con_data_idx, 0, 0, -1.0};
        vec_con_data.emplace_back(tmp_con);
    }
    // generate matrix: condition 3 mat_size * mat_size
    // generate non-diag elements
    for (unsigned int k = 1; k < (all_var_count - 1); k++)
    {
        for (unsigned int i = 0; i < mat_size; i++)
        {
            for (unsigned int j = i; j < mat_size; j++)
            {
                conData tmp_con = {2, k, i, j, 0.0};
                vec_con_data.emplace_back(tmp_con);
            }
        }
    }
    // generate 1 (diag)
    for (unsigned int i = 0; i < mat_size; i++)
    {
        conData tmp_con = {2, all_var_count - 1, i, i, 1.0};
        vec_con_data.emplace_back(tmp_con);
    }
    // transform data
    num_con_data = vec_con_data.size();
    con_data = new conData[num_con_data];
    for (size_t i = 0; i < num_con_data; i++)
    {
        con_data[i] = vec_con_data[i];
    }

    // generate con_vars
    int single_con_var_size = (mat_size + 1) * mat_size / 2;
    std::vector<conVar> vec_con_vars = {};
    for (int k = 1; k < (all_var_count - 1); k++)
    {
        for (int i = 0; i < single_con_var_size; i++)
        {
            conVar tmp_var = {2, k, i};
            vec_con_vars.emplace_back(tmp_var);
        }
    }
    // transform data
    num_con_vars = vec_con_vars.size();
    con_vars = new conVar[num_con_vars];
    for (size_t i = 0; i < num_con_vars; i++)
    {
        con_vars[i] = vec_con_vars[i];
    }
}

void StressMatrixProblem::updateAuxiliaries(const double *pars, double *aux)
{
    int n = _M_mat.size();
    int count = 0;
    int mat_size = _M_mat[0].cols();
    for (size_t i = 0; i < mat_size; i++)
    {
        for (size_t j = 0; j < mat_size; j++)
        {
            for (size_t k = 0; k < n; k++)
            {
                aux[count] = _M_mat[k](j, i);
                count++;
            }
        }
    }
}

double StressMatrixProblem::evalPolyObjective(const double *rpars, const double *ppars, const double *vars)
{
    double t0 = vars[_M_mat.size()];

    return t0;
}

void StressMatrixProblem::evalPolyConstraints(const double *vars, const double *pars, double *cons)
{
    int mat_size = _M_mat[0].cols();
    int n = _M_mat.size();
    double aux[1 + mat_size * mat_size * n];

    updateAuxiliaries(pars, aux);

    double t0 = 1.0;
    for (size_t i = 0; i < n; i++)
    {
        t0 -= vars[i] * vars[i];
    }
    cons[0] = t0;

    MatrixXd M;
    int count = 0;
    M.resize(mat_size, mat_size);
    // SUM(Mi * xi) + I
    for (size_t i = 0; i < mat_size; i++)
    {
        for (size_t j = 0; j < mat_size; j++)
        {
            if (j == i)
            {
                *(M.data() + count) = vars[n];
            }
            else
            {
                *(M.data() + count) = 0;
            }
            for (size_t k = 0; k < n; k++)
            {
                *(M.data() + count) += vars[k] * aux[n * count + k];
            }
            count++;
        }
    }
    cons[1] = M.eigenvalues().real().minCoeff();
}

void StressMatrixProblem::evalSdpObjectiveVariables(const double *rpars, const double *ppars, double **vars)
{
}

void StressMatrixProblem::evalSdpConstraintsVariables(const double *pars, const double *vmul, double **vars)
{
    int mat_size = _M_mat[0].cols();
    int n = _M_mat.size();
    double aux[1 + mat_size * mat_size * n];

    updateAuxiliaries(pars, aux);
    int count = 0;

    for (size_t k = 0; k < n; k++)
    {
        for (int i = 0; i < mat_size; i++)
        {
            for (int j = i; j < mat_size; j++)
            {
                *(vars[count]) = vmul[count] * (aux[n * i + n * mat_size * j + k]);
                count++;
            }
        }
    }
}