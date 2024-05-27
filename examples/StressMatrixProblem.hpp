/*
   This file was automatically generated on November 26, 2023 22:33:59
   by GpoSolver, (c) 2014--2015, Jan Heller, hellej1@cmp.felk.cvut.cz,
   Czech Technical University, Prague, Czech Republic

   Edited by 2023 Peng Yi, <yipeng3@mail2.sysu.edu.cn>
*/

#ifndef STRESS_MATRIX_PROBLEM_H
#define STRESS_MATRIX_PROBLEM_H

#include <gposolver/gpoproblem.h>
#include <vector>

namespace GpoSolver
{

    class StressMatrixProblem : public GpoProblem
    {
    private:
        int num_blocks;
        int num_cons;
        int num_pcons;
        int num_pvars;
        int num_rpars;
        int num_ppars;
        int problem_size;
        int relax_order;
        int max_degree;
        int rank_shift;
        int *block_sizes;

        int num_obj_data;
        int num_obj_vars;
        objData *obj_data;
        objVar *obj_vars;

        int num_con_data;
        int num_con_vars;
        conData *con_data;
        conVar *con_vars;

        std::vector<Eigen::MatrixXd> _M_mat = {};

    public:
        int getNumBlocks(void)
        {
            return num_blocks;
        }

        int getNumConstraints(void)
        {
            return num_cons;
        }

        int getNumProblemVariables(void)
        {
            return num_pvars;
        }

        int getNumPolyConstraints(void)
        {
            return num_pcons;
        }

        int getNumObjectiveParams(void)
        {
            return num_rpars;
        }

        int getNumConstraintsParams(void)
        {
            return num_ppars;
        }

        int getProblemSize(void)
        {
            return problem_size;
        }

        int getRelaxationOrder(void)
        {
            return relax_order;
        }

        int getMaxDegree(void)
        {
            return max_degree;
        }

        int getRankShift(void)
        {
            return rank_shift;
        }

        const int *getBlockSizes(void)
        {
            return block_sizes;
        }

        int getBlockSize(int i)
        {
            return block_sizes[i];
        }

        int getNumSdpObjectiveData(void)
        {
            return num_obj_data;
        }

        int getNumSdpObjectiveVariables(void)
        {
            return num_obj_vars;
        }

        const objData *getSdpObjectiveData(void)
        {
            return obj_data;
        }

        const objVar *getSdpObjectiveVariables(void)
        {
            return obj_vars;
        }

        const conData *getSdpConstraintsData(void)
        {
            return con_data;
        }

        const conVar *getSdpConstraintsVariables(void)
        {
            return con_vars;
        }
        int getNumSdpConstraintsData(void)
        {
            return num_con_data;
        }

        int getNumSdpConstraintsVariables(void)
        {
            return num_con_vars;
        }

        void updateAuxiliaries(const double *, double *);
        double evalPolyObjective(const double *, const double *, const double *);
        void evalPolyConstraints(const double *, const double *, double *);
        void evalSdpObjectiveVariables(const double *, const double *, double **);
        void evalSdpConstraintsVariables(const double *, const double *, double **);
        void initUserDefined(std::vector<Eigen::MatrixXd> M_mat);

    }; // StressMatrixProblem

} // GPOSolver namespace
#endif // STRESS_MATRIX_PROBLEM_H
