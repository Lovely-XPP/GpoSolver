#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>

#ifdef CSDP_SOLVER
#include <gposolver/gposolver_csdp.h>
#elif MOSEK_SOLVER
#include <gposolver/gposolver_mosek.h>
#elif SDPA_SOLVER
#include <gposolver/gposolver_sdpa.h>
#else
#error Define CSDP_SOLVER, MOSEK_SOLVER, or SDPA_SOLVER macros
#endif

#include "stress_matrix_problem.h"

using namespace GpoSolver;

int main(int argc, const char *argv[])
{
    #ifdef CSDP_SOLVER
        GpoSolverCsdp<StressMatrixProblem> gposolver(5);
    #elif MOSEK_SOLVER
        GpoSolverMosek<StressMatrixProblem> gposolver(5);

        gposolver.setParameter("restol", 0.1);
    #elif SDPA_SOLVER
        GpoSolverSdpa<StressMatrixProblem> gposolver(5);

        gposolver.setParameter("restol", 0.1);
        gposolver.setParameter("pivtol", 1e-6);
    #endif
    try
    {
        GpoSolverStatus status;
        GpoSolverSolutions sols;

        double rvals[] = {};
        double pvals[] = {-0.350710461573, 0.578570554153, -0.171793999573, -0.132228241552, -0.053096725511, 0.578570554153, -1.332395486659, -0.017560167271, -0.082037217330, -0.210991317447, -0.171793999573, -0.017560167271, -0.323840213592, -0.303825992990, -0.263797551786, -0.132228241552, -0.082037217330, -0.303825992990, -0.288277095959, -0.257179301896, -0.053096725511, -0.210991317447, -0.263797551786, -0.257179301896, -0.243942802115, -0.093419548058, -0.231131145200, -0.352564669783, -0.341215064528, -0.318515854017, -0.231131145200, 0.674306431035, 0.120127183214, 0.145586124648, 0.196504007517, -0.352564669783, 0.120127183214, -0.540234722100, -0.499488987898, -0.417997519495, -0.341215064528, 0.145586124648, -0.499488987898, -0.460116610945, -0.381371857038, -0.318515854017, 0.196504007517, -0.417997519495, -0.381371857038, -0.308120532124};
        // Call GpoSolver
        status = gposolver.solve(2, (double *)rvals, (double *)pvals, sols);
        cout << "GpoSolver: " << sols.size() << " solution(s) extracted" << endl;

        if (sols.size() == 0)
            return 0;
        else if (status == SUCCESS_GLOPT)
            cout << "Global optimality certified numerically" << endl;
        else if (status == SUCCESS)
            cout << "Alas, global optimality could not be certified" << endl;

        int c = 1;
        for (GpoSolverSolutions::iterator it = sols.begin(); it != sols.end(); it++)
        {
            std::vector<double> &sol = *(it);
            cout << "Solution " << c++ << ": ";
            for (int i = 0; i < sol.size(); i++)
                cout << sol[i] << " ";
            cout << endl;
        }
    }
    catch (exception &e)
    {
        cerr << "GpoSolver exception: " << e.what() << endl;
        return 1;
    }

    return 0;
}
