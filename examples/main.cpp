#include <StressMatrixProblem.hpp>
#include <gposolver/gposolver_sdpa.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>

int main()
{
    Eigen::MatrixXd formation(2, 9);
    Eigen::MatrixXd adjacency(9, 9);
    formation << 8, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    adjacency << 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0;
    // formation << 2.0000, 1.2000, 0, 0, 1.0000, 0, 1.6000, 2.0000, -2.0000, -2.0000;
    // adjacency << 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0;
    /* judge input */
    // Adjacency Matrix should be Square-Matrix
    if (adjacency.cols() != adjacency.rows())
    {
        throw std::logic_error("[Error] (calculate_stress) Adjacency Matrix should be Square-Matrix!");
    }
    // formation dimension should be equal to Adjacency Matrix
    if (formation.cols() != adjacency.cols())
    {
        throw std::logic_error("[Error] (calculate_stress) Formation Agents Count is not equal to Adjacency Matrix dimension!");
    }
    // Adjacency Matrix should be Symmetric Matrix
    if (!adjacency.isApprox(adjacency.transpose()))
    {
        throw std::logic_error("[Error] (calculate_stress) Adjacency Matrix should be Symmetric Matrix!");
    }
    if ((int)adjacency.sum() % 2 != 0)
    {
        throw std::logic_error("[Error] (calculate_stress) Adjacency Matrix should be Symmetric Matrix!");
    }

    // get parameters
    int n = formation.cols();
    int dim = formation.rows();
    int m = adjacency.sum() / 2;
    // init vars
    Eigen::MatrixXd stressMat = Eigen::MatrixXd::Zero(n, n);
    // generate H matrix
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, n);
    int pointer = 0;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            if (adjacency(i, j) != 0)
            {
                H(pointer, i) = -1;
                H(pointer, j) = 1;
                pointer++;
            }
        }
    }
    // svd Pbar
    Eigen::MatrixXd Pbar = Eigen::MatrixXd::Ones(dim + 1, n);
    Pbar.bottomRows(dim) = formation;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Pbar.transpose(), Eigen::ComputeFullU);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd U2 = U.rightCols(n - dim - 1);
    Eigen::MatrixXd Q = U2.transpose();
    // generate E matrix
    Eigen::MatrixXd E = Eigen::MatrixXd::Zero(dim * n, m);
    for (size_t i = 0; i < n; i++)
    {
        E.middleRows(i * dim, dim) = formation * H.transpose() * H.col(i).asDiagonal();
    }
    // SVD E
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(E, Eigen::ComputeFullV);
    Eigen::MatrixXd V = svd2.matrixV();
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(dim * n, m);
    S.topRows(std::min(dim * n, m)) = svd2.singularValues().asDiagonal();
    // get rank(S)
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(S);
    int S_rank = qr.rank();
    // generate B Matrrix
    Eigen::MatrixXd B = V.rightCols(V.cols() - S_rank);
    int num = B.cols();
    // generate M Matrixes
    int M_dim = n - dim - 1;
    std::vector<Eigen::MatrixXd> M = {};
    for (size_t i = 0; i < num; i++)
    {
        M.emplace_back(Q * H.transpose() * B.col(i).asDiagonal() * H * Q.transpose());
    }
    // init solution
    Eigen::MatrixXd x(num, 1);
    // Calculate Stress Matrix
    // setup solver
    Eigen::VectorXd stress;
    GpoSolver::GpoSolverSdpa<GpoSolver::StressMatrixProblem> gposolver(M);
    GpoSolver::GpoSolverStatus status;
    GpoSolver::GpoSolverSolutions sols;
    gposolver.setParameter("restol", 0.1);
    gposolver.setParameter("pivtol", 1e-6);
    // define problem
    double rvals[] = {};
    double pvals[] = {};
    int count = 0;
    // solve problem
    status = gposolver.solve(num, (double *)rvals, (double *)pvals, sols);
    // get solution
    for (GpoSolver::GpoSolverSolutions::iterator it = sols.begin(); it != sols.end(); it++)
    {
        std::vector<double> &sol = *(it);
        for (size_t i = 0; i < x.rows(); i++)
        {
            x(i, 0) = sol[i];
        }
        break;
    }
    // generate Stress
    stress = B * x;
    // generate Stress Matrix
    stress.normalize();
    stressMat = H.transpose() * stress.asDiagonal() * H;
    // if maximum eigen values < 0, then reverse Matrix
    Eigen::VectorXd eigen_values = stressMat.eigenvalues().real();
    std::cout << eigen_values << std::endl;
    double abs_max_eigen_value;
    if (abs(eigen_values.maxCoeff()) > abs(eigen_values.minCoeff()))
    {
        abs_max_eigen_value = eigen_values.maxCoeff();
    }
    else
    {
        abs_max_eigen_value = eigen_values.minCoeff();
    }
    if ((abs_max_eigen_value / abs(abs_max_eigen_value)) < 0)
    {
        stressMat = -stressMat;
    }
    eigen_values = stressMat.eigenvalues().real();
    /*
    if (eigen_values.minCoeff() < -1e-4)
    {
        throw std::runtime_error("Line Stress Matrix is not PSD!");
    }*/
    std::cout << stressMat << std::endl;
}