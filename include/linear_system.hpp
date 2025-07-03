#ifndef LINEAR_SYSTEM_HPP
#define LINEAR_SYSTEM_HPP

#include "matrix.hpp"
#include "exceptions.hpp"
#include <vector>

class LinearSystem
{
private:
    Matrix A;
    std::vector<double> b;

    void luDecompose(Matrix &L, Matrix &U, std::vector<int> &perm);
    std::vector<double> luSolve(const Matrix &L, const Matrix &U, const std::vector<int> &perm);

public:
    LinearSystem(const Matrix &A, const std::vector<double> &b);

    std::vector<double> solveGaussian(bool partialPivot = true);
    std::vector<double> solveGaussJordan();
    std::vector<double> solveLU();
    std::vector<double> solveJacobi(int maxIter = 1000, double tol = 1e-9);
    std::vector<double> solveGaussSeidel(int maxIter = 1000, double tol = 1e-9);
    Matrix invert();

};

#endif