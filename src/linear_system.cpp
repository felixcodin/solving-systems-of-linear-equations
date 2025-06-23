#include "linear_system.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

LinearSystem::LinearSystem(const Matrix &A, const std::vector<double> &b) : A(A), b(b) {}

std::vector<double> LinearSystem::solveGaussian(bool partialPivot)
{
    size_t n = A.n;
    Matrix M = A;
    std::vector<double> rhs = b;

    for (size_t k = 0; k < n; k++)
    {
        if (partialPivot)
        {
            size_t maxRow = k;
            for (size_t i = k + 1; i < n; i++)
            {
                if (std::abs(M(i, k)) > std::abs(M(maxRow, k)))
                    maxRow = i;
            }
            if (maxRow != k)
            {
                std::swap(M.matrix[k], M.matrix[maxRow]);
                std::swap(rhs[k], rhs[maxRow]);
            }
        }

        if (std::abs(M(k, k)) < 1e-12)
            throw std::runtime_error("Matrix is singular or nearly singular");

        for (size_t i = k + 1; i < n; i++)
        {
            double factor = M(i, k) / M(k, k);
            for (size_t j = k; j < n; j++)
                M(i, j) -= factor * M(k, j);
            rhs[i] -= factor * rhs[k];
        }
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = rhs[i];
        for (size_t j = i + 1; j < n; j++)
            x[i] -= M(i, j) * x[j];
        x[i] /= M(i, i);
    }
    return x;
}

void LinearSystem::luDecompose(Matrix &L, Matrix &U, std::vector<int> &perm)
{
    size_t n = A.n;
    L = Matrix(n);
    U = Matrix(n);
    perm.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        perm[i] = i;
        L(i, i) = 1.0;
    }

    Matrix M = A;
    for (size_t k = 0; k < n; k++)
    {
        for (size_t j = k; j < n; j++)
        {
            double sum = 0.0;
            for (size_t s = 0; s < k; s++)
                sum += L(k, s) * U(s, j);
            U(k, j) = M(k, j) - sum;
        }

        for (size_t i = k + 1; i < n; i++)
        {
            double sum = 0.0;
            for (size_t s = 0; s < k; s++)
                sum += L(i, s) * U(s, k);
            if (std::abs(U(k, k)) < 1e-12)
                throw std::runtime_error("Matrix is singular or nearly singular");
            L(i, k) = (M(i, k) - sum) / U(k, k);
        }
    }
}

std::vector<double> LinearSystem::luSolve(const Matrix &L, const Matrix &U, const std::vector<int> &perm)
{
    size_t n = A.n;
    std::vector<double> y(n), x(n);

    for (size_t i = 0; i < n; i++)
        y[i] = b[i];

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < i; j++)
            y[i] -= L(i, j) * y[j];
    }

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; j++)
            x[i] -= U(i, j) * x[j];
        if (std::abs(U(i, i)) < 1e-12)
            throw std::runtime_error("Matrix is singular or nearly singular");
        x[i] /= U(i, i);
    }
    return x;
}

std::vector<double> LinearSystem::solveLU()
{
    Matrix L(A.n), U(A.n);
    std::vector<int> perm;
    luDecompose(L, U, perm);
    return luSolve(L, U, perm);
}

std::vector<double> LinearSystem::solveJacobi(int maxIter, double tol)
{
    size_t n = A.n;
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for (int iter = 0; iter < maxIter; iter++)
    {
        for (size_t i = 0; i < n; i++)
        {
            double sigma = 0.0;
            for (size_t j = 0; j < n; j++)
            {
                if (j != i)
                    sigma += A(i, j) * x[j];
            }
            if (std::abs(A(i, i) < 1e-12))
                throw std::runtime_error("Zero diagonal element in Jacobi method");
            x_new[i] = (b[i] - sigma) / A(i, i);
        }

        double err = 0.0;
        for (size_t i = 0; i < n; i++)
            err += std::abs(x_new[i] - x[i]);
        if (err < tol)
            return x_new;

        x = x_new;
    }
    return x_new;
}