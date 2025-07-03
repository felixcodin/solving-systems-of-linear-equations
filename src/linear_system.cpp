#include "linear_system.hpp"
#include <cmath>
#include <algorithm>

LinearSystem::LinearSystem(const Matrix &A, const std::vector<double> &b) : A(A), b(b) {}

inline std::pair<size_t, size_t> computeRanks(const Matrix &A, const std::vector<double> &b)
{
    size_t n = A.n;
    Matrix R(n, n + 1);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            R(i, j) = A(i, j);
        R(i, n) = b[i];
    }
    
    size_t row = 0;
    for (size_t col = 0; col < n && row < n; col++)
    {
        size_t sel = row;
        for (size_t i = row; i < n; i++)
            if (std::abs(R(i, col)) > std::abs(R(sel, col)))
                sel = i;
        
        if (std::abs(R(sel, col)) < 1e-12)
            continue;
        
        std::swap(R.matrix[row], R.matrix[sel]);

        double pivot = R(row, col);
        for (size_t j = col; j < n; j++)
            R(row, j) /= pivot;
        
        for (size_t i = 0; i < n; i++)
        {
            if (i == row) continue;
            double f = R(i, col);
            for (size_t j = col; j <= n; j++)
                R(i, j) -= f*R(row, j);
        }
    }

    size_t rankA = 0, rankAb = 0;
    for (size_t i = 0; i < n; i++)
    {
        bool nonZeroA = false;
        for (size_t j = 0; j < n; j++)
            if (std::abs(R(i, j)) > 1e-12)
            {
                nonZeroA = true;
                break;
            }
        
        bool nonZeroAb = nonZeroA || (std::abs(R(i, n)) > 1e-12);
        if (nonZeroA) rankA++;
        if (nonZeroAb) rankAb++;
    }

    return {rankA, rankAb};
}

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
        {
            auto [rankA, rankAb] = computeRanks(A, b);
            if (rankAb > rankA) throw InconsistentSystem();
            if (rankA < A.n) throw InfiniteSolutions();
        }

        for (size_t i = k + 1; i < n; i++)
        {
            double factor = M(i, k) / M(k, k);
            for (size_t j = k; j < n; j++)
                M(i, j) -= factor * M(k, j);
            rhs[i] -= factor * rhs[k];
        }
    }

    std::vector<double> x(n);
    for (size_t idx = n; idx-- > 0;)
    {
        x[idx] = rhs[idx];
        for (size_t j = idx + 1; j < n; j++)
            x[idx] -= M(idx, j) * x[j];
        x[idx] /= M(idx, idx);
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
            {
                auto [rankA, rankAb] = computeRanks(A, b);
                if (rankAb > rankA) throw InconsistentSystem();
                if (rankA < A.n) throw InfiniteSolutions();
            }
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

    for (size_t i = n; i-- > 0;)
    {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; j++)
            x[i] -= U(i, j) * x[j];
        if (std::abs(U(i, i)) < 1e-12)
        {
            auto [rankA, rankAb] = computeRanks(A, b);
            if (rankAb > rankA) throw InconsistentSystem();
            if (rankA < A.n) throw InfiniteSolutions();
        }
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
            if (std::abs(A(i, i)) < 1e-12)
            {
                auto [rankA, rankAb] = computeRanks(A, b);
                if (rankAb > rankA) throw InconsistentSystem();
                if (rankA < A.n) throw InfiniteSolutions();
            }
            x_new[i] = (b[i] - sigma) / A(i, i);
        }

        double err = 0.0;
        for (size_t i = 0; i < n; i++)
            err += std::abs(x_new[i] - x[i]);
        if (err < tol)
            return x_new;

        x = x_new;
    }
    throw std::runtime_error("Jacobi method did not converge after " + std::to_string(maxIter) + " iterations");
}

std::vector<double> LinearSystem::solveGaussJordan()
{
    size_t n = A.n;
    Matrix M = A;
    std::vector<double> rhs = b;

    for (size_t i = 0; i < n; i++)
    {
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; k++)
        {
            if (std::abs(M(k, i)) > std::abs(M(maxRow, i)))
                maxRow = k;
        }

        if (std::abs(M(maxRow, i)) < 1e-12)
        {
            auto [rankA, rankAb] = computeRanks(A, b);
            if (rankAb > rankA) throw InconsistentSystem();
            if (rankA < A.n) throw InfiniteSolutions();
        }

        std::swap(M.matrix[i], M.matrix[maxRow]);
        std::swap(rhs[i], rhs[maxRow]);

        double pivot = M(i, i);
        for (size_t j = 0; j < n; j++)
            M(i, j) /= pivot;
        rhs[i] /= pivot;

        for (size_t k = 0; k < n; k++)
        {
            if (k == i) continue;
            double factor = M(k, i);
            for (size_t j = 0; j < n; j++)
                M(k, j) -= factor * M(i, j);
            rhs[k] -= factor * rhs[i];
        }
    }

    return rhs; 
}
