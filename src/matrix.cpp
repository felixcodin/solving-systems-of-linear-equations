#include "matrix.hpp"
#include <fstream>
#include <iostream>

Matrix::Matrix(size_t n) : n(n), m(n), matrix(n, std::vector<double>(n, 0.0)) {}
Matrix::Matrix(size_t n, size_t m) : n(n), m(m), matrix(n, std::vector<double>(m, 0.0)) {}

double &Matrix::operator()(size_t i, size_t j) { return matrix[i][j]; }
const double &Matrix::operator()(size_t i, size_t j) const { return matrix[i][j]; }

Matrix Matrix::loadFromFile(const std::string &filename)
{
    std::ifstream in(filename);
    size_t sz;
    in >> sz;
    Matrix M(sz);
    for (size_t i = 0; i < sz; i++)
        for (size_t j = 0; j < sz; j++)
            in >> M(i, j);
    return M;
}

void Matrix::print() const
{
    for (size_t i = 0; i < n; i++)
    {
        size_t j = 0;
        for (j; j < m - 1; j++)
            std::cout << matrix[i][j] << "x" << j + 1 << " ";
        std::cout << "= " << matrix[i][j] << "\n";
    }
}