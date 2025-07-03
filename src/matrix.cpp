#include "matrix.hpp"
#include <fstream>
#include <iostream>

Matrix::Matrix(size_t n) : n(n), matrix(n, std::vector<double>(n, 0.0)) {}

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
    for (const auto &row : matrix)
    {
        for (double value : row)
            std::cout << value << ' ';
        std::cout << '\n'; 
    }
}