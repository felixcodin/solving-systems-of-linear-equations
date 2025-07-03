#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <string>

class Matrix
{
public:
    size_t n;
    size_t m;
    std::vector<std::vector<double>> matrix;

    Matrix(size_t n);
    Matrix(size_t n, size_t m);

    double &operator()(size_t i, size_t j);
    const double &operator()(size_t i, size_t j) const;

    static Matrix loadFromFile(const std::string &filename);

    void print() const;
};


#endif