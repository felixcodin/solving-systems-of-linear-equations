#include "utils.hpp"
#include <random>
#include <stdexcept>
#include <iostream>
#include <fstream>

std::random_device rd;
std::mt19937 gen(rd());

void generateMatrixToFile(const std::string &filename, size_t n, double minVal, double maxVal)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Can't open file.");

    std::uniform_real_distribution<> distrib(minVal, maxVal);

    out << n << "\n";

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            out << distrib(gen) << " ";
        out << "\n";
    }
}

void generateVectorToFile(const std::string &filename, size_t n, double minVal, double maxVal)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Can't open file.");

    std::uniform_real_distribution<> distrib(minVal, maxVal);

    out << n << "\n";

    for (size_t i = 0; i < n; i++)
        out << distrib(gen) << " ";
    out << "\n";
}