#ifndef UTILS_HPP
#define UTILS_HPP

#include <chrono>
#include <utility>
#include <string>

template <typename Func, typename Arg>
double measureTime(Func func, Arg &&arg)
{
    auto start = std::chrono::high_resolution_clock::now();
    func(std::forward<Arg>(arg));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return duration.count() / 1000.0;
}

void generateMatrixToFile(const std::string &filename, size_t n, double minVal = -10.0, double maxVal = 10.0);
void generateVectorToFile(const std::string &filename, size_t n, double minVal = -100.0, double maxVal = 100.0);


#endif