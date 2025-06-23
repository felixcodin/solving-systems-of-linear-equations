#ifndef UTILS_HPP
#define UTILS_HPP

#include <chrono>
template <typename Func, typename Arg>
double measureTime(Func func, Arg &&arg)
{
    auto start = chrono::high_resolution_clock::now();
    func(forward<Arg>(arg));
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    return duration.count() / 1000.0;
}



#endif