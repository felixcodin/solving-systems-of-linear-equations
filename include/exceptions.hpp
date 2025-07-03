#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <stdexcept>

struct InconsistentSystem : std::runtime_error 
{
    InconsistentSystem() : std::runtime_error("No Solution.") {}
};

struct InfiniteSolutions : std::runtime_error
{
    InfiniteSolutions() : std::runtime_error("Infinite Soluiton.") {}
};


#endif