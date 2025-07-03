#include <iostream>
#include <fstream>
#include <functional>
#include "matrix.hpp"
#include "linear_system.hpp"
#include "utils.hpp"
using namespace std;

void printSolution(const std::string &methodName, std::function<std::vector<double>()> solver)
{
    std::vector<double> result;
    double elapsed = 0.0;

    cout << "===" << methodName << "===" << endl;
    try
    {
        elapsed = measureTime([&]() { result = solver(); });

        cout << "Solution: " << endl;
        for (size_t i = 0; i < result.size(); i++)
            cout << "x" << i + 1 << " = " << result[i] << endl;
        cout << "Time: " << elapsed << "ms" << endl;
    }
    catch(const InconsistentSystem &e)
    {
        cout << e.what() << endl;
    }
    catch(const InfiniteSolutions &e)
    {
        cout << e.what() << endl;
    }
    catch(const runtime_error &e)
    {
        cout << "Runtime error: " << e.what() << endl;
    }
    catch(...)
    {
        cout << "Unknown error occured." << endl;
    }
    cout << endl;
}

int main()
{
    generateMatrixToFile("inputMatrix.txt", 3);
    generateVectorToFile("inputVector.txt", 3);

    Matrix A = Matrix::loadFromFile("inputMatrix.txt");

    ifstream in("inputVector.txt");
    if (!in)
    {
        cout << "Cannot open the file 'inputVector.txt'. " << endl;
        return 1;
    }
    size_t n;
    in >> n;
    vector<double> b(n);
    for (size_t i = 0; i < n; i++)
        in >> b[i];

    LinearSystem system(A, b);
    
    Matrix extendedMatrix(n, n + 1);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            extendedMatrix(i, j) = A(i, j);
        extendedMatrix(i, n) = b[i];
    }
    extendedMatrix.print();

    printSolution("Gaussian Elimination", [&]() { return system.solveGaussian(); });
    printSolution("LU Decomposition", [&]() { return system.solveLU(); });
    printSolution("Gauss-Jordan", [&]() { return system.solveGaussJordan(); });
    printSolution("Jacobi", [&]() { return system.solveJacobi(); });

    return 0;
}