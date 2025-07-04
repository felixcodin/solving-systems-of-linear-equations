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
    cout << "======Choose one======" << endl;
    int key;
    cout << "1 = Random systems of linear equations." << endl;
    cout << "2 = Enter systems of linear equations." << endl;
    cin >> key;
    
    size_t sz;
    cout << "Enter simul equation number of unknowns (> 0): ";
    cin >> sz;
    if (sz <= 0)
    throw runtime_error("Number of unknown <= 0");
    
    
    vector<vector<double>> tempMatrix(sz, vector<double>(sz + 1, 0.0));
    switch(key)
    {
    case 1:
        cout << "Enter (n, m) of value to get generate: ";
        int minV, maxV;
        cin >> minV >> maxV;
        generateMatrixToFile("inputMatrix.txt", sz, minV, maxV);
        generateVectorToFile("inputVector.txt", sz, minV, maxV);
        break;
    
    case 2:
        for (size_t i = 0; i < sz; i++)
        {
            size_t j = 0;
            for (j; j < sz + 1; j++)
                cin >> tempMatrix[i][j];
        }
        
        {
            ofstream outMatrix("inputMatrix.txt");
            if (!outMatrix) throw runtime_error("Can't open inputMatrix.txt");
            ofstream outVector("inputVector.txt");
            if (!outVector) throw runtime_error("Can't open inputVector.txt");

            outMatrix << sz << endl; 
            outVector << sz << endl;
            for (size_t i = 0; i < sz; i++)
            {
                size_t j = 0;
                for (j; j < sz; j++)
                    outMatrix << tempMatrix[i][j] << " ";
                outMatrix << endl;
                outVector << tempMatrix[i][j] << " ";
            }
            outVector << endl;
        }
        break;

    default:
        cout << "Invalid Option." << endl;
        return 1;
    }

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
    cout << endl;
    extendedMatrix.print();

    printSolution("Gaussian Elimination", [&]() { return system.solveGaussian(); });
    printSolution("LU Decomposition", [&]() { return system.solveLU(); });
    printSolution("Gauss-Jordan", [&]() { return system.solveGaussJordan(); });
    printSolution("Jacobi", [&]() { return system.solveJacobi(); });
    printSolution("Gauss-Seidel", [&]() { return system.solveGaussSeidel(); });

    return 0;
}