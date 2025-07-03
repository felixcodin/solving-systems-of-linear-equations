#include <iostream>
#include <fstream>
#include "matrix.hpp"
#include "linear_system.hpp"
#include "utils.hpp"
using namespace std;


int main()
{
    // generateMatrixToFile("inputMatrix.txt", 3);
    // generateVectorToFile("inputVector.txt", 3);

    Matrix A = Matrix::loadFromFile("inputMatrix.txt");
    ifstream in("inputVector.txt");
    size_t n;
    in >> n;
    vector<double> b(n);
    for (size_t i = 0; i < n; i++)
        in >> b[i];

    LinearSystem hpt(A, b);

    vector<double> x = hpt.solveGaussian();
    cout << "Nghiem bang Gauss:\n";
    for (double v : x) cout << v << " ";
    cout << endl;

    return 0;
}