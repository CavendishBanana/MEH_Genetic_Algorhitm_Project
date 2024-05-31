// MEH_K_Gryko_main_project.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include "MatrixGeneticAlgorithm.h"
#include <random>
#include <limits>
using namespace std;


int getTargetFunValueNoHalf(int** matrix, unsigned int matrixSize)
{
    cout << "main: target function begin" << endl;
    int totalSum = 0;
    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            int localSum = 0;
            if (i - 1 >= 0)
            {
                localSum += matrix[i - 1][j];
            }
            if (j + 1 < matrixSize)
            {
                localSum += matrix[i][j + 1];
            }
            if (i + 1 < matrixSize)
            {
                localSum += matrix[i + 1][j];
            }
            if (j - 1 >= 0)
            {
                localSum += matrix[i][j - 1];
            }

            totalSum += (matrix[i][j] * localSum);
        }
        cout << "total sum for i = " << i << ": " << totalSum << endl;
    }
    cout << "main: target function done" << endl;
    return totalSum;
}

/*
if (i - 1u >= 0u && j - 1u >= 0u)
            {
                localSum += matrix[i - 1u][j - 1u];
            }
            if (i - 1u >= 0u && j + 1u < matrixSize)
            {
                localSum += matrix[i - 1u][j + 1u];
            }
            if (i + 1u < matrixSize && j - 1u >= 0u)
            {
                localSum += matrix[i + 1u][j - 1u];
            }
            if (i + 1u < matrixSize && j + 1u < matrixSize)
            {
                localSum += matrix[i + 1u][j + 1u];
            }

*/

double getTargetFunValueHalf(int** matrix, unsigned int matrixSize)
{
    int totalSum = 0;
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            int localSum = 0;
            if (i - 1 >= 0 && j - 1 >= 0)
            {
                localSum += matrix[i - 1][j - 1];
            }
            if (i - 1 >= 0 && j + 1 < matrixSize)
            {
                localSum += matrix[i - 1][j + 1];
            }
            if (i + 1 < matrixSize && j - 1 >= 0)
            {
                localSum += matrix[i + 1][j - 1];
            }
            if (i + 1 < matrixSize && j + 1 < matrixSize)
            {
                localSum += matrix[i + 1][j + 1];
            }

            totalSum += (matrix[i][j] * localSum);
        }

    }
    return static_cast<double>(totalSum)*0.5;
}

void printMatrix(int** matrix, unsigned int matrixSize)
{
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            cout << matrix[i][j];
            if (j < matrixSize - 1)
            {
                cout << "\t";
            }
        }
        cout << "\n";
    }
}

int main()
{
    unsigned int a = 8u;
    unsigned int b = a >> 1;
    unsigned int c = a >> 2;
    cout << "a: " << a << " b: "<<b<<" c: "<<c<<endl;
    unsigned int x = 1u;
    x = x - 10u;
    cout << "x: "<<x << endl;
    unsigned int x2 = numeric_limits<unsigned int>::max();
    x2 = x2 + 5u;
    cout << "x2: " << x2 << endl;
    

    default_random_engine engine;
    uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

    unsigned int matrixSize = 7;
    int** matrix = new int* [matrixSize];
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        matrix[i] = new int[matrixSize];
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            matrix[i][j] = distribution(engine);

        }
    }

    printMatrix(matrix, matrixSize);
    
    int startingTargetValue = getTargetFunValueNoHalf(matrix, matrixSize);
    cout << "Target value before optimization: " << startingTargetValue << endl;

    MatrixGeneticAlgorithm<unsigned short, int> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, int>(matrix, matrixSize, 2.0, 100.0, 3.0, 3.0, engine, 0.95, 0.07, 4);
    geneticAlgorithmSolver.solveWithNGenerations(10000, true);

    int** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
    int targetValueAfterOptimization = getTargetFunValueNoHalf(optimizedMatrix, matrixSize);
    cout << "target value after optimization: " << targetValueAfterOptimization << endl;
    cout << "Optimized matrix:" << endl;
    printMatrix(optimizedMatrix, matrixSize);

    for (unsigned int i = 0; i < matrixSize; i++)
    {
        delete[] matrix[i];
        delete[] optimizedMatrix[i];
    }
    delete[] optimizedMatrix;
    delete[] matrix;
}


