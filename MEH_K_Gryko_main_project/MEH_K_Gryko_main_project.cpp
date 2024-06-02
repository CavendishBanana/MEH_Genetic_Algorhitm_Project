// MEH_K_Gryko_main_project.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include "MatrixGeneticAlgorithm.h"
#include <random>
#include <limits>
using namespace std;


int getTargetFunValueNoHalf(long int** matrix, unsigned int matrixSize)
{
    //cout << "main: target function begin" << endl;
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
        //cout << "total sum for i = " << i << ": " << totalSum << endl;
    }
    //cout << "main: target function done" << endl;
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

double getTargetFunValueHalf(long int** matrix, unsigned int matrixSize)
{
    /*
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

    }*/
    return static_cast<double>(getTargetFunValueNoHalf(matrix, matrixSize))*0.5;
}

void printMatrix(long int** matrix, unsigned int matrixSize)
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

    default_random_engine engine;
    uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

    unsigned int matrixSize = 30;
    long int** matrix = new long int* [matrixSize];
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        matrix[i] = new long int[matrixSize];
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            matrix[i][j] = distribution(engine);

        }
    }

    printMatrix(matrix, matrixSize);
    cout << "===========================================================================================" << endl;
    MatrixGeneticAlgorithm<unsigned short, long int> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, long int>(matrix, matrixSize, 2.0, 15.0, 1.5, 2.0, engine, 0.95, 0.07, 4, 10, 0.2);
    
    double startingTargetValue = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(matrix);
    cout << "Target value before optimization: " << startingTargetValue << endl;
    geneticAlgorithmSolver.useNPointUniformMixedCrossoverInGeneratingChildrenPopulation();
    unsigned int generationsCount = 20000;
    geneticAlgorithmSolver.solveWithNGenerations(generationsCount, true);

    long int** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
    double targetValueForPassedOptimizedMatrix = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix);
    double targetValueForBestSolutionRememberedBySolver = geneticAlgorithmSolver.getCurrentBestSolutionTargetFunctionValue();
    cout << "target value for passed optimized matrix: " << targetValueForPassedOptimizedMatrix << endl;
    cout << "best target value remembered by solver " << targetValueForBestSolutionRememberedBySolver << endl;
    cout << "Their difference: " << targetValueForPassedOptimizedMatrix - targetValueForBestSolutionRememberedBySolver << endl;
    //double targetValueAfterOptimization = static_cast<double>(getTargetFunValueNoHalf(optimizedMatrix, matrixSize))*0.5;
    //cout << "solutions equal: " << (targetValueForPassedOptimizedMatrix == targetValueAfterOptimization) << " difference: " << targetValueAfterOptimization - targetValueForPassedOptimizedMatrix << endl;
    //cout << "both solver solutions equal: " << (targetValueForBestSolutionRememberedBySolver == targetValueForPassedOptimizedMatrix) << " diffrence: " << targetValueForBestSolutionRememberedBySolver - targetValueForPassedOptimizedMatrix << endl;
    //cout << "target value after optimization: " << targetValueAfterOptimization << endl;
    cout << "===========================================================================================" << endl;
    cout << "Optimized matrix:" << endl;
    printMatrix(optimizedMatrix, matrixSize);
    //cout << "optimized matrix - get target function value by passing array: " << geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix) << endl;
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        delete[] matrix[i];
        delete[] optimizedMatrix[i];
    }
    delete[] optimizedMatrix;
    delete[] matrix;
}


