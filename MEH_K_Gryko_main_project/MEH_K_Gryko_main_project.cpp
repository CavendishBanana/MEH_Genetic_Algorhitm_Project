// MEH_K_Gryko_main_project.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include "MatrixGeneticAlgorithm.h"
#include <random>
#include <string>
#include <limits>
#include <fstream>
using namespace std;

typedef int MatrixFieldType;


int getTargetFunValueNoHalf(MatrixFieldType** matrix, unsigned int matrixSize)
{
    //cout << "main: target function begin" << endl;
    MatrixFieldType totalSum = 0;
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

double getTargetFunValueHalf(MatrixFieldType** matrix, unsigned int matrixSize)
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

void printMatrix(MatrixFieldType** matrix, unsigned int matrixSize)
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

void saveMatrixToFile(MatrixFieldType** matrix, unsigned int matrixSize, std::string filePath)
{
    std::ofstream outdata;
    outdata.open(filePath);
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            outdata << matrix[i][j];
            if (j < matrixSize - 1)
            {
                outdata << "\t";
            }
        }
        outdata << "\n";
    }
    outdata.close();
}
void appendLineToFile(std::string filePath, std::string line, bool appendReturnAtTheEnd=true)
{
    std::ofstream outdata;
    outdata.open(filePath, std::ios::app);
    outdata << line;
    if (appendReturnAtTheEnd)
    {
        outdata << "\n";
    }
    outdata.close();
}

int main()
{
    
    std::string outputFolderPathEndingWithSlash = "C:\\Users\\krzys\\OneDrive\\Desktop\\studia_2_st\\semestr_4\\meh\\projekt\\solverOutput";

    default_random_engine engine;
    uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

    unsigned int matrixSize = 30;
    MatrixFieldType** matrix = new MatrixFieldType* [matrixSize];
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        matrix[i] = new  MatrixFieldType[matrixSize];
        for (unsigned int j = 0; j < matrixSize; j++)
        {
            matrix[i][j] = distribution(engine);

        }
    }

    printMatrix(matrix, matrixSize);
    cout << "===========================================================================================" << endl;
    MatrixGeneticAlgorithm<unsigned short, MatrixFieldType> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, MatrixFieldType>(matrix, matrixSize, 2.0, 15.0, 1.5, 2.0, engine, 0.95, 0.07, 4, 10, 0.2, 0.7);
    
    double startingTargetValue = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(matrix);
    cout << "Target value before optimization: " << startingTargetValue << endl;
    appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals.txt", "before: " + std::to_string(startingTargetValue));
    geneticAlgorithmSolver.useUniformCrossoverInGeneratingChildrenPopulation();
    unsigned int generationsCount = 20000;
    geneticAlgorithmSolver.solveWithNGenerations(generationsCount, true);

    MatrixFieldType** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
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
    saveMatrixToFile(optimizedMatrix, matrixSize, outputFolderPathEndingWithSlash + "matrix_optimized_mnpuc_avm.txt");
    appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals.txt", "after: " + std::to_string(targetValueForBestSolutionRememberedBySolver));
    //cout << "optimized matrix - get target function value by passing array: " << geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix) << endl;
    for (unsigned int i = 0; i < matrixSize; i++)
    {
        delete[] matrix[i];
        delete[] optimizedMatrix[i];
    }
    delete[] optimizedMatrix;
    delete[] matrix;
}


