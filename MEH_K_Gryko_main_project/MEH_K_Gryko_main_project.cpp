// Krzysztof Gryko - projekt na zaliczenie przedmiotu metaheurystyki - temat nr 20
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
            std::cout << matrix[i][j];
            if (j < matrixSize - 1)
            {
                std::cout << "\t";
            }
        }
        std::cout << "\n";
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

    /* 
    //------------------------------- EXPERIMENT 3 - MIXED CROSSOVER WITH DIFFERENT TYPES OF MUTATIONS, INCREASED MUTATION PROBABILITY AND INCREASED NO. OF SPECIMENS IN INVERSION MUTATION ----------

    std::string basePath = "C:\\Users\\krzys\\OneDrive\\Desktop\\studia_2_st\\semestr_4\\meh\\projekt\\solverOutput\\";
    std::cout<< "Experiment 3 - test mixed N Point altering crossover with different mutations, increased probability of mutations and increased no. of specimens in inversion mutation"<<std::endl;
    for (unsigned int ii = 0u; ii < 3u; ii++)
    {

        std::string crossoverTypeMsg = "mixed N point uniform";
        std::string mutationTypeMsg = "";
       
        std::string outputFolderPathEndingWithSlash = basePath;
        outputFolderPathEndingWithSlash += "exper_3_mixed_npoint_uniform_";
        switch (ii)
        {
        case 0:
            outputFolderPathEndingWithSlash += "altering\\";
            mutationTypeMsg = "altering";
            break;
        case 1:
            outputFolderPathEndingWithSlash += "inversion\\";
            mutationTypeMsg = "inversion";
            break;
        case 2:
            outputFolderPathEndingWithSlash += "mixed_altering_inversion\\";
            mutationTypeMsg = "mixed altering inversion";
            break;
        default:
            break;
        }

        

        default_random_engine engine;
        uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

        unsigned int matrixSize = 30;
        MatrixFieldType** matrix = new MatrixFieldType * [matrixSize];
        for (unsigned int i = 0; i < matrixSize; i++)
        {
            matrix[i] = new  MatrixFieldType[matrixSize];
            for (unsigned int j = 0; j < matrixSize; j++)
            {
                matrix[i][j] = distribution(engine);

            }
        }

        printMatrix(matrix, matrixSize);
        saveMatrixToFile(matrix, matrixSize, outputFolderPathEndingWithSlash + "unoptimized_matrix.txt");
        std::cout << "===========================================================================================" << endl;
        MatrixGeneticAlgorithm<unsigned short, MatrixFieldType> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, MatrixFieldType>(matrix, matrixSize, 2.0, 15.0, 1.5, 2.0, engine, 0.95, 0.11, 4, 10, 0.2, 0.7, 5);

        std::cout << "crossover type: " << crossoverTypeMsg << ", mutation type: " << mutationTypeMsg << std::endl;

        double startingTargetValue = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(matrix);
        std::cout << "Target value before optimization: " << startingTargetValue << endl;
        appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "before: " + std::to_string(startingTargetValue));

        geneticAlgorithmSolver.useNPointUniformMixedCrossoverInGeneratingChildrenPopulation();

        switch (ii)
        {
        case 0:
            geneticAlgorithmSolver.useAlteringValueMtuation();
            break;
        case 1:
            geneticAlgorithmSolver.useInversionMutation();
            break;
        case 2:
            geneticAlgorithmSolver.useMixedMutation();
            break;
        default:
            break;
        }

       

        unsigned int generationsCount = 20000;
        geneticAlgorithmSolver.solveWithNGenerations(generationsCount, true, true, outputFolderPathEndingWithSlash);

        MatrixFieldType** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
        double targetValueForPassedOptimizedMatrix = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix);
        double targetValueForBestSolutionRememberedBySolver = geneticAlgorithmSolver.getCurrentBestSolutionTargetFunctionValue();
        std::cout << "target value for passed optimized matrix: " << targetValueForPassedOptimizedMatrix << endl;
        std::cout << "best target value remembered by solver " << targetValueForBestSolutionRememberedBySolver << endl;
        std::cout << "Their difference: " << targetValueForPassedOptimizedMatrix - targetValueForBestSolutionRememberedBySolver << endl;

        std::cout << "===========================================================================================" << endl;
        std::cout << "Optimized matrix:" << endl;
        printMatrix(optimizedMatrix, matrixSize);
        saveMatrixToFile(optimizedMatrix, matrixSize, outputFolderPathEndingWithSlash + "optimized_matrix.txt");
        appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "after: " + std::to_string(targetValueForBestSolutionRememberedBySolver));

        for (unsigned int i = 0; i < matrixSize; i++)
        {
            delete[] matrix[i];
            delete[] optimizedMatrix[i];
        }
        delete[] optimizedMatrix;
        delete[] matrix;

    }

    // ---------------------------------------------------------------- EXPERIMENT 3 DONE -----------------------------------------------------------------------------
    */

    /*
    // ------------------------------------ EXPERIMENT 2 - MIXED MUTATION WITH INCREASED NO OF CPECIMENS IN INVERSION MUTATION --------------------------------- 
    std::cout<< "Experiment 2 - test mixed altering inversion mutation with different crossovers and increased no. of specimens in inversion mutation"<<std::endl;
    for (unsigned int ii = 0u; ii < 3u; ii++)
    {
            std::string crossoverTypeMsg = "";
            std::string mutationTypeMsg = "mixed altering inversion";
            
            std::string outputFolderPathEndingWithSlash = basePath;
            switch (ii)
            {
            case 0:
                outputFolderPathEndingWithSlash += "exper_2_npoint_";
                crossoverTypeMsg = "N point";
                break;
            case 1:
                outputFolderPathEndingWithSlash += "exper_2_uniform_";
                crossoverTypeMsg = "uniform";
                break;
            case 2:
                outputFolderPathEndingWithSlash += "exper_2_mixed_npoint_uniform_";
                crossoverTypeMsg = "mixed N point uniform";
                break;
            default:
                break;
            }

            outputFolderPathEndingWithSlash += "mixed_altering_inversion\\";

            default_random_engine engine;
            uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

            unsigned int matrixSize = 30;
            MatrixFieldType** matrix = new MatrixFieldType * [matrixSize];
            for (unsigned int i = 0; i < matrixSize; i++)
            {
                matrix[i] = new  MatrixFieldType[matrixSize];
                for (unsigned int j = 0; j < matrixSize; j++)
                {
                    matrix[i][j] = distribution(engine);

                }
            }

            printMatrix(matrix, matrixSize);
            saveMatrixToFile(matrix, matrixSize, outputFolderPathEndingWithSlash + "unoptimized_matrix.txt");
            std::cout << "===========================================================================================" << endl;
            MatrixGeneticAlgorithm<unsigned short, MatrixFieldType> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, MatrixFieldType>(matrix, matrixSize, 2.0, 15.0, 1.5, 2.0, engine, 0.95, 0.07, 4, 10, 0.2, 0.7, 5);

            double startingTargetValue = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(matrix);
            std::cout << "Target value before optimization: " << startingTargetValue << endl;

            std::cout << "crossover type: " << crossoverTypeMsg << ", mutation type: " << mutationTypeMsg << std::endl;

            appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "before: " + std::to_string(startingTargetValue));


            switch (ii)
            {
            case 0:
                geneticAlgorithmSolver.useNPointCrossoverInGeneratingChildrenPopulation();
                break;
            case 1:
                geneticAlgorithmSolver.useUniformCrossoverInGeneratingChildrenPopulation();
                break;
            case 2:
                geneticAlgorithmSolver.useNPointUniformMixedCrossoverInGeneratingChildrenPopulation();
                break;
            default:
                break;
            }

            geneticAlgorithmSolver.useMixedMutation();

            unsigned int generationsCount = 20000;
            geneticAlgorithmSolver.solveWithNGenerations(generationsCount, true, true, outputFolderPathEndingWithSlash);

            MatrixFieldType** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
            double targetValueForPassedOptimizedMatrix = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix);
            double targetValueForBestSolutionRememberedBySolver = geneticAlgorithmSolver.getCurrentBestSolutionTargetFunctionValue();
            std::cout << "target value for passed optimized matrix: " << targetValueForPassedOptimizedMatrix << endl;
            std::cout << "best target value remembered by solver " << targetValueForBestSolutionRememberedBySolver << endl;
            std::cout << "Their difference: " << targetValueForPassedOptimizedMatrix - targetValueForBestSolutionRememberedBySolver << endl;

            std::cout << "===========================================================================================" << endl;
            std::cout << "Optimized matrix:" << endl;
            printMatrix(optimizedMatrix, matrixSize);
            saveMatrixToFile(optimizedMatrix, matrixSize, outputFolderPathEndingWithSlash + "optimized_matrix.txt");
            appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "after: " + std::to_string(targetValueForBestSolutionRememberedBySolver));

            for (unsigned int i = 0; i < matrixSize; i++)
            {
                delete[] matrix[i];
                delete[] optimizedMatrix[i];
            }
            delete[] optimizedMatrix;
            delete[] matrix;

    }

    // ---------------------------------------------- experiment 2 done ------------------------------------------------------
    */


    

    // --------------------- DEFAULT EXPERIMENT ----------------------


    std::string basePath = "C:\\Users\\krzys\\OneDrive\\Desktop\\studia_2_st\\semestr_4\\meh\\projekt\\solverOutput\\";
    std::cout << "Default experiment - test all combinations of crossover and mutation" << std::endl;
    for (unsigned int ii = 0u; ii < 3u; ii++)
    {
        for (unsigned int jj = 0u; jj < 3u; jj++)
        {
            std::string crossoverTypeMsg = "";
            std::string mutationTypeMsg = "";
            std::string outputFolderPathEndingWithSlash = basePath;
            switch (ii)
            {
            case 0:
                outputFolderPathEndingWithSlash += "npoint_";
                crossoverTypeMsg = "N point";
                break;
            case 1:
                outputFolderPathEndingWithSlash += "uniform_";
                crossoverTypeMsg = "uniform";
                break;
            case 2:
                outputFolderPathEndingWithSlash += "mixed_npoint_uniform_";
                crossoverTypeMsg = "mixed N point uniform";
                break;
            default:
                break;
            }

            switch (jj)
            {
            case 0:
                outputFolderPathEndingWithSlash += "altering\\";
                mutationTypeMsg = "altering";
                break;
            case 1:
                outputFolderPathEndingWithSlash += "inversion\\";
                mutationTypeMsg = "inversion";
                break;
            case 2:
                outputFolderPathEndingWithSlash += "mixed_altering_inversion\\";
                mutationTypeMsg = "mixed altering inversion";
            default:
                break;
            }

            default_random_engine engine;
            uniform_int_distribution<int> distribution = uniform_int_distribution<int>(-100, 100);

            unsigned int matrixSize = 30;
            MatrixFieldType** matrix = new MatrixFieldType * [matrixSize];
            for (unsigned int i = 0; i < matrixSize; i++)
            {
                matrix[i] = new  MatrixFieldType[matrixSize];
                for (unsigned int j = 0; j < matrixSize; j++)
                {
                    matrix[i][j] = distribution(engine);

                }
            }

            printMatrix(matrix, matrixSize);
            saveMatrixToFile(matrix, matrixSize, outputFolderPathEndingWithSlash + "unoptimized_matrix.txt");
            std::cout << "===========================================================================================" << endl;
            MatrixGeneticAlgorithm<unsigned short, MatrixFieldType> geneticAlgorithmSolver = MatrixGeneticAlgorithm<unsigned short, MatrixFieldType>(matrix, matrixSize, 2.0, 15.0, 1.5, 2.0, engine, 0.95, 0.07, 4, 10, 0.2, 0.7, 2);

            double startingTargetValue = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(matrix);
            std::cout << "Target value before optimization: " << startingTargetValue << endl;

            std::cout << "crossover type: " << crossoverTypeMsg << ", mutation type: " << mutationTypeMsg << std::endl;
            appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "before: " + std::to_string(startingTargetValue));


            switch (ii)
            {
            case 0:
                geneticAlgorithmSolver.useNPointCrossoverInGeneratingChildrenPopulation();
                break;
            case 1:
                geneticAlgorithmSolver.useUniformCrossoverInGeneratingChildrenPopulation();
                break;
            case 2:
                geneticAlgorithmSolver.useNPointUniformMixedCrossoverInGeneratingChildrenPopulation();
                break;
            default:
                break;
            }

            switch (jj)
            {
            case 0:
                geneticAlgorithmSolver.useAlteringValueMtuation();
                break;
            case 1:
                geneticAlgorithmSolver.useInversionMutation();
                break;
            case 2:
                geneticAlgorithmSolver.useMixedMutation();
            default:
                break;
            }
            
            unsigned int generationsCount = 20000;
            geneticAlgorithmSolver.solveWithNGenerations(generationsCount, true, true, outputFolderPathEndingWithSlash);

            MatrixFieldType** optimizedMatrix = geneticAlgorithmSolver.getCurrentBestSolution();
            double targetValueForPassedOptimizedMatrix = geneticAlgorithmSolver.getTargetFunctionValueForPassedMatrix(optimizedMatrix);
            double targetValueForBestSolutionRememberedBySolver = geneticAlgorithmSolver.getCurrentBestSolutionTargetFunctionValue();
            std::cout << "target value for passed optimized matrix: " << targetValueForPassedOptimizedMatrix << endl;
            std::cout << "best target value remembered by solver " << targetValueForBestSolutionRememberedBySolver << endl;
            std::cout << "Their difference: " << targetValueForPassedOptimizedMatrix - targetValueForBestSolutionRememberedBySolver << endl;

            std::cout << "===========================================================================================" << endl;
            std::cout << "Optimized matrix:" << endl;
            printMatrix(optimizedMatrix, matrixSize);
            saveMatrixToFile(optimizedMatrix, matrixSize, outputFolderPathEndingWithSlash + "optimized_matrix.txt");
            appendLineToFile(outputFolderPathEndingWithSlash + "target_fun_vals_before_after.txt", "after: " + std::to_string(targetValueForBestSolutionRememberedBySolver));

            for (unsigned int i = 0; i < matrixSize; i++)
            {
                delete[] matrix[i];
                delete[] optimizedMatrix[i];
            }
            delete[] optimizedMatrix;
            delete[] matrix;
        }
    }

    // ----------------------------- DEFAULT EXPERIMENT DONE ----------------------------

    
}


