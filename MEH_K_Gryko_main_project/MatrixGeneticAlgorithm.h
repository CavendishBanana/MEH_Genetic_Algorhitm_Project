#pragma once

//#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <new>
#include <algorithm> 
#include <stdexcept>
#include <iomanip>
#include <string>
#include <fstream>
template <typename GeneField, typename MatrixField> class MatrixGeneticAlgorithm
{
private:
	MatrixField** matrixBeforeOptimization;
	MatrixField** matrixCopyForFitnessEvaluation;
	int* matrixCopyHashingRowsArray;
	int* matrixCopyHashingColsArray;
	int* matrixCopyHashingRowsArrayWithOffsets;
	int* matrixCopyHashingColsArrayWithOffsets;
	unsigned int inversionMutationGenesCount;
	unsigned int matrixSize;

	unsigned int nPointCrossoverRowsOffset;
	unsigned int nPointCrossoverColsOffset;
	unsigned int nPointCrossoverBaseSplitLength;
	unsigned int nPointCrossoverTotalSplitsCount;

	double mixedCrossoverNPointUniformThreshold;
	
	GeneField* bestSolution;
	MatrixField bestSolutionTargetFunctionValue;

	unsigned int populationSize;
	unsigned int chromosomeLength;
	unsigned int halfChromosomeLength;
	unsigned int* decodingArray;
	double crossoverProbability, mutationProbability;
	double alteringInversionMixedMutationAlteringPrefference;

	void (MatrixGeneticAlgorithm<GeneField, MatrixField>::*crossoverOperationFunctionPointer)(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2);
	void (MatrixGeneticAlgorithm<GeneField, MatrixField>::* mutationOperationFunctionPointer)(GeneField* specimen);

	GeneField decodingArrayRealSize;
	//const GeneField geneFieldOne;
	//const GeneField geneFieldTwo;
	//const GeneField geneFieldZero;
	static constexpr GeneField getGeneFieldZero();
	static constexpr GeneField getGeneFieldOne();
	static constexpr GeneField getGeneFieldTwo();
	static constexpr MatrixField getMatrixFieldZero();
	GeneField** parentPopulation;
	GeneField** childrenPopulation;
	MatrixField* parentPopulationTargetFunVals;
	MatrixField* childrenPopulationTargetFunVals;
	

	//selection process arrays and variables
	GeneField** parentPopulationSelectionHelperArray;
	MatrixField* parentPopulationSelectionTargetFunVals;
	GeneField** childrenPopulationSelectionHelperArray;
	MatrixField* childrenPopulationSelectionTargetFunVals;
	bool* selectionIsInNewParentGenerationArr;
	unsigned int tournamentGroupSize;
	unsigned int* indexesOfTournamentConteandants;
	unsigned int* populationsTournamentConteandants;
	double oneByAlteringMutationUnitSpan, alteringMutationMaxAcceptedValue;

	void fillDecodingArray();
	void generateStartingPopulation(double populationPower, double populationMultplier);
	inline GeneField getGenesCount();
	std::default_random_engine& randomEngine;
	std::uniform_int_distribution<GeneField> distribution;
	void copyMatrixToTargetFunctionArray(MatrixField** source);
	void initOptimizedMatrix(MatrixField** source);
	
	MatrixField targetFunctionValueNoHalf(GeneField* solution);
	void decodeMatrix(GeneField* chromosome);
	void crossoverUniform(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2);
	void crossoverNPoint(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2);
	void nPointUniformMixedCrossover(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2);
	void mutationByInversion(GeneField* specimen);
	void mutationByAlteringGeneValue(GeneField* specimen);
	void alteringValueAndInversionMixedMutation(GeneField* specimen);
	void findParamsForAlteringMutation(unsigned int runsCount = 10000000u);
	void findParamsForNPointCrossover();
	unsigned int findBestSolutionIndexInArray(MatrixField* targetFunctionValues, unsigned int targetFunValsArrayLength);
	void generateChildrenPopulation();
	void tournamentSelectionMiPlusLambda();
	GeneField getDecodingArrayRealSize();
	GeneField getDecodingArrayVisibleSize();
	void setDecodingArrayValue(GeneField geneValue, unsigned int from, unsigned int to);
	GeneField minimalColumnIndex();
	GeneField genesCount;
	void saveSolutionPackageToFile(const std::string& solutionFolder, const std::string& matrixFileName, const std::string& targetFunctionValuesFileName);

public:
	void solveWithNGenerations(unsigned int generationsCount, bool verbose, bool savePartialResults, std::string outputFolderPathEndingWithSlash);
	
	unsigned int* decodeWithArray(unsigned int* returnValueArray, GeneField geneValue);
	MatrixGeneticAlgorithm(MatrixField** staringMatrix, unsigned int startingMatrixSize, double populationMatrSizePower,
		double populationSizeMultiplier, double chromesomeLenMatrSizePower, double chromosomeSizeMultplier, 
		std::default_random_engine& randomEngine, double crossoverProbability, double mutationProbability, 
		unsigned int tournamentGroupSize, unsigned int nPointCrossoverSplitsCount, double mixedCrossoverNPointUniformThresholdForUniform, double mixedAltInvMutationAlterPrefference);

	MatrixGeneticAlgorithm(
		MatrixField** staringMatrix, unsigned int startingMatrixSize, double populationMatrSizePower,
		double populationSizeMultiplier, double chromesomeLenMatrSizePower, double chromosomeSizeMultplier,
		std::default_random_engine& randomEngine, double crossoverProbability, double mutationProbability,
		unsigned int tournamentGroupSize, unsigned int nPointCrossoverSplitsCount, double mixedCrossoverNPointUniformThresholdForUniform,
		double alteringInversionMixedMutationAlteringPrefference, unsigned int inversionMutationGenesCount);

	~MatrixGeneticAlgorithm();
	double targetFunctionValue(GeneField* solution);
	double getCurrentBestSolutionTargetFunctionValue();
	double getTargetFunctionValueForPassedMatrix(MatrixField** passedMatrix);
	MatrixField getCurrentBestSolutionTargetFunctionValueNoHalf();
	MatrixField** getCurrentBestSolution();
	void useNPointCrossoverInGeneratingChildrenPopulation();
	void useUniformCrossoverInGeneratingChildrenPopulation();
	void useNPointUniformMixedCrossoverInGeneratingChildrenPopulation();
	void useAlteringValueMtuation();
	void useInversionMutation();
	void useMixedMutation();
};

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverUniform(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2)
{
	//uniform crossover
	std::uniform_int_distribution<unsigned int> coinFlipping = std::uniform_int_distribution<unsigned int>(1, 10000);
	for (unsigned int i = 0u; i < chromosomeLength; i++)
	{
		unsigned int flip = coinFlipping(randomEngine);
		if (flip <= 5000)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
		}
		else
		{
			child1[i] = parent2[i];
			child2[i] = parent1[i];
		}
	}
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::findParamsForNPointCrossover()
{
	nPointCrossoverBaseSplitLength = (chromosomeLength - (chromosomeLength % nPointCrossoverTotalSplitsCount)) / nPointCrossoverTotalSplitsCount;
	nPointCrossoverRowsOffset = halfChromosomeLength % nPointCrossoverBaseSplitLength;
	nPointCrossoverColsOffset = (chromosomeLength - halfChromosomeLength) % nPointCrossoverBaseSplitLength;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useNPointCrossoverInGeneratingChildrenPopulation()
{
	crossoverOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverNPoint;
}
template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useUniformCrossoverInGeneratingChildrenPopulation()
{
	crossoverOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverNPoint;
}
template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useNPointUniformMixedCrossoverInGeneratingChildrenPopulation()
{
	crossoverOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::nPointUniformMixedCrossover;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useAlteringValueMtuation()
{
	mutationOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByAlteringGeneValue;
}
template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useInversionMutation()
{
	mutationOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByInversion;
}
template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::useMixedMutation()
{
	mutationOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::alteringValueAndInversionMixedMutation;
}


template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverNPoint(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2)
{
	GeneField* child1Backup = child1;
	GeneField* child2Backup = child2;
	unsigned int currentPartitionBeginIdx = 0u;
	for (unsigned int i = 0u; i < nPointCrossoverRowsOffset; i++)
	{
		child1[i] = parent1[i];
		child2[i] = parent2[i];
	}
	currentPartitionBeginIdx += nPointCrossoverRowsOffset;
	while (currentPartitionBeginIdx < halfChromosomeLength)
	{
		unsigned int thisSplitStop = currentPartitionBeginIdx + nPointCrossoverBaseSplitLength;
		thisSplitStop = (thisSplitStop < halfChromosomeLength ? thisSplitStop : halfChromosomeLength);
		unsigned int i = currentPartitionBeginIdx;
		while ( i < thisSplitStop)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
			i++;
		}

		//for (unsigned int i = currentPartitionBeginIdx; i < thisSplitStop; i++)
		//{
		//	child1[i] = parent1[i];
		//	child2[i] = parent2[i];
		//}
		currentPartitionBeginIdx = thisSplitStop;
		GeneField* tmp = child1;
		child1 = child2;
		child2 = tmp;
	}
	child1 = child1Backup;
	child2 = child2Backup;
	for (unsigned int i = halfChromosomeLength; i < halfChromosomeLength +  nPointCrossoverColsOffset; i++)
	{
		child1[i] = parent1[i];
		child2[i] = parent2[i];
	}
	currentPartitionBeginIdx = halfChromosomeLength + nPointCrossoverColsOffset;
	while (currentPartitionBeginIdx < chromosomeLength)
	{
		//unsigned int thisSplitStop = currentPartitionBeginIdx + nPointCrossoverBaseSplitLength;
		//for (unsigned int i = currentPartitionBeginIdx; i < thisSplitStop; i++)
		//{
		//	child1[i] = parent1[i];
		//	child2[i] = parent2[i];
		//}
		unsigned int thisSplitStop = currentPartitionBeginIdx + nPointCrossoverBaseSplitLength;
		thisSplitStop = (thisSplitStop < chromosomeLength ? thisSplitStop :chromosomeLength);
		unsigned int i = currentPartitionBeginIdx;
		while (i < thisSplitStop)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
			i++;
		}

		currentPartitionBeginIdx = thisSplitStop;
		GeneField* tmp = child1;
		child1 = child2;
		child2 = tmp;
	}
}
template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::nPointUniformMixedCrossover(GeneField* parent1, GeneField* parent2, GeneField* child1, GeneField* child2)
{
	crossoverNPoint(parent1, parent2, child1, child2);
	std::uniform_real_distribution<double> uniformCrossoverDistr = std::uniform_real_distribution<double>(0.0, 1.0);
	for (unsigned int i = 0u; i < chromosomeLength; i++)
	{
		if (uniformCrossoverDistr(randomEngine) <= mixedCrossoverNPointUniformThreshold)
		{
			GeneField tmp = child1[i];
			child1[i] = child2[i];
			child2[i] = tmp;
		}
	}

}


template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByInversion(GeneField* specimen)
{
	//inversion mutation
	std::uniform_int_distribution<unsigned int> mutationChooser = std::uniform_int_distribution<unsigned int>(0, chromosomeLength - inversionMutationGenesCount - 1);
	unsigned int mutationFirstGeneIdx = mutationChooser(randomEngine);
	unsigned int mutationLastGeneIdx = mutationFirstGeneIdx + inversionMutationGenesCount;
	unsigned int i = mutationFirstGeneIdx, j = mutationLastGeneIdx;
	while (i < j)
	{
		GeneField tmp = specimen[i];
		specimen[i] = specimen[j];
		specimen[j] = tmp;
		i++;
		j--;
	}
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::alteringValueAndInversionMixedMutation(GeneField* specimen)
{
	std::uniform_real_distribution<double> mutationTypeChoice(0.0, 1.0);
	if (mutationTypeChoice(randomEngine) < alteringInversionMixedMutationAlteringPrefference)
	{
		mutationByAlteringGeneValue(specimen);
	}
	else
	{
		mutationByInversion(specimen);
	}

}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::findParamsForAlteringMutation(unsigned int runsCount)
{
	std::normal_distribution<double> distrNormal = std::normal_distribution<double>(0.0, static_cast<double>(matrixSize));
	double smallest = 0.0, biggest = 0.0;
	for (unsigned int i = 0; i < runsCount; i++)
	{
		double randNum = distrNormal(randomEngine);
		smallest = std::min(smallest, randNum);
		biggest = std::max(biggest, randNum);
	}
	double bigFromSmallest = abs(smallest);
	alteringMutationMaxAcceptedValue = std::max(bigFromSmallest, smallest);
	oneByAlteringMutationUnitSpan = static_cast<double>(getGenesCount()) / alteringMutationMaxAcceptedValue;

}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByAlteringGeneValue(GeneField* specimen)
{
	std::uniform_int_distribution<unsigned int> mutationChooser = std::uniform_int_distribution<unsigned int>(0, chromosomeLength - 1);
	unsigned int mutationGene = mutationChooser(randomEngine);
	std::normal_distribution<double> distrNormal = std::normal_distribution<double>(0.0, static_cast<double>(matrixSize));
	//GeneField alterValue = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero();
	//bool mutationNotPerformed = true;
	while(true)
	{ 
		
		double valFromNormalDistr = distrNormal(randomEngine);
		valFromNormalDistr = std::min(std::max(-alteringMutationMaxAcceptedValue, valFromNormalDistr), alteringMutationMaxAcceptedValue);
		long long int valFromNormalDistrMappedToGenesSpan = static_cast<long long int>(round(valFromNormalDistr * oneByAlteringMutationUnitSpan));
		long long int geneToMutateValWithExtraSpace = static_cast<long long int>(specimen[mutationGene]);

		
		long long int newGene = geneToMutateValWithExtraSpace + valFromNormalDistrMappedToGenesSpan;
		if (newGene != geneToMutateValWithExtraSpace && newGene >= 0 && newGene < static_cast<long long int>(getGenesCount()))
		{
			specimen[mutationGene] = static_cast<GeneField>(newGene);
			//mutationNotPerformed = false;
			break;
		}
	}

}

template <typename GeneField, typename MatrixField> constexpr GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero()
{
	return static_cast<GeneField>(0u);
}

template <typename GeneField, typename MatrixField> constexpr GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne()
{
	return static_cast<GeneField>(1u);
}

template <typename GeneField, typename MatrixField>  constexpr GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo()
{
	return static_cast<GeneField>(2u);
}

template <typename GeneField, typename MatrixField>  constexpr MatrixField MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero()
{
	return static_cast<MatrixField>(0u);
}

template <typename GeneField, typename MatrixField> GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getDecodingArrayRealSize()
{
	return decodingArrayRealSize;
}

template <typename GeneField, typename MatrixField> GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getDecodingArrayVisibleSize()
{
	return decodingArrayRealSize >> MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne();
}

template <typename GeneField, typename MatrixField> GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::minimalColumnIndex()
{
	return decodingArrayRealSize >> MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo();
}

template <typename GeneField, typename MatrixField> unsigned int* MatrixGeneticAlgorithm<GeneField, MatrixField>::decodeWithArray(unsigned int* returnValueArray, GeneField geneValue)
{
	GeneField geneFieldFirst = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo()* geneValue;
	GeneField geneFieldSecond = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne() + geneFieldFirst;
	returnValueArray[0] = decodingArray[geneFieldFirst];
	returnValueArray[1] = decodingArray[geneFieldSecond];
	return returnValueArray;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::setDecodingArrayValue(GeneField geneValue, unsigned int from, unsigned int to)
{
	//GeneField geneFieldFirst = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo()* geneValue;
	//GeneField geneFieldSecond = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne() + geneFieldFirst;
	decodingArray[MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo() * geneValue] = from;
	decodingArray[MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo() * geneValue + MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne()] = to;
}

template <typename GeneField, typename MatrixField> MatrixGeneticAlgorithm<GeneField, MatrixField>::MatrixGeneticAlgorithm(
	MatrixField** staringMatrix, unsigned int startingMatrixSize, double populationMatrSizePower,
	double populationSizeMultiplier, double chromesomeLenMatrSizePower, double chromosomeSizeMultplier, 
	std::default_random_engine& randomEngine, double crossoverProbability, double mutationProbability, 
	unsigned int tournamentGroupSize, unsigned int nPointCrossoverSplitsCount, double mixedCrossoverNPointUniformThreshold, double alteringInversionMixedMutationAlteringPrefference) :
	//geneFieldOne(static_cast<GeneField>(1u)),
	//geneFieldTwo(static_cast<GeneField>(2u)),
	//geneFieldZero(static_cast<GeneField>(0u)),
	matrixSize(startingMatrixSize),
	chromosomeLength(static_cast<unsigned int>(ceil(pow(startingMatrixSize, chromesomeLenMatrSizePower)* static_cast<double>(chromosomeSizeMultplier)))),
	halfChromosomeLength(static_cast<unsigned int>(ceil(pow(startingMatrixSize, chromesomeLenMatrSizePower)* static_cast<double>(chromosomeSizeMultplier)))/2u),
	populationSize(static_cast<unsigned int>(ceil(pow(startingMatrixSize, populationMatrSizePower) * populationSizeMultiplier))),
	decodingArrayRealSize(static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())),
	randomEngine(randomEngine),
	distribution( std::uniform_int_distribution<GeneField>(MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero(), ( static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne()) / MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo() ) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())),
	inversionMutationGenesCount( static_cast<unsigned int>(std::max(floor(cbrt(static_cast<double>(startingMatrixSize)) * log2(chromesomeLenMatrSizePower)), 2.0)) ),
	crossoverProbability(crossoverProbability),
	mutationProbability(mutationProbability),
	tournamentGroupSize(tournamentGroupSize),
	genesCount((static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())) / MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo() ),
	nPointCrossoverTotalSplitsCount(nPointCrossoverSplitsCount),
	mixedCrossoverNPointUniformThreshold(mixedCrossoverNPointUniformThreshold),
	alteringInversionMixedMutationAlteringPrefference(alteringInversionMixedMutationAlteringPrefference)
{
	//constructor body start
	//std::cout << "Constructor starts working, chromosome length: "<<chromosomeLength << std::endl;
	if (populationSize % 2u == 1u)
	{
		populationSize++;
	}
	bestSolution = new GeneField[chromosomeLength];
	parentPopulation = new GeneField*[populationSize];
	childrenPopulation = new GeneField * [populationSize];
	parentPopulationSelectionHelperArray = new GeneField * [populationSize];
	childrenPopulationSelectionHelperArray = new GeneField * [populationSize];
	parentPopulationTargetFunVals = new MatrixField[populationSize];
	childrenPopulationTargetFunVals = new MatrixField[populationSize];
	parentPopulationSelectionTargetFunVals = new MatrixField[populationSize];
	childrenPopulationSelectionTargetFunVals = new MatrixField[populationSize];
	selectionIsInNewParentGenerationArr = new bool[populationSize * 2u];
	indexesOfTournamentConteandants = new unsigned int[tournamentGroupSize];
	populationsTournamentConteandants = new unsigned int[tournamentGroupSize];
	decodingArray = new unsigned int[getDecodingArrayRealSize()];
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		parentPopulation[i] = new GeneField[chromosomeLength];
		childrenPopulation[i] = new GeneField[chromosomeLength];
		
		selectionIsInNewParentGenerationArr[i] = false;
		selectionIsInNewParentGenerationArr[i + populationSize] = false;
	}
	findParamsForNPointCrossover();
	findParamsForAlteringMutation();
	//std::cout << "constructor: Arrays allocated" << std::endl;
	initOptimizedMatrix(staringMatrix);
	
	//std::cout << "Before decoding array filling" << std::endl;
	fillDecodingArray();
	//std::cout << "After decoding array filling" << std::endl;
	generateStartingPopulation(populationMatrSizePower, populationSizeMultiplier);
	
	crossoverOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverNPoint;
	mutationOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByAlteringGeneValue;
	//std::cout << "Constructor finished work" << std::endl;
}



template <typename GeneField, typename MatrixField> MatrixGeneticAlgorithm<GeneField, MatrixField>::MatrixGeneticAlgorithm(
	MatrixField** staringMatrix, unsigned int startingMatrixSize, double populationMatrSizePower,
	double populationSizeMultiplier, double chromesomeLenMatrSizePower, double chromosomeSizeMultplier,
	std::default_random_engine& randomEngine, double crossoverProbability, double mutationProbability,
	unsigned int tournamentGroupSize, unsigned int nPointCrossoverSplitsCount, double mixedCrossoverNPointUniformThreshold,
	double alteringInversionMixedMutationAlteringPrefference, unsigned int inversionMutationGenesCount) :
	//geneFieldOne(static_cast<GeneField>(1u)),
	//geneFieldTwo(static_cast<GeneField>(2u)),
	//geneFieldZero(static_cast<GeneField>(0u)),
	matrixSize(startingMatrixSize),
	chromosomeLength(static_cast<unsigned int>(ceil(pow(startingMatrixSize, chromesomeLenMatrSizePower)* static_cast<double>(chromosomeSizeMultplier)))),
	halfChromosomeLength(static_cast<unsigned int>(ceil(pow(startingMatrixSize, chromesomeLenMatrSizePower)* static_cast<double>(chromosomeSizeMultplier))) / 2u),
	populationSize(static_cast<unsigned int>(ceil(pow(startingMatrixSize, populationMatrSizePower)* populationSizeMultiplier))),
	decodingArrayRealSize(static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())),
	randomEngine(randomEngine),
	distribution(std::uniform_int_distribution<GeneField>(MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero(), (static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne()) / MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo()) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())),
	inversionMutationGenesCount(inversionMutationGenesCount),
	crossoverProbability(crossoverProbability),
	mutationProbability(mutationProbability),
	tournamentGroupSize(tournamentGroupSize),
	genesCount((static_cast<GeneField>(startingMatrixSize)* (static_cast<GeneField>(startingMatrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne())) / MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo()),
	nPointCrossoverTotalSplitsCount(nPointCrossoverSplitsCount),
	mixedCrossoverNPointUniformThreshold(mixedCrossoverNPointUniformThreshold),
	alteringInversionMixedMutationAlteringPrefference(alteringInversionMixedMutationAlteringPrefference)
{
	//constructor body start
	//std::cout << "Constructor starts working, chromosome length: "<<chromosomeLength << std::endl;
	if (populationSize % 2u == 1u)
	{
		populationSize++;
	}
	bestSolution = new GeneField[chromosomeLength];
	parentPopulation = new GeneField * [populationSize];
	childrenPopulation = new GeneField * [populationSize];
	parentPopulationSelectionHelperArray = new GeneField * [populationSize];
	childrenPopulationSelectionHelperArray = new GeneField * [populationSize];
	parentPopulationTargetFunVals = new MatrixField[populationSize];
	childrenPopulationTargetFunVals = new MatrixField[populationSize];
	parentPopulationSelectionTargetFunVals = new MatrixField[populationSize];
	childrenPopulationSelectionTargetFunVals = new MatrixField[populationSize];
	selectionIsInNewParentGenerationArr = new bool[populationSize * 2u];
	indexesOfTournamentConteandants = new unsigned int[tournamentGroupSize];
	populationsTournamentConteandants = new unsigned int[tournamentGroupSize];
	decodingArray = new unsigned int[getDecodingArrayRealSize()];
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		parentPopulation[i] = new GeneField[chromosomeLength];
		childrenPopulation[i] = new GeneField[chromosomeLength];

		selectionIsInNewParentGenerationArr[i] = false;
		selectionIsInNewParentGenerationArr[i + populationSize] = false;
	}
	findParamsForNPointCrossover();
	findParamsForAlteringMutation();
	//std::cout << "constructor: Arrays allocated" << std::endl;
	initOptimizedMatrix(staringMatrix);

	//std::cout << "Before decoding array filling" << std::endl;
	fillDecodingArray();
	//std::cout << "After decoding array filling" << std::endl;
	generateStartingPopulation(populationMatrSizePower, populationSizeMultiplier);

	crossoverOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::crossoverNPoint;
	mutationOperationFunctionPointer = &MatrixGeneticAlgorithm<GeneField, MatrixField>::mutationByAlteringGeneValue;
	//std::cout << "Constructor finished work" << std::endl;
}



template <typename GeneField, typename MatrixField> void  MatrixGeneticAlgorithm<GeneField, MatrixField>::copyMatrixToTargetFunctionArray(MatrixField** source)
{
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		for (unsigned int j = 0u; j < matrixSize; j++)
		{
			matrixCopyForFitnessEvaluation[i + 1u][j + 1u] = source[i][j];
		}
	}
}

template <typename GeneField, typename MatrixField> void  MatrixGeneticAlgorithm<GeneField, MatrixField>::decodeMatrix(GeneField* solution)
{
	unsigned int swapArray[2];
	//GeneField minColumnSwapIdx = minimalColumnIndex();
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		matrixCopyHashingRowsArray[i] = i;
		matrixCopyHashingColsArray[i] = i;
	}

	for (unsigned int i = 0u; i < halfChromosomeLength; i++)
	{
		decodeWithArray(swapArray, solution[i]);
		unsigned int tmp = matrixCopyHashingRowsArray[swapArray[0u]];
		matrixCopyHashingRowsArray[swapArray[0u]] = matrixCopyHashingRowsArray[swapArray[1u]];
		matrixCopyHashingRowsArray[swapArray[1u]] = tmp;
	}
	for (unsigned int i = halfChromosomeLength; i < chromosomeLength; i++)
	{
		decodeWithArray(swapArray, solution[i]);
		unsigned int tmp = matrixCopyHashingColsArray[swapArray[0u]];
		matrixCopyHashingColsArray[swapArray[0u]] = matrixCopyHashingColsArray[swapArray[1u]];
		matrixCopyHashingColsArray[swapArray[1u]] = tmp;
	}

}

template <typename GeneField, typename MatrixField> double MatrixGeneticAlgorithm<GeneField, MatrixField>::getCurrentBestSolutionTargetFunctionValue()
{
	return static_cast<double>(bestSolutionTargetFunctionValue) * 0.5;
}

template <typename GeneField, typename MatrixField> MatrixField MatrixGeneticAlgorithm<GeneField, MatrixField>::getCurrentBestSolutionTargetFunctionValueNoHalf()
{
	return bestSolutionTargetFunctionValue;
}

template <typename GeneField, typename MatrixField> double MatrixGeneticAlgorithm<GeneField, MatrixField>::targetFunctionValue(GeneField* solution)
{
	return static_cast<double>(targetFunctionValueNoHalf(solution)) * 0.5;
}

template <typename GeneField, typename MatrixField> double MatrixGeneticAlgorithm<GeneField, MatrixField>::getTargetFunctionValueForPassedMatrix(MatrixField** passedMatrix)
{

	copyMatrixToTargetFunctionArray(passedMatrix);
	unsigned int matrixSizeWithOffset = matrixSize + 1u;
	MatrixField targetFunValue = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
	for (unsigned int i = 1u; i < matrixSizeWithOffset; i++)
	{
		for (unsigned int j = 1u; j < matrixSizeWithOffset; j++)
		{
			targetFunValue += (matrixCopyForFitnessEvaluation[i][j] * (matrixCopyForFitnessEvaluation[i - 1][j] + matrixCopyForFitnessEvaluation[i + 1][j] +
				matrixCopyForFitnessEvaluation[i][j - 1] + matrixCopyForFitnessEvaluation[i][j + 1]));
		}
	}
	copyMatrixToTargetFunctionArray(matrixBeforeOptimization);
	return static_cast<double>(targetFunValue)*0.5;
}

template <typename GeneField, typename MatrixField> MatrixField MatrixGeneticAlgorithm<GeneField, MatrixField>::targetFunctionValueNoHalf(GeneField* solution)
{
	//copyMatrixToTargetFunctionArray(matrixBeforeOptimization);
	//std::cout << "targetFunctionValueNoHalf: after copying matrix to target function array" << std::endl;
	decodeMatrix(solution);
	//std::cout << "targetFunctionValueNoHalf: after decoding with decoding array" << std::endl;
	unsigned int matrixSizeWithOffset = matrixSize + 1u;
	MatrixField targetFunValue = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
	//std::cout << "targetFunctionValueNoHalf: before calculating target function value" << std::endl;
	
	
	for (int i = 0; i < static_cast<int>(matrixSize); i++)
	{
		for (int j = 0; j < static_cast<int>(matrixSize); j++)
		{
			targetFunValue += ( matrixCopyForFitnessEvaluation[ matrixCopyHashingRowsArray[i] + 1][ matrixCopyHashingColsArray[j] +1 ] * (
				matrixCopyForFitnessEvaluation[matrixCopyHashingRowsArray[i - 1] + 1][matrixCopyHashingColsArray[j] + 1] + 
				matrixCopyForFitnessEvaluation[matrixCopyHashingRowsArray[i + 1] + 1][matrixCopyHashingColsArray[j] + 1] + 
				matrixCopyForFitnessEvaluation[matrixCopyHashingRowsArray[i] + 1][matrixCopyHashingColsArray[j - 1] + 1] + 
				matrixCopyForFitnessEvaluation[matrixCopyHashingRowsArray[i] + 1][matrixCopyHashingColsArray[j + 1] + 1])  );
		}
	}

	//std::cout << "targetFunctionValueNoHalf: after calculating target function value" << std::endl;
	return targetFunValue;
}

template <typename GeneField, typename MatrixField> void  MatrixGeneticAlgorithm<GeneField, MatrixField>::initOptimizedMatrix(MatrixField** source)
{
	
	matrixBeforeOptimization = new MatrixField*[matrixSize];
	matrixCopyForFitnessEvaluation = new MatrixField*[matrixSize + 2u];
	matrixCopyHashingColsArrayWithOffsets = new int[matrixSize + 2u];
	matrixCopyHashingRowsArrayWithOffsets = new int[matrixSize + 2u];
	matrixCopyHashingColsArray = matrixCopyHashingColsArrayWithOffsets + 1;
	matrixCopyHashingRowsArray = matrixCopyHashingRowsArrayWithOffsets + 1;
	matrixCopyHashingColsArrayWithOffsets[0u] = -1;
	matrixCopyHashingColsArrayWithOffsets[matrixSize + 1u] = static_cast<int>(matrixSize);
	matrixCopyHashingRowsArrayWithOffsets[0u] = -1;
	matrixCopyHashingRowsArrayWithOffsets[matrixSize + 1u] = static_cast<int>(matrixSize);

	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		matrixCopyHashingColsArray[i] = i;
		matrixCopyHashingRowsArray[i] = i;
		matrixBeforeOptimization[i] = new MatrixField[matrixSize];
		for (unsigned int j = 0u; j < matrixSize; j++)
		{
			matrixBeforeOptimization[i][j] = source[i][j];
		}
	}
	for (unsigned int i = 0u; i < matrixSize + 2u; i++)
	{
		matrixCopyForFitnessEvaluation[i] = new MatrixField[matrixSize + 2u];
	}
	for (unsigned int i = 0u; i < matrixSize + 2u; i++)
	{
		matrixCopyForFitnessEvaluation[0u][i] = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
		matrixCopyForFitnessEvaluation[i][0u] = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
		matrixCopyForFitnessEvaluation[matrixSize + 1u][i] = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
		matrixCopyForFitnessEvaluation[i][matrixSize + 1u] = MatrixGeneticAlgorithm<GeneField, MatrixField>::getMatrixFieldZero();
	}
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		for (unsigned int j = 0u; j < matrixSize; j++)
		{
			matrixCopyForFitnessEvaluation[i + 1u][j + 1u] = source[i][j];
		}
	}
}

template <typename GeneField, typename MatrixField> MatrixGeneticAlgorithm<GeneField, MatrixField>::~MatrixGeneticAlgorithm()
{
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		delete[] matrixBeforeOptimization[i];
	}
	for (unsigned int i = 0u; i < matrixSize + 2u; i++)
	{
		delete[] matrixCopyForFitnessEvaluation[i];
	}
	delete[] bestSolution;
	delete[] matrixBeforeOptimization;
	delete[] matrixCopyForFitnessEvaluation;
	delete[] matrixCopyHashingColsArrayWithOffsets;
	delete[] matrixCopyHashingRowsArrayWithOffsets;
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		delete[] parentPopulation[i];
		delete[] childrenPopulation[i];
	}
	delete[] parentPopulation;
	delete[] childrenPopulation;
	delete[] decodingArray;

	delete[] parentPopulationSelectionHelperArray;
	delete[] childrenPopulationSelectionHelperArray;
	delete[] parentPopulationTargetFunVals;
	delete[] childrenPopulationTargetFunVals;
	delete[] parentPopulationSelectionTargetFunVals;
	delete[] childrenPopulationSelectionTargetFunVals;
	delete[] selectionIsInNewParentGenerationArr;
	delete[] indexesOfTournamentConteandants;
	delete[] populationsTournamentConteandants;

}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::fillDecodingArray()
{
	GeneField idx = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero();
	GeneField one = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne();
	GeneField two = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldTwo();
	//GeneField columsOffset = static_cast<GeneField>(matrixSize) * (static_cast<GeneField>(matrixSize) - MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldOne() );
	//GeneField columsOffset = minimalColumnIndex();
	for (unsigned int i = 0u; i < matrixSize - 1u; i++)
	{
		for (unsigned int j = i + 1u; j < matrixSize; j++)
		{
			
			decodingArray[idx] = i;
			decodingArray[idx + one] = j;
			//decodingArray[idx + columsOffset] = i;
			//decodingArray[idx + columsOffset + one] = j;
			idx += two;
			
			//setDecodingArrayValue(idx, i, j);
			//idx++;
		}
	}
	
	//std::cout << "fillDecodingArray: Decoding matrix: " << std::endl;
	//testingCode
	//for (GeneField i = static_cast<GeneField>(0u); i < getDecodingArrayRealSize() ; i+=two)
	//{
	//	std::cout << decodingArray[i] << " -> " << decodingArray[i + one] << std::endl;
	//}
	//std::cout << "fillDecodingArray: Decoding array printed" << std::endl;
	//testingCodeDone
}

template <typename GeneField, typename MatrixField> GeneField MatrixGeneticAlgorithm<GeneField, MatrixField>::getGenesCount()
{
	return genesCount;
}

template <typename GeneField, typename MatrixField> unsigned int MatrixGeneticAlgorithm<GeneField, MatrixField>::findBestSolutionIndexInArray(MatrixField* targetFunctionValues, unsigned int targetFunValsArrayLength)
{
	MatrixField bestValue = targetFunctionValues[0];
	unsigned int bestValueIdx = 0u;
	for (unsigned int i = 0u; i < targetFunValsArrayLength; i++)
	{
		if (targetFunctionValues[i] > bestValue)
		{
			bestValue = targetFunctionValues[i];
			bestValueIdx = i;
		}
	}
	return bestValueIdx;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::generateStartingPopulation(double populationPower, double populationMultplier)
{
	//fill every specimen with random genes
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		for (unsigned int j = 0u; j < chromosomeLength; j++)
		{
			parentPopulation[i][j] = distribution(randomEngine);
			//testingCode
			//if (parentPopulation[i][j] >= getGenesCount())
			//{
			//	std::cout << "\n\n\n\n\n\n\ngenerateStartingPopulation: too big gene!!! " << parentPopulation[i][j] << ", genesCount: " << genesCount <<"\n\n\n\n\n\n" << std::endl;
			//}
			//testingCodeDone
		}
	}
	//make sure that every gene is present in some specimen
	std::uniform_int_distribution<unsigned int> distributionPopulation = std::uniform_int_distribution<unsigned int>(0u, populationSize - 1);
	std::uniform_int_distribution<unsigned int> distributionChromosome = std::uniform_int_distribution<unsigned int>(0u, chromosomeLength - 1);
	//std::cout << "generate starting population: before loop with adding fixed amount of each gene for population" << std::endl;
	unsigned int singleGeneRepetitionsCount = static_cast<unsigned int>(floor(std::max(static_cast<double>(matrixSize), pow(static_cast<double>(matrixSize), std::max(cbrt(populationPower), 1.0)) * std::max(cbrt(populationMultplier), 1.0))));
	for (GeneField i = MatrixGeneticAlgorithm<GeneField, MatrixField>::getGeneFieldZero(); i < getGenesCount(); i++)
	{
		for (unsigned int j = 0; j < singleGeneRepetitionsCount; j++)
		{
			unsigned int specimenNo = distributionPopulation(randomEngine);
			unsigned int geneNo = distributionChromosome(randomEngine);
			parentPopulation[specimenNo][geneNo] = i;
		}
	}
	//std::cout << "generate starting population: after loop with adding fixed amount of each gene for population" << std::endl;
	//fill target function array with fitness values for starting population and find best solution
	//std::cout << "generate starting population: before calculating target fun vals" << std::endl;
	for (unsigned int i = 0; i < populationSize; i++)
	{
		parentPopulationTargetFunVals[i] = targetFunctionValueNoHalf(parentPopulation[i]);
	}
	//std::cout << "generate starting population: after calculating target fun vals" << std::endl;
	unsigned int bestSolutionIdx = findBestSolutionIndexInArray(parentPopulationTargetFunVals, populationSize);
	bestSolutionTargetFunctionValue = parentPopulationTargetFunVals[bestSolutionIdx];
	for (unsigned int i = 0; i < chromosomeLength; i++)
	{
		bestSolution[i] = parentPopulation[bestSolutionIdx][i];
	}
	//std::cout << "Starting population generated" << std::endl;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::generateChildrenPopulation()
{
	//std::cout << "generateChildrenPopulation: start generating population" << std::endl;
	std::uniform_real_distribution<double> geneticOperatorsProbabilityGenerator = std::uniform_real_distribution<double>(0.0, 1.0);
	std::uniform_int_distribution<unsigned int> parentsChoice = std::uniform_int_distribution<unsigned int>(0, populationSize - 1u);
	unsigned int pairsCount = populationSize >> 1u;
	unsigned int newChildrenCount = 0u;
	//go through all pairs of specimens and TRY to make them reproduce
	//shuffling of the parent population done in the selection stage
	for (unsigned int i = 0u; i < pairsCount; i += 2u)
	{
		if (geneticOperatorsProbabilityGenerator(randomEngine) <= crossoverProbability)
		{
			//crossoverUniform(parentPopulation[i], parentPopulation[i + 1u], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			//crossoverNPoint(parentPopulation[i], parentPopulation[i + 1u], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			(this->*crossoverOperationFunctionPointer)(parentPopulation[i], parentPopulation[i + 1u], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			newChildrenCount += 2u;
		}
	}
	//if some pair didnt reproduce fill missing specimens in children population with children of random specimens in parent population
	while (newChildrenCount < populationSize)
	{
		unsigned int parentOneIdx = parentsChoice(randomEngine);
		unsigned int parentTwoIdx = parentsChoice(randomEngine);
		if (parentOneIdx != parentTwoIdx)
		{
			//crossoverUniform(parentPopulation[parentOneIdx], parentPopulation[parentTwoIdx], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			//crossoverNPoint(parentPopulation[parentOneIdx], parentPopulation[parentTwoIdx], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			(this->*crossoverOperationFunctionPointer)(parentPopulation[parentOneIdx], parentPopulation[parentTwoIdx], childrenPopulation[newChildrenCount], childrenPopulation[newChildrenCount + 1u]);
			newChildrenCount += 2u;
		}
	}
	//try mutating each child specimen and find their target function values
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		if (geneticOperatorsProbabilityGenerator(randomEngine) <= mutationProbability)
		{
			//mutationByInversion(childrenPopulation[i]);
			//mutationByAlteringGeneValue(childrenPopulation[i]);
			(this->*mutationOperationFunctionPointer)(childrenPopulation[i]);
		}
		childrenPopulationTargetFunVals[i] = targetFunctionValueNoHalf(childrenPopulation[i]);
		//check if child is better than best solution
		if (childrenPopulationTargetFunVals[i] > bestSolutionTargetFunctionValue)
		{
			GeneField* newBestChildInChildrenPopulation = childrenPopulation[i];
			bestSolutionTargetFunctionValue = childrenPopulationTargetFunVals[i];
			for (unsigned int j = 0; j < chromosomeLength; j++)
			{
				bestSolution[j] = newBestChildInChildrenPopulation[j];
			}
		}
	}
	//std::cout << "generateChildrenPopulation: finish generating population" << std::endl;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::tournamentSelectionMiPlusLambda()
{
	//std::cout << "tournamentSelectionMiPlusLambda: begin selection" << std::endl;
	//for (unsigned int i = 0u; i < populationSize * 2u; i++)
	//{
	//	selectionIsInNewParentGenerationArr[i] = false;
	//}
	//reseting flags array not necessary because flags are reseted during moving of unselected specimens to helper array
	GeneField** parentChildrenPopulationsArr[] = { parentPopulation, childrenPopulation };
	MatrixField* parentChildrenTargetFunValsArr[] = { parentPopulationTargetFunVals, childrenPopulationTargetFunVals };
	unsigned int totalSelectedSpecimens = 0u;
	std::uniform_int_distribution<unsigned int> conteandantSelectioner = std::uniform_int_distribution<unsigned int>(0, populationSize - 1);
	std::uniform_int_distribution<unsigned int> populationSelectioner = std::uniform_int_distribution<unsigned int>(0, 1);
	while (totalSelectedSpecimens < populationSize)
	{
		//std::cout << "*";
		//find candidates
		unsigned int selectedForTournamentCount = 0u;
		while (selectedForTournamentCount < tournamentGroupSize)
		{
			unsigned int possibleConteandantPopulationIdx = populationSelectioner(randomEngine);
			unsigned int possibleConteandantIdx = conteandantSelectioner(randomEngine);
			if (!selectionIsInNewParentGenerationArr[possibleConteandantPopulationIdx * populationSize + possibleConteandantIdx])
			{
				selectionIsInNewParentGenerationArr[possibleConteandantPopulationIdx * populationSize + possibleConteandantIdx] = true;
				populationsTournamentConteandants[selectedForTournamentCount] = possibleConteandantPopulationIdx;
				indexesOfTournamentConteandants[selectedForTournamentCount] = possibleConteandantIdx;
				selectedForTournamentCount++;
			}
 		}
		//choose best conteandant
		MatrixField bestValue = parentChildrenTargetFunValsArr[populationsTournamentConteandants[0]][indexesOfTournamentConteandants[0]];
		unsigned int bestConteandantPopulationIdx = populationsTournamentConteandants[0];
		unsigned int bestConteandantIdx = indexesOfTournamentConteandants[0];
		for (unsigned int i = 0u; i < tournamentGroupSize; i++)
		{
			if (parentChildrenTargetFunValsArr[populationsTournamentConteandants[i]][indexesOfTournamentConteandants[i]] > bestValue)
			{
				bestValue = parentChildrenTargetFunValsArr[populationsTournamentConteandants[i]][indexesOfTournamentConteandants[i]];
				bestConteandantPopulationIdx = populationsTournamentConteandants[i];
				bestConteandantIdx = indexesOfTournamentConteandants[i];
			}
		}
		//fix the flags array - set selected flags of conteandants to false and then set the winner flag to true
		for (unsigned int i = 0u; i < tournamentGroupSize; i++)
		{
			selectionIsInNewParentGenerationArr[populationsTournamentConteandants[i] * populationSize + indexesOfTournamentConteandants[i] ] = false;
		}
		
		selectionIsInNewParentGenerationArr[bestConteandantPopulationIdx * populationSize + bestConteandantIdx] = true;

		//add tournament winner to new population parents helper array
		parentPopulationSelectionHelperArray[totalSelectedSpecimens] = parentChildrenPopulationsArr[bestConteandantPopulationIdx][bestConteandantIdx];
		parentPopulationSelectionTargetFunVals[totalSelectedSpecimens] = parentChildrenTargetFunValsArr[bestConteandantPopulationIdx][bestConteandantIdx];
		totalSelectedSpecimens++;
	}

	//move unselected specimens to helper children arrays
	unsigned int totalRejectedSpecimensMoved = 0u;
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		if (!selectionIsInNewParentGenerationArr[i])
		{
			childrenPopulationSelectionHelperArray[totalRejectedSpecimensMoved] = parentPopulation[i];
			//childrenPopulationSelectionTargetFunVals[totalRejectedSpecimensMoved] = parentPopulationTargetFunVals[i];
			//assignment not necessary as children from next generation have new target function assigned in this array
			totalRejectedSpecimensMoved++;
		}
		//reset flags array for next generation selection 
		selectionIsInNewParentGenerationArr[i] = false;
	}
	for (unsigned int i = 0u; i < populationSize; i++)
	{
		if (!selectionIsInNewParentGenerationArr[populationSize + i])
		{
			childrenPopulationSelectionHelperArray[totalRejectedSpecimensMoved] = childrenPopulation[i];
			//childrenPopulationSelectionTargetFunVals[totalRejectedSpecimensMoved] = childrenPopulationTargetFunVals[i];
			//assignment not necessary as children from next generation have new target function assigned in this array
			totalRejectedSpecimensMoved++;
		}
		selectionIsInNewParentGenerationArr[populationSize + i] = false;
		//reset flags array for next generation selection 
	}
	//swap current parent and children arrays with helper arrays
	GeneField** tmpParentPopulation = parentPopulation;
	parentPopulation = parentPopulationSelectionHelperArray;
	parentPopulationSelectionHelperArray = tmpParentPopulation;
	MatrixField* tmpParentFunVals = parentPopulationTargetFunVals;
	parentPopulationTargetFunVals = parentPopulationSelectionTargetFunVals;
	parentPopulationSelectionTargetFunVals = tmpParentFunVals;

	GeneField** tmpChildrenPopulation = childrenPopulation;
	childrenPopulation = childrenPopulationSelectionHelperArray;
	childrenPopulationSelectionHelperArray = tmpChildrenPopulation;
	MatrixField* tmpChildrenFunVals = childrenPopulationTargetFunVals;
	childrenPopulationTargetFunVals = childrenPopulationSelectionTargetFunVals;
	childrenPopulationSelectionTargetFunVals = tmpChildrenFunVals;
	//std::cout << "tournamentSelectionMiPlusLambda: selection done" << std::endl;
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::saveSolutionPackageToFile(const std::string& solutionFolder, const std::string& matrixFileName, const std::string& targetFunctionValuesFileName)
{
	MatrixField** bestSol = getCurrentBestSolution();
	std::ofstream outStream(solutionFolder + matrixFileName);
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		for (unsigned int j = 0u; j < matrixSize; j++)
		{
			outStream << bestSol[i][j];
			if (j < matrixSize - 1u)
			{
				outStream << "\t";
			}
			
		}
		outStream << "\n";
		delete[] bestSol[i];
	}
	delete[] bestSol;
	outStream.close();
	outStream.open(solutionFolder + targetFunctionValuesFileName, std::ios::app);
	outStream << std::fixed;
	//outStream << std::dec;
	outStream << std::setprecision(5)<<  getCurrentBestSolutionTargetFunctionValue() << "\n";
	outStream.close();
}

template <typename GeneField, typename MatrixField> void MatrixGeneticAlgorithm<GeneField, MatrixField>::solveWithNGenerations(unsigned int generationsCount, bool verbose, bool savePartialResults, std::string outputFolderPathEndingWithSlash)
{
	copyMatrixToTargetFunctionArray(matrixBeforeOptimization);
	if (!verbose && !savePartialResults)
	{
		for (unsigned int i = 0; i < generationsCount; i++)
		{
			generateChildrenPopulation();
			//std::cout << "makeGenerations: children population generated" << std::endl;
			tournamentSelectionMiPlusLambda();

		}
	}
	else
	{
		unsigned int intervalsCount = 10u;
		unsigned int comunicateCounter = 0u;
		unsigned int generationsInterval = generationsCount / intervalsCount;
		unsigned int generationsReminder = generationsCount - intervalsCount * generationsInterval;


		std::ios_base::fmtflags flagsBefore{ std::cout.flags() };
		std::cout << std::fixed;
		//std::cout << std::dec;
		
		if (verbose)
		{
			std::cout << "Starting evolution process" << std::endl;
		}
		unsigned int reminder_file_idx = 0u;
		for (unsigned int i = 0u; i < intervalsCount; i++)
		{
			for (unsigned int j = 0u; j < generationsInterval; j++)
			{
				generateChildrenPopulation();
				tournamentSelectionMiPlusLambda();
			}
			if (verbose)
			{
				std::cout << i + 1 << ":\t generated " << (i + 1) * generationsInterval << "\t\t generations(" << std::setprecision(2) << (10000.0 * static_cast<double>((i + 1) * generationsInterval) / static_cast<double>(generationsCount)) * 0.01 << "%) target function value: \t" << std::setprecision(5) << getCurrentBestSolutionTargetFunctionValue() << std::endl;
			}
			if (savePartialResults)
			{
				saveSolutionPackageToFile(outputFolderPathEndingWithSlash, std::string("solution_" + std::to_string(i) + ".txt"), "target_fun_vals_for_intermediate_solutions.txt");
			}
			reminder_file_idx++;
		}
		for (unsigned int i = 0u; i <generationsReminder; i++)
		{
			generateChildrenPopulation();
			tournamentSelectionMiPlusLambda();
		}
		if (generationsReminder > 0u && verbose)
		{
			std::cout << intervalsCount + 1 << ":\t generated " << intervalsCount * generationsInterval + generationsReminder<< "\t\t generations(" << std::setprecision(2) << (10000.0 * static_cast<double>(intervalsCount * generationsInterval + generationsReminder) / static_cast<double>(generationsCount)) * 0.01 << "%) target function value: \t" << std::setprecision(5) << getCurrentBestSolutionTargetFunctionValue() << std::endl;
			
		}
		if (generationsReminder > 0u && savePartialResults)
		{
			saveSolutionPackageToFile(outputFolderPathEndingWithSlash, std::string("solution_" + std::to_string(reminder_file_idx) + ".txt"), "target_fun_vals_for_intermediate_solutions.txt");
		}
		if (verbose)
		{
			std::cout << "Evolution process done" << std::endl;
		}
		std::cout.setf(flagsBefore, std::ios::floatfield);
	}
}

template <typename GeneField, typename MatrixField> MatrixField** MatrixGeneticAlgorithm<GeneField, MatrixField>::getCurrentBestSolution()
{
	decodeMatrix(bestSolution);
	MatrixField** decodedMatrix = new MatrixField * [matrixSize];
	for (unsigned int i = 0u; i < matrixSize; i++)
	{
		decodedMatrix[i] = new MatrixField[matrixSize];
		for (unsigned int j = 0u; j < matrixSize; j++)
		{
			decodedMatrix[i][j] = matrixCopyForFitnessEvaluation[ matrixCopyHashingRowsArray[i] + 1u][matrixCopyHashingColsArray[j] + 1u];
		}
	}
	return decodedMatrix;
}

;


/*
GeneField** parentPopulationSelectionHelperArray;
	MatrixField* parentPopulationSelectionTargetFunVals;
	GeneField** childrenPopulationSelectionHelperArray;
	MatrixField* childrenPopulationTargetFunVals;
	bool* selectionIsInNewParentGenerationArr;
	unsigned int tournamentGroupSize;
	unsigned int* indexesOfTournamentConteandants;

*/

