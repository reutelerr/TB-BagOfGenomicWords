import json
import os.path
import time
import sys
import re

from matplotlib import pyplot
import numpy as np

import kmer_Counting
import utils
import modelTrainingAndTesting
import random
import statistics
from src.constants import *
from src.testSequenceInjecter import needlemanWunschInjectedSequence, injectSequence


def displayHelp():
    print(
        "AVAILIBLE COMMANDS\n\n"
        " dictionary <action> <dictionaryName> <sourceFilepath>\n"
        "   - <action> values : 'new', 'add', 'filter'\n\n"
        " makeBOWs <dictionaryName> <sourceFilepath> <outputFilePath>\n\n"
        " convertCSVtoFASTA <csvFilePath> <fastaFilePath>\n\n"
        " findAlignments <sequenceFilePath>\n\n"
        " splitTrainingAndTesting <sequenceFilePath> <labelFilePath> <trainingFilePath> <trainingLabelPath> <testingFilePath> <testingLabelPath>\n\n"
        " injectSequence <seqFilePath> <newSeqFilePath> <labelFilePath>\n\n"
        " train <modelType> <trainingBOWsFilePath> <trainingLabelsFilePath> <testingBOWsFilePath> <testingLabelsFilePath>\n\n"
    )


def ReadSequenceFile(mode, dictionaryPath, sourcePath, outputPath=''):
    fileExtension = sourcePath.split('.')[-1]
    if fileExtension == "csv":
        sourceType = CSV
    if fileExtension == "fasta":
        sourceType = FASTA
    if fileExtension == "txt":
        sourceType = FASTA

    if mode == makeBOW:
        kmer_Counting.vectorize(dictionaryPath, sourcePath, outputPath, sourceType)
    else:
        kmer_Counting.buildDictionary(mode, dictionaryPath, sourcePath, sourceType)


def Filter(dicPath):
    kmer_Counting.filterByFrequency(dicPath)


def buildPhageDictionary():
    print("Building phage dictionary")
    ReadSequenceFile(NEW, phageKmerDictionaryPath, phageSeqPath)
    print("Filtering phage dictionary")
    Filter(phageKmerDictionaryPath)


def PlotDictionaryPageSizes():
    dictionaryPath = "phageKmerDictionaries/dic10"
    maxKmerLength = metaParameters['vectorization']['maxKmerLength']
    for i in range(10, maxKmerLength):
        print("Building dictionary")
        ReadSequenceFile(NEW, dictionaryPath, injectedPhageSeqPath)
        print("Filtering dictionary")
        Filter('phageKmerDictionaries/dic10')
        # TODO


def dictionaryStats():
    print("Building dictionary")
    ReadSequenceFile(NEW, phageKmerDictionaryPath, injectedPhageSeqPath)
    maxKmerLength = metaParameters['vectorization']['maxKmerLength']
    """
    for length in range(3, maxKmerLength + 1):
        dictionaryFileName = phageKmerDictionaryPath + "/kmerCounts" + str(length) + ".json"
        dictionaryFile = open(dictionaryFileName, "rb")
        (nucleotideCount, validKmers) = json.load(dictionaryFile)
        noOfPossibleKmers = 4 ** length
        noOfActualKmers = len(validKmers)
        noOf0OccurenceKmers = noOfPossibleKmers - noOfActualKmers
        print("No of kmers of length " + str(length) + " not in dataset : " + str(noOf0OccurenceKmers))
        print("Proportion of kmers not in dataset : " + str(noOf0OccurenceKmers/noOfPossibleKmers))
        ambiguousChars = ['n', 'm', 'r', 'w', 's', 'y', 'k', 'v', 'h', 'd', 'b']
        ambiguousKmers = [ambiguousKmer for ambiguousKmer in validKmers.values() if ambiguousChars in ambiguousKmer]
        print("Number of ambiguous kmers : " + sum(ambiguousKmers))
        print("Total number of kmers : " +sum(validKmers.values()))
        pyplot.boxplot(validKmer for validKmer in list(validKmers.values()))
        pyplot.title('Frequency distribution of ' + str(length) + '-mers')
        pyplot.show()
    #for length in range(3, maxKmerLength + 1):
    """

    print("Filtering dictionary")
    Filter('phageKmerDictionaries/dic10')
    dicSizes = []
    for length in range(3, maxKmerLength + 1):
        dictionaryFileName = phageKmerDictionaryPath + "/kmerCounts" + str(length) + ".json"
        dicSize = os.path.getsize(dictionaryFileName)
        dicSizes.append(dicSize)
        print("Dictionary size for length" + str(length) + " (in bytes): " + str(dicSize))
    pyplot.bar(np.arange(3, maxKmerLength + 1), dicSizes)
    frequencyThreshold = metaParameters['vectorization']['minimumFrequencyThreshold']
    frequencyThresholdModifier = metaParameters['vectorization']['filterThresholdLengthModifier']
    pyplot.title('Filtered Dictionary Sizes, with threshold = ' + str(frequencyThreshold) + '; modifier = ' + str(
        frequencyThresholdModifier))
    pyplot.xlabel('kmerLength')
    pyplot.ylabel('Size')
    pyplot.show()


def injectAndEvaluateModel(modelType, outputFilePath, buildDictionary=False, discard=0):
    variabilityType = metaParameters['sequenceInjection']['variabilityType']
    labelsPath = injectedPhageSeqFolder + "labels.txt"
    trainingPath = injectedPhageSeqFolder + "Training/"
    trainingSeqPath = trainingPath + "seq.csv"
    trainingLabelsPath = trainingPath + "labels.txt"
    testingPath = injectedPhageSeqFolder + "Testing/"
    testingSeqPath = testingPath + "seq.csv"
    testingLabelsPath = testingPath + "labels.txt"
    trainingBOWsPath = trainingPath + "Seq_BOW.json"
    testingBOWsPath = testingPath + "Seq_BOW.json"

    if buildDictionary:
        # Probably superfluous to rebuild dictionary every time
        print("Building dictionary")
        ReadSequenceFile(NEW, phageKmerDictionaryPath, injectedPhageSeqPath)
        print("Filtering dictionary")
        Filter(phageKmerDictionaryPath)

    injectSequence(phageSeqPath, injectedPhageSeqPath, labelsPath, variabilityType=variabilityType)
    print("Splitting training and testing data")
    utils.splitTrainingAndTesting(injectedPhageSeqPath, labelsPath, trainingPath, testingPath, discard=discard)
    print("Making BOWs (training)")
    ReadSequenceFile(makeBOW, phageKmerDictionaryPath, trainingSeqPath, trainingBOWsPath)
    print("Making BOWs (testing)")
    ReadSequenceFile(makeBOW, phageKmerDictionaryPath, testingSeqPath, testingBOWsPath)
    print("Building model")
    (score, refitTime, model) = modelTrainingAndTesting.trainAndTestModel(
        modelType,
        trainingBOWsPath,
        trainingLabelsPath,
        testingBOWsPath,
        testingLabelsPath,
        randomState=randomStateSeed,
        resultsFilePath=outputFilePath,
        verbose=4
    )
    return score, refitTime, model


def testInjectionRate():
    minInjectionRate = float(sys.argv[2])

    scores = []
    refitTimes = []
    start_time = time.time()
    injectionRates = np.arange(0.5, minInjectionRate-0.1, -0.1)
    for i in injectionRates:
        print("\nINJECTION RATE : "+str(i)+"\n")
        metaParameters['sequenceInjection']['injectionRate'] = i
        (score, refitTime, model) = injectAndEvaluateModel('MLP', 'testInjectionRate.txt', True, discard=0.75)
        scores.append(score)
        refitTimes.append(refitTimes)

    print("Scores : ")
    print(scores)
    print("RefitTimes : ")
    print(refitTimes)
    pyplot.plot(np.arange(0.5, minInjectionRate-0.1, -0.1), scores)
    pyplot.title('Score evolution in relation to injection rate')
    pyplot.xlabel('injection rate')
    pyplot.ylabel('Best score')
    pyplot.show()


    pyplot.plot(np.arange(0.5, minInjectionRate-0.1, -0.1), refitTimes)
    pyplot.title('Refit times in relation to injection rate')
    pyplot.xlabel('injection rate')
    pyplot.ylabel('Refit time')
    pyplot.show()


def testVariability(maxVariability):
    variabilityType = metaParameters['sequenceInjection']['variabilityType']
    injectedSequence = metaParameters['sequenceInjection']['injectedSequence']
    randomSate = random.randint(1, 100)

    scores = []
    refitTimes = []
    start_time = time.time()
    for i in range(0, maxVariability + 1):
        configurations = []
        configurationScores = []
        configurationRefitTimes = []
        variability = i
        metaParameters['sequenceInjection']['variability'] = variability
        noOfConfigurations = 1
        if variabilityType == 'fixedPosition' and variability != 0:
            noOfConfigurations = int(sys.argv[3])
        print("\nINJECTED SEQUENCE VARIABILITY : " + str(variability) + ', type : ' + str(sys.argv[2]) + ', number of trials per variability: ' + str(noOfConfigurations) + "\n")
        for j in range(0, int(noOfConfigurations)):
            if variabilityType == 'fixedPosition':
                configuration = random.sample(range(0, len(injectedSequence)), variability)
                metaParameters['sequenceInjection']['fixedVariabilityIndexes'] = configuration
                configurations.append(configuration)
                print("Variable positions configuration (" + str(j + 1) + "/" + str(noOfConfigurations) + ") : " + str(configuration))
            else:
                print("Variability " + str(variability) + " test (" + str(j + 1) + "/" + str(noOfConfigurations) + ")")

            (score, refitTime, model) = injectAndEvaluateModel('MLP', "testVariabilityResults.txt", (i == 0 and j == 0))

            configurationScores.append(score)
            configurationRefitTimes.append(refitTime)
            print("\n-------------------------------\n")
        if noOfConfigurations > 1 and variabilityType == 'fixedPosition':
            print("Configuration scores : " + str(zip(configurations, configurationScores)))
            pyplot.boxplot(configurationScores)
            pyplot.title('Configuration score distribution, variability ' + str(variability))
            pyplot.show()
            print("Configuration refit times : " + str(configurationRefitTimes) + "\n\n==============================\n")
        scores.append(statistics.mean(configurationScores))
        refitTimes.append(statistics.mean(configurationRefitTimes))
    print("%s executed in --- %s seconds ---" % (sys.argv[1], time.time() - start_time))

    print("Scores : ")
    print(scores)
    print("RefitTimes : ")
    print(refitTimes)
    pyplot.plot(np.arange(0, maxVariability + 1), scores)
    pyplot.title('Score evolution in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Best score')
    pyplot.show()

    pyplot.plot(np.arange(0, maxVariability + 1), refitTimes)
    pyplot.title('Refit times in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Refit time')
    pyplot.show()


def testPhageHostRecognition(bacteriaName, preprocess=False, negativeDiscardRate=0.80):
    scores = []
    refitTimes = []
    start_time = time.time()

    trainingFolder = phageSeqFolder + "Training/"
    testingFolder = phageSeqFolder + "Testing/"

    if preprocess:
        buildPhageDictionary()
        phageSeqFile = open(phageSeqPath)
        labelsPath = phageSeqFolder + "labels.txt"
        labelFile = open(labelsPath, 'w')

        print("Writing labels")
        headerLine = phageSeqFile.readline()
        line = phageSeqFile.readline()
        lineCounter = 0
        positiveLabelCounter = 0
        while line:
            lineCounter += 1
            lineElements = line.split(',')
            designation = lineElements[2]
            if bacteriaName in designation.lower():
                label = 1
                positiveLabelCounter += 1
                labelFile.write(str(label) + '\n')
            else:
                label = 0
                if random.random() > 0.80:
                    labelFile.write(str(label) + '\n')
                else:
                    labelFile.write("discard\n")
                    lineCounter -= 1
            line = phageSeqFile.readline()
        labelFile.close()
        phageSeqFile.close()

        print("Proportion of positive samples :" + str(positiveLabelCounter / lineCounter))

        print("Splitting training and testing data")
        utils.splitTrainingAndTesting(phageSeqPath, labelsPath, trainingFolder, testingFolder)

        print("Making BOWs (training)")
        ReadSequenceFile(makeBOW, phageKmerDictionaryPath, trainingFolder + "seq.csv", trainingFolder + "Seq_BOW.json")
        print("Making BOWs (testing)")
        ReadSequenceFile(makeBOW, phageKmerDictionaryPath, testingFolder + "seq.csv", testingFolder + "Seq_BOW.json")

    (score, refitTime, model) = modelTrainingAndTesting.trainAndTestModel(
        "RandomForest",
        trainingFolder + "Seq_BOW.json",
        trainingFolder + "labels.txt",
        testingFolder + "Seq_BOW.json",
        testingFolder + "labels.txt",
        resultsFilePath='testPhageHostRecognition.txt',
        verbose=4
    )


if sys.argv[1] == "help":
    displayHelp()

elif sys.argv[1] == "dictionary":
    if sys.argv[2] == "new":
        ReadSequenceFile(NEW, sys.argv[3], sys.argv[4])
    if sys.argv[2] == "add":
        ReadSequenceFile(ADD, sys.argv[3], sys.argv[4])
    if sys.argv[2] == "filter":
        Filter(sys.argv[3])
    if sys.argv[2] == "stats":
        dictionaryStats()
    if len(sys.argv) == 1:
        buildPhageDictionary()

elif sys.argv[1] == "makeBOWs":
    ReadSequenceFile(makeBOW, sys.argv[2], sys.argv[3], sys.argv[4])

elif sys.argv[1] == "convertCSVtoFASTA":
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

elif sys.argv[1] == "plotMeanSequenceLengths":
    utils.plotAverageFoldLength(utils.getSequenceLengths(sys.argv[2]), sys.argv[3])

elif sys.argv[1] == "findAlignments":
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")

elif sys.argv[1] == "splitTrainingAndTesting":
    utils.splitTrainingAndTesting(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

elif sys.argv[1] == "injectSequence":
    injectSequence(sys.argv[2], sys.argv[3], sys.argv[4])

elif sys.argv[1] == "train":
    start_time = time.time()
    modelTrainingAndTesting.trainAndTestModel(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    print("--- %s seconds ---" % (time.time() - start_time))

elif sys.argv[1] == "makePlot":
    data = metaParameters['plotData']
    data = data[0]


elif sys.argv[1] == "testVariability":
    maxVariability = int(sys.argv[2])
    testVariability(maxVariability)

elif sys.argv[1] == 'testInjectionRate':
    testInjectionRate()


elif sys.argv[1] == 'testPhageHostRecognition':
    preprocess = len(sys.argv) > 3 and sys.argv[3] == 'preprocess'
    testPhageHostRecognition(sys.argv[2], preprocess)

else:
    displayHelp()
