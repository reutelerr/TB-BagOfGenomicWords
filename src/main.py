import json
import os.path
import time
import sys

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
    if fileExtension == "csv": sourceType = CSV
    if fileExtension == "fasta": sourceType = FASTA
    if fileExtension == "txt": sourceType = FASTA

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
    for length in range(3, maxKmerLength):
        dictionaryFileName = phageKmerDictionaryPath + "/kmerCounts" + str(length) + ".json"
        dictionaryFile = open(dictionaryFileName, "rb")
        (nucleotideCount, validKmers) = json.load(dictionaryFile)
        noOfPossibleKmers = 4 ** maxKmerLength
        noOfActualKmers = validKmers
        noOf0OccurenceKmers = noOfPossibleKmers - noOfActualKmers
        print("No of kmers of length " + str(length) + "not in dataset : " + noOf0OccurenceKmers)
        pyplot.boxplot(list(validKmers.values))
        pyplot.title('Frequency distribution of ' + str(length) + '-mers')
        pyplot.show()

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


def injectAndEvaluateModel(modelType):
    variabilityType = metaParameters['sequenceInjection']['variabilityType']
    labelsPath = injectedPhageSeqFolder + "labels.txt"
    trainingPath = injectedPhageSeqFolder + "Training/"
    trainingSeqPath = trainingPath + "seq.csv"
    trainingLabelsPath = trainingPath + "labels.txt"
    testingPath = injectedPhageSeqFolder + "Testing/"
    testingSeqPath = testingPath + "seq.csv"
    testingLabelsPath = testingPath + "labels.txt"
    trainingBOWsPath = trainingPath + "/Seq_BOW.json"
    testingBOWsPath = testingPath + "/Seq_BOW.json"

    injectSequence(phageSeqPath, injectedPhageSeqPath, variabilityType=variabilityType)
    print("Splitting training and testing data")
    utils.splitTrainingAndTesting(injectedPhageSeqPath, labelsPath, trainingPath, testingPath)
    print("Making BOWs (training)")
    ReadSequenceFile(makeBOW, phageKmerDictionaryPath, trainingSeqPath, trainingBOWsPath)
    print("Making BOWs (testing)")
    ReadSequenceFile(makeBOW, phageKmerDictionaryPath, testingSeqPath, testingBOWsPath)
    print("Building model")
    (score, refitTime) = modelTrainingAndTesting.trainAndTestModel(
        modelType,
        trainingBOWsPath,
        trainingLabelsPath,
        testingBOWsPath,
        testingLabelsPath,
        randomState=randomSate,
        resultsFilePath='testVariabilityResults.txt',
        verbose=4
    )
    return score, refitTime


def testPhageHostRecognition(bacteriaName, preprocess=False):
    scores = []
    refitTimes = []
    start_time = time.time()

    buildPhageDictionary()

    trainingFolder = phageSeqFolder + "Training/"
    testingFolder = phageSeqFolder + "Testing/"

    if preprocess:
        phageSeqFile = open(phageSeqPath)
        labelsPath = phageSeqFolder + "labels.txt"
        labelFile = open(labelsPath, 'w')

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
            else:
                label = 0
            labelFile.write(str(label) + '\n')
            line = phageSeqFile.readline()

        print("Splitting training and testing data")
        utils.splitTrainingAndTesting(phageSeqPath, labelsPath, trainingFolder, testingFolder)

        print("Proportion of positive samples :" + str(positiveLabelCounter / lineCounter))
        print("Making BOWs (training)")
        ReadSequenceFile(makeBOW, phageKmerDictionaryPath, trainingFolder + "/seq.csv", trainingFolder + "/Seq_BOW.json")
        print("Making BOWs (testing)")
        ReadSequenceFile(makeBOW, phageKmerDictionaryPath, testingFolder + "/seq.csv", testingFolder + "/Seq_BOW.json")

    (score, refitTime) = modelTrainingAndTesting.trainAndTestModel(
        "MLP",
        trainingFolder + "/Seq_BOW.json",
        trainingFolder + "/labels.txt",
        testingFolder + "/Seq_BOW.json",
        testingFolder + "/labels.txt",
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
    if len(sys.argv) == 1:
        buildDictionary()

elif sys.argv[1] == "makeBOWs":
    ReadSequenceFile(makeBOW, sys.argv[2], sys.argv[3], sys.argv[4])

elif sys.argv[1] == "convertCSVtoFASTA":
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

elif sys.argv[1] == "plotMeanSequenceLengths":
    utils.plotAverageFoldLength(utils.getSequenceLengths(sys.argv[2]), sys.argv[3])

elif sys.argv[1] == "findAlignments":
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")

elif sys.argv[1] == "splitTrainingAndTesting":
    utils.splitTrainingAndTesting(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

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
        if i == 0:
            print("Building dictionary")
            ReadSequenceFile(NEW, phageKmerDictionaryPath, injectedPhageSeqPath)
            print("Filtering dictionary")
            Filter(phageKmerDictionaryPath)
        print("Injected sequence variability : " + str(variability) + ', type : ' + str(
            sys.argv[2]) + ', number of trials per variability: ' + str(noOfConfigurations))
        for j in range(0, int(noOfConfigurations)):
            if variabilityType == 'fixedPosition':
                configuration = random.sample(range(0, len(injectedSequence)), variability)
                metaParameters['sequenceInjection']['fixedVariabilityIndexes'] = configuration
                configurations.append(configuration)
                print("Variable positions configuration (" + str(j + 1) + "/" + str(noOfConfigurations) + ") : " + str(
                    configuration))
            else:
                print("Variability " + str(variability) + " test (" + str(j + 1) + "/" + str(noOfConfigurations) + ")")

            (score, refitTime) = injectAndEvaluateModel('MLP')

            configurationScores.append(score)
            configurationRefitTimes.append(refitTime)
            print("\n-------------------------------\n")
        if noOfConfigurations > 1 and variabilityType == 'fixedPosition':
            print("Configuration scores : " + str(zip(configurations, configurationScores)))
            pyplot.boxplot(configurationScores)
            pyplot.title('Configuration score distribution, variability ' + str(variability))
            pyplot.show()
            print(
                "Configuration refit times : " + str(configurationRefitTimes) + "\n\n==============================\n")
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

elif sys.argv[1] == 'testInjectionRate':
    injectedSequence = metaParameters['sequenceInjection']['injectedSequence']
    minInjectionRate = int(sys.argv[2])

    randomSate = random.randint(1, 100)

    scores = []
    refitTimes = []
    start_time = time.time()
    injectionRates = np.arange(0.5, minInjectionRate, -0.1)
    for i in injectionRates:
        metaParameters['sequenceInjection']['injectionRate'] = i
        (score, refitTime) = injectAndEvaluateModel('KNN')
        scores.append(score)
        refitTimes.append(refitTimes)

    print("Scores : ")
    print(scores)
    print("RefitTimes : ")
    print(refitTimes)
    pyplot.plot(np.arange(0.5, injectionRates), scores)
    pyplot.title('Score evolution in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Best score')
    pyplot.show()

    pyplot.plot(np.arange(0, injectionRates), refitTimes)
    pyplot.title('Refit times in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Refit time')
    pyplot.show()


elif sys.argv[1] == 'testPhageHostRecognition':
    preprocess = len(sys.argv) > 2 and sys.argv[2] == 'preprocess'
    testPhageHostRecognition(sys.argv[2], preprocess)

else:
    displayHelp()
