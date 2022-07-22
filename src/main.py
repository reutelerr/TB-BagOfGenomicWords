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


def ReadSequenceFile(mode, dictionaryPath, sourcePath, outputPath =''):
    fileExtension = sourcePath.split('.')[-1]
    if fileExtension == "csv" : sourceType = CSV
    if fileExtension == "fasta" : sourceType = FASTA
    if fileExtension == "txt" : sourceType = FASTA

    if mode == makeBOW:
        kmer_Counting.vectorize(dictionaryPath, sourcePath, outputPath, sourceType)
    else:
        kmer_Counting.buildDictionary(mode, dictionaryPath, sourcePath, sourceType)


def Filter(dicPath):
    kmer_Counting.filterByFrequency(dicPath)

def PlotDictionaryPageSizes():
    dictionaryPath = "kmerDictionaries/dic10"
    maxKmerLength = metaParameters['vectorization']['maxKmerLength']
    for i in range(10, maxKmerLength):
        print("Building dictionary")
        ReadSequenceFile(NEW, dictionaryPath, injectedSeqPath)
        print("Filtering dictionary")
        Filter('kmerDictionaries/dic10')
        #TODO

def dictionaryStats():
    print("Building dictionary")
    ReadSequenceFile(NEW, dictionaryPath, injectedSeqPath)
    maxKmerLength = metaParameters['vectorization']['maxKmerLength']
    for length in range(3, maxKmerLength):
        dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"
        dictionaryFile = open(dictionaryFileName, "rb")
        (nucleotideCount, validKmers) = json.load(dictionaryFile)
        noOfPossibleKmers = 4 ** maxKmerLength
        noOfActualKmers = validKmers
        noOf0OccurenceKmers = noOfPossibleKmers-noOfActualKmers
        print("No of kmers of length "+str(length)+"not in dataset : "+noOf0OccurenceKmers)
        pyplot.boxplot(list(validKmers.values))
        pyplot.title('Frequency distribution of '+str(length)+'-mers')
        pyplot.show()

    print("Filtering dictionary")
    Filter('kmerDictionaries/dic10')
    dicSizes = []
    for length in range(3, maxKmerLength+1):
        dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"
        dicSize = os.path.getsize(dictionaryFileName)
        dicSizes.append(dicSize)
        print("Dictionary size for length"+str(length)+" (in bytes): "+str(dicSize))
    pyplot.bar(np.arange(3, maxKmerLength+1), dicSizes)
    frequencyThreshold = metaParameters['vectorization']['minimumFrequencyThreshold']
    frequencyThresholdModifier = metaParameters['vectorization']['filterThresholdLengthModifier']
    pyplot.title('Filtered Dictionary Sizes, with threshold = '+str(frequencyThreshold)+'; modifier = '+str(frequencyThresholdModifier))
    pyplot.xlabel('kmerLength')
    pyplot.ylabel('Size')
    pyplot.show()





if (sys.argv[1] == "help"):
    displayHelp()

elif (sys.argv[1] == "dictionary"):
    if (sys.argv[2] == "new"):
        ReadSequenceFile(NEW, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "add"):
        ReadSequenceFile(ADD, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "filter"):
        Filter(sys.argv[3])

elif (sys.argv[1] == "makeBOWs"):
    ReadSequenceFile(makeBOW, sys.argv[2], sys.argv[3], sys.argv[4])

elif (sys.argv[1] == "convertCSVtoFASTA"):
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

elif (sys.argv[1] == "plotMeanSequenceLengths"):
    utils.plotAverageFoldLength(utils.getSequenceLengths(sys.argv[2]), sys.argv[3])

elif (sys.argv[1] == "findAlignments"):
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")

elif (sys.argv[1] == "splitTrainingAndTesting"):
    utils.splitTrainingAndTesting(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

elif (sys.argv[1] == "injectSequence"):
    injectSequence(sys.argv[2], sys.argv[3], sys.argv[4])

elif (sys.argv[1] == "train"):
    start_time = time.time()
    modelTrainingAndTesting.trainAndTestModel(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    print("--- %s seconds ---" % (time.time() - start_time))

elif (sys.argv[1] == "testVariability"):
    maxVariability = 10
    variabilityType = sys.argv[2]
    injectedSequence = metaParameters['sequenceInjection']['injectedSequence']
    randomSate = random.randint(1, 100)

    scores = []
    refitTimes = []
    start_time = time.time()
    for i in range(0, maxVariability+1):
        configurationScores = []
        variability = i
        configurationNo = 1
        if variabilityType == 'fixedPosition' and variability != 0:
            configurationNo = sys.argv[3]

        print("Injecting Sequence, variability : " + str(variability) + ', type : ' + str(sys.argv[2]) + ', number of trials per variability: ' + str(configurationNo))
        for i in range(0, int(configurationNo)):
            if variabilityType == 'fixedPosition':
                variabilityIndexes = random.sample(range(0, len(injectedSequence)), variability)
                metaParameters['sequenceInjection']['fixedVariabilityIndexes'] = variabilityIndexes
                print("Variable positions configuration (" + str(i+1) + "/" + str(6) + ") : " + str(variabilityIndexes))
            metaParameters['sequenceInjection']['variability'] = variability

            injectSequence(sourcePath, injectedSeqPath, labelsPath, variabilityType=variabilityType)
            if i==0:
                print("Building dictionary")
                ReadSequenceFile(NEW, dictionaryPath, injectedSeqPath)
                print("Filtering dictionary")
                Filter('kmerDictionaries/dic10')
            print("Splitting training and testing data")
            utils.splitTrainingAndTesting(injectedSeqPath, labelsPath, trainingSeqPath, trainingLabelsPath, testingSeqPath, testingLabelsPath)
            print("Making BOWs (training)")
            ReadSequenceFile(makeBOW, dictionaryPath, trainingSeqPath, trainingBOWsPath)
            print("Making BOWs (testing)")
            ReadSequenceFile(makeBOW, dictionaryPath, testingSeqPath, testingBOWsPath)
            print("Building model")
            (score, refitTime) = modelTrainingAndTesting.trainAndTestModel(
                "RandomForest",
                trainingBOWsPath,
                trainingLabelsPath,
                testingBOWsPath,
                testingLabelsPath,
                randomState=randomSate,
                resultsFilePath='testVariabilityResults.txt',
                verbose=4
            )
            configurationScores.append(score)
            refitTimes.append(refitTime)
            print("\n-------------------------------\n")
        scores.append(statistics.mean(configurationScores))
    print("%s executed in --- %s seconds ---" % (sys.argv[1], time.time() - start_time))

    print("Scores : ")
    print(scores)
    print("RefitTimes : ")
    print(refitTimes)
    pyplot.plot(np.arange(0, maxVariability+1), scores)
    pyplot.title('Score evolution in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Best score')
    pyplot.show()

    pyplot.plot(np.arange(0, maxVariability+1), refitTimes)
    pyplot.title('Refit times in relation to variability')
    pyplot.xlabel('Variability')
    pyplot.ylabel('Refit time')
    pyplot.show()

else:
    displayHelp()
