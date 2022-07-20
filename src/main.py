import time
import sys

from matplotlib import pyplot

import kmer_Counting
import utils
import modelTrainingAndTesting
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

if (sys.argv[1] == "help"):
    displayHelp()

if (sys.argv[1] == "dictionary"):
    if (sys.argv[2] == "new"):
        ReadSequenceFile(NEW, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "add"):
        ReadSequenceFile(ADD, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "filter"):
        Filter(sys.argv[3])

if (sys.argv[1] == "makeBOWs"):
    ReadSequenceFile(makeBOW, sys.argv[2], sys.argv[3], sys.argv[4])

if (sys.argv[1] == "convertCSVtoFASTA"):
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

if (sys.argv[1] == "plotMeanSequenceLengths"):
    utils.plotAverageFoldLength(utils.getSequenceLengths(sys.argv[2]), sys.argv[3])

if (sys.argv[1] == "findAlignments"):
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")

if (sys.argv[1] == "splitTrainingAndTesting"):
    utils.splitTrainingAndTesting(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

if (sys.argv[1] == "injectSequence"):
    injectSequence(sys.argv[2], sys.argv[3], sys.argv[4])

if (sys.argv[1] == "train"):
    start_time = time.time()
    modelTrainingAndTesting.trainAndTestModel(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    print("--- %s seconds ---" % (time.time() - start_time))

if (sys.argv[1] == "testVariability"):
    sourcePath = "sequences/PhageWholeDNASeq/PhageWholeDNASeq.csv"
    injectedSeqPath = "sequences/InjectedPhageWholeDNASeq/InjectedPhageWholeDNASeq.csv"
    labelsPath = "sequences/InjectedPhageWholeDNASeq/labels.txt"
    trainingSeqPath = "sequences/InjectedPhageWholeDNASeq/Training/seq.csv"
    trainingLabelsPath = "sequences/InjectedPhageWholeDNASeq/Training/labels.txt"
    testingSeqPath = "sequences/InjectedPhageWholeDNASeq/Testing/seq.csv"
    testingLabelsPath = "sequences/InjectedPhageWholeDNASeq/Testing/labels.txt"
    trainingBOWsPath = "sequences/InjectedPhageWholeDNASeq/Training/Seq_BOW.json"
    testingBOWsPath = "sequences/InjectedPhageWholeDNASeq/Testing/Seq_BOW.json"
    dictionaryPath = "kmerDictionaries/dic10"

    scores = []
    refitTimes = []
    start_time = time.time()
    for i in range(0, 7):
        variability = i*0.05
        print("Injecting Sequence, variability : "+str(variability))
        metaParameters['sequenceInjection']['variability'] = variability
        injectSequence(sourcePath, injectedSeqPath, labelsPath)
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
        (score, refitTime) = modelTrainingAndTesting.trainAndTestModel("Perceptron", trainingBOWsPath, trainingLabelsPath, testingBOWsPath, testingLabelsPath, verbose=4)
        scores.append(score)
        refitTimes.append(refitTime)
        print("\n-------------------------------\n")
    print("%s executed in --- %s seconds ---" % (sys.argv[1], time.time() - start_time))

    pyplot.bar(0.05*np.arange(0, 7), scores)
    pyplot.xlabel('Variability')
    pyplot.ylabel('Best score')
    pyplot.show()

    pyplot.bar(0.05 * np.arange(0, 7), refitTimes)
    pyplot.xlabel('Variability')
    pyplot.ylabel('Refit time')
    pyplot.show()
