import time
import sys

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

    start_time = time.time()
    if mode == makeBOW:
        kmer_Counting.vectorize(dictionaryPath, sourcePath, outputPath, sourceType)
    else:
        kmer_Counting.buildDictionary(mode, dictionaryPath, sourcePath, sourceType)
    print("--- %s seconds ---" % (time.time() - start_time))

def Filter(dicPath):
    start_time = time.time()
    kmer_Counting.filterByFrequency(dicPath)
    print("--- %s seconds ---" % (time.time() - start_time))

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
    sourcePath = "sequences/PhageContigsSeq/PhageContigsSeq.csv"
    injectedSeqPath = "sequences/InjectedPhageContigsSeq/InjectedPhageContigsSeq.csv"
    labelsPath = "sequences/InjectedPhageContigsSeq/labels.txt"
    trainingSeqPath = "sequences/InjectedPhageContigsSeq/Training/seq.csv"
    trainingLabelsPath = "sequences/InjectedPhageContigsSeq/Training/labels.txt"
    testingSeqPath = "sequences/InjectedPhageContigsSeq/Testing/seq.csv"
    testingLabelsPath = "sequences/InjectedPhageContigsSeq/Testing/labels.txt"
    trainingBOWsPath = "sequences/InjectedPhageContigsSeq/Training/Seq_BOW.json"
    testingBOWsPath = "sequences/InjectedPhageContigsSeq/Testing/Seq_BOW.json"
    dictionaryPath = "kmerDictionaries/dic10"

    start_time = time.time()
    for i in range(0, 10):
        variability = i*0.05
        print("Injecting Sequence, variability : "+str(variability))
        injectSequence(sourcePath, injectedSeqPath, labelsPath, variability=variability)
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
        modelTrainingAndTesting.trainAndTestModel("Perceptron", trainingBOWsPath, trainingLabelsPath, testingBOWsPath, testingLabelsPath, verbose=3)
        print("\n-------------------------------\n")
