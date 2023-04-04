import os
import random
import time

from tqdm import tqdm
from statistics import mean
from matplotlib import pyplot

from src.constants import metaParameters


def blocks(file, size=65536):
    while True:
        b = file.read(size)
        if not b: break
        yield b


def CSVtoFASTA(csvfilePath, fastaFilePath):
    csvFile = open(csvfilePath)
    fastaFile = open(fastaFilePath, 'w')

    noLines = sum(bl.count("\n") for bl in blocks(csvFile))
    csvFile.seek(0)
    with tqdm(total=noLines, position=0, leave=True) as pbar:
        headers = csvFile.readline().split(',')
        line = csvFile.readline()
        while line:
            pbar.update(1)
            lineElements = line.split(',')
            head = lineElements[3]
            startIndex = head.index('>')
            head = head[startIndex:] + ' ' + lineElements[4]
            head = head.replace('\"', '')
            sequence = lineElements[-1].lower()

            fastaFile.write(head + '\n')
            fastaFile.write(sequence)
            fastaFile.write('\n')
            line = csvFile.readline()

    fastaFile.close()
    csvFile.close()


def FASTAtoCSV(fastaFilePath, csvFilePath):
    csvFile = open(csvFilePath, 'w')
    fastaFile = open(fastaFilePath)
    noLines = sum(bl.count("\n") for bl in blocks(fastaFile))
    fastaFile.seek(0)
    with tqdm(total=noLines, position=0, leave=True) as pbar:
        headers = "fasta_head, sequence"
        line = csvFile.readline()
        while line:
            pbar.update(1)
            lineElements = line.split(',')
            head = lineElements[3]
            startIndex = head.index('>')
            head = head[startIndex:] + ' ' + lineElements[4]
            head = head.replace('\"', '')
            sequence = lineElements[-1].lower()

            fastaFile.write(head + '\n')
            fastaFile.write(sequence)
            fastaFile.write('\n')
            line = csvFile.readline()

    fastaFile.close()
    csvFile.close()


# Randomly separates training and testing data
def splitTrainingAndTesting(seqFilePath, labelFilePath, trainingFolder, testingFolder, split=0.9, discard=0):
    seqFile = open(seqFilePath)
    labelFile = open(labelFilePath)

    os.makedirs(os.path.dirname(trainingFolder), exist_ok=True)
    trainingSeqFile = open(trainingFolder + "seq.csv", "w")
    trainingLabelFile = open(trainingFolder + "labels.txt", "w")
    os.makedirs(os.path.dirname(testingFolder), exist_ok=True)
    testingSeqFile = open(testingFolder + "seq.csv", "w")
    testingLabelFile = open(testingFolder + "labels.txt", "w")

    random.seed()

    headerLine = seqFile.readline()
    headers = headerLine.split(',')
    testingSeqFile.write(headerLine)
    trainingSeqFile.write(headerLine)
    line = seqFile.readline()
    label = labelFile.readline()

    while line:
        if random.random() > discard:
            if label != "discard\n":
                if random.random() < split:
                    trainingSeqFile.write(line)
                    trainingLabelFile.write(label)
                else:
                    testingSeqFile.write(line)
                    testingLabelFile.write(label)
        line = seqFile.readline()
        label = labelFile.readline()

    seqFile.close()
    labelFile.close()
    trainingSeqFile.close()
    trainingLabelFile.close()
    testingSeqFile.close()
    testingLabelFile.close()


def getSequenceLengths(seqFilePath):
    seqLengths = []
    seqFile = open(seqFilePath)

    headerLine = seqFile.readline()
    headers = headerLine.split(',')
    line = seqFile.readline()
    while line:
        lineElements = line.split(',')
        id = lineElements[0]
        sequence = lineElements[-1].lower()
        seqLengths.append(len(sequence))
        line = seqFile.readline()
    print("Number of sequences :" + str(len(seqLengths)))
    return seqLengths


def saveToFile(data, resultsFilePath='results.txt'):
    resultsFile = open(resultsFilePath, "a")
    resultsFile.write(data)
    resultsFile.close()


def plotAverageFoldLength(sequenceLengths, cvFolds=metaParameters['modelTraining']['cvFolds']):
    meanSeqLengths = []
    cvFolds = int(cvFolds)
    for i in range(0, cvFolds):
        foldStart = int(i * len(sequenceLengths) / cvFolds)
        foldEnd = int((i + 1) * len(sequenceLengths) / cvFolds)  # Not included
        meanSeqLength = mean(sequenceLengths[foldStart:foldEnd])
        meanSeqLengths.append(meanSeqLength)
    pyplot.bar(range(1, cvFolds + 1), meanSeqLengths)
    pyplot.xlabel('CV-fold')
    pyplot.ylabel('Average Sequence Length')
    pyplot.show()


def timerWrapper(func):
    def wrap(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print("%s executed in --- %s seconds ---" % (func.__name__, time.time() - start))
        return result

    return wrap
