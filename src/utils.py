import os
import random
import time

from tqdm import tqdm

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
            head = head[startIndex:]+' '+lineElements[4]
            head = head.replace('\"', '')
            sequence = lineElements[-1].lower()

            fastaFile.write(head+'\n')
            fastaFile.write(sequence)
            fastaFile.write('\n')
            line = csvFile.readline()

    fastaFile.close()
    csvFile.close()

#Randomly separates training and testing data
def splitTrainingAndTesting(seqFilePath, labelFilePath, trainingSeqFilePath, trainingLabelFilePath, testingSeqFilePath, testingLabelFilePath, split = 0.9):
    seqFile = open(seqFilePath)
    labelFile = open(labelFilePath)

    os.makedirs(os.path.dirname(trainingSeqFilePath), exist_ok=True)
    trainingSeqFile = open(trainingSeqFilePath, "w")
    trainingLabelFile = open(trainingLabelFilePath, "w")
    os.makedirs(os.path.dirname(testingSeqFilePath), exist_ok=True)
    testingSeqFile = open(testingSeqFilePath, "w")
    testingLabelFile = open(testingLabelFilePath, "w")

    random.seed()

    line = seqFile.readline()
    label = labelFile.readline()

    while line:
        if(random.random()<split):
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

def timerWrapper(func):
    def wrap(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print("--- %s seconds ---" % (time.time() - start))
        return result

    return wrap