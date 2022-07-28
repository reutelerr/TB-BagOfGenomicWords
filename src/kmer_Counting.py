import json
import os
import threading
import time

import numpy as np
from tqdm import tqdm

from src.constants import CSV, FASTA, metaParameters
from src.utils import blocks, timerWrapper

seqLimiter = 100000

metaParams = metaParameters.get('vectorization')
maxKmerLength = metaParams.get('maxKmerLength')
minimumFrequencyThreshold = metaParams.get('minimumFrequencyThreshold')
filterThresholdLengthModifier = metaParams.get('filterThresholdLengthModifier')

#Gets the maximum frequency of this sequence as a subsequence of longer ones
def occurenceFreqAsSubSeq(seq, frequent_seq_indexes):
    maxOccurenceFreq = 0
    for seqMap in frequent_seq_indexes:
        for frequentSeq in seqMap:

            # Is seq a subsequence of frequentSeq, and does it occur more than the current max ?
            if seq in frequentSeq and seqMap[frequentSeq] > maxOccurenceFreq:
                maxOccurenceFreq = seqMap[frequentSeq]
    return maxOccurenceFreq


def countKmers(sequence, kmerCounts, nucleotideCounter, kmerLength, newCounts) :
    for i in range(len(sequence) - kmerLength + 1):
        nucleotideCounter += 1
        kmer = sequence[i:i + kmerLength]
        if kmerCounts.__contains__(kmer):
            kmerCounts[kmer] += 1
        else:
            if not newCounts:
                kmerCounts[kmer] = 1
    return nucleotideCounter + kmerLength


def writeDictionaryPage(sequences, length, dictionaryPath):
    sequenceCount = 0
    sequence = ""
    kmerCounts = {}
    nucleotideCounter = 0

    dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"

    with tqdm(total=len(sequences), position=0, leave=True) as pbar:
        for sequence in sequences:
            nucleotideCounter = countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
            if length == maxKmerLength:
                pbar.update(1)
    os.makedirs(os.path.dirname(dictionaryFileName), exist_ok=True)
    dictionaryFile = open(dictionaryFileName, "w")
    dictionaryFile.write(json.dumps((nucleotideCounter, kmerCounts), indent=4))
    dictionaryFile.close()


@timerWrapper
def vectorize(dictionaryPath, sourcePath, outputPath, sourceType):
    seqFile = open(sourcePath)

    noLines = sum(bl.count("\n") for bl in blocks(seqFile))
    print("Total number of lines in file : "+str(noLines))
    os.makedirs(os.path.dirname(outputPath), exist_ok=True)
    output = open(outputPath, "w")

    seqFile.seek(0)
    time.sleep(0.01)

    sequenceCount = 0
    sequence = ""
    validKmers = {}
    kmerCounts = {}
    nucleotideCounter = 0

    BOWs = []

    with tqdm(total=noLines, position=0, leave=True) as pbar:

        if sourceType == CSV:
            headerLine = seqFile.readline()
            line = seqFile.readline()
            while line:
                pbar.update(1)
                lineElements = line.split(',')
                id = lineElements[0]
                sequence = lineElements[-1].lower()

                nucleotideCounter = 0
                for length in range(maxKmerLength, 2, -1):

                    dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"
                    dictionaryFile = open(dictionaryFileName, "rb")
                    (nucleotideCount, size, validKmers) = json.load(dictionaryFile)
                    for entry in validKmers:
                        validKmers[entry] = 0
                    nucleotideCounter = countKmers(sequence, validKmers, nucleotideCounter, length, True)
                    kmerCounts.update(validKmers)
                BOWs = ((id, nucleotideCounter, len(kmerCounts), kmerCounts))
                output.write(json.dumps(BOWs, indent=4))
                output.write("\n===\n")
                line = seqFile.readline()

            dictionaryFile.close()
            seqFile.close()
            output.close()

        if sourceType == FASTA:
            for line in seqFile:
                pbar.update(1)
                if line[0] == '>':
                    if sequence != "":
                        nucleotideCounter = countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
                    sequence = ""
                    sequenceCount += 1
                else:
                    rawLine = line[:len(line) - 1]
                    sequence = sequence + rawLine  # Stitching the sequence back together line by line>
            if sequence != "":
                countKmers(sequence, kmerCounts, nucleotideCounter, length, False)



@timerWrapper
def buildDictionary(mode, dictionaryPath, sourcePath, sourceType=CSV):
    seqFile = open(sourcePath)

    noLines = sum(bl.count("\n") for bl in blocks(seqFile))
    print("Total number of lines in file : "+str(noLines))

    for length in range(maxKmerLength, 2, -1):
        print("Counting K-mers of length : "+str(length))
        seqFile.seek(0)
        time.sleep(0.01)

        sequenceCount = 0
        sequence = ""
        kmerCounts = {}
        nucleotideCounter = 0

        dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"

        with tqdm(total=noLines, position=0, leave=True) as pbar:

            #Reading the file
            if sourceType == FASTA:
                for line in seqFile:
                    pbar.update(1)
                    if sequenceCount < seqLimiter:
                        if line[0] == '>':
                            if sequence != "":
                                nucleotideCounter = countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
                            sequence = ""
                            sequenceCount+=1
                        else:
                            rawLine = line[:len(line)-1]
                            sequence = sequence + rawLine#Stitching the sequence back together line by line
                if sequence != "":
                    countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
            if sourceType == CSV:
                headers = seqFile.readline().split(',')
                line = seqFile.readline()
                while line:
                    pbar.update(1)
                    lineElements = line.split(',')
                    id = lineElements[0]
                    sequence = lineElements[-1].lower()[:len(lineElements[-1])-1]
                    nucleotideCounter = countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
                    line = seqFile.readline()

            os.makedirs(os.path.dirname(dictionaryFileName), exist_ok=True)
            dictionaryFile = open(dictionaryFileName, "w")
            dictionaryFile.write(json.dumps((nucleotideCounter, kmerCounts), indent=4))
        dictionaryFile.close()
    seqFile.close()



@timerWrapper
def filterByFrequency(dictionaryFilePath):
    frequent_kmers = [None] * maxKmerLength
    for length in range(maxKmerLength, 2, -1):
        kmerCountFile = open(dictionaryFilePath + "/kmerCounts" + str(length) + ".json", "rb")
        (nucleotideCount, kmerCounts) = json.load(kmerCountFile)
        frequent_kmers[length - 1] = {}

        for kmer in kmerCounts:
            noOfOccurences = kmerCounts[kmer]
            if noOfOccurences > nucleotideCount/pow(4, len(kmer)/filterThresholdLengthModifier) * minimumFrequencyThreshold:
                dic = frequent_kmers[length - 1]
                dic[kmer] = kmerCounts[kmer]

        kmerCountFile.close()
        kmerCountFile = open(dictionaryFilePath + "/kmerCounts" + str(length) + ".json", "w")
        kmerCountFile.write(json.dumps((nucleotideCount, len(frequent_kmers[length-1]), frequent_kmers[length-1]), indent=4))
        kmerCountFile.close()
















