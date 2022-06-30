import json
import os
import sys
import time

import numpy as np
from tqdm import tqdm
import pickle

from src.constants import CSV, FASTA, makeBOW
from src.utils import blocks

minimumFrequencyThreshold = 2.0
seqLimiter = 100000
maxLength = 10



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
            if not newCounts: #onlyExisting determines wether we are counting kmers only already belonging to the dictionary or adding new ones
                kmerCounts[kmer] = 1
    return nucleotideCounter + kmerLength

def vectorize(dictionaryPath, sourcePath, outputPath, sourceType):
    sequence_file = open(sourcePath)

    noLines = sum(bl.count("\n") for bl in blocks(sequence_file))
    print("Total number of lines in file : "+str(noLines))
    os.makedirs(os.path.dirname(outputPath), exist_ok=True)
    output = open(outputPath, "w")

    sequence_file.seek(0)
    time.sleep(0.01)

    sequenceCount = 0
    sequence = ""
    validKmers = {}
    kmerCounts = {}
    nucleotideCounter = 0

    with tqdm(total=noLines, position=0, leave=True) as pbar:

        if sourceType == CSV:
            headers = sequence_file.readline().split(',')
            line = sequence_file.readline()
            while line:
                pbar.update(1)
                lineElements = line.split(',')
                id = lineElements[0]
                sequence = lineElements[-1].lower()

                nucleotideCounter = 0
                for length in range(maxLength, 2, -1):

                    dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"
                    dictionaryFile = open(dictionaryFileName, "rb")
                    (nucleotideCount, validKmers) = json.load(dictionaryFile)

                    for entry in validKmers:
                        validKmers[entry] = 0
                    nucleotideCounter = countKmers(sequence, validKmers, nucleotideCounter, length, True)
                    kmerCounts.update(validKmers)
                output.write(json.dumps((id, nucleotideCounter, kmerCounts), indent=4))
                line = sequence_file.readline()
            dictionaryFile.close()
            sequence_file.close()
            output.close()




def buildDictionary(mode, dictionaryPath, sourcePath, sourceType = CSV):
    sequence_file = open(sourcePath)

    noLines = sum(bl.count("\n") for bl in blocks(sequence_file))
    print("Total number of lines in file : "+str(noLines))

    for length in range(maxLength, 2, -1):
        print("Counting K-mers of length : "+str(length))
        sequence_file.seek(0)
        time.sleep(0.01)

        sequenceCount = 0
        sequence = ""
        kmerCounts = {}
        nucleotideCounter = 0

        dictionaryFileName = dictionaryPath + "/kmerCounts" + str(length) + ".json"

        with tqdm(total=noLines, position=0, leave=True) as pbar:

            #Reading the file
            if sourceType == FASTA:
                for line in sequence_file:
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
                            # TODO find a way to deal with very long sequences
                if sequence != "":
                    countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
            if sourceType == CSV:
                headers = sequence_file.readline().split(',')
                line = sequence_file.readline()
                while line:
                    pbar.update(1)
                    lineElements = line.split(',')
                    id = lineElements[0]
                    sequence = lineElements[-1].lower()[:len(lineElements[-1])-1]
                    nucleotideCounter = countKmers(sequence, kmerCounts, nucleotideCounter, length, False)
                    line = sequence_file.readline()

            os.makedirs(os.path.dirname(dictionaryFileName), exist_ok=True)
            dictionaryFile = open(dictionaryFileName, "w")
            dictionaryFile.write(json.dumps((nucleotideCounter, kmerCounts), indent=4))
        dictionaryFile.close()
    sequence_file.close()


def filterByFrequency(dictionaryFilePath):
    frequent_kmers = [None]*maxLength
    for length in range(maxLength, 2, -1):
        kmerCountFile = open(dictionaryFilePath + "/kmerCounts" + str(length) + ".json", "rb")
        (nucleotideCount, kmerCounts) = json.load(kmerCountFile)
        frequent_kmers[length - 1] = {}


        for kmer in tqdm(kmerCounts):
            noOfOccurences = kmerCounts[kmer]
            if noOfOccurences > nucleotideCount/pow(4, len(kmer)) * minimumFrequencyThreshold:
                # occurences independent of longer frequent seq
                #noOfEffectiveOccurences = kmerCounts[kmer] - occurenceFreqAsSubSeq(kmer, frequent_kmers[length:maxLength - 1])
                #if noOfEffectiveOccurences > nucleotideCount/pow(4, len(kmer)) * minimumFrequencyThreshold :  # Removing any non-frequent kner
                    dic = frequent_kmers[length - 1]
                    dic[kmer] = kmerCounts[kmer]

        kmerCountFile.close()
        kmerCountFile = open(dictionaryFilePath + "/kmerCounts" + str(length) + ".json", "w")
        kmerCountFile.write(json.dumps((nucleotideCount, frequent_kmers[length-1]), indent=4))
        kmerCountFile.close()
















