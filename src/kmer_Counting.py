import json
import os
import sys
import time

import numpy as np
from tqdm import tqdm
import pickle

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

def countKmers(sequence, sequenceCounts, nucleotideCounter, kmerLength) :
    for i in range(len(sequence) - kmerLength + 1):
        nucleotideCounter += 1
        kmer = sequence[i:i + kmerLength]
        if sequenceCounts.__contains__(kmer):
            sequenceCounts[kmer] += 1
        else:
            sequenceCounts[kmer] = 1
    return nucleotideCounter + kmerLength


def readSequenceFile(mode, indexFilePath, sourcePath, sourceType):
    sequence_file = open("../"+sourcePath)
    def blocks(file, size=65536):
        while True:
            b = file.read(size)
            if not b: break
            yield b

    noLines = sum(bl.count("\n") for bl in blocks(sequence_file))
    print("Total number of lines in file : "+str(noLines))

    for length in range(maxLength, 2, -1):
        print("Counting K-mers of length : "+str(length))
        sequence_file.seek(0)
        time.sleep(0.01)

        with tqdm(total=noLines, position=0, leave=True) as pbar:

            sequenceCount = 0
            sequence = ""
            sequenceIndexes = {}
            nucleotideCounter = 0

            #Reading the file
            if sourceType == 0 :
                for line in sequence_file:
                    pbar.update(1)
                    if sequenceCount < seqLimiter:
                        if line[0] == '>':
                            if sequence != "":
                                nucleotideCounter = countKmers(sequence, sequenceIndexes, nucleotideCounter, length)
                            sequence = ""
                            sequenceCount+=1
                        else:
                            rawLine = line[:len(line)-1]
                            sequence = sequence + rawLine#Stitching the sequence back together line by line
                            # TODO find a way to deal with very long sequences
                if sequence != "":
                    countKmers(sequence, sequenceIndexes, nucleotideCounter, length)
            if sourceType == 1:
                headers = sequence_file.readline().split(',')
                line = sequence_file.readline()
                while line:
                    pbar.update(1)
                    lineElements = line.split(',')
                    sequence = lineElements[-1].lower()
                    nucleotideCounter = countKmers(sequence, sequenceIndexes, nucleotideCounter, length)
                    line = sequence_file.readline()



            sequenceIndexesFilename = indexFilePath + "/sequences" +str(length) + ".json"
            os.makedirs(os.path.dirname(sequenceIndexesFilename), exist_ok=True)
            sequenceIndexesFile = open(sequenceIndexesFilename, "w")
            sequenceIndexesFile.write(json.dumps((nucleotideCounter, sequenceIndexes), indent=4))

    sequence_file.close()


def filterByFrequency(indexFilePath):
    frequent_sequence_indexes = [None]*maxLength
    for length in range(maxLength, 2, -1):
        sequenceCountFile = open(indexFilePath+"/sequences"+str(length)+".json", "rb")
        (nucleotideCount, sequence_indexes) = json.load(sequenceCountFile)
        frequent_sequence_indexes[length - 1] = {}


        for kmer in tqdm(sequence_indexes):
            noOfOccurences = sequence_indexes[kmer]
            if noOfOccurences > nucleotideCount/pow(4, len(kmer)) * minimumFrequencyThreshold:
                # occurences independent of longer frequent seq
                #noOfEffectiveOccurences = sequence_indexes[kmer] - occurenceFreqAsSubSeq(kmer, frequent_sequence_indexes[length:maxLength - 1])
                #if noOfEffectiveOccurences > nucleotideCount/pow(4, len(kmer)) * minimumFrequencyThreshold :  # Removing any non-frequent kner
                    dic = frequent_sequence_indexes[length - 1]
                    dic[kmer] = sequence_indexes[kmer]

        sequenceCountFile.close()
        sequenceCountFile = open(indexFilePath + "/sequences" + str(length) + ".json", "w")
        sequenceCountFile.write(json.dumps((nucleotideCount, frequent_sequence_indexes[length-1]), indent=4))
        sequenceCountFile.close()
















