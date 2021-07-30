import sys
import numpy as np
from tqdm import tqdm
import pickle

minimumFrequencyThreshold = 1.0
seqLimiter = 100000
maxLength = 10

#Gets the maximum frequency of this sequence as a subsequence of longer ones
def occurenceFreqAsSubSeq(seq, frequent_seq_indexes):
    maxOccurenceFreq = 0
    for seqMap in frequent_seq_indexes:
        for frequentSeq in seqMap:

            # Is seq a subsequence of frequentSeq, and does it occur more than the current max ?
            if seq in frequentSeq and len(seqMap[frequentSeq]) > maxOccurenceFreq:
                maxOccurenceFreq = len(seqMap[frequentSeq])
    return maxOccurenceFreq

def countKmers(mode, indexFilePath, sourcePath):
    for length in tqdm(range(maxLength, 0, -1)):
        sequence_file = open("../"+sourcePath)

        sequenceCount = 0
        sequence = ""
        sequenceIndexes = {}
        nucleotideCounter = 0

        for line in sequence_file:
            if sequenceCount < seqLimiter:
                if line[0] == '>':
                    if sequence != "":
                        for i in range(len(sequence)-length+1):
                            nucleotideCounter+=1
                            kmer = sequence[i:i+length]
                            if sequenceIndexes.__contains__(kmer):
                                sequenceIndexes[kmer].append(nucleotideCounter)
                            else:
                                sequenceIndexes[kmer] = [nucleotideCounter]
                        #print(len(sequence))
                    #print(line)
                    sequence = ""
                    sequenceCount+=1
                else:
                    rawLine = line[:len(line)-1]
                    sequence = sequence + rawLine#Stitching the sequence back together line by line
                    # TODO find a way to deal with very long sequences
        if sequence != "":
            for i in range(len(sequence) - length + 1):
                nucleotideCounter += 1
                kmer = sequence[i:i + length]
                if sequenceIndexes.__contains__(kmer):
                    sequenceIndexes[kmer].append(nucleotideCounter)
                else:
                    sequenceIndexes[kmer] = [nucleotideCounter]

        sequenceIndexesFilename = indexFilePath + "/sequences" +str(length)
        sequenceIndexesFile = open(sequenceIndexesFilename, "wb")
        pickle.dump((nucleotideCounter, sequenceIndexes), sequenceIndexesFile)
        print(sequenceCount)


def filterByFrequency(indexFilePath):
    frequent_sequence_indexes = [{}]*maxLength
    for length in range(maxLength, 0, -1):
        (nucleotideCount, sequence_indexes) = pickle.load(open(indexFilePath+"/sequences"+str(length), "rb"))

        for kmer in tqdm(sequence_indexes):
            noOfOccurences = len(sequence_indexes[kmer])
            if noOfOccurences/nucleotideCount > pow(0.25, len(kmer)) * minimumFrequencyThreshold:
                noOfEffectiveOccurences = len(sequence_indexes[kmer]) - occurenceFreqAsSubSeq(kmer, frequent_sequence_indexes[length:maxLength - 1])
                if noOfEffectiveOccurences / nucleotideCount > pow(0.25, len(kmer)) * minimumFrequencyThreshold:  # Removing any non-frequent kner
                    frequent_sequence_indexes[length - 1][kmer] = sequence_indexes[kmer]
















