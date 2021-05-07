import sys
import numpy as np

minimumOccurenceThreshold = 1.0
seqLimiter = 10


def occurenceFreqAsSubSeq(seq, freq_seq_indexes):
    maxOccurenceFreq = 0
    for seq_len in freq_seq_indexes:
        for freq_seq in seq_len:
            if seq in freq_seq and len(freq_seq_indexes[freq_seq]) > maxOccurenceFreq:
                maxOccurenceFreq = len(freq_seq_indexes[seq])
    return maxOccurenceFreq

def countKmers(maxLength):
    sequence_file = open("../sequences.fasta")
    frequent_sequence_indexes = [{}]*maxLength
    for length in range(maxLength, 0, -1):
        nucleotideCounter = 0
        sequenceCount = 0
        sequence = ""
        sequence_indexes = {}
        for line in sequence_file:
            if sequenceCount < seqLimiter:
                if line[0] == '>':
                    if sequence != "":
                        for i in range(len(sequence)-length+1):
                            nucleotideCounter+=1
                            kmer = sequence[i:i+length]
                            if sequence_indexes.__contains__(kmer):
                                sequence_indexes[kmer].append(nucleotideCounter)
                            else:
                                sequence_indexes[kmer] = [nucleotideCounter]
                        print(len(sequence))
                    print(line)
                    sequence = ""
                    sequenceCount+=1
                else:
                    sequence = sequence + line#Stitching the sequence together line by line
                    # TODO find a way to deal with very long sequences

        for kmer in sequence_indexes:
            noOfEffectiveOccurences = len(sequence_indexes[kmer]) - occurenceFreqAsSubSeq(kmer, frequent_sequence_indexes[length:maxLength-1])
            if noOfEffectiveOccurences/nucleotideCounter > pow(0.25, len(kmer)) * minimumOccurenceThreshold:  # Removing any non-frequent kner
                frequent_sequence_indexes[length-1][kmer] = sequence_indexes[kmer]


countKmers(10)











