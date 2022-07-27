import math
import os
import random

from Bio.Blast import NCBIWWW
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tqdm import tqdm
from matplotlib import pyplot

from src.constants import *
from src.utils import blocks, timerWrapper

maxKmerLength = metaParameters['vectorization']['maxKmerLength']
minKmerLength = metaParameters['vectorization']['minKmerLength']

@timerWrapper
def needlemanWunschInjectedSequence(seqFilePath, outputPath=None, sourceType=CSV):
    injectedSequence = metaParameters['sequenceInjection']['injectedSequence']
    seqFile = open(seqFilePath)
    noLines = sum(bl.count("\n") for bl in blocks(seqFile))
    seqFile.seek(0)

    maxScore = 0
    bestSequence = ''

    line = seqFile.readline()
    scores = []
    sequenceCount = 0
    limiter = -1

    if sourceType == CSV:
        headers = seqFile.readline().split(',')
        line = seqFile.readline()

    with tqdm(total=noLines, position=0, leave=True) as pbar:
        while line and sequenceCount != limiter:
            pbar.update(1)
            if sourceType == CSV:
                lineElements = line.split(',')
                id = lineElements[0]
                sequence = lineElements[-1].lower()[:len(lineElements[-1]) - 1]
            if sourceType == FASTA:
                if line != '\n' and line[0] != '>':
                    rawLine = line[:len(line) - 1]
                    sequence = sequence + rawLine  # Stitching the sequence back together line by line
                    continue
                sequenceCount += 1
            score = pairwise2.align.localxs(injectedSequence, sequence, -3.0, -0.0, score_only=True, penalize_end_gaps=False)
            scores.append(score)
            if score>maxScore:
                maxScore = score
                bestSequence = sequence
            line = seqFile.readline()

    print('Best score: ' + str(maxScore) + '/' + str(len(injectedSequence)))
    alignment = pairwise2.align.localxs(injectedSequence, bestSequence, -3.0, -0.0, one_alignment_only=True, penalize_end_gaps=False)
    print('Best alignment:' + str(alignment))
    if outputPath != None:
        outputFile = open(outputPath, "w")
        outputFile.write('Proximity scores to test sequence "' + injectedSequence + '"\n')
        print("Sorting and saving scores...")
        scores.sort()
        outputFile.write(str(scores))
    pyplot.boxplot(scores)
    pyplot.title('Best Needleman-Wunsch alignement scores')
    pyplot.show()


def kmerLoss(substitutionIndexes, sequence=metaParameters['sequenceInjection']['injectedSequence']):
    kmerTotal = ((len(sequence)+1)*(maxKmerLength+1-minKmerLength))+minKmerLength-((maxKmerLength*(maxKmerLength+1))/2)#Number of kmers the unaltered sequence is composed of
    #TODO


@timerWrapper
def injectSequence(seqFilePath, newFilePath, labelFilePath, variabilityType='randomPosition', sourceType=CSV):

    metaParams = metaParameters.get('sequenceInjection')
    injectedSequence = metaParams.get('injectedSequence')
    injectionRate = metaParams.get('injectionRate')
    variability = metaParams.get('variability')
    fixedVariabilityIndexes = metaParams.get('fixedVariabilityIndexes')

    random.seed()
    seqFile = open(seqFilePath, "r")
    os.makedirs(os.path.dirname(newFilePath), exist_ok=True)
    newFile = open(newFilePath, "w")
    labelFile = open(labelFilePath, "w")

    if sourceType == CSV:
        headerLine = seqFile.readline()
        headers = headerLine.split(',')
        newFile.write(headerLine)
        line = seqFile.readline()
        while line:
            if random.random()<injectionRate :
                modifiedSequence = ""
                if variabilityType == 'randomPosition':
                    if variability>0:
                        for i in range(len(injectedSequence)):
                            if random.random()<variability/len(injectedSequence):
                                otherNucleotideValues = nucleotideValues.copy()
                                otherNucleotideValues.remove(injectedSequence[i])
                                modifiedSequence += otherNucleotideValues[random.randint(0, len(otherNucleotideValues)-1)]
                            else:
                                modifiedSequence += injectedSequence[i]
                    else:
                        modifiedSequence = injectedSequence
                if variabilityType == 'fixedPosition':
                    for i in range(len(injectedSequence)):
                        if i in fixedVariabilityIndexes:
                            otherNucleotideValues = nucleotideValues.copy()
                            otherNucleotideValues.remove(injectedSequence[i])
                            modifiedSequence += nucleotideValues[random.randint(0, len(nucleotideValues)-1)]
                        else:
                            modifiedSequence += injectedSequence[i]



                label = 1
                seqStart = line.rindex(',')+1
                insertIndex = random.randint(seqStart, len(line)-1)
                if not line[0:insertIndex]:
                    raise Exception('must not insert at beginning of line')
                if not line[insertIndex:len(line)]:
                    raise Exception('must not insert AFTER line')
                newFile.write(line[0:insertIndex]+modifiedSequence+line[insertIndex:len(line)])

            else:
                label = 0
                newFile.write(line)
            labelFile.write(str(label)+'\n')
            line = seqFile.readline()

    seqFile.close()
    newFile.close()
    labelFile.close()
