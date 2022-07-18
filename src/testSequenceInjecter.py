import math
import os
import random

from Bio.Blast import NCBIWWW
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tqdm import tqdm
from matplotlib import pyplot

from src.constants import CSV, nucleotideValues, metaParameters
from src.utils import blocks

def needlemanWunschInjectedSequence(fastaFilePath, outputPath  = None):
    fastafile = open(fastaFilePath)
    noLines = sum(bl.count("\n") for bl in blocks(fastafile))
    fastafile.seek(0)

    maxScore = 0
    bestSequence = ''

    line = fastafile.readline()
    scores = []
    sequenceCount = 0
    limiter = -1

    with tqdm(total=noLines, position=0, leave=True) as pbar:
        while(line and sequenceCount != limiter):
            pbar.update(1)
            if(line == '\n' or line[0] == '>'):
                line = fastafile.readline()
                continue
            sequenceCount += 1
            sequence = line
            score = pairwise2.align.localxs(injectedSequence, sequence, -0.5, -0.25, score_only=True)
            scores.append(score)
            if(score>maxScore):
                maxScore = score
                bestSequence = sequence
            line = fastafile.readline()

    print(maxScore)
    alignment = pairwise2.align.localxs(injectedSequence, bestSequence, -0.5, -0.25, one_alignment_only=True)
    print(alignment)
    if(outputPath != None):
        outputFile = open(outputPath, "w")
        outputFile.write("Proximity scores to test sequence 'gcaattagatctaatgggacggaggcct'\n")
        scores.sort()
        outputFile.write(str(scores))
        pyplot.boxplot(scores)
        pyplot.show()
    #Best alignement to date : gcaatta--gatctaat-gg-gacggaggc--ct with a score of 22.5 / 28


def injectSequence(seqFilePath, newFilePath, labelFilePath, sourceType=CSV):

    metaParams = metaParameters.get('sequenceInjection')
    injectedSequence = metaParams.get('injectedSequence')
    injectionRate = metaParams.get('injectionRate')
    variability = metaParams.get('variability')

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
                if variability>0:
                    for i in range(len(injectedSequence)):
                        if random.random()<variability:
                            modifiedSequence += nucleotideValues[random.randint(0, len(nucleotideValues)-1)]
                        else:
                            modifiedSequence += injectedSequence[i]
                else:
                    modifiedSequence = injectedSequence

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
