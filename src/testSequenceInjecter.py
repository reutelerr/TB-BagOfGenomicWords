from Bio.Blast import NCBIWWW
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tqdm import tqdm
from matplotlib import pyplot

from src.constants import CSV
from src.utils import blocks

injectedSequence = 'gcaattagatctaatgggacggaggcct'

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

def injectSequence(seq, seqFile, rate = 0.5, sourceType = CSV):
    if sourceType == CSV:
        for line in seqFile:
            return None
            #if rand()<rate :
                #seqStart = line.lastindexof(',')+1
                #seqLength = line.length - seqStart
                #insertIndex = seqStart + Math.floor(rand() * seqLength)
