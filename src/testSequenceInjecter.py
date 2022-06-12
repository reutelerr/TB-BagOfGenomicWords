from Bio.Blast import NCBIWWW
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tqdm import tqdm

from src.utils import blocks

injectedSequence = 'gcaattagatctaatgggacggaggcct'

def needlemanWunschInjectedSequence(fastaFilePath):
    fastafile = open(fastaFilePath)
    noLines = sum(bl.count("\n") for bl in blocks(fastafile))
    fastafile.seek(0)

    maxScore = 0
    bestSequence = ''

    line = fastafile.readline()

    with tqdm(total=noLines, position=0, leave=True) as pbar:
        while(line):
            pbar.update(1)
            if(line == '\n' or line[0] == '>'):
                line = fastafile.readline()
                continue
            sequence = line
            score = pairwise2.align.localxs(injectedSequence, sequence, -0.5, -0.25, score_only=True)
            if(score>maxScore):
                maxScore = score
                bestSequence = sequence
            line = fastafile.readline()

    print(maxScore)
    alignment = pairwise2.align.localxs(injectedSequence, bestSequence, -0.5, -0.25, one_alignment_only=True)
    print(alignment)
    #Best alignement to date : gcaatta--gatctaat-gg-gacggaggc--ct
