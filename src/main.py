import time
import sys

import kmer_Counting
import utils
from constants import *
from src.testSequenceInjecter import needlemanWunschInjectedSequence, injectSequence


def displayHelp():
    print(
        "AVAILIBLE COMMANDS\n\n"
        " dictionary <action> <dictionaryName> <sourceFilepath>\n"
        "   - <action> values : 'new', 'add', 'filter'\n\n"
        " makeBOWs <dictionaryName> <sourceFilepath> <outputFilePath>\n\n"
        " convertCSVtoFASTA <csvFilePath> <fastaFilePath>\n\n"
        " findAlignments <sequenceFilePath>\n\n"
        " injectSequence <seqFilePath> <newSeqFilePath> <labelFilePath>"
    )


def ReadSequenceFile(mode, dictionaryPath, sourcePath, outputPath =''):
    fileExtension = sourcePath.split('.')[-1]
    if fileExtension == "csv" : sourceType = CSV
    if fileExtension == "fasta" : sourceType = FASTA
    if fileExtension == "txt" : sourceType = FASTA

    start_time = time.time()
    if mode == makeBOW:
        kmer_Counting.vectorize(dictionaryPath, sourcePath, outputPath, sourceType)
    else:
        kmer_Counting.buildDictionary(mode, dictionaryPath, sourcePath, sourceType)
    print("--- %s seconds ---" % (time.time() - start_time))

def Filter(dicPath):
    start_time = time.time()
    kmer_Counting.filterByFrequency(dicPath)
    print("--- %s seconds ---" % (time.time() - start_time))

if (sys.argv[1] == "help"):
    displayHelp()

if (sys.argv[1] == "dictionary"):
    if (sys.argv[2] == "new"):
        ReadSequenceFile(NEW, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "add"):
        ReadSequenceFile(ADD, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "filter"):
        Filter(sys.argv[3])

if (sys.argv[1] == "makeBOWs"):
    ReadSequenceFile(makeBOW, sys.argv[2], sys.argv[3], sys.argv[4])

if (sys.argv[1] == "convertCSVtoFASTA"):
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

if (sys.argv[1] == "findAlignments"):
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")

if (sys.argv[1] == "injectSequence"):
    injectSequence(sys.argv[2], sys.argv[3], sys.argv[4])
