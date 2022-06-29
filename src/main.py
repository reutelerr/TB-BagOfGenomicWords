import time
import sys

import kmer_Counting
import utils
from constants import *
from src.testSequenceInjecter import needlemanWunschInjectedSequence



def displayHelp():
    print(
        "AVAILIBLE COMMANDS\n\n"
        " dictionary <action> <dictionaryName> <sourceFilepath>\n"
        "   - <action> values : 'new', 'add', 'filter'\n\n"
        " vectorize <dictionaryName> <sourceFilepath> <outputFilePath>\n\n"


    )


def ReadSequenceFile(mode, dictionaryName, sourcePath, outputPath =''):
    dictionaryPath = dictionaryName
    fileExtension = sourcePath.split('.')[-1]
    if fileExtension == "csv" : sourceType = CSV
    if fileExtension == "fasta" : sourceType = FASTA
    if fileExtension == "txt" : sourceType = FASTA

    start_time = time.time()
    if mode == VECTORIZE:
        kmer_Counting.vectorize(dictionaryName, sourcePath, outputPath, sourceType)
    else:
        kmer_Counting.buildDictionary(mode, dictionaryPath, sourcePath, sourceType)
    print("--- %s seconds ---" % (time.time() - start_time))

def Filter(name):
    indexFilePath = "sequenceIndexFiles/"+name
    start_time = time.time()
    kmer_Counting.filterByFrequency(indexFilePath)
    print("--- %s seconds ---" % (time.time() - start_time))

#print(f"Arguments count: {len(sys.argv)}")
#for i, arg in enumerate(sys.argv):
    #print(f"Argument {i:>6}: {arg}")

if (sys.argv[1] == "help"):
    displayHelp()

if (sys.argv[1] == "dictionary"):
    if (sys.argv[2] == "new"):
        ReadSequenceFile(NEW, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "add"):
        ReadSequenceFile(ADD, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "filter"):
        Filter(sys.argv[3])

if (sys.argv[1] == "vectorize"):
    ReadSequenceFile(VECTORIZE, sys.argv[2], sys.argv[3], sys.argv[4])

if (sys.argv[1] == "convert"):
    utils.CSVtoFASTA(sys.argv[2], sys.argv[3])

if (sys.argv[1] == "verifySequence"):
    needlemanWunschInjectedSequence(sys.argv[2], "proximateSequenceScores.txt")
