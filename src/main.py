import time
import sys

import kmer_Counting

NEW = 0
ADD = 1

def displayHelp():
    print(
        "AVAILIBLE COMMANDS\n"
        " dictionary <action> <dictionaryName> <filepath>\n"
        " <action> values : 'new', 'add', 'filter'\n"
    )


def Index(mode, name, sourcePath):
    destinationPath = "../sequenceIndexFiles/"+name

    start_time = time.time()
    kmer_Counting.countKmers(mode, destinationPath, sourcePath)
    print("--- %s seconds ---" % (time.time() - start_time))

def Filter(name):
    indexFilePath = "../sequenceIndexFiles/"+name
    start_time = time.time()
    kmer_Counting.filterByFrequency()
    print("--- %s seconds ---" % (time.time() - start_time))

#print(f"Arguments count: {len(sys.argv)}")
#for i, arg in enumerate(sys.argv):
    #print(f"Argument {i:>6}: {arg}")

if (sys.argv[1] == "help"):
    displayHelp()

if (sys.argv[1] == "dictionary" and len(sys.argv)>4):
    if (sys.argv[2] == "new"):
        Index(NEW, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "add"):
        Index(ADD, sys.argv[3], sys.argv[4])
    if (sys.argv[2] == "filter"):
        Filter(sys.argv[3])
