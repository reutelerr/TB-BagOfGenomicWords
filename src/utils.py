import os

from tqdm import tqdm

def blocks(file, size=65536):
    while True:
        b = file.read(size)
        if not b: break
        yield b

def CSVtoFASTA(csvfilePath, fastaFilePath):
    csvFile = open(csvfilePath)
    fastaFile = open(fastaFilePath, 'w')

    noLines = sum(bl.count("\n") for bl in blocks(csvFile))
    csvFile.seek(0)
    with tqdm(total=noLines, position=0, leave=True) as pbar:

        headers = csvFile.readline().split(',')
        line = csvFile.readline()
        while line:
            pbar.update(1)
            lineElements = line.split(',')
            head = lineElements[3]
            startIndex = head.index('>')
            head = head[startIndex:]+' '+lineElements[4]
            head = head.replace('\"', '')
            sequence = lineElements[-1].lower()

            fastaFile.write(head+'\n')
            fastaFile.write(sequence)
            fastaFile.write('\n')
            line = csvFile.readline()

    fastaFile.close()
    csvFile.close()
