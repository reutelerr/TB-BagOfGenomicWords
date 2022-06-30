import json
from sklearn.neighbors import KNeighborsClassifier

def readBOWFile(seqFile):
    data = []
    str=""
    line = seqFile.readline()
    sequencesLoaded = 0
    while line:
        while line != "===\n":
            rawLine = line[:len(line) - 1]
            str = str + rawLine
            line = seqFile.readline()
        data.append(json.loads(str))
        sequencesLoaded += 1
        print('\rSequences loaded : [%d]\n' % sequencesLoaded, end="")
        str=""
        line = seqFile.readline()
    return data


def prepareData(seqFilePath, labelsFilePath):
    seqFile = open(seqFilePath)
    labelsFile = open(labelsFilePath)
    data = readBOWFile(seqFile)
    seq_BOWs = [counts for (index, nulceotideCount, vecLength, counts) in data]
    countValues = [list(seq_BOW.values()) for seq_BOW in seq_BOWs]
    labels = labelsFile.readlines()
    return (countValues, labels)

def trainAndTestKNN(trainingFilePath, trainingLabelsFilePath, testingFilePath, testingLabelsFilePath):
    print("Reading training data")
    (trainingBOWs, trainingLabels) = prepareData(trainingFilePath, trainingLabelsFilePath)

    print("Training model")
    neigh = KNeighborsClassifier(n_neighbors=10)
    neigh.fit(trainingBOWs, trainingLabels)

    print("Reading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)

    print("Testing model")
    scores = neigh.score(testingBOWs, testingLabels)
    print(scores)







