import json
import sys

import sklearn
from sklearn.model_selection import cross_val_score
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
        print('\rSequences loaded : [%d]' % sequencesLoaded, end="")
        sys.stdout.flush()
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
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)

    print("Training model (unscaled)")
    neigh = KNeighborsClassifier(n_neighbors=100)
    scores = cross_val_score(neigh, trainingBOWs, trainingLabels)
    neigh.fit(trainingBOWs, trainingLabels)

    print("Training model (scaled)")
    neighScaled = KNeighborsClassifier(n_neighbors=100)
    scoresScaled = cross_val_score(neighScaled, scaledTrainingBOWs, trainingLabels)
    neighScaled.fit(scaledTrainingBOWs, trainingLabels)

    print("Reading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.normalize(testingBOWs)



    print("Testing model")
    testing_score = neigh.score(testingBOWs, testingLabels)
    testing_scaled_score = neighScaled.score(scaledTestingBOWs, testingLabels)
    print("Unscaled scores :"+str(scores))
    print("-"+str(testing_score))
    print("Scaled scores :"+str(scoresScaled))
    print("-"+str(testing_scaled_score))







