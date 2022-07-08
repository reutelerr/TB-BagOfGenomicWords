import json
import sys

import numpy
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier

from src.utils import timerWrapper


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
    countValues = numpy.array([list(seq_BOW.values()) for seq_BOW in seq_BOWs])
    labels = labelsFile.readlines()
    return (countValues, labels)

@timerWrapper
def trainModel(model, dataset, labels):
    scores = cross_val_score(model, dataset, labels)
    model.fit(dataset, labels)
    return (model, scores)

def trainAndTestKNN(trainingFilePath, trainingLabelsFilePath, testingFilePath, testingLabelsFilePath):
    print("Reading training data")
    (trainingBOWs, trainingLabels) = prepareData(trainingFilePath, trainingLabelsFilePath)
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)
    scaledAndNormalizedTrainingBOWs = sklearn.preprocessing.normalize(scaledTrainingBOWs)

    print("Training model (unscaled)")
    neigh = KNeighborsClassifier(n_neighbors=5)
    (neigh, scores) = trainModel(neigh, trainingBOWs, trainingLabels)


    print("Training model (scaled)")
    neighScaled = KNeighborsClassifier(n_neighbors=5)
    (neighScaled, scoresScaled) = trainModel(neighScaled, scaledTrainingBOWs,trainingLabels)

    print("Training model (scaled and normalized)")
    neighScaledAndNormalized = KNeighborsClassifier(n_neighbors=5)
    (neighScaledAndNormalized, scoresScaledAndNormalized) = trainModel(neighScaledAndNormalized, scaledAndNormalizedTrainingBOWs, trainingLabels)

    print("Reading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.normalize(testingBOWs)
    scaledAndNormalizedTestingBOWs = sklearn.preprocessing.normalize(scaledTestingBOWs)

    print("Testing model")
    testing_score = neigh.score(testingBOWs, testingLabels)
    testing_scaled_score = neighScaled.score(scaledTestingBOWs, testingLabels)
    testing_scaled_normalized_score = neighScaledAndNormalized.score(scaledAndNormalizedTestingBOWs, testingLabels)
    print("Unscaled scores :"+str(scores))
    print("-"+str(testing_score))
    print("Scaled scores :"+str(scoresScaled))
    print("-"+str(testing_scaled_score))
    print("Scaled and Normalized scores :" + str(scoresScaledAndNormalized))
    print("-" + str(testing_scaled_normalized_score))







