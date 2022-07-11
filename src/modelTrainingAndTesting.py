import json
import sys
import random

import numpy
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier

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
    print('\n')
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

def trainAndTestKNN(modelType, trainingFilePath, trainingLabelsFilePath, testingFilePath, testingLabelsFilePath):
    random.seed()

    print("Reading training data")
    (trainingBOWs, trainingLabels) = prepareData(trainingFilePath, trainingLabelsFilePath)
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)
    scaledAndNormalizedTrainingBOWs = sklearn.preprocessing.normalize(scaledTrainingBOWs)
    model = None
    modelScaled = None
    modelScaledAndNormalized = None


    print("Model type : "+modelType)
    if(modelType=="KNN"):
        n_neighbors = 5
        model = KNeighborsClassifier(n_neighbors=n_neighbors)
        modelScaled = KNeighborsClassifier(n_neighbors=n_neighbors)
        modelScaledAndNormalized = KNeighborsClassifier(n_neighbors=n_neighbors)
        print("n_neighbors : "+str(n_neighbors))
    if(modelType=="Perceptron"):
        hidden_layer_sizes = (10)
        random_state = random.randint(1, 100)
        model = MLPClassifier(solver='sgd', hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        modelScaled = MLPClassifier(solver='sgd', hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        modelScaledAndNormalized = MLPClassifier(solver='sgd', hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        print("hidden_layer_sizes : "+str(hidden_layer_sizes))
        print("random_state : "+str(random_state))

    print("Training model (unscaled)")
    (model, scores) = trainModel(model, trainingBOWs, trainingLabels)

    print("Training model (scaled)")
    (modelScaled, scoresScaled) = trainModel(modelScaled, scaledTrainingBOWs, trainingLabels)

    print("Training model (scaled and normalized)")
    (modelScaledAndNormalized, scoresScaledAndNormalized) = trainModel(modelScaledAndNormalized, scaledAndNormalizedTrainingBOWs, trainingLabels)

    print("Reading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.scale(testingBOWs)
    scaledAndNormalizedTestingBOWs = sklearn.preprocessing.normalize(scaledTestingBOWs)

    print("Testing model")
    testing_score = model.score(testingBOWs, testingLabels)
    testing_scaled_score = modelScaled.score(scaledTestingBOWs, testingLabels)
    testing_scaled_normalized_score = modelScaledAndNormalized.score(scaledAndNormalizedTestingBOWs, testingLabels)
    print("Unscaled scores :"+str(scores))
    print('->'+str(testing_score))
    print("Scaled scores :"+str(scoresScaled))
    print('->'+str(testing_scaled_score))
    print("Scaled and Normalized scores :"+str(scoresScaledAndNormalized))
    print('->'+str(testing_scaled_normalized_score))







