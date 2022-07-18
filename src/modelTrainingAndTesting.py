import json
import sys
import random

import numpy
import sklearn
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
import warnings


from src.utils import timerWrapper
from src.constants import metaParameters

metaParams = metaParameters.get('modelTraining')
cvFolds = metaParams.get('cvFolds')
gridSearchParams = metaParams.get('gridSearchParams')
scoring = metaParams.get('scoring')


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
    seqFile.close()
    labelsFile.close()
    return (countValues, labels)

@timerWrapper
def trainModel(model, dataset, labels):
    model.fit(dataset, labels)
    return model


def trainAndTestModel(modelType, trainingFilePath, trainingLabelsFilePath, testingFilePath, testingLabelsFilePath, verbose=1):
    random.seed()

    print("Reading training data")
    (trainingBOWs, trainingLabels) = prepareData(trainingFilePath, trainingLabelsFilePath)
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)
    scaledAndNormalizedTrainingBOWs = sklearn.preprocessing.normalize(scaledTrainingBOWs)
    model = None
    gridSearcher = None

    print("Model type : "+modelType)
    if modelType== "KNN":
        n_neighbors = 5
        model = KNeighborsClassifier(n_neighbors=n_neighbors)
        print("n_neighbors : "+str(n_neighbors))
    if modelType== "Perceptron":
        random_state = random.randint(1, 100)
        model = MLPClassifier(max_iter=500)
        gridSearcher = GridSearchCV(model, gridSearchParams, cv=cvFolds, scoring=scoring, verbose=verbose)
        print(metaParams)
        print("random_state : "+str(random_state))

    print("\nTraining model (scaled)")
    (gridSearcher) = trainModel(gridSearcher, scaledTrainingBOWs, trainingLabels)
    print("Scores : "+str(gridSearcher.cv_results_))
    print("Best params : " + str(gridSearcher.best_params_))
    print("Best mean cross-validation score : " + str(max(gridSearcher.cv_results_['mean_test_score'])))

    print("\nReading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.scale(testingBOWs)
    print("Testing model")
    testing_scaled_score = gridSearcher.score(scaledTestingBOWs, testingLabels)
    print("Testing data score : " + str(testing_scaled_score))







