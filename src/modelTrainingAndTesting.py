import json
import os
import sys
import random

import numpy
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from statistics import mean
import warnings

from src import utils
from src.utils import timerWrapper
from src.constants import metaParameters

metaParams = metaParameters.get('modelTraining')
cvFolds = metaParams.get('cvFolds')
scoring = metaParams.get('scoring')

@timerWrapper
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

@timerWrapper
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

@timerWrapper
def trainAndTestModel(modelType, trainingFilePath, trainingLabelsFilePath, testingFilePath, testingLabelsFilePath, randomState=None, resultsFilePath=None, verbose=3):
    random.seed()

    print("Reading training data")
    (trainingBOWs, trainingLabels) = prepareData(trainingFilePath, trainingLabelsFilePath)
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)
    scaledAndNormalizedTrainingBOWs = sklearn.preprocessing.normalize(scaledTrainingBOWs)
    model = None
    gridSearcher = None
    gridSearchParams = metaParams[modelType]['gridSearchParams']

    print("Model type : "+modelType)
    if modelType == "KNN":
        model = KNeighborsClassifier()
    if modelType == "MLP":
        model = MLPClassifier(max_iter=500, random_state=randomState)
    if modelType == "RandomForest":
        model = RandomForestClassifier()
    # TODO
    print(str(gridSearchParams))
    print("random_state : " + str(randomState))
    gridSearcher = GridSearchCV(model, metaParams['RandomForest']['gridSearchParams'], cv=cvFolds, scoring=scoring, verbose=verbose, n_jobs=os.cpu_count())



    print("\nTraining model (scaled)")
    (gridSearcher) = trainModel(gridSearcher, scaledTrainingBOWs, trainingLabels)

    # Getting mean times and scores for each parameter value
    results = gridSearcher.cv_results_
    paramSets = results.get('params')
    meanScores = results.get('mean_test_score')
    meanFitTimes = results.get('mean_fit_time')
    paramStats = {}
    for name, values in gridSearchParams.items():
        paramMeanScores = {}
        if len(values) > 1:
            for value in values:
                scoreValues = []
                fitTimeValues = []
                for i in range(len(paramSets)):
                    if(paramSets[i].get(name) == value):
                        scoreValues.append(meanScores[i])
                        fitTimeValues.append(meanFitTimes[i])
                paramMeanScores[value] = {'meanScore': mean(scoreValues), 'meanFitTime': mean(fitTimeValues)}
        paramStats[name] = paramMeanScores

    #Getting mean scores for each fold
    meanFoldScores = []
    for i in range(0, cvFolds):
        foldScores = results.get('split'+str(i)+'_test_score')
        meanFoldScores.append(mean(foldScores))

    shortResults = {
        'variability': metaParams['sequenceInjection']['variability'],
        'params': metaParams[modelType],
        'paramStats': paramStats,
        'meanFoldScores': meanFoldScores
    }

    if resultsFilePath is not None:
        resultsFile = open(resultsFilePath, 'a')
        resultsFile.write(json.dumps(shortResults, indent=4))

    print("Raw results : "+str(gridSearcher.cv_results_))
    print("\n")
    print("Parameter mean stats: ")
    print(paramStats)
    print("Best params : " + str(gridSearcher.best_params_))
    print("Best mean cross-validation score : " + str(max(gridSearcher.cv_results_['mean_test_score'])))
    print("Refit time for best score : " + str(gridSearcher.refit_time_))

    print("\nReading testing data")
    (testingBOWs, testingLabels) = prepareData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.scale(testingBOWs)
    print("Testing model")
    testing_scaled_score = gridSearcher.score(scaledTestingBOWs, testingLabels)
    print("Testing data score : " + str(testing_scaled_score))
    return testing_scaled_score, gridSearcher.refit_time_







