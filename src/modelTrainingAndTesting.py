import json
import os
import sys
import random

import numpy
import sklearn
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix, roc_curve, auc
from statistics import mean
import warnings

from src import utils
from src.utils import timerWrapper
from src.constants import metaParameters

modelParams = metaParameters.get('modelTraining')
cvFolds = modelParams.get('cvFolds')
scoring = modelParams.get('scoring')

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
def readData(seqFilePath, labelsFilePath):
    seqFile = open(seqFilePath)
    labelsFile = open(labelsFilePath)
    data = readBOWFile(seqFile)
    seq_BOWs = [counts for (index, nulceotideCount, vecLength, counts) in data]
    countValues = numpy.array([list(seq_BOW.values()) for seq_BOW in seq_BOWs])
    labels = [int(labelLine[0]) for labelLine in labelsFile.readlines()]
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
    (trainingBOWs, trainingLabels) = readData(trainingFilePath, trainingLabelsFilePath)
    scaler = sklearn.preprocessing.StandardScaler().fit(trainingBOWs)
    scaledTrainingBOWs = scaler.transform(trainingBOWs)
    scaledAndNormalizedTrainingBOWs = sklearn.preprocessing.normalize(scaledTrainingBOWs)
    model = None
    gridSearcher = None
    gridSearchParams = modelParams[modelType]['gridSearchParams']

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
    gridSearcher = GridSearchCV(model, modelParams[modelType]['gridSearchParams'], cv=cvFolds, scoring=scoring, verbose=verbose, n_jobs=os.cpu_count())

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
                paramMeanScores[str(value)] = {'meanScore': mean(scoreValues), 'meanFitTime': mean(fitTimeValues)}
        paramStats[name] = paramMeanScores

    #Getting mean scores for each fold
    meanFoldScores = []
    for i in range(0, cvFolds):
        foldScores = results.get('split'+str(i)+'_test_score')
        meanFoldScores.append(mean(foldScores))



    print("Raw results : "+str(gridSearcher.cv_results_))
    print("\n")
    print("Parameter mean stats: ")
    print(paramStats)
    print("Best params : " + str(gridSearcher.best_params_))
    print("Best mean cross-validation score : " + str(max(gridSearcher.cv_results_['mean_test_score'])))
    print("Refit time for best score : " + str(gridSearcher.refit_time_))

    print("\nReading testing data")
    (testingBOWs, testingLabels) = readData(testingFilePath, testingLabelsFilePath)
    scaledTestingBOWs = sklearn.preprocessing.scale(testingBOWs)
    print("Testing model")
    testing_scaled_score = gridSearcher.score(scaledTestingBOWs, testingLabels)
    confusionMatrix = confusion_matrix(testingLabels, gridSearcher.predict(scaledTestingBOWs))
    y_pred_proba = gridSearcher.predict_proba(testingBOWs)

    fpr, tpr, _ = roc_curve(testingLabels, y_pred_proba[:, 1])
    roc_auc = auc(fpr, tpr)
    print("roc_auc score = "+str(roc_auc))

    plt.figure()
    lw = 2
    plt.plot(
        fpr,
        tpr,
        color="darkorange",
        lw=lw,
        label="ROC curve (area = %0.2f)" % roc_auc,
    )
    plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Receiver operating characteristic")
    plt.legend(loc="lower right")
    plt.show()

    print("Testing data score : " + str(testing_scaled_score))
    print("Testing data confusion matrix : " + str(confusionMatrix))

    shortResults = {
            'injectionrate': metaParameters['sequenceInjection']['injectionRate'],
            'variability': metaParameters['sequenceInjection']['variability'],
            'params': modelParams[modelType],
            'paramStats': paramStats,
            'meanFoldScores': meanFoldScores,
            'bestConfusionMatrix': str(confusion_matrix),
        }

    if resultsFilePath is not None:
        resultsFile = open(resultsFilePath, 'a')
        resultsFile.write(json.dumps(shortResults, indent=4))

    return testing_scaled_score, gridSearcher.refit_time_, gridSearcher







