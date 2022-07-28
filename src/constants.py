import random
import numpy as np

NEW = 0
ADD = 1
makeBOW = 2

FASTA = 0
CSV = 1

nucleotideValues = ['a', 'c', 'g', 't']

randomStateSeed = random.randint(0, 100)

metaParameters = {
    'vectorization': {
        'minimumFrequencyThreshold': 0.5,
        'minKmerLength': 3,
        'maxKmerLength': 10,
        'filterThresholdLengthModifier': 1.1,
    },
    'sequenceInjection': {
        'injectedSequence': 'gcaattagatctaatgggacggaggcct',
        'injectionRate': 0.5,
        'variabilityType': 'fixedPosition',
        'variability': 2, #Number of substituted nucleotides in the sequence (average in the case of random variability)
        'fixedVariabilityIndexes': [3, 17],
        'noOfConfigurations': 10
    },
    'modelTraining': {
        'KNN': {
            'gridSearchParams': {
                'n_neighbors': [5, 10, 20, 50, 100, 500, 1000],
                'weights': ['uniform', 'distance']
            }
        },
        'MLP': {
            'gridSearchParams': {
                'hidden_layer_sizes': [5, 10, 20, 50, 100, (50, 5)],
                'solver': ['sgd'],
                'activation': ['logistic'],
                'alpha': list(0.1 ** np.arange(2, 5)),
                'learning_rate': ['adaptive'],
                'learning_rate_init': [0.01]
            }
        },
        'RandomForest': {
            'gridSearchParams': {
                'n_estimators': [400],
                'min_samples_split': [5],
                'min_samples_leaf': [1],
                'max_features': ['sqrt']
            }
        },
        'scoring': 'f1',
        'cvFolds': 5,
    },
    'plotData': [6.6546831130981445, 16.577001333236694, 10.10099983215332, 19.050995588302612, 8.281992197036743, 17.59199833869934, 16.979514122009277, 16.942002773284912, 17.967002153396606, 9.079991579055786, 18.29199743270874, 16.3249990940094, 10.277997732162476, 10.356998920440674, 21.10798740386963, 17.6329984664917, 9.689000368118286, 11.276999473571777, 9.921990633010864, 24.95999813079834, 23.384995937347412, 12.328008890151978, 18.807042360305786, 19.929998636245728, 12.997000217437744, 12.612975597381592, 10.054999828338623, 22.43999934196472, 12.776999473571777, 12.452004671096802, 12.332996606826782]
}

phageSeqFolder = "sequences/PhageWholeDNASeq/"
injectedPhageSeqFolder = "sequences/InjectedPhageWholeDNASeq/"
#bacteriaSeqPath = "sequences/PseudomonasWholeDNASeq/PseudomonasWholeDNASeq.csv"

phageSeqPath = phageSeqFolder + "PhageWholeDNASeq.csv"
injectedPhageSeqPath = injectedPhageSeqFolder + "InjectedPhageWholeDNASeq.csv"

phageKmerDictionaryPath = "phageKmerDictionaries/dic10"
bacteriaKmerDictionaryPath = "bacteriaKmerDictionaries/dic10"


