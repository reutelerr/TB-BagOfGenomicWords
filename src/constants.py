NEW = 0
ADD = 1
makeBOW = 2

FASTA = 0
CSV = 1

import numpy as np

nucleotideValues = ['a', 'c', 'g', 't']

metaParameters = {
    'vectorization': {
        'minimumFrequencyThreshold': 0.5,
        'maxKmerLength': 10,
        'filterThresholdLengthModifier': 1.2,
    },
    'sequenceInjection': {
        'injectedSequence': 'gcaattagatctaatgggacggaggcct',
        'injectionRate': 0.5,
        'variability': 0.0,
        'fixedVariabilityIndexes': [3, 15]
    },
    'modelTraining': {
        'gridSearchParams': {'hidden_layer_sizes': [8, 10, 12], 'solver': ['sgd'], 'activation': ['logistic'], 'alpha': 0.1 ** np.arange(2, 3), 'learning_rate': ['adaptive'], 'learning_rate_init': [0.01]},
        'scoring': 'accuracy',
        'cvFolds': 10,
    }

}