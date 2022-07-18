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
    },
    'modelTraining': {
        'gridSearchParams': {'hidden_layer_sizes': [5, 10, 15], 'solver': ['sgd'], 'activation': ['logistic'], 'alpha': 0.1 ** np.arange(1, 5), 'learning_rate': ['adaptive'], 'learning_rate_init': [0.01, 0.001, 0.0001]},
        'scoring': 'accuracy',
        'cvFolds': 5,
    }

}