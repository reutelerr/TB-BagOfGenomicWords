import numpy as np

NEW = 0
ADD = 1
makeBOW = 2

FASTA = 0
CSV = 1

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
        'variability': 1, #Number of substituted nucleotides in the sequence (average in the case of random variability)
        'fixedVariabilityIndexes': [3, 15]
    },
    'modelTraining': {
        'MLP': {
            'gridSearchParams': {
                'hidden_layer_sizes': [8, 10, 12],
                'solver': ['sgd'],
                'activation': ['logistic'],
                'alpha': list(0.1 ** np.arange(2, 3)),
                'learning_rate': ['adaptive'],
                'learning_rate_init': [0.01]
            }
        },
        'RandomForest': {
            'gridSearchParams': {
                'n_estimators': [200, 400],
                'min_samples_split': [2, 3, 5],
                'min_samples_leaf': [3, 5, 10],
                'max_features': ['sqrt']
            }
        },
        'scoring': 'accuracy',
        'cvFolds': 5,
    }
}

sourcePath = "sequences/PhageWholeDNASeq/PhageWholeDNASeq.csv"
injectedSeqPath = "sequences/InjectedPhageWholeDNASeq/InjectedPhageWholeDNASeq.csv"
labelsPath = "sequences/InjectedPhageWholeDNASeq/labels.txt"
trainingSeqPath = "sequences/InjectedPhageWholeDNASeq/Training/seq.csv"
trainingLabelsPath = "sequences/InjectedPhageWholeDNASeq/Training/labels.txt"
testingSeqPath = "sequences/InjectedPhageWholeDNASeq/Testing/seq.csv"
testingLabelsPath = "sequences/InjectedPhageWholeDNASeq/Testing/labels.txt"
trainingBOWsPath = "sequences/InjectedPhageWholeDNASeq/Training/Seq_BOW.json"
testingBOWsPath = "sequences/InjectedPhageWholeDNASeq/Testing/Seq_BOW.json"
dictionaryPath = "kmerDictionaries/dic10"

