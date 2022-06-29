# TB-BagOfGenomicWords

## K-mer Counting

To use a "bag of words" approach, we build a dictionary of all the K-mers we find, 
as long as they have a given minimum frequency. Thus, only kmers that appear often will be considered

We can then convert sequences into "Bags of Words" usable in ML models, counting the occurences of each "word" of the previously built dictionary

###Commands :

Use "dictionary new <dictionaryName> <sourceFile>" to create a new dictionary

Use "dictionary filter <dictionaryName>" to filter out less frequent K-mers

Use "makeBOWs <dictionaryName> <sourceFilepath> <outputFilePath>" to convert sequences in text form (FASTA, CSV, ...) into Bags of Words in json format


Use "findAlignments <sequenceFilePath>" to get the scores of the best alignments in a sequence file. Useful to verify if a sequence is naturally occuring in our data. Said sequence is for now simply a constant in testSequenceInjecter.py


