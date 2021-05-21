import time

import kmer_Counting

def Index():
    start_time = time.time()
    kmer_Counting.countKmers()
    print("--- %s seconds ---" % (time.time() - start_time))

def Filter():
    start_time = time.time()
    kmer_Counting.filterByFrequency()
    print("--- %s seconds ---" % (time.time() - start_time))

Filter()