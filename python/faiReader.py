from itertools import *

def readFAI(file):
# 1       247249719       3       60      61
# 2       242951149       251370554       60      61
# 3       199501827       498370892       60      61
    return [line.split() for line in open(file)]


def readFAIContigOrdering(file):
# 1       247249719       3       60      61
# 2       242951149       251370554       60      61
# 3       199501827       498370892       60      61
    return dict([[rec[0], i] for rec, i in izip(readFAI(file), count())])
