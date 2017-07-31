import argparse

import numpy as np
from sklearn.metrics import *

parser = argparse.ArgumentParser(description="train neural net")
parser.add_argument('-p', '--predict', help="give directories with files")
parser.add_argument('-t', '--truth', help="give directories with files")
paths = parser.parse_args()
paths = vars(paths)
predictinput = paths['predict']
truthinput = paths['truth']

predictarray = np.load(predictinput)
trutharray = np.load(truthinput)
newpredictarray = []

print list(predictarray)[:1000]
print list(trutharray)[:1000]

print predictarray.shape
print trutharray.shape

print recall_score(predictarray, trutharray)
print precision_score(predictarray, trutharray)
print f1_score(predictarray, trutharray)
