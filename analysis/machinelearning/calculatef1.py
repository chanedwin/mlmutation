import argparse

import numpy as np
from sklearn.metrics import *

#python wrapper script that takes in 2 vcf files, one vcf files of predicted mutations and another ground truth
#and calculates the f1 score

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

print recall_score(predictarray, trutharray)
print precision_score(predictarray, trutharray)
print f1_score(predictarray, trutharray)
