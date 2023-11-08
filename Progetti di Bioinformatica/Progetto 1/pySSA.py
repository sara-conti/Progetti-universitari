# to use this script:
# iPython pySSA.py modelFileName tEnd

import modelParser as mP
import numpy as np
import sys


def binCoeff(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def propensities(props, xIn, cIn, leftSideIn):
    leftDim = leftSideIn.shape
    for row in range(leftDim[0]):
        props[row] = cIn[row]
        for col in range(leftDim[1]):
            props[row] *= binCoeff(xIn[col], leftSideIn[row, col])
            # print( " -- ", row, col, " -- ", xIn[col],)
            # print(  leftSideIn[row, col],)
            # print(  binCoeff(xin[col], leftsidein[row, col]))
    # print(props)


def tossTau(aZero):
    tau = (1.0 / aZero) * np.log(1.0 / np.random.random())
    return tau


def tossRule(propensities, aZero):
    aProbs = propensities / aZero
    # print(aProbs)
    rulesSep = aProbs.cumsum()
    mu = np.where(np.random.random() < rulesSep)[0][0]
    return mu


def applyRule(deltaMatrix, iRule, aX):
    aX += deltaMatrix[iRule][:]


def feedStep(aX, dictFeed):
    for spec in dictFeed:
        aX[spec] = dictFeed[spec]


fileModelName = sys.argv[1]
tEnd = float(sys.argv[2])

print("modelling approach: |", mP.modellingApproach(fileModelName), "|")

iRulesLines, iSpecsLines, iFeedsLines = mP.findRulesAndSpecs(fileModelName)
aRules, aConstants = mP.rulesAndConstants(fileModelName, iRulesLines)
aConstants = aConstants.astype(float)
dictAlphabet = mP.rulesToAlphabet(aRules)
print("alphabet:index \t\t", dictAlphabet)
print("constants: ", aConstants)
leftSide, rightSide, delta = mP.rulesToMatrices(aRules, dictAlphabet)
print(leftSide)
print(rightSide)
print(delta)

fileName = "dynamics"
dynFile = open(fileName, 'w')

X = mP.speciesInit(fileModelName, iSpecsLines, dictAlphabet).astype(int)
print(X)
dictFeedSpecs = mP.speciesInFeed(fileModelName, iFeedsLines, dictAlphabet, X)
print("dictFeed:\n", dictFeedSpecs)

aPropensities = np.zeros(len(aRules), dtype=float)
t = 0.0
time = [0.0]
# tEnd = 100.0
while (t < tEnd):
    propensities(aPropensities, X, aConstants, leftSide)
    # print(aPropensities)
    a0 = aPropensities.sum()
    if a0 > 0:
        deltaT = tossTau(a0)
        rule = tossRule(aPropensities, a0)
        # print(deltaT, rule)
        t += deltaT
        time.append(t)
        applyRule(delta, rule, X)
        feedStep(X, dictFeedSpecs)
        # print(t, X)
        dynFile.write(str(t) + '\t')
        # np.savetxt(dynFile, X, newline='\t')
        X.tofile(dynFile, sep="\t")
        dynFile.write('\n')
    else:
        t = tEnd
        time.append(t)
        dynFile.write(str(t) + '\t')
        # np.savetxt(dynFile, X, newline='\t')
        X.tofile(dynFile, sep="\t")
        dynFile.write('\n')

dynFile.close()
