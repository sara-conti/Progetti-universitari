# to use this script:
# iPython pyODE.py modelFileName Volume tEnd

import modelParser as mP
import numpy as np
import scipy.integrate as spi
import sys
import matplotlib.pyplot as plt


fileModelName = sys.argv[1]
print(fileModelName)
Vol = float(sys.argv[2])
tEnd = float(sys.argv[3])

iRulesLines, iSpecsLines, iFeedsLines = mP.findRulesAndSpecs(fileModelName)
aRules, aConstants = mP.rulesAndConstants(fileModelName, iRulesLines)
dictAlphabet = mP.rulesToAlphabet(aRules)
print("alphabet:index \t\t", dictAlphabet)
leftSide, rightSide, delta = mP.rulesToMatrices(aRules, dictAlphabet)
leftSide = leftSide.astype(float)
rightSide = rightSide.astype(float)
delta = delta.astype(float)
print("L\n", leftSide, "\n")
print("R\n", rightSide, "\n")
print("D\n", delta, "\n")
aConstants = aConstants.astype(float)
print("c ", aConstants)

X = mP.speciesInit(fileModelName, iSpecsLines, dictAlphabet).astype(float)
print("X_in: ", X)
dictFeedSpecs = mP.speciesInFeed(fileModelName, iFeedsLines, dictAlphabet, X)
print("dictFeed:\n", dictFeedSpecs)

approachIn = mP.modellingApproach(fileModelName)
if 'Stoc' in approachIn:
    print("Stoch->Det")
    X, aConstants = mP.stoch2Det(Vol, X, leftSide, aConstants)
print("X: ", X)
print("K: ", aConstants)


def myCauchy(x, t):
    # dx = np.zeros(len(dictAlphabet) - 1, dtype=float)
    terms = np.product(np.power(x, leftSide), axis=1)
    teta = delta.T * aConstants
    dx = teta.dot(terms)
    for spec in dictFeedSpecs:
        dx[spec] = 0.0
    # print(dx)
    return dx


print("\n")
# terms = np.product(np.power(X, leftSide), axis=1)
# print("terms: ", terms)
# teta = delta.T * aConstants
# print("teta:\n", teta)
# xPunto = teta.dot(terms)
# for spec in dictFeedSpecs:
#     xPunto[spec] = 0.0

time = np.linspace(0.0, tEnd, 10000)
yInit = X
y = spi.odeint(myCauchy, yInit, time, atol=1e-9)

print(y[-1] * Vol * 6.022e23)

plt.grid()
plt.xlabel('t', fontsize=22)
plt.ylabel('concentrazione (mg)', fontsize=22)
for spec in dictAlphabet:
    lista = [dictAlphabet['lAPAP'], dictAlphabet['gAPAP'], dictAlphabet['pAPAP'], dictAlphabet['tAPAP'], dictAlphabet['uAPAP']]
    plt.plot(time, y[:, lista] * Vol * 6.022e23, '-', label=spec)
import pandas as pd
#l = dictAlphabet['lAPAP'] * Vol * 6.022e23
#g = dictAlphabet['gAPAP'] * Vol * 6.022e23
#p = dictAlphabet['pAPAP'] * Vol * 6.022e23
#t = dictAlphabet['tAPAP'] * Vol * 6.022e23
df=pd.DataFrame(list(zip(time, y[:, dictAlphabet['lAPAP']] * Vol * 6.022e23,
y[:, dictAlphabet['gAPAP']] * Vol * 6.022e23,
y[:, dictAlphabet['pAPAP']] * Vol * 6.022e23,
y[:, dictAlphabet['tAPAP']] * Vol * 6.022e23,
y[:, dictAlphabet['uAPAP']] * Vol * 6.022e23))
,columns=['time','lAPAP','gAPAP', 'pAPAP','tAPAP','uAPAP'])
#y[:, dictAlphabet['gAPAP']] * Vol * 6.022e23,
#y[:, dictAlphabet['pAPAP']] * Vol * 6.022e23,
#y[:, dictAlphabet['tAPAP']] * Vol * 6.022e23),
#columns=['time','lAPAP', 'gAPAP', 'pAPAP','tAPAP']))

#print(df)
df.to_csv("600.csv",index=False)
plt.legend(['lAPAP','gAPAP','pAPAP','tAPAP', "uAPAP"])
fileName = fileModelName + 'ODEpy.pdf'
plt.savefig(fileName, format='pdf')
plt.show()
