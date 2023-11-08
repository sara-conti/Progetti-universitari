import numpy as np
import re
import scipy.special as sps


def uniquefy(extList):
    seen = set()
    seen_add = seen.add
    return [x for x in extList if not (x in seen or seen_add(x))]


def termsToRule(terms, matrix, rule, dictAlphabetIn):
    for term in terms:
        # print( 'term |', term, '|' )
        if term is not '':
            parti = re.split('\*', term)
            # print( "--- ", parti )
            if re.match('\d', parti[0]) is not None:
                # print( rule, " : ", dictAlphabetIn[parti[1]] )
                matrix[rule][dictAlphabetIn[parti[1]]] = int(parti[0])
            else:
                # print( rule, dictAlphabetIn[parti[0]] )
                matrix[rule][dictAlphabetIn[parti[0]]] = 1


def findRulesAndSpecs(fileModelName):
    with open(fileModelName, 'r') as fileModel:
        linesToParse = fileModel.readlines()
    fileModel.close

    i = 0
    iLineRulesStart = 0
    iLineRulesStop = 0
    iLineSpecsStart = 0
    iLineSpecsStop = 0
    iLineFeedsStart = 0
    iLineFeedsStop = 0
    for line in linesToParse:
        if 'Rules>' in line:
            iLineRulesStart = i
        if 'X0>' in line:
            iLineSpecsStart = i
        if 'XFeed>' in line:
            iLineFeedsStart = i
        if '<Rules' in line:
            iLineRulesStop = i
        if '<X0' in line:
            iLineSpecsStop = i
        if '<XFeed' in line:
            iLineFeedsStop = i
        i += 1
    # print( iLineRulesStart, iLineRulesStop )
    rulesLines = np.array([iLineRulesStart, iLineRulesStop], dtype=int)
    specsLines = np.array([iLineSpecsStart, iLineSpecsStop], dtype=int)
    feedsLines = np.array([iLineFeedsStart, iLineFeedsStop], dtype=int)
    return rulesLines, specsLines, feedsLines


def rulesAndConstants(fileModelName, rulesLines):
    with open(fileModelName, 'r') as fileModel:
        linesToParse = fileModel.readlines()
    fileModel.close

    rules = []
    constants = []
    for line in linesToParse[rulesLines[0] + 1:rulesLines[1]]:
        tokens = line.strip().split(';')
        rules.append(tokens[0].replace(' ', ''))
        constants.append(tokens[1].replace(' ', ''))

    # nRules = len(rules)
    print( "rules \t\t\t", rules )
    fConstants = np.array(constants)
    print( "constants \t\t", fConstants )
    return rules, fConstants


def rulesToAlphabet(rules):
    alphabet = []
    for line in rules:
        elementi = re.split('[\+(\-\>)]', line)
        for elemento in elementi:
            parti = re.split('\*', elemento)
            for el in parti:
                # print( "el -> ", el )
                el = re.split('^\d', el)[0]
                # print( "el <- ", el )
                if el is not '' and el not in alphabet:
                    alphabet.append(el)
                # print
    # print( "alphabet \t\t", alphabet )
    dictAlphabetIn = {}
    index = 0
    for symbol in alphabet:
        dictAlphabetIn[symbol] = index
        index += 1
    # print( "alphabet:index \t\t", dictAlphabetIn )
    return dictAlphabetIn


def rulesToMatrices(rules, dictAlphabetIn):
    nRules = len(rules)
    nAlphabet = len(dictAlphabetIn)
    lftSide = np.zeros((nRules, nAlphabet), int)
    rgtSide = np.zeros((nRules, nAlphabet), int)
    iRule = 0
    for line in rules:
        leftRight = line.split('->')
        # print( 'Left  >>>>\n', )
        terms = leftRight[0].split('+')
        termsToRule(terms, lftSide, iRule, dictAlphabetIn)
        # print( 'Right >>>>\n', )
        terms = leftRight[1].split('+')
        termsToRule(terms, rgtSide, iRule, dictAlphabetIn)
        # print( '------------' )
        iRule += 1
    dlt = rgtSide - lftSide
    return lftSide, rgtSide, dlt


def speciesInit(fileModelName, rulesLines, dictAlphabetIn):
    with open(fileModelName, 'r') as fileModel:
        linesToParse = fileModel.readlines()
    fileModel.close
    nAlphabetIn = len(dictAlphabetIn)
    X = np.zeros(nAlphabetIn, float)
    i = 0
    for line in linesToParse[rulesLines[0] + 1:rulesLines[1]]:
        tokens = line.strip().split('=')
        spec = tokens[0].replace(' ', '')
        amount = float(tokens[1].replace(' ', ''))
        X[dictAlphabetIn[spec]] = amount
        # print( dictAlphabetIn[spec], spec, X[dictAlphabetIn[spec]] )
        i += 1
    return X


def speciesInFeed(fileModelName, feedLines, dictAlphabetIn, X):
    with open(fileModelName, 'r') as fileModel:
        linesToParse = fileModel.readlines()
    fileModel.close
    dictFeed = {}
    for line in linesToParse[feedLines[0] + 1:feedLines[1]]:
        spec = line.strip()
        if spec in dictAlphabetIn:
            dictFeed[dictAlphabetIn[spec]] = X[dictAlphabetIn[spec]]
    return dictFeed


def modellingApproach(fileModelName):
    with open(fileModelName, 'r') as fileModel:
        linesToParse = fileModel.readlines()
    fileModel.close
    i = 0
    iLineModelling = 0
    for line in linesToParse:
        if 'Kind>' in line:
            iLineModelling = i
        i += 1
    tokens = linesToParse[iLineModelling].strip().split('Kind>')
    modelling = tokens[1].replace(' ', '')
    return modelling


def stoch2Det(Vol, xIN, leftSideIn, aStochConsts):
    NAV = Vol * 6.02214129e23
    xOut = xIN / NAV
    # print( "xIN [molecules] ", xIN )
    Kmolec = np.sum(leftSideIn, axis=1)
    # print( Kmolec )
    prods = np.prod(sps.factorial(leftSideIn, exact=False), axis=1)
    # print( prods )
    kRates = aStochConsts * (NAV**(Kmolec - 1)) / prods
    return xOut, kRates


#############
# usage example
# fileModelName = "model"
# iRulesLines, iSpecsLines, iFeedLines = findRulesAndSpecs(fileModelName)
# aRules, aConstants = rulesAndConstants(fileModelName, iRulesLines)
# dictAlphabet = rulesToAlphabet(aRules)
# print( "alphabet:index \t\t", dictAlphabet )
# leftSide, rightSide, delta = rulesToMatrices(aRules, dictAlphabet)
# print( leftSide )
# print( rightSide )
# print( delta )

# X0 = speciesInit(fileModelName, iSpecsLines, dictAlphabet)
# print( X0 )
# iRule = 0
# iSpec = 1
##############

# produttoria = np.prod(X**leftSide[iRule, :])
# print( produttoria )
# print( delta[iRule, iSpec] * produttoria )
