import csv
import numpy as np
import matplotlib.pyplot as plt
import random

MAX_SEQUENCE = 7
MIN_SEQUENCE = 2
MAX_COMMON = 10
SEQUENCE_LENGTH = 8


def openCSV(csvFile):
    array = []
    with open(csvFile, 'r') as breastCancer:
        reader = csv.reader(breastCancer)
        for element, row in enumerate(reader):
            if element > 0:
                array.append([row[1],row[2]]) 
    return array     


def main():
    
    table = extractAnalysis('ACPs_Breast_cancer.csv', 'output_breast_cancer.csv')
    outputCSV('output_breast_cancer_genetic.csv', table)
    
    table = extractAnalysis('ACPs_Lung_cancer.csv', 'output_lung_cancer.csv')
    outputCSV('output_lung_cancer_genetic.csv', table)
    
    #createBellCurve(sequenceValuesBreastCancer, 'Breast Cancer')
    #sequenceValuesLungCancer = extractAnalysis('ACPs_Lung_cancer.csv', 'output_lung_cancer.csv')
    
def extractAnalysis(csvFileRead, csvFileWrite):
    db = openCSV(csvFileRead)
    sequenceValues = createsequenceValues(db)
    return findPopularSequences(sequenceValues)
    
    #print(sequenceValues)
    
def createsequenceValues(array):
    sequenceValues = {}

    for i in range (MIN_SEQUENCE,MAX_SEQUENCE): #for each sequence length
        for j in range (0,len(array)):    #for each amino acid
            for k in range (0, len(array[j][0]) - i + 1): #for each sequence in the amino acid 
                #print(len(array[j][0]))
                
                sequence = array[j][0][0+k:i+k] #sequence is value of length max-min
                if sequence in sequenceValues: #if the sequence is already in the dictionary, increment frequency
                    classValue = array[j][1]
                    if classValue == "very active":
                        sequenceValues[sequence][0] += 1
                        sequenceValues[sequence][4] += 1
                    elif classValue == "mod. active":
                        sequenceValues[sequence][1] += 1
                        sequenceValues[sequence][4] += 1
                    elif classValue == "inactive - exp":
                        sequenceValues[sequence][2] += 1
                        sequenceValues[sequence][4] += 1
                    elif classValue == "inactive - virtual":
                        sequenceValues[sequence][3] += 1
                        sequenceValues[sequence][4] += 1
                else: #if the sequence is not in the dictionary, add it
                    classValue = array[j][1]
                    if classValue == "very active":
                        sequenceValues[sequence] = [1, 0, 0, 0, 1, i]
                    elif classValue == "mod. active":
                        sequenceValues[sequence] = [0, 1, 0, 0, 1, i]
                    elif classValue == "inactive - exp":
                        sequenceValues[sequence] = [0, 0, 1, 0, 1, i]
                    elif classValue == "inactive - virtual":
                        sequenceValues[sequence] = [0, 0, 0, 1, 1, i]
    return sequenceValues

def findPopularSequences(sequenceValues):
    popularSequences = []
    feilds = ['Sequence', 'Frequency Very Active', 'Frequency Mod. Active', 'Frequency Inactive - Exp', 'Frequency Inactive - Virtual', 'Frequency Total', 'Sequence Length']       

    max10VeryActive = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    max10VeryActiveSequences = ["", "", "", "", "", "", "", "", "", ""]
    max10ModActive = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    max10ModActiveSequences = ["", "", "", "", "", "", "", "", "", ""]
    max10InactiveExp = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    max10InactiveExpSequences = ["", "", "", "", "", "", "", "", "", ""]
    max10InactiveVirtual = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    max10InactiveVirtualSequences = ["", "", "", "", "", "", "", "", "", ""]
    
    allMax = {}
    
    for key in sequenceValues:
        
        if sequenceValues[key][5] == 2:
        
            for i in range(0,MAX_COMMON):
                if sequenceValues[key][0] > max10VeryActive[i]:
                    max10VeryActive[i] = sequenceValues[key][0]
                    max10VeryActiveSequences[i] = key
                    break
            
            for i in range(0,MAX_COMMON):        
                if sequenceValues[key][1] > max10ModActive[i]:
                    max10ModActive[i] = sequenceValues[key][1]
                    max10ModActiveSequences[i] = key
                    break
                    
            for i in range(0,MAX_COMMON):
                if sequenceValues[key][2] > max10InactiveExp[i]:
                    max10InactiveExp[i] = sequenceValues[key][2]
                    max10InactiveExpSequences[i] = key
                    break
                    
            for i in range(0,MAX_COMMON):
                if sequenceValues[key][3] > max10InactiveVirtual[i]:
                    max10InactiveVirtual[i] = sequenceValues[key][3]
                    max10InactiveVirtualSequences[i] = key
                    break
                
              
    allMax['Very Active'] = [max10VeryActive, max10VeryActiveSequences]
    allMax['Mod. Active'] = [max10ModActive, max10ModActiveSequences]
    allMax['Inactive - Exp'] = [max10InactiveExp, max10InactiveExpSequences]
    allMax['Inactive - Virtual'] = [max10InactiveVirtual, max10InactiveVirtualSequences]
    
    return createTable(allMax)
        
def createTable(allMax):
    
    #for array in allMax:
    #    print(allMax[array])
    table = []
    temp = []
    countVeryActive = 0
    countModActive = 0
    countInactiveExp = 0
    countInactiveVirtual = 0
    maxVal = 0
        
    for i in range (0,MAX_COMMON*MAX_COMMON):
        
        for j in range (0,SEQUENCE_LENGTH):
            randomRow = random.randrange(0, 4)
            randomColumn = random.randrange(0, MAX_COMMON)

            if randomRow == 0:
                randomRow = 'Very Active'
                countVeryActive += allMax[randomRow][0][randomColumn]
            if randomRow == 1:
                randomRow = 'Mod. Active'
                countModActive += allMax[randomRow][0][randomColumn]
            if randomRow == 2:
                randomRow = 'Inactive - Exp'
                countInactiveExp += allMax[randomRow][0][randomColumn]
            if randomRow == 3:  
                randomRow = 'Inactive - Virtual'
                countInactiveVirtual += allMax[randomRow][0][randomColumn]
            
            temp.append(allMax[randomRow][1][randomColumn])
            
        maxVal = max(countVeryActive, countModActive, countInactiveExp, countInactiveVirtual)
        
        if countVeryActive == maxVal:
            temp.append('Very Active')
        if countModActive == maxVal:
            temp.append('Mod. Active')
        if countInactiveExp == maxVal:
            temp.append('Inactive - Exp')
        if countInactiveVirtual == maxVal:
            temp.append('Inactive - Virtual')    
        table.append(temp)
        temp = []
        countVeryActive = 0
        countModActive = 0
        countInactiveExp = 0
        countInactiveVirtual = 0
    
    return table
    
    #print(table)

def outputCSV(csvFile, table):
    #this is adding each sequence to the csv file rows
    rows = []
    temp = []
    for key in table:
        temp = []
        for item in range(0,SEQUENCE_LENGTH//2):
            temp.append(key[item])
        temp.append(key[SEQUENCE_LENGTH])
        rows.append(temp)
     
    with open(csvFile, "w") as csvfile:
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerows(rows)

                    
            
        
main()