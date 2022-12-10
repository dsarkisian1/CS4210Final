import csv
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt

NUMBEROFROWS = 50
CROSSOVER = 0.5

def main():

    accuracy = 0

    table = openCSV("output_breast_cancer_genetic.csv")
    #while(accuracy < 50):
    decimalValues, binaryTable = createDecimalBinaryCode(table)
    geneticAlgorithm(binaryTable, decimalValues)
    accuracy = max(fitness(binaryTable))
    #print(accuracy)
    #print(fitness(binaryTable))
    
    
    
def openCSV(csvFile):
    array = []
    with open(csvFile, 'r') as breastCancer:
        reader = csv.reader(breastCancer)
        for element, row in enumerate(reader):
            if element > 0 and element < NUMBEROFROWS:
                array.append([row[0],row[1],row[2],row[3],row[4]]) 
    return array   
                

def createDecimalBinaryCode(table):
    
    decimalValues = {}
    tempDec = []
    tempBinary = []
    decimalTable = []
    binaryTable = []
    count = 0
    
    for row in table:
        for item in row:
            if item in decimalValues:
                tempDec.append(decimalValues[item])
                tempBinary.append((bin(decimalValues[item])[2:]).zfill(8))
            elif row.index(item) != 4:
                decimalValues[item] = count
                tempDec.append(decimalValues[item])
                tempBinary.append((bin(decimalValues[item])[2:]).zfill(8))
                count += 1
            else:
                if item == "Very Active":
                    #tempDec.append(0)
                    tempBinary.append("00")
                elif item == "Mod. Active":
                    #tempDec.append(1)
                    tempBinary.append("01")
                elif item == "Inactive - Exp":
                    #tempDec.append(2)
                    tempBinary.append("10")
                elif item == "Inactive - Virtual":
                    #tempDec.append(3)
                    tempBinary.append("11")
                tempDec.append(item)
        
        decimalTable.append(tempDec)
        binaryTable.append(tempBinary)
        #print(tempBinary)
        tempDec = []
        tempBinary = []
        
    print(binaryTable)
    #print(decimalValues)


    return decimalValues, binaryTable
    #print(decimalTable)
        
def geneticAlgorithm(binaryTable, decimalValues):
    fitnessTable = fitness(binaryTable)
    
    mutate(binaryTable, decimalValues)
    
    for i in range(0,int(NUMBEROFROWS*CROSSOVER)):
        gene1Index = findNextHighestFitness(fitnessTable)
        gene2Index = findNextHighestFitness(fitnessTable)
        crossover(binaryTable[gene1Index], binaryTable[gene2Index])
        #print(gene1Index)
        #print(gene2Index)
    
    
def mutate(binaryTable, decimalValues):
    numberOfMutations = random.randrange(0, NUMBEROFROWS/10)
    
    
    for i in range (0, numberOfMutations):
        randomReplacementIndex = random.randrange(0, len(decimalValues))
        randomReplacement = (bin(randomReplacementIndex)[2:]).zfill(8)
        randomRow = random.randrange(0, NUMBEROFROWS)
        randomColumn = random.randrange(0, 4)
        binaryTable[randomRow][randomColumn] = randomReplacement
    
    
            
    

    

def crossover(gene1, gene2):

    randomColumn = random.randrange(0, 4)
        
    temp = gene1[randomColumn]
    gene1[randomColumn] = gene2[randomColumn]
    gene2[randomColumn] = temp
  
    
def findNextHighestFitness(fitnessTable):
    index = fitnessTable.index(max(fitnessTable))
    fitnessTable[index] = 0
    return index
    
def fitness(binaryTable):
    fitnessTable = []
    accuracyTemp = 0
    for row in binaryTable:
        #print(row)
        accuracyTemp = 0
        for row2 in binaryTable:
            if row[0:3] == row2[0:3]:
                if row[4] == row2[4]:
                    accuracyTemp += 1
            else: 
                if row[4] != row2[4]:
                    accuracyTemp += 1
        fitnessTable.append(accuracyTemp)
        
    return fitnessTable
    

                
    
    
main()