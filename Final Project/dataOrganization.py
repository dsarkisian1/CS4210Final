#####
#This program will read amino acid codes, sort them into 2,3,4,5... sequences, 
#and then calculate the bell curve for each sequence
#Once the sequence with the best distribution is found, the ID3 machine learning will 
#be used to make a predictive model of which sequences are most likely to be active
#####

import csv
import numpy as np
import matplotlib.pyplot as plt

MAX_SEQUENCE = 7
MIN_SEQUENCE = 2
#This will be the header of the csv file
feilds = ['Sequence', 'Frequency Very Active', 'Frequency Mod. Active', 'Frequency Inactive - Exp', 'Frequency Inactive - Virtual', 'Frequency Total', 'Sequence Length']       
    

def main():
    
    sequenceValuesBreastCancer = extractAnalysis('ACPs_Breast_cancer.csv', 'output_breast_cancer.csv')
    createBellCurve(sequenceValuesBreastCancer, 'Breast Cancer')
    sequenceValuesLungCancer = extractAnalysis('ACPs_Lung_cancer.csv', 'output_lung_cancer.csv')
    createBellCurve(sequenceValuesLungCancer, 'Lung Cancer')
    outputPreProcessed(sequenceValuesBreastCancer, 'ACPs_Breast_cancer_preprocessed.csv', feilds, 4)
    outputPreProcessed(sequenceValuesLungCancer, 'ACPs_Lung_cancer_preprocessed.csv', feilds, 4)
   
#####extractAnalysis#####
#This function is automating the process of extracting the data from the csv file
#then running the data through the preprocessing to prepare for the machine learning algorithms
#preprocessing includes dividing sequences into small segments, tracking their frequency
#then using statistics to perform bell curves to determine the best sequence length (based off distribution)
#########################
def extractAnalysis(csvFileRead, csvFileWrite): 
    array = openCSV(csvFileRead) #opens the csv file and reads the data, writes to array
    
    sequenceValues = createsequenceValues(array) #dictionary that holds each list of sequences

   #Creation of output.csv
    outputCSV(csvFileWrite, feilds, sequenceValues)  
    return sequenceValues

#######openCSV######
#opens the csv file and reads the data, writes to array
####################
def openCSV(csvFile):
    array = []
    with open(csvFile, 'r') as breastCancer:
        reader = csv.reader(breastCancer)
        for element, row in enumerate(reader):
            if element > 0:
                array.append([row[1],row[2]]) 
    return array      
            
######createsequenceValues######
#This function will take the array of sequences and divide them into smaller sequences
#It will also track the frequency of their occurances
################################
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
                 
######outputCSV######
#this function will take the headers for each column, and the dictionary of sequences,
#and will write the data to a csv file
#####################   
def outputCSV(csvFile, feilds, sequenceValues):
    #this is adding each sequence to the csv file rows
    rows = []
    for key in sequenceValues:
        rows.append([key, sequenceValues[key][0], sequenceValues[key][1], sequenceValues[key][2], sequenceValues[key][3], sequenceValues[key][4], sequenceValues[key][5]])
     
    with open(csvFile, "w") as csvfile:
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerow(feilds)
        csvwriter.writerows(rows)
 
    
#####BELL CURVE CALCULATION#####
def bellCurve(sequence, mu, sigma):
    y = []
    for i in range (0, len(sequence)):
        y.append(1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-((sequence[i] - mu)**2)/(2 * sigma**2)))
        
    return y

def createBellCurve(sequenceValues, cancerType):
    statistics = {} #a dictionary that holds the statistics for each sequence length
    sequence = [] #an array that holds each sequence frequency
    count = 0
    sum = 0
    sequenceLength = MIN_SEQUENCE
    
    #calculating the number of occurances for each length of sequence to find the mean
    for key in sequenceValues:    
        tempSequenceLength = sequenceValues[key][5]
        
        #this increments a counter that resets the index of the 2d array column to 0 when the 
        #row is incremented
        if sequenceLength == tempSequenceLength:
            count += 1
            sum += sequenceValues[key][4] #calculating the total number of occurances for each sequence length
            sequence.append(sequenceValues[key][4]) #the frequency of the sequence
        else:
            sequenceLength = tempSequenceLength
            statistics[tempSequenceLength] = [sum,count] #total number of all occurances
            sum = 0
            count = 0
        
            #calculating the mean for each length of sequence
            for key in statistics:
                mean = statistics[key][0]/statistics[key][1]
        
                #calculating the value of each occurance count - the mean squared
                stDevNumerator = 0
                for value in sequence:
                    stDevNumerator += (value - mean)**2
            
                #calculating the standard deviation for each length of sequence
                stDev = (stDevNumerator/statistics[key][1])**.5

                sequence.sort()

                yValues = bellCurve(sequence, mean, stDev)
                plt.style.use('seaborn')
                plt.figure(figsize = (6, 6))
                plt.plot(sequence, yValues, color = 'black', linestyle = 'dashed')
                plt.scatter(sequence, yValues, marker = 'o', s = 25, color = 'red')
                plt.title(cancerType + " Sequence Length: " + str(sequenceLength - 1))
                plt.show()
                #plt.savefig(cancerType + " Sequence Length: " + str(sequenceLength - 1) + ".png")

                sequence.clear()
            
def outputPreProcessed(sequenceValues, csvFile, feilds, sequenceLength):
    rows = []
    for sequence in sequenceValues:
        if sequenceValues[sequence][5] == sequenceLength:
            rows.append([sequence, sequenceValues[sequence][0], sequenceValues[sequence][1], sequenceValues[sequence][2], sequenceValues[sequence][3], sequenceValues[sequence][4], sequenceValues[sequence][5]])

    with open(csvFile, "w") as csvfile:
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerow(feilds)
        csvwriter.writerows(rows)     


    
#####MAIN#####
#Calling of main 
main()