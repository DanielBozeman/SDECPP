import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import csv


def standardPlotter():
    print("starting")
    
    df = pd.read_csv("build/output.csv", index_col=False)

    df.drop(df.tail(1).index,inplace=True) 

    print(df.shape)

    print(df.head)

    times = np.arange(df.shape[0])

    for column in df.columns:
        values = df[column]

        plt.plot(values)
        #plt.hist(values)

    plt.show()

    return

def multiGrayLastFull():
    df = pd.read_csv("build/output.csv", index_col=False)

    df.drop(df.tail(1).index,inplace=True) 

    print(df.shape)

    print(df.head)

    times = np.arange(df.shape[0])

    

    for column in df.iloc[:,:-1]:
         values = df[column]

         plt.plot(values, color='gray', alpha=0.5)

    value = df[df.columns[-2]]

    plt.plot(value, color = 'blue', zorder=2)

    plt.show()

    return

def multiGrayAveragePlot():
    df = pd.read_csv("build/output.csv", index_col=False)

    df.drop(df.tail(1).index,inplace=True) 


    averagePath = np.zeros(df.shape[0])

    for index,row in enumerate(df.to_numpy()):
        row = row[:-1]

        average = np.median(row, axis=0)

        averagePath[index] = average

    plt.plot(averagePath, color = 'red', zorder=3)

    for column in df.iloc[:,:-1]:
         values = df[column]

         plt.plot(values, color='gray', alpha=0.5)

    value = df[df.columns[-2]]

    plt.plot(value, color = 'blue', zorder=2)

    plt.show()

    return

def twoVectorHist():
    df = pd.read_csv("build/output.csv", index_col=False)

    times = df.columns[0]

    for i in range(len(df.columns)-1):
    #for i in range(5):
        values = df.iloc[:,i]

        plt.hist(values, (40*(i+1)), alpha=0.5)

    plt.show()

    return

def twoVectorAgainst():
    df = pd.read_csv("build/output.csv", index_col=False)

    print(df.shape)

    print(df.head)

    times = df.iloc[0]

    print(times)
    
    plt.plot(times)

    plt.show()

    return

def normalPlotter():
    x = np.linspace(-200, 200, 100)

    df = pd.read_csv("build/output.csv", index_col=False, header = None)

    for i in range(len(df.columns)):

        values = df.iloc[:,i]
        
        print(values)

        mean = values[0]
        variance = values[1]

        # Calculate y values using the probability density function (PDF)
        y = norm.pdf(x, mean, variance)  # Mean = 0, Standard Deviation = 1

        # Plot the curve
        plt.plot(x, y)

    plt.title("Normal Distribution Curves")
    plt.xlabel("x")
    plt.ylabel("Probability Density")
    plt.show()

def trueVsEst():

    #Reading data
    dfData = pd.read_csv("build/dataOutput.csv", index_col=False)

    dfEstimation = pd.read_csv("build/estimationOutput.csv", index_col=False)
    
    dfTrue = pd.read_csv("build/trueOutput.csv", index_col=False)

    print(dfData.shape)
    print(dfEstimation.shape)
    print(dfTrue.shape)

    estimateAveragePath = np.zeros(dfEstimation.shape[0])

    for index,row in enumerate(dfEstimation.to_numpy()):
        row = row[:-1]

        average = np.median(row, axis=0)

        estimateAveragePath[index] = average

    trueAveragePath = np.zeros(dfTrue.shape[0])

    for index,row in enumerate(dfTrue.to_numpy()):
        row = row[:-1]

        average = np.median(row, axis=0)

        trueAveragePath[index] = average
    
    for column in dfEstimation.iloc[:,:-1]:
        values = dfEstimation[column]

        plt.plot(values, color='green', alpha=0.25, linewidth = 0.5)

    plt.plot(estimateAveragePath, color = 'yellow')

    # for column in dfTrue.iloc[:,:-1]:
    #       values = dfTrue[column]

    #       plt.plot(values, color='green', alpha=0.25, linewidth = 0.5)

    # plt.plot(trueAveragePath, color = 'yellow')

    value = dfData[dfData.columns[-2]]

    #print(dfData.head)

    plt.plot(value, color = 'red', zorder=2, linewidth = 1)

    #plt.xlim(-5,1200)

    #plt.ylim(-5,15)

    #plt.yticks(np.arange(-5,16,2))

    plt.show()

    return

if __name__ == "__main__":
    #multiGrayLastFull()
    standardPlotter()
    #multiGrayAveragePlot()
    #twoVectorAgainst()
    #twoVectorHist()
    #normalPlotter()
    #trueVsEst()

    
    