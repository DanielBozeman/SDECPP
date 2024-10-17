import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv


def standardPlotter():
    df = pd.read_csv("build/output.csv", index_col=False)

    df.drop(df.tail(1).index,inplace=True) 

    print(df.shape)

    print(df.head)

    times = np.arange(df.shape[0])

    for column in df.columns:
        values = df[column]

        plt.plot(values)

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

if __name__ == "__main__":
    multiGrayLastFull()
    #standardPlotter()
    #multiGrayAveragePlot()

    
    