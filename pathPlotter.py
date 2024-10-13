import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":

    df = pd.read_csv("build/output.csv", index_col=False)

    df.drop(df.tail(1).index,inplace=True) 

    print(df.shape)

    print(df.head)

    times = np.arange(df.shape[0])

    for column in df.columns:
        values = df[column]

        plt.plot(values)

    plt.show()
    