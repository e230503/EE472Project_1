How To run the script:

-> Project has only one script named "main.py".
-> In main.py there is code section at the lines 989 and 990:
if __name__ == "__main__":
    Run("ieee57cdf.txt", 100, 5e-13)
Run basically executes all the code. There are 3 parameters for that function:
first one -> directory of the data that is going to be read.
second one -> maximum amounth of iteration for the Newton Raphson method if it is exceeded simulation ends
without checking anything.
third one -> max value of deltaP and deltaQ can take for clarify the convergence and end the simulation.
-> One can change parameters in "Run" method and run the code.
-> PyCharm is used while constructing the script "main.py".

In this code:
import time
import matplotlib.pyplot as plt
import cmath as mat
import numpy as np
import math as mat
from enum import Enum

are imported one can download required libraries

