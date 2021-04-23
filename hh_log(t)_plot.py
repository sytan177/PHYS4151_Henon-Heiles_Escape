import numpy as np
from mpl_toolkits import mplot3d
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import random
import csv
import sys

h = [0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.35, 0.40, 0.45, 0.50]
t_xy = []
t_ydy = []

for E0 in h:
    exit_file_ydy = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_500_phase_%s_10000s.txt"%(E0)
    exit_file_xy = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_%s_10000s.txt"%(E0)
    final_state_xy = []
    final_state_ydy = []

    file = open(exit_file_ydy,'r')
    lines = file.readlines()
    tokens_column_number = 1
    result_table=[]
    for x in lines:
        result_table.append(x.split())
    file.close()
    result_table = pd.DataFrame(result_table, columns=['ind', 'x0', 'y0', 'dy0', 'exit', 't'])
    result_table = result_table.astype(float)

    l1 = random.sample(range(len(result_table)), 30000)

    t_ydy.append(np.average(result_table['t'][:]))

    file1 = open(exit_file_xy, 'r')
    lines1 = file1.readlines()
    tokens_column_number = 1
    result_table1 = []
    for x1 in lines1:
        result_table1.append(x1.split())
    file1.close()

    result_table1 = pd.DataFrame(result_table1, columns=['ind', 'x0', 'y0', 'dy0', 'exit', 't'])
    result_table1 = result_table1.astype(float)

    t_xy.append(np.average(result_table1['t']))


plt.figure(figsize=(5., 5))
plt.plot(h, np.log10(t_ydy))
plt.scatter(h, np.log10(t_ydy), s=10.)
plt.plot(h, np.log10(t_xy))
plt.scatter(h, np.log10(t_xy), s=10.)
plt.ylim([0.5,3.0])
plt.xlabel('h')
plt.ylabel(r'$log_{10}(<t_{esc}>)$')

#plt.scatter(result_table['y0'][-60000:], result_table['dy0'][-60000:])
plt.show()