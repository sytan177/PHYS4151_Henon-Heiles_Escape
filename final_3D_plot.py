import numpy as np
from mpl_toolkits import mplot3d
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import csv
import sys
from sympy import var, plot_implicit
from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list('basin', ['black', 'green', 'orange', 'blue'], N=4)

r = np.sqrt(2)

# Define Henon-Heiles potential
def V_hh(x,y):
    v = 1 / 2 * (x ** 2 + y ** 2) + (x ** 2 * y - y ** 3 / 3)
    return v

def H(x, y, xdot, ydot):
    h = V_hh(x,y) - 1/2*(xdot**2+ydot**2)
    return h

h = [0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3]

def get_plot_E(E):
    initial_path = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/3D/hh_initial_300_%s.txt'%(E)
    exit_file = "/Users/shuyut/Downloads/PHYS4151_Project/3D/escape_list_3D_%s_10000s.txt"%(E)

    ini_states = []
    with open(initial_path) as inifile:
        for line in inifile:
            ini_states.append(float(line))
    # for x-y plane
    ini_states = np.array(ini_states).reshape(int(len(ini_states)/4), 4)
    x0 = ini_states[:,0]
    y0 = ini_states[:,1]
    dy0 = ini_states[:,3]

    file = open(exit_file,'r')
    lines = file.readlines()
    tokens_column_number = 1
    result_table=[]
    for x in lines:
        result_table.append(x.split())
    file.close()

    result_table = pd.DataFrame(result_table, columns=['ind', 'x0', 'y0', 'dy0', 'exit', 't'])
    result_table = result_table.astype(float)

    ind1 = result_table.index[result_table['exit'] == 1.0 ].tolist()
    escape_1 = result_table['ind'].astype(int)[ind1]-1
    ind2 = result_table.index[result_table['exit'] == 2.0 ].tolist()
    escape_2 = result_table['ind'].astype(int)[ind2]-1
    ind3 = result_table.index[result_table['exit'] == 3.0 ].tolist()
    escape_3 = result_table['ind'].astype(int)[ind3]-1
    tol_escape = np.concatenate((np.array(escape_1), np.array(escape_2), np.array(escape_3))).tolist()

    setT = set(range(len(x0)))
    non_escape = list(setT.difference(tol_escape))

    return result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape


def get_plot_slice(dyr):
    initial_path = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/3D/hh_initial_300_%s.txt'%(0.3)
    exit_file = "/Users/shuyut/Downloads/PHYS4151_Project/3D/escape_list_3D_%s_10000s.txt"%(0.3)

    ini_states = []
    with open(initial_path) as inifile:
        for line in inifile:
            ini_states.append(float(line))
    # for x-y plane
    ini_states = np.array(ini_states).reshape(int(len(ini_states)/4), 4)
    x0 = ini_states[:,0]
    y0 = ini_states[:,1]
    dy0 = ini_states[:,3]

    file = open(exit_file,'r')
    lines = file.readlines()
    tokens_column_number = 1
    result_table=[]
    for x in lines:
        result_table.append(x.split())
    file.close()

    result_table = pd.DataFrame(result_table, columns=['ind', 'x0', 'y0', 'dy0', 'exit', 't'])
    result_table = result_table.astype(float)
    result_table = result_table[(round(result_table['dy0'], 2) == dyr)]
    ind1 = result_table.index[result_table['exit'] == 1.0 ].tolist()
    escape_1 = result_table['ind'].astype(int)[ind1]-1
    ind2 = result_table.index[result_table['exit'] == 2.0 ].tolist()
    escape_2 = result_table['ind'].astype(int)[ind2]-1
    ind3 = result_table.index[result_table['exit'] == 3.0 ].tolist()
    escape_3 = result_table['ind'].astype(int)[ind3]-1
    tol_escape = np.concatenate((np.array(escape_1), np.array(escape_2), np.array(escape_3))).tolist()

    setT = set(range(len(x0)))
    non_escape = list(setT.difference(tol_escape))

    return result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape


# exit channels plot in x-y plane
def threeD_channel():
    fig = plt.figure(figsize=(15, 13.8))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1, projection='3d')
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        print(len(escape_1), len(escape_2), len(escape_3))
        color = np.zeros(len(x0))
        color[list(escape_1)] = 1.
        color[list(escape_2)] = 2.
        color[list(escape_3)] = 3.
        ax.scatter(x0, y0, dy0, s=0.05, c=color, cmap = cmap, vmin=0., vmax=3.)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\dot{y}$');
        #ax.set_xlim(-1.5, 1.5)
        #ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.4, 1.4)
        ax.set_title("h=%s" % (h[i]), y=1.07, pad=-14)
    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return

dy = [-0.59, -0.45, -0.31, -0.11, 0, 0.11, 0.31, 0.45, 0.59]
def threeD_slice():
    fig = plt.figure(figsize=(15, 13.8))
    for i in range(len(dy)):
        ax = fig.add_subplot(3,3,i+1, projection='3d')
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_slice(dy[i])
        a = dy0
        color = np.zeros(len(x0))
        color[list(escape_1)] = 1.
        color[list(escape_2)] = 2.
        color[list(escape_3)] = 3.
        ax.scatter(x0, y0, dy0, s=0.05, c=color, cmap='rainbow', vmin=0., vmax=3.)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\dot{y}$');
        #ax.set_xlim(-1.5, 1.5)
        #ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.4, 1.4)
        title = r"$\dot{y}"+"=$%s" % (dy[i])
        ax.set_title(title, y=1.05, pad=-14)

    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return a

def threeD_t():
    fig = plt.figure(figsize=(17, 13.8))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1, projection='3d')
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        t_escape = result_table['t']
        x_e = result_table['x0']
        y_e = result_table['y0']
        dy_e = result_table['dy0']
        points = ax.scatter(x_e, y_e, dy_e, s=0.03, c=np.log10(t_escape), cmap="nipy_spectral", vmin=0., vmax=4.)
        cb = plt.colorbar(points)
        cb.set_label(r'log($t_esp$)')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\dot{y}$');
        #ax.set_xlim(-1.5, 1.5)
        #ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.4, 1.4)
        ax.set_title("h=%s" % (h[i]), y=1.07, pad=-14)

    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return

def threeD_slice_t():
    fig = plt.figure(figsize=(18, 13.8))
    for i in range(len(dy)):
        ax = fig.add_subplot(3,3,i+1, projection='3d')
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_slice(dy[i])
        a = dy0
        t_escape = result_table['t']
        x_e = result_table['x0']
        y_e = result_table['y0']
        dy_e = result_table['dy0']
        points = ax.scatter(x_e, y_e, dy_e, s=0.03, c=np.log10(t_escape), cmap="rainbow", vmin=-2, vmax=4.)
        cb = plt.colorbar(points)
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\dot{y}$');
        #ax.set_xlim(-1.5, 1.5)
        #ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.4, 1.4)
        cb.set_label(r'log($t_esp$)')
        title = r"$\dot{y}"+"=$%s" % (dy[i])
        ax.set_title(title, y=1.05, pad=-14)

    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return a

threeD_channel()