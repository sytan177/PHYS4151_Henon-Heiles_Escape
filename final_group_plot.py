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

def ZVC(E):
    esp = 1e-5
    r = np.sqrt(2)
    delta = 0.001
    x = np.arange(-r, r, delta)
    y = np.arange(-r, r, delta)
    coord = []
    for x_i in x:
        for y_i in y:
            if x_i**2+y_i**2 <= r**2 and abs(V_hh(x_i,y_i) - E) < esp:
                coord.append([x_i,y_i])
    return coord

def draw_circle(r):
    theta = np.linspace(0, 2 * np.pi, 100)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    plt.plot(x,y,ls='dotted')
    return

h = [0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30]

def get_plot_E(E):
    #initial_path = "/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_300_%s.txt"%(E)
    initial_path = "/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_500_phase_%s.txt"%(E)
    #exit_file = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_%s_10000s.txt"%(E)
    exit_file = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_500_phase_%s_10000s.txt"%(E)

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


# exit channels plot in x-y plane
def xy_channel():
    fig = plt.figure(figsize=(10.2, 10.0))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1)
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        color = np.zeros(len(x0))
        color[list(escape_1)] = 1.
        color[list(escape_2)] = 2.
        color[list(escape_3)] = 3.
        ax.scatter(x0, y0, s=0.01, c=color, cmap=cmap, vmin=0., vmax=3.)
        draw_circle(np.sqrt(2.))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_title("h=%s" % (h[i]))
        # plot ZVC
        x = np.linspace(min(x0), max(x0), 100)
        y = np.linspace(min(y0), max(y0)-0.0085*(-min(y0)+max(y0)), 100)

        X, Y = np.meshgrid(x, y)
        pot = V_hh(X, Y)
        cp = ax.contour(x, y, pot, levels=[h[i]], colors='black', linewidths=1.3)

    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return

# escape period plot in x-y plane
def xt_t():
    fig = plt.figure(figsize=(12.5, 10.0))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1)
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        t_escape = result_table['t']
        x_e = result_table['x0']
        y_e = result_table['y0']
        points = plt.scatter(x_e, y_e, s = 0.5, c = np.log10(t_escape), cmap = "YlGn")
        cb = plt.colorbar(points)
        draw_circle(np.sqrt(2.))
        plt.clim(0.,4)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title("h=%s"%(h[i]))
        cb.set_label(r'log($t_{esp}$)')
        x = np.linspace(min(x0), max(x0), 100)
        y = np.linspace(min(y0), max(y0) - 0.0085 * (-min(y0) + max(y0)), 100)

        X, Y = np.meshgrid(x, y)
        pot = V_hh(X, Y)
        cp = ax.contour(x, y, pot, levels=[h[i]], colors='black', linewidths=1.3)
    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_t.png')
    plt.show()
    return

def ydy_channel():
    fig = plt.figure(figsize=(10.2, 10.0))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1)
        result_table, x0, y0, dy0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        color = np.zeros(len(x0))
        color[list(escape_1)] = 1.
        color[list(escape_2)] = 2.
        color[list(escape_3)] = 3.
        ax.scatter(y0, dy0,  s=0.005, c=color, cmap=cmap, vmin=0., vmax=3.)
        ax.set_xlabel('y')
        ax.set_ylabel(r'$\dot{y}$')
        #ax.set_xlim(-1.5, 1.5)
        #ax.set_ylim(-1.5, 1.5)
        ax.set_title("h=%s" % (h[i]))

    fig.set_facecolor('w')
    plt.tight_layout()
    plt.savefig('xy_channel.png')
    plt.show()
    return

# escape period plot in x-y plane
def yt_t():
    fig = plt.figure(figsize=(12.5, 10.0))
    for i in range(len(h)):
        ax = fig.add_subplot(3,3,i+1)
        result_table, x0, y0, escape_1, escape_2, escape_3, non_escape = get_plot_E(h[i])
        t_escape = result_table['t']
        x_e = result_table['y0']
        y_e = result_table['dy0']
        points = plt.scatter(x_e, y_e, s = 0.5, c = np.log10(t_escape), cmap = "YlGn")
        cb = plt.colorbar(points)
        ax.set_title("h=%s"%(h[i]))
        cb.set_label(r'log($t_{esp}$)')
        plt.clim(0.,4)
    fig.set_facecolor('w')
    plt.tight_layout()
    #plt.savefig('xy_t.png')
    plt.show()
    return

ydy_channel()
'''
# for x-y plane plot
plt.figure(figsize=(6.5,6.5))
plt.scatter(x0[escape_1], y0[escape_1], s = 0.1, c = 'green')
plt.scatter(x0[escape_3], y0[escape_3], s = 0.1, c = 'blue')
plt.scatter(x0[escape_2], y0[escape_2], s = 0.1, c = 'orange')
plt.scatter(x0[non_escape], y0[non_escape], s = 0.1, c = "black")
plt.xlabel('y')
plt.ylabel(r'$dy$')
plt.title("h=%s"%(E0))
plt.show()
'''

