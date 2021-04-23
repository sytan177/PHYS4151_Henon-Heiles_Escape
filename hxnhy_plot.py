import numpy as np
from mpl_toolkits import mplot3d
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list('basin', ['black', 'green', 'orange', 'blue'], N=4)

def HHH(x,y,dx,dy):
    E =1./2.*(dx**2.+dy**2.+x**2.+y**2.)+x**2.*y-1./3.*y**3.
    return E.round(decimals=9)

def list_TF(escape_list, time_list):
    tb = [1]*len(epsilon)
    for i in range(len(epsilon)):
        ind = 2*i
        if str(escape_list[ind+1]) == str(escape_list[0]) and str(escape_list[ind+2]) == str(escape_list[0]):
            if abs(time_list[ind+1] - time_list[0]) < 1.0e-2 and abs(time_list[ind+2] - time_list[0]) < 1.0e-2:
                tb[i] = 0
    # return uncertain states
    return tb

def power_law(x, a, b, c):
    return a*x**b+c

def linear_fit(x, a, b):
    return a*x+b

n = 500
r_sca_reg = np.sqrt(2.)
'''
x0 = np.linspace(-r_sca_reg, r_sca_reg, n )
y0 = np.linspace(-r_sca_reg, r_sca_reg, n )
E0 = np.linspace(1./6., 0.4, n )
'''

xh_init = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_500_coord_xh.txt'
yh_init = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_500_coord_yh.txt'

xh_file = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_500_xh_10000s.txt"
yh_file = "/Users/shuyut/Downloads/PHYS4151_Project/escape_list_500_yh_10000s.txt"

fractal_init_xh = '/Users/shuyut/Downloads/PHYS4151_Project/HH_Fractal/fractal_init_coord/hh_initial_fractal_xh.txt'
fractal_final_xh = '/Users/shuyut/Downloads/PHYS4151_Project/HH_Fractal/escape_list_fractal_xh_10000s.txt'

fractal_init_yh = '/Users/shuyut/Downloads/PHYS4151_Project/HH_Fractal/fractal_init_coord/hh_initial_fractal_yh.txt'
fractal_final_yh = '/Users/shuyut/Downloads/PHYS4151_Project/HH_Fractal/escape_list_fractal_yh_10000s.txt'

initial_path = xh_init
exit_file = xh_file

ini_states = []
with open(initial_path) as inifile:
    for line in inifile:
        ini_states.append(float(line))

# for x-y plane
ini_states = np.array(ini_states).reshape(int(len(ini_states)/4), 4)
x0 = ini_states[:,0]
y0 = ini_states[:,1]
dx0 = ini_states[:,2]
dy0 = ini_states[:,3]
h0 = HHH(x0, y0, dx0, dy0)

file = open(exit_file,'r')
lines = file.readlines()
tokens_column_number = 1
result_table=[]
for x in lines:
    x.rstrip('\x00')
    result_table.append(x.split()[-6:])
file.close()

result_table = pd.DataFrame(result_table, columns=['ind', 'x0', 'y0', 'dy0', 'exit', 't'])
result_table = result_table.astype(float)

initial_coords = pd.DataFrame(zip(h0, x0, y0, dx0, dy0), columns=["h", "x", "y", "v_x", "v_y"])


ind1 = result_table.index[result_table['exit'] == 1.0 ].tolist()
escape_1 = result_table['ind'].astype(int)[ind1]-1
ind2 = result_table.index[result_table['exit'] == 2.0 ].tolist()
escape_2 = result_table['ind'].astype(int)[ind2]-1
ind3 = result_table.index[result_table['exit'] == 3.0 ].tolist()
escape_3 = result_table['ind'].astype(int)[ind3]-1
tol_escape = np.concatenate((np.array(escape_1), np.array(escape_2), np.array(escape_3))).tolist()

#setE = set(tol_escape)
setT = set(range(len(x0)))
non_escape = list(setT.difference(tol_escape))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={'width_ratios': [1, 1.2], 'wspace' : 0.3})

color = np.zeros(len(x0))
color[list(escape_1)] = 1.
color[list(escape_2)] = 2.
color[list(escape_3)] = 3.
ax1.scatter(x0, h0, s=0.01, c=color, cmap = cmap, vmin=0., vmax=3.)
ax1.set_xlabel('X')
ax1.set_ylabel('h')

t_escape = result_table['t']
x_e = result_table['y0']
ind_e = result_table['ind'].astype(int)-1
y_e = h0[ind_e]
points = ax2.scatter(x_e, y_e, s = 0.5, c = np.log10(t_escape), cmap = "YlGn")
ax2.set_xlabel('y')
ax2.set_ylabel('h')
cb = plt.colorbar(points)
cb.set_label(r'log($t_{esp}$)')
plt.show()

'''
final_states = np.zeros(len(h0))
final_time = np.zeros(len(h0))

final_states[escape_1] = "1"
final_states[escape_2] = "2"
final_states[escape_3] = "3"
final_states[non_escape] = "0"

final_time[escape_1] = result_table['t'].astype(float)[ind1]
final_time[escape_2] = result_table['t'].astype(float)[ind2]
final_time[escape_3] = result_table['t'].astype(float)[ind3]
final_time[non_escape] = 0.

h_values = list(dict.fromkeys(h0))
epsilon = [1.0e-1, 5.0e-2, 1.0e-2, 5.0e-3, 1.0e-3, 5.0e-4, 1.0e-4, 5.0e-5, 1.0e-5, 5.0e-6, 1.0e-6]
len_o = 2*len(epsilon)+1
#D = list(range(len(h_values)))
D = [0]*len(h_values)
for i in range(len(h_values)):
#for i in range(1):
    inds = [j for j, x in enumerate(h0) if x == h_values[i]]
    num_saved = int(len(inds)/len_o)
    r_tb = [0]*len(epsilon)
    k = 1
    while k < num_saved:
        ind = inds[len_o*(k-1):len_o*k]
        escape_con1 = final_states[ind]
        escape_con2 = final_time[ind]
        r_tb = [x + y for x, y in zip(r_tb, list_TF(escape_con1, escape_con2))]
        #print(list_TF(escape_con1, escape_con2))
        k += 1
    ucty_D = np.array(r_tb)/float(num_saved)
    #ucty_D = [ucty_D[i]/epsilon[i] for i in range(len(epsilon))]
    #popt, pcov = curve_fit(power_law, epsilon, ucty_D)
    try:
        x_to_fit = np.array(np.log10(epsilon))
        y_to_fit = np.array(np.log10(ucty_D))
        m, b = np.polyfit(x_to_fit, y_to_fit, 1)
        D[i] = 1 - m
    except:
        D[i] = 1 - 0
    #plt.scatter(x_to_fit, y_to_fit)
    #plt.plot(x_to_fit, m*x_to_fit+b)

plt.figure(figsize=(5.2,4.8))
plt.plot(h_values, D, color='green')
plt.xlabel("E")
plt.ylabel(r"$D_{0}$")
plt.ylim([0.6,1.])
plt.show()
'''