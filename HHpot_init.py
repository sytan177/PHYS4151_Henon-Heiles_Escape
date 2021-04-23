import numpy as np
from mpl_toolkits import mplot3d
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import warnings

# Initial value of Hamiltonian

n = 100

h = [0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30]
h_extra = [0.35, 0.40, 0.45, 0.50]
#epsilon = [7.0e-3, 5.0e-3, 3.0e-3, 1.0e-3, 7.0e-4, 5.0e-4, 3.0e-4]
epsilon = [1.0e-1, 5.0e-2, 1.0e-2, 5.0e-3, 1.0e-3, 5.0e-4, 1.0e-4, 5.0e-5, 1.0e-5, 5.0e-6, 1.0e-6]
r_sca_reg = np.sqrt(2.)

def polar_V(r, phi):
    pot = 1./2.*r**2 + 1./3.*r**3*np.sin(3*phi)
    return pot

def find_init(r, phi, dr, dphi, E):
    V = polar_V(r, phi)
    to_solve = V + 1./2.*( dr**2. + r**2.*dphi**2. ) - E
    return to_solve

def find_dphi(E, r, phi):
    KE = E - polar_V(r, phi)
    #print("KE = ", KE)
    d_phi = np.sqrt(2*KE)/r
    return d_phi

def cartesian_V(x, y):
    return 1./2.*(x**2+y**2)+x**2*y-1./3.*y**3

def find_E_phase(y, dy):
    return 1./2.*(dy**2+y**2)-1./3.*y**3.

def find_dx(y, dy, E):
    dx = np.sqrt(2*E-dy**2-y**2+2./3.*y**3)
    return dx

def find_dy_xh(x, h):
    return np.sqrt(2.*h-x**2.)

def find_dx_yh(y, h):
    return np.sqrt(2./3.*y**3+2.*h-y**2.)


def check_y(y, h):
    check = 0
    if (2. / 3. * y ** 3 + 2. * h - y ** 2.) >= 0.:
        for i in range(len(epsilon)):
            if (2. / 3. * (y+epsilon[i]) ** 3 + 2. * h - (y+epsilon[i]) ** 2.) >= 0.:
                if (2. / 3. * (y-epsilon[i]) ** 3 + 2. * h - (y-epsilon[i]) ** 2.) >= 0.:
                    check += 0.
                else:
                    check += 1.
            else:
                check += 1.
    else:
        check += 1.
    if check != 0.:
        return 0
    else:
        return 1

# 1. For 2D x-y plane, ( 1024 X 1024 ) grids, then solve d_phi with E given.
# Initial total energy (hamiltonian)
E0 = 0.5
'''
# Initial configuration plane in polar coordinate
r0 = np.linspace(0., r_sca_reg, n)[1:]
phi0 = np.linspace(0., 2.*np.pi, n)[1:]
# Initial phase space distribution
dr0 = 0.
plane = np.array(np.meshgrid(r0, phi0))
pos = []
for i in range(len(r0)):
    for j in range(len(phi0)):
        coord  = plane[:, i, j]
        rad = coord[0]
        ang = coord[1]
        if polar_V(rad, ang) <= E0:
            dphi = find_dphi(E0, rad, ang)
            #print(dphi)
            pos.append([rad, ang, dr0, dphi])
        else:
            pass

# From polar to Cartesian coordinate:
R0, PHI0, DR0, DPHI0 = np.array(pos).T

X0 = R0*np.cos(PHI0)
Y0 = R0*np.sin(PHI0)
DX0 = -R0*np.sin(PHI0)*DPHI0
DY0 = R0*np.cos(PHI0)*DPHI0
'''

def get_coord(E):
    # OR 1. generate initial uniform coordinates in x-y plane
    x0 = np.linspace(-r_sca_reg, r_sca_reg, n+1)
    y0 = np.linspace(-r_sca_reg, r_sca_reg, n+1)

    plane = np.array(np.meshgrid(x0, y0))
    pos = []
    for i in range(len(x0)):
        for j in range(len(y0)):
            coord  = plane[:, i, j]
            x_coord = coord[0]
            y_coord = coord[1]
            r = np.sqrt(x_coord**2.+y_coord**2.)
            phi = np.arctan2(y_coord, x_coord)
                  #+np.pi/2.
            if polar_V(r, phi) <= E and r <= r_sca_reg and r != 0.:
                dphi = find_dphi(E, r, phi)
            #print(dphi)
                dx = -r * np.sin(phi) * dphi
                dy = r * np.cos(phi) * dphi
                pos.append([x_coord, y_coord, dx, dy])
            else:
                pass

    return pos

def get_phase_coord(E):
    # 2. generate initial uniform coordinates in y-dy plane
    y0 = np.linspace(-1., 2., n )
    dy0 = np.linspace(-2., 2., n )
    plane = np.array(np.meshgrid(y0, dy0))
    pos = []
    for i in range(len(y0)):
        for j in range(len(dy0)):
            coord = plane[:, i, j]
            y_coord = coord[0]
            dy_coord = coord[1]
            if find_E_phase(y_coord, dy_coord) <= E:
                dx = find_dx(y_coord, dy_coord, E)
                pos.append([0., y_coord, dx, dy_coord])
            else:
                pass
    return pos

def get_xh_coord():
    # 3. generate initial uniform coordinates in x-h plane
    x0 = np.linspace(-r_sca_reg, r_sca_reg, n )
    E0 = np.linspace(1./6., 0.4, n )
    plane = np.array(np.meshgrid(x0, E0))
    pos = []
    for i in range(len(x0)):
        for j in range(len(E0)):
            coord = plane[:, i, j]
            x_coord = coord[0]
            h_coord = coord[1]
            if (2.*h_coord-x_coord**2.) >= 0:
                dy = find_dy_xh(x_coord, h_coord)
                pos.append([x_coord, 0., 0., dy])
    return pos

def get_yh_coord():
    # 3. generate initial uniform coordinates in x-h plane
    y0 = np.linspace(-r_sca_reg, r_sca_reg, n )
    E0 = np.linspace(1./6., 0.4, n )
    plane = np.array(np.meshgrid(y0, E0))
    pos = []
    for i in range(len(y0)):
        for j in range(len(E0)):
            coord = plane[:, i, j]
            y_coord = coord[0]
            h_coord = coord[1]
            if 2./3.*y_coord**3+2.*h_coord-y_coord**2. >= 0:
                dx = find_dx_yh(y_coord, h_coord)
                pos.append([ 0., y_coord, dx, 0.])
    return pos

def ucty_coord(coord):
    if str(coord) == 'y':
        y0 = np.linspace(-0.75, 1.5, 10000)
        E0 = np.linspace(0.17, 0.4, 50)
        pos = []
        for i in range(len(y0)):
            for j in range(len(E0)):
                y_coord = y0[i]
                h_coord = E0[j]
                if check_y(y_coord, h_coord) == 1:
                    dx = find_dx_yh(y_coord, h_coord)
                    pos.append([0., y_coord, dx, 0.])
                    for k in range(len(epsilon)):
                            dx1 = find_dx_yh(y_coord+epsilon[k], h_coord)
                            pos.append([0., y_coord+epsilon[k], dx1, 0.])
                            dx2 = find_dx_yh(y_coord-epsilon[k], h_coord)
                            pos.append([0., y_coord-epsilon[k], dx2, 0.])
    elif str(coord) == 'x':
        x0 = np.linspace(-1.0, 0, 5000)
        #E0 = np.linspace(1. / 6., 0.4, 30)
        E0 = np.linspace(0.17, 0.4, 50)
        pos = []
        for i in range(len(x0)):
            for j in range(len(E0)):
                x_coord = x0[i]
                h_coord = E0[j]
                if (2. * h_coord - (abs(x_coord)+1.0e-1) ** 2.) >= 0.:
                    dy = find_dy_xh(x_coord, h_coord)
                    pos.append([x_coord, 0., 0., dy])
                    try:
                        for k in range(len(epsilon)):
                            dy1 = find_dy_xh(x_coord+epsilon[k], h_coord)
                            pos.append([x_coord+epsilon[k], 0., 0., dy1])
                            dy2 = find_dy_xh(x_coord-epsilon[k], h_coord)
                            pos.append([x_coord-epsilon[k], 0., 0., dy2])
                    except:
                        print(epsilon[k])
                        pass
    return pos

def threeD_coord(E):
    x0 = np.linspace(-r_sca_reg, r_sca_reg, n + 1)
    y0 = np.linspace(-r_sca_reg, r_sca_reg, n + 1)
    dy0 = np.linspace(-r_sca_reg, r_sca_reg, n + 1)
    pos = []
    for i in range(len(x0)):
        for j in range(len(y0)):
            for k in range (len(dy0)):
                x_coord = x0[i]
                y_coord = y0[j]
                dy_coord = dy0[k]
                check = 2.*(E+1./3.*y_coord**3.-x_coord**2*y_coord)-\
                        (x_coord**2+y_coord**2+dy_coord**2)
                if check >= 0. and x_coord**2+y_coord**2 <= 2. :
                    dx = np.sqrt(check)
                    pos.append([x_coord, y_coord, dx, dy_coord])
                else:
                    pass

    return pos


pos = threeD_coord(0.17)
#check_e = checl_conE(X0,Y0,DX0,DY0)
X0, Y0, DX0, DY0 = np.array(pos).T
# Initial distribution in Cartesian
#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.scatter(X0, Y0, DY0, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#plt.scatter(X0, Y0, s = 5.)
#plt.show()
'''
E = 1./2.*(X0**2.+DX0**2.+Y0**2.+DY0**2)+X0**2.*Y0-1./3.*Y0**3
coord = np.array([X0, Y0, DX0, DY0]).T
coord_1D = coord.ravel()

# Output the initial conditions
df = pd.DataFrame(zip(X0, Y0, DX0, DY0), columns = ["x", "y", "v_x", "v_y"])
df_1D = pd.DataFrame(coord_1D, columns=None)

#fn = "/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_300_%s.txt"%E0
fn = '/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/extra_energy/hh_initial_500_phase_%s.txt'%(E0)
#fn = "/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/hh_initial_500_coord_yh.txt"
#fn = '/Users/shuyut/Downloads/PHYS4151_Project/HH_Fractal/fractal_init_coord/hh_initial_fractal_yh.txt'
try:
    os.remove(fn)
except:
    pass

df_1D.to_csv(fn, index=False, header=False)
'''
'''
group = 100
i = 1
plt.figure(figsize=(5,5))
while i < 100:
    plt.scatter(X0[int((i-1)*group):int(i*group)],Y0[int((i-1)*group):int(i*group)], marker=1)
    i += 1
plt.xlabel('x')
plt.ylabel('y')
'''

'''
plt.scatter(Y0, E, s = 0.5)
plt.show()
'''

'''
plt.clf()

test_ind = np.arange(0,len(X0),5000)
fig = plt.figure()
ax = Axes3D(fig)

ax.scatter(X0[test_ind], Y0[test_ind], DY0[test_ind])
ax.scatter(X0[test_ind], Y0[test_ind], -DY0[test_ind])

ax.set_xlabel('X0')
ax.set_ylabel('Y0')
ax.set_zlabel('DY0')
plt.show()
'''


for tol_e in h:
    pos = threeD_coord(tol_e)
    X0, Y0, DX0, DY0 = np.array(pos).T
    coord = np.array([X0, Y0, DX0, DY0]).T
    coord_1D = coord.ravel()
    # Output the initial conditions
    df = pd.DataFrame(zip(X0, Y0, DX0, DY0), columns=["x", "y", "v_x", "v_y"])
    df_1D = pd.DataFrame(coord_1D, columns=None)
    fn = "/Users/shuyut/Downloads/PHYS4151_Project/Initial_coords/3D/hh_initial_300_%s.txt" % tol_e
    try:
        os.remove(fn)
    except:
        pass
    print("E=", tol_e, " and number of particle is", len(df))
    df_1D.to_csv(fn, index=False, header=False)
'''

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(X0, Y0, DY0)
plt.show()
'''
