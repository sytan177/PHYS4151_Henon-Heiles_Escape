import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.integrate import odeint
from scipy.optimize import fsolve


# Define Henon-Heiles potential
def V_hh(x,y):
    v = 1 / 2 * (x ** 2 + y ** 2) + (x ** 2 * y - y ** 3 / 3)
    return v

def H(x, y, xdot, ydot):
    h = V_hh(x,y) - 1/2*(xdot**2+ydot**2)
    return h

# Henon-Heiles Potantial Plot
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

x = np.linspace(-1.3,1.3, 100)
y = np.linspace(-1.1,1.4, 100)

X, Y = np.meshgrid(x, y)
pot = V_hh(X, Y)

level = [*np.arange(min(pot.flatten())-0.2, 1./6., 0.03),
         *np.arange(1./6., max(pot.flatten()), 0.10)]
cp = ax.contour(x, y, pot, levels = list(level),
                colors = 'black',  linewidths = 0.5)
#cp = ax.contour(x, y, pot, levels = list([0.17, 0.18, 0.19, 0.20, 0.21]),
#                colors = 'black',  linewidths = 0.5)
# Energy of Escape and Saddle Point
x_saddle = [0, -np.sqrt(3)/2, np.sqrt(3)/2]
y_saddle = [1., -1./2., -1./2.]
cp = ax.contour(x, y, pot, levels = [1./6.],
                colors = 'red',  linewidths = 2.)

plt.scatter(x_saddle, y_saddle, color = 'blue', s = 60, zorder=2)
plt.scatter(0., 0., s = 60, zorder=2)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.show()


'''
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, pot, c=pot, cmap='Blues')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V(x,y)')
plt.show()
'''
