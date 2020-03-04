import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
from mpl_toolkits.mplot3d import Axes3D
import random as rnd

def iniconf(box_size, dist, n_b):
    limit = 0.5 * (box_size - dist)
    x = np.linspace(-limit, limit, n_b)
    y = np.linspace(-limit, limit, n_b)
    flatX, Y = np.meshgrid(x, y)
    flatX = np.matrix.flatten(flatX)
    Y = np.matrix.flatten(Y)
    
    X, Z = np.meshgrid(x, Y)
    Z = np.matrix.flatten(Z)
    X, Y = np.meshgrid(x, flatX)
    X = np.matrix.flatten(X)
    Y = np.matrix.flatten(Y)
    
    return X, Y, Z

def hardsphere(rij, xij, yij, zij):
    dlr = 50.0
    dla = 49.0
    a2 = (dlr / (dlr - dla)) * (dlr / dla)**(dla / (dlr - dla))
    tem = 1.47370
    
    if rij < (dlr / dla)**(1.0 / (dlr - dla)):
        uij = (a2 / tem) * ((1.0 / rij)**dlr - (1.0 / rij)**dla) + 1.0 / tem
        fij = dlr * (1.0 / rij)**(dlr + 1.0) - dla * (1.0 / rij)**(dla + 1.0)
        fij = (a2 / tem) * fij
    else:
        uij = 0.0
        fij = 0.0
        
    fxij = fij * xij / rij
    fyij = fij * yij / rij
    fzij = fij * zij / rij
    
    return fxij, fyij, fzij, uij

def force(x, y, z, n, box_l, rc):
    ener = 0.0
    fx = np.zeros(n) 
    fy = np.zeros(n)
    fz = np.zeros(n)
    
    for i in range(n - 1):
        for j in range(i + 1, n):
            xij = x[i] - x[j]
            yij = y[i] - y[j]
            zij = z[i] - z[j]
            xij = xij - box_l * np.rint(xij / box_l)
            yij = yij - box_l * np.rint(yij / box_l)
            zij = zij - box_l * np.rint(zij / box_l)
            rij2 = xij**2 + yij**2 + zij**2
            rij = sqrt(rij2)
            
            if rij < rc:
                fxij, fyij, fzij, uij = hardsphere(rij, xij, yij, zij)
                ener = ener + uij
                fx[i] = fx[i] + fxij 
                fy[i] = fy[i] + fyij
                fz[i] = fz[i] + fzij
                fx[j] = fx[j] - fxij 
                fy[j] = fy[j] - fyij
                fz[j] = fz[j] - fzij
    
    return fx, fy, fz, ener

def position(x, y, z, fx, fy, fz, pbc, boxl, N):
    deltat = 0.00001
    sigma = sqrt(2.0 * deltat)
    dx = np.ones(N)
    dy = np.ones(N)
    dz = np.ones(N)
    rng = np.random.default_rng()

    for i in range(N):
        dx[i] = sigma * rng.standard_normal()
        dy[i] = sigma * rng.standard_normal()
        dz[i] = sigma * rng.standard_normal()

    xf = x + (fx * deltat) + dx
    yf = y + (fy * deltat) + dy
    zf = z + (fz * deltat) + dz
    
    if pbc > 0.0:
        xf = xf - boxl * np.rint(xf / boxl)
        yf = yf - boxl * np.rint(yf / boxl)
        zf = zf - boxl * np.rint(zf / boxl)
    
    return xf, yf, zf

def graficar(X, Y, Z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, c = 'b', marker = 'o') 
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')