# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:15:09 2015

@author: Binoy Pilakkat
"""
import time
import numpy as np
from scipy.sparse import linalg
from scipy import sparse
import matplotlib.pyplot as plt
from numba import autojit


def matSolve(matA, matB):
    """
    Solves the matrixes using sparse matrix implementation
    """
    spA = sparse.csr_matrix(matA)
    curSol = linalg.spsolve(spA, matB)
    return curSol
numSol = autojit()(matSolve)
numSol.func_name = "numSol"


def blsolve(nx, ny, dx, dy, U0, V0):
    """
    nx - Number of points along x - axis
    ny - Number of points along y - axis
    dx - Grid spacing along x
    dy - Grid spacing along y
    U0,V0 - Initial velocity conditions

    """
    matU = np.zeros((ny, nx), dtype='float64')
    matV = np.zeros((ny, nx))
    matU[-1, :] = U0
    ################################################
    # bottom boundary
    matU[0, :] = 0
    matV[0, :] = 0
    ################################################
    # Left wall
    matU[:, 0] = U0
    matV[:, 0] = V0

    progress = 0.0
    for coln_u in np.arange(1, nx):                 # iterate through columns
        progress += 1
        matA = np.zeros((ny, ny), dtype='float64')

        matA[0, 0] = 1
        matA[-1, -1] = 1
        matB = np.zeros((ny, 1))
        matB[-1, :] = matU[-1, coln_u]
        for rown_u in np.arange(1, ny-1):           # iterate through rows
            matA[rown_u, rown_u-1] = -(matV[rown_u][coln_u-1]/dy) - (nu/dy**2)
            matA[rown_u, rown_u] = 1 * ((matU[rown_u, coln_u-1]/dx) + \
                                        matV[rown_u, coln_u-1]/dy + \
                                        (2*nu)/dy**2)
            matA[rown_u, rown_u+1] = -nu/dy**2
            matB[rown_u] = 1.0*(matU[rown_u, coln_u-1])**2/dx
        matU[:, coln_u] = numSol(matA, matB)
        for rown_u in np.arange(1, ny): # v calculation
            matV[rown_u, coln_u] = matV[rown_u-1, coln_u] \
                        + ((dy/dx)*(matU[rown_u, coln_u-1] - \
                        matU[rown_u, coln_u]))
    return [matU,matV]

if __name__ == '__main__':
    #############################################################################
    ## input Parameter
    nu = 1.1886*10**-3                  # kinematic viscosity
    Re = 500                            # REynold's number
    L = 1.50                            # Length of theplate
    bl = 5*L/np.sqrt(Re)                # Height of the domain
    H = 2.0 * bl                        # DOmain height
    nx = 500                            # Discretization along x
    ny = 500                            # Discretization along y
    dx = L/(nx-1)
    dy = H/(ny-1)
    U0 = (Re*nu)/L   # ms-1
    V0 = 0
    ###########################################################################
    start = time.time()
    x = np.linspace(0, L, nx)
    y = np.linspace(0, H, ny)
    [matU, matV] = blsolve(nx, ny, dx, dy, U0, V0)
    it = np.nditer(matU, flags=['multi_index'])
    u99 = []
    for val in it:  # Locating boundary layer
        if not val/U0 > 0.99:
            u99.append(it.multi_index)
    ###########################################################################
    # plt.quiver(x,y,matU,matV,width = 0.001)
    layer = {}
    for (r, c) in u99:
        layer[c] = r
    #    plt.plot(x[indx[1]],y[indx[0]],'^')
    ly = [y[c] for (r, c) in layer.items()]
    lx = [x[r] for (r, c) in layer.items()]
    yblass = 5 / np.sqrt(U0*x[1:]/nu) * x[1:]
    print 'total timn:', round(time.time()-start, 2), 'seconds'
    fig, ax = plt.subplots()
    ax.plot(x[1:], yblass, label='Blassius Solution')
    print "sd:", np.std(ly-yblass)
    plt.plot(lx, ly, label='Numerical Solution')
    plt.legend()
    plt.show()
