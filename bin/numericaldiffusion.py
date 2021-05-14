# -*- coding: utf-8 -*-
import numpy as np
import os
from scipy.constants import N_A

class NumericalDiffusion:
    """
    Simulates diffusion with constant D in one dimension from start to stop with the timestep deltat and xstep number of bins.
    Implemented for spherical symetry but will be extended to linear symmetry. 

    Variables:
        float start
        float stop
        int xsteps
        int deltat
        float D
        np.array((n,1)) U
        strinng how
    """
    A = np.array([])
    def __init__(self, start, stop, xsteps, deltat, D, U = None, how='spherical'):
        self.start = start
        self.stop = stop
        self.xsteps = xsteps
        self.deltat = deltat
        self.D = D
        self.deltax = (stop - start)/xsteps
        self.time = 0
        self.path = "../data/DiffusionArrays/Array_" + '_'.join([str(x) for x in [self.deltax, deltat, xsteps, D]])
        self.U = U
        if self.U == None:
            self.U = np.zeros((xsteps, 1))
        
        self.makeArray(xsteps, deltat, self.deltax, D)
    
    def makeArray(self, xsteps, deltat, deltax, D):
            K = D*deltat / deltax**2
            r = lambda i : i*deltax + self.start
            
            xn = lambda i : 1 + 2*K/r(i)**2 *(r(i)**2 + deltax**2)
            xp = lambda i : -(r(i) + deltax)**2 * K/r(i)**2
            xm = lambda i : -(r(i) - deltax)**2 *K/r(i)**2

            A = np.diag([xn(i) for i in range(xsteps)], 0)
            A += np.diag([xp(i) for i in range(xsteps - 1)], 1)
            A += np.diag([xm(i) for i in range(xsteps - 1)], -1)
            A[0,0] = 1 + K/r(0)**2 * (r(0)**2 + deltax**2)
            A[-1,-1] = 1 + K/r(xsteps)**2 *(r(xsteps)**2 + deltax**2)
            A = np.linalg.inv(A)
            self.A = A

    def timeStep(self):
        #Steps one timme step forward 
        self.U = np.matmul(self.A, self.U)
        self.time += self.deltat


class MonomericDiffusion(NumericalDiffusion):
    """
    Simulates diffusion where more particles are generated at start.
    """
    def __init__(self, start, stop, xsteps, deltat, D, ke):
        super().__init__(start, stop, xsteps, deltat, D, U = None, how='spherical')
        self.ke = ke
    

    def timeStep(self):
        self.U[0,0] += self.ke*self.deltat/N_A*3/4/np.pi/(self.start**2*self.deltax)
        super().timeStep()

class EcoliMonomericDiffusion(MonomericDiffusion):
    """
    Siulates diffusion for CsgA monomers. All values in nm
    """
    ke = 1e-10 * N_A * 1e-12
    R0 = 380
    D = 82114
    
    def __init__(self, dist, xsteps, deltat):
        super().__init__(self.R0, dist + self.R0, xsteps, deltat, self.D, self.ke)