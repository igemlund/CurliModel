import numpy as np
import randomComplexAddition
from numericaldiffusion import EcoliDiffusion, UniformEcoliDiffusion
from curliutil import ListDict
import cmath
from scipy.constants import N_A
import copy

class CsgAFibril:
    SIGMA = 3.47
    KPLUS = 1.4*10**6*1e24 #mole/nm3/s
    UNITL = 4.

    def __init__(self, index):
        self.pos = 0
        self.alpha = np.random.normal(0,self.SIGMA*np.pi/180)
        self.size = 1
        self.index = index
    
    def __lt__(self, other):
        if self.pos < other.pos:
            return True
        return False

class FibrilFormation(object):
    """
    Simulates Curli fibril formation. 
    
    Attributes:
        float dist : Maximum distance from cell membrane to simulate
        int xsteps : Number of bins to devide dist into
        float deltat : length of one timestep
        EcoliMonomeriDiffusion diffusion : Simulates monomeric diffusion
        float deltax : length of one bin
        list<ListDict> endpointSet : Keeps track of all fibril endpoints
        np.array massProfile : Keeps track of mass distribution in every bin
        int index : Number of fibrils in simulation, used for creating unique fibril objects. 
    
    Limits:
        deltax = dist / xsteps > CsgAFibril.UNITL 
        dist > 1000 nm (at least)
        t > 1e3
    """
    CSGBRATE = 1.3e-13*N_A*1e-12
    def __init__(self, dist, xsteps, deltat, how='uniform', what = None, fibril0 = None):
        self.dist = dist
        self.xsteps = xsteps
        self.deltat = deltat
        self.deltax = dist / xsteps
        if how == 'uniform':
            self.timeStep = self.__uniformTimeStep
            if what == None:
                self.C = UniformEcoliDiffusion(dist, xsteps, deltat)
            else:
                self.C = what
        else:
            self.timeStep = self.__nonUniformTimeStep
            self.C = EcoliDiffusion(dist, xsteps, deltat)
        
        if fibril0 == None:
            self.endpointSets = [ListDict() for _ in range(xsteps)]
            self.massProfile = np.zeros(xsteps)
            self.index = 0
            self.totalMass = 0
        else:
            self.endpointSets = copy.deepcopy(fibril0.endpointSets)
            self.massProfile = copy.deepcopy(fibril0.massProfile)
            self.index = copy.deepcopy(fibril0.index)
            self.totalMass = copy.deepcopy(fibril0.totalMass)
            self.C.U = copy.deepcopy(fibril0.C.U)
    
    def __elongate(self,fibril):
        pos0 = fibril.pos
        fibril.alpha += np.random.normal(0,fibril.SIGMA*np.pi/180)
        fibril.size += 1
        fibril.pos += fibril.UNITL*np.cos(fibril.alpha)
        diff = fibril.pos % self.deltax
        fi = lambda i : int(np.floor(i))
        if fi(fibril.pos / self.deltax) != fi(pos0 / self.deltax):
            self.massProfile[fi(pos0 / self.deltax)] += ((fibril.pos -pos0)-diff) / (fibril.pos -pos0)
            self.massProfile[fi(fibril.pos / self.deltax)] += diff / (fibril.pos -pos0)
            return
        self.massProfile[fi(pos0 / self.deltax)] += 1
        self.totalMass += 1
        return
         
    def __nonUniformTimeStep(self):
        for x in range(self.xsteps):
            if len(self.endpointSets[x]) > 0:
                mC = self.C.U[x,0]
                fN = len(self.endpointSets[x])
                xPos = x*self.deltax + self.C.R0
                dV = 4/3*np.pi*((xPos+self.deltax)**3 - xPos**3)
                mN = mC * N_A * dV
                nElongations = int(np.random.poisson(max(CsgAFibril.KPLUS * fN * mC*self.deltat,0)))
                toElongate = np.random.choice(fN, size=nElongations)
                toElongate = list(map(lambda f : self.endpointSets[x].items[f], toElongate))
                list(map(self.__elongate, toElongate))
                list(map(lambda f: self.endpointSets[x].remove(f), toElongate))
                try:
                    list(map(lambda f : self.endpointSets[int(f.pos // self.deltax)].add(f), set(toElongate)))
                except IndexError:
                    raise "Dist to small. Fibril out of bounds."
                
                mN -= nElongations
                self.C.U[x,0] = mN /dV /N_A
        self.C.timeStep()
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CsgAFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils

    def __uniformTimeStep(self):
        diff = 0
        for x in range(self.xsteps):
            if len(self.endpointSets[x]) > 0 and self.C.U > 0:
                mC = self.C.U
                fN = len(self.endpointSets[x])
                xPos = x*self.deltax + self.C.R0
                dV = 4/3*np.pi*((xPos+self.deltax)**3 - xPos**3)
                mN = mC * N_A * dV
                nElongations = int(np.random.poisson(max(CsgAFibril.KPLUS* fN * mC*self.deltat,0)))
                toElongate = np.random.choice(fN, size=nElongations)
                toElongate = list(map(lambda f : self.endpointSets[x].items[f], toElongate))
                list(map(self.__elongate, toElongate))
                list(map(lambda f: self.endpointSets[x].remove(f), toElongate))
                try:
                    list(map(lambda f : self.endpointSets[int(f.pos // self.deltax)].add(f), set(toElongate)))
                except IndexError:
                    raise "Dist to small. Fibril out of bounds."

                diff -= nElongations
        self.C.timeStep(diff)
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CsgAFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils
