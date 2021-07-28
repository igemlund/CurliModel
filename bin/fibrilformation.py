import numpy as np
import randomComplexAddition
from numericaldiffusion import inhibitedCsgAC
from curliutil import ListDict
import cmath
from scipy.constants import N_A
import copy

class CsgAFibril:
    SIGMA = 3.47 #degrees
    KPLUS = 1.4*10**6 #/mol s
    UNITL = 4e-9

    def __init__(self, index):
        self.pos = 0
        self.alpha = np.random.normal(0,self.SIGMA*np.pi/180)
        self.size = 1
        self.index = index
    
    def __lt__(self, other):
        if self.pos < other.pos:
            return True
        return False

class UniformFibrilFormation(object):
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
        dt > 1e3
    """
    
    def __init__(self, dist, xsteps, deltat,  inhibitors = [], cBacteria=1e12, concentrationProfile = None, fibril0 = None, how='uniform'):
        self.dist = dist
        self.CSGBRATE = 1.3e-13*N_A/cBacteria #n / bac / s
        self.xsteps = xsteps
        self.deltat = deltat
        self.deltax = dist / xsteps
        self.cBacteria = cBacteria
        if concentrationProfile == None:
            self.C = inhibitedCsgAC(dist, xsteps, deltat, cBacteria, how = how,U0=None, inhibitors=inhibitors)
        else:
            self.C = concentrationProfile
        
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
    
    def elongate(self,fibril):
        pos0 = fibril.pos
        fibril.alpha += np.random.normal(0,fibril.SIGMA*np.pi/180)
        fibril.size += 1
        fibril.pos += fibril.UNITL*np.cos(fibril.alpha)
        diff = fibril.pos % self.deltax
        """
        while fibril.pos / self.deltax >= self.xsteps:
            self.endpointSets = copy.deepcopy(self.endpointSets) + [ListDict() for _ in range(self.xsteps)]
            massProfiletmp = np.zeros(self.xsteps*2)
            massProfiletmp[:self.xsteps] = self.massProfile
            self.massProfile = massProfiletmp
            self.xsteps *= 2
            self.C = inhibitedCsgAC(self.dist, self.xsteps, self.deltat, self.cBacteria, self.C.how, self.C.U, self.C.inhibitors)
            """
        fi = lambda i : int(np.floor(i))
        self.totalMass += 1
        if fi(fibril.pos / self.deltax) != fi(pos0 / self.deltax):
            self.massProfile[fi(pos0 / self.deltax)] += (abs(fibril.pos -pos0)-diff) / abs(fibril.pos -pos0)
            self.massProfile[fi(fibril.pos / self.deltax)] += diff / abs(fibril.pos -pos0)
            return
        else:
            self.massProfile[fi(pos0 / self.deltax)] += 1
        return
         
    def timeStep(self):
        for x in np.random.choice([*range(self.xsteps)], size=self.xsteps, replace=False):
            if len(self.endpointSets[x]) > 0:
                kwrates = {'kplus':CsgAFibril.KPLUS}
                for inh in self.C.inhibitors:
                    kwrates = inh.rateFunc(self.C, kwrates)
                
                mC = self.C.U
                fN = len(self.endpointSets[x])
                dV = 1/self.cBacteria #The volume occupied by one bacteria in 1e12 bac/dm3 sol.
                mN = mC * N_A * dV
                nElongations = int(np.random.poisson(max(kwrates['kplus'] * fN* mC*self.deltat,0)))
                toElongate = np.random.choice(fN, size=nElongations)
                toElongate = list(map(lambda f : self.endpointSets[x].items[f], toElongate))
                list(map(self.elongate, toElongate))
                list(map(lambda f: self.endpointSets[x].remove(f), toElongate))
                try:
                    for f in set(toElongate):
                        if f.pos > 0:
                            self.endpointSets[int(f.pos // self.deltax)].add(f)
                except IndexError:
                    raise "Dist to small. Fibril out of bounds."
                
                mN -= nElongations
                self.C.U = mN / dV /N_A
        self.C.timeStep()
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CsgAFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils
        
        if self.C.U > 10*self.C.time*self.C.ke:
            raise ValueError('Total fibril concentration much greater than production')

class DiffusiveFibrilFormation(UniformFibrilFormation):
    def __init__(self, dist, xsteps, deltat,  inhibitors=[], cBacteria=1e12, concentrationProfile = None, fibril0 = None):
        super().__init__(dist, xsteps, deltat, inhibitors, cBacteria,concentrationProfile, fibril0, how='spherical')

    def timeStep(self):
        for x in np.random.choice([*range(self.xsteps)], size=self.xsteps, replace=False):
            if len(self.endpointSets[x]) > 0 and self.C.U[x,0] > 0:
                mC = self.C.U[x,0]
                fN = len(self.endpointSets[x])
                xPos = x*self.deltax + self.C.R0
                dV = 4/3*np.pi*((xPos+self.deltax)**3 - xPos**3)*1e3
                mN = mC * N_A * dV
                nElongations = int(np.random.poisson(max(CsgAFibril.KPLUS* fN * mC*self.deltat,0)))
                if nElongations > mN:
                    nElongations = int(mN)
                toElongate = np.random.choice(fN, size=nElongations)
                toElongate = list(map(lambda f : self.endpointSets[x].items[f], toElongate))
                list(map(self.elongate, toElongate))
                list(map(lambda f: self.endpointSets[x].remove(f), toElongate))
                try:
                    for f in set(toElongate):
                        if f.pos > 0:
                            self.endpointSets[int(f.pos // self.deltax)].add(f)
                except IndexError:
                    raise "Dist to small. Fibril out of bounds."

                mN -= nElongations
                self.C.U[x,0] = mN /dV /N_A
        self.C.timeStep()
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CsgAFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils
        
        if sum(self.C.U)[0] > 10*self.C.time*self.C.ke:
            raise ValueError('Total fibril concentration much greater than production')