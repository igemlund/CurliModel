import numpy as np
import randomComplexAddition
from numericaldiffusion import inhibitedCsgAC
from curliutil import ListDict
import cmath
from scipy.constants import N_A
import copy

class CurliFibril(object):
    """
    An object describing one Curli fibril.
    
    Attributes
    ----------
    SIGMA: constant 3.47
        average angle difference between two Curli subunits (/degrees)
    KPLUS: constant 1.4e6
        Elongation rate (/mol/s)
    UNITL: constant 4e-9
        length of one subunit
    pos: float
        distance from cell membrane (/m)
    alpha: float
        angle deviation from original orientation
    size: int
        number of subunits
    index: int
        number of fibrils created before this one
    
    Methods
    -------
    None
    
    """
    SIGMA = 3.47 #degrees
    KPLUS = 21000 #/mol s
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
    Simulates Curli fibril formation given environment and inhibitors.

    ...

    Attributes
    ----------
    dist: float
        total simulated distance from bacteria cell membrane
    CSGBRATE: constant 1.3e-13
        production rate of CsgB protein (nbr/bacteria/s)
    xsteps: int
        number of bins to keep track of fibrils in
    deltat: float
        length of one timestep
    deltax: float
        length of one bin
    cBacteria: float
        bacteria concentration

    Methods
    ----------
    timeStep:
        Steps time forward one timestep
        Args:
            None
        Returns:
            None
    getTime:
        Returns the total time simulated this far
        Args:
            None
        Returns:
            total time (float)
    getMass:
        Returns the total mass of Curli fibrils (/CsgA mass)
        Args:
            None
        Returns:
            total mass (int)
    getMassProfile:
        Returns mass distribution of Curli at different distances from the membrane
        Args:
            None
        Returns:
            1d array of mass distribution (np.array)
    getFibrils:
        Returns list of all fibrils
        Args:
            None
        Returns:
            list of CurliFibril objects (list)

    """
    
    
    def __init__(self, dist, xsteps, deltat,  inhibitors = [], cBacteria=1e12, initial_concentration = None, fibril0 = None):
        """
        Initialises the class.
        
        Args:
            dist (float): Maximum distance from cell membrane to simulate (m)
            xsteps (int): Number of bins to devide dist into
            deltat (float): Length of one timestep (s)
            inhibitors (list): List of inhibitor objects. Defaults to [].
            cBacteria (float): Bacteria concentration in environment. Defaults to 1e12.
            inital_concentration (float): Inital CsgA concentration. Defaults to None.
            fibril0 (UniformFibrilFormation): Inital fibrils. Defaults to None.
        
        Returns:
            None
        """
        self.dist = dist
        self.CSGBRATE = 1.3e-13*N_A/cBacteria #n / bac / s
        self.xsteps = xsteps
        self.deltat = deltat
        self.deltax = dist / xsteps
        self.cBacteria = cBacteria
        if initial_concentration == None:
            self.C = inhibitedCsgAC(dist, xsteps, deltat, cBacteria, how = 'uniform',c0=None, inhibitors=inhibitors)
        else:
            self.C = initial_concentration
        
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
            self.C.c = copy.deepcopy(fibril0.C.c)
    
    def getTime(self):
        return self.C.time
    def getMass(self):
        return self.totalMass
    def getFibrils(self):
        fibrils = []
        for i in self.endpointSets:
            fibrils += i.items
        return fibrils
    def getConcentration(self):
        return self.C.c
    def getMassProfile(self):
        return self.massProfile
    
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
            self.C = inhibitedCsgAC(self.dist, self.xsteps, self.deltat, self.cBacteria, self.C.how, self.C.c, self.C.inhibitors)
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
                kwrates = {'kplus':CurliFibril.KPLUS}
                for inh in self.C.inhibitors:
                    kwrates = inh.rateFunc(self.C, kwrates)
                
                mC = self.C.c
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
                self.C.c = mN / dV /N_A
        self.C.timeStep()
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CurliFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils
        
        if self.C.c > 10*self.C.time*self.C.ke:
            raise ValueError('Total fibril concentration much greater than production')

class DiffusiveFibrilFormation(UniformFibrilFormation):
    def __init__(self, dist, xsteps, deltat,  inhibitors=[], cBacteria=1e12, concentrationProfile = None, fibril0 = None):
        super().__init__(dist, xsteps, deltat, inhibitors, cBacteria,concentrationProfile, fibril0)
        if concentrationProfile == None:
            self.C = inhibitedCsgAC(dist, xsteps, deltat, cBacteria, how = 'spherical',U0=None, inhibitors=inhibitors)
        else:
            self.C = concentrationProfile

    def timeStep(self):
        for x in np.random.choice([*range(self.xsteps)], size=self.xsteps, replace=False):
            if len(self.endpointSets[x]) > 0 and self.C.U[x,0] > 0:
                mC = self.C.U[x,0]
                fN = len(self.endpointSets[x])
                xPos = x*self.deltax + self.C.R0
                dV = 4/3*np.pi*((xPos+self.deltax)**3 - xPos**3)*1e3
                mN = mC * N_A * dV
                nElongations = int(np.random.poisson(max(CsgAFibril.KPLUS* fN * mC*self.deltat,0)))
                toElongate = np.random.choice(fN, size=nElongations)
                toElongate = list(map(lambda f : self.endpointSets[x].items[f], toElongate))
                list(map(self.elongate, toElongate))
                list(map(lambda f: self.endpointSets[x].remove(f), toElongate))
                try:
                    for f in set(toElongate):
                        if f.pos > 0:
                            self.endpointSets[int(f.pos // self.deltax)].add(f)
                except IndexError:
                    raise "Dist too small. Fibril out of bounds."

                mN -= nElongations
                self.C.U[x,0] = mN /dV /N_A
        self.C.timeStep()
        nNewFibrils = np.random.poisson(self.CSGBRATE * self.deltat)
        list(map(lambda f : self.endpointSets[0].add(CsgAFibril(f)), range(self.index, self.index + nNewFibrils)))
        self.index += nNewFibrils
        
        fxp = lambda x: x*self.deltax + self.C.R0
        f = lambda x: self.C.U[x,0]*(4/3*np.pi*(fxp(x)  + self.deltax)**3 - fxp(x)**3)*1e3*N_A
        if sum(map(f, range(self.xsteps))) > 10*self.C.time*self.C.ke*N_A*self.cBacteria:
            raise ValueError('Total fibril concentration much greater than production')