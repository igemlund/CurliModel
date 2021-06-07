import numpy as np
from scipy.constants import N_A
from numericaldiffusion import Inh

class Ainh(Inh):
    KS = 10**24.5
    def __init__(self, pAinh, dist, xsteps, deltat, how, U0):
        super().__init__(dist, xsteps, deltat, pAinh, how, U0)

    def bindingfunction(self, monC):
        if self.how == 'uniform' and monC.how == 'uniform':
            self.U -= max(self.KS*self.U*monC.U*self.deltat,0)
            monC.U -= max(self.KS*self.U*monC.U*self.deltat,0)

class CsgC(Inh):
    def __init__(self, pAinh, dist, xsteps, deltat, how, U0):
        super().__init__(dist, xsteps, deltat, pAinh, how, U0)
    
    def rateFunc(self, other, kwrates):
        kwrates['kplus'] = np.sqrt(np.exp(-88*self.U*1e+30))*kwrates['kplus']
        return kwrates
