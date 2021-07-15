import numpy as np
from scipy.constants import N_A
from numericaldiffusion import Inh

class Ainh_fixed:
    U = 0.3*1e-3*1e24 #Concentration in mol/nm3
    def __init__(self, deltat, KS):
        self.KS = KS
        self.deltat = deltat
    
    def timeStep(self, other):
        other.U -= max(self.KS*self.U*other.U*self.deltat, 0)
    def rateFunc(self,other, kwrates):
        return kwrates
    
class Ainh(Inh):
    """The inhibitor binding to CsgA. 
    The binding constant KS was estimated by assuming first order binding to the CsgA monomers \
        and then computing the equilibrium constant from previous published results. 
    """
    #KS = 10**24.5
    def __init__(self, pAinh, dist, xsteps, deltat, how, U0, KS):
        super().__init__(dist, xsteps, deltat, pAinh, how, U0)
        self.KS = KS

    def bindingfunction(self, monC):
        """The function declaring how the inhibitor will affect the CsgA concentration. \
            The function edits the concentration inplace. 
            Note that the function only is implemented for uniform concentrations. 

        Args:
            monC (CsgADiffusion): An object describing the CsgA monomers.
        """
        if self.how == 'uniform' and monC.how == 'uniform':
            self.U -= max(self.KS*self.U*monC.U*self.deltat,0)
            monC.U -= max(self.KS*self.U*monC.U*self.deltat,0)

class CsgC(Inh):
    """The CsgC chaperone. 
    """
    def __init__(self, pAinh, dist, xsteps, deltat, how, U0):
        super().__init__(dist, xsteps, deltat, pAinh, how, U0)
    
    def rateFunc(self, other, kwrates):
        """A function describing the affect the CsgC concentration will have on the reaction rates.\
            Values computed by fitting an exponential function to previously published data. 
            The CsgC is assumed to not degrade or diffuse away. -_('H')_-


        Args:
            other (CsgADiffusion): Has no affect on the CsgC 
            kwrates (dict<string, float>): Rate constants

        Returns:
            dict<string,float>: New rate constants
        """
        kwrates['kplus'] = np.sqrt(np.exp(-88*self.U*1e+30))*kwrates['kplus']
        return kwrates

