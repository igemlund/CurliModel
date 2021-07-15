from os import times
from numpy.lib.function_base import diff
from inhibitors import Ainh_fixed
from fibrilformation import UniformFibrilFormation
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
from time import time
np.random.RandomState(41)

totime = 3*3600
dist = 10000
xst = 150
timesteps = int(1e2)
deltat =  totime / timesteps

def plot(inhibitors, mainplot=None, label=None):
    A = []
    t = []
    time0 = time()
    time1 = time()
    diffusion = UniformFibrilFormation(dist, xst, deltat,  inhibitors = inhibitors)
    for i in range(timesteps):
        
        if i % (timesteps // 10) == 0:
            print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time() - time0), dt = str(time() - time1)))
            time1 = time()
        if i % (timesteps / 10):
            A.append(diffusion.totalMass)
            t.append(i*totime / timesteps)
        diffusion.timeStep()
    A.append(diffusion.totalMass)
    t.append(totime)
    plt.plot(t,A, linewidth=5, alpha = 0.5, label=label)
    if not mainplot == None:
        plt.plot(mainplot[0], mainplot[1], alpha = 0.2, label = 'Main')
    plt.legend()
    return t, A

t0, A0 = plot([])
M0 = A0[-1]
M = 0
failstop = 20
upper = 1e30
lower = 1e20

while failstop > 0 and abs(M0 - M) / M0 > 0.1:
    KS = (upper + lower) / 2
    print(f"\n################ KS = {KS} ################\n")
    t, A = plot([Ainh_fixed(totime/timesteps, KS)], [t0, A0], label=str(KS))
    if A[-1] - M0*0.7 > 0:
        lower = KS
    elif A[-1] - M0*0.7 < 0:
        upper = KS
    else:
        break
    failstop -= 1

print(f'Estimated KS: {KS}')
    

    
