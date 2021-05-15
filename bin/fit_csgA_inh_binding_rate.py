from os import times
from numpy.lib.function_base import diff
from numericaldiffusion import UnifCurliFormWAInh
from fibrilformation import FibrilFormation
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
import time
np.random.RandomState(41)

time0 = time.time()
time1 = time.time()
timesteps = int(1e2)
dist = 2500
xst = 150
totim = 5*3600
for j in np.linspace(-1, 1, 6):
    diffusion = FibrilFormation(dist, xst, totim/timesteps, how='uniform', what=UnifCurliFormWAInh(dist, xst, totim / timesteps, \
        10**j*1e24, 0.3*1e-3*1e-24))
    
    A = []
    t = []
    for i in range(timesteps):
        if i % (timesteps // 10) == 0:
            print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
            time1 = time.time()
        if i % (timesteps / 1000):
            A.append(diffusion.totalMass)
            t.append(i*totim / timesteps)
        diffusion.timeStep()

    plt.plot(t,A, label=str(j))

diffusion = FibrilFormation(dist, xst, totim/timesteps, how='uniform')
    
A = []
t = []
for i in range(timesteps):
    if i % (timesteps // 10) == 0:
        print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
        time1 = time.time()
    if i % (timesteps / 1000):
        A.append(diffusion.totalMass*0.3)
        t.append(i*totim / timesteps)
    diffusion.timeStep()

plt.plot(t,A, linewidth=5, alpha = 0.2, c = 'blue')

plt.legend()
plt.savefig('tmp.png')
    
    

    
