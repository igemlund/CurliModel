from numericaldiffusion import CurliAinhSecretion
from fibrilformation import FibrilFormation
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
np.random.RandomState(41)

time0 = time.time()
time1 = time.time()
timesteps = int(1e4)
dist = 10000
xst = 1000
totim = 1*3600

diffusion = FibrilFormation(dist, xst, totim/timesteps, how='uniform')
    
Ae = []
te = []
for i in range(timesteps):
    if i % (timesteps // 10) == 0:
        print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
        time1 = time.time()
    if (i + 1) % (timesteps / 1000):
        Ae.append(diffusion.totalMass)
        te.append(i*totim / timesteps)
    diffusion.timeStep()

for j in np.linspace(-10, -5, 5):
    diffusionInh =  FibrilFormation(dist, xst, totim/timesteps, how='uniform', \
        what=CurliAinhSecretion(dist, xst, totim/timesteps, \
        10**j), fibril0 = diffusion)
    A = []
    t = []
    for i in range(timesteps):
        if i % (timesteps // 10) == 0:
            print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
            print(diffusionInh.C.U*10**24)
            time1 = time.time()
        if (i+1) % (timesteps / 1000):
            A.append(diffusionInh.totalMass)
            t.append(i*totim / timesteps + totim)
        diffusionInh.timeStep()

    plt.plot(t,A, alpha = 0.8, label="{} /s".format(str(10**j*10**24 /N_A *10**12)))

for i in range(timesteps):
    if i % (timesteps // 10) == 0:
        print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
        time1 = time.time()
    if (i+1) % (timesteps / 1000):
        Ae.append(diffusion.totalMass)
        te.append(i*totim / timesteps + totim)
    diffusion.timeStep()

plt.plot(te,Ae, linewidth=5, alpha = 0.2, c = 'blue', label ='0')

plt.legend()
plt.title('Effect of CsgA Inhibitor Secretion Rate on Curli Growth. \n \
    Secretion starting after half the time')
plt.xlabel('Time /s')
plt.ylabel('Mass /nbr monomer masses')
plt.savefig('tmp.png')
    
    

    
