import fibrilformation as ff

from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
import time
from curliutil import checkMakeDir
np.random.RandomState(41)

time0 = time.time()
time1 = time.time()
timesteps = int(1e3)
dist = 1e-5
xst = 500
diffusion = ff.UniformFibrilFormation(dist, xst, 10*3600/timesteps)

for i in range(timesteps):
    if i % (timesteps // 10) == 0:
        print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
        print(f'Total mass: {diffusion.totalMass}')
        time1 = time.time()
    diffusion.timeStep()

fibrils = []
for i in diffusion.endpointSets:
    fibrils += i.items

xstag = np.linspace(0., dist, xst + 1)
x = (xstag[:-1] + diffusion.deltax/2) / 1e-9
mu = chr(956)
x, diffusion.massProfile = x[:-200], diffusion.massProfile[:-200]
fig, axs = plt.subplots(2,2, figsize=(12,9))
fig.suptitle("After 10 h of growth")
axs[0,0].scatter(x, diffusion.massProfile)
axs[0,0].set_title('Mass distribution')
axs[0,0].set_xlabel('Distance from membrane (/nm)')
axs[0,0].set_ylabel('Mass /CsgA mass')
#axs[0,1].plot(x, diffusion.C.U[:-200,0] / 1e-6)
axs[0,1].set_title('Monomer concentration')
axs[0,1].set_ylabel(f"Concentration /{mu}M")
axs[0,1].set_yscale('log')
axs[1,0].hist(list(map(lambda i : i.pos*1e9, fibrils)))
axs[1,0].set_title('Endpoint distrubution')
axs[1,0].set_yscale("log")
axs[1,0].set_xlabel('Distance from membrane /nm')
axs[1,1].hist(list(map(lambda i : i.size, fibrils)))
axs[1,1].set_title('Size distribution')
axs[1,1].set_yscale("log")
axs[1,1].set_xlabel('Fibril size /CsgA masses')

checkMakeDir('../figures/')
plt.savefig("../figures/Curli_After_10h_spherical.png")

