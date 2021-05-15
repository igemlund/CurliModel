from fibrilformation import FibrilFormation
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
import time
from curliutil import checkMakeDir
np.random.RandomState(41)

time0 = time.time()
time1 = time.time()
timesteps = int(1e2)
dist = 2500
xst = 150
diffusion = FibrilFormation(dist, xst, 10.*3600/timesteps, how='uniform')

for i in range(timesteps):
    if i % (timesteps // 10) == 0:
        print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
        time1 = time.time()
    diffusion.timeStep()

fibrils = []
for i in diffusion.endpointSets:
    fibrils += i.items

xstag = np.linspace(0., dist, xst + 1)
x = xstag[:-1] + diffusion.deltax/2

fig, axs = plt.subplots(2,2, figsize=(12*2,9*2))
fig.suptitle("After 15 h of growth")
axs[0,0].scatter(x, diffusion.massProfile)
axs[0,0].set_title('Mass distribution')
#axs[0,1].plot(x, diffusion.diffusion.U[:,0]*1e30)
#axs[0,1].set_title('Monomer concentration')
#axs[0,1].set_ylabel("Concentration /mu M")
print(diffusion.C.U)
axs[1,0].hist(list(map(lambda i : i.pos, fibrils)))
axs[1,0].set_title('Endpoint distrubution')
axs[1,0].set_yscale("log")
axs[1,1].hist(list(map(lambda i : i.size, fibrils)))
axs[1,1].set_title('Size distribution')
axs[1,1].set_yscale("log")
list(map(lambda i : i.set_xlabel("Distance from membrane / nm"), axs.flatten()))

checkMakeDir('../figures/')
plt.savefig("../figures/Curli_After_10h_uniform.png")

