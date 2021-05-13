#%%
from numericaldiffusion import EcoliMonomericDiffusion
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt

timesteps = int(1e2)
diffusion = EcoliMonomericDiffusion(500, 50, 5*3600/timesteps)

for i in range(timesteps):
    diffusion.timeStep()

xstag = np.linspace(0., 500, 51)
x = xstag[:-1] + diffusion.deltax/2

plt.plot(x, diffusion.U[:,0]*10**18)
plt.show()
# %%
