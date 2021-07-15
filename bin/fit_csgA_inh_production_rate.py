
from fibrilformation import UniformFibrilFormation
from scipy.constants import N_A 
import numpy as np
import matplotlib.pyplot as plt
import time
from inhibitors import Ainh, CsgC
np.random.RandomState(41)

def main():
    time0 = time.time()
    time1 = time.time()
    timesteps = int(1e4)
    dist = 10000
    xst = 100
    totim = 3600

    for j in [0] + [i for i in np.linspace(-13,-9,8)]:
        args = (10**j,  dist, xst, totim/timesteps,'uniform', 0)
        if j == 0:
            args = (0,  dist, xst, totim/timesteps,'uniform', 0)

        diffusionInh =  UniformFibrilFormation(dist, xst, totim/timesteps,[CsgC(*args)])
        A = []
        t = []
        for i in range(timesteps):
            if i % (timesteps // 10) == 0:
                print("{prc}% completed in {t}s. Delta T = {dt}".format(prc = str(i / timesteps *100), t = str(time.time() - time0), dt = str(time.time() - time1)))
                print(diffusionInh.C.U*10**24, diffusionInh.C.inhibitors[0].U*10**24)
                time1 = time.time()
            A.append(diffusionInh.totalMass)
            t.append(i*totim / timesteps )
            diffusionInh.timeStep()

        if j == 0:
            plt.plot(t,A, alpha = 0.2, lw=5)
        else:
            plt.plot(t,A, alpha = 0.8, label="{} /s".format(str(round(10**j*10**24/N_A *10**12, 2))))

    plt.legend()
    plt.title('Effect of CsgC Chaperone Secretion Rate on Curli Growth. \n \
        Secretion')
    plt.xlabel('Time /s')
    plt.ylabel('Mass /nbr monomer masses')
    plt.savefig('tmp4.png')
        

if __name__ == '__main__':
    main()

        
