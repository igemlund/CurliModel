import numpy as np
import cmath
import os

class randomComplexAddition:

    def __init__(self, sigma):
        self.sigma = sigma*np.pi / 180
        self.path = "../data/randomComplexNumbers"
        if not os.path.exists('../data'):
            os.makedirs('../data')
        if not os.path.exists(self.path):
            self.randomArray = np.array([[cmath.rect(1, np.random.normal(self.sigma, 1)) for _ in range(1000)]])
            self.__addRandomDist(self.randomArray[0], 1)
        else:
            self.randomArray = np.loadtxt(self.path, dtype=complex)
    
    def __addRandomDist(self, lis, ind):
        n = 1 << ind - 1
        for i in range(n):
            lis += np.array([cmath.rect(1, np.random.normal(self.sigma, 1)) for _ in range(len(lis))])
        if self.randomArray.shape[0] >= ind:
            self.randomArray = np.vstack((self.randomArray, lis))
        else:
            self.randomArray[ind] = lis
        np.savetxt(self.path, self.randomArray)

    def __randomNumber(self, n):
        while self.randomArray.shape[0] <= n:
            self.__addRandomDist(self.randomArray[-1], self.randomArray.shape[0])
        
        return self.randomArray[n, np.random.randint(0,self.randomArray.shape[1])]

    def getRandomNumber(self, n, size, alhpa0):
        if 0 == n:
            return 0
        index = n.bit_length() - 1
        x = self.__randomNumber(index) * size 
        return x + self.getRandomNumber(n - (1 << index), size, cmath.phase(x)) 
