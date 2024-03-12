# classes:
#Dirchlet:
import numpy as np
import scipy as sp 

class solucaoPressPress:
    def __init__(self, p0, pw, phi, mi, k, L, c):
        self.po = p0
        self.pw = pw
        self.phi = phi 
        self.mi = mi
        self.k = k 
        self.L = L
        self.c = c

    def PressPress(self, x, t):

        p = np.zeros((len(x), len(t)))
        sum = 0 
        sum_old = 100
        err = 100
        eppara = 10**-3
        n = 1 
        while err > eppara:
            sum = sum + np.exp(- (n*np.pi/self.L)**2*(self.k/(self.phi*self.mi*self.c))*t)*np.sin(n*np.pi*x/self.L)
            n = n + 1
            err = abs((sum - sum_old)/sum) * 100
            sum_old = sum 
        p = (self.po - self.pw) * (x/self.L + 2/np.pi * sum) 

        return p