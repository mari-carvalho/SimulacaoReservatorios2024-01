# Fluxo Linear com Condição de Contorno Pressão-Pressão (Dirchlet)

# Importando as Bibliotecas:

import numpy as np 
import math as mt
import matplotlib.pyplot as plt

class solucaoPressPress:
    def __init__(self, p0, pw, phi, mi, k, L, c):
        self.po = p0
        self.pw = pw
        self.phi = phi 
        self.mi = mi
        self.k = k 
        self.l = L
        self.ct = c

    def PressPress(self, x, t):

        if x == 0:
            p = self.pw
            return p # o que esta abaixo esta dentro do "else"
        elif x == self.l:
            p = self.po
            return  p 
        else:
            sum = 0 
            sum_old = 100
            err = 1000
            eppara = 1e-6
            n = 0
            while err >= eppara:
                n += 1
                i = n-1
                sum += (np.exp(-((n*np.pi/self.l)**2)*(self.k/(self.phi*self.mi*self.ct))*t)/n*np.sin(n*np.pi*x/self.l))
                err = abs((sum-sum_old)/sum)*100
                sum_old = sum
            p = (self.po - self.pw)*((x/self.l)+(2/np.pi)*sum)+self.pw
            return p
    
