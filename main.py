from test_aula import solucaoPressPress
import numpy as np 

k = 100 # md
mi = 3 # cp
phi = 0.20  
c = 130*10**-6 # (kgf/cm²)-1
p0 = 100 # kgf/cm²
pw = 70 # kgf/cm²
qw = 35 # m³std/d
L = 100 # m 
A = 2000 # m²
x = np.linspace(0,L,1000)
t = np.linspace(0,100,10)

solucaoPP = solucaoPressPress(p0,pw,phi,mi,k,L,c)
solucaoFP = solucaoFluxoPressao(po,qw,phi,mi,k,L,c)
pPP = np.zeros((len(t), len(x)))
pFP = np.zeros((len(t), len(x)))
for i in range(1,len())

p = solucao.PressPress(x,t)
print(p)