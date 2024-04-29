import numpy as np
import matplotlib.pyplot as plt
from solver_gauss_seidel import gauss_seidel
from FTCS import FTCS
from BTCS import BTCS
from CN import CN 

# Definindo as Variáveis de Entrada:

p0 = 19000000
pw = 9000000
qw = 0.01
q0 = 100
cc = 'pp'
mi = 0.001
k = 9.869233e-14
h = 10 
phi = 0.2
c = 2.04e-9
L = 20
A = 30  
x0 = 0
xf = L
t0 = 0
tf = 100
h_t = 0.25
h_x = 0.5

#calc_FTCS = FTCS.calculate_FTCS_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
#calc_BTCS = BTCS.calculate_BTCS_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
calc_BTCS = BTCS.calculate_BTCS_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
calc_BTCS = BTCS.calculate_BTCS_fp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
calc_BTCS = BTCS.calculate_BTCS_ff(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
calc_CN = CN.calculate_CN_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)