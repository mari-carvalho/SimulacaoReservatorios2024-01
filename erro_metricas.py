import numpy as np
import matplotlib.pyplot as plt
from solver_gauss_seidel import gauss_seidel
from fluxo_linerar_dirchlet_aula import solucaoPressPress
from FTCS import FTCS
from BTCS import BTCS
from CN import CN 

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

n_x = (xf-x0)/(h_x)
n_t = (tf-t0)/(h_t)

x = np.linspace(0,L,int(n_x)+1) # de 0 a 10, com 1000 elementos
t = np.linspace(0,tf,int(n_t)+1) # de 0 a 1000, com 10 elementos
tam1 = len(x)
print(tam1)
# não importa o tamanho de cada vetor, são independentes e formam uma matriz de i linhas e j colunas 

p_dirchlet_matriz = np.zeros((len(t), len(x))) # matriz de 10 linhas e 1000 colunas 
tam2 = len(p_dirchlet_matriz)
print(tam2)

solucao_dirchlet = solucaoPressPress(p0, pw, phi, mi, k, L, c) # chamando a classe 

for i in range(len(t)): # vai percorrer todas as linhas, todos os tempos
    p_dirchlet_matriz[i, :] = solucao_dirchlet.PressPress(x, t[i]) # chamando a função # para cada valor de tempo (linhas), vai calcular um vetor de pressões com os valores de posição x, a linha da vez (tempo) por todas as colunas (posições x)
print(p_dirchlet_matriz)

for i in range(len(t)):
    plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-') # vai plotar todo o vetor de x por cada um dos vetores de pressão gerados, na linha da vez (tempo) vai plotar todas as colunas (vetor de posições x) por todas as linhas (vetor de pressões gerado)
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()

x_ex, y_ex = x, p_dirchlet_matriz

x_calc, y_calc = FTCS.calculate_FTCS_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)

# Cálculo do Erro:

#lembrar de pegar uma linha da matriz











#calc_FTCS = FTCS.calculate_FTCS_fp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
#calc_FTCS = FTCS.calculate_FTCS_ff(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)