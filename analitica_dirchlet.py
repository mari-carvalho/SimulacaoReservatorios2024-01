import numpy as np 
import math as mt
import matplotlib.pyplot as plt
from fluxo_linerar_dirchlet import solucaoPressPress

def calculate_analitica_dirchlet(h_t):
        
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
    h_x = 0.5

    n_x = (xf-x0)/(h_x)
    n_t = (tf-t0)/(h_t)
    print('nx', n_x)

    x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
    t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos
    p = np.zeros((int(n_t)+1,int(n_x)+1))

    # Alimentando os vetores:
    for i in range(len(x)):
        if i == 0:
            x[i] = x0
        elif i == 1:
            x[i] = i*(h_x/2)
        elif i == len(x):
            x[i] = L 
        elif i == len(x)-1:
            print(i)
            x[i] = x[i-1] + (h_x/2)
        else:
            x[i] = x[i-1] + h_x
                
    for i in range(len(t)):
        if i == 0:
            t[i] = t0
        elif i == len(t):
            t[i] = tf
        else:
            t[i] = i*h_t      
    print('x', x)
    print('t', t)
    print('p', p)
    tam1 = len(x)
    print(tam1)
    # não importa o tamanho de cada vetor, são independentes e formam uma matriz de i linhas e j colunas 

    p_dirchlet_matriz = np.zeros((len(t), len(x))) # matriz de 10 linhas e 1000 colunas 
    tam2 = len(p_dirchlet_matriz)
    print(tam2)

    solucao_dirchlet = solucaoPressPress(p0, pw, phi, mi, k, L, c) # chamando a classe 

    for i in range(len(t)): # vai percorrer todas as linhas, todos os tempos
        for j in range(len(x)):
            if i == 0:
                p_dirchlet_matriz[i,j] = p0
            else:
                p_dirchlet_matriz[i, j] = solucao_dirchlet.PressPress(x[j], t[i]) # chamando a função # para cada valor de tempo (linhas), vai calcular um vetor de pressões com os valores de posição x, a linha da vez (tempo) por todas as colunas (posições x)
    print(p_dirchlet_matriz)

    for i in range(len(t)):
        plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-') # vai plotar todo o vetor de x por cada um dos vetores de pressão gerados, na linha da vez (tempo) vai plotar todas as colunas (vetor de posições x) por todas as linhas (vetor de pressões gerado)
    plt.xlabel('x [m]')
    plt.ylabel('Pressão [kgf/cm²]')
    plt.grid(True)
    plt.show()

    return x, t, p_dirchlet_matriz