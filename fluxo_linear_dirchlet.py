# Fluxo Linear com Condição de Contorno Pressão-Pressão (Dirchlet)

# Importando as Bibliotecas:

import numpy as np 
import math as mt 
import matplotlib.pyplot as plt 

# Dados de Entrada:

p0 = 19000000
pw = 9000000
qw = 0.01
mi = 0.001
k = 9.869233e-14
h = 10 
phi = 0.2
c = 2.04e-9
L = 20
A = 30 

x = np.linspace(0,L,1000) # de 0 a 10, com 1000 elementos
t = np.linspace(0,1000,10) # de 0 a 1000, com 10 elementos

p_dirchlet_matriz = np.zeros((len(t), len(x)))

def calculate_dirchlet(p0:float, pw:float, x:float, L:float, k:float, phi:float, mi:float, c:float, t:float) -> float:

    p_dirchlet = np.zeros(len(x))

    for i in range(len(x)): # um valor de pressão é calculado para cada valor de x, resultando em um vetor no final 
        # usa-se o while porque o somatório tende a ir de n a infinito, precisaria de muitas interações
        erro = 1000 # chutamos valor alto no início para ir diminuindo progressivamente 
        eppara = 1e-3 # critério de parada 
        n = 0 # o somatório começa de 1 e vai ao infinito (até o critério de parada nesse caso), primeiro n vai ser 0 + 1, ou seja, 1 
        sum = 0 
        sum_old = 100 
        while erro >= eppara:
            n += 1 # atualiza o valor de n (somatório)
            k = n - 1 # em relação aos somatórios, para diferenciar as somas de cada interação do while 
            sum = (np.exp(-((n*np.pi/L)**2)*(k/(phi*mi*c))*t)/n * np.sin(n*np.pi*x[i]/L)) # vai usar o valor de x da interação (for) e o valor de tempo t passado na chamada da função
            erro = abs((sum[k] - sum_old)/sum[k])*100 # calcula o erro para o critério de parada 
            sum_old = sum[k] # depois de calcular o novo valor do somatório, adota esse como antigo para a próxima interação, acumulando os resultados, o valor recém calculado deve ser incluído na soma 
        p_dirchlet[i] = (p0-pw) * ((x[i]/L) + (2/np.pi) * sum) + pw # não usa [i] no sum, porque já é o resultado de somatório da última interação, último sum calculado. Vai adicionar um valor de pressão ao vetor de pressões que depois vai ser posicionado na linha de t, cada valor de pressão vai se encaixar em uma posição de x (coluna)
    return p_dirchlet # devolve o vetor de pressões que vai entrar na matriz, cada valor de pressão do vetor entra em uma posição de x, ou seja, em uma coluna

for i in range(len(t)):
    p_dirchlet_matriz[i,:] = calculate_dirchlet(p0, pw, x, L, k, phi, mi, c, t[i]) # vai preencher a matriz com um vetor de pressões que é calculado dentro da função def e que depende das posições de x, cada valor de pressão entra em uma posição de x (colunas)
print(p_dirchlet_matriz)

for i in range(len(t)):
    plt.plot(x, p_dirchlet_matriz[i, :], linestyle='-') # vai plotar todo o vetor de x por cada um dos vetores de pressão gerados, na linha da vez (tempo) vai plotar todas as colunas (vetor de posições x) por todas as linhas (vetor de pressões gerado)
plt.xlabel('x [m]')
plt.ylabel('Pressão [kgf/cm²]')
plt.grid(True)
plt.show()