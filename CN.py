from scipy.linalg import solve

class CN():

    def calculate_CN_pp():

        import numpy as np
        import matplotlib.pyplot as plt 
        import sympy as sp 
        from solver_gauss_seidel import gauss_seidel

        n_x = (xf-x0)/(h_x)
        n_t = (tf-t0)/(h_t)
        print('n_x',n_x)
        print('n_t',n_t)

        x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
        t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos

        # Alimentando os vetores:
        for i in range(len(x)):
            if i == 0:
                x[i] = x0
            elif i == 1:
                x[i] = i*(h_x/2)
            elif i == len(x):
                x[i] = L 
            elif i == len(x)-1:
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

        def calculate_eta(k:float, phi:float, mi:float, c:float) -> float:
            eta = k/(phi*mi*c)

            return eta 

        eta = calculate_eta(k,phi,mi,c)

        def calculate_rx(h_t:float, h_x:float) -> float:
            rx = (h_t)/(h_x**2)

            return rx 

        rx = calculate_rx(h_t, h_x)

        # Criando o método MDF de BTCS:

        # Critérios de Scarvorought, 1966:
        n = 12 
        Eppara = 0.5*10**-n

        # Número Máximo de Interações:
        maxit = 1000

        ai = -1/2*eta*rx 
        bi = 1 + rx*eta
        an = -(2/3)*rx*eta
        b1 = 1 + 2*rx*eta

        p_coeficientes = np.zeros((int(n_x)-1, int(n_x)-1)) # a matriz de coeficientes deve ser quadrada
        p_old = np.ones(int(n_x)-1)*p0 # vai atualizar cada linha da matriz 
        p_solucoes = np.zeros((int(n_t)+2, int(n_x)+1)) # matriz de soluções pode seguir a mesma lógica da FTCS, não quadrada. a matriz de soluções deve ter dimensão de 402 em linhas, porque o vetor de h começa em 0 + 1 = 1 (linha 0 da matriz de soluções é a condição inicial p0), para ir até 401 deve ter uma dimensão a mais (21 linhas, preenche até 20 começando de 0; 401 linhas, preenche até 400 começando de 0; 402 linhas, preenche até 400, começando de 1 ) 
        d = np.zeros(int(n_x)-1) # vai guardar os valores de p no tempo anterior mais 8/3*eta*rx*(Pw ou P0)
        h = 0 # para acompanhar o tamanho do vetor de tempo (0 a 9, 10 elementos), p_soluções deve ter uma posição a frente (1 a 9, 9 elementos)
        p_solucoes[h, :]  = p0 # primeira linhas todas as colunas, tempo = 0 

        for j in range(len(t)): # 0 a 4 (tamanho de t), tempo 4 elemento 5; 1 a 4 (tamanho de t, mesmo for), tempo 4 elemento 4; precisa de mais um elemento
            h = h + 1 
            p_coeficientes = np.zeros((int(n_x)-1, int(n_x)-1))
            for i in range(len(p_coeficientes)): # variando a linha
                if i == 0:
                    p_coeficientes[i,0] = b1
                    p_coeficientes[i,1] = an
                    d[i] = (1 - 2*rx*eta)*p_old[i] + (2/3*rx*eta)*p_old[i+1] + 8/3*rx*eta*pw
                elif i == len(p_coeficientes)-1: # o último, N
                    p_coeficientes[i,len(p_coeficientes)-2] = an
                    p_coeficientes[i,len(p_coeficientes)-1] = b1
                    d[i] = (1 - 2*rx*eta)*p_old[i] + (2/3*rx*eta)*p_old[i-1] + 8/3*rx*eta*p0
                else:
                    p_coeficientes[i,i-1] = ai # linha 1, coluna 0 (i-1)
                    p_coeficientes[i,i] = bi
                    p_coeficientes[i,i+1] = ai
                    d[i] = (1 - rx*eta)*p_old[i] + (1/2*rx*eta)*p_old[i+1] + (1/2*rx*eta)*p_old[i-1]# condição central é 0

            x0 = p_old # os primeiros valores de chute inicial vão ser os valores de p calculadas no tempo anterior 
            p_new = gauss_seidel(p_coeficientes,d,x0,Eppara,maxit)
            #p_new = solve(p_coeficientes,d)
            p_old = p_new # atualiza a matriz, coloca o vetor de pressões calculado no tempo anterior (p_new) em p_old 
            p_new = np.insert(p_new, 0, pw) # inserindo colunas
            p_new = np.append(p_new, p0) # append sempre no final 
            p_solucoes[h, :] = p_new # vai guardar na matriz de solucoes todos os vetores de pressão calculados nos tempos 
        
        print(p_solucoes)
        print(p_coeficientes)
        tam1 = len(p_solucoes[0])
        tam2 = len(p_coeficientes)
        print('tam', tam1)
        print('tam2', tam2)

        # Plotagem:
        time = [0,10,20,30,40,50,60,70,80,90,100]
        for i in range(len(t)):
            if t[i] in time:
                plt.plot(x, p_solucoes[i, :], linestyle='-', label=f't = {t[i]}')

        plt.legend()
        plt.title('Formulação CN - Dirchlet')
        plt.xlabel('Comprimento (m)')
        plt.ylabel('Pressão (psia)')
        plt.grid()
        plt.show()

    def calculate_CN_pp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x):

        import numpy as np
        import matplotlib.pyplot as plt 
        import sympy as sp 
        from solver_gauss_seidel import gauss_seidel

        n_x = (xf-x0)/(h_x)
        n_t = (tf-t0)/(h_t)
        print('n_x',n_x)
        print('n_t',n_t)

        x = np.zeros(int(n_x)+1) # de 0 ao tamanho do reservatório com 10 elementos na malha 
        t = np.zeros(int(n_t)+1) # de 0 a 10 segundos com 10 elementos

        # Alimentando os vetores:
        for i in range(len(x)):
            if i == 0:
                x[i] = x0
            elif i == 1:
                x[i] = i*(h_x/2)
            elif i == len(x):
                x[i] = L 
            elif i == len(x)-1:
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

        def calculate_eta(k:float, phi:float, mi:float, c:float) -> float:
            eta = k/(phi*mi*c)

            return eta 

        eta = calculate_eta(k,phi,mi,c)

        def calculate_rx(h_t:float, h_x:float) -> float:
            rx = (h_t)/(h_x**2)

            return rx 

        rx = calculate_rx(h_t, h_x)

        # Criando o método MDF de BTCS:

        # Critérios de Scarvorought, 1966:
        n = 12 
        Eppara = 0.5*10**-n

        # Número Máximo de Interações:
        maxit = 1000

        ai = -1/2*eta*rx 
        bi = 1 + rx*eta
        an = -(2/3)*rx*eta
        b1 = 1 + 2*rx*eta

        p_coeficientes = np.zeros((int(n_x)-1, int(n_x)-1)) # a matriz de coeficientes deve ser quadrada
        p_old = np.ones(int(n_x)-1)*p0 # vai atualizar cada linha da matriz 
        p_solucoes = np.zeros((int(n_t)+2, int(n_x)+1)) # matriz de soluções pode seguir a mesma lógica da FTCS, não quadrada. a matriz de soluções deve ter dimensão de 402 em linhas, porque o vetor de h começa em 0 + 1 = 1 (linha 0 da matriz de soluções é a condição inicial p0), para ir até 401 deve ter uma dimensão a mais (21 linhas, preenche até 20 começando de 0; 401 linhas, preenche até 400 começando de 0; 402 linhas, preenche até 400, começando de 1 ) 
        d = np.zeros(int(n_x)-1) # vai guardar os valores de p no tempo anterior mais 8/3*eta*rx*(Pw ou P0)
        h = 0 # para acompanhar o tamanho do vetor de tempo (0 a 9, 10 elementos), p_soluções deve ter uma posição a frente (1 a 9, 9 elementos)
        p_solucoes[h, :]  = p0 # primeira linhas todas as colunas, tempo = 0 

        for j in range(len(t)): # 0 a 4 (tamanho de t), tempo 4 elemento 5; 1 a 4 (tamanho de t, mesmo for), tempo 4 elemento 4; precisa de mais um elemento
            h = h + 1 
            p_coeficientes = np.zeros((int(n_x)-1, int(n_x)-1))
            for i in range(len(p_coeficientes)): # variando a linha
                if i == 0:
                    p_coeficientes[i,0] = b1
                    p_coeficientes[i,1] = an
                    d[i] = (1 - 2*rx*eta)*p_old[i] + (2/3*rx*eta)*p_old[i+1] + 8/3*rx*eta*pw
                elif i == len(p_coeficientes)-1: # o último, N
                    p_coeficientes[i,len(p_coeficientes)-2] = an
                    p_coeficientes[i,len(p_coeficientes)-1] = b1
                    d[i] = (1 - 2*rx*eta)*p_old[i] + (2/3*rx*eta)*p_old[i-1] + 8/3*rx*eta*p0
                else:
                    p_coeficientes[i,i-1] = ai # linha 1, coluna 0 (i-1)
                    p_coeficientes[i,i] = bi
                    p_coeficientes[i,i+1] = ai
                    d[i] = (1 - rx*eta)*p_old[i] + (1/2*rx*eta)*p_old[i+1] + (1/2*rx*eta)*p_old[i-1]# condição central é 0

            x0 = p_old # os primeiros valores de chute inicial vão ser os valores de p calculadas no tempo anterior 
            p_new = gauss_seidel(p_coeficientes,d,x0,Eppara,maxit)
            #p_new = solve(p_coeficientes,d)
            p_old = p_new # atualiza a matriz, coloca o vetor de pressões calculado no tempo anterior (p_new) em p_old 
            p_new = np.insert(p_new, 0, pw) # inserindo colunas
            p_new = np.append(p_new, p0) # append sempre no final 
            p_solucoes[h, :] = p_new # vai guardar na matriz de solucoes todos os vetores de pressão calculados nos tempos 
        
        print(p_solucoes)
        print(p_coeficientes)
        tam1 = len(p_solucoes[0])
        tam2 = len(p_coeficientes)
        print('tam', tam1)
        print('tam2', tam2)

        # Plotagem:
        time = [0,10,20,30,40,50,60,70,80,90,100]
        for i in range(len(t)):
            if t[i] in time:
                plt.plot(x, p_solucoes[i, :], linestyle='-', label=f't = {t[i]}')

        plt.legend()
        plt.title('Formulação CN - Dirchlet')
        plt.xlabel('Comprimento (m)')
        plt.ylabel('Pressão (psia)')
        plt.grid()
        plt.show()

