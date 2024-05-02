import numpy as np
import matplotlib.pyplot as plt
from solver_gauss_seidel import gauss_seidel
from analitica_dirchlet import calculate_analitica_dirchlet
from FTCS import FTCS
from BTCS import BTCS
from CN import CN 

n_x = 10
h_t = [0.1, 0.2, 0.25]

def calculate_h_t_ex():
    x_ex_array = []
    t_ex_array = []
    p_ex_array = []
    for i in h_t:

        x_ex, t_ex, p_ex = calculate_analitica_dirchlet(i)

        x_ex_array.append(x_ex)
        p_ex_array.append(p_ex)
        t_ex_array.append(t_ex)

    return x_ex_array, t_ex_array, p_ex_array

x_ex_array, t_ex_array, p_ex_array = calculate_h_t_ex()

def calculate_h_t_calc():
    x_calc_array = []
    t_calc_array = []
    p_calc_array = []
    for i in h_t:
        
        x_calc, t_calc, p_calc = FTCS.calculate_FTCS_pp(i)

        x_calc_array.append(x_calc)
        p_calc_array.append(p_calc)
        t_calc_array.append(t_calc)

    return x_calc_array, t_calc_array, p_calc_array
    

x_calc_array, t_calc_array, p_calc_array = calculate_h_t_calc()

# Cálculo do Erro:

#L2_list = []
'''
for i in range(len(p_calc_array)): # acesso a matriz menor
    p_calc = p_calc_array[i]
    p_ex = p_ex_array[i]
    sum = 0 
    y_calc = p_calc[6]
    y_ex = p_ex[6]
    for k in range(len(y_ex)): # acesso a cada linha
        sum = sum + ((y_ex[k]-y_calc[k])/(y_ex[k]))**2
    L2 = np.sqrt((1/(n_x**2))*sum)
    L2_list.append(L2)
print('L2', L2_list)
'''
E_inf_list = []
E_inf_antes_list = []
E_inf_depois_list  = []
for i in range(len(p_calc_array)): # acesso a matriz menor
    p_calc = p_calc_array[i]
    p_ex = p_ex_array[i]
    sum = 0 
    y_calc = p_calc[6]
    y_ex = p_ex[6]
    E_inf_antes = (y_ex[i] - y_calc[i])
    E_inf_antes_list.append(E_inf_antes_list)
    E_inf_depois = max(E_inf_antes_list)
    E_inf_depois_list.append(E_inf_depois)

print(E_inf_depois_list)


# Plotagem:

plt.plot(h_t, E_inf_depois_list, linestyle='-', label='Erro Analítica/Explícita')

plt.title('Norma E$n \infty$')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('L2')
plt.show()











#calc_FTCS = FTCS.calculate_FTCS_fp(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)
#calc_FTCS = FTCS.calculate_FTCS_ff(p0,pw,qw,q0,cc,mi,k,h,phi,c,L,A,x0,xf,t0,tf,h_t,h_x)