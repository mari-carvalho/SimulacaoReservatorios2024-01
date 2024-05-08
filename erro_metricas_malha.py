import numpy as np
import matplotlib.pyplot as plt
from solver_gauss_seidel import gauss_seidel
from analitica_dirchlet import calculate_analitica_dirchlet
from FTCS import FTCS
from BTCS import BTCS
from CN import CN

h_x = [0.5, 0.6, 0.7]
h_t = 0.2

def calculate_h_t_ex():
    x_ex_array = []
    t_ex_array = []
    p_ex_array = []
    for i in h_x:
        x_ex, t_ex, p_ex = calculate_analitica_dirchlet(h_t, i)

        x_ex_array.append(x_ex)
        p_ex_array.append(p_ex)
        t_ex_array.append(t_ex)

    return x_ex_array, t_ex_array, p_ex_array


x_ex_array, t_ex_array, p_ex_array = calculate_h_t_ex()


def calculate_h_t_calc():
    x_calc_array = []
    t_calc_array = []
    p_calc_array = []
    for i in h_x:
        x_calc, t_calc, p_calc = BTCS.calculate_BTCS_pp(h_t, i)

        x_calc_array.append(x_calc)
        p_calc_array.append(p_calc)
        t_calc_array.append(t_calc)

    return x_calc_array, t_calc_array, p_calc_array


x_calc_array, t_calc_array, p_calc_array = calculate_h_t_calc()

# Cálculo do Erro:

'''
# Norma L2
L2_list = []
for i in range(len(p_calc_array)):  # acesso a matriz menor
    p_calc = p_calc_array[i]
    p_ex = p_ex_array[i]
    sum = 0
    y_calc = p_calc[6]
    y_ex = p_ex[6]
    for k in range(len(y_ex)):  # acesso a cada linha
        sum = sum + ((abs((y_ex[k] - y_calc[k]) / (y_ex[k]))) ** 2)
    L2 = np.sqrt((1 / (n_x ** 2)) * sum)
    L2_list.append(L2)
print('L2', L2_list)

L2_log_list = np.log(L2_list)
for i in range(len(L2_list)):
    L2_log = np.log(L2_list[i])
    L2_log_list.append(L2_log)
'''

# Norma E_inf
E_inf_depois_list = []

for i in range(len(p_calc_array)):  # acesso a matriz menor
    p_calc = p_calc_array[i]
    p_ex = p_ex_array[i]
    y_calc = p_calc[6]
    y_ex = p_ex[6]
    E_inf_antes_list = []
    for k in range(len(y_calc)):
        E_inf_antes = abs((y_ex[k] - y_calc[k]))
        E_inf_antes_list.append(E_inf_antes)
    E_inf_depois = max(E_inf_antes_list)
    E_inf_depois_list.append(E_inf_depois)

print(E_inf_depois_list)

E_inf_depois_log_list = np.log(E_inf_depois_list)

'''
for i in range(len(E_inf_depois_list)):
    E_inf_depois_log = np.log(E_inf_depois_list[i])
    E_inf_depois_log_list.append(E_inf_depois_log)

'''

'''
# Norma E_abs

err_abs_total_list = []
# test

for i in range(len(p_calc_array)):  # acesso a matriz menor
    p_calc = p_calc_array[i]
    p_ex = p_ex_array[i]
    y_calc = p_calc[6]
    y_ex = p_ex[6]
    err_abs_list = []
    sum = 0
    for k in range(len(y_ex)):
        err_abs = abs(y_ex[k] - y_calc[k])
        err_abs_list.append(err_abs)
    for j in range(len(err_abs_list)):
        sum = sum + err_abs_list[j]
    err_abs_total = 1 / n_x * sum
    err_abs_total_list.append(err_abs_total)

err_abs_total_log_list = np.log(err_abs_total_list)
for i in range(len(err_abs_total_list)):
    err_abs_total_log = np.log(err_abs_total_list[i])
    err_abs_total_log_list.append(err_abs_total_log)

'''
# Log de h_t

h_x_log_list = np.log(h_x)

'''
for i in range(len(h_t)):
    h_t_log = np.log(h_t[i])
    h_t_log_list.append(h_t_log)
'''
print('h_t_log', h_x_log_list)

# Plotagem:
'''
plt.plot(h_t_log_list, L2_log_list, linestyle='-', label='Erro Analítica/Explícita')

plt.title('Norma L2')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('L2')
plt.show()

'''

# Plotagem:

plt.plot(h_x_log_list, E_inf_depois_log_list, linestyle='-', label='Erro Analítica/Explícita')

plt.title('Norma E$ \infty$')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('L2')
plt.show()

# Plotagem:
'''
plt.plot(h_t_log_list, err_abs_total_log_list, linestyle='-', label='Erro Analítica/Explícita')

plt.title('Norma Erro Absoluto')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('L2')
plt.show()
'''
