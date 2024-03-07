# Fluxo Linear com Condição de Contorno Pressão-Pressão (Dirchlet)

# Importando Bibliotecas:

import numpy as np 
import math as mt 

# Criando a Função para a Solução Analítica:

def calculate_Dirchlet(pe:float, pw:float, x:float, L:float, k:float, phi:float, mi:float, c:float, t:float) -> float:

    p_dirchlet = (pe-pw) * ((x/L) + (2/mt.pi))