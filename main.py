import numpy as np
import matplotlib.pyplot as plt
from solvers import solvers
from FTCS import FTCS
from BTCS import BTCS
from CN import CN
from erro_metricas_L2 import erros
from solucao_solvers import tempo_computacional
from estabilidade import estabilidade

# Definindo as Variáveis de Entrada:

#calc_FTCS = FTCS.calculate_FTCS_pp()
#calc_erros_tempo = erros.calculate_erros_tempo()
calc_erros_malha = erros.calculate_erros_malha()
#calc_tempo_comp = tempo_computacional.calculate_tempo_computacional()
#calc_estabilidade = estabilidade.calculate_estabilidade()