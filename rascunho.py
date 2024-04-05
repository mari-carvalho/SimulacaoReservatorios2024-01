            if entrada == qw and saida == p0:
            if i == 0:
                p[i,j] = qw
            elif i == 1: # bloco 1
                p[i,j] = eta*rx*p[i-1,j+1] + (1-eta*rx)*p[i-1,j] + eta*rx*qw*((mi*delta_x)/(k*A))
            elif i == len(x): #bloco N
                p[i,j] = (4/3)*eta*rx*p[i-1,j-1] + (1-4*eta*rx)*p[i-1,j] + (8/3)*eta*rx*p0
            elif i == len(x) + 1: # N+1
                p[i,j] = p0
            else: # blocos interiores
                p[i,j] = eta*rx*p[i-1,j-1] + (1-2*eta*rx)*p[i-1,j] + eta*rx*p[i-1,j+1]
        if entrada == qw and saida == q0:
            if i == 0:
                p[i,j] = qw
            elif i == 1: # bloco 1
                p[i,j] = eta*rx*p[i-1,j+1] + (1-eta*rx)*p[i-1,j] + eta*rx*qw*((mi*delta_x)/(k*A))
            elif i == len(x): #bloco N
                p[i,j] = -eta*rx*p[i-1,j-1] + (1+eta*rx)*p[i-1,j] - eta*rx*q0*((mi*delta_x)/(k*A))
            elif i == len(x) + 1: # N+1
                p[i,j] = p0
            else: # blocos interiores
                p[i,j] = eta*rx*p[i-1,j-1] + (1-2*eta*rx)*p[i-1,j] + eta*rx*p[i-1,j+1]