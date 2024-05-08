import numpy as np

def TDMA(A, D):
    a = np.diagonal(A, offset = -1)
    b = np.diagonal(A, offset = 0)
    c = np.diagonal(A, offset = +1)
    d = D

    n = len(d)
    c_ = np.zeros(n-1)
    d_ = np.zeros(n)
    x = np.zeros(n)

    # Forward elimination
    c_[0] = c[0]/b[0]
    d_[0] = d[0]/b[0]
    for i in range(1, n-1):
        c_[i] = c[i] / (b[i] - a[i-1] * c_[i-1])
    for i in range(1, n):
        d_[i] = (d[i] - a[i-1] * d[i-1]) / (b[i] - a[i-1] * c_[i-1])

    # Backward substitution
    x[-1] = d_[-1]
    for i in range(n-2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i+1]

    return x