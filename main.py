from sympy import *
from math import *
import cmath
import matplotlib.pyplot as plt
import numpy as np

def DFT(x_n, char):
    N = len(x_n)
    omega = float(2*pi/N)
    X_terms = []
    complex_term = 0

    for k in range(N):
        sum_row = 0
        print(f'{char}[{k}] = ')
        for n in range(N):
            complex_term = complex(cos(k*n*omega), -sin(k*n*omega))
            sum_row += complex_term*x_n[n]

            polar_term = cmath.polar(complex_term)
            if n == N-1:
                print(f'{x_n[n]}*{polar_term[0]}e^{polar_term[1]}j = {round_complex(sum_row,3)}')
            else:
                print(f'{x_n[n]}*{polar_term[0]}e^{polar_term[1]}j + ')
        print('\n')
        X_terms.append(round_complex(sum_row, 3))

    print(f'{char}[k] = {X_terms}')
    return X_terms

def inverse_DFT(X_n, char):
    N = len(X_n)
    omega = float(2*pi/N)
    x_terms = []
    complex_term = 0

    for k in range(N):
        sum_row = 0
        print(f'{char}[{k}] = {1/N}*(')
        for n in range(N):
            complex_term = complex(cos(k * n * omega), sin(k * n * omega))
            sum_row += (1/N) * complex_term * X_n[n]
            polar_term = cmath.polar(complex_term)
            if n == N - 1:
                print(f'{X_n[n]}*{polar_term[0]}e^{polar_term[1]}j) =  {round_complex(sum_row,3)}')
            else:
                print(f'{X_n[n]}*{polar_term[0]}e^{polar_term[1]}j + ')
        print('\n')
        x_terms.append(round_complex(sum_row, 3))

    print(f'{char}[n] = {x_terms}')
    return x_terms

def round_complex(x,n):
    return (round(x.real, n) + round(x.imag,n)*1j)

def padding(h, y):
    if len(h) != len(y):
        for i in range(len(y)-len(h)):
            h.append(0)
    return h

def deconvolution(h_n, y_n, char_h, char_y, char_x):
    X = []
    X_duration = len(y_n)-len(h_n)+1

    h_padded = padding(h_n,y_n)
    H = DFT(h_padded, char_h)
    Y = DFT(y_n, char_y)

    for i in range(len(y_n)):
        X.append(round_complex(Y[i]/H[i],3))

    return unpadding(inverse_DFT(X, char_x.lower()), X_duration)

def unpadding(x, duration):
    i = len(x)
    boolean = False

    while i-1 > duration:
        if x[i-1] == (0+0j):
            boolean = True
        i-=1

    return x[:duration] if (boolean == True) else x

def plotter(list):
    magnitude = [abs(i) for i in list]
    x = np.arange(len(magnitude))
    plt.stem(x, magnitude)
    plt.xticks(np.arange(min(x), max(x)+1,1.0))
    plt.xlabel('Index')
    plt.ylabel('Magnitude')
    plt.show()

def unit_step(x):
    if x >= 0:
        return 1
    else:
        return 0

def function(x):
    return (4*unit_step(x))-((0.25**x)*unit_step(x))

def pole_locator(H_denom, var):
    return solve(factor(H_denom, extension=[I]),var)

def zero_locator(H_num, var):
    return solve(factor(H_num, extension=[I]),var)