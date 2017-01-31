
from scipy.optimize import newton
import numpy as np
import scipy as sc
import slab_functions as sf


def solver_forwards(func, W_init, W_fin, ntries):
    W = W_init
    for i in range(ntries):
        try:
            W = newton(func, W, tol=1e-5, maxiter=50)
            if W < W_init or W > W_fin:
                W = W_init + i * (W_fin - W_init)/ntries
            else:
                break
        except RuntimeError:
            W = W_init + i * (W_fin - W_init)/ntries
        if i == ntries-1:
            print('Did not converge after ' + str(ntries) + ' tries')
    return np.real(W)
    
def solver_backwards(func, W_init, W_fin, ntries):
    W = W_init
    for i in range(ntries):
        try:
            W = newton(func, W, tol=1e-5, maxiter=50)
            if W > W_init or W < W_fin:
                W = W_init + i * (W_fin - W_init)/ntries
            else:
                break
        except RuntimeError:
            W = W_init + i * (W_fin - W_init)/ntries
        if i == ntries-1:
            print('Did not converge after ' + str(ntries) + ' tries')
    return np.real(W)


def what_mode_k(W, K, R1):
    vx_b1 = sf.constA_kink(W, K, R1)*(sc.cosh(sf.m1(W, R1)*(-2*K)) + 
                                      sc.sinh(sf.m1(W, R1)*(-2*K)))
    vx_b2 = sf.constD_kink(W, K, R1)*(sc.cosh(sf.m2(W)*(2*K)) - 
                                      sc.sinh(sf.m2(W)*(2*K)))
    if vx_b1 == 0 or vx_b2 == 0:
        return ['ERROR: divide by zero in what_mode_k','','']
    else:
        if vx_b2 / vx_b1 > 0:
            return ['KINK','asym','k']
        else:
            return ['SAUSAGE','asym','k']

def what_mode_s(W, K, R1):
    vx_b1 = sf.constA_saus(W, K, R1)*(sc.cosh(sf.m1(W, R1)*(-2*K)) + 
                                      sc.sinh(sf.m1(W, R1)*(-2*K)))
    vx_b2 = sf.constD_saus(W, K, R1)*(sc.cosh(sf.m2(W)*(2*K)) - 
                                      sc.sinh(sf.m2(W)*(2*K)))
    if vx_b1 == 0 or vx_b2 == 0:
        return ['ERROR: divide by zero in what_mode_s','','']
    else:
        if vx_b2 / vx_b1 > 0:
            return ['KINK','asym','s']
        else:
            return ['SAUSAGE','asym','s']
            
            
            
def what_mode_sym_k(W, K):
    vx_b1 = sf.constA_kink_sym(W, K)*(sc.cosh(sf.me(W)*(-2*K)) + 
                                      sc.sinh(sf.me(W)*(-2*K)))
    vx_b2 = sf.constD_kink_sym(W, K)*(sc.cosh(sf.me(W)*(2*K)) - 
                                      sc.sinh(sf.me(W)*(2*K)))
    if vx_b1 == 0 or vx_b2 == 0:
        return ['ERROR: divide by zero in what_mode_sym_k','','']
    else:
        if vx_b2 / vx_b1 > 0:
            return ['KINK','sym','k']
        else:
            return ['SAUSAGE','sym','k']

def what_mode_sym_s(W, K):
    vx_b1 = sf.constA_saus_sym(W, K)*(sc.cosh(sf.me(W)*(-2*K)) + 
                                      sc.sinh(sf.me(W)*(-2*K)))
    vx_b2 = sf.constD_saus_sym(W, K)*(sc.cosh(sf.me(W)*(2*K)) - 
                                      sc.sinh(sf.me(W)*(2*K)))
    if vx_b1 == 0 or vx_b2 == 0:
        return ['ERROR: divide by zero in what_mode_sym_s','','']
    else:
        if vx_b2 / vx_b1 > 0:
            return ['KINK','sym','s']
        else:
            return ['SAUSAGE','sym','s']
            