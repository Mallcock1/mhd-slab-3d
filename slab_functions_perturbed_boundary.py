
import numpy as np
import scipy as sc
from scipy.optimize import fsolve

# SBB
# Define the sound speeds and alfven speeds.
c2 = 1.2
c0 = 1.
vA = 0.9
cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))
ce = 1.2 #0.7


## SBS
## Define the sound speeds and alfven speeds.
#c2 = 1.2
#c0 = 1.
#vA = 1.3
#cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))
#ce = 1.2

R2 = 2.
Re = 2.

def c1(R1):
    return c2 * sc.sqrt(R2 / R1)

B0 = 10.

# variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

def m0(W):
    return sc.sqrt((c0**2 - W**2)*(vA**2 - W**2) / 
                   ((c0**2 + vA**2)*(cT**2 - W**2)))
                   
def m00(W):
    return sc.sqrt(1 - W**2/c0**2)
    
def m1(W, R1):
    return sc.sqrt(1 - W**2*(c2*sc.sqrt(R2/R1))**(-2))
    
def m2(W):
    return sc.sqrt(1 - W**2/c2**2)

def me(W):
    return sc.sqrt(1 - W**2*ce**(-2))

#def const(W):
#    if m0(W)**2 < 0:
#        constant = 0.3 # for bosy modes
#    else:
#        constant = 0.05 # for surface modes
#    return constant

const = 0.05 #for surface
#const = 0.3 # Body

# MAKE THESE GUYS FUCTIONS!
constC_s = const
constB_kink = const
constC_s_sym = const
constB_kink_sym = const
    
    
def disp_rel_asym(W, K, R1):
    return ((W**4*m0(W)**2*R1*R2 + (vA**2 - W**2)**2*m1(W, R1)*m2(W) -
            0.5*m0(W)*W**2*(vA**2 - W**2)*(R2*m1(W, R1) + R1*m2(W))*
            (sc.tanh(m0(W)*K) + (sc.tanh(m0(W)*K))**(-1))) /
            (vA**2 - W**2)*(c0**2 - W**2)*(cT**2 - W**2))

def disp_rel_saus_sym(W, K):
    return (sc.sqrt(ce**2-W**2) -
            ce*Re*W**2*m0(W)*sc.tanh(m0(W)*K) / (vA**2-W**2))
            
def disp_rel_kink_sym(W, K):
    return (sc.sqrt(ce**2-W**2) -
            ce*Re*W**2*m0(W) / ((vA**2-W**2)*sc.tanh(m0(W)*K)))
    

def lamb0(W):
    return -(vA**2-W**2)*1.j/(m0(W)*W)

def lamb00(W):
    return W*1.j/m00(W)

def lamb1(W, R1):
    return R1*W*1.j/m1(W, R1)
    
def lamb2(W):
    return R2*W*1.j/m2(W)

def lambe(W):
    return Re*W*1.j/me(W)
    
def bz_eq_1d(x, z, t, W, K, R1):
    truth = np.array((x <= (K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                 (x >= (-K + xi_boundary_kink(z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bz_eqfunction = np.zeros(len(x))
    for i in indices:
        bz_eqfunction[i] = B0
    return bz_eqfunction
    
def bz_eq(x, z, t, W, K, R1):
    bz_eq_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_eq_array[:,i] = bz_eq_1d(x, z[i], t, W, K, R1)
    return bz_eq_array


def constC_saus(W):
    if m0(W)**2 < 0:
        return constC_s*1j
    else:
        return constC_s
        
def constC_saus_sym(W):
    if m0(W)**2 < 0:
        return constC_s_sym*1j
    else:
        return constC_s_sym

##############################################
# symmetric sausage solution

constB_saus_sym = 0

def constA_saus_sym(W, K):
    return (-constC_saus_sym(W)*sc.sinh(m0(W)*K) /
            (sc.cosh(me(W)*K) - sc.sinh(me(W)*K)))

def constD_saus_sym(W, K):
    return ((constC_saus_sym(W)*sc.sinh(m0(W)*K)) /
            (sc.cosh(me(W)*K) - sc.sinh(me(W)*K)))

def vxhat_saus_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhatfunction[i] = (constB_saus_sym*sc.cosh(m0(W)*x[i]) + 
                            constC_saus_sym(W)*sc.sinh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhatfunction[i] = constA_saus_sym(W, K)*(sc.cosh(me(W)*x[i]) + 
                                                  sc.sinh(me(W)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhatfunction[i] = constD_saus_sym(W, K)*(sc.cosh(me(W)*x[i]) - 
                                                  sc.sinh(me(W)*x[i]))
    return vxhatfunction
    
def vxhat_dash_saus_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhat_dashfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhat_dashfunction[i] = m0(W)*(constB_saus_sym*sc.sinh(m0(W)*x[i]) + 
                                       constC_saus_sym(W)*sc.cosh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhat_dashfunction[i] = me(W)*constA_saus_sym(W, K)*(sc.sinh(me(W)*x[i]) + 
                                                                 sc.cosh(me(W)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhat_dashfunction[i] = me(W)*constD_saus_sym(W, K)*(sc.sinh(me(W)*x[i]) - 
                                                             sc.cosh(me(W)*x[i]))
    return vxhat_dashfunction
    
def p_tothat_saus_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    p_tothatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        p_tothatfunction[i] = (constB_saus_sym*sc.sinh(m0(W)*x[i]) + 
                                     constC_saus_sym(W)*sc.cosh(m0(W)*x[i])) * lamb0(W)
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        p_tothatfunction[i] = constA_saus_sym(W, K)*(sc.sinh(me(W)*x[i]) + 
                                                           sc.cosh(me(W)*x[i])) * lambe(W)
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        p_tothatfunction[i] = constD_saus_sym(W, K)*(sc.sinh(me(W)*x[i]) - 
                                                           sc.cosh(me(W)*x[i])) * lambe(W)
    return p_tothatfunction

def vx_saus_sym(x, z, t, W, K):
    return np.outer(vxhat_saus_sym(x, W, K), np.exp(1j*(z-t)))

def vx_dash_saus_sym(x, z, t, W, K):
    return np.outer(vxhat_dash_saus_sym(x, W, K), np.exp(1j*(z-t)))
    
def p_tot_saus_sym(x, z, t, W, K):
    return np.outer(p_tothat_saus_sym(x, W, K), np.exp(1j*(z-t)))


##############################################
# symmetric kink solution

constC_kink_sym = 0

def constA_kink_sym(W, K):
    return (constB_kink_sym*sc.cosh(m0(W)*K) /
            (sc.cosh(me(W)*K) - sc.sinh(me(W)*K)))

def constD_kink_sym(W, K):
    return ((constB_kink_sym*sc.cosh(m0(W)*K)) /
            (sc.cosh(me(W)*K) - sc.sinh(me(W)*K)))

def vxhat_kink_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhatfunction[i] = (constB_kink_sym*sc.cosh(m0(W)*x[i]) + 
                            constC_kink_sym*sc.sinh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhatfunction[i] = constA_kink_sym(W, K)*(sc.cosh(me(W)*x[i]) + 
                                                  sc.sinh(me(W)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhatfunction[i] = constD_kink_sym(W, K)*(sc.cosh(me(W)*x[i]) - 
                                                  sc.sinh(me(W)*x[i]))
    return vxhatfunction
    
def vxhat_dash_kink_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhat_dashfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhat_dashfunction[i] = m0(W)*(constB_kink_sym*sc.sinh(m0(W)*x[i]) + 
                                       constC_kink_sym*sc.cosh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhat_dashfunction[i] = me(W)*constA_kink_sym(W, K)*(sc.sinh(me(W)*x[i]) + 
                                                                 sc.cosh(me(W)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhat_dashfunction[i] = me(W)*constD_kink_sym(W, K)*(sc.sinh(me(W)*x[i]) - 
                                                             sc.cosh(me(W)*x[i]))
    return vxhat_dashfunction
    
def p_tothat_kink_sym(x, W, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    p_tothatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        p_tothatfunction[i] = (constB_kink_sym*sc.sinh(m0(W)*x[i]) + 
                                     constC_kink_sym*sc.cosh(m0(W)*x[i])) * lamb0(W)
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        p_tothatfunction[i] = constA_kink_sym(W, K)*(sc.sinh(me(W)*x[i]) + 
                                                           sc.cosh(me(W)*x[i])) * lambe(W)
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        p_tothatfunction[i] = constD_kink_sym(W, K)*(sc.sinh(me(W)*x[i]) - 
                                                           sc.cosh(me(W)*x[i])) * lambe(W)
    return p_tothatfunction

def vx_kink_sym(x, z, t, W, K):
    return np.outer(vxhat_kink_sym(x, W, K), np.exp(1j*(z-t)))

def vx_dash_kink_sym(x, z, t, W, K):
    return np.outer(vxhat_dash_kink_sym(x, W, K), np.exp(1j*(z-t)))
    
def p_tot_kink_sym(x, z, t, W, K):
    return np.outer(p_tothat_kink_sym(x, W, K), np.exp(1j*(z-t)))


##############################################
# sausage solution

def constB_saus(W, K, R1):
    return constC_saus(W)*((lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)) /
                 (lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)))

def constA_saus(W, K, R1):
    return ((constB_saus(W, K, R1)*sc.cosh(m0(W)*K) - constC_saus(W)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m1(W, R1)*K) - sc.sinh(m1(W, R1)*K)))

def constD_saus(W, K, R1):
    return ((constB_saus(W, K, R1)*sc.cosh(m0(W)*K) + constC_saus(W)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m2(W)*K) - sc.sinh(m2(W)*K)))

def vxhat_saus(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhatfunction[i] = (constB_saus(W, K, R1)*sc.cosh(m0(W)*x[i]) + 
                           constC_saus(W)*sc.sinh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhatfunction[i] = constA_saus(W, K, R1)*(sc.cosh(m1(W, R1)*x[i]) + 
                                                  sc.sinh(m1(W, R1)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhatfunction[i] = constD_saus(W, K, R1)*(sc.cosh(m2(W)*x[i]) - 
                                                  sc.sinh(m2(W)*x[i]))
    return vxhatfunction
    
def vxhat_dash_saus(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhat_dashfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhat_dashfunction[i] = m0(W)*(constB_saus(W, K, R1)*sc.sinh(m0(W)*x[i]) + 
                                       constC_saus(W)*sc.cosh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhat_dashfunction[i] = m1(W, R1)*constA_saus(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                                 sc.cosh(m1(W, R1)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhat_dashfunction[i] = m2(W)*constD_saus(W, K, R1)*(sc.sinh(m2(W)*x[i]) - 
                                                             sc.cosh(m2(W)*x[i]))
    return vxhat_dashfunction
    
def p_tothat_saus(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    p_tothatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        p_tothatfunction[i] = (constB_saus(W, K, R1)*sc.sinh(m0(W)*x[i]) + 
                                     constC_saus(W)*sc.cosh(m0(W)*x[i])) * lamb0(W)
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        p_tothatfunction[i] = constA_saus(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                               sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1)
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        p_tothatfunction[i] = constD_saus(W, K, R1)*(sc.sinh(m2(W)*x[i]) - 
                                                           sc.cosh(m2(W)*x[i])) * lamb2(W)
    return p_tothatfunction

def vx_saus(x, z, t, W, K, R1):
    return np.outer(vxhat_saus(x, W, K, R1), np.exp(1j*(z-t)))

def vx_dash_saus(x, z, t, W, K, R1):
    return np.outer(vxhat_dash_saus(x, W, K, R1), np.exp(1j*(z-t)))
    
def p_tot_saus(x, z, t, W, K, R1):
    return np.outer(p_tothat_saus(x, W, K, R1), np.exp(1j*(z-t)))
    
    
##############################################
# kink solution

def constC_kink(W, K, R1):
    return constB_kink*((lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)) /
                 (lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)))

def constA_kink(W, K, R1):
    return ((constB_kink*sc.cosh(m0(W)*K) - constC_kink(W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m1(W, R1)*K) - sc.sinh(m1(W, R1)*K)))

def constD_kink(W, K, R1):
    return ((constB_kink*sc.cosh(m0(W)*K) + constC_kink(W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m2(W)*K) - sc.sinh(m2(W)*K)))

def vxhat_kink(x, W, K, R1):
    if type(x) == np.float64:
        if np.abs(x) <= K:
            vxhatfunction = (constB_kink*sc.cosh(m0(W)*x) + 
                            constC_kink(W, K, R1)*sc.sinh(m0(W)*x))
        elif x < -K:
            vxhatfunction = constA_kink(W, K, R1)*(sc.cosh(m1(W, R1)*x) + 
                                                  sc.sinh(m1(W, R1)*x))
        elif x > K:
            vxhatfunction = constD_kink(W, K, R1)*(sc.cosh(m2(W)*x) - 
                                                  sc.sinh(m2(W)*x))
    else:
        truth = np.array(np.abs(x) <= K*np.ones(len(x)))
        indices = np.where(truth == True)
        vxhatfunction = np.zeros(len(x), dtype=complex)
        for i in indices:
            vxhatfunction[i] = (constB_kink*sc.cosh(m0(W)*x[i]) + 
                                constC_kink(W, K, R1)*sc.sinh(m0(W)*x[i]))
        truth2 = np.array(x < -K*np.ones(len(x)))
        indices2 = np.where(truth2 == True)
        for i in indices2:
            vxhatfunction[i] = constA_kink(W, K, R1)*(sc.cosh(m1(W, R1)*x[i]) + 
                                                      sc.sinh(m1(W, R1)*x[i]))
        truth3 = np.array(x > K*np.ones(len(x)))
        indices3 = np.where(truth3 == True)
        for i in indices3:
            vxhatfunction[i] = constD_kink(W, K, R1)*(sc.cosh(m2(W)*x[i]) - 
                                                      sc.sinh(m2(W)*x[i]))
    return vxhatfunction
    
def vxhat_dash_kink(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhat_dashfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        vxhat_dashfunction[i] = m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
                                       constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhat_dashfunction[i] = m1(W, R1)*constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                                 sc.cosh(m1(W, R1)*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhat_dashfunction[i] = m2(W)*constD_kink(W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                                                             sc.cosh(m2(W)*x[i]))
    return vxhat_dashfunction
    
def vzhat_kink(x, W, K, R1):
    if type(x) == np.float64:
        if np.abs(x) <= K:
            vzhat_function = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB_kink*sc.sinh(m0(W)*x) +
                                       constC_kink(W, K, R1)*sc.cosh(m0(W)*x))
        elif x < -K:
            vzhat_function = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x) + 
                                                                 sc.cosh(m1(W, R1)*x))
        elif x > K:
            vzhat_function = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD_kink(W, K, R1)*(sc.sinh(m2(W)*x) -
                                                             sc.cosh(m2(W)*x))
    else:
        truth = np.array(np.abs(x) <= K*np.ones(len(x)))
        indices = np.where(truth == True)
        vzhat_function = np.zeros(len(x), dtype=complex)
        for i in indices:
            vzhat_function[i] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
                                           constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i]))
        truth2 = np.array(x < -K*np.ones(len(x)))
        indices2 = np.where(truth2 == True)
        for i in indices2:
            vzhat_function[i] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                                     sc.cosh(m1(W, R1)*x[i]))
        truth3 = np.array(x > K*np.ones(len(x)))
        indices3 = np.where(truth3 == True)
        for i in indices3:
            vzhat_function[i] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD_kink(W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                                                                 sc.cosh(m2(W)*x[i]))
    return vzhat_function
    
def p_tothat_kink(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    p_tothatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        p_tothatfunction[i] = (constB_kink*sc.sinh(m0(W)*x[i]) +
                                     constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i])) * lamb0(W)
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        p_tothatfunction[i] = constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) +
                                                               sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1)
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        p_tothatfunction[i] = constD_kink(W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                                                           sc.cosh(m2(W)*x[i])) * lamb2(W)
    return p_tothatfunction

def vx_kink(x, z, t, W, K, R1):
    if type(vxhat_kink(x, W, K, R1)) == np.complex128:
        return vxhat_kink(x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vxhat_kink(x, W, K, R1), np.exp(1j*(z-t)))

def vx_dash_kink(x, z, t, W, K, R1):
    return np.outer(vxhat_dash_kink(x, W, K, R1), np.exp(1j*(z-t)))

def vz_kink(x, z, t, W, K, R1):
    if type(vzhat_kink(x, W, K, R1)) == np.complex128:
        return vzhat_kink(x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vzhat_kink(x, W, K, R1), np.exp(1j*(z-t)))
    
def xix_kink(x, z, t, W, K, R1):
    if type(vxhat_kink(x, W, K, R1)) == np.complex128:
        return (1j * vxhat_kink(x, W, K, R1) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vxhat_kink(x, W, K, R1) / W, np.exp(1j*(z-t)))
    
def xiz_kink(x, z, t, W, K, R1):
    if type(vzhat_kink(x, W, K, R1)) == np.complex128:
        return (1j * vzhat_kink(x, W, K, R1) /W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vzhat_kink(x, W, K, R1) / W, np.exp(1j*(z-t)))
    
def p_tot_kink(x, z, t, W, K, R1):
    return np.outer(p_tothat_kink(x, W, K, R1), np.exp(1j*(z-t)))
    
    
def xihat_boundary_kink(W, K, R1, boundary='r'):
    if boundary == 'r' or boundary == 'right':
        xihat_kink = (1j / W) * (constB_kink*sc.cosh(m0(W)*K) + 
                                 constC_kink(W, K, R1)*sc.sinh(m0(W)*K))
    if boundary == 'l' or boundary == 'left':
        xihat_kink = (1j / W) * (constB_kink*sc.cosh(m0(W)*-K) + 
                                 constC_kink(W, K, R1)*sc.sinh(m0(W)*-K))   
    return xihat_kink

def xi_boundary_kink(z, t, W, K, R1, boundary='r'):           
    return xihat_boundary_kink(W, K, R1, boundary) * np.exp(1j*(z-t))
    

def bxhat_kink(x, z, t, W, K, R1):
    truth = np.array((x <= (K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                     (x >= (-K + xi_boundary_kink(z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bxhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bxhat_function[i] = (-B0/W)*(constB_kink*sc.cosh(m0(W)*x[i]) +
                                   constC_kink(W, K, R1)*sc.sinh(m0(W)*x[i]))
#    truth2 = np.array(x < (-K + xi_boundary_kink(z, t, W, R1, K, boundary='l'))*np.ones(len(x)))
#    indices2 = np.where(truth2 == True)
#    for i in indices2:
#        bxhat_function[i] = 0.
#    truth3 = np.array(x > (K + xi_boundary_kink(z, t, W, R1, K, boundary='r'))*np.ones(len(x)))
#    indices3 = np.where(truth3 == True)
#    for i in indices3:
#        bxhat_function[i] = 0.
    return bxhat_function
    
def bzhat_kink(x, z, t, W, K, R1):
    truth = np.array((x <= (K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                     (x >= (-K + xi_boundary_kink(z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bzhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bzhat_function[i] = (-1j*B0/W)*m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
                                   constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i]))
#    truth2 = np.array(x < (-K + xi_boundary_kink(z, t, W, R1, K, boundary='l'))*np.ones(len(x)))
#    indices2 = np.where(truth2 == True)
#    for i in indices2:
#        bzhat_function[i] = 0.
#    truth3 = np.array(x > (K + xi_boundary_kink(z, t, W, R1, K, boundary='r'))*np.ones(len(x)))
#    indices3 = np.where(truth3 == True)
#    for i in indices3:
#        bzhat_function[i] = 0.
    return bzhat_function
    
    
def bx_kink(x, z, t, W, K, R1):
    bx_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bx_array[:,i] = bxhat_kink(x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bx_array
    
def bz_kink(x, z, t, W, K, R1):
    bz_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_array[:,i] = bzhat_kink(x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bz_array
    
# In the slab has an extra -1 factor to make it work. Not sure why it needs this.
def rho_hat_kink(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    rho_hatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        rho_hatfunction[i] = m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
                             constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i])) * lamb00(W) / (c0**2 * m00(W))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        rho_hatfunction[i] = constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) +
                             sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1) / c1(R1)**2
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        rho_hatfunction[i] = constD_kink(W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                             sc.cosh(m2(W)*x[i])) * lamb2(W) / c2**2
    return rho_hatfunction
    
def rho_kink(x, z, t, W, K, R1):
    rho_kinkfunc = np.outer(rho_hat_kink(x, W, K, R1), np.exp(1j*(z-t)))
    return rho_kinkfunc
    
    
# A more general function, that can find the lagrangian density.
def rho_kink_lagrang(x, z, t, W, K, R1):
    rho_kink_xz = np.real(np.outer(rho_hat_kink(x, W, K, R1), np.exp(1j*(z-t))))
    NUG = np.real(np.array(np.meshgrid(x,z)) + np.array([xix_kink(x, z, t, W, K, R1), 
                                                 xiz_kink(x, z, t, W, K, R1)]))
    f = sc.interpolate.interp2d(NUG[0,:,:], NUG[1,:,:], rho_kink_xz, kind='cubic')
    return f(x,z)
                
                
# for plotting rho(x + xi,t) at (x,t)
#def rho_hat_kink_pert(x, z, t, W, K, R1):
#    truth = np.array((x <= (K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
#                     (x >= (-K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x))))
#    indices = np.where(truth == True)
#    rho_hatfunction = np.zeros(len(x), dtype=complex)
#    for i in indices:
#        rho_hatfunction[i] = -m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
#                             constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i])) * lamb00(W) / (c0**2 * m00(W))
#    truth2 = np.array(x < (-K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x)))
#    indices2 = np.where(truth2 == True)
#    for i in indices2:
#        rho_hatfunction[i] = constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) +
#                             sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1) / c1(R1)**2
#    truth3 = np.array(x > (K + xi_boundary_kink(z, t, W, K, R1, boundary='r'))*np.ones(len(x)))
#    indices3 = np.where(truth3 == True)
#    for i in indices3:
#        rho_hatfunction[i] = constD_kink(W, K, R1)*(sc.sinh(m2(W)*x[i]) -
#                             sc.cosh(m2(W)*x[i])) * lamb2(W) / c2**2
#    return rho_hatfunction
#
#def rho_kink_pert(x, z, t, W, K, R1):
#    rho_kink_array = np.zeros((len(x), len(z)), dtype=complex)
#    for i in range(len(z)):
#        rho_kink_array[:,i] = rho_hat_kink_pert(x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
#    return rho_kink_array
#    
def rho_kink_pert(x, z, t, W, K, R1):
    rho_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    rho_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            pert = [x[i] + xix_kink(x[i], z[j], t, W, K, R1),
                    z[j] + xiz_kink(x[i], z[j], t, W, K, R1)]
            if abs(x[i]) <= (K + xix_kink(x[i], z[j], t, W, K, R1)):
                rho_hat_vals[i,j] = -m0(W)*(constB_kink*sc.sinh(m0(W)*pert[0]) +
                                constC_kink(W, K, R1)*sc.cosh(m0(W)*pert[0])) * lamb00(W) / (c0**2 * m00(W))
            elif x[i] < (-K + xix_kink(x[i], z[j], t, W, K, R1)):
                rho_hat_vals[i,j] = constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*pert[0]) +
                                sc.cosh(m1(W, R1)*pert[0])) * lamb1(W, R1) / c1(R1)**2
            elif x[i] > (K + xix_kink(x[i], z[j], t, W, K, R1)):
                rho_hat_vals[i,j] = constD_kink(W, K, R1)*(sc.sinh(m2(W)*pert[0]) -
                                     sc.cosh(m2(W)*pert[0])) * lamb2(W) / c2**2
            
            rho_vals[i,j] = rho_hat_vals[i,j] * np.exp(1j*(pert[1]-t))
    return rho_vals
    

# for plotting rho(x,t) at (x + xi,t)
def rho_kink_pert2(x, z, t, W, K, R1):
    rho_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    rho_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix_kink(r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz_kink(r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                rho_hat_vals[i,j] = -m0(W)*(constB_kink*sc.sinh(m0(W)*sol[0]) +
                                constC_kink(W, K, R1)*sc.cosh(m0(W)*sol[0])) * lamb00(W) / (c0**2 * m00(W))
            elif sol[0] < -K:
                rho_hat_vals[i,j] = constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) +
                                sc.cosh(m1(W, R1)*sol[0])) * lamb1(W, R1) / c1(R1)**2
            elif sol[0] > K:
                rho_hat_vals[i,j] = constD_kink(W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                     sc.cosh(m2(W)*sol[0])) * lamb2(W) / c2**2
            
            rho_vals[i,j] = rho_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return rho_vals
    
def vx_kink_pert(x, z, t, W, K, R1):
    vx_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vx_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix_kink(r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz_kink(r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vx_hat_vals[i,j] = (constB_kink*sc.cosh(m0(W)*sol[0]) + 
                            constC_kink(W, K, R1)*sc.sinh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vx_hat_vals[i,j] = constA_kink(W, K, R1)*(sc.cosh(m1(W, R1)*sol[0]) + 
                                                  sc.sinh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vx_hat_vals[i,j] = constD_kink(W, K, R1)*(sc.cosh(m2(W)*sol[0]) - 
                                                  sc.sinh(m2(W)*sol[0]))
            
            vx_vals[i,j] = vx_hat_vals[i,j] * np.exp(1j*(z[j]-t))

    return vx_vals
    
def vz_kink_pert(x, z, t, W, K, R1):
    vz_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vz_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix_kink(r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz_kink(r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vz_hat_vals[i,j] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB_kink*sc.sinh(m0(W)*sol[0]) +
                                       constC_kink(W, K, R1)*sc.cosh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vz_hat_vals[i,j] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA_kink(W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) + 
                                                                 sc.cosh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vz_hat_vals[i,j] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD_kink(W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                                             sc.cosh(m2(W)*sol[0]))
            
            vz_vals[i,j] = vz_hat_vals[i,j] * np.exp(1j*(z[j]-t))

    return vz_vals    
    