
import numpy as np
import scipy as sc

# SBB
# Define the sound speeds and alfven speeds.
c2 = 1.2
c0 = 1.
vA = 0.9
cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))
ce = 1.2 #0.7


## for xi of x slow surface sf GS and maybe others
#c2 = 0.7
#c0 = 1.
#vA = 0.4
#cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))


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

#constC_s = 0.08
#constB_kink = 0.08
#constC_s_sym = 0.08
#constB_kink_sym = 0.08

# for xi of x slow surface sf GS and maybe others
constC_s = 0.05
constB_kink = 0.05
constC_s_sym = 0.05
constB_kink_sym = 0.05

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
def m1(W, R1):
    return sc.sqrt(1 - W**2*(c2*sc.sqrt(R2/R1))**(-2))
    
def m2(W):
    return sc.sqrt(1 - W**2/c2**2)

def me(W):
    return sc.sqrt(1 - W**2*ce**(-2))
    
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

def lamb1(W, R1):
    return R1*W*1.j/m1(W, R1)
    
def lamb2(W):
    return R2*W*1.j/m2(W)

def lambe(W):
    return Re*W*1j/me(W)
    
def bz_eq(x, K):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    bz_eqfunction = np.zeros(len(x))
    for i in indices:
        bz_eqfunction[i] = B0
    return bz_eqfunction

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

def bxhat_kink(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    bxhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bxhat_function[i] = (-B0/W)*m0(W)*(constB_kink*sc.cosh(m0(W)*x[i]) +
                                   constC_kink(W, K, R1)*sc.sinh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        bxhat_function[i] = 0.
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        bxhat_function[i] = 0.
    return bxhat_function
    
def bzhat_kink(x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    bzhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bzhat_function[i] = (-1j*B0/W)*m0(W)*(constB_kink*sc.sinh(m0(W)*x[i]) +
                                   constC_kink(W, K, R1)*sc.cosh(m0(W)*x[i]))
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        bzhat_function[i] = 0.
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        bzhat_function[i] = 0.
    return bzhat_function
    
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
    return np.outer(vxhat_kink(x, W, K, R1), np.exp(1j*(z-t)))

def vx_dash_kink(x, z, t, W, K, R1):
    return np.outer(vxhat_dash_kink(x, W, K, R1), np.exp(1j*(z-t)))

def bx_kink(x, z, t, W, K, R1):
    return np.outer(bxhat_kink(x, W, K, R1), np.exp(1j*(z-t)))
    
def bz_kink(x, z, t, W, K, R1):
    return np.outer(bzhat_kink(x, W, K, R1), np.exp(1j*(z-t)))
    
def p_tot_kink(x, z, t, W, K, R1):
    return np.outer(p_tothat_kink(x, W, K, R1), np.exp(1j*(z-t)))