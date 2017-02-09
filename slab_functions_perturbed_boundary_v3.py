
import numpy as np
import scipy as sc
from scipy.optimize import fsolve

# SBB
# Define the sound speeds and alfven speeds.
c2 = 1.2
c0 = 1.
vA = 0.9
cT = sc.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))

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


mode_options = ['slow kink surf', 'slow saus surf', 'slow saus body 3',
                'slow kink body 3', 'slow saus body 2', 'slow kink body 2', 
                'slow saus body 1', 'slow kink body 1', 'fast saus body 1',
                'fast kink body 1', 'fast saus body 2', 'fast kink body 1',
                'fast saus body 3', 'fast saus body 3', 'fast kink surf',
                'fast saus surf']
                
kink_mode_options = ['slow kink surf', 'slow kink body 1', 'slow kink body 2',
                     'slow kink body 3', 'fast kink body 1', 'fast kink body 2',
                     'fast kink body 3', 'fast kink surf']
saus_mode_options = ['slow saus surf', 'slow saus body 1', 'slow saus body 2',
                     'slow saus body 3', 'fast saus body 1', 'fast saus body 2',
                     'fast saus body 3', 'fast saus surf']
slow_surf_mode_options = ['slow kink surf', 'slow saus surf']
fast_surf_mode_options = ['fast kink surf', 'fast saus surf']
slow_body_1_mode_options = ['slow kink body 1', 'slow saus body 1']
slow_body_2_mode_options = ['slow kink body 2', 'slow saus body 2']
slow_body_3_mode_options = ['slow kink body 3', 'slow saus body 3']
fast_body_1_mode_options = ['fast kink body 1', 'fast saus body 1']
fast_body_2_mode_options = ['fast kink body 2', 'fast saus body 2']
fast_body_3_mode_options = ['fast kink body 3', 'fast saus body 3']

R2 = 2.


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

def lamb0(W):
    return -(vA**2-W**2)*1.j/(m0(W)*W)

def lamb00(W):
    return W*1.j/m00(W)

def lamb1(W, R1):
    return R1*W*1.j/m1(W, R1)
    
def lamb2(W):
    return R2*W*1.j/m2(W)

def required_xi(mode, K):
    if mode in slow_surf_mode_options:
        return K / 3.
    if mode in slow_body_1_mode_options:
        return K / 90. #30.
    if mode in slow_body_2_mode_options:
        return K / 110. #90.
    if mode in slow_body_3_mode_options:
        return K / 250.
    if mode in fast_body_1_mode_options:
        return K / 50.

def const(mode, W, K, R1):
    const_val_r = (W * required_xi(mode, K) / (constB_dash(mode, W, K, R1)*sc.cosh(m0(W)*K) +
                                     constC_dash(mode, W, K, R1)*sc.sinh(m0(W)*K)))
    const_val_l = (W * required_xi(mode, K) / (constB_dash(mode, W, K, R1)*sc.cosh(m0(W)*-K) +
                                     constC_dash(mode, W, K, R1)*sc.sinh(m0(W)*-K)))
    return min(const_val_r, const_val_l)

#def const(mode):
#    if mode in slow_surf_mode_options:
#        return 0.05
#    elif mode in slow_body_1_mode_options:
#        return 0.23
#    elif mode in slow_body_2_mode_options:
#        return 0.1
#    elif mode in fast_body_1_mode_options:
#        return 0.9
#    
def disp_rel_asym(W, K, R1):
    return ((W**4*m0(W)**2*R1*R2 + (vA**2 - W**2)**2*m1(W, R1)*m2(W) -
            0.5*m0(W)*W**2*(vA**2 - W**2)*(R2*m1(W, R1) + R1*m2(W))*
            (sc.tanh(m0(W)*K) + (sc.tanh(m0(W)*K))**(-1))) /
            (vA**2 - W**2)*(c0**2 - W**2)*(cT**2 - W**2))

    
def bz_eq_1d(mode, x, z, t, W, K, R1):
    truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                 (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bz_eqfunction = np.zeros(len(x))
    for i in indices:
        bz_eqfunction[i] = B0
    return bz_eqfunction
    
def bz_eq(mode, x, z, t, W, K, R1):
    bz_eq_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_eq_array[:,i] = bz_eq_1d(mode, x, z[i], t, W, K, R1)
    return bz_eq_array


##############################################

def constB(mode, W, K, R1):
    if mode in kink_mode_options:
        return const(mode, W, K, R1)
    elif mode in saus_mode_options:
        return const(mode, W, K, R1)*((lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)) /
                               (lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)))

def constC(mode, W, K, R1):
    if mode in kink_mode_options:
        return const(mode, W, K, R1)*((lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)) /
                            (lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)))
    elif mode in saus_mode_options:
        return const(mode, W, K, R1)
    
def constA(mode, W, K, R1):
    return ((constB(mode, W, K, R1)*sc.cosh(m0(W)*K) - constC(mode, W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m1(W, R1)*K) - sc.sinh(m1(W, R1)*K)))

def constD(mode, W, K, R1):
    return ((constB(mode, W, K, R1)*sc.cosh(m0(W)*K) + constC(mode, W, K, R1)*sc.sinh(m0(W)*K)) /
            (sc.cosh(m2(W)*K) - sc.sinh(m2(W)*K)))

# constants without const(mode). Used for finding appropriate value for const(mode)
def constB_dash(mode, W, K, R1):
    if mode in kink_mode_options:
        return 1.
    elif mode in saus_mode_options:
        return ((lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)) /
                               (lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)))

def constC_dash(mode, W, K, R1):
    if mode in kink_mode_options:
        return ((lamb1(W, R1)*sc.cosh(m0(W)*K)+lamb0(W)*sc.sinh(m0(W)*K)) /
                            (lamb0(W)*sc.cosh(m0(W)*K)+lamb1(W, R1)*sc.sinh(m0(W)*K)))
    elif mode in saus_mode_options:
        return 1.

def vxhat(mode, x, W, K, R1):
    if type(x) == np.float64:
        if np.abs(x) <= K:
            vxhatfunction = (constB(mode, W, K, R1)*sc.cosh(m0(W)*x) + 
                            constC(mode, W, K, R1)*sc.sinh(m0(W)*x))
        elif x < -K:
            vxhatfunction = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*x) + 
                                                  sc.sinh(m1(W, R1)*x))
        elif x > K:
            vxhatfunction = constD(mode, W, K, R1)*(sc.cosh(m2(W)*x) - 
                                                  sc.sinh(m2(W)*x))
    else:
        truth = np.array(np.abs(x) <= K*np.ones(len(x)))
        indices = np.where(truth == True)
        vxhatfunction = np.zeros(len(x), dtype=complex)
        for i in indices:
            vxhatfunction[i] = (constB(mode, W, K, R1)*sc.cosh(m0(W)*x[i]) + 
                                constC(mode, W, K, R1)*sc.sinh(m0(W)*x[i]))
        truth2 = np.array(x < -K*np.ones(len(x)))
        indices2 = np.where(truth2 == True)
        for i in indices2:
            vxhatfunction[i] = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*x[i]) + 
                                                      sc.sinh(m1(W, R1)*x[i]))
        truth3 = np.array(x > K*np.ones(len(x)))
        indices3 = np.where(truth3 == True)
        for i in indices3:
            vxhatfunction[i] = constD(mode, W, K, R1)*(sc.cosh(m2(W)*x[i]) - 
                                                      sc.sinh(m2(W)*x[i]))
    return vxhatfunction
    
def vzhat(mode, x, W, K, R1):
    if type(x) == np.float64:
        if np.abs(x) <= K:
            vzhat_function = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x) +
                                       constC(mode, W, K, R1)*sc.cosh(m0(W)*x))
        elif x < -K:
            vzhat_function = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x) + 
                                                                 sc.cosh(m1(W, R1)*x))
        elif x > K:
            vzhat_function = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*x) -
                                                             sc.cosh(m2(W)*x))
    else:
        truth = np.array(np.abs(x) <= K*np.ones(len(x)))
        indices = np.where(truth == True)
        vzhat_function = np.zeros(len(x), dtype=complex)
        for i in indices:
            vzhat_function[i] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                                           constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i]))
        truth2 = np.array(x < -K*np.ones(len(x)))
        indices2 = np.where(truth2 == True)
        for i in indices2:
            vzhat_function[i] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) + 
                                                                     sc.cosh(m1(W, R1)*x[i]))
        truth3 = np.array(x > K*np.ones(len(x)))
        indices3 = np.where(truth3 == True)
        for i in indices3:
            vzhat_function[i] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                                                                 sc.cosh(m2(W)*x[i]))
    return vzhat_function
    
def vx(mode, x, z, t, W, K, R1):
    if type(vxhat(mode, x, W, K, R1)) == np.complex128:
        return vxhat(mode, x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vxhat(mode, x, W, K, R1), np.exp(1j*(z-t)))

def vz(mode, x, z, t, W, K, R1):
    if type(vzhat(mode, x, W, K, R1)) == np.complex128:
        return vzhat(mode, x, W, K, R1) * np.exp(1j*(z-t))
    else:
        return np.outer(vzhat(mode, x, W, K, R1), np.exp(1j*(z-t)))
        
def vx_pert(mode, x, z, t, W, K, R1):
    vx_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vx_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vx_hat_vals[i,j] = (constB(mode, W, K, R1)*sc.cosh(m0(W)*sol[0]) + 
                            constC(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vx_hat_vals[i,j] = constA(mode, W, K, R1)*(sc.cosh(m1(W, R1)*sol[0]) + 
                                                  sc.sinh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vx_hat_vals[i,j] = constD(mode, W, K, R1)*(sc.cosh(m2(W)*sol[0]) - 
                                                  sc.sinh(m2(W)*sol[0]))
            vx_vals[i,j] = vx_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return vx_vals
    
def vz_pert(mode, x, z, t, W, K, R1):
    vz_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    vz_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                vz_hat_vals[i,j] = (1j * c0**2 / (c0**2 - W**2)) * m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]) +
                                       constC(mode, W, K, R1)*sc.cosh(m0(W)*sol[0]))
            elif sol[0] < -K:
                vz_hat_vals[i,j] = (1j * c1(R1)**2 / (c1(R1)**2 - W**2)) * m1(W, R1)*constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) + 
                                                                 sc.cosh(m1(W, R1)*sol[0]))
            elif sol[0] > K:
                vz_hat_vals[i,j] = (1j * c2**2 / (c2**2 - W**2)) * m2(W)*constD(mode, W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                                             sc.cosh(m2(W)*sol[0]))
            vz_vals[i,j] = vz_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return vz_vals    
    
def xix_hat(mode, x, W, K, R1):
    return vxhat(mode, x, W, K, R1) / W
    
def xiz_hat(mode, x, W, K, R1):
    return vzhat(mode, x, W, K, R1) / W

def xix(mode, x, z, t, W, K, R1):
    if type(vxhat(mode, x, W, K, R1)) == np.complex128:
        return (1j * vxhat(mode, x, W, K, R1) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vxhat(mode, x, W, K, R1) / W, np.exp(1j*(z-t)))
    
def xiz(mode, x, z, t, W, K, R1):
    if type(vzhat(mode, x, W, K, R1)) == np.complex128:
        return (1j * vzhat(mode, x, W, K, R1) / W) * np.exp(1j*(z-t))
    else:
        return np.outer(1j * vzhat(mode, x, W, K, R1) / W, np.exp(1j*(z-t)))

def xixhat_boundary(mode, W, K, R1, boundary='r'):
    if boundary == 'r' or boundary == 'right':
        xixhat = (1j / W) * (constB(mode, W, K, R1)*sc.cosh(m0(W)*K) + 
                                 constC(mode, W, K, R1)*sc.sinh(m0(W)*K))
    if boundary == 'l' or boundary == 'left':
        xixhat = (1j / W) * (constB(mode, W, K, R1)*sc.cosh(m0(W)*-K) + 
                                 constC(mode, W, K, R1)*sc.sinh(m0(W)*-K))
    return xixhat

def xix_boundary(mode, z, t, W, K, R1, boundary='r'):           
    return xixhat_boundary(mode, W, K, R1, boundary) * np.exp(1j*(z-t))

def bxhat(mode, x, z, t, W, K, R1):
    truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                     (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bxhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bxhat_function[i] = (-B0/W)*(constB(mode, W, K, R1)*sc.cosh(m0(W)*x[i]) +
                                   constC(mode, W, K, R1)*sc.sinh(m0(W)*x[i]))
    return bxhat_function
    
def bzhat(mode, x, z, t, W, K, R1):
    truth = np.array((x <= (K + xix_boundary(mode, z, t, W, K, R1, boundary='r'))*np.ones(len(x))) &
                     (x >= (-K + xix_boundary(mode, z, t, W, K, R1, boundary='l'))*np.ones(len(x))))
    indices = np.where(truth == True)
    bzhat_function = np.zeros(len(x), dtype=complex)
    for i in indices:
        bzhat_function[i] = (-1j*B0/W)*m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                                   constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i]))
    return bzhat_function
    
def bx(mode, x, z, t, W, K, R1):
    bx_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bx_array[:,i] = bxhat(mode, x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bx_array
    
def bz(mode, x, z, t, W, K, R1):
    bz_array = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(z)):
        bz_array[:,i] = bzhat(mode, x, z[i], t, W, K, R1) * np.exp(1j*(z[i]-t))
    return bz_array

# In the slab has an extra -1 factor to make it work. Not sure why it needs this.
def rho_hat(mode, x, W, K, R1):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    rho_hatfunction = np.zeros(len(x), dtype=complex)
    for i in indices:
        rho_hatfunction[i] = m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*x[i]) +
                             constC(mode, W, K, R1)*sc.cosh(m0(W)*x[i])) * lamb00(W) / (c0**2 * m00(W))
#        if mode in slow_surf_mode_options + slow_body_1_mode_options + slow_body_2_mode_options:
#            rho_hatfunction[i] = -rho_hatfunction[i]
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        rho_hatfunction[i] = constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*x[i]) +
                             sc.cosh(m1(W, R1)*x[i])) * lamb1(W, R1) / c1(R1)**2
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        rho_hatfunction[i] = constD(mode, W, K, R1)*(sc.sinh(m2(W)*x[i]) -
                             sc.cosh(m2(W)*x[i])) * lamb2(W) / c2**2
    return rho_hatfunction
    
def rho(mode, x, z, t, W, K, R1):
    rho_func = np.outer(rho_hat(mode, x, W, K, R1), np.exp(1j*(z-t)))
    return rho_func

# for plotting rho(x,t) at (x + xi,t)
def rho_pert(mode, x, z, t, W, K, R1):
    rho_hat_vals = np.zeros((len(x), len(z)), dtype=complex)
    rho_vals = np.zeros((len(x), len(z)), dtype=complex)
    for i in range(len(x)):
        for j in range(len(z)):
            def func(r):
                return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
                        r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
            sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
            if abs(sol[0]) <= K:
                rho_hat_vals[i,j] = m0(W)*(constB(mode, W, K, R1)*sc.sinh(m0(W)*sol[0]) +
                                constC(mode, W, K, R1)*sc.cosh(m0(W)*sol[0])) * lamb00(W) / (c0**2 * m00(W))
                if mode in slow_surf_mode_options + slow_body_1_mode_options + slow_body_2_mode_options:
                    rho_hat_vals[i,j] = -rho_hat_vals[i,j]
            elif sol[0] < -K:
                rho_hat_vals[i,j] = constA(mode, W, K, R1)*(sc.sinh(m1(W, R1)*sol[0]) +
                                sc.cosh(m1(W, R1)*sol[0])) * lamb1(W, R1) / c1(R1)**2
            elif sol[0] > K:
                rho_hat_vals[i,j] = constD(mode, W, K, R1)*(sc.sinh(m2(W)*sol[0]) -
                                     sc.cosh(m2(W)*sol[0])) * lamb2(W) / c2**2
            
            rho_vals[i,j] = rho_hat_vals[i,j] * np.exp(1j*(z[j]-t))
    return rho_vals
    