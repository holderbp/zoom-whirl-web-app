import numpy as np
import scipy.integrate as spi
import scipy.optimize as spo

#
#--- Module parameters
#
G = 1  # using "geometrized" units where G=c=1 but left "G" in for clarity
M = 1  # set to unity such that all dimensional quantities in units of M (leave "M" in for clarity)
M_over_m_default = 100000  # sets scale of GW amplitude
ti = 0
tf_default = 1000
very_large_r = 100000

# create potential graph data
def Veff(r, ell, getderiv=False):
    if getderiv:
        return G*M/r**2 - ell**2/r**3 + 3*G*M*ell**2/r**4
    else:
        return -G*M/r + ell**2/(2*r**2) - G*M*ell**2 / r**3

def VNeff(r, ell):
    return -G*M/r + ell**2/(2*r**2)

def Veff_minus_E(r, ell, E):
    return Veff(r, ell) - E

def get_Veff_maxmin_r_values(ell):
    if ell**2 > 12*(G*M)**2:
        radical = np.sqrt(ell**2 - 12*(G*M)**2)
        # inner-unstable and outer-stable circular orbits
        return [
            ell**2 / (2*G*M) - ell / (2*G*M) * radical,
            ell**2 / (2*G*M) + ell / (2*G*M) * radical
            ]
    elif ell**2 == 12*(G*M)**2:
        # at the minimum value of angular momentum
        return [ell**2/(2*G*M), ell**2/(2*G*M)]
    else:
        # no solutions
        return [None, None]

def derivs(t, X, ell):
    r = X[0]
    vr = X[1]
    phi = X[2]
    drdt = vr
    dvdt = -Veff(r, ell, getderiv=True)
    dphidt = ell/r**2
    return [drdt, dvdt, dphidt]

def get_peri_and_apoapsis(ell, E):
    # find r values of circular orbits (rp is in between)
    [IUCO, ISCO] = get_Veff_maxmin_r_values(ell)
    rp = spo.bisect(Veff_minus_E, a=IUCO, b=ISCO, disp=True, args=(ell, E))
    # then bisect with very large value
    ra = spo.bisect(Veff_minus_E, a=ISCO, b=very_large_r, disp=True, args=(ell, E))
    return rp, ra

def get_rwell(Vpeak, ell):
    if (ell**2 < 12*(G*M)**2):
        return None
    elif (ell > 4):
        E = 0.0
        [IUCO, ISCO] = get_Veff_maxmin_r_values(ell)
        rwell = spo.bisect(Veff_minus_E, a=IUCO, b=ISCO, disp=True, args=(ell, E))
        return rwell
    else:
        E = Vpeak
        [IUCO, ISCO] = get_Veff_maxmin_r_values(ell)
        rwell = spo.bisect(Veff_minus_E, a=ISCO, b=very_large_r, disp=True, args=(ell, E))
        return rwell
    
def get_eccentricity(rp, ra):
    return (ra - rp) / (rp + ra)
    
def get_orbit(ell, E, tmax):
    #
    # Note: User must set zwoc.ell and zwoc.E first
    #
    #--- Get the value of r at periapsis
    #
    [rp, ra] = get_peri_and_apoapsis(ell, E)
    #
    #--- set initial conditions, always at periapsis
    #
    X0 = [rp, 0.0, 0.0]
    #
    #--- Times to evaluate function
    #
    Ntsteps = int(tmax-ti)+1
    evaltimes = np.linspace(ti, tmax, num=Ntsteps)
    #
    #--- evolve
    #
    sol = spi.solve_ivp(derivs, [ti, tmax], X0, t_eval = evaltimes,
                        rtol = 1e-8, atol = 1e-8, args = (ell,))
    #
    #
    #  Fixme: I removed event tracking of apoapsis... do we need it?
    #
    # get eccentricity
    ecc = get_eccentricity(rp, ra)
    # grab solution vectors
    [r_t, vr_t, phi_t] = sol.y
    # and return all info
    return evaltimes, r_t, phi_t, rp, ra, ecc

# get Iddot for gravitational wave plotting
def get_Iddot(t, r, phi, M_over_m):
    m = M/M_over_m
    N = len(r)
    XY = np.zeros(shape=(2, N))
    I = np.zeros(shape=(2, 2, N))
    I_dot = np.zeros(shape=(2, 2, N - 2))
    I_ddot = np.zeros(shape=(2, 2, N - 4))
    for i in range(N):
        XY[0][i] = r[i]*np.cos(phi[i])
        XY[1][i] = r[i]*np.sin(phi[i])
        I[0][0][i] = m*((XY[0][i]**2) - (1/3)*(r[i]**2))
        I[1][1][i] = m*((XY[1][i]**2) - (1/3)*(r[i]**2))
        I[0][1][i] = 2*m*(XY[0][i]*XY[1][i])
        I[1][0][i] = 2*m*(XY[0][i]*XY[1][i])
    for i in range(N - 2):
        I_dot[0][0][i] = (I[0][0][i + 2] - I[0][0][i])/(t[i + 2] - t[i])
        I_dot[1][1][i] = (I[1][1][i + 2] - I[1][1][i]) / (t[i + 2] - t[i])
        I_dot[1][0][i] = (I[1][0][i + 2] - I[1][0][i]) / (t[i + 2] - t[i])
        I_dot[0][1][i] = (I[0][1][i + 2] - I[0][1][i]) / (t[i + 2] - t[i])
    for i in range(N - 4):
        I_ddot[0][0][i] = (I_dot[0][0][i + 2] - I_dot[0][0][i])/(t[i + 2] - t[i])
        I_ddot[1][1][i] = (I_dot[1][1][i + 2] - I_dot[1][1][i]) / (t[i + 2] - t[i])
        I_ddot[1][0][i] = (I_dot[1][0][i + 2] - I_dot[1][0][i]) / (t[i + 2] - t[i])
        I_ddot[0][1][i] = (I_dot[0][1][i + 2] - I_dot[0][1][i]) / (t[i + 2] - t[i])
    return I_ddot
