import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt

#
# See my notes 2023-02-18 (actually this update from 2023-03-05)
#
#   I wanted to get the inverse transform (E, ell) -> (p, e)
#   where p is the semilatus rectum and e is eccentricity. The
#   forward transform gives equations E(p,e) and ell(p,e) using
#   (r-, r+) and E - V(r) = 0.  I inverted these to get this
#   cubic equation for p(E, ell).  It turns out there are
#   multiple roots in the relevant region, so it's not so
#   straightforward to find the right one (have to provide
#   random guesses and then check the found root against the
#   known periapsis, r-).
#
#   BUT, all of this is nonsense, because the periapsis and
#   apoapsis can be found easily by rootfinding E - V(r) = 0
#   (see below), and then the eccentricity can be found by the
#   defining equations of (p, e):
#
#      r- = p/(1+e) ; r+ = p/(1-e)  -> e = (r+ - r-)/(r+ + r-)
#
#   so there is no necessity to find p on its own.  Duh.
#

G = 1
M = 1
very_large_r = 1e4

# starting point of simulation
ell = 3.7
E = -0.0355

def get_Veff_maxmin_r_values():
    if ell**2 > 12*(G*M)**2:
        radical = np.sqrt(ell**2 - 12*(G*M)**2)
        return [ell**2 / (2*G*M) - ell / (2*G*M) * radical,
                ell**2 / (2*G*M) + ell / (2*G*M) * radical]
    elif ell**2 == 12*(G*M)**2:
        # at the minimum value of angular momentum
        return [ell**2/(2*G*M), ell**2/(2*G*M)]
    else:
    	# no solutions
        return [None, None]

def Veff(r, getderiv=False):
    if getderiv:
        return G*M/r**2 - ell**2/r**3 + 3*G*M*ell**2/r**4
    else:
        return -G*M/r + ell**2/(2*r**2) - G*M*ell**2 / r**3

def Veff_minus_E(r):
    return Veff(r) - E

def func(p):
    """
    one of the roots of this equation is the value of
    the semilatus rectum, p, given (ell, E)
    """
    A = (2*E + 1)
    B = -(ell**2 + 4)
    C = 8*ell**2
    D = -16*ell**2
    return A * p**3 + B * p**2 + C * p + D

def ecc(p):
    return np.sqrt(p - p**2/ell**2 - 3)

while True:
    # get energy and angular momentum values
    ell = input("Enter angular momentum: ")
    E = input("Enter energy: ")
    ell = float(ell)
    E = float(E)
    print("-----------------------------------------------")
    print("Using ell =", ell, " and E =", E)
    # get peri- and apo-apsis
    [IUCO, ISCO] = get_Veff_maxmin_r_values()
    rp = spo.bisect(Veff_minus_E, a=IUCO, b=ISCO, disp=True)
    # bisect with very large value to get apoapsis
    ra = spo.bisect(Veff_minus_E, a=ISCO, b=very_large_r, disp=True)
    print("IUCO =", IUCO, "ISCO =", ISCO)
    print("rp =",rp, " ra =", ra)
    #
    # get p
    #
    #  p = (1+e)*rp, so if e=0, p=rp and if e=1, p=2*rp
    #
    while True:
        guess = np.random.rand() * (2*rp - rp) + rp
        print("* Trying to find root with guess =", guess)
        p = spo.fsolve(func, guess)
        p = p[0]
        if (p < rp) | (p > 2*rp):
            print("Root found, but outside of [rp, 2*rp]")
            print("p =", p)
        else:
            e = ecc(p)
            rpp = p/(1+e)
            if (np.abs(rp - rpp)/np.average([rp, rpp]) < 0.01):
                break
            else:
                e = ecc(p)
                rpp = p/(1+e)
                print("Root found, but it is wrong:")
                print("p =", p, "e =", e, "rpp =", rpp)
    print("Correct root found:")
    print("p =", p, " e =", e, " rperi =", rpp)
    print("-----------------------------------------------")    
    # plot the p function
    dr = 0.01
    rs = np.linspace(rp, 2*rp, num=1000)
    Vs = Veff(rs)
    pfunc = func(rs)
    plt.plot(rs, pfunc)
    plt.plot(rs, 0*rs)
    plt.show()
