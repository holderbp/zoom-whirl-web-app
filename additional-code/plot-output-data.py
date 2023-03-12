import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np

# import parameters
X = np.loadtxt('pars.csv', delimiter=',', skiprows=1)
ell = X[0]
E = X[1]
ecc = X[2]
periap = X[3]
M_over_m = X[4]

# import orbit data
X = np.loadtxt('orbit.csv', delimiter=',', skiprows=1)
orb_t = X[:,0]
orb_r = X[:,1]
orb_phi = X[:,2]

# import grav wave data
X = np.loadtxt('grav-wave.csv', delimiter=',', skiprows=1)
gw_t = X[:,0]
gw_Hplus = X[:,1]
gw_Hcross = X[:,2]

# import effective potential data
X = np.loadtxt('eff-pot.csv', delimiter=',', skiprows=1)
ep_r = X[:,0]
ep_V = X[:,1]
ep_E = X[:,2]

# use gridspec to make 4 figures
gs = gs.GridSpec(7, 2) # nrows, ncols
fig = plt.figure(figsize=(8, 10.5))
# effective potential
ax1 = fig.add_subplot(gs[0:3,0])
ax1.plot(ep_r, ep_V)
ax1.plot(ep_r, ep_E)
ax1.set_xlabel(r'$r/M$')
ax1.set_ylabel(r'$V_{\rm eff}$')
ax1.text(0.8, 0.2, r"$\ell$" + f" = {ell:.6f}\n" + r"$E$" + f" = {E:.6f}",
         horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)
# orbit
ax3 = fig.add_subplot(gs[0:3,1], polar=True)
ax3.plot(orb_phi, orb_r)
ax3.text(0.8, 0, r"$e$" + f" = {ecc:.3f}\n" + r"$r_p$" + f" = {periap:.2f}",
         horizontalalignment='center',
         verticalalignment='top', transform=ax3.transAxes)
# Hplus
ax2 = fig.add_subplot(gs[3:5,:])
ax2.plot(gw_t, gw_Hplus)
ax2.set_xlabel(r'$\tau/M$')
ax2.set_ylabel(r'$H_+$')
ax2.set_title(f"M/m = {round(M_over_m):d}")
# Hcross
ax4 = fig.add_subplot(gs[5:7,:])
ax4.plot(gw_t, gw_Hcross)
ax4.set_xlabel(r'$\tau/M$')
ax4.set_ylabel(r'$H_x$')

fig.tight_layout()

# show plots
#plt.show()
# or save
plt.savefig('zw.pdf')

            
