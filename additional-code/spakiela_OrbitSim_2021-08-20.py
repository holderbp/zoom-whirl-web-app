import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import fsolve, newton, bisect
from matplotlib.widgets import Slider, TextBox, Button
import matplotlib.gridspec as gridspec
from fractions import Fraction as frac
from matplotlib.animation import FuncAnimation
import time


# Known issues:
# When scrolling to zoom in or out at the end of an animation, the project closes
# Gravitational wave plots do not scale correctly when updating to a new orbit
# If pause button is hit twice, plots are cleared


# Start Time
t_i = 0
# End Time
t_f = 1000
# Step Interval = h
Step = 1
# Time array
Time = np.arange(t_i, t_f, Step)
# Mass of massive body
M = 1
# Mass of particle (For radiation and waveform)
m = M/100000
# Gravitational Constant
G = 1
# Angular Momentum
L = 4
# Initial radius
r_i = 4
r_p = r_i
# Initial radial velocity
rdot_i = 0
# Initial Phi
phi_i = 0
# Radius of innermost stable circular orbit (Moore 10.11)
ISCO = (6 * G * M) / (1 - np.sqrt(1 - 12 * (G * M / L) ** 2))
# Radius of innermost unstable circular orbit (Moore 10.11)
IUCO = (6 * G * M) / (1 + np.sqrt(1 - 12 * (G * M / L) ** 2))


# Creating the effective potential function
r = sp.Symbol('r')
U_Eff = (-G*M/r + L**2/(2*r**2) - G*M*(L**2)/r**3)
U_Eff_Func = sp.lambdify(r, U_Eff)
Eff_Force = -sp.diff(U_Eff, r)
Eff_Force_Func = sp.lambdify(r, Eff_Force)
Phi_dot = L/r**2
Phi_dot_Func = sp.lambdify(r, Phi_dot)
U_n_Eff = -G*M/r + L**2/(2*r**2)
U_n_Eff_Func = sp.lambdify(r, U_n_Eff)


# Instantiating E
E = U_Eff_Func(r_i)


# Plotting Effective Potential
Ueff_Array = np.linspace(1, 50, 1000)
Zeros_Array = np.zeros_like(Ueff_Array)
for i in range(Zeros_Array.size):
    Zeros_Array[i] = U_Eff_Func(Ueff_Array[i])


# Creating a line to represent particles energy
def energy_line(Radius):
    y1 = np.zeros_like(Ueff_Array)
    for i in range(Ueff_Array.size):
        y1[i] = U_Eff_Func(Radius)
        # y1[i] = U_n_Eff_Func(Radius)
    return y1


# Python function to find the root of the effective potential
root1 = sp.lambdify(r,  E - G*M/r + L**2/(2*r**2) - G*M*(L**2)/r**3)


# Method to calculate eccentricity
def get_e(p, a):
    return (a - p) / (p + a)


# Integration method for solve_ivp
def deriv(t, y):
    return [y[1], Eff_Force_Func(y[0]), Phi_dot_Func(y[0])]


# event tracking method to track apoapsis terminal
def apoapsis(t, y):
    return y[1]

# event tracking method to track apoapsis nonterminal
def apoapsis_nt(t, y):
    return y[1]


apoapsis.terminal = True
apoapsis.direction = -1
apoapsis_nt.terminal = False
apoapsis_nt.direction = -1


# Array of initial conditions
y0 = [r_i, rdot_i, phi_i]
# Array of solutions
sol = solve_ivp(deriv, y0=y0, t_span=[t_i, t_f], t_eval=Time, rtol=1e-8, atol=1e-8, events=apoapsis_nt)

Ecc = get_e(r_i, sol.y_events[0][0][0])


# Creates I_ddot tensor
XY = np.zeros(shape=(2, sol.y[0].size))
I = np.zeros(shape=(2, 2, sol.y[0].size))
I_dot = np.zeros(shape=(2, 2, sol.y[0].size - 2))
I_ddot = np.zeros(shape=(2, 2, sol.y[0].size - 4))


# Method to calculate gravitational waves
def get_H():
    global XY, I, I_dot, I_ddot
    for i in range(sol.y[0].size):
        XY[0][i] = sol.y[0][i]*np.cos(sol.y[2][i])
        XY[1][i] = sol.y[0][i]*np.sin(sol.y[2][i])
        I[0][0][i] = m*((XY[0][i]**2) - (1/3)*(sol.y[0][i]**2))
        I[1][1][i] = m*((XY[1][i]**2) - (1/3)*(sol.y[0][i]**2))
        I[0][1][i] = 2*m*(XY[0][i]*XY[1][i])
        I[1][0][i] = 2*m*(XY[0][i]*XY[1][i])
    for i in range(sol.y[0].size - 2):
        I_dot[0][0][i] = (I[0][0][i + 2] - I[0][0][i])/(sol.t[i + 2] - sol.t[i])
        I_dot[1][1][i] = (I[1][1][i + 2] - I[1][1][i]) / (sol.t[i + 2] - sol.t[i])
        I_dot[1][0][i] = (I[1][0][i + 2] - I[1][0][i]) / (sol.t[i + 2] - sol.t[i])
        I_dot[0][1][i] = (I[0][1][i + 2] - I[0][1][i]) / (sol.t[i + 2] - sol.t[i])
    for i in range(sol.y[0].size - 4):
        I_ddot[0][0][i] = (I_dot[0][0][i + 2] - I_dot[0][0][i])/(sol.t[i + 2] - sol.t[i])
        I_ddot[1][1][i] = (I_dot[1][1][i + 2] - I_dot[1][1][i]) / (sol.t[i + 2] - sol.t[i])
        I_ddot[1][0][i] = (I_dot[1][0][i + 2] - I_dot[1][0][i]) / (sol.t[i + 2] - sol.t[i])
        I_ddot[0][1][i] = (I_dot[0][1][i + 2] - I_dot[0][1][i]) / (sol.t[i + 2] - sol.t[i])


# Calling gravitational wave method to fill array to be plotted
get_H()


# Creates Plotting window
fig1 = plt.figure(num='Orbit Applet', figsize=(16, 9))
plt.subplots_adjust(bottom=.25, left=.2)
spec1 = gridspec.GridSpec(nrows=2, ncols=2, figure=fig1)

# Plotting orbit
ax1 = fig1.add_subplot(spec1[0, 0], projection='polar')
# ax1.set_title('Orbit')
orbit, = plt.polar(sol.y[2], sol.y[0])
orbit_trace, = plt.polar(sol.y[2], sol.y[0], alpha=.5, linestyle="--", color="b")
orbit_dot, = plt.polar(sol.y[2][0], sol.y[0][0], "or")
plt.ylim(0, 20)

# Plotting GW
ax3 = fig1.add_subplot(spec1[1, 0])
ax3.set_title('H+')
h_xx, = plt.plot(sol.t[10: - 4], I_ddot[0][0][10:])
h_xx_dot, = plt.plot(sol.t[10], I_ddot[0][0][10], "or")
plt.xlabel('t')
plt.xlim(0, t_f)

ax4 = fig1.add_subplot(spec1[1, 1])
ax4.set_title('Hx')
h_xy, = plt.plot(sol.t[10: - 4], I_ddot[1][0][10:])
h_xy_dot, = plt.plot(sol.t[10], I_ddot[1][0][10], "or")
plt.xlabel('t')
plt.xlim(0, t_f)

# Plotting Effective Potential
ax2 = fig1.add_subplot(spec1[0, 1])
ax2.set_title('Effective Potential')
a1, = plt.plot(Ueff_Array, Zeros_Array)
plt.xlabel('r/M')
plt.ylabel('E')
plt.ylim(-.05, .05)

# Adds periapsis on Effective Potential
a2, = plt.plot(r_i, U_Eff_Func(r_i), "or")

# Adds Energy Line on Effective Potential
a3, = plt.plot(Ueff_Array, energy_line(r_i), "--r", alpha=.3)

# Add slider for energy
ax_E = plt.axes([.45, .15, .4, .02], facecolor='lightgoldenrodyellow')
s_E = Slider(ax_E, 'Energy', -.05, .05, valinit=U_Eff_Func(r_i), valstep=.00001)

# Add slider for angular momentum
ax_L = plt.axes([.45, .1, .4, .02], facecolor='lightgoldenrodyellow')
s_L = Slider(ax_L, 'Ang Momentum', 0, 10, valinit=L, valstep=.01)

# Add Slider for eccentricity
ax_e = plt.axes([.05, .25, .012, .5], facecolor='lightgoldenrodyellow')
s_e = Slider(ax_e, 'Eccentricity', 0, 1, valinit=Ecc, valstep=.001, orientation='vertical')

# Add Slider for periapsis
ax_rp = plt.axes([.1, .25, .012, .5], facecolor='lightgoldenrodyellow')
s_rp = Slider(ax_rp, 'Periapsis', 0, 10, valinit=r_p, valstep=.01, orientation='vertical')

# Add text box to manually input energy
ax_text_L = plt.axes([.25, .15, .1, .025])
text_bot_E = TextBox(ax_text_L, ' ', initial=str(int(U_Eff_Func(r_i) * (10 ** 5)) / (10 ** 5)))

# Add text box to manually input angular momentum
ax_text_E = plt.axes([.25, .1, .1, .025])
text_bot_L = TextBox(ax_text_E, ' ', initial=str(L))

# Add Pause and Play button
ax_pause = plt.axes([.05, .1, .05, .05])
ax_play = plt.axes([.11, .1, .05, .05])
pause = Button(ax_pause, '||')
play = Button(ax_play, '|>')


# Function to zoom in and out of orbit plot w/ scroll wheel
#   Need to fix scrolling ending program
def zoom(ax, scale):
    def zoom_event(event):
        cur_ylim = ax.get_ylim()
        if event.button == 'up':
            scale_fac = 1/scale
        elif event.button == 'down':
            scale_fac = scale
        else:
            scale_fac = 1
        ax.set_ylim([0, cur_ylim[1]*scale_fac])
        animation.event_source.stop()
        plt.draw()
        animation.event_source.start()
    fig = ax.get_figure()
    fig.canvas.mpl_connect('scroll_event', zoom_event)
    return zoom_event


f = zoom(ax1, 1.5)


# Method to update all global variables used for recalculating orbits from slider values
def update_globals():
    global ISCO, IUCO, root1, r_i,  U_Eff, U_Eff_Func, Eff_Force, Eff_Force_Func, Phi_dot, Phi_dot_Func
    ISCO = (6 * G * M) / (1 - np.sqrt(1 - 12 * (G * M / L) ** 2))
    IUCO = (6 * G * M) / (1 + np.sqrt(1 - 12 * (G * M / L) ** 2))
    U_Eff = (-(G * M / r) + ((L ** 2) / (2 * (r ** 2))) - ((G * M * (L ** 2)) / (r ** 3)))
    U_Eff_Func = sp.lambdify(r, U_Eff)
    Eff_Force = -sp.diff(U_Eff, r)
    Eff_Force_Func = sp.lambdify(r, Eff_Force)
    Phi_dot = L / r ** 2
    Phi_dot_Func = sp.lambdify(r, Phi_dot)


# Function to update plots when energy is changed
def update_e(val):
    update_figures(s_E.val, L, r_p, Ecc)


# Function to change angular momentum
def update_l(val):
    update_figures(E, s_L.val, r_p, Ecc)


# Method to update periapsis
def update_rp(val):
    update_figures(E, L, s_rp.val, Ecc)


# Method to update eccentricity
def update_ecc(val):
    update_figures(E, L, r_p, s_e.val)


# Method to update all figures with new orbital information
def update_figures(Eval, Lval, rpval, eccval):
    global E, L, ISCO, IUCO, root1, r_i, r_p, Ecc, U_Eff, U_Eff_Func, Eff_Force, Eff_Force_Func, Phi_dot, Phi_dot_Func, y0, sol, animation
    start = time.time()
    animation.event_source.stop()
    # Check if call originated from periapsis slider
    if rpval != r_p:
        print("Rp slider was updated")
        r_p = rpval
        E = (M * (1 - Ecc) * (4 * M - (1 + Ecc) * r_p)) / (2 * r_p * ((1 + Ecc) * r_p - (3 + Ecc ** 2) * M))
        L = ((1 + Ecc) * r_p) / (np.sqrt((1 + Ecc) * (r_p / M) - (3 + (Ecc ** 2))))
        update_globals()
        # Handles problem with finding circular orbits due to floating point values
        if np.round(E, 5) == np.round(U_Eff_Func(IUCO), 5):
            r_i = IUCO
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        elif np.round(E, 5) == np.round(U_Eff_Func(ISCO), 5):
            r_i = ISCO
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        else:
            root1 = sp.lambdify(r, E + (G * M / r) - ((L ** 2) / (2 * (r ** 2))) + ((G * M * (L ** 2)) / (r ** 3)))
            r_i = bisect(root1, a=IUCO, b=ISCO, disp=True)
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        s_L.set_val(L)
    # Check if call was made from eccentricity slider
    elif eccval != Ecc:
        print("Ecc slider was updated")
        Ecc = eccval
        E = (M * (1 - Ecc) * (4 * M - (1 + Ecc) * r_p)) / (2 * r_p * ((1 + Ecc) * r_p - (3 + Ecc ** 2) * M))
        L = ((1 + Ecc) * r_p) / (np.sqrt((1 + Ecc) * (r_p / M) - (3 + (Ecc ** 2))))
        update_globals()
        # Handles problem with finding circular orbits due to floating point values
        if np.round(E, 5) == np.round(U_Eff_Func(IUCO), 5):
            r_i = IUCO
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        elif np.round(E, 5) == np.round(U_Eff_Func(ISCO), 5):
            r_i = ISCO
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        else:
            root1 = sp.lambdify(r, E + (G * M / r) - ((L ** 2) / (2 * (r ** 2))) + ((G * M * (L ** 2)) / (r ** 3)))
            r_i = bisect(root1, a=IUCO, b=ISCO, disp=True)
            s_E.set_val(E)
            animation = FuncAnimation(fig1, update, blit=True, interval=10)
            animation.event_source.start()
        s_L.set_val(L)
    # Check if call was made from angular momentum slider
    elif Lval != L:
        print("L slider was updated")
        L = Lval
        update_globals()
        root1 = sp.lambdify(r, E + (G * M / r) - ((L ** 2) / (2 * (r ** 2))) + ((G * M * (L ** 2)) / (r ** 3)))
        r_i = bisect(root1, a=IUCO, b=ISCO, disp=True)
        r_a = bisect(root1, a=ISCO, b=500000, disp=True)
        Ecc = get_e(r_i, r_a)
        r_p = r_i
        s_rp.set_val(r_p)
        s_e.set_val(Ecc)
        animation = FuncAnimation(fig1, update, blit=True, interval=10)
        animation.event_source.start()
    # Check if call was made from energy slider
    elif Eval != E:
        print("E slider was updated")
        E = Eval
        root1 = sp.lambdify(r, E + (G * M / r) - (L ** 2 / (2 * (r ** 2))) + ((G * M * (L ** 2)) / (r ** 3)))
        r_i = bisect(root1, a=IUCO, b=ISCO, disp=True)
        r_a = bisect(root1, a=ISCO, b=500000, disp=True)
        Ecc = get_e(r_i, r_a)
        r_p = r_i
        s_rp.set_val(r_p)
        s_e.set_val(Ecc)
        animation = FuncAnimation(fig1, update, blit=True, interval=10)
        animation.event_source.start()
    # When all slider values are equal to global values that agree with each other, all plotting information is updated
    else:
        print("Sliders are now updated")
        y0 = [r_i, rdot_i, phi_i]
        sol = solve_ivp(deriv, y0=y0, t_eval=Time, t_span=[t_i, t_f], rtol=1e-8, atol=1e-8, events=apoapsis_nt)
        get_H()
        orbit_trace.set_data(sol.y[2], sol.y[0])
        y = np.zeros_like(Ueff_Array)
        for i in range(y.size):
            y[i] = U_Eff_Func(Ueff_Array[i])
        a1.set_data(Ueff_Array, y)
        a2.set_data(r_i, U_Eff_Func(r_i))
        a3.set_data(Ueff_Array, energy_line(r_i))
        text_bot_E.set_val(E)
        text_bot_L.set_val(L)
        # Need to fix auto scaling for gravitational wave plots
        ax3.relim()
        ax3.autoscale_view()
        ax4.relim()
        ax4.autoscale_view()
        print(time.time() - start)


# Function to update Energy from text box
def submit_E(text):
    if float(text) != E:
        s_E.set_val(eval(text))


# Method to update angular momentum from text box
def submit_L(text):
    if float(text) != L:
        s_L.set_val(eval(text))


# Method to start the animation from the start button
def play_ani(id):
    animation.event_source.start()


# Method to stop the animation from the stop button
def pause_ani(id):
    animation.event_source.stop()


# calls updating functions for sliders and text box
s_E.on_changed(update_e)
s_L.on_changed(update_l)
s_rp.on_changed(update_rp)
s_e.on_changed(update_ecc)
text_bot_E.on_text_change(submit_E)
text_bot_L.on_text_change(submit_L)
pause.on_clicked(pause_ani)
play.on_clicked(play_ani)


# Method defining clear plots used for blit mapping in animation
def init_func():
    a2.set_data([], [])
    orbit_dot.set_data([], [])
    orbit.set_data([], [])
    h_xx.set_data([], [])
    h_xx_dot.set_data([], [])
    h_xy_dot.set_data([], [])
    h_xy.set_data([], [])
    return a2, orbit, orbit_dot, \
           h_xx, h_xx_dot, h_xy, h_xy_dot


# Method used to create the animation
def update(frame_number):
    start = time.time()
    a2.set_data(sol.y[0][frame_number], U_Eff_Func(r_i))
    orbit_dot.set_data(sol.y[2][frame_number], sol.y[0][frame_number])
    orbit.set_data(sol.y[2][:frame_number], sol.y[0][:frame_number])
    h_xx.set_data(sol.t[:frame_number], I_ddot[0][0][:frame_number])
    h_xx_dot.set_data(sol.t[frame_number], I_ddot[0][0][frame_number - 1])
    h_xy_dot.set_data(sol.t[frame_number], I_ddot[1][0][frame_number - 1])
    h_xy.set_data(sol.t[:frame_number], I_ddot[1][0][:frame_number])
    print(time.time() - start)
    if frame_number == (Step * t_f) - 5:
        animation.event_source.stop()
    return a2, orbit, orbit_dot, \
           h_xx, h_xx_dot, h_xy, h_xy_dot


animation = FuncAnimation(fig1, update, blit=True, interval=10)


plt.show()
