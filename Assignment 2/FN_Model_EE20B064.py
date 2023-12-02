#Importing required Libraries
import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 0.5
b = 0.1
r = 0.1

# Defining the time array for the v, w vs time plots
tmax = 100
dt = 0.1
t = np.arange(0, tmax, dt)

# Defining the dynamics of the system and the nullclines
def f(v):
    return v*(a-v)*(v-1)
def dvdt(v, w, Iext):
    return f(v) - w + Iext
def dwdt(v, w):
    return b*v - r*w
def v_nc(v):
    return f(v) + Iext
def w_nc(v):
    return b*v/r

# Phase portrait plotter
def phase_portrait(Iext, init_conds, traj, ylims):
    # x and y axes limits for the phase portrait
    xlim1 = -0.5
    xlim2 = 1.5
    ylim1 = ylims[0]
    ylim2 = ylims[1]
    # Defining points where the arrows will be plotted
    v = np.linspace(xlim1, xlim2, 100)
    w = np.linspace(xlim1, xlim2, 100)
    # Making the meshed grid and defining the derivatives (ie. arrows)
    v_grid, w_grid = np.meshgrid(v, w)
    vdot = dvdt(v_grid, w_grid, Iext)
    wdot = dwdt(v_grid, w_grid)
    # If initial points given, draw the trajectory
    if traj:
        for i in range(len(init_conds)):
            plt.streamplot(v_grid, w_grid, vdot, wdot, start_points = [init_conds[i]], integration_direction='forward')
    # Else, plot the phase portrait
    else:
        plt.streamplot(v_grid, w_grid, vdot, wdot, density=2)
    # Plotting the v and w nullclines for reference
    plt.plot(v, v_nc(v), label='v nullcline')
    plt.plot(v, w_nc(v), label='w nullcline')
    plt.ylim(ylim1, ylim2)
    plt.xlabel('v')
    plt.ylabel('w')
    plt.legend()
    plt.title('Phase Plot for $I_{ext}$ =' + str(round(Iext, 2)))
    plt.grid()
    plt.show()


# Time plot
def vs_time(Iext, v0, w0, plt_title):
    # Defining v, w arrays to store the variables over time
    v = np.zeros_like(t)
    w = np.zeros_like(t)
    v[0] = v0
    w[0] = w0
    # Updating the variables (ie. Forward Euler Integration)
    for time in range(len(t)-1):
        v[time+1] = v[time] + (dvdt(v[time], w[time], Iext))*dt
        w[time+1] = w[time] + (dwdt(v[time], w[time]))*dt

    plt.xlabel('Time (t)')
    plt.ylabel('Voltage (V)')
    plt.plot(t, v, label='v(t)')
    plt.plot(t, w, label='w(t)')
    plt.title(plt_title)
    plt.legend()
    plt.grid()
    plt.show()

# [Time plots: (i) V(0) = 0.4 < a and (ii) V(0) = 0.6 > a]
# Case 1 : Iext = 0
Iext = 0
phase_portrait(Iext, [], False, [-0.5, 0.5])
phase_portrait(Iext, [[0,0], [0.4,0], [0.6, 0], [1, 0]], True, [-0.5, 0.5])
vs_time(Iext, 0.4, 0, "$I_{ext}$ = 0, v(0) = 0.4, w(0) = 0")
vs_time(Iext, 0.6, 0, "$I_{ext}$ = 0, v(0) = 0.6, w(0) = 0")

# Case 2 : Iext = 0.6
Iext = 0.6
phase_portrait(Iext, [], False, [0, 1])
phase_portrait(Iext, [[-0.25, 0.2], [0.2, 1], [0.8, 0], [1.3, 0.7], [0.6, 0.6]], True, [0, 1])
vs_time(Iext, 0.4, 0, "$I_{ext}$ = 0.6, v(0) = 0.4, w(0) = 0")
vs_time(Iext, 0.6, 0, "$I_{ext}$ = 0.6, v(0) = 0.6, w(0) = 0")

# Case 3 : Iext = 1.2
Iext = 1.2
phase_portrait(Iext, [], False, [0.7, 1.5])
phase_portrait(Iext, [[0.3,0.8], [1.5, 1], [1.2, 1.4], [0, 1.3]], True, [0.7, 1.5])
vs_time(Iext, 0.4, 0, "$I_{ext}$ = 1.2, v(0) = 0.4, w(0) = 0")
vs_time(Iext, 0.6, 0, "$I_{ext}$ = 1.2, v(0) = 0.6, w(0) = 0")

# Case 4 : Iext = 0.02, b,r = 0.01, 0.8
Iext = 0.01
b = 0.02
r = 0.6
phase_portrait(Iext, [], False, [-0.25, 0.25])
phase_portrait(Iext, [[0,-0.3],[-0.15, 0.3], [0.4,-0.07],[0.55,-0.05],[1,-0.3], [1.2, 0.3]], True, [-0.25, 0.25])
vs_time(Iext, 0.4, 0, "$I_{ext}$ = 0.01, v(0) = 0.4, w(0) = 0")
vs_time(Iext, 0.6, 0, "$I_{ext}$ = 0.01, v(0) = 0.6, w(0) = 0")