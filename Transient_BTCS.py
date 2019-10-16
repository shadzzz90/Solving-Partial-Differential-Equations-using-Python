
import numpy as np
import matplotlib.pyplot as plt
from math import *

DENSITY = 7800   # Kg/m^3
THERMAL_COND = 55 # W/m-K
SPECIFIC_HEAT = 450 # J/Kg-K
LENGTH = 0.6 # m
TIME =1000# secs

h = 1e5# W/m-K
NUM_PTS_X = 100
NUM_PTS_T = 1000
DELTA_X = LENGTH/NUM_PTS_X
DELTA_T = TIME/NUM_PTS_T

SIGMA = (THERMAL_COND*DELTA_T)/(DENSITY*SPECIFIC_HEAT*DELTA_X**2)
BETA = DELTA_T/(DENSITY*SPECIFIC_HEAT)
T_amb = 293.00 # K

ALPHA = (h*DELTA_X)/THERMAL_COND



def initialTemp (x = 0):

    f= 293.00
    return f

def heatSource (x):

    q = 2e7*exp(-1000*(x-0.2)**2)

    return q


def tridiag_solver(a, b, c, f):
    """Implementing Thomas Algorithm"""

    n = len(f)
    a_dash = np.zeros(n)
    c_dash = np.zeros(n-1)
    y = np.zeros(n)
    u = np.zeros(n)
    a_dash[0] = a[0]

    for i in range(1,n):
        a_dash[i] = a[i] - (c[i-1]*b[i-1]/a_dash[i-1])

    for i in range(0,n-1):
        c_dash[i] = c[i]/a_dash[i]

    y[0] = f[0]

    for i in range(1,n):
        y[i] = f[i] - c_dash[i-1]*y[i-1]

    u[n-1] = y[n-1]/a_dash[n-1]

    for i in range(n-2,-1,-1):
        u[i] = (y[i]-b[i]*u[i+1])/a_dash[i]

    return u




def meshing():

    spaceMeshLocations = np.zeros(NUM_PTS_X+1)
    timeMeshLocations = np.zeros(NUM_PTS_T+1)

    for i in range(0,NUM_PTS_X+1):
        spaceMeshLocations[i] = i*DELTA_X

    for i in range(0,NUM_PTS_T+1):
        timeMeshLocations[i] = i*DELTA_T

    return spaceMeshLocations, timeMeshLocations


def plottemp(u,spaceLocation, u_time, timeLocation, flux):

    plt.figure(1)
    plt.plot(spaceLocation,u[0], 'r')
    plt.plot(spaceLocation,u[1], 'g')
    plt.plot(spaceLocation, u[2], 'b')
    plt.plot(spaceLocation, u[3], '--r')
    plt.plot(spaceLocation, u[4], '--g')

    plt.xlabel('Length (m)')
    plt.ylabel('Tempreature (K)')
    plt.title('Length ($x$)  vs. Tempreature ($u$)')
    label0 = 't = 0 secs'
    label1 = 't = 10 secs'
    label2 = 't = 100 secs'
    label3 = 't = 500 secs'
    label4 = 't = 1000 secs'
    label5 = 'x = 0.2 m'

    plt.gca().legend((label0,label1,label2, label3, label4), loc = 'upper right')
    plt.ylim(280,1800)
    plt.xlim(0,0.6)
    plt.show()

    plt.figure(2)
    plt.plot(timeLocation,u_time,'-r')
    plt.xlabel('Time (secs)')
    plt.ylabel('Tempreature (K)')
    plt.title('Time ($t$)  vs. Tempreature ($u$) at $x = 0.2$ m')
    plt.xlim(0, 1000)

    plt.show()

    plt.figure(3)
    plt.plot(spaceLocation, flux[0,:], 'r')
    plt.plot(spaceLocation, flux[1,:], 'g')
    plt.plot(spaceLocation, flux[2,:], 'b')
    plt.plot(spaceLocation, flux[3,:], '--r')
    plt.plot(spaceLocation, flux[4,:], '--g')

    plt.xlabel('Length (m)')
    plt.ylabel('Flux ($W/m^2$)')
    plt.title('Flux  vs. Tempreature')
    label0 = 't = 0 secs'
    label1 = 't = 10 secs'
    label2 = 't = 100 secs'
    label3 = 't = 500 secs'
    label4 = 't = 1000 secs'


    plt.gca().legend((label0, label1, label2, label3, label4), loc='upper right')
    # plt.ylim(280, 1800)
    plt.xlim(0, 0.6)
    plt.show()

def flux(u_hist, spaceMeshLocations):
    """Calculating the flux"""

    FLUX = np.zeros((len(u_time),len(spaceMeshLocations)))
    n = len(spaceMeshLocations)
    flux = np.zeros(n)

    for i in range(0,len(u_hist)):

        u_temp = u_hist[i]

        flux[0] = - THERMAL_COND*((u_temp[1]-(1-ALPHA)*u_temp[0]-ALPHA*T_amb)/(2*DELTA_X))

        for j in range(1,n-1):

            flux[j] = - THERMAL_COND*((u_temp[j+1]-u_temp[j-1])/(2*DELTA_X))

        flux[n-1] = - THERMAL_COND*(((1-ALPHA)*u_temp[n-1] + ALPHA*T_amb-u_temp[n-2])/(2*DELTA_X))

        FLUX[i,:] = flux

    return FLUX



def solver_BTCS(spaceMeshLocations, timeMeshLocations):

    """Solving using BTCS"""

    n = len(spaceMeshLocations)

    b = np.full(n-1,-SIGMA)
    c = np.full(n-1, -SIGMA)
    a = np.zeros(n)

    a[0] = (1+SIGMA+SIGMA*ALPHA)
    a[n-1] = (1+SIGMA+SIGMA*ALPHA)
    a[1:n-1] = (1+2*SIGMA)

    heatSource_term = np.zeros(n)

    for i in range(0, len(spaceMeshLocations)):
        heatSource_term[i] = heatSource(spaceMeshLocations[i])*BETA


    BC_term = np.zeros(n)

    BC_term[0] = SIGMA*ALPHA*T_amb
    BC_term[n-1] = SIGMA*ALPHA*T_amb

    u_hist = []

    const_value = np.full(n,0.2)

    difference = spaceMeshLocations-const_value

    loc_index = np.argwhere(difference >= 0)
    index = loc_index[0]

    u_initial = np.full(n,initialTemp())
    u_initial[0] = (u_initial[1]+ALPHA*T_amb)/(1+ALPHA)
    u_initial[n-1] = (u_initial[n-2] + ALPHA * T_amb) / (1 + ALPHA)

    u_time = np.zeros(len(timeMeshLocations))

    u_time[0] = u_initial[index]

    u_hist.append(u_initial)

    u_curr = u_initial+heatSource_term+BC_term

    for i in range(1,len(timeMeshLocations)):

        u_next = tridiag_solver(a,b,c, u_curr)
        u_time[i] = u_next[index]

        u_curr = u_next+heatSource_term+BC_term

        if timeMeshLocations[i] == 10:

            u_hist.append(u_next)

        elif timeMeshLocations[i] == 100:

            u_hist.append(u_next)

        elif timeMeshLocations[i] == 500:

            u_hist.append(u_next)

        elif timeMeshLocations[i] == 1000:

            u_hist.append(u_next)



    return u_hist, u_time


spaceMeshLocations, timeMeshLocations = meshing()

u_history, u_time = solver_BTCS(spaceMeshLocations, timeMeshLocations)

heat_flux = flux(u_history,spaceMeshLocations)

plottemp(u_history,spaceMeshLocations, u_time, timeMeshLocations, heat_flux)










