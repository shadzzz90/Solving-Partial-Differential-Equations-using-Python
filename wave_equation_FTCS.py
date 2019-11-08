
import numpy as np
import matplotlib.pyplot as plt
from math import *


# BETA = 0.02 # /sec

BETA = 0 # /sec

NUM_PTS_X = 100
TIME = 120# secs
L = 0.6 # m
c = 1 # m/s
a = 0.1
DELTA_X = L/NUM_PTS_X

# using the CFL stability condition

CFL= 0.99
DELTA_T = (CFL*DELTA_X)/c

NUM_PTS_T = int(TIME/ DELTA_T)

# calculating constants

ALPHA = (BETA*DELTA_T)/2
GAMMA = (c*DELTA_T/DELTA_X)**2

A = (ALPHA-1)/(ALPHA+1)
B = GAMMA/(ALPHA+1)
C = (2*(1-GAMMA))/(ALPHA+1)

#
u_beta_0 =  np.zeros(NUM_PTS_T+1)
u_beta_002 =  np.zeros(NUM_PTS_T+1)


def meshing():

    spaceMeshLocation = np.zeros(NUM_PTS_X+1)
    timeMeshLocations = np.zeros(NUM_PTS_T+1)

    for i in range(0,NUM_PTS_X+1):
        spaceMeshLocation[i] = i*DELTA_X

    for i in range(0,NUM_PTS_T+1):
        timeMeshLocations[i] = i*DELTA_T

    return spaceMeshLocation, timeMeshLocations


def initialCondition(spaceLocations):

    u_initial = np.zeros(len(spaceLocations))

    for i in range(0,len(spaceLocations)):

        if spaceLocations[i] <= 0.3:

            u_initial[i] = (2*a*spaceLocations[i])/L

        else:

            u_initial[i] = 2 * a * (1-spaceLocations[i] / L)


    return u_initial

def plotdisp(u,spaceLocation,timeLocations):

    plt.plot(spaceLocation,u[0,:], 'y')
    plt.plot(spaceLocation,u[1,:], 'r')
    plt.plot(spaceLocation, u[2, :], 'g')
    plt.plot(spaceLocation, u[3, :], 'b')
    plt.plot(spaceLocation, u[4, :], '--r')
    plt.plot(spaceLocation, u[5, :], '--g')
    plt.plot(spaceLocation, u[6, :], '--b')

    plt.xlabel('Length (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Length ($x$)  vs. Displacement ($u$) for CFL = 0.99')
    # plt.title('Length ($x$)  vs. Displacement ($u$) for Beta = {0:0.2f}'.format(BETA))

    label0 = 't = 0 secs'
    label1 = 't = 1 secs'
    label2 = 't = 5 secs'
    label3 = 't = 10 secs'
    label4 = 't = 25 secs'
    label5 = 't = 50 secs'
    label6 = 't = 100 secs'



    plt.gca().legend((label0,label1,label2,label3,label4,label5, label6), loc = 'upper right')
    # plt.ylim(0,0.12)
    plt.xlim(0,0.6)
    plt.show()

    u_beta_0 = np.load('u_beta_0_0.3.npy')
    u_beta_002 = np.load('u_beta_002_0.3.npy')


    plt.plot(timeLocations,u_beta_0)
    plt.xlabel('t (sec)')
    plt.ylabel('u (m)')
    plt.title('Time (t)  vs. Displacement (u) at x = 0.3 m  (Beta = 0)')
    plt.show()
    plt.plot(timeLocations, u_beta_002)
    plt.xlabel('t(sec)')
    plt.ylabel('u (m)')
    plt.title('Time (t)  vs. Displacement (u) at x = 0.3 m (Beta = 0.02)$')
    plt.show()



def solver(u_initial, spaceLocations, timeLocations):



        index = np.argwhere(spaceLocations == 0.3)


        u_prev = np.zeros(len(spaceLocations))
        u_prev_prev = np.zeros(len(spaceLocations))

        u_prev[:] = u_initial[:]
        u_prev_prev[:] = u_initial[:]

        u_curr = np.zeros(len(spaceLocations))

        u = np.zeros((7, len(spaceLocations)))
        u[0,:] = u_initial

        u_beta_0[0] = u_initial[index[0][0]]
        u_beta_002[0]= u_initial[index[0][0]]




        for m in range(1, len(timeLocations)):

            u_curr[0] = 0
            u_curr[len(spaceLocations)-1] = 0

            for i in range(1,len(spaceLocations)-1):
                u_curr[i] = (A*u_prev_prev[i])+(B*u_prev[i+1])+(C*u_prev[i])+(B*u_prev[i-1])


            # if BETA == 0:
            #
            #     u_beta_0[m] = u_curr[index[0][0]]
            #
            # elif BETA == 0.02:
            #
            #     u_beta_002[m] = u_curr[index[0][0]]



            if round(round(timeLocations[m],2),1) == 1.00:
                u[1, :] = u_curr[:]


            elif round(round(timeLocations[m],2),1)  == 5.00:
                u[2, :] = u_curr[:]


            elif round(round(timeLocations[m],2),1) == 10.00:
                u[3, :] = u_curr[:]


            elif round(round(timeLocations[m],2),1) == 25.00:
                u[4, :] = u_curr[:]


            elif round(round(timeLocations[m],2),1)== 50.00:
                u[5, :] = u_curr[:]


            elif round(round(timeLocations[m],2),1) == 100.00:
                u[6, :] = u_curr[:]


            u_prev_prev[:] = u_prev[:]
            u_prev[:] = u_curr[:]


        # if BETA == 0:
        #
        #     np.save('u_beta_0_0.3', u_beta_0)
        #
        # elif BETA == 0.02:
        #
        #     np.save('u_beta_002_0.3', u_beta_002)



        return u




spaceLocations, timeLocations = meshing()

u_initial = initialCondition(spaceLocations)

u= solver(u_initial,spaceLocations, timeLocations)

plotdisp(u, spaceLocations, timeLocations)