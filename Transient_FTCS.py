
import numpy as np
import matplotlib.pyplot as plt

DENSITY = 2.7889*1e3   # Kg/m^3
THERMAL_COND = 0.2473*1e3 # W/m-K
SPECIFIC_HEAT = 0.7639*1e3 # J/Kg-K
LENGTH = 0.12 # m
TIME = 2000  # secs

N1 = 1000
N2 = 2000
N3 = 3000
N4 = 5000
N5 = 10000
N6 = 20000

h = 3000 # W/m-K
NUM_PTS_X = 100
DELTA_X = LENGTH/NUM_PTS_X
LAMBDA = THERMAL_COND/(DENSITY*SPECIFIC_HEAT)
DELTA_T = DELTA_X ** 2 / (2* LAMBDA)
SIGMA = (LAMBDA * DELTA_T) / DELTA_X ** 2
NUM_PTS_T = int(TIME/ DELTA_T)
print(NUM_PTS_T)
BETA = (DELTA_X*h)/THERMAL_COND

T_amb = 298 # K

assert (N1 + N2 + N3 + N4 + N5 + N6) <= NUM_PTS_T, "Change NUM_PTS_T"

def initialTemp (x = 0):

    f= 673
    return f

def heatSource (x = 0 , t = 0):

    q = 0

    return q

def initialCondition(spaceLocation):

    up = np.full(len(spaceLocation), initialTemp())

    return up


def meshing():

    spaceMeshLocation = np.zeros(NUM_PTS_X+1)
    timeMeshLocations = np.zeros(NUM_PTS_T+1)

    for i in range(0,NUM_PTS_X+1):
        spaceMeshLocation[i] = i*DELTA_X

    for i in range(0,NUM_PTS_T+1):
        timeMeshLocations[i] = i*DELTA_T

    return spaceMeshLocation, timeMeshLocations


def plottemp(u,spaceLocation, time):

    plt.plot(spaceLocation,u[0,:], 'r')
    plt.plot(spaceLocation,u[1,:], 'g')
    plt.plot(spaceLocation, u[2, :], 'b')
    plt.plot(spaceLocation, u[3, :], '--r')
    plt.plot(spaceLocation, u[4, :], '--g')
    plt.plot(spaceLocation, u[5, :], '--b')

    plt.xlabel('Length (m)')
    plt.ylabel('Tempreature (K)')
    plt.title('Length ($x$)  vs. Tempreature ($u$)')
    label0 = 't ={0:0.2f} secs'.format(time[0])
    label1 = 't ={0:0.2f} secs'.format(time[1])
    label2 = 't ={0:0.2f} secs'.format(time[2])
    label3 = 't ={0:0.2f} secs'.format(time[3])
    label4 = 't ={0:0.2f} secs'.format(time[4])
    label5 = 't ={0:0.2f} secs'.format(time[5])

    plt.gca().legend((label0,label1,label2, label3, label4, label5), loc = 'upper right')
    plt.ylim(298,678)
    plt.xlim(0,0.12)
    plt.show()


def solver_FTCS():


    spaceLocation, timeLocation = meshing()
    up = initialCondition(spaceLocation)
    uc = np.zeros(len(spaceLocation))
    u = np.zeros((6, len(spaceLocation)))
    time = []
    u_ss = np.full(len(spaceLocation), T_amb)

    for t in range(1, len(timeLocation)):

        for i in range(1, len(spaceLocation)-1):

            uc[i] = up[i]+SIGMA*(up[i+1]-2*up[i]+up[i-1])+((heatSource()*DELTA_T)/(DENSITY*THERMAL_COND))

        uc[0] = (uc[1] + BETA*T_amb)/(1+BETA)
        uc[len(spaceLocation)-1] = (uc[len(spaceLocation)-2]+BETA*T_amb)/(1+BETA)
        up = uc

        if timeLocation[t] == N1*DELTA_T:
            u[0,:] = up
            time.append(N1*DELTA_T)

        elif timeLocation[t] ==  N2*DELTA_T:
            u[1,:] = up
            time.append(N2*DELTA_T)

        elif timeLocation[t] ==  N3*DELTA_T:
            u[2,:] = up
            time.append(N3 * DELTA_T)

        elif timeLocation[t] == N4 * DELTA_T:
            u[3, :] = up
            time.append(N4 * DELTA_T)

        elif timeLocation[t] == N5 * DELTA_T:
            u[4, :] = up
            time.append(N5 * DELTA_T)

        elif timeLocation[t] == N6 * DELTA_T:
            u[5, :] = up
            time.append(N6 * DELTA_T)

        u_check = np.subtract(up, u_ss)

        if (u_check < 0.001).all():
            print("Steady state reached at {0:0.2f} secs".format(timeLocation[t]))
            break

    plottemp(u, spaceLocation, time)




solver_FTCS()










