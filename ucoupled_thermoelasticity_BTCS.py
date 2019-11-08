
import numpy as np
import matplotlib.pyplot as plt
from math import *

NUM_PTS_R_TC = 10
NUM_PTS_R_STL =10
NUM_PTS_T = 10
TIME = 5  # secs
T_inital = 293 # K
T_amb = 303 # K

h = 200 # W/m^2-K
RHO_TC = 15.88*1e3 # Tungsten Carbide Kg/m^3
Cp_TC = 0.292 # KJ/Kg-K
K_TC = 88 # W/m-K

RHO_STL = 8000 # Stainless Steel AISI 304 Kg/m^3
Cp_STL = 0.4 # KJ/Kg-K
K_STL = 13.8 # W/m-K

DIFF_TC = K_TC/(RHO_TC*Cp_TC)
DIFF_STL = K_STL/(RHO_STL*Cp_STL)

inDia_TC = 70*1e-3 # m
outDia_TC = 80*1e-3 # m

inDia_STL = 80*1e-3 # m
outDia_STL = 200*1e-3 # m

thickness_TC = outDia_TC - inDia_TC
thickness_STL = outDia_STL -inDia_STL

DELTA_R_TC = thickness_TC/NUM_PTS_R_TC
DELTA_R_STL = thickness_STL/NUM_PTS_R_STL

DELTA_T = TIME / NUM_PTS_T

betaTC = (DIFF_TC * DELTA_T) / (DELTA_R_TC**2)
betaSTL = (DIFF_STL * DELTA_T) / (DELTA_R_STL**2)

GAMMA = (K_STL*DELTA_R_TC)/(K_TC*DELTA_R_STL)

LAMBDA = (h*DELTA_R_STL)/K_STL


def matrix_builder_STL():


    n_STL = len(spaceLocations_STL)-2

    A_STL = np.zeros((n_STL,n_STL))

    for i in range(0,n_STL-1):
        A_STL[i][i] = (1+2*betaSTL)

    temp = (betaSTL/(1+LAMBDA))*(1+(1/n_STL))

    A_STL[n_STL-1][n_STL-1] = 1+2*betaSTL-temp


    for i in range(0,n_STL-1):
        A_STL[i][i+1] = -betaSTL*(1+(1/(i+1)))
        A_STL[i+1][i] = -betaSTL*(1-(1/(i+2)))


    b_STL = np.full(n_STL, T_inital)

    b_STL[n_STL-1] = b_STL[n_STL-1] + ((LAMBDA*betaSTL)/(1+LAMBDA))*(1+(1/n_STL))*T_amb


    return A_STL, b_STL

def matrix_builder_TC(T_last):

    n_TC = len(spaceLocations_TC) - 2

    A_TC = np.zeros((n_TC, n_TC))

    for i in range(0, n_TC - 1):
        A_TC[i][i] = (1 + 2 * betaTC)

    temp = (betaTC / (1 + GAMMA)) * (1 + (1 / n_TC))

    A_TC[n_TC - 1][n_TC - 1] = 1 + 2 * betaTC - temp

    for i in range(0,n_TC-1):
        A_TC[i][i+1] = -betaTC*(1+(1/(i+1)))
        A_TC[i+1][i] = -betaTC*(1-(1/(i+2)))

    b_TC = np.full(n_TC, T_inital)

    b_TC[n_TC - 1] = b_TC[n_TC - 1] + ((GAMMA * betaTC) / (1 + GAMMA)) * (1 + (1 / n_TC)) * T_last

    return A_TC, b_TC




def meshing():

    spaceMeshLocations_TC = np.zeros(NUM_PTS_R_TC+1)
    spaceMeshLocations_STL = np.zeros(NUM_PTS_R_STL + 1)
    timeMeshLocations = np.zeros(NUM_PTS_T+1)

    for i in range(0,NUM_PTS_R_TC+1):
        spaceMeshLocations_TC[i] = i*DELTA_R_TC

    for i in range(0, NUM_PTS_R_STL + 1):
        spaceMeshLocations_STL[i] = i * DELTA_R_STL

    for i in range(0,NUM_PTS_T+1):
        timeMeshLocations[i] = i*DELTA_T

    return spaceMeshLocations_TC, spaceMeshLocations_STL, timeMeshLocations


def SR_solver(A,b):
    """Sucessive Relaxation Solver"""

    n=len(b)
    relax_factor = 1
    epsilon = 1e-8

    # A = np.random.randn(n,n)+4*np.eye(n)
    # b = np.random.randn(n)
    #
    # x = np.linalg.solve(A,b.T)

    x_inital = np.zeros(n)
    x_next = np.zeros(n)
    iteration_number =1


    while True:

        sum =0
        iteration_number = iteration_number+1

        for j in range(0,n):
            sum = sum + A[0][j]*x_inital[j]

        x_next[0] = x_inital[0]+(b[0]-sum)*relax_factor/A[0][0]  # First point

        sum = 0

        for j in range(1,n):
            sum = sum + A[1][j]*x_inital[j]

        x_next[1] = x_inital[1] + (relax_factor * (b[1] -A[1][0]*x_next[0]-sum)) / A[1][1]

        for i in range(2,n):

            sum1 = 0

            for j in range(0,i):

                sum1 = sum1 + A[i][j]*x_next[j]

            sum2 = 0

            for j in range(i,n):
                sum2 = sum2+A[i][j]*x_inital[j]

            x_next[i] = x_inital[i]+(relax_factor*(b[i]-sum1-sum2))/A[i][i]

        x_inital = x_next

        if np.linalg.norm(np.matmul(A,x_next)-b) <= epsilon*np.linalg.norm(b):
            print('Total Iteration : {}'.format(iteration_number))
            break

    return x_next



spaceLocations_TC, spaceLocations_STL, timeLocations = meshing()

A_STL, b_STL = matrix_builder_STL()

x_temp = SR_solver(A_STL, b_STL)

A_TC, b_TC = matrix_builder_TC(x_temp[0])

print(SR_solver(A_TC, b_TC), x_temp)

