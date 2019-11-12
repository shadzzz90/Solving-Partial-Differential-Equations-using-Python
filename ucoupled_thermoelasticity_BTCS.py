
import numpy as np
import matplotlib.pyplot as plt
# from math import *

NUM_PTS_R_TC = 5
NUM_PTS_R_STL = 5
NUM_PTS_T = 10
TIME = 2 # secs
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

betaTC = (DIFF_TC * DELTA_T) / (DELTA_R_TC)
betaSTL = (DIFF_STL * DELTA_T) / (DELTA_R_STL)

GAMMA = (K_STL*DELTA_R_TC)/(K_TC*DELTA_R_STL)

LAMBDA = (h*DELTA_R_STL)/K_STL


# def matrix_builder_STL():
#
#
#     n_STL = len(spaceLocations_STL)-2
#
#     A_STL = np.zeros((n_STL,n_STL))
#
#     for i in range(0,n_STL-1):
#         A_STL[i][i] = (1+2*betaSTL)
#
#     temp = (betaSTL/(1+LAMBDA))*(1+(1/n_STL))
#
#     A_STL[n_STL-1][n_STL-1] = 1+2*betaSTL-temp
#
#
#     for i in range(0,n_STL-1):
#         A_STL[i][i+1] = -betaSTL*(1+(1/(i+1)))
#         A_STL[i+1][i] = -betaSTL*(1-(1/(i+2)))
#
#
#     b_STL = np.full(n_STL, T_inital)
#
#     b_STL[n_STL-1] = b_STL[n_STL-1] + ((LAMBDA*betaSTL)/(1+LAMBDA))*(1+(1/n_STL))*T_amb
#
#
#     return A_STL, b_STL
#
# def matrix_builder_TC(T_first_STL):
#
#     n_TC = len(spaceLocations_TC) - 2
#
#     A_TC = np.zeros((n_TC, n_TC))
#
#     for i in range(0, n_TC - 1):
#         A_TC[i][i] = (1 + 2 * betaTC)
#
#     temp = (betaTC / (1 + GAMMA)) * (1 + (1 / n_TC))
#
#     A_TC[n_TC - 1][n_TC - 1] = 1 + 2 * betaTC - temp
#
#     for i in range(0,n_TC-1):
#         A_TC[i][i+1] = -betaTC*(1+(1/(i+1)))
#         A_TC[i+1][i] = -betaTC*(1-(1/(i+2)))
#
#     b_TC = np.full(n_TC, T_inital)
#
#     b_TC[n_TC - 1] = b_TC[n_TC - 1] + ((GAMMA * betaTC) / (1 + GAMMA)) * (1 + (1 / n_TC)) * T_first_STL
#
#     return A_TC, b_TC

def matrix_builder(spaceLocations_STL, spaceLocations_TC):

    n = (len(spaceLocations_STL)-2)+(len(spaceLocations_TC)-2)
    A_TC = np.zeros(len(spaceLocations_TC)-3)
    A_STL = np.zeros(len(spaceLocations_STL)-2)
    C_TC = np.zeros(len(spaceLocations_TC)-2)
    C_STL = np.zeros(len(spaceLocations_STL)-3)
    B_TC = np.zeros(len(spaceLocations_TC)-2)
    B_STL = np.zeros(len(spaceLocations_STL)-2)
    A = np.zeros((n,n))


    for i in range(0, len(A_STL)):

        A_STL[i] = - betaSTL*(1/DELTA_R_STL-1/(inDia_STL+(i+1)*DELTA_R_STL))

    for i in range(0, len(A_TC)):
        A_TC[i] = - betaTC * (1 / DELTA_R_TC - 1 / (inDia_TC + (i+1) * DELTA_R_TC))

    A_STL[0] = A_STL[0]/(1+GAMMA)

    for i in range(0, len(C_TC)):

        C_TC[i] = - betaTC*(1/DELTA_R_TC+1/(inDia_TC+(i+1)*DELTA_R_TC))

    for i in range(0, len(C_STL)):
        C_STL[i] = - betaSTL * (1 / DELTA_R_STL + 1 / (inDia_STL + (i+1) * DELTA_R_STL))

    C_TC[-1] = (GAMMA*C_TC[-1])/(1+GAMMA)

    for i in range(0,len(B_TC)):
        B_TC[i] = 1+(2*betaTC)/DELTA_R_TC

    for i in range(0,len(B_STL)):
        B_STL[i] = 1+(2*betaSTL)/DELTA_R_STL

    B_TC[-1] = B_TC[-1]+C_TC[-1]/(1+GAMMA)
    B_STL[0] = B_STL[0]+(A_STL[0]*GAMMA)/(1+GAMMA)
    B_STL[-1] = B_STL[-1]+ C_STL[-1]/(1+LAMBDA)


    for i in range(0,len(spaceLocations_TC)-2):

        A[i][i] = B_TC[i]

    for i in range((len(spaceLocations_TC)-2), n):

        A[i][i] = B_STL[i-(len(spaceLocations_TC)-2)]


    for i in range(0,len(spaceLocations_TC)-2):

        A[i][i+1] = C_TC[i]

    for i in range(0, len(spaceLocations_TC) - 3):
        A[i+1][i] = A_TC[i]


    for i in range((len(spaceLocations_TC)-2), n-1):

        A[i][i+1] = C_STL[i-(len(spaceLocations_TC)-2)]


    for i in range((len(spaceLocations_TC)-3), n-1):

        A[i+1][i] = A_STL[i-(len(spaceLocations_TC)-3)]


    return A





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

# def mainSolver(timeLocations):
#
#
#
#     A_STL, b_STL = matrix_builder_STL()
#
#     x_STL = SR_solver(A_STL, b_STL)
#
#     A_TC, b_TC = matrix_builder_TC(x_STL[0])
#
#     x_TC = SR_solver(A_TC, b_TC)
#
#     for i in range(1,len(timeLocations)):
#
#         x_curr_STL = SR_solver(A_STL,x_STL)
#
#         x_curr_TC =  SR_solver(A_TC,x_TC)
#
#         x_STL[:] = x_curr_STL[:]
#         x_TC[:] = x_curr_TC[:]
#     # x_STL2 = SR_solver(A_STL, x_STL)
#     # x_TC2 = SR_solver(A_TC,x_TC)
#
#
#     return  x_TC, x_STL

# def test(timeLocations):
#
#     T = np.zeros(len(timeLocations))
#
#     for i in range(0, len(timeLocations)):
#
#         T[i] = 293+553*(1 - exp(-10*i*DELTA_T))
#
#
#     plt.plot(timeLocations, T)
#     plt.show()
#





spaceLocations_TC, spaceLocations_STL, timeLocations = meshing()

print(matrix_builder(spaceLocations_STL,spaceLocations_TC))

#
# print(mainSolver(timeLocations))

# test(timeLocations)