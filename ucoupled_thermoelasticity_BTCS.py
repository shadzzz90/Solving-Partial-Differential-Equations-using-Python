
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter
from math import *
from matplotlib import colors as mcolors
from matplotlib import cm

NUM_PTS_R_TC = 50
NUM_PTS_R_STL = 50
NUM_PTS_T = 100
TIME = 100 # secs
T_inital = 293 # K
T_amb = 303 # K
T_final = 573 # K

h = 200 # W/m^2-K
RHO_TC = 15.25*1e3 # Tungsten Carbide Kg/m^3
Cp_TC = 0.184 # KJ/Kg-K
K_TC = 28 # W/m-K

RHO_STL = 7833 # Steel 0.5% C Kg/m^3
Cp_STL = 0.465 # KJ/Kg-K
K_STL = 54 # W/m-K

DIFF_TC = K_TC/(RHO_TC*Cp_TC)
DIFF_STL = K_STL/(RHO_STL*Cp_STL)

inR_TC = 70*1e-3 # m
outR_TC = 80*1e-3 # m

inR_STL = 80*1e-3 # m
outR_STL = 200*1e-3 # m

thickness_TC = outR_TC - inR_TC
thickness_STL = outR_STL -inR_STL

DELTA_R_TC = thickness_TC/NUM_PTS_R_TC
DELTA_R_STL = thickness_STL/NUM_PTS_R_STL

DELTA_T = TIME / NUM_PTS_T

betaTC = (DIFF_TC * DELTA_T) / DELTA_R_TC
betaSTL = (DIFF_STL * DELTA_T) / DELTA_R_STL

GAMMA = (K_STL*DELTA_R_TC)/(K_TC*DELTA_R_STL)

LAMBDA = (h*DELTA_R_STL)/K_STL

nu_TC = 0.2
nu_STL = 0.3

Coeff_Thermal_Exp_TC = 4.5 * 1e-6 # /K
Coeff_Thermal_Exp_STL = 1.24 * 1e-5 # /K

E_TC = 620*1e9 # GPa
E_STL = 210*1e9 # GPa

Phi_TC = (1-nu_TC)/(1+nu_TC)
Phi_STL = (1-nu_STL)/(1+nu_STL)

Xi_TC = 2*nu_TC/(1+nu_TC)
Xi_STL = 2*nu_STL/(1+nu_STL)

Psi_TC = (1+nu_TC)*Coeff_Thermal_Exp_TC/(1-nu_TC)
Psi_STL = (1+nu_STL)*Coeff_Thermal_Exp_STL/(1-nu_STL)

Eta_TC = E_TC*(1-nu_TC)/(DELTA_R_TC*(1+nu_TC)*(1-2*nu_TC))
Eta_STL = E_STL*(1-nu_STL)/(DELTA_R_STL*(1+nu_STL)*(1-2*nu_STL))

Mu_TC = 2*E_TC*nu_TC/(outR_TC*(1+nu_TC)*(1-2*nu_TC))
Mu_STL = 2*E_STL*nu_STL/(inR_STL*(1+nu_STL)*(1-2*nu_STL))

Kappa_TC = E_TC*DIFF_TC/(1-2*nu_TC)
Kappa_STL = E_STL*DIFF_STL/(1-2*nu_STL)



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
#
# def matrix_builder(spaceLocations_STL, spaceLocations_TC):
#
#     n = (len(spaceLocations_STL)-2)+(len(spaceLocations_TC)-2)
#     A_TC = np.zeros(len(spaceLocations_TC)-3)
#     A_STL = np.zeros(len(spaceLocations_STL)-2)
#     C_TC = np.zeros(len(spaceLocations_TC)-2)
#     C_STL = np.zeros(len(spaceLocations_STL)-3)
#     B_TC = np.zeros(len(spaceLocations_TC)-2)
#     B_STL = np.zeros(len(spaceLocations_STL)-2)
#     A = np.zeros((n,n))
#
#
#     for i in range(0, len(A_STL)):
#
#         A_STL[i] = - betaSTL*(1/DELTA_R_STL-1/(inDia_STL+(i+1)*DELTA_R_STL))
#
#     for i in range(0, len(A_TC)):
#         A_TC[i] = - betaTC * (1 / DELTA_R_TC - 1 / (inDia_TC + (i+1) * DELTA_R_TC))
#
#     A_STL[0] = A_STL[0]/(1+GAMMA)
#
#     for i in range(0, len(C_TC)):
#
#         C_TC[i] = - betaTC*(1/DELTA_R_TC+1/(inDia_TC+(i+1)*DELTA_R_TC))
#
#     for i in range(0, len(C_STL)):
#         C_STL[i] = - betaSTL * (1 / DELTA_R_STL + 1 / (inDia_STL + (i+1) * DELTA_R_STL))
#
#     C_TC[-1] = (GAMMA*C_TC[-1])/(1+GAMMA)
#
#     for i in range(0,len(B_TC)):
#         B_TC[i] = 1+(2*betaTC)/DELTA_R_TC
#
#     for i in range(0,len(B_STL)):
#         B_STL[i] = 1+(2*betaSTL)/DELTA_R_STL
#
#     B_TC[-1] = B_TC[-1]+C_TC[-1]/(1+GAMMA)
#     B_STL[0] = B_STL[0]+(A_STL[0]*GAMMA)/(1+GAMMA)
#     B_STL[-1] = B_STL[-1]+ C_STL[-1]/(1+LAMBDA)
#
#
#     for i in range(0,len(spaceLocations_TC)-2):
#
#         A[i][i] = B_TC[i]
#
#     for i in range((len(spaceLocations_TC)-2), n):
#
#         A[i][i] = B_STL[i-(len(spaceLocations_TC)-2)]
#
#
#     for i in range(0,len(spaceLocations_TC)-2):
#
#         A[i][i+1] = C_TC[i]
#
#     for i in range(0, len(spaceLocations_TC) - 3):
#         A[i+1][i] = A_TC[i]
#
#
#     for i in range((len(spaceLocations_TC)-2), n-1):
#
#         A[i][i+1] = C_STL[i-(len(spaceLocations_TC)-2)]
#
#
#     for i in range((len(spaceLocations_TC)-3), n-1):
#
#         A[i+1][i] = A_STL[i-(len(spaceLocations_TC)-3)]
#
#
#     b = np.full(n, T_inital)
#
#
#
#     return A, b, n, C_STL
#
#
# def boundary_vector(t,n, C_STL):
#
#     BoundaryTerm = np.zeros(n)
#     T = T_inital + (T_final-T_inital)*(1-exp(-10*t))
#
#     BoundaryTerm[0] = betaTC*T*(1/DELTA_R_TC-1/(inDia_TC+DELTA_R_TC))
#
#     BoundaryTerm[-1] = (-LAMBDA*C_STL[-1]*T_amb)/(1+LAMBDA)
#
#     return BoundaryTerm

def matrix_builder(spaceLocations, timeLocations):


    A = np.zeros((len(spaceLocations), len(spaceLocations)))

    A[0][0] = 1

    for i in range(1, NUM_PTS_R_TC):
        A[i][i+1] = -betaTC*(1/DELTA_R_TC + 1/spaceLocations[i])
        A[i][i-1] = -betaTC*(1/DELTA_R_TC - 1/spaceLocations[i])
        A[i][i] = 1+(2*betaTC)/DELTA_R_TC

    A[NUM_PTS_R_TC][NUM_PTS_R_TC+1] = -GAMMA
    A[NUM_PTS_R_TC][NUM_PTS_R_TC] = 1+GAMMA
    A[NUM_PTS_R_TC][NUM_PTS_R_TC-1] = -1

    for i in range(NUM_PTS_R_TC+1, len(spaceLocations)-1):
        A[i][i+1] = -betaSTL*(1/DELTA_R_STL + 1/spaceLocations[i])
        A[i][i-1] = -betaSTL*(1/DELTA_R_STL - 1/spaceLocations[i])
        A[i][i] = 1+(2*betaSTL)/DELTA_R_STL

    A[-1][-1] = (1+LAMBDA)
    A[-1][-2] = -1

    b = np.full(len(spaceLocations), T_inital)
    b[0] =  T_inital+(T_final-T_inital)*(1 - exp(-10*DELTA_T))
    b[NUM_PTS_R_TC] = 0
    b[-1] = LAMBDA*T_amb


    return A,b

def matrix_builder_DISP(spaceLocations, timeLocations, T_history):


    A = np.zeros((len(spaceLocations), len(spaceLocations)))

    A[0][0] = Xi_TC/spaceLocations[0] - Phi_TC/DELTA_R_TC
    A[0][1] = Phi_TC / DELTA_R_TC


    for i in range(1, NUM_PTS_R_TC):
        A[i][i+1] = (2/Psi_TC)*(1/DELTA_R_TC + 1/spaceLocations[i])
        A[i][i-1] = (2/Psi_TC)*(1/DELTA_R_TC - 1/spaceLocations[i])
        A[i][i] = - (4/Psi_TC)*(1/DELTA_R_TC + DELTA_R_TC/(spaceLocations[i]**2))

    A[NUM_PTS_R_TC][NUM_PTS_R_TC+1] = -Eta_STL
    A[NUM_PTS_R_TC][NUM_PTS_R_TC] = Eta_TC+Eta_STL+Mu_TC-Mu_STL
    A[NUM_PTS_R_TC][NUM_PTS_R_TC-1] = -Eta_TC

    for i in range(NUM_PTS_R_TC+1, len(spaceLocations)-1):
        A[i][i + 1] = (2 / Psi_STL) * (1 / DELTA_R_STL + 1 / spaceLocations[i])
        A[i][i - 1] = (2 / Psi_STL) * (1 / DELTA_R_STL - 1 / spaceLocations[i])
        A[i][i] = - (4 / Psi_STL) * (1 / DELTA_R_STL + DELTA_R_STL / (spaceLocations[i] ** 2))

    A[-1][-1] = Xi_STL/spaceLocations[-1] + Phi_STL/DELTA_R_STL
    A[-1][-2] = - Phi_STL / DELTA_R_STL

    b = np.zeros((len(timeLocations), len(spaceLocations)))

    for i in range(0,b.shape[0]):
        b[i][0] = DIFF_TC*(T_history[i][0]-T_amb)

        for j in range(1, b.shape[1]-1):
            b[i][j] = T_history[i][j+1] - T_history[i][j-1]

        b[i][NUM_PTS_R_TC] = (Kappa_TC - K_STL)*T_history[i][NUM_PTS_R_TC] - (Kappa_TC - K_STL) * T_amb
        b[i][-1] = DIFF_STL * (T_history[i][-1] - T_amb)
    # b[0] =  D
    # b[NUM_PTS_R_TC] = 0
    # b[-1] = LAMBDA*T_amb


    return A,b







def meshing():

    spaceMeshLocations_TC = np.zeros(NUM_PTS_R_TC+1)
    spaceMeshLocations_STL = np.zeros(NUM_PTS_R_STL + 1)
    timeMeshLocations = np.zeros(NUM_PTS_T+1)

    for i in range(0,NUM_PTS_R_TC+1):
        spaceMeshLocations_TC[i] = inR_TC+i*DELTA_R_TC

    for i in range(0, NUM_PTS_R_STL + 1):
        spaceMeshLocations_STL[i] = inR_STL+i * DELTA_R_STL

    for i in range(0,NUM_PTS_T+1):
        timeMeshLocations[i] = i*DELTA_T

    # n = len(spaceMeshLocations_STL)+len(spaceMeshLocations_TC)-1
    # spaceMeshLocations = np.zeros(n)

    spaceMeshLocations_STL= np.delete(spaceMeshLocations_STL,0)

    spaceMeshLocations = np.concatenate((spaceMeshLocations_TC, spaceMeshLocations_STL))


    return spaceMeshLocations, timeMeshLocations


def SR_solver(A,b):
    """Sucessive Relaxation Solver"""

    n=len(b)
    relax_factor = 1.76
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
#         T[i] = T_inital+(T_final-T_inital)*(1 - exp(-10*i*DELTA_T))

    #
    # plt.plot(timeLocations, T)
    # plt.show()

# def b_matrix_builder(b,BoundaryTerm ):
#
#     b_total = b+BoundaryTerm
#
#     return b_total

# def cc(arg):
#     return mcolors.to_rgba(arg, alpha=0.6)

def plotter(T_history,u_history, spaceLocations, timeLocations):

    plt.contour(spaceLocations,timeLocations,T_history)
    plt.xlabel('Length (m)')
    plt.xlim(0.07,0.2)
    plt.ylim(0, 100)
    plt.ylabel('Time (s)')

    plt.show()

    plt.contour(spaceLocations, timeLocations, u_history)
    plt.xlabel('Length (m)')
    plt.ylabel('Time (s)')
    plt.xlim(0.07, 0.2)
    plt.ylim(0, 100)

    plt.show()

    # time = timeLocations.tolist()
    #
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # verts = []
    #
    # for i in range(0, T_history.shape[0]):
    #
    #     verts.append(list(zip(spaceLocations, T_history[i,:])))
    #
    # zs = [0.0, 1.0, 2.0]
    # verts = [list(zip(spaceLocations,T_history[0,:])),list(zip(spaceLocations,T_history[30,:])), list(zip(spaceLocations,T_history[80,:]))]
    #
    # poly = PolyCollection(verts, facecolors = ['r','g','b'])
    # poly.set_alpha(0.7)
    # ax.add_collection3d(poly,zs=zs, zdir='x')
    #
    # plt.show()
    #
    # fig = plt.figure()
    # ax = fig.gca(projection = '3d')
    # zs = range(0,T_history.shape[0])
    # verts = []
    #
    # for z in zs:
    #     ys = T_history[z,:]
    #     verts.append(list(zip(spaceLocations,ys)))
    #
    # poly = LineCollection(verts, linewidths=5)
    #
    # # ax.plot_trisurf(poly, zs=zs)
    # ax.add_collection3d(poly, zs= zs, zdir='y')
    # ax.set_xlim3d(0.07, 0.20)
    # ax.set_ylim3d(0,100)
    # ax.set_zlim3d(0,600)
    #
    #
    #
    # #
    # # plt.plot(spaceLocations, T_history[0,:], 'r')
    # # plt.plot(spaceLocations, T_history[20, :], 'g')
    # # plt.plot(spaceLocations, T_history[30, :],'b')
    # # plt.plot(spaceLocations, T_history[50, :], '-ro')
    # # plt.plot(spaceLocations, T_history[80, :], '-go')
    # # plt.xlabel('Length (m)')
    # # plt.ylabel('Temp')
    # plt.show()
    #
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # zs = range(0, u_history.shape[0])
    # verts = []
    #
    # for z in zs:
    #     ys = u_history[z, :]
    #     verts.append(list(zip(spaceLocations, ys)))
    #
    # poly = LineCollection(verts, linewidths=5)
    #
    # # ax.plot_trisurf(poly, zs=zs)
    # ax.add_collection3d(poly, zs=zs, zdir='y')
    # ax.set_xlim3d(0.07, 0.20)
    # ax.set_ylim3d(0, 100)
    # ax.set_zlim3d(0, 1)
    # plt.show()


    plt.plot(spaceLocations, u_history[0,:], 'r')
    plt.plot(spaceLocations, u_history[2, :], 'g')
    plt.plot(spaceLocations, u_history[3, :],'b')
    plt.plot(spaceLocations, u_history[5, :], '-ro')
    plt.plot(spaceLocations, u_history[100, :], '-go')
    plt.xlabel('Length (m)')
    plt.ylabel('displacement (m')
    plt.show()



# def main_solver():
#
#
#     spaceLocations_TC, spaceLocations_STL, timeLocations = meshing()
#
#
#
#     A, b, n, C_STL = matrix_builder(spaceLocations_STL, spaceLocations_TC)
#
#     BoundaryTerm = boundary_vector(DELTA_T, n, C_STL)
#     # BoundaryTerm[-1] = (-LAMBDA*C_STL[-1]*T_amb)/(1+LAMBDA)
#
#     b_total = b_matrix_builder(b, BoundaryTerm)
#
#     x_intial = np.linalg.solve(A, b_total)
#     x_curr = np.zeros(n)
#     T_history = np.zeros((len(timeLocations), n))
#     T_history_new = np.zeros((len(timeLocations), n+3))
#     T_history[0,:] = x_intial[:]
#
#     for i in range(1, len(timeLocations)):
#
#         BoundaryTerm = boundary_vector((i+1)*DELTA_T,n, C_STL)
#
#         x_total = b_matrix_builder(x_intial, BoundaryTerm)
#
#         x_curr = np.linalg.solve(A, x_total)
#
#         x_intial[:] = x_curr[:]
#
#         T_history[i,:] = x_curr[:]
#
#
#     for i in range(0, T_history.shape[0]):
#
#         temp = T_history[i,:]
#
#         index = int(n / 2)
#
#         T_interface = (temp[index-1]+ GAMMA*temp[index])/(1+GAMMA)
#
#         T_n = (temp[-1] + (LAMBDA*T_amb))/(1+LAMBDA)
#
#         T_0 = T_inital+(T_final-T_inital)*(1 - exp(-10*(i+1)*DELTA_T))
#
#
#         temp = np.insert(temp,index,T_interface)
#
#         temp = np.insert(temp,0,T_0)
#
#         temp = np.append(temp,T_n)
#
#         T_history_new[i,:] = temp
#
#         temp[:] = 0
#
#
#     spaceLocations_TC = spaceLocations_TC+inDia_TC
#     spaceLocations_STL = spaceLocations_STL+inDia_STL
#
#     spaceLocations_STL = np.delete(spaceLocations_STL,0)
#
#     spaceLocation_actual = np.concatenate((spaceLocations_TC,spaceLocations_STL))
#
#
#
#     return T_history_new, spaceLocation_actual, timeLocations
#

def temp_Solver(A, b, spaceLocations, timeLocations):

    x_inital =  np.linalg.solve(A, b)
    T_history = np.zeros((len(timeLocations), len(spaceLocations)))
    T_history[0, :] = T_inital
    T_history[0, -1] = T_amb
    T_history[1,:] = x_inital[:]

    for i in range(2, len(timeLocations)):

        x_inital[0] = T_inital+(T_final-T_inital)*(1 - exp(-10*(i+1)*DELTA_T))
        x_inital[NUM_PTS_R_TC] = 0
        x_inital[-1]= LAMBDA*T_amb

        x_curr = np.linalg.solve(A,x_inital)

        x_inital[:] = x_curr[:]
        T_history[i,:] = x_curr[:]
        x_curr[:] = 0

    return T_history

def disp_Solver(A, b, spaceLocations, timeLocations):


    u_history = np.zeros((len(timeLocations), len(spaceLocations)))

    for i in range(0, u_history.shape[0]):

        u_history[i,:] = np.linalg.solve(A,b[i,:])



    return u_history








spaceLocations, timeLocations = meshing()
#
A, b = matrix_builder(spaceLocations, timeLocations)

T_history = temp_Solver(A,b,spaceLocations, timeLocations)

# plotter(T_history, spaceLocations, timeLocations)

#
# _,_, timeLocations = meshing()
#
# test(timeLocations)
A_disp, b_disp = matrix_builder_DISP(spaceLocations,timeLocations, T_history)

u_history = disp_Solver(A_disp, b_disp, spaceLocations, timeLocations)

plotter(T_history, u_history, spaceLocations,timeLocations)




