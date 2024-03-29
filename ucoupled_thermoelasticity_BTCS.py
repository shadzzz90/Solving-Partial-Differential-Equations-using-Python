
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter
from math import *
from matplotlib import colors as mcolors
from matplotlib import cm

NUM_PTS_R_TC = 100
NUM_PTS_R_STL = 4000
NUM_PTS_T = 100
TIME = 100 # secs
T_inital = 293 # K
T_amb = 303 # K
T_final = 573 # K
# T_final = 303 # K

h = 200 # W/m^2-K
RHO_TC = 15.25*1e3 # Tungsten Carbide Kg/m^3
Cp_TC = 0.184 *1e3 # KJ/Kg-K
K_TC = 28 # W/m-K

# RHO_STL = 7833 # Steel 0.5% C Kg/m^3
# Cp_STL = 0.465*1e3 # KJ/Kg-K
# K_STL = 54 # W/m-K

RHO_STL = RHO_TC # Steel 0.5% C Kg/m^3
Cp_STL = Cp_TC # KJ/Kg-K
K_STL = K_TC # W/m-K

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
# nu_STL = 0.3
nu_STL = nu_TC

Coeff_Thermal_Exp_TC = 4.5 * 1e-6 # /K
# Coeff_Thermal_Exp_STL = 1.25 * 1e-5 # /K

Coeff_Thermal_Exp_STL = Coeff_Thermal_Exp_TC

E_TC = 620*1e9 # GPa
# E_STL = 210*1e9 # GPa

E_STL = E_TC

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

tau_TC = E_TC*(1-nu_TC)/((1+nu_TC)*(1-2*nu_TC))
tau_STL = E_STL*(1-nu_STL)/((1+nu_STL)*(1-2*nu_STL))

Theta_TC = 2*E_TC*nu_TC/((1+nu_TC)*(1-2*nu_TC))
Theta_STL = 2*E_STL*nu_STL/((1+nu_STL)*(1-2*nu_STL))





def matrix_builder_TEMP(spaceLocations, timeLocations):


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
            if j == NUM_PTS_R_TC:
                b[i][j] = (Kappa_TC - Kappa_STL) * (T_history[i][j] - T_amb)
            else:
                b[i][j] = T_history[i][j+1] - T_history[i][j-1]


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


def SR_solver(A,b, relax_factor):
    """Sucessive Relaxation Solver"""

    n=len(b)
    # relax_factor = 1.76
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

    return x_next, iteration_number


def plotter(T_history,u_history, stress_history,  spaceLocations, timeLocations):

    plt.contourf(spaceLocations,timeLocations,T_history, cmap = 'rainbow')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Tempreature (K)')
    plt.xlabel('Length (m)')
    plt.xlim(0.07,0.2)
    plt.ylim(0, 100)
    plt.ylabel('Time (s)')

    plt.show()

    plt.contourf(spaceLocations, timeLocations, u_history, cmap = 'rainbow')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Displacement (m)')
    plt.xlabel('Length (m)')
    plt.ylabel('Time (s)')
    plt.xlim(0.07, 0.2)
    plt.ylim(0, 100)

    plt.show()

    plt.contourf(spaceLocations, timeLocations, stress_history, cmap ='rainbow')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Stress (Pa)')
    plt.xlabel('Length (m)')
    plt.ylabel('Time (s)')
    plt.xlim(0.07, 0.2)
    plt.ylim(0, 100)

    plt.show()


    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    zs = range(0,T_history.shape[0])
    verts = []

    for z in zs:
        ys = T_history[z,:]
        verts.append(list(zip(spaceLocations,ys)))

    poly = LineCollection(verts, linewidths=1)

    # ax.plot_trisurf(poly, zs=zs)
    ax.add_collection3d(poly, zs= zs, zdir='y')
    ax.set_xlim3d(0.07, 0.20)
    ax.set_xlabel('Length (m)')
    ax.set_ylim3d(0,100)
    ax.set_ylabel('Time (secs)')
    ax.set_zlim3d(0,600)
    ax.set_zlabel('Tempreature (K)')
    plt.show()



    fig = plt.figure()
    ax = fig.gca(projection='3d')
    zs = range(0, u_history.shape[0])
    verts = []

    for z in zs:
        ys = u_history[z, :]
        verts.append(list(zip(spaceLocations, ys)))

    poly = LineCollection(verts, linewidths=1)

    # ax.plot_trisurf(poly, zs=zs)
    ax.add_collection3d(poly, zs=zs, zdir='y')
    ax.set_xlim3d(0.07, 0.20)
    ax.set_xlabel('Length (m)')
    ax.set_ylim3d(0, 100)
    ax.set_ylabel('Time (secs)')
    ax.set_zlim3d(-3e-4, 1e-4 )
    ax.set_zlabel('Displacement (m)')
    plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    zs = range(0, stress_history.shape[0])
    verts = []

    for z in zs:
        ys = stress_history[z, :]
        verts.append(list(zip(spaceLocations, ys)))

    poly = LineCollection(verts, linewidths=1)

    # ax.plot_trisurf(poly, zs=zs)
    ax.add_collection3d(poly, zs=zs, zdir='y')
    ax.set_xlim3d(0.07, 0.20)
    ax.set_xlabel('Length (m)')
    ax.set_ylim3d(0, 100)
    ax.set_ylabel('Time (secs)')
    ax.set_zlim3d(-1e8, 9e8)
    ax.set_zlabel('Stress (Pa)')

    plt.show()


    plt.plot(spaceLocations, T_history[0,:], 'r')
    plt.plot(spaceLocations, T_history[10, :], 'g')
    plt.plot(spaceLocations, T_history[25, :],'b')
    plt.plot(spaceLocations, T_history[50, :], '-ro')
    plt.plot(spaceLocations, T_history[100, :], '-go')
    plt.xlabel('Length (m)')
    plt.title('Tempreature (K) vs Length(m)')
    plt.ylabel('Tempearture (K)')
    label0 = 't = 0 secs'
    label1 = 't = 10 secs'
    label2 = 't = 25 secs'
    label3 = 't = 50 secs'
    label4 = 't = 100 secs'

    plt.gca().legend((label0, label1, label2, label3, label4), loc='upper right')
    plt.show()


    plt.plot(spaceLocations, u_history[0,:], 'r')
    plt.plot(spaceLocations, u_history[10, :], 'g')
    plt.plot(spaceLocations, u_history[25, :],'b')
    plt.plot(spaceLocations, u_history[50, :], '-ro')
    plt.plot(spaceLocations, u_history[100, :], '-go')
    plt.xlabel('Length (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacement (m) vs Length(m)')
    label0 = 't = 0 secs'
    label1 = 't = 10 secs'
    label2 = 't = 25 secs'
    label3 = 't = 50 secs'
    label4 = 't = 100 secs'

    plt.gca().legend((label0, label1, label2, label3, label4), loc='upper right')
    plt.show()

    plt.plot(spaceLocations, stress_history[0, :], 'r')
    plt.plot(spaceLocations, stress_history[2, :], 'g')
    plt.plot(spaceLocations, stress_history[3, :], 'b')
    plt.plot(spaceLocations, stress_history[5, :], '-ro')
    plt.plot(spaceLocations, stress_history[100, :], '-go')
    plt.xlabel('Length (m)')
    plt.ylabel('Stress (Pa)')
    plt.title('Stress (Pa) vs Length(m)')
    label0 = 't = 0 secs'
    label1 = 't = 10 secs'
    label2 = 't = 25 secs'
    label3 = 't = 50 secs'
    label4 = 't = 100 secs'

    plt.gca().legend((label0, label1, label2, label3, label4), loc='upper right')
    plt.show()




def temp_Solver(A, b, spaceLocations, timeLocations, relaxFactor):

    # iteration = 0
    # x_inital, itr =  SR_solver(A, b, relaxFactor)
    x_inital= np.linalg.solve(A, b)
    # iteration = iteration+itr

    T_history = np.zeros((len(timeLocations), len(spaceLocations)))
    T_history[0, :] = T_inital
    # T_history[0, -1] = T_amb
    T_history[1,:] = x_inital[:]

    for i in range(2, len(timeLocations)):

        x_inital[0] = T_inital+(T_final-T_inital)*(1 - exp(-10*(i+1)*DELTA_T))
        x_inital[NUM_PTS_R_TC] = 0
        x_inital[-1]= LAMBDA*T_amb

        # x_curr, itr = SR_solver(A,x_inital, relaxFactor)
        x_curr = np.linalg.solve(A,x_inital)
        # iteration = iteration+itr

        x_inital[:] = x_curr[:]
        T_history[i,:] = x_curr[:]
        x_curr[:] = 0

    # return T_history, iteration
    return T_history

def disp_Solver(A, b, spaceLocations, timeLocations, relaxFactor):

    iterations = 0

    u_history = np.zeros((len(timeLocations), len(spaceLocations)))

    for i in range(0, u_history.shape[0]):

        # u_history[i,:], itr = SR_solver(A,b[i,:], relaxFactor)
        u_history[i, :] = np.linalg.solve(A, b[i, :])
        # iterations = iterations + itr



    # return u_history, iterations
    return u_history

def stress_Solver(T_history, u_history, spaceLocations, timeLocations):

    stress_history = np.zeros_like(T_history)

    for i in range(0, stress_history.shape[0]):
        stress_history[i][0] = 0

        for j in range(1, NUM_PTS_R_TC+1):
            stress_history[i][j] = (tau_TC/DELTA_R_TC + Theta_TC/spaceLocations[j])*u_history[i][j] - (tau_TC/DELTA_R_TC)*u_history[i][j-1]-Kappa_TC*(T_history[i][j] - T_amb)

        for j in range(NUM_PTS_R_TC+1,stress_history.shape[1]-1):
            stress_history[i][j] = (tau_STL /DELTA_R_STL + Theta_STL / spaceLocations[j]) * u_history[i][j] - (tau_STL /DELTA_R_STL) * u_history[i][j - 1] - Kappa_STL * (T_history[i][j] - T_amb)

        # stress_history[i][NUM_PTS_R_TC] = stress_history[i][NUM_PTS_R_TC+1]
        stress_history[i][-1] = 0


    return stress_history



spaceLocations, timeLocations = meshing()

# relaxFactors = np.arange(0.1,1.98,0.01)
relaxFactors = np.full(1,1.97)
# iterations = np.zeros_like(relaxFactors)
index = 0

for relaxFactor in relaxFactors:

    A, b = matrix_builder_TEMP(spaceLocations, timeLocations)
    # T_history, temp_iteration = temp_Solver(A,b,spaceLocations, timeLocations, relaxFactor)
    T_history = temp_Solver(A, b, spaceLocations, timeLocations, relaxFactor)
    A_disp, b_disp = matrix_builder_DISP(spaceLocations,timeLocations, T_history)
    # u_history, disp_iteration = disp_Solver(A_disp, b_disp, spaceLocations, timeLocations, relaxFactor)
    u_history = disp_Solver(A_disp, b_disp, spaceLocations, timeLocations, relaxFactor)
    stress_history = stress_Solver(T_history,u_history,spaceLocations, timeLocations)
    plotter(T_history, u_history, stress_history, spaceLocations,timeLocations)
    # total_iterations = (temp_iteration+disp_iteration)
    # print('Total Iteration for relax factor {0:0f} is {0:0}'.format(relaxFactor,total_iterations))
    # iterations[index] = total_iterations
    # index = index+1

# np.save('relaxation_factor', relaxFactors)
# np.save('Total_iterations', iterations)
#
# plt.plot(relaxFactors, iterations)
# plt.title("Parametric Study for optimal relaxation factor")
# plt.xlabel("Relaxation Factor")
# plt.ylabel("Total Iterations")
# plt.show()






