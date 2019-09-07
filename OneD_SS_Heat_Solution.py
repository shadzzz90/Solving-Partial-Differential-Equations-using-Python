
# Problem 2,3,4
  # Take L = 0.08, K = 50
  # q = 1.5e6,  alpha =30, beta = 100

import numpy as np
from math import exp, sin, pi
import matplotlib.pyplot as plt

LENGTH = 0.08 # Length in meters

CONDUCTIVITY = 50

ALPHA = 30.0 # Temprature at x = 0

BETA = 100.0 # Temprature at x = L

NUMPOINTS = 10 # Number of Points
# NUMPOINTS2 = 10
# NUMPOINTS3 = 50

GRIDSPACING = LENGTH/NUMPOINTS # Grid spacing

# GRIDSPACING2 = LENGTH/NUMPOINTS2
#
# GRIDSPACING3 = LENGTH/NUMPOINTS3


def mesh(gridSpacing, numPoints):

    """Finds the mesh coordinates for the problems """

    meshLocations = []

    for pt in range(0, numPoints+1):

        meshLocations.append(pt*gridSpacing)

    return meshLocations

def heatsource(L, x):

    """Defines the internal heat source"""

    halfL = L/2
    #q = x**4
    #q = 900*exp(-75*x)
    # q = sin(15*pi*x)

    q = 1.5e6

    return q

def solver(A, b):

    """Solves the matrix equation"""
    A = np.matrix(A)
    A_inv = np.linalg.inv(A)
    x = np.dot(A_inv,b)

    return x

def formAF_fourthorder(numPoints, gridSpacing, meshLocations, conductivity, Length, uO, uL):

    """Generate equation for A and F for fourth order scheme"""
    A = np.zeros((numPoints + 1, numPoints + 1))
    A[0][0] = 1
    A[1][0] = 1
    A[1][1] = -2
    A[1][2] = 1
    A[numPoints - 1][numPoints - 2] = 1
    A[numPoints - 1][numPoints-1] = -2
    A[numPoints-1][numPoints] = 1
    A[numPoints][numPoints] = 1

    for index in range(2, numPoints-1):
       A[index][index] = -30
       A[index][index - 1] = 16.0
       A[index][index + 1] = 16.0
       A[index][index - 2] = -1.0
       A[index][index + 2] = -1.0


    F = np.zeros(numPoints + 1)
    F[1] = -heatsource(Length, meshLocations[1])*gridSpacing**2/conductivity
    F[numPoints-1] = -heatsource(Length, meshLocations[numPoints-1])*gridSpacing**2/conductivity

    for index in range(2, numPoints-1):
         F[index] = -12*heatsource(Length, meshLocations[index]) * gridSpacing ** 2 / conductivity

    F[0] = uO
    F[numPoints] = uL

    return A, F



def formAF(numPoints, gridSpacing, meshLocations, conductivity, Length, uO, uL ):

    """Generates the equation for A and F """

    A = np.zeros((numPoints+1, numPoints+1))
    A[0][0] = 1
    A[numPoints][numPoints] = 1

    for index in range(1, numPoints):
        A[index][index] = -2.0
        A[index][index-1] = 1.0
        A[index][index+1] = 1.0

    F = np.zeros(numPoints+1)

    for index in range(1,numPoints):
        F[index] = -heatsource(Length, meshLocations[index])*gridSpacing**2/conductivity

    F[0] = uO
    F[numPoints] = uL

    return A, F

def plottemp(meshLocations, u, u2, xmin, xmax): #, umin, umax):  #, meshLocations2, u2, meshLocations3, u3,  xmin, xmax):

    """This function plots tempreature"""
    meshLocations = np.array([meshLocations])
    u = np.array(u)
    # meshLocations2 = np.array([meshLocations2])
    u2 = np.array(u2)
    # meshLocations3 = np.array([meshLocations3])
    # u3 = np.array(u3)

    plt.plot(meshLocations.flatten(), u.flatten(),'-rs')
    plt.plot(meshLocations.flatten(), u2.flatten(), '-bo')
    # plt.plot(meshLocations3.flatten(), u3.flatten(), '-go')
    plt.gca().legend(('Second Order Accurate', 'Fourth Order Accurate', '$n=50$'))

    plt.title('Steady State Tempreature Distribution in the Plane Wall  \n ($q  = 1.5 x 10^6$, $\\alpha = 30, \\beta = 100$) \n Fourth Order Scheme vs Second Order Scheme')
    plt.xlabel('x(m)')
    plt.ylabel('Tempreature (C)')
    plt.xlim([xmin, xmax])
    # plt.ylim([umin, umax])
    plt.show()


# #1
# meshLocations = mesh(GRIDSPACING, NUMPOINTS)  # Find the mesh coordinates
#
# A, F = formAF(numPoints=NUMPOINTS, meshLocations=meshLocations, gridSpacing=GRIDSPACING,
#              conductivity= CONDUCTIVITY, Length=LENGTH, uO=ALPHA, uL= BETA)             # Find the values of A and F
#
# u = solver(A, F) # Solve the Matrix
#
# xmin = 0.0
# xmax = LENGTH
#
# umin = 0.8*u.min()
# umax = 1.2*u.max()
#
# plottemp(meshLocations, u, xmax=xmax, xmin=xmin, umax=umax, umin=umin)
#
# #2
#
# meshLocations2 = mesh(GRIDSPACING2, NUMPOINTS2)  # Find the mesh coordinates
#
# A, F = formAF(numPoints=NUMPOINTS2, meshLocations=meshLocations2, gridSpacing=GRIDSPACING2,
#              conductivity= CONDUCTIVITY, Length=LENGTH, uO=ALPHA, uL= BETA)             # Find the values of A and F
#
# u2 = solver(A, F) # Solve the Matrix
#
# xmin = 0.0
# xmax = LENGTH
# # umin = 0.8*u.min()
# # umax = 1.2*u.max()
#
# #3
#
# meshLocations3 = mesh(GRIDSPACING3, NUMPOINTS3)  # Find the mesh coordinates
#
# A, F = formAF(numPoints=NUMPOINTS3, meshLocations=meshLocations3, gridSpacing=GRIDSPACING3,
#               conductivity=CONDUCTIVITY, Length=LENGTH, uO=ALPHA, uL=BETA)  # Find the values of A and F
#
# u3 = solver(A, F)  # Solve the Matrix
#
# xmin = 0.0
# xmax = LENGTH
# # umin = 0.8*u.min()
# # umax = 1.2*u.max()
#
# plottemp(meshLocations=meshLocations, u=u, meshLocations2=meshLocations2, u2=u2,
#          meshLocations3= meshLocations3, u3= u3,  xmin=xmin, xmax=xmax)
#
#

#
# meshLocations = mesh(GRIDSPACING, NUMPOINTS)
#
# A, F = formAF_fourthorder(numPoints=NUMPOINTS,gridSpacing=GRIDSPACING,meshLocations= meshLocations,
#                          conductivity=CONDUCTIVITY,Length=LENGTH, uO=ALPHA, uL=BETA)
# u_4th_order = solver(A,F)
#
# xmin = 0.0
# xmax = LENGTH
# umin = 0.8*u_4th_order.min()
# umax = 1.2*u_4th_order.max()
#
# plottemp(meshLocations,u_4th_order, xmin=xmin, xmax=xmax, umin=umin, umax=umax)


meshLocations = mesh(GRIDSPACING, NUMPOINTS)  # Find the mesh coordinates

A2, F2 = formAF(numPoints=NUMPOINTS, meshLocations=meshLocations, gridSpacing=GRIDSPACING,
             conductivity= CONDUCTIVITY, Length=LENGTH, uO=ALPHA, uL= BETA)             # Find the values of A and F

u_2nd_order = solver(A2, F2) # Solve the Matrix

A4, F4 = formAF_fourthorder(numPoints=NUMPOINTS,gridSpacing=GRIDSPACING,meshLocations= meshLocations,
                         conductivity=CONDUCTIVITY,Length=LENGTH, uO=ALPHA, uL=BETA)
u_4th_order = solver(A4, F4)

xmin = 0.0
xmax = LENGTH

plottemp(meshLocations,u_2nd_order, u_4th_order, xmin=xmin, xmax=xmax)