__author__ = "Robin Goix <robin.goix.rg@gmail.com>"
__date__ = "2014-04-22"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from Square.py originally developed 
#   by Erich L Foster <erichlf@gmail.com>
#

'''
This problem allows to assess the Peregrine Code by simulating its 
solitons solutions.
It need to be run using the following parameters:
- Nx = 2047
- Ny = 1
- lambda0 = 1.
- a0 = 1.
- h0 = 1.
'''

from problembase import *
from numpy import *
from scipy.interpolate import interp1d
import numpy as np

x0 = -40.
x1 = 40.
y0 = -0.5
y1 = 0.5

# No-slip boundary
class Y_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS)
           
class X_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS)
           
# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        global x0, x1, y0, y1

        # Create mesh
        Nx = 2047
        Ny = 1
        lambda0 = float(options["lambda0"])
        x0 = x0/lambda0
        x1 = x1/lambda0
        y0 = y0/lambda0
        y1 = y1/lambda0
        mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)
        self.mesh = mesh
    def initial_conditions(self, V, Q):
       
        eta0 = genfromtxt('eta.txt')[np.newaxis] #Get an array of array with the height solution from the Matlab code
        eta0 = eta0[0] #Get the array with the height solution from the Matlab code
        eta00 = np.zeros(4096) #Create a new array twice as long as eta0

        u0 = genfromtxt('u.txt')[np.newaxis]#Get an array with the x-velocity solution from the Matlab code
        u0 = u0[0]  #Get the array with the x-velocity solution from the Matlab code
        u00 = np.zeros(4096) #Create a new array twice as long as u0
        
        #Fill the longer arrays to get the right values on the 2D-mesh coordinates
        i = 0
        while(i<=4094):
            i += 1
            eta00[i] = eta0[floor(i/2.)]
            u00[i] = u0[floor(i/2.)]

        eta_initial = Function(Q) #Create an empty function defined on the Q Space
        eta_initial.vector()[:] = eta00 #Fill the function's nodes with the values of the Matlab solution
        eta_0 = Expression("eta_initial",eta_initial=eta_initial)

        u_initial = Function(Q)
        u_initial.vector()[:]=u00
        u_0=Function(V)
        u_0=Expression(("u_initial","0.0"),u_initial=u_initial)
        
        return u_0, eta_0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bc_X = DirichletBC(V.sub(0), 0.0, X_SlipBoundary())
        bc_Y = DirichletBC(V.sub(1), 0.0, Y_SlipBoundary())
        
        bcs = [bc_X, bc_Y]
        
        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'PeregrineSoliton'
