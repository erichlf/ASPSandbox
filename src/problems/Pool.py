__author__ = "Robin Goix <robin.goix.rg@gmail.com>"
__date__ = "2014-04-22"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from Square.py originally developed
#   by Erich L Foster <erichlf@gmail.com>
#

'''
This is a simple problem with no forcing and on a rectangle.
This is basically a blank problem that we can adapt with optional inputs.
'''

from problembase import *
from numpy import array

x0 = -6.
x1 = 60.
y0 = -25.
y1 = 25.

#Physical Parameters
hd = 2. #Depth of the pool in the central lane [m]
hb = 0.3 #Depth at the boundaries [m]
ad = 0.8 #height of the moving object [m]

# Slip boundary
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

        global x0, x1, y0, y1, hd, hb, ad

        #Scaling Parameters
        g = 9.8 #Gravity
        self.lambda0 = options['lambda0'] #typical wavelength
        a0 = options['a0'] #Typical wave height
        h0 = options['h0'] #Typical depth
        self.sigma = h0/self.lambda0
        self.c0 = (h0*g)**(0.5)
        self.epsilon = a0/h0

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        x0 = x0/self.lambda0
        x1 = x1/self.lambda0
        y0 = y0/self.lambda0
        y1 = y1/self.lambda0
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        #Scaled Parameters
        self.t0 = 0.
        self.T = options['T']*self.c0/self.lambda0 #Final time
        self.k = options['dt']*self.c0/self.lambda0 #time step

        #DEFINITION OF THE OBJECT
        #Scaled Parameters
        self.hd = hd/h0 #depth
        self.ad = ad/a0 #height of the moving object
        self.hb = hb/h0 #depth at the boundary

        #Defintion of the shape of the seabed
        seabed = 'hd - (hd-hb)/21.*(x[1]>4./lambda0 ? 1. : 0.)*(lambda0*x[1]-4.)'\
                + ' + (hd-hb)/21.*(x[1]<(-4./lambda0) ? 1. : 0.)*(lambda0*x[1]+4.)'

        #Definition of the wave_object
        self.vmax = ((self.hd*h0+self.ad*a0)*g)**(0.5) #Max Speed of the moving object [m.s^(-1)]
        traj = 'vmax*lambda0/c0*t*exp(-4./(lambda0/c0*t+0.05))'
        movingObject = '-(x[1]<3/lambda0 ? 1. : 0.)*(x[1]>0 ? 1. : 0.)*(lambda0*x[0]-'+traj+'>-6 ? 1. : 0.)'\
                    +'*ad*0.5*0.5*(1.-tanh(0.5*lambda0*x[1]-2.))*(tanh(10*(1.-(lambda0*x[0]-' + traj + ')'\
                    +'-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0] - ' + traj + ')+pow(lambda0*x[1],2)/5 + 0.5)))'\
                    +'-(x[1]>-3/lambda0 ? 1. : 0.)*(x[1]<=0 ? 1. : 0.)*(lambda0*x[0]-'+traj+'>-6 ? 1. : 0.)'\
                    +'*ad*0.5*0.5*(1.+tanh(0.5*lambda0*x[1]+2.))*(tanh(10*(1.-(lambda0*x[0]-' + traj + ')'\
                    +'-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0] - ' + traj + ')+pow(lambda0*x[1],2)/5 + 0.5)))'
        self.D = seabed
        self.zeta0 = movingObject

    def initial_conditions(self, V, Q):
        u0 = Expression(("0.0", "0.0"))
        eta0 = Expression("0.0")

        return u0, eta0

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
        return 'Pool'
