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

#start location of our object
x2 = -3.
x3 = 3.
y2 = -3.
y3 = 3.

#Physical Parameters
hd = 2. #Depth of the pool in the central lane [m]
hb = 0.3 #Depth at the boundaries [m]
ad = 0.8 #height of the moving object [m]
g = 9.8 #Gravity
lambda0 = 20. #typical wavelength
a0 = 0.8 #Typical wave height
h0 = 2. #Typical depth

#scaling our domain
x0 = x0/lambda0
x1 = x1/lambda0
y0 = y0/lambda0
y1 = y1/lambda0
x2 = x2/lambda0
x3 = x3/lambda0
y2 = y2/lambda0
y3 = y3/lambda0

#Scaled Parameters
hd = hd/h0 #depth
ad = ad/a0 #height of the moving object
hb = hb/h0 #depth at the boundary
c0 = (h0*g)**(0.5)

#maximum velocity of wave object
vmax = ((hd*h0+ad*a0)*g)**(0.5) #Max Speed of the moving object [m.s^(-1)]

# Slip boundary
class Y_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS)

class X_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS)

class Object(Expression):
    '''
        This expression takes the CST parameters and outputs the y cooresponding to
        x for the wave object. We use 3rd order Bernstein polynomials and thus we
        require 4 w-parameter, while the N parameters define the type of shape.
        The following code has been adapted from the code
        written by Pramudita Satria Palar for matlab which can be found at
        http://www.mathworks.com/matlabcentral/fileexchange/42239-airfoil-generation-using-cst-parameterization-method
    '''
    def __init__(self, N, w, dz, t):
        global x0, x1, y0, y1, x2, x3, y2, y3, hd, hb, ad, a0, h0, c0, g, vmax

        self.N1 = N[0]
        self.N2 = N[1]
        self.w = w
        self.t = t
        self.dz = dz

    def eval(self, value, x):
        N1 = self.N1
        N2 = self.N2
        w = self.w
        dz = self.dz
        t = self.t

        u = vmax*lambda0/c0*t*exp(-4./(lambda0/c0*t+0.05))

        X = ((x[0] - x2 - u*t)/(x3 - x2), (x[1] - y2)/(y3 - y2))

        if(X[0]>=0. and X[0]<=1. and x[1]>=-0.25 and x[1]<=0.25):
            # Class function; taking input of N1 and N2
            C = X[0]**N1*(1. - X[0])**N2

            # Shape function; using Bernstein Polynomials
            n = len(w) - 1 # Order of Bernstein polynomials

            S = 0;
            for i in range(0,n):
                K = float(factorial(n)/(factorial(i)*factorial(n-i)));
                S = S + w[i]*K*X[0]**i*(1. - X[0])**(n-i)

            value[0] = C*S + X[0]*dz
        else:
            value[0] = 0.

    def value_shape(self):
        return ()

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        global x0, x1, y0, y1, x2, x3, y2, y3, hd, hb, ad, a0, h0, c0, g, vmax

        #Scaling Parameters
        self.sigma = h0/lambda0
        self.epsilon = a0/h0

        self.N = (0.005, 5E-6)
        self.w = (1.,1.,1.,1.)
        self.dz = 0.1

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        #Scaled Parameters
        self.t0 = 0.
        self.T = options['T']*c0/lambda0 #Final time
        self.k = options['dt']*c0/lambda0 #time step

        Q = FunctionSpace(self.mesh, 'CG', 1)
        self.zeta0 = Object(N=self.N,w=self.w,dz=self.dz,t=self.t0)

        #Defintion of the shape of the seabed
        self.D = str(hd) + ' - ' + str((hd-hb)/21.) + '*(x[1]>4./' \
            + str(lambda0) + ' ? 1. : 0.)*(' + str(lambda0) + '*x[1]-4.)' \
            + ' +  ' + str((hd-hb)/21.) + '*(x[1]<(-4./' + str(lambda0) \
            + ') ? 1. : 0.)*(' + str(lambda0) + '*x[1]+4.)'

        #Define the bounds
        self.ub = '(-' + str(ad) + '*(x[1]<0. ? 1. : 0.)' \
                        + '*(x[1]>' + str(-3./lambda0) \
                        + '? 1. : 0.)*(x[0]>' + str(-1.5/lambda0) + '? 1. : 0.)' \
                        + '*(x[0]<' + str(-1.5/lambda0) + '? 1. : 0.) < -0.5 ? -' \
                        + str(ad) + '*(x[1]<0. ? 1. : 0.)' \
                        + '*(x[1]>' + str(-3./lambda0) \
                        + '? 1. : 0.)*(x[0]>' + str(-1.5/lambda0) + '? 1. : 0.)'\
                        + '*(x[0]<' + str(2./lambda0) + '? 1. : 0.) : 0 )'
        self.lb = '(-' + str(ad) + '*(x[1]<' + str(3./lambda0) + ' ? 1. : 0.)' \
                        + '*(x[1]>' + str(-3./lambda0) \
                        + '? 1. : 0.)*(x[0]>' + str(-3./lambda0) + '? 1. : 0.)'\
                        + '*(x[0]<' + str(1.5/lambda0) + '? 1. : 0.) < -0.5 ? -' \
                        + str(ad) + '*(x[1]<' + str(3./lambda0) + '? 1. : 0.)'\
                        + '*(x[1]>' + str(-3./lambda0) \
                        + '? 1. : 0.)*(x[0]>' + str(-3./lambda0) + ' ? 1. : 0.)'\
                        + '*(x[0]<' + str(1.5/lambda0) + ' ? 1. : 0.) : 0 )'

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
