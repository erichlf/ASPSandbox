__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from problembase import *
import numpy as np

#outer square dimensions
Xmin = 0.
Xmax = 1.
Ymin = 0.
Ymax = 1.
#inner square dimensions
xmin = 4./9.
xmax = 5./9.
ymin = 4./9.
ymax = 5./9.

kappa1 = 1000.;  kappa2 = 1.
theta = pi/4.
rho = 1.; c = 1.
TR = 0.5; TA = 0.5; omega = 3.*pi

class InitialConditions(Expression):
    def eval(self,values,x):
        values[0] = 0.

class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0],Xmin) or near(x[0],Xmax) \
                or near(x[1],Ymin) or near(x[1],Ymax))

class InnerBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0],xmin) or near(x[0],xmax) \
                or near(x[1],ymin) or near(x[1],ymax))


# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        self.Nx = options["Nx"]

        outerRect = Rectangle(Xmin, Ymin, Xmax, Ymax)
        innerRect = Rectangle(xmin, ymin, xmax, ymax)
        domain = outerRect - innerRect
        self.mesh = Mesh(domain, self.Nx)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']


        K = np.array([[kappa1, 0],[0, kappa2]])
        Theta = np.array([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])
        kappa = np.dot(Theta,np.dot(K,np.transpose(Theta)))
        self.kappa = Expression((('k1','k2'), ('k3','k4')), \
                k1=kappa[0,0], k2=kappa[0,1], k3=kappa[1,0], k4=kappa[1,1])
        self.rho = rho
        self.c = c

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        g = Expression('TR - TA*cos(omega*t)', TR=TR, TA=TA, omega=omega, t=t)
        bc0 = DirichletBC(W, Constant(0.0), OuterBoundary())
        bc1 = DirichletBC(W, g, InnerBoundary())

        return [bc0, bc1]

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'Plate'
