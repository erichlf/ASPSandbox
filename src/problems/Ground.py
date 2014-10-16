__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from AFES import *
from AFES import Problem as ProblemBase

D = 1.5
W = D / 2.
n = 5.
TR = 10.
TA = 10.
omega = 7.27E-5
kappa_0 = 2.3
kappa_1 = 100
rho = 1500.
c = 1480


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = TR


class Surface2D(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.)


class Surface3D(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 0.)


class Problem(ProblemBase):

    '''
        This is an implementation of the problem described in the FEniCS
        documentation
        http://fenicsproject.org/documentation/tutorial/timedep.html#tut-timedep-diffusion2-sin-fig1
        where we simulate the effect of heating an inhomogeneous medium.
    '''

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.dim = options['dim']

        # Create mesh
        Nx = options['Nx']
        Ny = options['Ny']

        if self.dim == 2:
            self.mesh = RectangleMesh(-W / 2., -D, W / 2., 0., Nx, Ny)
        else:
            Nz = options['Nz']
            self.mesh = BoxMesh(- W / 2., - W / 2., - D, W / 2., W / 2., 0.,
                                Nx, Ny, Nz)

        self.t0 = 0.
        period = 2. * pi / omega
        self.T = 10. * period
        self.k = period / n

        if self.dim == 2:
            kappa = 'x[0] > -W/4 && x[0] < W/4 '\
                '&& x[1] > -D/2 && x[1] < -D/2 + D/4 ? '\
                'kappa_1 : kappa_0'
        else:
            kappa = 'x[0] > -W/4 && x[0] < W/4 '\
                '&& x[1] > -W/4 && x[1] < W/4 ' \
                '&& x[2] > -D/2 && x[2] < -D/2 + D/4 ? '\
                'kappa_1 : kappa_0'
        self.kappa = Expression(kappa, D=D, W=W,
                                kappa_0=kappa_0, kappa_1=kappa_1)
        self.rho = rho
        self.c = c

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        T0 = Expression('TR + TA*sin(omega*t)', TR=TR, TA=TA, omega=omega, t=t)
        if self.dim == 2:
            bcs = DirichletBC(W, T0, Surface2D())
        else:
            bcs = DirichletBC(W, T0, Surface3D())

        return bcs

    def F(self, t):
        # Mass source/sink
        return Constant(0)

    def update(self, W, t):

        return self.boundary_conditions(W, t)

    def functional(self, mesh, u):

        M = u * dx

        return M

    def __str__(self):
        return 'Ground'
