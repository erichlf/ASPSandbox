__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Problem as ProblemBase
from dolfin_adjoint import *


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.opt_control = None

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]

        try:
            x0 = float(options["x0"])
            x1 = float(options["x1"])
            y0 = float(options["y0"])
            y1 = float(options["y1"])
            self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)
        except:
            self.mesh = UnitSquareMesh(Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        try:
            self.kappa = Constant(options['kappa'])
        except:
            self.kappa = Constant(1E-2)  # heat Coefficient
        try:
            self.rho = Constant(options['rho'])
        except:
            self.rho = Constant(1.)  # density

        try:
            self.c = Constant(options['c'])
        except:
            self.c = Constant(1.)  # heat Coefficient

    def initial_conditions(self, V):
        u0 = project(Expression('sin(pi*x[0])*sin(pi*x[1])'), V)

        return u0

    def boundary_conditions(self, V, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, Constant(0.0), 'on_boundary')

        return bcs

    def F(self, t):
        if self.opt_control is None:
            W = FunctionSpace(self.mesh, "DG", 0)
            m = Function(W, name='Control')
        else:
            m = self.opt_control

        return m

    def Optimize(self, solver, V, u):
        # file for solution to optimization
        optfile = File(solver.s + '_Opt.pvd')

        x = SpatialCoordinate(V.mesh())
        d = sin(pi*x[0])*sin(pi*x[1])
        dp = project(d, V)

        # Functionnal to be minimized
        J = Functional((0.5 * inner(u - d, u - d)) * dx * dt[FINISH_TIME])

        Jhat = ReducedFunctional(J, Control(solver.f))  # Reduced Functional
        self.opt_control = minimize(Jhat, method = "L-BFGS-B",
                bounds = (-10, 10), options = {"gtol": 2e-8, "disp": True})
        self.opt_control = project(self.opt_control, V)
        if solver.plotSolution:
            plot(self.opt_control, title='Optimization result.')
            interactive()
        else:
            optfile << self.opt_control

    def __str__(self):
        return 'Square'
