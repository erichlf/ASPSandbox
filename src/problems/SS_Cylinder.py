from Cylinder import *
from Cylinder import Problem as Cylinder


class Problem(Cylinder):

    '''
        Purely a container to call Cylinder. This
        convert Cylinder.py to a steady state problem.
    '''

    def __init__(self, options, cube=False):
        Cylinder.__init__(self, options, cube=cube)

        H = ymax
        if self.mesh.topology().dim() == 2:
            self.noSlip = Constant((0, 0))
            self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)', '0.0'),
                                Um=Um, H=ymax)
            self.Ubar = 4. / 3. * Um * ymax * (H - ymax / 2.) / (H * H)
        else:
            self.noSlip = Constant((0, 0, 0))
            self.U = Expression(('16*Um*x[1]*x[2]*(H - x[1])*(H - x[2])' +
                                 '/pow(H,4)', '0.0', '0.0'),
                                Um=Um ** 2, H=ymax)
            self.Ubar = 16. / 9. * Um * ymax * zmax * (H - ymax / 2.) \
                * (H - zmax / 2.) / pow(H, 4)

        self.Re = self.Ubar*Diameter/self.nu

    def boundary_conditions(self, W):
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip, NoSlipBoundary(self.dim))

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def F(self):
        return Constant((0.0, 0.0))

    def functional(self, W, w):

        # nu = self.nu
        # mesh = W.mesh()
        # psimarker = PsiMarker()

        if W.mesh().topology().dim == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
        else:
            print
            print 'This problem can only be used in 2D.'
            sys.exit(1)

        # n = FacetNormal(mesh)

        # I = Identity(2)
        # sigma = p*I - nu*self.epsilon(u)
        # theta = Constant((1.0, 0.0))

        # g = Expression(
        #     ("200.0*exp(-200.0*(pow(x[0] - 0.5, 2) + pow(x[1] - 0.3, 2)))",
        #      "0.0"))

        # M1 = psimarker*p*n[0]*ds  # Drag (only pressure)
        # M2 = psimarker*p*n[1]*ds  # Lift (only pressure)
        # M3 = inner(g, u)*dx  # Mean of the velocity in a region
        # M4 = psimarker*dot(dot(sigma, n), theta)*ds  # Drag (full stress)
        M5 = u[0]*dx  # Mean of the x-velocity in the whole domain

        return M5

    def __str__(self):
        return 'SS_Cylinder'
