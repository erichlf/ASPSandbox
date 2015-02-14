from SS_Cylinder import *
from SS_Cylinder import Problem as SS_Cylinder
import sys


class Problem(SS_Cylinder):

    '''
        Purely a container to call SS_Cylinder. This
        convert SS_Cylinder.py to a viscoelestic problem.
    '''

    def __init__(self, options, cube=False):
        SS_Cylinder.__init__(self, options, cube=cube)

    def functional(self, W, w):

        # nu = self.nu
        # mesh = W.mesh()
        # psimarker = PsiMarker()

        if W.mesh().topology().dim() == 2:
            (u, p, tau) = (as_vector((w[0], w[1])), w[2],
                           as_tensor(((w[3], w[4]), (w[5], w[6]))))
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
        return 'SS_3F_Cylinder'
