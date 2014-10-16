__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from dolfin import *
from AFES.solverbase import SolverBase as Solver
from AFES.problembase import ProblemBase as Problem

# Default options
OPTIONS = {
    'dim': 2,  # number of dimensions
    'Nx': 20,  # number of elements along x axis
    'Ny': 20,  # number of elements along y axis
    'Nz': 20,  # number of elements along z axis
    'k': 0.01,  # time-step
    'T': 1.0,  # final-time
    'theta': 0.5,  # time-stepping method
    'adaptive': False,  # mesh adaptivity
    'adapt_ratio': 0.1,  # percent of mesh to refine
    'max_adaptations': 30,  # max number of times to adapt mesh
    'adaptive_TOL': 1E-20,  # tolerance for terminating adaptivity
    'optimize': False,  # optimize as defined in solver
    'onDisk': 0.,  # percent of steps on disk
    'save_solution': False,
    'save_frequency': 1,
    'plot_solution': True,
    'debug': False,
    'check_mem_usage': False,
    'absolute_tolerance': 1e-25,
    'relative_tolerance': 1e-12,
    'monitor_convergence': False,
    'velocity_order': 1,  # order of velocity element
    'pressure_order': 1,  # order of height/pressure element
    'stabilize': False,  # stabilize?
    'initialMesh': None  # to use for initial computation
}
