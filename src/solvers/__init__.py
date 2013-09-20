__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

# List of solvers
solvers = ['NSE', 'SteadyStokes', 'Stokes', 'LinearSWE', 'SWE', 'SWEAdvection',
    'DensityNSE']

# Wrapper for solver classes
def Solver(name, options):
#   Return solver instance for given solver name"
    exec('from %s import Solver as NamedSolver' % name)
    return NamedSolver(options)
