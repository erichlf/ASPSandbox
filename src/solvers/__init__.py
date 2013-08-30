__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

# List of solvers
solvers = ['NSE', 'Stokes', 'LinearSWE', 'SWE']

# Wrapper for solver classes
def Solver(name, options):
#   Return solver instance for given solver name"
    exec('from %s import Solver as NamedSolver' % name)
    return NamedSolver(options)
