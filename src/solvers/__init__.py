__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

# List of solvers
solvers = [
    'NSE',
    'SWE',
    'Peregrine',
    'QGE',
    'DensityNSE',
    'Heat',
    'ADR',
    'MixedWaveALE',
    'SS_NSE',
]


def Solver(name, options):  # Wrapper for solver classes
    # Return solver instance for given solver name"
    exec('from %s import Solver as NamedSolver' % name)
    return NamedSolver(options)
