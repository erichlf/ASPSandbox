__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

# List of problems
problems = ['Channel', 'Drop', 'Cavity', 'Cylinder']

# Wrapper for problem classes
def Problem(name, options):
#   Return problem instance for given problem name"
    exec('from %s import Problem as NamedProblem' % name)
    return NamedProblem(options)
