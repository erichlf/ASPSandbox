__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

# List of problems
problems = ['Channel', 'Drop', 'Cavity', 'Cylinder','Square', 'Point', 'Pool', 'PoolInMotion', 'PeregrineSoliton']

# Wrapper for problem classes
def Problem(name, options):
#   Return problem instance for given problem name"
    exec('from %s import Problem as NamedProblem' % name)
    return NamedProblem(options)
