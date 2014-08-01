__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-10-01"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Modified by Erich L Foster, 2013

from problembase import *
from Cylinder import Problem as Cylinder
from numpy import array

# this is basically the cylinder problem but the shape is a square/box
# Problem definition
class Problem(Cylinder):

    def __init__(self, options):
        self = Cylinder.__init__(self, options, cube=True)

    def __str__(self):
        return 'Cube'
