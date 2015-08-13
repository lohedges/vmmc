"""
@package demo
@author Lester Hedges
@brief Class for a simple particle container.
"""

class Particle:
    """ Class for a simple particle container. """

    def __init__(self):
        self.index = None       # Particle index.
        self.position = []      # Position vector.
        self.orientation = []   # Orientation unit vector.
        self.cell = None        # Cell index.
        self.pos_cell = None    # Position in the cell list.
