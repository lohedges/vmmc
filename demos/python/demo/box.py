"""
@package demo
@author Lester Hedges
@brief Class for a periodic, cuboidal simulation box object.
"""

class Box:
    """ Class for periodic, cuboidal simulation boxes. """

    def __init__(self, size):
        """ Initialise with box size list. """

        # Check dimensionality.
        if len(size) != 2 and len(size) != 3:
            raise ValueError('Invalid dimensionality (must be 2 or 3).')
            raise SystemExit

        self.dimension = len(size)
        self.box_size = size
        self.pos_min_image = [0.5*x for x in size]
        self.neg_min_image = [-0.5*x for x in size]

    def periodic_boundaries(self, vec):
        """ Enforce periodic boundary conditions. """
        for i, x in enumerate(vec):
            if x < 0:
                vec[i] += self.box_size[i]
            elif x >= self.box_size[i]:
                vec[i] -= self.box_size[i]

    def minimum_image(self, vec):
        """ Enforce minumum image separation. """
        for i, x in enumerate(vec):
            if x < self.neg_min_image[i]:
                vec[i] += self.box_size[i]
            elif x >= self.pos_min_image[i]:
                vec[i] -= self.box_size[i]
