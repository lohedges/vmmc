"""
@package demo
@author Lester Hedges
@brief Derived class for implementing the square-well potential.
"""

from model import Model
from operator import mul, sub

class SquareWellium(Model):
    """ Derived class for implementing the square-well potential. """

    def __init__(self,
                 box,
                 particles,
                 cells,
                 max_interactions,
                 interaction_energy,
                 interaction_range):

        # Call base class constructor.
        Model.__init__(self, box, particles, cells, max_interactions, interaction_energy, interaction_range)

    def compute_pair_energy(self, particle1, position1, orientation1, particle2, position2, orientation2):
        """ Compute pair interaction energy between two particles. """

        # Calculate particle separation.
        sep = map(sub, position1, position2)

        # Compute minimum image.
        self.box.minimum_image(sep)

        # Compute squared norm of separation vector.
        norm_sqd = sum(map(mul, sep, sep))

        if norm_sqd < 1:
            return float("inf")
        if norm_sqd < self.squared_cutoff_distance:
            return -self.interaction_energy
        return 0
