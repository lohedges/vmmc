"""
@package demo
@author Lester Hedges
@brief Derived class for implementing the Lennard-Jones potential.
"""

import math
from model import Model
from operator import mul, sub

class LennardJonesium(Model):
    """ Derived class for implementing the Lennard-Jones potential. """

    def __init__(self,
                 box,
                 particles,
                 cells,
                 max_interactions,
                 interaction_energy,
                 interaction_range):

        # Call base class constructor.
        Model.__init__(self, box, particles, cells, max_interactions, interaction_energy, interaction_range)

        # Work out potential shift.
        self.potential_shift = math.pow(1/interaction_range, 12) - math.pow(1/interaction_range, 6)

    def compute_pair_energy(self, particle1, position1, orientation1, particle2, position2, orientation2):
        """ Compute pair interaction energy between two particles. """

        # Calculate particle separation.
        sep = map(sub, position1, position2)

        # Compute minimum image.
        self.box.minimum_image(sep)

        # Compute squared norm of separation vector.
        norm_sqd = sum(map(mul, sep, sep))

        if norm_sqd < self.squared_cutoff_distance:
            r2_inv = 1 / norm_sqd
            r6_inv = r2_inv*r2_inv*r2_inv
            return 4*self.interaction_energy*((r6_inv*r6_inv) - r6_inv - self.potential_shift)
        else:
            return 0
