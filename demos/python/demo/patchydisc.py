# Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
@package demo
@author Lester Hedges
@brief Derived class for implementing a simple two-dimension patchy disc potential.
"""

import math
from model import Model
from operator import mul, sub

class PatchyDisc(Model):
    """ Derived class for implementing a simple two-dimension patchy disc potential. """

    def __init__(self,
                 box,
                 particles,
                 cells,
                 max_interactions,
                 interaction_energy,
                 interaction_range):

        # Call base class constructor.
        Model.__init__(self, box, particles, cells, max_interactions, interaction_energy, interaction_range)

        # Work out the angle between patches.
        self.patch_separation = 2*math.pi / max_interactions

        # Create rotation matrix lookup tables.
        self.cos_theta = [0] * max_interactions
        self.sin_theta = [0] * max_interactions

        for x in range(0, max_interactions):
            self.cos_theta[x] = math.cos(x*self.patch_separation)
            self.sin_theta[x] = math.sin(x*self.patch_separation)

    def compute_pair_energy(self, particle1, position1, orientation1, particle2, position2, orientation2):
        """ Compute pair interaction energy between two discs. """

        # Calculate particle separation.
        sep = map(sub, position1, position2)

        # Compute minimum image.
        self.box.minimum_image(sep)

        # Compute squared norm of separation vector.
        norm_sqd = sum(map(mul, sep, sep))

        # Discs overlap.
        if norm_sqd < 1:
            return float("inf")

        # Zero energy sum.
        energy = 0

        # Initialise disc coordinate arrays.
        coord1 = [0] * 2
        coord2 = [0] * 2

        # Test interactions between all patch pairs.
        for i in range(0, self.max_interactions):

            # Compute position of patch i on first disc
            coord1[0] = position1[0] + 0.5*(orientation1[0]*self.cos_theta[i] - orientation1[1]*self.sin_theta[i])
            coord1[1] = position1[1] + 0.5*(orientation1[0]*self.sin_theta[i] + orientation1[1]*self.cos_theta[i])

            # Enforce periodic boundaries.
            self.box.periodic_boundaries(coord1)

            for j in range(0, self.max_interactions):

                # Compute position of patch j on second disc
                coord2[0] = position2[0] + 0.5*(orientation2[0]*self.cos_theta[i] - orientation2[1]*self.sin_theta[i])
                coord2[1] = position2[1] + 0.5*(orientation2[0]*self.sin_theta[i] + orientation2[1]*self.cos_theta[i])

                # Enforce periodic boundaries.
                self.box.periodic_boundaries(coord1)

                # Calculate patch separation.
                sep = map(sub, coord1, coord2)

                # Compute minimum image.
                self.box.minimum_image(sep)

                # Compute squared norm of separation vector.
                norm_sqd = sum(map(mul, sep, sep))

                # Patches interact.
                if norm_sqd < self.squared_cutoff_distance:
                    energy -= self.interaction_energy

        return energy

    def compute_interactions(self, particle, position, orientation, interactions):
        """ Compute the number of interactions between two discs. """

        # Zero number of interactions.
        num_interactions = 0

        # Loop over all nearest neighbour cells.
        for cell in self.cells.neighbours[self.particles[particle].cell]:

            # Loop over all particles in cell.
            for x in range(0, self.cells.tally[cell]):

                neighbour = self.cells.particles[cell][x]

                # Make sure particles are different.
                if neighbour != particle:

                    # Calculate pair energy.
                    energy = self.compute_pair_energy(particle, position, orientation, neighbour,
                        self.particles[neighbour].position, self.particles[neighbour].orientation)

                    # Discs interact.
                    if energy < 0:

                        # Check maxmium overlaps.
                        if num_interactions == self.max_interactions:
                            raise ValueError('Maximum number of interactions exceeded.')
                            raise SystemExit

                        interactions[num_interactions] = neighbour
                        num_interactions += 1

        return num_interactions
