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
@brief Class for initialising particle configurations.
"""

import math
import random
from operator import mul, sub

class Initialise:
    """ Class initialising particle configurations. """

    def __init__(self):
        """ Constructor. """
        self.max_trials = 100000000

    def check_overlap(self, particle, particles, cells, box):
        """ Helper function to check particle overlap following trial insertions. """

        # Loop over all nearest neighbour cells.
        for cell in cells.neighbours[particle.cell]:

            # Loop over all particles in cell.
            for x in range(0, cells.tally[cell]):

                neighbour = cells.particles[cell][x]

                # Make sure particles are different.
                if neighbour != particle.index:

                    # Calculate particle separation.
                    sep = map(sub, particle.position, particles[neighbour].position)

                    # Compute minimum image.
                    box.minimum_image(sep)

                    # Compute squared norm of separation vector.
                    norm_sqd = sum(map(mul, sep, sep))

                    # Overlap if norm_sqd is less than particle diameter (box is scaled in diameter units).
                    if norm_sqd < 1:
                        return True

        # If we get this far, no overlaps.
        return False

    def random(self, particles, cells, box):
        """ Generate a random, non-overlapping, particle configuration. """

        for x in range(0, len(particles)):

            # Zero number of trial insertions.
            num_trials = 0

            # Whether particle overlaps.
            is_overlap = True

            # Set the particle index.
            particles[x].index = x

            # Keep trying to insert a particle until there is no overlap.
            while is_overlap:

                # Increment trials.
                num_trials += 1

                # Normalisation factor
                norm = 0

                # Resize position and orientation lists.
                particles[x].position = [0] * box.dimension
                particles[x].orientation = [0] * box.dimension

                # Generate a random position and orientation.
                for d in range(0, box.dimension):
                    particles[x].position[d] = random.random()*box.box_size[d]
                    particles[x].orientation[d] = random.gauss(0, 1)
                    norm += particles[x].orientation[d]*particles[x].orientation[d]

                # Normalise orientation.
                norm = math.sqrt(norm)
                particles[x].orientation[:] = [i/norm for i in particles[x].orientation]

                # Calculate cell index.
                particles[x].cell = cells.get_cell(particles[x])

                # See if there is any overlap between particles.
                is_overlap = self.check_overlap(particles[x], particles, cells, box)

                if num_trials == self.max_trials:
                    raise ValueError('Maximum number of trial insertions reached.')
                    raise SystemExit

            # Update cell list.
            cells.init_cell(particles[x].cell, particles[x]);
