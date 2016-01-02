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
@brief Base class for implementing model specific potentials.
"""

from operator import mul, sub

class Model:
    """ Base class for implementing model specific potentials. """

    def __init__(self,
                 box,
                 particles,
                 cells,
                 max_interactions,
                 interaction_energy,
                 interaction_range):
        """ Initialise Model object. """

        self.box = box
        self.particles = particles
        self.cells = cells
        self.max_interactions = max_interactions
        self.interaction_energy = interaction_energy
        self.interaction_range = interaction_range

        # Work out squared cut-off distance.
        self.squared_cutoff_distance = interaction_range * interaction_range

    def compute_energy(self, particle, position, orientation):
        """ Compute the total interaction energy felt by a particle. """

        # Zero interaction energy.
        energy = 0

        # Check all neighbouring cells, including the same cell
        for cell in self.cells.neighbours[self.particles[particle].cell]:

            # Loop over all particles in cell.
            for x in range(0, self.cells.tally[cell]):

                neighbour = self.cells.particles[cell][x]

                # Make sure particles are different.
                if neighbour != particle:

                    # Calculate the model specific pair energy.
                    energy += self.compute_pair_energy(particle, position, orientation, neighbour,
                        self.particles[neighbour].position, self.particles[neighbour].orientation)

                    # Early exit for hard core overlaps and large finite energy repulsions.
                    if energy > 1e6:
                        return float("inf")

        return energy

    def compute_interactions(self, particle, position, orientation, interactions):
        """ Compute the number of interactions energy felt by a particle. """

        # Zero number of interactions.
        num_interactions = 0

        # Loop over all nearest neighbour cells.
        for cell in self.cells.neighbours[self.particles[particle].cell]:

            # Loop over all particles in cell.
            for x in range(0, self.cells.tally[cell]):

                neighbour = self.cells.particles[cell][x]

                # Make sure particles are different.
                if neighbour != particle:

                    # Calculate particle separation.
                    sep = map(sub, position, self.particles[neighbour].position)

                    # Compute minimum image.
                    self.box.minimum_image(sep)

                    # Compute squared norm of separation vector.
                    norm_sqd = sum(map(mul, sep, sep))

                    # Overlap if norm_sqd is less than particle diameter (box is scaled in diameter units).
                    if norm_sqd < self.squared_cutoff_distance:

                        # Check maxmium overlaps.
                        if num_interactions == self.max_interactions:
                            raise ValueError('Maximum number of interactions exceeded.')
                            raise SystemExit

                        interactions[num_interactions] = neighbour
                        num_interactions += 1

        return num_interactions

    def apply_post_move_updates(self, particle, position, orientation):
        """ Apply updates following a trial move. """

        # Copy coordinates and orientations.
        self.particles[particle].position = position
        self.particles[particle].orientation = orientation

        # Calculate the particle's cell index.
        new_cell = self.cells.get_cell(self.particles[particle])

        # Update cell lists if necessary.
        if self.particles[particle].cell != new_cell:
            self.cells.update_cell(new_cell, self.particles[particle], self.particles)

    def get_energy(self):
        """ Get the total system energy. """

        # Zero energy.
        energy = 0

        for x in range(0, len(self.particles)):
            energy += self.compute_energy(x, self.particles[x].position, self.particles[x].orientation);

        return energy/(2*len(self.particles))
