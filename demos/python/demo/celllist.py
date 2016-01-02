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
@brief Class for the cell list object.
"""

import math

class CellList:
    """ Class for the cell list. """

    def __init__(self, dimension = None, box_size = None, interaction_range = None):
        """ Cell list constructor. """

        # Set default dimension.
        self.dimension = 3

        if dimension is None and box_size is None and interaction_range is None:
            pass
        else:
            self.dimension = dimension
            self.initialise(box_size, interaction_range)

    def initialise(self, box_size, interaction_range):
        """ Initalise the cell list. """

        # Initialise data structures.

        self.cells_per_axis = [0] * self.dimension  # Number of cells per axis.
        self.cell_spacing = [0] * self.dimension    # Cell spacing along each axis.

        # Determine minimal cell spacing.
        for x in range(0, self.dimension):

            # Initialise number of cells.
            self.cells_per_axis[x] = 1

            while ((box_size[x] / self.cells_per_axis[x]) > interaction_range):
                self.cells_per_axis[x] += 1

            self.cells_per_axis[x] -=1
            self.cell_spacing[x] = box_size[x] / self.cells_per_axis[x]

            # Check that the number of cells per axis is large enough.
            if self.cells_per_axis[x] < 3:
                raise ValueError('Simulation box is too small for cell lists.')
                raise SystemExit

        # Estimate the number of particles per cell from the interaction range.
        # (Assumes the diameter is one.)
        if self.dimension is 3:
            self.max_particles = int((self.cell_spacing[0]*self.cell_spacing[1]*self.cell_spacing[2]) / ((4/3)*math.pi*0.5*0.5*0.5))
        else:
            self.max_particles = int((self.cell_spacing[0]*self.cell_spacing[1]) / (math.pi*0.5*0.5))

        # Add a buffer, e.g. if particles can overlap.
        self.max_particles += 10

        # Work out the total number of cells.
        self.num_cells = self.cells_per_axis[0]*self.cells_per_axis[1]
        if self.dimension is 3:
            self.num_cells *= self.cells_per_axis[2]

        # Resize data structures.

        self.index = [0] * self.num_cells       # Index of each cell.
        self.tally = [0] * self.num_cells       # Number of particles in each cell.
        self.particles = [0] * self.num_cells   # Indices of particles in each cell.
        self.neighbours = [0] * self.num_cells  # Indices of nearest neighbour cells.

        if self.dimension is 3:

            # Set number of nearest neighbour cells.
            self.num_neighbours = 27

            # Loop over all cells in x direction.
            for i in range(0, self.cells_per_axis[0]):

                # Loop over all cells in y direction.
                for j in range(0, self.cells_per_axis[1]):

                    # Loop over all cells in z direction.
                    for k in range(0, self.cells_per_axis[2]):

                        # Evaluate cell index.
                        m = i + self.cells_per_axis[0]*j + self.cells_per_axis[0]*self.cells_per_axis[1]*k

                        # Resize particle and nearest neighbour lists.
                        self.particles[m] = [0] * self.max_particles
                        self.neighbours[m] = [0] * self.num_neighbours

                        # Zero particle tally.
                        self.tally[m] = 0

                        # Set cell index
                        self.index[m] = m

                        # Initalise nearest neighbour counter
                        nn_count = 0

                        # x loop for nearest neighbours
                        for a in range(0, 3):

                            x = (i + (a-1) + self.cells_per_axis[0]) % self.cells_per_axis[0]

                            # y loop for nearest neighbours
                            for b in range(0, 3):

                                y = (j + (b-1) + self.cells_per_axis[1]) % self.cells_per_axis[1]

                                # z loop for nearest neighbours
                                for c in range(0, 3):

                                    z = (k + (c-1) + self.cells_per_axis[2]) % self.cells_per_axis[2]

                                    # Nearest neighbour cell index.
                                    nn = x + y*self.cells_per_axis[0] + z*self.cells_per_axis[0]*self.cells_per_axis[1]

                                    self.neighbours[m][nn_count] = nn
                                    nn_count += 1

        else:

            # Set number of nearest neighbour cells.
            self.num_neighbours = 9

            # Loop over all cells in x direction.
            for i in range(0, self.cells_per_axis[0]):

                # Loop over all cells in y direction.
                for j in range(0, self.cells_per_axis[1]):

                    # Evaluate cell index.
                    m = i + self.cells_per_axis[0]*j

                    # Resize particle and nearest neighbour lists.
                    self.particles[m] = [0] * self.num_neighbours
                    self.neighbours[m] = [0] * self.num_neighbours

                    # Zero particle tally.
                    self.tally[m] = 0

                    # Set cell index
                    self.index[m] = m

                    # Initalise nearest neighbour counter
                    nn_count = 0

                    # x loop for nearest neighbours
                    for a in range(0, 3):

                        x = (i + (a-1) + self.cells_per_axis[0]) % self.cells_per_axis[0]

                        # y loop for nearest neighbours
                        for b in range(0, 3):

                            y = (j + (b-1) + self.cells_per_axis[1]) % self.cells_per_axis[1]

                            # Nearest neighbour cell index.
                            nn = x + y*self.cells_per_axis[0]

                            self.neighbours[m][nn_count] = nn
                            nn_count += 1

    def reset(self):
        """ Reset cell tally counters. """
        self.tally = [0] * self.num_cells

    def get_cell(self, particle):
        """ Get a particle's cell index. """

        # Work out x and y indices.
        cellx = int(particle.position[0] / self.cell_spacing[0])
        celly = int(particle.position[1] / self.cell_spacing[1])

        # Compute the cell index.
        cell = cellx + celly * self.cells_per_axis[0]

        if self.dimension is 3:
            cellz = int(particle.position[2] / self.cell_spacing[2])
            cell += cellz * self.cells_per_axis[0] * self.cells_per_axis[1]

        return cell

    def init_cell(self, new_cell, particle):
        """ Insert a particle into the cell list. """

        # Add to new list.
        self.particles[new_cell][self.tally[new_cell]] = particle.index
        particle.cell = new_cell
        particle.pos_cell = self.tally[new_cell]
        self.tally[new_cell] += 1

        if self.tally[new_cell] == self.max_particles:
            raise ValueError('Maximum number of particles per cell exceeded.')
            raise SystemExit

    def init_cell_list(self, particles):
        """ Insert all particles into the cell list. """
        [init_cell(get_cell(particle), particle) for particle in particles]

    def update_cell(self, new_cell, particle, particles):
        """ Update the cell list. """

        # Remove particle from the old list.
        self.tally[particle.cell] -= 1
        self.particles[particle.cell][particle.pos_cell] = self.particles[particle.cell][self.tally[particle.cell]]
        particles[self.particles[particle.cell][self.tally[particle.cell]]].pos_cell = particle.pos_cell

        # Add particle to the new cell list
        self.particles[new_cell][self.tally[new_cell]] = particle.index
        particle.cell = new_cell
        particle.pos_cell = self.tally[new_cell]
        self.tally[new_cell] += 1

        if self.tally[new_cell] == self.max_particles:
            raise ValueError('Maximum number of particles per cell exceeded.')
            raise SystemExit
