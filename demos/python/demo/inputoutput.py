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
@brief Class for reading and writing particle configurations.
"""

import math
import sys

class InputOutput:
    """ Class for reading and writing particle configurations. """

    def load_configuration(self, file_name, box, particles, cells, is_isotropic):
        """ Load particle configuration from file. """

        # Try to open the restart file.
        try:
            f = open(file_name, 'w')
        except IOError:
            print >> sys.stderr, "I/O Error:", file_name, "doesn't exist!"
            raise SystemExit

        # Read data.
        data = [line.strip() for line in f]

        # Check configuration size is consistent.
        if len(data) != len(particles):
            raise ValueError('Particle configuration is too large.')
            raise SystemExit

        # Reset the cell list.
        cells.reset()

        # Process each particle from the configuration file.
        for x, particle in enumerate(data):

            # Resize particle position and orientation vectors.
            particles[x].position = [0] * box.dimension
            particles[x].orientation = [0] * box.dimension

            # Split data string.
            p = particle.split()

            # Read data for each dimension.
            for d in range(0, box.dimension):

                # Position
                particles[x].position[d] = float(p[d])

                if is_isotropic:
                    # Assign dummy orientation.
                    particles[x].orientation[d] = 1 / math.sqrt(box.dimension)
                else:
                    particles[x].orientation[d] = float(p[d + box.dimension])

            # Enforce periodic boundary conditions.
            box.periodic_boundaries(particles[x].position)

            # Calculate the particle's cell index.
            particles[x].cell = cells.get_cell(particles[x])

            # Update the cell list.
            cells.init_cell(particles[x].cell, particles[x])

        # Close the input file.
        f.close()

    def save_configuration(self, file_name, box, particles, is_isotropic):
        """ Save particle configuration to file. """

        # Open file for writing.
        f = open(file_name, 'w')

        # Loop over all particles.
        for p in particles:

            # Write particle position.
            f.write("%5.4f %5.4f" % (p.position[0], p.position[1]))
            if box.dimension == 3:
                f.write(" %5.4f" % p.position[2])

            # Write particle orientation.
            if not is_isotropic:
                f.write(" %5.4f %5.4f" % (p.orientation[0], p.orientation[1]))
                if box.dimension == 3:
                    f.write(" %5.4f" % p.orientation[2])

            # Terminate line.
            f.write('\n')

        # Close the output file.
        f.close()

    def append_xyz_trajectory(self, dimension, particles, clear_file):
        """ Append particle configuration to a VMD xyz trajectory file. """

        # Wipe existing trajectory file.
        if clear_file is True:
            f = open("trajectory.xyz", 'w')
            f.close()

        f = open("trajectory.xyz", 'a')
        f.write("%d\n\n" % len(particles));

        # Write configuration.
        for p in particles:
            f.write("0 %5.4f %5.4f %5.4f\n" % (p.position[0], p.position[1], p.position[2] if dimension == 3 else 0))

        # Close the output file.
        f.close()

    def vmd_script(self, box_size):
        """ Write a VMD TcL script for setting view and particle attributes. """

        f = open("vmd.tcl", 'w')

        # Turn on lights 0 and 1
        f.write("light 0 on\n")
        f.write("light 1 on\n")
        f.write("light 2 off\n")
        f.write("light 3 off\n")

        # Position the stage and axes.
        f.write("axes location off\n")
        f.write("stage location off\n")

        # Set orthographic projection.
        f.write("display projection orthographic\n")

        # Set drawing method to van der Waals radius.
        f.write("mol modstyle 0 0 VDW 1 30\n")

        # Set sensible atom radius.
        f.write("set sel [atomselect top \"name X\"]\n")
        f.write("atomselect0 set radius 0.5\n")

        # Set default particle to blue.
        f.write("color Name X blue\n")

        # Turn off depth cue.
        f.write("display depthcue off\n")

        # Define box boundaries.
        f.write("set minx 0\n")
        f.write("set maxx %5.4f\n" % box_size[0])
        f.write("set miny 0\n")
        f.write("set maxy %5.4f\n" % box_size[1])
        if len(box_size) == 3:
            f.write("set minz 0\n")
            f.write("set maxz %5.4f\n" % box_size[2])
        else:
            f.write("set minz 0\n")
            f.write("set maxz 0\n")

        # Set colours.
        f.write("draw materials off\n")
        f.write("draw color white\n")

        # Draw cube edges.
        f.write("draw line \"$minx $miny $minz\" \"$maxx $miny $minz\"\n")
        f.write("draw line \"$minx $miny $minz\" \"$minx $maxy $minz\"\n")
        f.write("draw line \"$minx $miny $minz\" \"$minx $miny $maxz\"\n")
        f.write("draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n")
        f.write("draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\"\n")
        f.write("draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\"\n")
        f.write("draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\"\n")
        f.write("draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\"\n")
        f.write("draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\"\n")
        f.write("draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\"\n")
        f.write("draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\"\n")
        f.write("draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\"\n")

        # Rotate box.
        if len(box_size) == 3:
            f.write("rotate x by -60\n")
            f.write("rotate y by -30\n")
            f.write("rotate z by -15\n")

        f.close();
