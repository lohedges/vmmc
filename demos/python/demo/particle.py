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
