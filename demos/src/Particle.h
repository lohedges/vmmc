/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>

/*! \file Particle.h
    \brief A simple particle data type.
*/

//! Structure containing attributes for an individual particle.
struct Particle
{
public:
    //! Default constructor.
    Particle();

    unsigned int index;                 //!< Particle index.
    std::vector<double> position;       //!< The x,y,z coordinates of particle.
    std::vector<double> orientation;    //!< The orientation of the particle (unit vector).

    unsigned int cell;                  //!< Index of the cell in which the particle is located.
    unsigned int posCell;               //!< Position of particle in the corresponding cell list.
};

#endif  /* _PARTICLE_H */
