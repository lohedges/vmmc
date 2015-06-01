/*
  Copyright (c) 2015 Lester Hedges <lester.hedges+vmmc@gmail.com>

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

#ifndef _LENNARDJONESIUM_H
#define _LENNARDJONESIUM_H

#include <cstdlib>
#include <limits>

#include "Box.h"
#include "CellList.h"
#include "Potential.h"

/*! \file LennardJonesium.h
*/

//! Class defining the square-well potential.
class LennardJonesium : public Potential
{
public:
    //! Constructor.
    /*! \param box_
            A reference to the simulation box object.

        \param particles_
            A reference to the particle list.

        \param cells_
            A reference to the cell list object.

        \param maxInteractions_
            The maximum number of interactions per particle.

        \param interactionEnergy_
            The square well interaction energy (in units of kBT).

        \param interactionRange_
            The square well interaction range (in units of the particle diameter).
     */
	LennardJonesium(Box&, std::vector <Particle>&, CellList&, unsigned int, double, double);

    //! Calculate the total interaction energy felt by a particle.
    /*! \param index
            The particle index.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the first particle.

        \return
            The total interaction energy.
     */
    double computeEnergy(unsigned int, double[], double[]);

    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \return
            The pair energy between particles 1 and 2.
     */
    double computePairEnergy(unsigned int, double[], double[], unsigned int, double[], double[]);

private:
    double potentialShift;  //!> shift factor to zero potential at cut-off
};

#endif	/* _LENNARDJONESIUM_H */
