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

#ifndef _SQUAREWELLIUMWALL_H
#define _SQUAREWELLIUMWALL_H

#include "Model.h"

/*! \file SquareWelliumWall.h
*/

//! Class defining a square-well potential plus a wall.
/*!
    The wall is placed at a distance of -0.5 below the bottom of the box,
    i.e. half a particle diameter. The wall is in the y dimension in 2D,
    and the z dimension in 3D.
 */
class SquareWelliumWall : public Model
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

        \param wallInteractionEnergy_
            The interaction energy between particles and the wall (in units of kBT).

        \param wallInteractionRange_
            The interaction range between particles and the wall (in units of the particle diameter).
     */
    SquareWelliumWall(Box&, std::vector<Particle>&, CellList&, unsigned int, double, double, double, double);

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
#ifndef ISOTROPIC
    double computePairEnergy(unsigned int, double[], double[], unsigned int, double[], double[]);
#else
    double computePairEnergy(unsigned int, double[], unsigned int, double[]);
#endif

    //! Calculate the interaction energy between a particle and the wall.
    /*! \param particle
            The index of the particle.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the particle.

        \return
            The pair interaction energy between the particle and the wall.
     */
#ifndef ISOTROPIC
    double computeWallEnergy(unsigned int, double[], double[]);
#else
    double computeWallEnergy(unsigned int, double[]);
#endif

    //! Test whether a particle moves outside of the non-periodic boundaries.
    /*! \param particle
            The index of the particle.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the particle.

        \return
            The pair interaction energy between the particle and the wall.
     */
#ifndef ISOTROPIC
    bool isOutsideBoundary(unsigned int, double[], double[]);
#else
    bool isOutsideBoundary(unsigned int, double[]);
#endif

private:
    /// The interaction energy between particles and the wall.
    double wallInteractionEnergy;

    /// The range of interaction between particles and the wall.
    double wallInteractionRange;
};

#endif  /* _SQUAREWELLIUMWALL_H */
