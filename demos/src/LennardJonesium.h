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

#include "CellList.h"
#include "Model.h"

/*! \file LennardJonesium.h
*/

// Global infinity constant for hard core repulsions.
extern double INF;

//! Class defining the Lennard-Jones potential.
class LennardJonesium : public vmmc::Model
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
    LennardJonesium(Box&, std::vector<Particle>&, CellList&, unsigned int, double, double);

    //! Calculate the total interaction energy felt by a particle.
    /*! \param particle
            The particle index.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the first particle.

        \return
            The total interaction energy.
    */
#ifndef ISOTROPIC
    double energyCallback(unsigned int, double[], double[]);
#else
    double energyCallback(unsigned int, double[]);
#endif

    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param particle2
            The index of the first particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \return
            The pair energy between particles 1 and 2.
    */
#ifndef ISOTROPIC
    double pairEnergyCallback(unsigned int, double[], double[], unsigned int, double[], double[]);
#else
    double pairEnergyCallback(unsigned int, double[], unsigned int, double[]);
#endif

    //! Determine the interactions for a given particle.
    /*! \param particle
            The particle index.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the particle.

        \param interactions
            An array to store the indices of neighbours with which the particle interacts.

        \return
            The number of interactions.
    */
#ifndef ISOTROPIC
    unsigned int interactionsCallback(unsigned int, double[], double[], unsigned int[]);
#else
    unsigned int interactionsCallback(unsigned int, double[], unsigned int[]);
#endif

    //! Apply any post-move updates for a given particle.
    /*! \param particle
            The particle index.

        \param position
            The position of the particle following the  move.

        \param orientation
            The orientation of the particle following the  move.
    */
#ifndef ISOTROPIC
    void postMoveCallback(unsigned int, double[], double[]);
#else
    void postMoveCallback(unsigned int, double[]);
#endif

    //! Check for non-pairwise energy contributions.
    /*! \param particle
            The particle index.

        \param position
            The position of the particle.

        \param orientation
            The orientation of the particle.

        \param return
            The total non-pairwise energy felt by the particle.
    */
#ifndef ISOTROPIC
    double nonPairwiseCallback(unsigned int, double[], double[]);
#else
    double nonPairwiseCallback(unsigned int, double[]);
#endif

    //! Check custom boundary condition.
    /*! \param particle
            The particle index.

        \param position
            The position of the particle following the virtual move.

        \param orientation
            The orientation of the particle following the virtual move.

        \param return
            Whether the particle is outside the boundary.
    */
#ifndef ISOTROPIC
    bool boundaryCallback(unsigned int, double[], double[]);
#else
    bool boundaryCallback(unsigned int, double[]);
#endif

    //! Calculate global energy.
    double getEnergy();

private:
    Box& box;                           //!< A reference to the simulation box.
    std::vector<Particle>& particles;   //!< A reference to the particle list.
    CellList& cells;                    //!< A reference to the cell list.

    unsigned int maxInteractions;       //!< The maximum number of interactions per particle.
    double interactionEnergy;           //!< Interaction energy scale (in units of kBT).
    double interactionRange;            //!< Size of interaction range (in units of particle diameter).
    double squaredCutOffDistance;       //!< Squared cut-off distance.
    double potentialShift;              //!< Shift factor to zero potential at cut-off..
};

#endif  /* _LENNARDJONESIUM_H */
