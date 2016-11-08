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

#ifndef _MODEL_H
#define _MODEL_H

/*! \file Model.h
*/

namespace vmmc
{
    //! Pure abstract base class defining general access functions for the model potential.
    class Model
    {
    public:
        //! Constructor.
        Model();

        //! Destructor.
        ~Model();

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
        virtual double energyCallback(unsigned int, const double*, const double*);
#else
        virtual double energyCallback(unsigned int, const double*);
#endif

        //! Calculate the pair energy between two particles.
        /*! \param particle1
                The index of the first particle.

            \param position1
                The position vector of the first particle.

            \param orientation1
                The orientation vector of the first particle.

            \param particle2
                The index of the second particle.

            \param position2
                The position vector of the second particle.

            \param orientation2
                The orientation vector of the second particle.

            \return
                The pair energy between particles 1 and 2.
        */
#ifndef ISOTROPIC
        virtual double pairEnergyCallback(unsigned int, const double*, const double*, unsigned int, const double*, const double*);
#else
        virtual double pairEnergyCallback(unsigned int, const double*, unsigned int, const double*);
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
        virtual unsigned int interactionsCallback(unsigned int, const double*, const double*, unsigned int*);
#else
        virtual unsigned int interactionsCallback(unsigned int, const double*, unsigned int*);
#endif

        //! Apply any post-move updates for a given particle.
        /*! \param particle
                The particle index.

            \param position
                The position of the particle following the virtual move.

            \param orientation
                The orientation of the particle following the virtual move.
        */
#ifndef ISOTROPIC
        virtual void postMoveCallback(unsigned int, const double*, const double*);
#else
        virtual void postMoveCallback(unsigned int, const double*);
#endif

        //! Check for non-pairwise energy contributions.
        /*! \param particle
                The particle index.

            \param position
                The position of the particle.

            \param orientation
                The orientation of the particle.

            \return
                The total non-pairwise energy felt by the particle.
        */
#ifndef ISOTROPIC
        virtual double nonPairwiseCallback(unsigned int, const double*, const double*);
#else
        virtual double nonPairwiseCallback(unsigned int, const double*);
#endif

        //! Check custom boundary condition.
        /*! \param particle
                The particle index.

            \param position
                The position of the particle following the virtual move.

            \param orientation
                The orientation of the particle following the virtual move.

            \return
                Whether the particle is outside the boundary.
        */
#ifndef ISOTROPIC
        virtual bool boundaryCallback(unsigned int, const double*, const double*);
#else
        virtual bool boundaryCallback(unsigned int, const double*);
#endif
    };
}

#endif  /* _MODEL_H */
