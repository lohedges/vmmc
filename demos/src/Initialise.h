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

#ifndef _INITIALISE_H
#define _INITIALISE_H

/*! \file Initialise.h
*/

// FORWARD DECLARATIONS

class  Box;
class  CellList;
struct Particle;
class  MersenneTwister;

//! Class for the initialisation of particle configurations.
class Initialise
{
public:
    //! Default constructor.
    Initialise();

    //! Initialise a random particle configuration.
    /*! \param particles
            A reference to a vector of particles.

        \param cells
            A reference to the cell list container.

        \param box
            A reference to the simulation box.

        \param rng
            A reference to the random number generator.

        \param isSpherocylinder
            Whether particles are confined to a sphereocyliner.
     */
    void random(std::vector<Particle>&, CellList&, Box&, MersenneTwister&, bool);

    //! Check whether particle is within spherocylinder.
    /*! \param index
            The particle index.

        \param position
            The position of the particle following the virtual move.

        \param orientation
            The orientation of the particle following the virtual move.

        \return
            Whether the particle lies outside of the spherocylinder.
    */
#ifndef ISOTROPIC
    bool outsideSpherocylinder(unsigned int, const double*, const double*);
#else
    bool outsideSpherocylinder(unsigned int, const double*);
#endif

private:
    /// Copy of the simulation box size.
    std::vector<double> boxSize;

    //! Helper function for testing particle insertions.
    /*! \param particle
            A reference to the trial particle.

        \param particles
            A reference to the particle list.

        \param cells
            A refernce to the cell list.

        \param box
            A reference to the simulation box.
     */
    bool checkOverlap(Particle&, std::vector<Particle>&, CellList&, Box&);

    /// Maximum number of trial particle insertions (per particle).
    static const unsigned int MAX_TRIALS = 100000000;
};

#endif  /* _INITIALISE_H */
