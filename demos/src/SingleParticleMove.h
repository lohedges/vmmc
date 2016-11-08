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

#ifndef _SINGLEPARTICLEMOVE_H
#define _SINGLEPARTICLEMOVE_H

#include <vector>

#include "MersenneTwister.h"

/*! \file SingleParticleMove.h
    \brief A class for executing single particle Move Monte Carlo moves (translations and rotations).
*/

// FORWARD DECLARATIONS

class  Model;
struct Particle;

//! Container for storing move parameters.
struct MoveParams
{
    unsigned int seed;                          //!< Index of the seed particle.
    bool isRotation;                            //!< Whether the move is a rotation.
    double stepSize;                            //!< The magnitude of the trial move.
    std::vector<double> trialVector;            //!< Vector for trial move.
    Particle preMoveParticle;                   //!< Particle state before the trial move.
};

class SingleParticleMove
{
public:
    //! Constructor.
    /*! \param model_
            A pointer to the model object.

        \param maxTrialTranslation_
            The maximum trial translation (in units of the reference particle diameter).

        \param maxTrialRotation_
            The maximum trial rotation.

        \param probTranslate_
            The probability of performing a translation move (versus a rotation).

        \param isIsotropic_
            Whether the potential is isotropic.
     */
    SingleParticleMove(Model*, double, double, double, bool);

    //! Overloaded ++ operator. Perform a single step.
    void operator ++ (const int);

    //! Overloaded += operator. Perform "n" steps.
    void operator += (const int);

    //! Perform a single trial move.
    void step();

    //! Perform a specified number of trial moves.
    /*! \param nSteps
            The number of attempted trial moves.
     */
    void step(const int);

    //! Get the number of attempted moves.
    /*! \return
            The number of attempted moves.
     */
    unsigned long long getAttempts() const;

    //! Get the number of accepted moves.
    /*! \return
            The number of accepted virtual moves.
     */
    unsigned long long getAccepts() const;

    //! Get the number of accepted rotation moves.
    /*! \return
            The number of accepted rotation moves.
     */
    unsigned long long getRotations() const;

    //! Reset statistics.
    void reset();

    MersenneTwister rng;                        //!< Random number generator.

private:
    MoveParams moveParams;                      //!< Parameters for the trial move.
    Model* model;                               //!< A pointer to the model object.
    unsigned long long nAttempts;               //!< Number of attempted moves.
    unsigned long long nAccepts;                //!< Number of accepted moves.
    unsigned long long nRotations;              //!< Number of accepted rotations.

    double maxTrialTranslation;                 //!< The maximum trial translation (in units of the reference diameter).
    double maxTrialRotation;                    //!< The maximum trial rotation.
    double probTranslate;                       //!< The relative probability of translational moves (vs rotations).
    bool is3D;                                  //!< Whether the simulation is three-dimensional.
    bool isIsotropic;                           //!< Whether the potential is isotropic.
    double energyChange;                        //!< Energy change resulting from trial move.

    //! Propose a trial particle translation/rotation.
    void proposeMove();

    //! Determine whether move is accepted.
    bool accept();

    //! Calculate an unbiased rotation vector in 3D (Beard & Schlick, BJ 85 2973 (2003)).
    /*! \param v1
            The vector about which to rotate (either the position or orientation).

        \param v2
            Rotation unit vector.

        \param v3
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate3D(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);

    //! Calculate a simple in plane rotatation vector.
    /*! \param v1
            The vector to rotate (either the position or orientation).

        \param v2
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate2D(std::vector<double>&, std::vector<double>&, double);

    //! Compute the norm of a vector.
    /*! \param vec
            A reference to the vector

        \return
            The norm of the vector.
     */
    double computeNorm(std::vector<double>&);
};

#endif  /* _SINGLEPARTICLEMOVE_H */
