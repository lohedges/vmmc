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

#ifndef _VMMC_H
#define _VMMC_H

#include <functional>
#include <vector>

#include "MersenneTwister.h"

/*! \file VMMC.h
    \brief A simple class for executing Virtual Move Monte Carlo moves (cluster translations and rotations).
    For algorithmic details, see:
    -# Avoiding unphysical kinetic traps in Monte Carlo simulations of strongly attractive particles,\n
    S. Whitelam and P.L. Geissler,
    <a href="http://dx.doi.org/10.1063/1.2790421"> Journal of Chemical Physics, 127, 154101 (2007)</a>.
    -# Approximating the dynamical evolution of systems of strongly interacting overdamped particles,\n
    S. Whitelam, <a href="http://dx.doi.org/10.1080/08927022.2011.565758">
    Molecular Simulation, 37 (7) (2011)</a>. (Preprint version available
    <a href= "http://arxiv.org/abs/1009.2008">here</a>.)
*/

namespace vmmc
{
    // CALLBACK FUNCTION PROTOTYPES

    //! Calculate the energy for a given particle.
    /*! \param index
            The particle index.

        \param position
            The position of the particle.

        \param orientation
            The orientation of the particle.

        \return
            The total interaction energy felt by the particle.
    */
#ifndef ISOTROPIC
    typedef std::function<double (unsigned int, const double*, const double*)> EnergyCallback;
#else
    typedef std::function<double (unsigned int, const double*)> EnergyCallback;
#endif

    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position of the first particle.

        \param orientation1
            The orientation of the first particle.

        \param particle2
            The index of the second particle.

        \param position2
            The position of the second particle.

        \param orientation2
            The orientation of the second particle.

        \return
            The pair interaction energy between particles 1 and 2.
    */
#ifndef ISOTROPIC
    typedef std::function<double (unsigned int, const double*, const double*, unsigned int, const double*, const double*)> PairEnergyCallback;
#else
    typedef std::function<double (unsigned int, const double*, unsigned int, const double*)> PairEnergyCallback;
#endif

    //! Determine the interactions for a particle.
    /*! \param index
            The particle index.

        \param position
            The position of the particle.

        \param orientation
            The orientation of the particle.

        \param interactions
            An array to store the indices of neighbours with which the particle interacts.

        \return
            The number of interactions.
    */
#ifndef ISOTROPIC
    typedef std::function<unsigned int (unsigned int, const double*, const double*, unsigned int[])> InteractionsCallback;
#else
    typedef std::function<unsigned int (unsigned int, const double*, unsigned int[])> InteractionsCallback;
#endif

    //! Apply any post-move updates for a given particle.
    /*! \param index
            The particle index.

        \param position
            The position of the particle following the virtual move.

        \param orientation
            The orientation of the particle following the virtual move.
    */
#ifndef ISOTROPIC
    typedef std::function<void (unsigned int, const double*, const double*)> PostMoveCallback;
#else
    typedef std::function<void (unsigned int, const double*)> PostMoveCallback;
#endif

    //! Calculate the non-pairwise energy felt by a particle.
    /*! \param index
            The particle index.

        \param position
            The position of the particle following.

        \param orientation
            The orientation of the particle.

        \return
            The total non-pairwise energy felt by the particle.
    */
#ifndef ISOTROPIC
    typedef std::function<double (unsigned int, const double*, const double*)> NonPairwiseCallback;
#else
    typedef std::function<double (unsigned int, const double*)> NonPairwiseCallback;
#endif

    //! Check custom boundary condition.
    /*! \param index
            The particle index.

        \param position
            The position of the particle following the virtual move.

        \param orientation
            The orientation of the particle following the virtual move.

        \return
            Whether the particle lies outside the custom boundary.
    */
#ifndef ISOTROPIC
    typedef std::function<bool (unsigned int, const double*, const double*)> BoundaryCallback;
#else
    typedef std::function<bool (unsigned int, const double*)> BoundaryCallback;
#endif

    // DATA TYPES

    //! Container for storing virtual move parameters.
    struct Parameters
    {
        unsigned int seed;                          //!< Index of the seed particle.
        bool isRotation;                            //!< Whether the move is a rotation.
        double stepSize;                            //!< The magnitude of the trial move.
        std::vector<double> trialVector;            //!< Vector for trial move.
    };

    //! Container for storing particle attributes during the virtual move.
    class Particle
    {
    public:
        //! Default constructor.
        Particle();

        //! Constructor.
        /*! \param dimension
                The number of dimensions.
        */
        Particle(unsigned int);

        unsigned int index;                         //!< Particle index.
        bool isMoving;                              //!< Whether the particle is part of the virtual move.
        bool isFrustrated;                          //!< Whether the particle is involved in a frustrated link.
        unsigned int posFrustated;                  //!< Index in the frustrated links array.
        std::vector<double> preMovePosition;        //!< Particle position before the virtual move.
        std::vector<double> postMovePosition;       //!< Particle position following the virtual move.
        std::vector<double> clusterPosition;        //!< Position of the particle in the moving cluster (relative to seed).
#ifndef ISOTROPIC
        std::vector<double> preMoveOrientation;     //!< Particle orientation before the virtual move.
        std::vector<double> postMoveOrientation;    //!< Particle orientation following the virtual move.
#endif
    };

    //! Container for storing callback functions
    struct CallbackFunctions
    {
        EnergyCallback energyCallback;              //!< Callback function to calculate particle energies.
        PairEnergyCallback pairEnergyCallback;      //!< Callback function to calculate pair energies.
        InteractionsCallback interactionsCallback;  //!< Callback function to determine particle interactions.
        PostMoveCallback postMoveCallback;          //!< Callback function to apply any post-move updates.
        NonPairwiseCallback nonPairwiseCallback;    //!< Callback function to calculate non-pairwise interaction energies.
        BoundaryCallback boundaryCallback;          //!< Callback function to apply custom boundary conditions.

        bool isNonPairwise;                         //!< Whether the non-pairwise energy callback is defined.
        bool isCustomBoundary;                      //!< Whether the boundary callback is defined.
    };

    //! Main VMMC class.
    class VMMC
    {
    public:
        //! Constructor.
        /*! \param nParticles_
                The number of particles in the simulation box.

            \param dimension_
                The dimension of the simulation box.

            \param coordinates
                The coordinates of all particles in the system.

            \param orientations
                The orientations of all particle in the system.

            \param maxTrialTranslation_
                The maximum trial translation (in units of the reference particle diameter).

            \param maxTrialRotation_
                The maximum trial rotation.

            \param probTranslate_
                The probability of performing a translation move (versus a rotation).

            \param referenceRadius_
                Reference particle radius (for Stokes scaling).

            \param maxInteractions_
                Maximum number of interactions per particle.

            \param boxSize_
                The size of the periodic simulation box in each dimension.

            \param isIsotropic_
                Whether the potential of each particle is isotropic.

            \param isRepusive_
                Whether there are finite repulsive interactions.

            \param callbacks_
                Callback function container.
        */
#ifndef ISOTROPIC
        VMMC(unsigned int, unsigned int, double*, double*, double, double, double, double, unsigned int, double*, bool*, bool,
#else
        VMMC(unsigned int, unsigned int, double*, double, double, double, double, unsigned int, double*, bool,
#endif
            const CallbackFunctions&);

        //! Overloaded ++ operator. Perform a single VMMC step.
        void operator ++ (const int);

        //! Overloaded += operator. Perform "n" VMMC steps.
        void operator += (const int);

        //! Perform a single VMMC trial move.
        void step();

        //! Perform a specified number of VMMC trial moves.
        /*! \param nSteps
                The number of attempted VMMC trial moves.
        */
        void step(const int);

        //! Get the number of attempted moves.
        /*! \return
                The number of attempted virtual moves.
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

        //! Get the number of accepted translation moves for each cluster size.
        /*! \param clusterStatistics
                An array into which the cluster statistics will be copied.
        */
        void getClusterTranslations(unsigned long long[]) const;

        //! Get the number of accepted translation moves for each cluster size.
        /*! \return
                A const reference to the cluster statistics vector.
        */
        const std::vector<unsigned long long>& getClusterTranslations() const;

        //! Get the number of accepted rotation moves for each cluster size.
        /*! \param clusterStatistics
                An array into which the cluster statistics will be copied.
        */
        void getClusterRotations(unsigned long long[]) const;

        //! Get the number of accepted rotation moves for each cluster size.
        /*! \return
                A const reference to the cluster statistics vector.
        */
        const std::vector<unsigned long long>& getClusterRotations() const;

        //! Reset statistics.
        void reset();

        MersenneTwister rng;                        //!< Random number generator.

    private:
        Parameters moveParams;                      //!< Parameters for the trial move.
        unsigned long long nAttempts;               //!< Number of attempted moves.
        unsigned long long nAccepts;                //!< Number of accepted moves.
        unsigned long long nRotations;              //!< Number of accepted rotations.

        unsigned int nParticles;                    //!< The number of particles in the simulation box.
        unsigned int dimension;                     //!< The dimension of the simulation box.
        double maxTrialTranslation;                 //!< The maximum trial translation (in units of the reference diameter).
        double maxTrialRotation;                    //!< The maximum trial rotation.
        double probTranslate;                       //!< The relative probability of translational moves (vs rotations).
        double referenceRadius;                     //!< Reference particle radius (for Stokes scaling).
        unsigned int maxInteractions;               //!< Maximum number of interactions per particle.
        std::vector<double> boxSize;                //!< The size of the simulation box in each dimension.
#ifndef ISOTROPIC
        std::vector<bool> isIsotropic;              //!< Whether the potential of each particle is isotropic.
#endif
        bool isRepusive;                            //!< Whether there are finite repulsive interactions.
        bool is3D;                                  //!< Whether the simulation is three-dimensional.

        CallbackFunctions callbacks;                //!< Callback functions.

        std::vector<Particle> particles;            //!< Vector of particles.

        unsigned int nMoving;                                   //!< The number of particles in the cluster.
        std::vector<unsigned int> moveList;                     //!< the indices of particles in the cluster.
        std::vector<unsigned long long> clusterTranslations;    //!< Array for storing the number of translations for each cluster size.
        std::vector<unsigned long long> clusterRotations;       //!< Array for storing the number of rotations for each cluster size

        unsigned int nFrustrated;                               //!< The number of frustrated links.
        std::vector<unsigned int> frustratedLinks;              //!< Array of particles involved in frustrated links.

        unsigned int nInteractions;                             //!< The number of pair interactions for particles in the cluster.
        std::vector<std::vector<unsigned int> > interactions;   //!< Indices of particle pairs that interact in the cluster.
        std::vector<std::vector<double> > pairEnergyMatrix;     //!< Pair energies for particle interactions in the cluster.

        unsigned int cutOff;                        //!< The cut-off cluster size for the trial move.
        bool isEarlyExit;                           //!< Whether trial move aborted early.

        //! Propose a trial particle translation/rotation.
        void proposeMove();

        //! Determine whether move is accepted.
        bool accept();

        //! Compute the hydrodynamic radius of the moving cluster.
        double computeHydrodynamicRadius() const;

        //! Compute particle's position and orientation following the trial move.
        /*! \param particle
                Index of the particle.

            \param direction
                Whether move is forward (1) or reverse (-1).

            \param postMoveParticle
                The particle data structure.
        */
        void computePostMoveParticle(unsigned int, int, Particle&);

        //! Initiate a particle ready for the virtual move.
        /*! \param particle
                Index of the particle.

            \param linker
                A reference to the linking particle.
        */
        void initiateParticle(unsigned int, Particle&);

        //! Recursively assign additional particles to the moving cluster.
        /*! \param particle
                Index of the trial particle.
        */
        void recursiveMoveAssignment(unsigned int);

        //! Apply/unnapply the virtual move.
        void swapMoveStatus();

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

        //! Calculate the minimum image separation between two coordinates (from v1 to v2).
        /*! \param v1
                The coordinate vector of the first particle.

            \param v2
                The coordinate vector of the second particle.

            \param sep
                The minimum image separation vector.
        */
        void computeSeparation(std::vector<double>&, std::vector<double>&, std::vector<double>&);

        //! Enforce periodic boundary conditions.
        /*! \param vec
                The coordinate vector.
        */
        void applyPeriodicBoundaryConditions(std::vector<double>&);

        //! Compute the norm of a vector.
        /*! \param vec
                A reference to the vector

            \return
                The norm of the vector.
        */
        double computeNorm(std::vector<double>&);
    };
}

#endif /* _VMMC_H */
