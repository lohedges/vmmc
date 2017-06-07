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

#include <cstdlib>
#include <iostream>

#include "VMMC.h"

namespace vmmc
{
    Particle::Particle() {}

    Particle::Particle(unsigned int dimension)
    {
        // Resize position/orientation vectors.
        preMovePosition.resize(dimension);
        postMovePosition.resize(dimension);
        clusterPosition.resize(dimension);
#ifndef ISOTROPIC
        preMoveOrientation.resize(dimension);
        postMoveOrientation.resize(dimension);
#endif
    }

    VMMC::VMMC(
        unsigned int nParticles_,
        unsigned int dimension_,
        double* coordinates,
#ifndef ISOTROPIC
        double* orientations,
#endif
        double maxTrialTranslation_,
        double maxTrialRotation_,
        double probTranslate_,
        double referenceRadius_,
        unsigned int maxInteractions_,
        double* boxSize_,
#ifndef ISOTROPIC
        bool* isIsotropic_,
#endif
        bool isRepusive_,
        const CallbackFunctions& callbacks_) :

        nAttempts(0),
        nAccepts(0),
        nRotations(0),
        nParticles(nParticles_),
        dimension(dimension_),
        maxTrialTranslation(maxTrialTranslation_),
        maxTrialRotation(maxTrialRotation_),
        probTranslate(probTranslate_),
        referenceRadius(referenceRadius_),
        maxInteractions(maxInteractions_),
        isRepusive(isRepusive_),
        callbacks(callbacks_)
    {
        // Check number of particles.
        if ((nParticles == 0) ||
            (nParticles > (1 + std::numeric_limits<unsigned int>::max() - nParticles)))
        {
            std::cerr << "[ERROR] VMMC: Number of particle must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check dimensionality.
        if (dimension == 3) is3D = true;
        else if (dimension == 2) is3D = false;
        else
        {
            std::cerr << "[ERROR] VMMC: Invalid dimensionality!\n";
            exit(EXIT_FAILURE);
        }

        // Check maximum trial translation.
        if (maxTrialTranslation < 0)
        {
            std::cerr << "[ERROR] VMMC: Maximum trial translation must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check maximum trial rotation.
        if (maxTrialRotation < 0)
        {
            std::cerr << "[ERROR] VMMC: Maximum trial rotation must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check reference radius.
        if (referenceRadius < 0)
        {
            std::cerr << "[ERROR] VMMC: Reference radius must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // N.B. There's no need to check probTranslate since anything less than zero
        // will be treated as zero, and anything greater than one will be treated as one.

        // Store simulation box size.
        boxSize.resize(dimension);
        for (unsigned int i=0;i<dimension;i++)
        {
            boxSize[i] = boxSize_[i];

            // Check box size.
            if (boxSize[i] < 0)
            {
                std::cerr << "[ERROR] VMMC: Box length must be > 0!\n";
                exit(EXIT_FAILURE);
            }
        }

        // Allocate memory.
        moveParams.trialVector.resize(dimension);
        particles.resize(nParticles);
        moveList.resize(nParticles);
        clusterTranslations.resize(nParticles);
        clusterRotations.resize(nParticles);
        frustratedLinks.resize(nParticles);
#ifndef ISOTROPIC
        isIsotropic.resize(nParticles);
#endif

        // Create particle container.
        for (unsigned int i=0;i<nParticles;i++)
        {
            // Resize vectors.
            particles[i].preMovePosition.resize(dimension);
            particles[i].postMovePosition.resize(dimension);
            particles[i].clusterPosition.resize(dimension);
#ifndef ISOTROPIC
            particles[i].preMoveOrientation.resize(dimension);
            particles[i].postMoveOrientation.resize(dimension);
#endif

            // Initialise moving boolean flag.
            particles[i].isMoving = false;

            // Initialise frustrated boolean flag.
            particles[i].isFrustrated = false;

            // Copy particle coordinates and orientations.
            for (unsigned int j=0;j<dimension;j++)
            {
                particles[i].preMovePosition[j] = coordinates[dimension*i + j];
#ifndef ISOTROPIC
                particles[i].preMoveOrientation[j] = orientations[dimension*i + j];
#endif

                // Check coordinate.
                if ((particles[i].preMovePosition[j] < 0) ||
                    (particles[i].preMovePosition[j] > boxSize[j]))
                {
                    std::cerr << "[ERROR] VMMC: Coordinates must run from 0 to the box size!\n";
                    exit(EXIT_FAILURE);
                }
            }

#ifndef ISOTROPIC
            // Check that orientation is a unit vector.
            if (std::abs(1.0 - computeNorm(particles[i].preMoveOrientation)) > 1e-6)
            {
                std::cerr << "[ERROR] VMMC: Particle orientations must be unit vectors!\n";
                exit(EXIT_FAILURE);
            }
#endif

#ifndef ISOTROPIC
            // Store particle potential style.
            isIsotropic[i] = isIsotropic_[i];
#endif
        }

        // Allocate memory for pair interaction matrix (finite repulsions only).
        if (isRepusive)
        {
            // Maximum number of pair interactions.
            unsigned int nPairs = (nParticles*maxInteractions)/2;

            interactions.resize(nPairs);
            for (unsigned int i=0;i<nPairs;i++)
                interactions[i].resize(2);

            // Construct a triangular matrix to save memory.
            pairEnergyMatrix.resize(nParticles);
            for (unsigned int i=0;i<nParticles;i++)
                pairEnergyMatrix[i].resize(i);
        }

        // Check for non-pairwise energy callback function.
        if (callbacks.nonPairwiseCallback == nullptr) callbacks.isNonPairwise = false;
        else callbacks.isNonPairwise = true;

        // Check for custom boundary callback function.
        if (callbacks.boundaryCallback == nullptr) callbacks.isCustomBoundary = false;
        else callbacks.isCustomBoundary = true;

        std::cout << "Initialised VMMC";
#ifdef ISOTROPIC
        std::cout << " (isotropic)";
#endif
        std::cout << ".\nseed\t" << rng.getSeed() << '\n';
        // Print version info.
#ifdef COMMIT
        std::cout << "commit\t" << COMMIT << '\n';
#endif
        // Print branch info.
#ifdef BRANCH
        std::cout << "branch\t" << BRANCH << '\n';
#endif
    }

    void VMMC::step(const int nSteps)
    {
        for (int i=0;i<nSteps;i++)
            step();
    }

    void VMMC::operator ++ (const int)
    {
        step();
    }

    void VMMC::operator += (const int nSteps)
    {
        step(nSteps);
    }

    void VMMC::step()
    {
        // Increment number of attempted moves.
        nAttempts++;

        // Reset number of moving particles.
        nMoving = 0;

        // Reset number of frustrated links.
        nFrustrated = 0;

        // Reset number of pair interactions.
        nInteractions = 0;

        // Reset early exit flag.
        isEarlyExit = false;

        // Propose a move for the cluster.
        proposeMove();

        // Move hasn't been aborted.
        if (!isEarlyExit)
        {
            // Check for acceptance and apply move.
            if (accept())
            {
                // Increment number of accepted moves.
                nAccepts++;

                // Increment number of rotations.
                nRotations += moveParams.isRotation;

                // Tally cluster size.
                if (moveParams.isRotation) clusterRotations[nMoving-1]++;
                else clusterTranslations[nMoving-1]++;
            }
            else
            {
                // Undo move.
                if (!isEarlyExit) swapMoveStatus();
            }
        }

        // Reset the move list.
        for (unsigned int i=0;i<nMoving;i++) particles[moveList[i]].isMoving = false;

        // Reset frustrated links.
        for (unsigned int i=0;i<nFrustrated;i++) particles[frustratedLinks[i]].isFrustrated = false;

        // Reset pair interaction matrix.
        if (isRepusive)
        {
            for (unsigned int i=0;i<nInteractions;i++)
                pairEnergyMatrix[interactions[i][0]][interactions[i][1]] = 0;
        }
    }

    unsigned long long VMMC::getAttempts() const
    {
        return nAttempts;
    }

    unsigned long long VMMC::getAccepts() const
    {
        return nAccepts;
    }

    unsigned long long VMMC::getRotations() const
    {
        return nRotations;
    }

    void VMMC::getClusterTranslations(unsigned long long clusterStatistics[]) const
    {
        for (unsigned int i=0;i<nParticles;i++)
            clusterStatistics[i] = clusterTranslations[i];
    }

    const std::vector<unsigned long long>& VMMC::getClusterTranslations() const
    {
        return clusterTranslations;
    }

    void VMMC::getClusterRotations(unsigned long long clusterStatistics[]) const
    {
        for (unsigned int i=0;i<nParticles;i++)
            clusterStatistics[i] = clusterRotations[i];
    }

    const std::vector<unsigned long long>& VMMC::getClusterRotations() const
    {
        return clusterRotations;
    }

    void VMMC::reset()
    {
        nAttempts = nAccepts = nRotations = 0;
        std::fill(clusterTranslations.begin(), clusterTranslations.end(), 0);
        std::fill(clusterRotations.begin(), clusterRotations.end(), 0);
    }

    void VMMC::proposeMove()
    {
        // Choose a seed particle.
        moveParams.seed = rng.integer(0, nParticles-1);

        // Get a uniform random number in range [0-1].
        double r = rng();

        // Make sure the divisor doesn't blow things up.
        while (r == 0) r = rng();

        // Cluster size cut-off.
        cutOff = int(1.0/r);

        // Choose a random point on the surface of the unit sphere/circle.
        for (unsigned int i=0;i<dimension;i++)
            moveParams.trialVector[i] = rng.normal();

        // Normalise the trial vector.
        double norm = computeNorm(moveParams.trialVector);
        for (unsigned int i=0;i<dimension;i++)
            moveParams.trialVector[i] /= norm;

        // Neighbour index (for isotropic rotations).
        unsigned int neighbour;

        // Choose the move type.
        if (rng() < probTranslate)
        {
            // Translation.
            moveParams.isRotation = false;

            // Scale step-size to uniformly sample unit sphere/circle.
            if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
            else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
        }
        else
        {
            // Rotation.
            moveParams.isRotation = true;
            moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);

            // Check whether seed particle is isotropic.
#ifndef ISOTROPIC
            if (isIsotropic[moveParams.seed])
#endif
            {
                // Cluster size cut-off (minimum size is two).
                cutOff = int(2.0/r);

                unsigned int pairInteractions[maxInteractions];

                // Get a list of pair interactions.
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(moveParams.seed, &particles[moveParams.seed].preMovePosition[0],
                    &particles[moveParams.seed].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(moveParams.seed,
                    &particles[moveParams.seed].preMovePosition[0], pairInteractions);
#endif

                // Abort move if there are no neighbours, else choose one at random.
                if (nPairs == 0) isEarlyExit = true;
                else neighbour = pairInteractions[rng.integer(0, nPairs-1)];
            }
        }

        if (!isEarlyExit)
        {
            // Initialise the seed particle.
            particles[moveParams.seed].clusterPosition = particles[moveParams.seed].preMovePosition;
            initiateParticle(moveParams.seed, particles[moveParams.seed]);

            // Check that trial move of seed hasn't triggered early exit condition.
            if (!isEarlyExit)
            {
#ifndef ISOTROPIC
                if (isIsotropic[moveParams.seed] && moveParams.isRotation)
#else
                if (moveParams.isRotation)
#endif
                {
                    // Initialise neighbouring particle.
                    initiateParticle(neighbour, particles[moveParams.seed]);

                    // Recursively recruit neighbours to the cluster.
                    recursiveMoveAssignment(neighbour);
                }
                else
                {
                    // Recursively recruit neighbours to the cluster.
                    recursiveMoveAssignment(moveParams.seed);
                }

                // Check whether the cluster is too large.
                if (nMoving > cutOff) isEarlyExit = true;
            }
        }
    }

    bool VMMC::accept()
    {
        // Abort if early exit condition has been triggered.
        if (isEarlyExit) return false;

        // Any remaining frustrated links must be external to the cluster.
        if (nFrustrated > 0)
        {
            isEarlyExit = true;
            return false;
        }

        // Calculate the approximate Stokes scaling factor.
        double scaleFactor = (nMoving > 1) ? computeHydrodynamicRadius() : 1.0;

        // Stokes drag rejection.
        if (rng() > scaleFactor)
        {
            isEarlyExit = true;
            return false;
        }

        // Energy variables.
        double energy;
        double excessEnergy = 0;

        // Construct pair interaction matrix (finite repulsions only).
        if (isRepusive)
        {
            unsigned int x, y;
            unsigned int nPairs;
            unsigned int pairInteractions[maxInteractions];

            // Check all particles in the moving cluster.
            for (unsigned int i=0;i<nMoving;i++)
            {
                // Get a list of pair interactions.
#ifndef ISOTROPIC
                nPairs = callbacks.interactionsCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);
#else
                nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0], pairInteractions);
#endif

                // Test all pair interactions.
                for (unsigned int j=0;j<nPairs;j++)
                {
#ifndef ISOTROPIC
                    energy = callbacks.pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0], &particles[moveList[i]].preMoveOrientation[0],
                        pairInteractions[j], &particles[pairInteractions[j]].preMovePosition[0],
                        &particles[pairInteractions[j]].preMoveOrientation[0]);
#else
                    energy = callbacks.pairEnergyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                        pairInteractions[j], &particles[pairInteractions[j]].preMovePosition[0]);
#endif

                    x = moveList[i];
                    y = pairInteractions[j];

                    // Make sure leading index is larger.
                    if (x < y)
                    {
                        x = y;
                        y = moveList[i];
                    }

                    // Check to see if pair interaction has already been logged.
                    if (pairEnergyMatrix[x][y] == 0)
                    {
                        interactions[nInteractions][0] = x;
                        interactions[nInteractions][1] = y;
                        nInteractions++;

                        // Store pair energy.
                        pairEnergyMatrix[x][y] = energy;
                    }
                }
            }
        }

        // Check for non-pairwise energy contributions.
        if (callbacks.isNonPairwise)
        {
            // Check all particles in the moving cluster.
            for (unsigned int i=0;i<nMoving;i++)
            {
#ifndef ISOTROPIC
                excessEnergy -= callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif
            }
        }

        // Apply the move.
        swapMoveStatus();

        // Check for overlaps (or finite repulsions).
        for (unsigned int i=0;i<nMoving;i++)
        {
            // Check for non-pairwise energy contributions.
            if (callbacks.isNonPairwise)
            {
#ifndef ISOTROPIC
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif

                // Early exit for large non-pairwise energies.
                if (excessEnergy > 1e6) return false;
            }

            if (!isRepusive)
            {
#ifndef ISOTROPIC
                energy = callbacks.energyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                energy = callbacks.energyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif

                // Overlap.
                if (energy > 1e6) return false;
            }
            else
            {
                double x, y;
                double pairEnergy;
                unsigned int pairInteractions[maxInteractions];

#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0], pairInteractions);
#endif

                for (unsigned int j=0;j<nPairs;j++)
                {
#ifndef ISOTROPIC
                    energy = callbacks.pairEnergyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                        &particles[moveList[i]].preMoveOrientation[0], pairInteractions[j],
                        &particles[pairInteractions[j]].preMovePosition[0], &particles[pairInteractions[j]].preMoveOrientation[0]);
#else
                    energy = callbacks.pairEnergyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                        pairInteractions[j], &particles[pairInteractions[j]].preMovePosition[0]);
#endif

                    // Early exit test for hard core overlaps and large finite energy repulsions.
                    if (energy > 1e6) return false;

                    x = moveList[i];
                    y = pairInteractions[j];

                    if (x < y)
                    {
                        x = y;
                        y = moveList[i];
                    }

                    // Repulsive interaction.
                    if (energy > 0)
                    {
                        // Check that particles didn't previously interact.
                        if (pairEnergyMatrix[x][y] == 0)
                            excessEnergy += energy;
                    }
                    else
                    {
                        // Neighbour isn't part of the moving cluster.
                        if (!particles[pairInteractions[j]].isMoving)
                        {
                            // Particles no longer interact.
                            if (energy == 0)
                            {
                                pairEnergy = pairEnergyMatrix[x][y];

                                // Particles previously felt a repulsive interaction.
                                if (pairEnergy > 0)
                                    excessEnergy -= pairEnergy;
                            }
                        }
                    }
                }
            }
        }

        if (isRepusive || callbacks.isNonPairwise)
        {
            if (rng() > exp(-excessEnergy)) return false;
        }

        // Move successful.
        return true;
    }

    double VMMC::computeHydrodynamicRadius() const
    {
        std::vector<double> centerOfMass(dimension);
        std::vector<double> delta(dimension);

        double hydroRadius;

        // Calculate center of mass of the moving cluster (translations only).
        if (!moveParams.isRotation)
        {
            for (unsigned int i=0;i<nMoving;i++)
            {
                for (unsigned int j=0;j<dimension;j++)
                    centerOfMass[j] += particles[moveList[i]].clusterPosition[j];
            }
        }

        // Second pass to calculate the mean square extent perpendicular to motion.
        for (unsigned int i=0;i<nMoving;i++)
        {
            if (!moveParams.isRotation)
            {
                for (unsigned int j=0;j<dimension;j++)
                    delta[j] = particles[moveList[i]].clusterPosition[j] - centerOfMass[j] / (double) nMoving;
            }
            else
            {
                for (unsigned int j=0;j<dimension;j++)
                    delta[j] = particles[moveList[i]].clusterPosition[j] - particles[moveParams.seed].preMovePosition[j];
            }

            double a1 = delta[0]*moveParams.trialVector[1] - delta[1]*moveParams.trialVector[0];
            hydroRadius = a1*a1;

            if (is3D)
            {
                double a2 = delta[1]*moveParams.trialVector[2] - delta[2]*moveParams.trialVector[1];
                double a3 = delta[2]*moveParams.trialVector[0] - delta[0]*moveParams.trialVector[2];

                hydroRadius += a2*a2 + a3*a3;
            }
        }

        // Calculate scale factor from Stokes' law.
        double rEff = referenceRadius + sqrt(hydroRadius / (double) nMoving);
        double scaleFactor = referenceRadius / rEff;

        // For rotations.
        if (moveParams.isRotation) scaleFactor *= scaleFactor*scaleFactor;

        return scaleFactor;
    }

    void VMMC::computePostMoveParticle(unsigned int particle, int direction, Particle& postMoveParticle)
    {
        // Initialise post-move position and orientation.
        postMoveParticle.postMovePosition = particles[particle].preMovePosition;
#ifndef ISOTROPIC
        postMoveParticle.postMoveOrientation = particles[particle].preMoveOrientation;
#endif

        if (!moveParams.isRotation) // Translation.
        {
            for (unsigned int i=0;i<dimension;i++)
                postMoveParticle.postMovePosition[i] += direction*moveParams.stepSize*moveParams.trialVector[i];
        }
        else                        // Rotation.
        {
            std::vector<double> v1(dimension);
            std::vector<double> v2(dimension);

            // Calculate coordinates relative to the global rotation point.
            for (unsigned int i=0;i<dimension;i++)
                v1[i] = particles[particle].clusterPosition[i] - particles[moveParams.seed].clusterPosition[i];

            // Calculate position rotation vector.
            if (is3D) rotate3D(v1, moveParams.trialVector, v2, direction*moveParams.stepSize);
            else rotate2D(v1, v2, direction*moveParams.stepSize);

            // Update position.
            for (unsigned int i=0;i<dimension;i++)
                postMoveParticle.postMovePosition[i] += v2[i];

#ifndef ISOTROPIC
            // Only update orientations for anisotropic particles.
            if (!isIsotropic[particle])
            {
                // Calculate orientation rotation vector.
                if (is3D) rotate3D(postMoveParticle.postMoveOrientation, moveParams.trialVector, v2, direction*moveParams.stepSize);
                else rotate2D(postMoveParticle.postMoveOrientation, v2, direction*moveParams.stepSize);

                // Update orientation.
                for (unsigned int i=0;i<dimension;i++)
                    postMoveParticle.postMoveOrientation[i] += v2[i];
            }
#endif
        }

        // Only check forward move.
        if (direction == 1)
        {
            // Check custom boundary condition.
            if (callbacks.isCustomBoundary)
            {
#ifndef ISOTROPIC
                bool isOutsideBoundary = callbacks.boundaryCallback(particle,
                    &postMoveParticle.postMovePosition[0], &postMoveParticle.postMoveOrientation[0]);
#else
                bool isOutsideBoundary = callbacks.boundaryCallback(particle, &postMoveParticle.postMovePosition[0]);
#endif
                // Particle has moved outside boundary. Abort move!
                if (isOutsideBoundary) isEarlyExit = true;
            }
        }

        // Apply periodic boundary conditions.
        applyPeriodicBoundaryConditions(postMoveParticle.postMovePosition);
    }

    void VMMC::initiateParticle(unsigned int particle, Particle& linker)
    {
        std::vector<double> delta(dimension);

        // Calculate minumum image separation.
        computeSeparation(linker.clusterPosition, particles[particle].preMovePosition, delta);

        // Assign cluster position based on minumum image separation.
        for (unsigned int i=0;i<dimension;i++)
            particles[particle].clusterPosition[i] = linker.clusterPosition[i] + delta[i];

        // Update move list.
        particles[particle].isMoving = true;
        moveList[nMoving] = particle;
        nMoving++;

        // See if particle was previously participating in a frustrated link.
        if (particles[particle].isFrustrated)
        {
            // Decrement number of frustated links.
            nFrustrated--;
            particles[particle].isFrustrated = false;
            frustratedLinks[particles[particle].posFrustated] = frustratedLinks[nFrustrated];
            particles[frustratedLinks[nFrustrated]].posFrustated = particles[particle].posFrustated;
        }

        // Calculate updated position and orientation.
        computePostMoveParticle(particle, 1, particles[particle]);
    }

    void VMMC::recursiveMoveAssignment(unsigned int particle)
    {
        // Abort if any early exit conditions have been triggered.
        if (!isEarlyExit)
        {
            // Abort if the cluster size cut-off is exceeded.
            if (nMoving <= cutOff)
            {
                Particle reverseMoveParticle(dimension);

                // Calculate coordinates under reverse trial move.
                computePostMoveParticle(particle, -1, reverseMoveParticle);

                unsigned int pairInteractions[maxInteractions];

                // Get list of interactions.
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(particle, &particles[particle].preMovePosition[0],
                    &particles[particle].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(particle,
                    &particles[particle].preMovePosition[0], pairInteractions);
#endif

                // Loop over all interactions.
                for (unsigned int i=0;i<nPairs;i++)
                {
                    unsigned int neighbour = pairInteractions[i];

                    // Make sure link hasn't been tested already.
                    if (!particles[neighbour].isMoving)
                    {
                        // Pre-move pair energy.
#ifndef ISOTROPIC
                        double initialEnergy = callbacks.pairEnergyCallback(particle,
                            &particles[particle].preMovePosition[0], &particles[particle].preMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double initialEnergy = callbacks.pairEnergyCallback(particle, &particles[particle].preMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Post-move pair energy.
#ifndef ISOTROPIC
                        double finalEnergy = callbacks.pairEnergyCallback(particle,
                            &particles[particle].postMovePosition[0], &particles[particle].postMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double finalEnergy = callbacks.pairEnergyCallback(particle, &particles[particle].postMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Pair energy following the reverse virtual move.
#ifndef ISOTROPIC
                        double reverseMoveEnergy = callbacks.pairEnergyCallback(particle,
                            &reverseMoveParticle.postMovePosition[0], &reverseMoveParticle.postMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double reverseMoveEnergy = callbacks.pairEnergyCallback(particle, &reverseMoveParticle.postMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Forward link weight.
                        double linkWeight = std::max(1.0-exp(initialEnergy-finalEnergy),0.0);

                        // Reverse link weight.
                        double reverseLinkWeight = std::max(1.0-exp(initialEnergy-reverseMoveEnergy),0.0);

                        // Test links.
                        if (rng() <= linkWeight)
                        {
                            if (rng() > reverseLinkWeight/linkWeight)
                            {
                                // Particle isn't already participating in a frustrated link.
                                if (!particles[neighbour].isFrustrated)
                                {
                                    particles[neighbour].isFrustrated = true;
                                    particles[neighbour].posFrustated = nFrustrated;
                                    frustratedLinks[nFrustrated] = neighbour;
                                    nFrustrated++;
                                }
                            }
                            else
                            {
                                // Prepare neighbour for virtual move.
                                initiateParticle(neighbour, particles[particle]);

                                // Continue search from neighbour.
                                recursiveMoveAssignment(neighbour);
                            }
                        }
                    }
                }
            }
        }
    }

    void VMMC::swapMoveStatus()
    {
        // Swap the pre- and post-move positions and orientations.
        for (unsigned int i=0;i<nMoving;i++)
        {
            particles[moveList[i]].preMovePosition.swap(particles[moveList[i]].postMovePosition);
#ifndef ISOTROPIC
            particles[moveList[i]].preMoveOrientation.swap(particles[moveList[i]].postMoveOrientation);
#endif
        }

        // Apply any post-move updates.
        for (unsigned int i=0;i<nMoving;i++)
#ifndef ISOTROPIC
            callbacks.postMoveCallback(moveList[i], &particles[moveList[i]].preMovePosition[0], &particles[moveList[i]].preMoveOrientation[0]);
#else
            callbacks.postMoveCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif
    }

    void VMMC::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
    {
        double c = cos(angle);
        double s = sin(angle);

        double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

        v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
        v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
        v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
    }

    void VMMC::rotate2D(std::vector<double>& v1, std::vector<double>& v2, double angle)
    {
        double c = cos(angle);
        double s = sin(angle);

        v2[0] = (v1[0]*c - v1[1]*s) - v1[0];
        v2[1] = (v1[0]*s + v1[1]*c) - v1[1];
    }

    void VMMC::computeSeparation(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& sep)
    {
        for (unsigned int i=0;i<dimension;i++)
        {
            sep[i] = v2[i] - v1[i];

            if (sep[i] < -0.5*boxSize[i])
            {
                sep[i] += boxSize[i];
            }
            else
            {
                if (sep[i] >= 0.5*boxSize[i])
                {
                    sep[i] -= boxSize[i];
                }
            }
        }
    }

    void VMMC::applyPeriodicBoundaryConditions(std::vector<double>& vec)
    {
        for (unsigned int i=0;i<vec.size();i++)
        {
            if (vec[i] < 0)
            {
                vec[i] += boxSize[i];
            }
            else
            {
                if (vec[i] >= boxSize[i])
                {
                    vec[i] -= boxSize[i];
                }
            }
        }
    }

    double VMMC::computeNorm(std::vector<double>& vec)
    {
        double normSquared = 0;

        for (unsigned int i=0;i<vec.size();i++)
            normSquared += vec[i]*vec[i];

        return sqrt(normSquared);
    }
}
