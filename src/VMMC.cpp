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

#include "VMMC.h"

VMMC_Particle::VMMC_Particle() {}

VMMC_Particle::VMMC_Particle(unsigned int dimension)
{
    // Resize position/orientation vectors.
    preMovePosition.resize(dimension);
    preMoveOrientation.resize(dimension);
    postMovePosition.resize(dimension);
    postMoveOrientation.resize(dimension);
    pseudoPosition.resize(dimension);
}

VMMC::VMMC(unsigned int nParticles_,
           unsigned int dimension_,
           double coordinates[],
           double orientations[],
           double maxTrialTranslation_,
           double maxTrialRotation_,
           double probTranslate_,
           double referenceRadius_,
           unsigned int maxInteractions_,
           double boxSize_[],
           bool isRepusive_,
           const VMMC_energyCallback& energyCallback_,
           const VMMC_pairEnergyCallback& pairEnergyCallback_,
           const VMMC_interactionsCallback& interactionsCallback_,
           const VMMC_postMoveCallback& postMoveCallback_) :

           nParticles(nParticles_),
           dimension(dimension_),
           maxTrialTranslation(maxTrialTranslation_),
           maxTrialRotation(maxTrialRotation_),
           probTranslate(probTranslate_),
           referenceRadius(referenceRadius_),
           maxInteractions(maxInteractions_),
           isRepusive(isRepusive_),
           energyCallback(energyCallback_),
           pairEnergyCallback(pairEnergyCallback_),
           interactionsCallback(interactionsCallback_),
           postMoveCallback(postMoveCallback_)
{
    // Check dimensionality is valid.
    if (dimension == 3) is3D = true;
    else if (dimension == 2) is3D = false;
    else
    {
        std::cerr << "[ERROR]: Invalid dimensionality!\n";
        exit(EXIT_FAILURE);
    }

    // Store simulation box size.
    boxSize.resize(dimension);
    for (unsigned int i=0;i<dimension;i++)
        boxSize[i] = boxSize_[i];

    // Allocate memory.
    moveParams.trialVector.resize(dimension);
    particles.resize(nParticles);
    moveList.resize(nParticles);
    clusterTranslations.resize(nParticles);
    clusterRotations.resize(nParticles);
    frustratedLinks.resize(nParticles);

    // Create particle container.
    for (unsigned int i=0;i<nParticles;i++)
    {
        // Resize vectors.
        particles[i].preMovePosition.resize(dimension);
        particles[i].preMoveOrientation.resize(dimension);
        particles[i].postMovePosition.resize(dimension);
        particles[i].postMoveOrientation.resize(dimension);
        particles[i].pseudoPosition.resize(dimension);

        // Copy particle coordinates and orientations.
        for (unsigned int j=0;j<dimension;j++)
        {
            particles[i].preMovePosition[j] = coordinates[dimension*i + j];
            particles[i].preMoveOrientation[j] = orientations[dimension*i + j];
        }
    }

    // Allocate memory for pair interaction matrix (finite repulsions only).
    if (isRepusive)
    {
        interactions.resize(nParticles);
        for (unsigned int i=0;i<nParticles;i++)
            interactions[i].resize(2);

        // Construct a triangular matrix to save memory.
        pairEnergyMatrix.resize(nParticles);
        for (unsigned int i=0;i<nParticles;i++)
            pairEnergyMatrix[i].resize(i);
    }

    // Print version info.
#ifdef COMMIT
    std::cout << "Initialised VMMC: commit " <<  COMMIT << '\n';
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

    // Propose a move for the new pseudocluster.
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
            if (moveParams.isRotation) clusterRotations[nMoving]++;
            else clusterTranslations[nMoving]++;
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

const std::vector <unsigned long long>& VMMC::getClusterTranslations() const
{
    return clusterTranslations;
}

void VMMC::getClusterRotations(unsigned long long clusterStatistics[]) const
{
    for (unsigned int i=0;i<nParticles;i++)
        clusterStatistics[i] = clusterRotations[i];
}

const std::vector <unsigned long long>& VMMC::getClusterRotations() const
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
    moveParams.seed = rng.randInt(nParticles-1);

    // Get a uniform random number in range [0-1].
    double r = rng();

    // Make sure the divisor doesn't blow things up.
    while (r == 0) r = rng();

    // Cluster size cut-off.
    cutOff = int(1.0/r);

    // Choose a random point on the surface of the unit sphere/circle.
    for (unsigned int i=0;i<dimension;i++)
        moveParams.trialVector[i] = rng.randNorm(0,1);

    // Normalise the trial vector.
    double norm = computeNorm(moveParams.trialVector);
    for (unsigned int i=0;i<dimension;i++)
        moveParams.trialVector[i] /= norm;

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
    }

    // Initialise the seed particle.
    particles[moveParams.seed].pseudoPosition = particles[moveParams.seed].preMovePosition;
    initiateParticle(moveParams.seed, particles[moveParams.seed]);

    // Recursively recruit neighbours to the cluster.
    recursiveMoveAssignment(moveParams.seed);

    // Check whether the cluster is too large.
    if (nMoving > cutOff) isEarlyExit = true;
}

bool VMMC::accept()
{
    // Any remaining frustrated links must be external to the cluster.
    if (nFrustrated > 0)
    {
        isEarlyExit = true;
        return false;
    }

    // Calculate the approximate Stoke's scaling factor.
    double scaleFactor = (nMoving > 1) ? computeHydrodynamicRadius() : 1.0;

    // Stoke's drag rejection.
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

        // Check all particles in the pseudo-cluster.
        for (unsigned int i=0;i<nMoving;i++)
        {
            // Get a list of pair interactions.
            nPairs = interactionsCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);

            // Test all pair interactions.
            for (unsigned int j=0;j<nPairs;j++)
            {
                energy = pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0], &particles[moveList[i]].preMoveOrientation[0],
                        pairInteractions[j], &particles[pairInteractions[j]].preMovePosition[0],
                        &particles[pairInteractions[j]].preMoveOrientation[0]);

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
                }

                // Store pair energy.
                pairEnergyMatrix[x][y] = energy;
            }
        }
    }

    // Apply the move.
    swapMoveStatus();

    // Check for overlaps (or finite repulsions).
    for (unsigned int i=0;i<nMoving;i++)
    {
        if (!isRepusive)
        {
            energy = energyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);

            // Overlap.
            if (energy > 1e6) return false;
        }
        else
        {
            double x, y;
            double pairEnergy;
            unsigned int pairInteractions[maxInteractions];

            unsigned int nPairs = interactionsCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);

            for (unsigned int j=0;j<nPairs;j++)
            {
                energy = pairEnergyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                        &particles[moveList[i]].preMoveOrientation[0], pairInteractions[j],
                        &particles[pairInteractions[j]].preMovePosition[0], &particles[pairInteractions[j]].preMoveOrientation[0]);

                x = moveList[i];
                y = pairInteractions[j];

                if (x < y)
                {
                    x = y;
                    y = moveList[i];
                }

                if (energy > 0)             // Repulsive interaction.
                {
                    // Check that particles didn't previously interact.
                    if (pairEnergyMatrix[x][y] == 0)
                        excessEnergy += energy;
                }
                else
                {
                    if (energy == 0)       // Particles no longer interact.
                    {
                        pairEnergy = pairEnergyMatrix[x][y];

                        // Particles previously felt a repulsive interaction.
                        if (pairEnergy > 0) excessEnergy -= energy;
                    }
                }
            }
        }
    }

    if (isRepusive)
    {
        if (rng() > exp(-excessEnergy)) return false;
    }

    // Move succesful.
    return true;
}

double VMMC::computeHydrodynamicRadius() const
{
    std::vector <double> centerOfMass(dimension);
    std::vector <double> delta(dimension);

    double hydroRadius;

    // Calculate center of mass of the moving cluster (translations only).
    if (!moveParams.isRotation)
    {
        for (unsigned int i=0;i<nMoving;i++)
        {
            for (unsigned int j=0;j<dimension;j++)
                centerOfMass[j] += particles[moveList[i]].pseudoPosition[j];
        }
    }

    // Second pass to calculate the mean square extent perpendicular to motion.
    for (unsigned int i=0;i<nMoving;i++)
    {
        if (!moveParams.isRotation)
        {
            for (unsigned int j=0;j<dimension;j++)
                delta[j] = particles[moveList[i]].pseudoPosition[j] - centerOfMass[j] / (double) nMoving;
        }
        else
        {
            for (unsigned int j=0;j<dimension;j++)
                delta[j] = particles[moveList[i]].pseudoPosition[j] - particles[moveParams.seed].preMovePosition[j];
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

void VMMC::computeCoords(unsigned int particle, VMMC_Particle& postMoveParticle)
{
    // Initialise post-move position and orientation.
    postMoveParticle.postMovePosition = particles[particle].preMovePosition;
    postMoveParticle.postMoveOrientation = particles[particle].preMoveOrientation;

    if (!moveParams.isRotation) // Translation.
    {
        for (unsigned int i=0;i<dimension;i++)
            postMoveParticle.postMovePosition[i] += moveParams.stepSize*moveParams.trialVector[i];
    }
    else                        // Rotation.
    {
        std::vector <double> v1(dimension);
        std::vector <double> v2(dimension);

        // Calculate coordinates relative to the global rotation point.
        for (unsigned int i=0;i<dimension;i++)
            v1[i] = postMoveParticle.pseudoPosition[i] - particles[moveParams.seed].pseudoPosition[i];

        // Calculate position rotation vector.
        if (is3D) rotate3D(v1, moveParams.trialVector, v2, moveParams.stepSize);
        else rotate2D(v1, v2, moveParams.stepSize);

        // Update position.
        for (unsigned int i=0;i<dimension;i++)
            postMoveParticle.postMovePosition[i] += v2[i];

        // Calculate orientation rotation vector.
        if (is3D) rotate3D(postMoveParticle.preMoveOrientation, moveParams.trialVector, v2, moveParams.stepSize);
        else rotate2D(postMoveParticle.preMoveOrientation, v2, moveParams.stepSize);

        // Update orientation.
        for (unsigned int i=0;i<dimension;i++)
            postMoveParticle.postMoveOrientation[i] += v2[i];
    }

    // Apply periodic boundary conditions.
    applyPeriodicBoundaryConditions(postMoveParticle.postMovePosition);
}

void VMMC::initiateParticle(unsigned int particle, VMMC_Particle& linker)
{
    std::vector <double> delta(dimension);

    // Calculate minumum image separation.
    computeSeparation(linker.pseudoPosition, particles[particle].preMovePosition, delta);

    // Assign pseudo-coordinate based on minumum image separation.
    for (unsigned int i=0;i<dimension;i++)
        particles[particle].pseudoPosition[i] = linker.pseudoPosition[i] + delta[i];

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

    // Calculate final coordinates.
    computeCoords(particle, particles[particle]);
}

void VMMC::recursiveMoveAssignment(unsigned int particle)
{
    // Abort if the cluster size cut-off is exceeded.
    if (nMoving <= cutOff)
    {
        VMMC_Particle reverseMoveParticle(dimension);

        // Calculate coordinates under reverse trial move.
        moveParams.stepSize = -moveParams.stepSize;
        computeCoords(particle, reverseMoveParticle);
        moveParams.stepSize = -moveParams.stepSize;

        unsigned int pairInteractions[maxInteractions];

        // Get list of interactions.
        unsigned int nPairs = interactionsCallback(particle, &particles[particle].preMovePosition[0],
                &particles[particle].preMoveOrientation[0], pairInteractions);

        // Loop over all interactions.
        for (unsigned int i=0;i<nPairs;i++)
        {
            unsigned int neighbour = pairInteractions[i];

            // Make sure link hasn't been tested already.
            if (!particles[neighbour].isMoving)
            {
                // Pre-move pair energy.
                double initialEnergy = pairEnergyCallback(particle,
                        &particles[particle].preMovePosition[0], &particles[particle].preMoveOrientation[0],
                        neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);

                // Post-move pair energy.
                double finalEnergy = pairEnergyCallback(particle,
                        &particles[particle].postMovePosition[0], &particles[particle].postMovePosition[0],
                        neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);

                // Pair energy following the reverse virtual move.
                double reverseMoveEnergy = pairEnergyCallback(particle,
                        &reverseMoveParticle.postMovePosition[0], &reverseMoveParticle.postMoveOrientation[0],
                        neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);

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

                        // Continue search from neighbor.
                        recursiveMoveAssignment(neighbour);
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
        particles[moveList[i]].preMoveOrientation.swap(particles[moveList[i]].postMoveOrientation);
    }

    // Apply any post-move updates.
    for (unsigned int i=0;i<nMoving;i++)
        postMoveCallback(moveList[i], &particles[moveList[i]].preMovePosition[0], &particles[moveList[i]].preMoveOrientation[0]);
}

void VMMC::rotate3D(std::vector <double>& v1, std::vector <double>& v2, std::vector <double>& v3, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] *v1[2]*v2[2];

    v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
    v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
    v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
}

void VMMC::rotate2D(std::vector <double>& v1, std::vector <double>& v2, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    v2[0] = v1[0] - (v1[0]*c - v1[1]*s);
    v2[1] = v1[1] - (v1[0]*s + v1[1]*c);
}

void VMMC::computeSeparation(std::vector <double>& v1, std::vector <double>& v2, std::vector <double>& sep)
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

void VMMC::applyPeriodicBoundaryConditions(std::vector <double>& vec)
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

double VMMC::computeNorm(std::vector <double>& vec)
{
    double normSquared = 0;

    for (unsigned int i=0;i<vec.size();i++)
        normSquared += vec[i]*vec[i];

    return sqrt(normSquared);
}
