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

#include "Demo.h"
#include "VMMC.h"

int main(int argc, char** argv)
{
    // Simulation parameters.
    unsigned int dimension = 3;                     // dimension of simulation box
    unsigned int nParticles = 1000;                 // number of particles
    double interactionEnergy = 2.3;                 // pair interaction energy scale (in units of kBT)
    double interactionRange = 1.1;                  // size of interaction range (in units of particle diameter)
    double wallInteractionEnergy = 3.5;             // particle-wall interaction energy scale (in units of kBT)
    double wallInteractionRange = 0.1;              // particle-wall interaction range (particle surface to wall distance)
    double density = 0.05;                          // particle density
    double baseLength;                              // base length of simulation box
    unsigned int maxInteractions = 15;              // maximum number of interactions per particle

    // Data structures.
    std::vector<Particle> particles(nParticles);    // particle container
    CellList cells;                                 // cell list
#ifndef ISOTROPIC
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
#endif

    // Work out base length of simulation box (particle diameter is one).
    if (dimension == 2) baseLength = std::pow((nParticles*M_PI)/(2.0*density), 1.0/2.0);
    else baseLength = std::pow((nParticles*M_PI)/(6.0*density), 1.0/3.0);

    std::vector<double> boxSize;
    std::vector<bool> isPeriodic;
    for (unsigned int i=0;i<dimension;i++)
    {
        boxSize.push_back(baseLength);

        // Final dimension is non-periodic.
        if (i < (dimension - 1)) isPeriodic.push_back(true);
        else isPeriodic.push_back(false);
    }

    // Initialise simulation box object.
    Box box(boxSize, isPeriodic);

    // Initialise input/output class,
    InputOutput io;

    // Create VMD script.
    io.vmdScript(boxSize);

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise the square well potential model.
    SquareWelliumWall squareWelliumWall(box, particles, cells, maxInteractions,
        interactionEnergy, interactionRange, wallInteractionEnergy, wallInteractionRange);

    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise particle initialisation object.
    Initialise initialise;

    // Generate a random particle configuration.
    initialise.random(particles, cells, box, rng, false);

    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
#ifndef ISOTROPIC
    double orientations[dimension*nParticles];
#endif

    // Copy particle coordinates and orientations into C-style arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = particles[i].position[j];
#ifndef ISOTROPIC
            orientations[dimension*i + j] = particles[i].orientation[j];
#endif
        }

#ifndef ISOTROPIC
        // Set all particles as isotropic.
        isIsotropic[i] = true;
#endif
    }

    // Initialise the VMMC callback functions.
    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
#ifndef ISOTROPIC
    callbacks.energyCallback =
        std::bind(&SquareWelliumWall::computeEnergy, squareWelliumWall, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&SquareWelliumWall::computePairEnergy, squareWelliumWall, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&SquareWelliumWall::computeInteractions, squareWelliumWall, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&SquareWelliumWall::applyPostMoveUpdates, squareWelliumWall, _1, _2, _3);
    callbacks.nonPairwiseCallback =
        std::bind(&SquareWelliumWall::computeWallEnergy, squareWelliumWall, _1, _2, _3);
    callbacks.boundaryCallback =
        std::bind(&SquareWelliumWall::isOutsideBoundary, squareWelliumWall, _1, _2, _3);
#else
    callbacks.energyCallback =
        std::bind(&SquareWelliumWall::computeEnergy, squareWelliumWall, _1, _2);
    callbacks.pairEnergyCallback =
        std::bind(&SquareWelliumWall::computePairEnergy, squareWelliumWall, _1, _2, _3, _4);
    callbacks.interactionsCallback =
        std::bind(&SquareWelliumWall::computeInteractions, squareWelliumWall, _1, _2, _3);
    callbacks.postMoveCallback =
        std::bind(&SquareWelliumWall::applyPostMoveUpdates, squareWelliumWall, _1, _2);
    callbacks.nonPairwiseCallback =
        std::bind(&SquareWelliumWall::computeWallEnergy, squareWelliumWall, _1, _2);
    callbacks.boundaryCallback =
        std::bind(&SquareWelliumWall::isOutsideBoundary, squareWelliumWall, _1, _2);
#endif

    // Initialise VMMC object.
#ifndef ISOTROPIC
    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], isIsotropic, false, callbacks);
#else
    vmmc::VMMC vmmc(nParticles, dimension, coordinates,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], false, callbacks);
#endif

    // Execute the simulation.
    for (unsigned int i=0;i<1000;i++)
    {
        // Increment simulation by 1000 Monte Carlo Sweeps.
        vmmc += 1000*nParticles;

        // Append particle coordinates to an xyz trajectory.
        if (i == 0) io.appendXyzTrajectory(dimension, particles, true);
        else io.appendXyzTrajectory(dimension, particles, false);

        // Report.
        printf("sweeps = %9.4e, energy = %5.4f\n", ((double) (i+1)*1000), squareWelliumWall.getEnergy());
    }

    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
