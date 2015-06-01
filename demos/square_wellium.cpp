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

// CALLBACK PROTOTYPES

void applyPostMoveUpdates(unsigned int, double[], double[]);

// FUNCTION PROTOTYPES

double getEnergy(SquareWellium&);

// GLOBALS

unsigned int dimension = 3;             // dimension of simulation box
unsigned int nParticles = 1000;         // number of particles
std::vector <Particle> particles;       // particle container
CellList cells;                         // cell list

// MAIN FUNCTION

int main(int argc, char** argv)
{
    // Simulation parameters...

    double interactionEnergy = 2.6;         // pair interaction energy scale (in units of kBT)
    double interactionRange = 0.1;          // size of interaction range (in units of particle diameter)
    double density = 0.05;                  // particle density
    double baseLength;                      // base length of simulation box
    unsigned int maxInteractions = 15;      // maximum number of interactions per particle

    // Resize particle container.
    particles.resize(nParticles);

    // Work out base length of simulation box (particle diameter is one).
    if (dimension == 2) baseLength = std::pow((nParticles*M_PI)/(2.0*density), 1.0/2.0);
    else baseLength = std::pow((nParticles*M_PI)/(6.0*density), 1.0/3.0);

    std::vector <double> boxSize;
    for (unsigned int i=0;i<dimension;i++)
        boxSize.push_back(baseLength);

    // Initialise simulation box object.
    Box box(boxSize);

    // Initialise input/output class,
    InputOutput io;

    // Create VMD script.
    io.vmdScript(boxSize);

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, 1 + interactionRange);

    // Initialise square well potential object.
    SquareWellium squareWellium(box, particles, cells,
            maxInteractions, interactionEnergy, interactionRange);

    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise particle initialisation object.
    Initialise initialise;

    // Generate a random particle configuration.
    initialise.random(particles, cells, box, rng);

    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    // Copy particle coordinates/orientations into C-style arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }
    }

    // Initialise the VMMC callback functions...

    // Bind model specific member functions for potential callbacks.
    using namespace std::placeholders;
    VMMC_energyCallback energyCallback =
        std::bind(&SquareWellium::computeEnergy, squareWellium, _1, _2, _3);
    VMMC_pairEnergyCallback pairEnergyCallback =
        std::bind(&SquareWellium::computePairEnergy, squareWellium, _1, _2, _3, _4, _5, _6);
    VMMC_interactionsCallback interactionsCallback =
        std::bind(&SquareWellium::computeInteractions, squareWellium, _1, _2, _3, _4);

    // Wrap a free function for applying the post-move updates.
    VMMC_postMoveCallback postMoveCallback = applyPostMoveUpdates;

    // Initalise VMMC object.
    VMMC vmmc(nParticles, dimension, coordinates, orientations, 0.15, 0.2, 0.5, 0.5, maxInteractions,
            &boxSize[0], true, false, energyCallback, pairEnergyCallback, interactionsCallback, postMoveCallback);

    // Execute the simulation.
    for (unsigned int i=0;i<1000;i++)
    {
        // Increment simulation by 1000 Monte Carlo Sweeps.
        vmmc += 1000*nParticles;

        // Append particle coordinates to an xyz trajectory.
        if (i == 0) io.appendXyzTrajectory(dimension, particles, true);
        else io.appendXyzTrajectory(dimension, particles, false);

        // Report.
        printf("sweeps = %9.4e, energy = %5.4f\n", ((double) (i+1)*1000), getEnergy(squareWellium));
    }

    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}

// CALLBACK DEFINITIONS

void applyPostMoveUpdates(unsigned int particle, double position[], double orientation[])
{
    // Copy coordinates/orientations.
    for (unsigned int i=0;i<dimension;i++)
    {
        particles[particle].position[i] = position[i];
        particles[particle].orientation[i] = orientation[i];
    }

    // Calculate the particle's cell index.
    unsigned int newCell = cells.getCell(particles[particle]);

    // Update cell lists if necessary.
    if (particles[particle].cell != newCell)
        cells.updateCell(newCell, particles[particle], particles);
}

// FUNCTION DEFINITIONS

double getEnergy(SquareWellium& squareWellium)
{
    double energy = 0;

    for (unsigned int i=0;i<nParticles;i++)
        energy += squareWellium.computeEnergy(i, &particles[i].position[0], &particles[i].orientation[0]);

    return energy/(2*nParticles);
}
