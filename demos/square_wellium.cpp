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

#include <limits>

#include "Demo.h"
#include "VMMC.h"

// CALLBACK PROTOTYPES

double computeEnergy(unsigned int, double[], double[]);
double computePairEnergy(unsigned int, double[], double[], unsigned int, double[], double[]);
unsigned int computeInteractions(unsigned int, double[], double[], unsigned int[]);
void applyPostMoveUpdates(unsigned int, double[], double[]);

// FUNCTION PROTOTYPES

double getEnergy();
void appendXyzTrajectory(std::vector <Particle>&, bool);
void minimumImage(std::vector <double>&);

// GLOBALS

std::vector <Particle> particles;       // particle container
CellList cells;                         // cell list
unsigned int dimension = 3;             // dimension of simulation box
unsigned int nParticles = 1000;         // number of particles
double interactionEnergy = 2.6;         // pair interaction energy scale (in units of kBT)
double interactionRange = 0.1;          // size of interaction range (in units of particle diameter)
double density = 0.05;                  // particle density
double squaredCutOffDistance;           // squared interaction cut-off distance
double baseLength;                      // base length of simulation box
unsigned int maxInteractions = 15;      // maximum number of interactions per particle

// PORTABLE NUMERIC CONSTANTS

double INF = std::numeric_limits<double>::infinity();

// MAIN FUNCTION

int main(int argc, char** argv)
{
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

    // Work out cut-off distance.
    double cutOffDistance = 1 + interactionRange;

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, cutOffDistance);

    // Squared cut-off distance.
    squaredCutOffDistance = cutOffDistance*cutOffDistance;

    // Initialise random number generator.
    MTRand rng;

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

    // Initialise VMMC callback functions.
    VMMC_energyCallback energyCallback = computeEnergy;
    VMMC_pairEnergyCallback pairEnergyCallback = computePairEnergy;
    VMMC_interactionsCallback interactionsCallback = computeInteractions;
    VMMC_postMoveCallback postMoveCallback = applyPostMoveUpdates;

    // Initalise VMMC object.
    VMMC vmmc(nParticles, dimension, coordinates, orientations, 0.15, 0.2, 0.5, 0.5, maxInteractions,
            &boxSize[0], false, energyCallback, pairEnergyCallback, interactionsCallback, postMoveCallback);

    // Execute the simulation.
    for (unsigned int i=0;i<1000;i++)
    {
        // Increment simulation by 1000 Monte Carlo Sweeps.
        vmmc += 1000*nParticles;

        // Append particle coordinates to a xyz trajectory.
        if (i == 0) appendXyzTrajectory(particles, true);
        else appendXyzTrajectory(particles, false);

        // Report.
        printf("sweeps = %9.4e, energy = %5.4f\n", ((double) (i+1)*1000), getEnergy());
    }

    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}

// CALLBACK DEFINITIONS

double computeEnergy(unsigned int particle, double position[], double orientation[])
{
    unsigned int cell;      // cell index
    unsigned int neighbour; // index of neighbouring particle
    double energy = 0;      // energy counter

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector <double> sep(dimension);

                // Calculate separation.
                for (unsigned int k=0;k<dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<dimension;k++)
                    normSqd += sep[k]*sep[k];

                if (normSqd < 1) return INF;
                if (normSqd < squaredCutOffDistance) energy -= interactionEnergy;
            }
        }
    }

    return energy;
}

double computePairEnergy(unsigned int particle1, double position1[],
        double orientation1[], unsigned int particle2, double position2[], double orientation2[])
{
    // Separation vector.
    std::vector <double> sep(dimension);

    // Calculate separation.
    for (unsigned int i=0;i<dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    minimumImage(sep);

    double normSqd = 0;

    // Calculate squared norm of vector.
    for (unsigned int i=0;i<dimension;i++)
        normSqd += sep[i]*sep[i];

    if (normSqd < 1) return INF;
    if (normSqd < squaredCutOffDistance) return -interactionEnergy;
    return 0;
}

unsigned int computeInteractions(unsigned int particle,
        double position[], double orientation[], unsigned int interactions[])
{
    unsigned int cell;              // cell index
    unsigned int neighbour;         // index of neighbouring particle
    unsigned int nInteractions = 0; // interaction counter

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector <double> sep(dimension);

                // Compute separation.
                for (unsigned int k=0;k<dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Particles interact.
                if (normSqd < squaredCutOffDistance)
                {
                    interactions[nInteractions] = neighbour;
                    nInteractions++;

                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR]: Maximum number of interactions exceeded!\n";
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    return nInteractions;
}

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

double getEnergy()
{
    double energy = 0;

    for (unsigned int i=0;i<nParticles;i++)
        energy += computeEnergy(i, &particles[i].position[0], &particles[i].orientation[0]);

    return energy/(2*nParticles);
}

void appendXyzTrajectory(std::vector <Particle>& particles, bool clearFile)
{
    FILE* pFile;

    // Wipe existing trajectory file.
    if (clearFile)
    {
        pFile = fopen("trajectory.xyz", "w");
        fclose(pFile);
    }

    pFile = fopen("trajectory.xyz", "a");
    fprintf(pFile, "%lu\n\n", particles.size());

    for (unsigned int i=0;i<particles.size();i++)
    {
        fprintf(pFile, "0 %5.4f %5.4f %5.4f\n",
            particles[i].position[0], particles[i].position[1], (dimension == 3) ? particles[i].position[2] : 0);
    }

    fclose(pFile);
}

void minimumImage(std::vector <double>& separation)
{
    for (unsigned int i=0;i<dimension;i++)
    {
        if (separation[i] < -0.5*baseLength)
        {
            separation[i] += baseLength;
        }
        else
        {
            if (separation[i] >= 0.5*baseLength)
            {
                separation[i] -= baseLength;
            }
        }
    }
}
