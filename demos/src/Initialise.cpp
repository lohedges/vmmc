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

#include "Initialise.h"

Initialise::Initialise()
{
}

void Initialise::random(std::vector <Particle>& particles, CellList& cells, Box& box, MTRand& rng)
{
    for (unsigned i=0;i<particles.size();i++)
    {
        // Current number of attempted particle insertions.
        unsigned int nTrials = 0;

        // Whether particle overlaps.
        bool isOverlap = true;

        // Whether this is the first attempted insertion.
        bool isFirstAttempt = true;

        // Temporary vector.
        std::vector <double> vec(box.dimension);

        // Set particle index.
        particles[i].index = i;

        // Keep trying to insert particle until there is no overlap.
        while (isOverlap)
        {
            nTrials++;

            // Generate a random position.
            vec[0] = rng()*box.boxSize[0];
            vec[1] = rng()*box.boxSize[1];
            if (box.dimension == 3) vec[2] = rng()*box.boxSize[2];

            particles[i].position = vec;

            // Generate a random orientation.
            vec[0] = rng.randNorm(0,1);
            vec[1] = rng.randNorm(0,1);
            if (box.dimension == 3) vec[2] = rng.randNorm(0,1);

            // Calculate vector norm.
            double norm = vec[0]*vec[0] + vec[1]*vec[1];
            if (box.dimension == 3) norm += vec[2]*vec[2];
            norm = sqrt(norm);

            // Convert orientation to a unit vector.
            vec[0] /= norm;
            vec[1] /= norm;
            if (box.dimension == 3) vec[2] /= norm;

            particles[i].orientation = vec;

            // Calculate the particle's cell index.
            unsigned int newCell = cells.getCell(particles[i]);

            // Update the cell list.
            if (isFirstAttempt)
            {
                particles[i].cell = newCell;
                cells.initCell(newCell, particles[i]);
                isFirstAttempt = false;
            }
            else
            {
                unsigned int oldCell = particles[i].cell;

                if (oldCell != newCell)
                    cells.updateCell(newCell, particles[i], particles);
            }

            // See if there is any overlap between particles.
            isOverlap = checkOverlap(particles[i], particles, cells, box);

            // Check trial limit isn't exceeded.
            if (nTrials == MAX_TRIALS)
            {
                std::cerr << "[ERROR]: Maximum number of trial insertions reached.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
}

bool Initialise::checkOverlap(Particle& particle, std::vector <Particle>& particles, CellList& cells, Box& box)
{
    unsigned int cell, neighbour;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells[particle.cell].neighbours.size();i++)
    {
        cell = cells[particle.cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure particles are different.
            if (neighbour != particle.index)
            {
                // Particle separtion vector.
                std::vector <double> sep(box.dimension);

                sep[0] = particle.position[0] - particles[neighbour].position[0];
                sep[1] = particle.position[1] - particles[neighbour].position[1];
                if (box.dimension == 3) sep[2] = particle.position[2] - particles[neighbour].position[2];

                // Compute minimum image.
                box.minimumImage(sep);

                // Squared norm of vector.
                double normSqd = sep[0]*sep[0] + sep[1]*sep[1];
                if (box.dimension == 3) normSqd += sep[2]*sep[2];

                // Overlap if normSqd is less than particle diameter (box is scaled in diameter units).
                if (normSqd < 1) return true;
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
