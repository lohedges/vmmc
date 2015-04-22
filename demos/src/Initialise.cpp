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

void Initialise::random(std::vector <Particle>& particles, CellList& cells, Box& box, MersenneTwister& rng)
{
    for (unsigned i=0;i<particles.size();i++)
    {
        // Current number of attempted particle insertions.
        unsigned int nTrials = 0;

        // Whether particle overlaps.
        bool isOverlap = true;

        // Temporary vector.
        std::vector <double> vec(box.dimension);

        // Set particle index.
        particles[i].index = i;

        // Keep trying to insert particle until there is no overlap.
        while (isOverlap)
        {
            nTrials++;

            // Generate a random position.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng()*box.boxSize[j];

            particles[i].position = vec;

            // Generate a random orientation.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng.normal();

            // Calculate vector norm.
            double norm = 0;
            for (unsigned int j=0;j<box.dimension;j++)
                norm += vec[j]*vec[j];
            norm = sqrt(norm);

            // Convert orientation to a unit vector.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] /= norm;

            particles[i].orientation = vec;

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // See if there is any overlap between particles.
            isOverlap = checkOverlap(particles[i], particles, cells, box);

            // Check trial limit isn't exceeded.
            if (nTrials == MAX_TRIALS)
            {
                std::cerr << "[ERROR]: Maximum number of trial insertions reached.\n";
                exit(EXIT_FAILURE);
            }
        }

        // Update cell list.
        cells.initCell(particles[i].cell, particles[i]);
    }
}

bool Initialise::checkOverlap(Particle& particle, std::vector <Particle>& particles, CellList& cells, Box& box)
{
    unsigned int cell, neighbour;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
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

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = particle.position[k] - particles[neighbour].position[k];

                // Compute minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Overlap if normSqd is less than particle diameter (box is scaled in diameter units).
                if (normSqd < 1) return true;
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
