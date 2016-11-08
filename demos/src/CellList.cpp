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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "CellList.h"
#include "Particle.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

CellList::CellList() : dimension(3)
{
}

CellList::CellList(unsigned int dimension_, const std::vector<double>& boxSize, double range) : dimension(dimension_)
{
    this->initialise(boxSize, range);
}

CellList& CellList::operator = (const CellList& cells)
{
    // resize cell list
    (*this).resize(cells.size());

    // copy across individual cells
    for (unsigned int i=0;i<cells.size();i++)
    {
        at(i) = cells[i];
    }

    return *this;
}

CellList& CellList::operator = (const std::vector<Cell>& cells)
{
    // resize cell list
    (*this).resize(cells.size());

    // copy across individual cells
    for (unsigned int i=0;i<cells.size();i++)
    {
        at(i) = cells[i];
    }

    return *this;
}

void CellList::initialise(const std::vector<double>& boxSize, double range)
{
    unsigned int i,j,k,m;
    unsigned int a,b,c;
    unsigned int x,y,z;
    unsigned int nn,nnCount;

    cellsPerAxis.resize(dimension);
    cellSpacing.resize(dimension);

    for (i=0;i<dimension;i++)
    {
        cellsPerAxis[i] = 1;

        while ((boxSize[i] / (double) cellsPerAxis[i]) > range)
        {
            cellsPerAxis[i]++;
        }
        cellsPerAxis[i]--;
        cellSpacing[i] = boxSize[i] / (double) cellsPerAxis[i];

        // check that number of cells per axis is large enough
        if (cellsPerAxis[i] < 3)
        {
            std::cerr << "[ERROR] CellList: Simulation box is too small (min cells per axis is 3)\n";
            exit(EXIT_FAILURE);
        }
    }

    // Estimate maximum number of particles per cell from interaction range.
    // (Assumes particle diameter is one.)
    if (dimension == 3) maxParticles = (cellSpacing[0]*cellSpacing[1]*cellSpacing[2]) / ((4.0/3.0)*M_PI*0.5*0.5*0.5);
    else maxParticles = (cellSpacing[0]*cellSpacing[1]) / (M_PI*0.5*0.5);

    // Add a buffer, e.g. if particles can overlap.
    maxParticles += 10;

    nCells = cellsPerAxis[0]*cellsPerAxis[1];
    if (dimension == 3) nCells *= cellsPerAxis[2];

    // resize cell list array
    (*this).resize(nCells);

    if (dimension == 3)
    {
        nNeighbours = 27;

        // loop over all cells x direction
        for (i=0;i<cellsPerAxis[0];i++)
        {
            // loop over all cells y direction
            for (j=0;j<cellsPerAxis[1];j++)
            {
                // loop over all cells z direction
                for (k=0;k<cellsPerAxis[2];k++)
                {
                    // cell index
                    m = i + cellsPerAxis[0]*j + cellsPerAxis[0]*cellsPerAxis[1]*k;

                    // resize neighbour list and particle arrays
                    at(m).neighbours.resize(nNeighbours);
                    at(m).particles.resize(maxParticles);

                    nnCount = 0;

                    // x loop for nearest neighbours
                    for (a=0;a<3;a++)
                    {
                        x = (i+(a-1)+cellsPerAxis[0])%cellsPerAxis[0];

                        // y loop for nearest neighbours
                        for (b=0;b<3;b++)
                        {
                            y = (j+(b-1)+cellsPerAxis[1])%cellsPerAxis[1];

                            // z loop for nearest neighbours
                            for (c=0;c<3;c++)
                            {
                                z = (k+(c-1)+cellsPerAxis[2])%cellsPerAxis[2];

                                // nn cell index
                                nn = x + y*cellsPerAxis[0] + z*cellsPerAxis[0]*cellsPerAxis[1];

                                at(m).neighbours[nnCount] = nn;
                                nnCount++;
                            }
                        }
                    }

                    // zero tally for each cell
                    at(m).tally = 0;
                    at(m).index = m;
                }
            }
        }
    }
    else
    {
        nNeighbours = 9;

        // loop over all cells x direction
        for (i=0;i<cellsPerAxis[0];i++)
        {
            // loop over all cells y direction
            for (j=0;j<cellsPerAxis[1];j++)
            {
                // cell index
                m = i + cellsPerAxis[0]*j;

                // resize neighbour list and particle arrays
                at(m).neighbours.resize(nNeighbours);
                at(m).particles.resize(maxParticles);

                nnCount = 0;

                // x loop for nearest neighbours
                for (a=0;a<3;a++)
                {
                    x = (i+(a-1)+cellsPerAxis[0])%cellsPerAxis[0];

                    // y loop for nearest neighbours
                    for (b=0;b<3;b++)
                    {
                        y = (j+(b-1)+cellsPerAxis[1])%cellsPerAxis[1];

                        // nn cell index
                        nn = x + y*cellsPerAxis[0];

                        at(m).neighbours[nnCount] = nn;
                        nnCount++;
                    }
                }

                // zero tally for each cell
                at(m).tally = 0;
                at(m).index = m;
            }
        }
    }
}

void CellList::reset()
{
    for (unsigned int i=0;i<nCells;i++) at(i).tally = 0;
}

int CellList::getCell(const Particle& particle)
{
    int cell,cellx,celly;

    cellx = int(particle.position[0]/cellSpacing[0]);
    celly = int(particle.position[1]/cellSpacing[1]);

    cell = cellx + celly*cellsPerAxis[0];

    if (dimension == 3)
    {
        int cellz = int(particle.position[2]/cellSpacing[2]);
        cell += cellz*cellsPerAxis[0]*cellsPerAxis[1];
    }

    return cell;
}

void CellList::initCell(int newCell, Particle& particle)
{
    // Add to new list
    at(newCell).particles[at(newCell).tally] = particle.index;
    particle.cell = newCell;
    particle.posCell = at(newCell).tally;
    at(newCell).tally++;

    if (at(newCell).tally == maxParticles)
    {
        std::cerr << "[ERROR] CellList: Maximum number of particles per cell exceeded!\n";
        exit(EXIT_FAILURE);
    }
}

void CellList::initCellList(std::vector<Particle>& particles)
{
    for (unsigned int i=0;i<particles.size();i++)
    {
        initCell(getCell(particles[i]), particles[i]);
    }
}

void CellList::updateCell(int newCell, Particle& particle, std::vector<Particle>& particles)
{
    // Remove from old list
    at(particle.cell).tally--;
    at(particle.cell).particles[particle.posCell] = at(particle.cell).particles[at(particle.cell).tally];
    particles[at(particle.cell).particles[at(particle.cell).tally]].posCell = particle.posCell;

    // Add to new list
    at(newCell).particles[at(newCell).tally] = particle.index;
    particle.cell = newCell;
    particle.posCell = at(newCell).tally;
    at(newCell).tally++;

    if (at(newCell).tally == maxParticles)
    {
        std::cerr << "[ERROR] CellList: Maximum number of particles per cell exceeded!\n";
        exit(EXIT_FAILURE);
    }
}

void CellList::setDimension(unsigned int dimension_)
{
    dimension = dimension_;
}

unsigned int CellList::getNeighbours() const
{
    return nNeighbours;
}
