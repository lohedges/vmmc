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

#ifndef _CELLLIST_H
#define _CELLLIST_H

#include <cmath>

#include "Box.h"
#include "Particle.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

/*! \file CellList.h
*/

//! Structure containing attributes for an individual cell.
struct Cell
{
    unsigned int index;                     //!< cell index
    unsigned int tally;                     //!< number of particles in the cell
    std::vector<unsigned int> particles;    //!< indices of particles in the cell
    std::vector<unsigned int> neighbours;   //!< indices of nearest neighbour cells
};

//! Container class for storing a list of cells.
//! This class contains the main cell list that is manipulated by the simulation.
class CellList : public std::vector<Cell>
{
public:
    //! Default constructor.
    CellList();

    //! Constructor.
    /*! \param dimension
            The dimensionality of the simulation box.

        \param boxSize
            The size of the simulation box in each dimension.

        \param range
            Maximum interaction range.
     */
    CellList(unsigned int, const std::vector<double>&, double);

    //! Copy constructor. \param cells A reference to an existing CellList object.
    CellList& operator = (const CellList&);

    //! Copy constructor. \param cells A vector of existing Cell data structures.
    CellList& operator = (const std::vector<Cell>&);

    //! Initialise cell lists.
    /*! \param boxSize
            The size of the simulation box in each dimension.

        \param range
            Maximum interaction range.
     */
    void initialise(const std::vector<double>&, double);

    //! Reset cell lists (zero cell tallys).
    void reset();

    //! Get cell index for a particle.
    /*! \param particle
            Reference to a particle.

        \return
            The cell index.
     */
    int getCell(const Particle&);

    //! Initialise cell list for an individual particle.
    /*! \param newCell
            The index of the cell in which the particle is located.

        \param particle
            Reference to a particle.
     */
    void initCell(int, Particle&);

    //! Initialise cell list for all particles.
    /*! \param particles Reference to a vector of particles.
     */
    void initCellList(std::vector<Particle>&);

    //! Update cell list for an individual particle.
    /*! \param newCell
            The index of the cell in which the particle is located.

        \param particle
            Reference to a particle.

        \param particles
            Reference to a vector of particles.
     */
    void updateCell(int, Particle&, std::vector<Particle>&);

    //! Set the dimensionality of the cell list.
    /*! \param dimension_
            The dimensionality of the simulation.
     */
    void setDimension(unsigned int);

    //! Get the number of neighbours per cell.
    unsigned int getNeighbours() const;

private:
    unsigned int dimension;                     //!< dimension of the simulation box
    unsigned int nCells;                        //!< total number of cells
    unsigned int nNeighbours;                   //!< number of neighbours per cell
    unsigned int maxParticles;                  //!< maximum number of particles per cell
    std::vector<unsigned int> cellsPerAxis;     //!< number of cells per axis
    std::vector<double> cellSpacing;            //!< spacing between cells
};

#endif  /* _CELLLIST_H */
