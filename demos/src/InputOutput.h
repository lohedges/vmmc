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

#ifndef _INPUTOUTPUT_H
#define _INPUTOUTPUT_H

#include <string>

/*! \file InputOutput.h
    \brief A class for reading/writing data.
*/

// FORWARD DECLARATIONS

class  Box;
class  CellList;
struct Particle;

class InputOutput
{
public:
    //! Default constructor.
    InputOutput();

    //! Load a restart configuration from a plain text file.
    /*! \param fileName
            The path to the restart file.

        \param box
            A reference to the simulation box object.

        \param particles
            A reference to the particle container.

        \param cells
            A refence to the cell list.

        \param isIsotropic
            Whether the potential is isotropic (no orientation data).
     */
    void loadConfiguration(std::string, Box&, std::vector<Particle>&, CellList&, bool);

    //! Save a restart configuration to a plain text file.
    /*! \param fileName
            The path to the restart file.

        \param box
            A reference to the simulation box object.

        \param particles
            A reference to the particle container.

        \param isIsotropic
            Whether the potential is isotropic (no orientation data).
     */
    void saveConfiguration(std::string, Box&, std::vector<Particle>&, bool);

    //! Append a particle configuration to an existing xyz trajectory.
    /*! \param dimension
            The dimension of the simulation box.

        \param particles
            A vector of particles.

        \param clearFile
            Whether to clear the trajectory file before writing.
     */
    void appendXyzTrajectory(unsigned int, const std::vector<Particle>&, bool);

    //! Create a VMD TcL script to set the particle view and draw a bounding box.
    /*! \param boxSize
            The size of the simulation box in each dimension.
     */
    void vmdScript(const std::vector<double>&);

    //! Create a VMD TcL script to set the particle view and draw sphereocylinrical boundary.
    /*! \param boxSize
            The size of the simulation box in each dimension.
     */
    void vmdSpherocylinder(const std::vector<double>&);
};

#endif  /* _INPUTOUTPUT_H */
