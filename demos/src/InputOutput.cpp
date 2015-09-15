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

#include "InputOutput.h"

InputOutput::InputOutput() {}

void InputOutput::loadConfiguration(std::string fileName, Box& box,
    std::vector<Particle>& particles, CellList& cells, bool isIsotropic)
{
    std::ifstream dataFile;

    // Reset cell list.
    cells.reset();

    // Attempt to read data file.
    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        for (unsigned int i=0;i<particles.size();i++)
        {
            // Set particle index.
            particles[i].index = i;

            // Resize position and orientation vectors.
            particles[i].position.resize(box.dimension);
            particles[i].orientation.resize(box.dimension);

            // Load position.
            for (unsigned int j=0;j<box.dimension;j++)
                dataFile >> particles[i].position[j];

            if (!isIsotropic)
            {
                // Load orientation.
                for (unsigned int j=0;j<box.dimension;j++)
                    dataFile >> particles[i].orientation[j];
            }
            else
            {
                // Assign dummy orientation.
                for (unsigned int j=0;j<box.dimension;j++)
                    particles[i].orientation[j] = 1.0/sqrt(box.dimension);
            }

            // Enforce periodic boundary conditions.
            box.periodicBoundaries(particles[i].position);

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // Update cell list.
            cells.initCell(particles[i].cell, particles[i]);
        }
    }
    else
    {
        std::cerr << "[ERROR] InputOutput: Invalid restart file!\n";
        exit(EXIT_FAILURE);
    }

    // Close file stream.
    dataFile.close();
}

void InputOutput::saveConfiguration(std::string fileName, Box& box,
    std::vector<Particle>& particles, bool isIsotropic)
{
    // Create file pointer.
    FILE *pFile = fopen(fileName.c_str(), "w");

    for (unsigned int i=0;i<particles.size();i++)
    {
        // Write particle position.
        fprintf(pFile, "%5.4f %5.4f", particles[i].position[0], particles[i].position[1]);
        if (box.dimension == 3) fprintf(pFile, " %5.4f", particles[i].position[2]);

        // Write particle orientation.
        if (!isIsotropic)
        {
            fprintf(pFile, " %5.4f %5.4f", particles[i].orientation[0], particles[i].orientation[1]);
            if (box.dimension == 3) fprintf(pFile, " %5.4f", particles[i].orientation[2]);
        }

        // Terminate line.
        fprintf(pFile, "\n");
    }

    // Close file pointer.
    fclose(pFile);
}

void InputOutput::appendXyzTrajectory(unsigned int dimension, const std::vector<Particle>& particles, bool clearFile)
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

void InputOutput::vmdScript(const std::vector<double>& boxSize)
{
    FILE *pFile;

    pFile = fopen("vmd.tcl", "w");

    // Turn on lights 0 and 1.
    fprintf(pFile, "light 0 on\n");
    fprintf(pFile, "light 1 on\n");
    fprintf(pFile, "light 2 off\n");
    fprintf(pFile, "light 3 off\n");

    // Position the stage and axes.
    fprintf(pFile, "axes location off\n");
    fprintf(pFile, "stage location off\n");

    // Set orthographic projection.
    fprintf(pFile, "display projection orthographic\n");

    // Set drawing method to van der Waals radius.
    fprintf(pFile, "mol modstyle 0 0 VDW 1 30\n");

    // Set sensible atom radius.
    fprintf(pFile, "set sel [atomselect top \"name X\"]\n");
    fprintf(pFile, "atomselect0 set radius 0.4\n");

    // Set default particle to blue.
    fprintf(pFile, "color Name X blue\n");

    // Turn off depth cue.
    fprintf(pFile, "display depthcue off\n");

    // Define box boundaries.
    fprintf(pFile, "set minx 0\n");
    fprintf(pFile, "set maxx %5.4f\n", boxSize[0]);
    fprintf(pFile, "set miny 0\n");
    fprintf(pFile, "set maxy %5.4f\n", boxSize[1]);
    if (boxSize.size() == 3)
    {
        fprintf(pFile, "set minz 0\n");
        fprintf(pFile, "set maxz %5.4f\n", boxSize[2]);
    }
    else
    {
        fprintf(pFile, "set minz 0\n");
        fprintf(pFile, "set maxz 0\n");
    }

    // Set colours.
    fprintf(pFile, "draw materials off\n");
    fprintf(pFile, "draw color white\n");

    // Draw cube edges.
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$maxx $miny $minz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$minx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$minx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\"\n");

    // Rotate box.
    if (boxSize.size() == 3)
    {
        fprintf(pFile, "rotate x by -60\n");
        fprintf(pFile, "rotate y by -30\n");
        fprintf(pFile, "rotate z by -15\n");
    }

    fclose(pFile);
}

void InputOutput::vmdSpherocylinder(const std::vector<double>& boxSize)
{
    FILE *pFile;

    // Radius of spherocylinder cap.
    double r = 0.5*boxSize[0];

    pFile = fopen("vmd.tcl", "w");

    // Turn on lights 0 and 1.
    fprintf(pFile, "light 0 on\n");
    fprintf(pFile, "light 1 on\n");
    fprintf(pFile, "light 2 off\n");
    fprintf(pFile, "light 3 off\n");

    // Position the stage and axes.
    fprintf(pFile, "axes location off\n");
    fprintf(pFile, "stage location off\n");

    // Set orthographic projection.
    fprintf(pFile, "display projection orthographic\n");

    // Set drawing method to van der Waals radius.
    fprintf(pFile, "mol modstyle 0 0 VDW 1 30\n");

    // Set sensible atom radius.
    fprintf(pFile, "set sel [atomselect top \"name X\"]\n");
    fprintf(pFile, "atomselect0 set radius 0.4\n");

    // Set default particle to blue.
    fprintf(pFile, "color Name X blue\n");

    // Turn off depth cue.
    fprintf(pFile, "display depthcue off\n");

    // Set colours.
    fprintf(pFile, "draw materials off\n");
    fprintf(pFile, "draw color white\n");

    // Draw half spherocylinder.
    fprintf(pFile, "draw sphere { %3.2f %3.2f %3.2f } radius %3.2f resolution 25\n", r, r, r, r);
    fprintf(pFile, "draw cylinder { %3.2f %3.2f %3.2f } { %3.2f %3.2f %3.2f } radius %3.2f resolution 25 filled no\n", r, r, r, r, r, 0.5*boxSize[2], r);
    fprintf(pFile, "draw material Transparent\n");

    fprintf(pFile, "rotate x by -60\n");
    fprintf(pFile, "rotate y by -30\n");
    fprintf(pFile, "rotate z by -15\n");

    fclose(pFile);
}
