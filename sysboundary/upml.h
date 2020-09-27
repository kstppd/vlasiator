/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef UPML_H
#define UPML_H

#include "sysboundarycondition.h"
#include<vector>
#include "fsgrid.hpp"
#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"

bool assingFSGridCells(FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,
                       FsGrid< fsgrids::technical, 2> & technicalGrid,
                       dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                       int width,int  side, int comp);


bool buildPMLGrid(FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,FsGrid< fsgrids::technical, 2> & technicalGrid, dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);


#endif
