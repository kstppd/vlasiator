#pragma once
/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "../definitions.h"
#include "../logger.h"
#include "../mpiconversion.h"
#include "../object_wrapper.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../spatial_cell_wrapper.hpp"

// External API for Vlasiator
namespace ASTERIX {
void compress_vdfs_fourier_mlp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                               size_t number_of_spatial_cells);

std::size_t probe_network_size_in_bytes(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                        size_t number_of_spatial_cells);

void compress_vdfs_zfp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, size_t number_of_spatial_cells);
} // namespace ASTERIX
