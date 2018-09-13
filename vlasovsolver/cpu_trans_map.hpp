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
#ifndef CPU_TRANS_MAP_H
#define CPU_TRANS_MAP_H

#include <vector>

#include "vec.h"
#include "../common.h"
#include "../spatial_cell.hpp"

struct setOfPencils {

   uint N; // Number of pencils in the set
   uint sumOfLengths;
   std::vector<uint> lengthOfPencils; // Lengths of pencils
   std::vector<CellID> ids; // List of cells
   std::vector<Realv> x,y; // x,y - position
   std::vector<bool> periodic;

   setOfPencils() {
      N = 0;
      sumOfLengths = 0;
   }

   void addPencil(std::vector<CellID> idsIn, Real xIn, Real yIn, bool periodicIn) {

      N += 1;
      sumOfLengths += idsIn.size();
      lengthOfPencils.push_back(idsIn.size());
      ids.insert(ids.end(),idsIn.begin(),idsIn.end());
      x.push_back(xIn);
      y.push_back(yIn);
      periodic.push_back(periodicIn);
   }

   std::vector<CellID> getIds(const uint pencilId) {

      std::vector<CellID> idsOut;
      
      if (pencilId > N) {
         return idsOut;
      }

      CellID ibeg = 0;
      for (uint i = 0; i < pencilId; i++) {
         ibeg += lengthOfPencils[i];
      }
      CellID iend = ibeg + lengthOfPencils[pencilId];
    
      for (uint i = ibeg; i <= iend; i++) {
         idsOut.push_back(ids[i]);
      }

      return idsOut;
   }

};

bool do_translate_cell(spatial_cell::SpatialCell* SC);
bool trans_map_1d(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& localPropagatedCells,
                  const std::vector<CellID>& remoteTargetCells,
                  const uint dimension,
                  const Realv dt,
                  const uint popID);
void update_remote_mapping_contribution(dccrg::Dccrg<spatial_cell::SpatialCell,
                                        dccrg::Cartesian_Geometry>& mpiGrid,
                                        const uint dimension,
                                        int direction,
                                        const uint popID);

void compute_spatial_source_neighbors(const dccrg::Dccrg<SpatialCell,
                                      dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      SpatialCell **neighbors);

void compute_spatial_target_neighbors(const dccrg::Dccrg<SpatialCell,
                                      dccrg::Cartesian_Geometry>& mpiGrid,
                                      const CellID& cellID,
                                      const uint dimension,
                                      SpatialCell **neighbors);
void copy_trans_block_data(SpatialCell** source_neighbors,
                           const vmesh::GlobalID blockGID,
                           Vec* values,
                           const unsigned char* const cellid_transpose,
                           const uint popID);

bool trans_map_1d_amr(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& localPropagatedCells,
                  const std::vector<CellID>& remoteTargetCells,
                  const uint dimension,
                  const Realv dt,
                  const uint popID);

setOfPencils buildPencilsWithNeighbors( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid, 
					setOfPencils &pencils, CellID startingId,
					std::vector<CellID> ids, uint dimension, 
					std::vector<uint> path);

void get_seed_ids(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID> &localPropagatedCells,
                  const uint dimension,
                  std::vector<CellID> &seedIds);

#endif
