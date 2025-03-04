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

// These tools  are fwd declared here and implemented at the end of the file for
// better clarity. They are not for external usage and as such they do not go
// into the header file

#include "../definitions.h"
#include "../logger.h"
#include "../mpiconversion.h"
#include "../object_wrapper.h"
#include "../spatial_cell_wrapper.hpp"
#include "../velocity_blocks.h"
#include <fsgrid.hpp>
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <span>
#include <stdexcept>
#include <unordered_map>

#define MLP_KEY 42

namespace ASTERIX {
struct VCoords {
   Real vx, vy, vz;
   VCoords operator+(const VCoords& other) { return {vx + other.vx, vy + other.vy, vz + other.vz}; }
   VCoords operator-(const VCoords& other) { return {vx - other.vx, vy - other.vy, vz - other.vz}; }
};

struct OrderedVDF {

   struct OrderedVDFHeader{
      std::size_t ignore_bytes;
      std::size_t octree_bytes;
   };
   
   std::vector<vmesh::GlobalID> blocks_to_ignore;
   std::size_t sparse_vdf_bytes = {0};
   std::vector<Realf> vdf_vals;
   std::array<Real, 6> v_limits;     // vx_min,vy_min,vz_min,vx_max,vy_max,vz_max
   std::array<std::size_t, 3> shape; // x,y,z
   std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
      return i * (shape[1] * shape[2]) + j * shape[2] + k;
   }

   Realf& at(std::size_t i, std::size_t j, std::size_t k) noexcept { return vdf_vals.at(index(i, j, k)); }

   const Realf& at(std::size_t i, std::size_t j, std::size_t k) const noexcept { return vdf_vals.at(index(i, j, k)); }

   bool save_to_file(const char* filename) const noexcept {
      std::ofstream file(filename, std::ios::out | std::ios::binary);
      if (!file) {
         std::cerr << "Could not open file for writting! Exiting!" << std::endl;
         return false;
      }
      file.write((char*)shape.data(), 3 * sizeof(size_t));
      if (!file) {
         std::cerr << "Error writing shape data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.data(), vdf_vals.size() * sizeof(Realf));
      if (!file) {
         std::cerr << "Error writing vdf_vals data to file!" << std::endl;
         return false;
      }
      return true;
   }
};

struct UnorderedVDF {
   std::vector<Realf> vdf_vals;
   std::vector<std::array<Real, 3>> vdf_coords;
   std::array<Real, 6> v_limits{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                                std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                                std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   bool save_to_file(const char* filename) const noexcept {
      std::ofstream file(filename, std::ios::out | std::ios::binary);
      if (!file) {
         std::cerr << "Could not open file for writting! Exiting!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.size(), sizeof(size_t));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)v_limits.data(), 6 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_coords.data(), vdf_coords.size() * 3 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing vdf_coords data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.data(), vdf_vals.size() * sizeof(Realf));
      if (!file) {
         std::cerr << "Error writing vdf_vals data to file!" << std::endl;
         return false;
      }
      return true;
   }
};

template<typename T>
class PhaseSpace6D{
   public:
      PhaseSpace6D(FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                   dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, uint popID) {

         const auto local_start=technicalGrid.getLocalStart();
         const auto local_size=technicalGrid.getLocalSize();
         const auto coords0 = technicalGrid.getPhysicalCoords(0,0,0);
         const auto coords1 = technicalGrid.getPhysicalCoords(local_size[0],local_size[1],local_size[2]);
         
         rv_lims[0] =coords0[0]; 
         rv_lims[1] =coords0[1]; 
         rv_lims[2] =coords0[2]; 
         rv_lims[6] =coords1[0]; 
         rv_lims[7] =coords1[1]; 
         rv_lims[8] =coords1[2]; 
         T fsum = 0.0;
         Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;
         const auto gridDims(technicalGrid.getLocalSize());
         for (FsGridTools::FsIndex_t k = 0; k < gridDims[2]; k+=1) {
            for (FsGridTools::FsIndex_t j = 0; j < gridDims[1]; j+=1) {
               for (FsGridTools::FsIndex_t i = 0; i < gridDims[0]; i+=1) {
                  const std::array<FsGridTools::FsIndex_t, 3> globalIndices = technicalGrid.getGlobalIndices(i, j, k);
                  const dccrg::Types<3>::indices_t indices = {
                      {(uint64_t)globalIndices[0], (uint64_t)globalIndices[1], (uint64_t)globalIndices[2]}};
                  CellID dccrgCell =
                      mpiGrid.get_existing_cell(indices, 0, mpiGrid.mapping.get_maximum_refinement_level());
                  SpatialCell* sc = mpiGrid[dccrgCell];
                  vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = sc->get_velocity_blocks(popID);
                  const std::size_t total_blocks = blockContainer.size();
                  const Realf* data = blockContainer.getData();
                  const Real* blockParams = sc->get_block_parameters(popID);
                  const auto coords = technicalGrid.getPhysicalCoords(i, j, k);
                  for (std::size_t n = 0; n < total_blocks; ++n) {
                     const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
                     const Realf* vdf_data = &data[n * WID3];
                     for (uint kk = 0; kk < WID; kk+=1) {
                        for (uint jj = 0; jj < WID; jj+=1) {
                           for (uint ii = 0; ii < WID; ii+=1) {
                              const T vx = (bp[BlockParams::VXCRD] + (ii + 0.5) * bp[BlockParams::DVX]);
                              const T vy = (bp[BlockParams::VYCRD] + (jj + 0.5) * bp[BlockParams::DVY]);
                              const T vz = (bp[BlockParams::VZCRD] + (kk + 0.5) * bp[BlockParams::DVZ]);
                              const T value = static_cast<T>(vdf_data[cellIndex(ii, jj, kk)]);
                              const T scaled_value =
                                  std::abs(std::log10(std::max(value, static_cast<T>(0.1f * sparse))));
                               
                              // const T scaled_value =
                              //     std::abs(std::max(value, static_cast<T>(0.1f * sparse)));                               
                              this->space.push_back(coords[0]);
                              this->space.push_back(coords[1]);
                              this->space.push_back(coords[2]);
                              this->space.push_back(vx);
                              this->space.push_back(vy);
                              this->space.push_back(vz);
                              this->f.push_back(scaled_value);
                              fsum += scaled_value;
                              this->nrows++;

                              rv_lims[3] = std::min(rv_lims[3], vx);
                              rv_lims[4] = std::min(rv_lims[4], vy);
                              rv_lims[5] = std::min(rv_lims[5], vz);
                              rv_lims[9] = std::max(rv_lims[9], vx);
                              rv_lims[10] = std::max(rv_lims[10], vy);
                              rv_lims[11] = std::max(rv_lims[11], vz);
                           }
                        }
                     } // over vspace
                  }    // over blocks
               }
            }
         }// over real space         
         //Mean and std
         norms.mean=fsum/this->f.size();
         T diff2={0.0};
         for (auto val : f) {
            diff2 += (val - norms.mean) * (val - norms.mean);
         }
         norms.std=std::sqrt(diff2 / f.size());

         // auto min_f=*std::min_element(f.begin(),f.end());
         // auto max_f=*std::max_element(f.begin(),f.end());
         // auto min_v=*std::min_element(space.begin(),space.end());
         // auto max_v=*std::max_element(space.begin(),space.end());
         // std::cerr<<"Normed-> "<<min_f<<" "<<max_f<<" "<<min_v<<" "<<max_v<<std::endl;
         // for (int i=0;i<12;++i){
         //    std::cerr<<rv_lims[i]<<", ";
         // }
         // std::cerr<<std::endl;
         normalize();
      }
      
   PhaseSpace6D(const std::vector<char>&serialized_state,FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, uint popID) {

      const auto local_start=technicalGrid.getLocalStart();
      const auto local_size=technicalGrid.getLocalSize();
      const auto coords0 = technicalGrid.getPhysicalCoords(0,0,0);
      const auto coords1 = technicalGrid.getPhysicalCoords(local_size[0],local_size[1],local_size[2]);
      deserialize_from((unsigned char*)serialized_state.data());
      
      const auto vxmin = this->rv_lims[3] + (T(0.5) * (this->rv_lims[9] - this->rv_lims[3]));
      const auto vymin = this->rv_lims[4] + (T(0.5) * (this->rv_lims[10] - this->rv_lims[4]));
      const auto vzmin = this->rv_lims[5] + (T(0.5) * (this->rv_lims[11] - this->rv_lims[5]));

      const auto vxmax = this->rv_lims[3] + ((T(1.5)) * (this->rv_lims[9] - this->rv_lims[3]));
      const auto vymax = this->rv_lims[4] + ((T(1.5)) * (this->rv_lims[10] - this->rv_lims[4]));
      const auto vzmax = this->rv_lims[5] + ((T(1.5)) * (this->rv_lims[11] - this->rv_lims[5]));

      const auto DV = 6.0e4;
      const auto gridDims(technicalGrid.getLocalSize());
      for (FsGridTools::FsIndex_t k = 0; k < gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j = 0; j < gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i = 0; i < gridDims[0]; i++) {
               const auto rcoords = technicalGrid.getPhysicalCoords(i, j, k);
               for (float vz = vzmin; vz < vzmax; vz += DV) {
                  for (float vy = vymin; vy < vymax; vy += DV) {
                     for (float vx = vxmin; vx < vxmax; vx += DV) {
                        this->space.push_back(rcoords[0]);
                        this->space.push_back(rcoords[1]);
                        this->space.push_back(rcoords[2]);
                        this->space.push_back(vx);
                        this->space.push_back(vy);
                        this->space.push_back(vz);
                        this->f.push_back(0.0);
                        this->nrows++;
                     }
                  }
               } // over vspace
            }
         }
      } // over real space
         // auto min_f=*std::min_element(f.begin(),f.end());
         // auto max_f=*std::max_element(f.begin(),f.end());
         // auto min_v=*std::min_element(space.begin(),space.end());
         // auto max_v=*std::max_element(space.begin(),space.end());
         // std::cerr<<"Normed-> "<<min_f<<" "<<max_f<<" "<<min_v<<" "<<max_v<<std::endl;
         // for (int i=0;i<12;++i){
         //    std::cerr<<rv_lims[i]<<", ";
         // }
         // std::cerr<<std::endl;
      normalize();
   }
      
   struct Norms{
      T mean={T(0)};
      T std={T(0)};
   };
   


   void normalize()noexcept{
      for (auto& val:f){
         val=(val-norms.mean)/norms.std;
      }
      
      //Normalize rv coords too
      for (std::size_t row =0; row<this->nrows; ++row){
         for (std::size_t col=0; col<this->ncols; ++col){
            std::size_t index=row*ncols+col; 
            this->space[index]=(this->space[index]-this->rv_lims[col])/(this->rv_lims[col+ncols]-this->rv_lims[col]) - T(0.5);
         }
      }
   }

   void denormalize(T sparse) noexcept {
      for (auto& val : f) {
         val = val * norms.std + norms.mean;
         if (val<sparse){
            val=0.0;
         }
      }

      for (std::size_t row = 0; row < this->nrows; ++row) {
         for (std::size_t col = 0; col < this->ncols; ++col) {
            std::size_t index = row * ncols + col;
            this->space[index] =
                (this->space[index] + T(0.5)) * (this->rv_lims[col + ncols] - this->rv_lims[col]) + this->rv_lims[col];
         }
      }
   }

   std::size_t serialized_bytes_size()const noexcept{
      return sizeof(Norms)+12*sizeof(T)+mlp_representation_nbytes;
   }
   
   void serialize_into(unsigned char* buffer) const {
      std::size_t write_index = 0;
      std::memcpy(&buffer[write_index], &norms.mean, sizeof(Norms));
      write_index += sizeof(Norms);
      std::memcpy(&buffer[write_index], &rv_lims[0], 12 * sizeof(T));
      write_index += 12 * sizeof(T);
      std::memcpy(&buffer[write_index], mlp_representation, mlp_representation_nbytes);
      write_index += mlp_representation_nbytes;
   }
   
   void deserialize_from(unsigned char* buffer) {
      std::size_t read_index = 0;
      std::memcpy(&norms.mean,&buffer[read_index], sizeof(Norms));
      read_index += sizeof(Norms);
      std::memcpy(&rv_lims[0],&buffer[read_index], 12 * sizeof(T));
      read_index += 12 * sizeof(T);
   }

   //This holds a 2D representation of the 6D hypercube coordinates where each rows is [x,y,z,vx,vy,vz]
   std::vector<T> space;
   //This holds a 1D representation of the 6D hypercube's phase space value where each rows is [f]
   std::vector<T> f;
   std::size_t nrows={0};
   std::size_t ncols = {6};
   std::size_t cols = {0};
   T* mlp_representation=nullptr;
   std::size_t mlp_representation_nbytes=0;
   std::array<T, 12> rv_lims = {
       std::numeric_limits<T>::max(),    std::numeric_limits<T>::max(),    std::numeric_limits<T>::max(),
       std::numeric_limits<T>::max(),    std::numeric_limits<T>::max(),    std::numeric_limits<T>::max(),
       std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(),
       std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};

   Norms norms = {};
};

struct VDFUnion {

   struct MinMaxValues {
      Realf min = std::numeric_limits<Real>::lowest();
      Realf max = std::numeric_limits<Real>::max();
      Realf mean = 0.0;
   };

   struct SerializedVDFUnionHeader {
      std::size_t key;
      std::size_t total_size;
      std::size_t rows;
      std::size_t cols;
      std::size_t n_weights;
   };

   std::size_t nrows = 0, ncols = 0;
   std::vector<MinMaxValues> norms;
   std::vector<CellID> cids;
   std::vector<std::array<Real, 3>> vcoords_union;
   std::vector<VCoords> vbulk_union;
   std::vector<Realf> vspace_union;
   std::unordered_map<vmesh::LocalID, std::size_t> map;
   std::size_t size_in_bytes;
   double* network_weights = nullptr;
   std::size_t n_weights;
   std::array<Real, 6> v_limits{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                                std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                                std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   std::size_t total_serialized_size_bytes() const {
      return sizeof(SerializedVDFUnionHeader) + cids.size() * sizeof(CellID) + norms.size() * sizeof(MinMaxValues) +
             vbulk_union.size() * sizeof(VCoords) + vcoords_union.size() * 3 * sizeof(Real) + 6*sizeof(Real)+
             n_weights * sizeof(double) + map.size() * sizeof(std::pair<vmesh::LocalID, std::size_t>);
      ;
   }

   void print(){
         int myRank;
   int mpiProcs;
   MPI_Comm_size(MPI_COMM_WORLD, &mpiProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   std::cout<<"============"<<myRank<<"================"<<std::endl;
      printf("Size %zu %zu\n",nrows,ncols);
      printf("Network size %zu\n",n_weights);
      printf("First weight %f\n ",network_weights[0]);
      printf("Last weight %f\n ",network_weights[n_weights-1]);
      
   for (int i=0;i<6;++i){
      std::cout<<v_limits[i]<<", ";
   }
   std::cout<<std::endl;
   std::cout<<"=========================="<<std::endl;
      
      // for (std::size_t i=0; i< nrows;++i){
      //    printf("%zu %f %f %f\n",i,vcoords_union[i][0],vcoords_union[i][1],vcoords_union[i][2]);
      // }
   }


   
   void serialize_into(unsigned char* buffer) const {
      SerializedVDFUnionHeader header;
      header.key = MLP_KEY;
      header.total_size = total_serialized_size_bytes();
      header.rows = nrows;
      header.cols = ncols;
      header.n_weights = n_weights;
      std::size_t write_index = 0;

      std::memcpy(&buffer[write_index], &header, sizeof(SerializedVDFUnionHeader));
      write_index += sizeof(SerializedVDFUnionHeader);

      std::memcpy(&buffer[write_index], &cids[0], cids.size() * sizeof(CellID));
      write_index += cids.size() * sizeof(CellID);

      std::memcpy(&buffer[write_index], &norms[0], norms.size() * sizeof(MinMaxValues));
      write_index += norms.size() * sizeof(MinMaxValues);

      std::memcpy(&buffer[write_index], &vbulk_union[0], vbulk_union.size() * sizeof(VCoords));
      write_index += vbulk_union.size() * sizeof(VCoords);
      
      std::memcpy(&buffer[write_index], &v_limits[0],6 * sizeof(Real));
      write_index += 6 * sizeof(Real);

      std::memcpy(&buffer[write_index], &vcoords_union[0], vcoords_union.size() * 3 * sizeof(Real));
      write_index += vcoords_union.size() * 3 * sizeof(Real);

      std::memcpy(&buffer[write_index], &network_weights[0], n_weights * sizeof(double));
      write_index += n_weights * sizeof(double);

      for (const auto& kval : map) {
         std::memcpy(&buffer[write_index], &kval, sizeof(std::pair<vmesh::LocalID, std::size_t>));
         write_index += sizeof(std::pair<vmesh::LocalID, std::size_t>);
      }

      assert(header.total_size == write_index);
   }

   void deserialize_from(const unsigned char* buffer) {
      const SerializedVDFUnionHeader* const header = reinterpret_cast<const SerializedVDFUnionHeader*>(&buffer[0]);
      assert(header->key = MLP_KEY && "Blame Kostis Papadakis for this!");

      // Inflate vspave union
      vspace_union.resize(header->cols * header->rows);
      ncols = header->cols;
      nrows = header->rows;

      // Recover cids in this union;
      std::size_t read_index = sizeof(SerializedVDFUnionHeader);
      std::size_t cids_size = header->cols;
      cids.resize(cids_size);

      std::memcpy(cids.data(), &buffer[read_index], cids_size * sizeof(CellID));
      read_index += cids_size * sizeof(CellID);

      std::size_t norms_size = cids_size;
      norms.resize(cids_size);
      std::memcpy(norms.data(), &buffer[read_index], norms_size * sizeof(MinMaxValues));
      read_index += norms_size * sizeof(MinMaxValues);

      std::size_t vbulk_size = cids_size;
      vbulk_union.resize(vbulk_size);
      std::memcpy(vbulk_union.data(), &buffer[read_index], vbulk_size * sizeof(VCoords));
      read_index += vbulk_size * sizeof(VCoords);
      
      std::memcpy(&v_limits[0], &buffer[read_index], 6 * sizeof(Real));
      read_index +=  6 * sizeof(Real);

      std::size_t vcoords_size = header->rows;
      vcoords_union.resize(vcoords_size);
      std::memcpy(vcoords_union.data(), &buffer[read_index], vcoords_size * 3 * sizeof(Real));
      read_index += vcoords_size * 3 * sizeof(Real);

      if (network_weights != nullptr) {
         free(network_weights);
      }

      network_weights = (double*)malloc(header->n_weights * sizeof(double));
      n_weights = header->n_weights;
      std::memcpy(network_weights, &buffer[read_index], header->n_weights * sizeof(double));
      read_index += n_weights * sizeof(double);

      while (read_index < header->total_size) {
         const std::pair<vmesh::LocalID, std::size_t>* kval =
             reinterpret_cast<const std::pair<vmesh::LocalID, std::size_t>*>(&buffer[read_index]);
         map[kval->first] = kval->second;
         read_index += sizeof(std::pair<vmesh::LocalID, std::size_t>);
      }

      assert(read_index == total_size && "Size mismatch while reading in serialized VDF Union!");
   }

   std::size_t index_2d(std::size_t row, std::size_t col) { return row * ncols + col; };

   void sparsify(Realf sparse) {
      std::for_each(vspace_union.begin(), vspace_union.end(), [sparse](Realf& x) {
         if (x - sparse < 0.0) {
            x = 0.0;
         }
      });
   }

   void normalize_union() {
      const std::size_t nVDFS = ncols;
      norms.resize(nVDFS);

      for (std::size_t v = 0; v < nVDFS; ++v) {
         Realf sum = 0;
         for (std::size_t i = 0; i < nrows; ++i) {
            sum += vspace_union[index_2d(i, v)];
         }
         Realf mean_val = sum / nrows;

         for (std::size_t i = 0; i < nrows; ++i) {
            vspace_union[index_2d(i, v)] -= mean_val;
         }

         Realf min_val = std::numeric_limits<Realf>::max();
         Realf max_val = std::numeric_limits<Realf>::lowest();
         for (std::size_t i = 0; i < nrows; ++i) {
            min_val = std::min(min_val, vspace_union[index_2d(i, v)]);
            max_val = std::max(max_val, vspace_union[index_2d(i, v)]);
         }
         Realf range = max_val - min_val;
         for (std::size_t i = 0; i < nrows; ++i) {
            vspace_union[index_2d(i, v)] = (vspace_union[index_2d(i, v)] - min_val) / range;
         }
         norms[v] = MinMaxValues{.min = min_val, .max = max_val, .mean = mean_val};
      }
   }

   void unormalize_union() {
      const std::size_t nVDFS = ncols;
      for (std::size_t v = 0; v < nVDFS; ++v) {
         const Realf max_val = norms[v].max;
         const Realf min_val = norms[v].min;
         const Realf mean_val = norms[v].mean;
         const Realf range = max_val - min_val;
         for (std::size_t i = 0; i < nrows; ++i) {
            vspace_union[index_2d(i, v)] = vspace_union[index_2d(i, v)] * range + min_val + mean_val;
         }
      }
   }

   void scale(Realf sparse) {
      std::for_each(vspace_union.begin(), vspace_union.end(),
                    [sparse](Realf& value) { value = std::abs(std::log10(std::max(value, 0.1f*sparse))); });
   }

   void unscale(Realf sparse) {
      std::for_each(vspace_union.begin(), vspace_union.end(),
                    [](Realf& value) { value = std::pow(10.0, -1.0 * value); });
   }

   bool save_to_file(const char* filename) const noexcept {
      std::ofstream file(filename, std::ios::out | std::ios::binary);
      if (!file) {
         std::cerr << "Could not open file for writting! Exiting!" << std::endl;
         return false;
      }

      file.write((char*)vspace_union.size(), sizeof(size_t));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)v_limits.data(), 6 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)vcoords_union.data(), vcoords_union.size() * 3 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing vdf_coords data to file!" << std::endl;
         return false;
      }

      file.write((char*)vspace_union.data(), vspace_union.size() * sizeof(Realf));
      if (!file) {
         std::cerr << "Error writing vdf_vals data to file!" << std::endl;
         return false;
      }
      return true;
   }
};

auto extract_pop_vdf_from_spatial_cell(SpatialCell* sc, uint popID) -> UnorderedVDF;

auto extract_union_pop_vdfs_from_cids(const std::span<const CellID> cids, uint popID,
                                      const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                      bool center_vdfs = false) -> VDFUnion;

auto extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(SpatialCell* sc, uint popID, int zoom) -> OrderedVDF;

constexpr auto isPow2(std::unsigned_integral auto val) -> bool { return (val & (val - 1)) == 0; };

auto overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const std::vector<Realf>& new_vspace) -> void;

auto overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const OrderedVDF& vdf) -> void;

auto overwrite_cellids_vdfs(const std::span<const CellID> cids, uint popID,
                            dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                            const std::vector<std::array<Real, 3>>& vcoords, const std::vector<Realf>& vspace_union,
                            const std::unordered_map<vmesh::LocalID, std::size_t>& map_exists_id) -> void;

auto overwrite_cellids_vdf_single_cell(const std::span<const CellID> cids, uint popID, SpatialCell* sc, size_t cc,
                                     const std::vector<std::array<Real, 3>>& vcoords,
                                     const std::vector<Realf>& vspace_union,
                                     const std::unordered_map<vmesh::LocalID, std::size_t>& map_exists_id)->void; 

auto dump_vdf_to_binary_file(const char* filename, CellID cid) -> void;

auto dump_vdf_to_binary_file(const char* filename, CellID cid,
                             dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) -> void;

// https://en.wikipedia.org/wiki/Entropy_(information_theory)
template <typename T>
requires(std::is_same_v<T, float> || std::is_same_v<T, double>) auto shannon_entropy(const std::vector<T>& data) -> T {
   const std::size_t sz = data.size();
   if (sz == 0) {
      return 0.0;
   }

   using key_t = std::conditional_t<std::is_same_v<T, float>, uint32_t, uint64_t>;
   std::unordered_map<key_t, int> frequency;
   for (std::size_t i = 0; i < sz; ++i) {
      frequency[*(reinterpret_cast<const key_t*>(&data[i]))]++;
   }
   T entropy = 0.0;
   for (const auto& [byte, count] : frequency) {
      T pk = static_cast<T>(count) / sz;
      entropy -= pk * std::log2(pk);
   }
   return entropy;
}

template <typename T>
requires(std::is_same_v<T, float> ||
         std::is_same_v<T, double>) auto theoritical_lossless_compression_ratio(const std::vector<T>& data,
                                                                                std::size_t bits) -> T {
   T entorpy = shannon_entropy(data);
   return static_cast<T>(bits) / entorpy;
}

template <typename NetworkType>
requires(std::is_same_v<NetworkType, float> || std::is_same_v<NetworkType, double>) auto calculate_total_size_bytes(
    const std::vector<std::size_t>& architecture, std::size_t fourier_order, std::size_t output_dim) -> std::size_t {
   if (architecture.empty()) {
      throw std::runtime_error("Architecture cannot be empty.");
   }
   std::size_t input_dim = 2 * fourier_order;
   std::size_t total_size = 0;
   total_size += input_dim * architecture[0];
   total_size += architecture[0];

   for (std::size_t i = 1; i < architecture.size(); ++i) {
      total_size += architecture[i - 1] * architecture[i];
      total_size += architecture[i];
   }

   total_size += architecture.back() * output_dim;
   total_size += output_dim;

   return total_size * sizeof(NetworkType);
}

template <typename NetworkType>
requires(std::is_same_v<NetworkType, float> || std::is_same_v<NetworkType, double>) auto calculate_hidden_neurons(
    std::size_t N_input, std::size_t N_output, std::size_t num_hidden_layers, std::size_t target_size)
    -> std::vector<std::size_t> {
   std::vector<std::size_t> neurons(num_hidden_layers + 2); // 2 input and output
   neurons[0] = N_input;
   neurons[num_hidden_layers + 1] = N_output;

   // We guess this heyuristically
   std::size_t initial_hidden_size = 1;
   for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
      neurons[i] = initial_hidden_size;
   }
   std::size_t current_size = calculate_total_size_bytes<NetworkType>(neurons);

   while (current_size < target_size) {
      for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
         neurons[i]++;
      }
      current_size = calculate_total_size_bytes<NetworkType>(neurons);
   }

   while (current_size > target_size) {
      for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
         if (neurons[i] > 1) {
            neurons[i]--;
         }
      }
      current_size = calculate_total_size_bytes<NetworkType>(neurons);
   }
   return neurons;
}
Real get_Non_MaxWellianity(const SpatialCell* cell, uint popID);

} // namespace ASTERIX
