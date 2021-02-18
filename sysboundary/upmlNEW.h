#pragma once
#include <iostream>
#include <fsgrid.hpp>
#include <cmath>
#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace ABC{
   
   
   
   
   class UPML{

      public:
         UPML(FsGrid< std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,FsGrid< fsgrids::technical, 2> & technicalGrid,dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
         bool resetPML(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt);
         int widthXP,widthXM,widthYP,widthYM,widthZP,widthZM,start;
         Real alpha;
         bool logcells;
      private:
         bool calculateParameters(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt);
         bool buildConductivity(FsGrid< std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt);
         bool getParameters();
         bool classifyCells(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,
                               FsGrid< fsgrids::technical, 2> & technicalGrid,
                               dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

   };









};
