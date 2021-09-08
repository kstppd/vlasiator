#pragma once
#include "../definitions.h"
#include "../common.h"
#include <fsgrid.hpp>
#include "../readparameters.h"
#include "../parameters.h"

namespace PerfectlyMatchedLayer{


   class History{

      public:

         History();
         bool push(const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid);
         bool getAvg(FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGridAvg);
         bool getDiffB(FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                        FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGridAvg,
                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid);
      private:

         int maxLength = Parameters::upmlHistoryLength;
         int active = Parameters::upmlHistory;
         std::vector<FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>> fgVec;



   };

   class UPML{

      public:

         UPML();
         void update(FsGrid<std::array<Real, fsgrids::upml::N_UPML>, FS_STENCIL_WIDTH>& fsUpml,
               FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, const Real dt);

      private:
         // void outputVarsPML(std::map<int, std::string>& output); 
         void getParameters();
         void addParameters();
         Real upmlWidth;
         bool pml_Xp,pml_Xm;
         bool pml_Yp,pml_Ym;
         bool pml_Zp,pml_Zm;
   };

};
