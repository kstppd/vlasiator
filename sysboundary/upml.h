#pragma once
#include "../definitions.h"
#include "../common.h"
#include <fsgrid.hpp>

namespace PerfectlyMatchedLayer{


   class UPML{

      public:

         UPML();
         void update(FsGrid<std::array<Real, fsgrids::upml::N_UPML>, FS_STENCIL_WIDTH>& fsUpml,
              const FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, const Real dt);

      private:
         void outputVarsPML(std::map<int, std::string>& output); 
         bool getParameters();
         void addParameters();
         int upmlWidth;
         bool pml_Xp,pml_Xm;
         bool pml_Yp,pml_Ym;
         bool pml_Zp,pml_Zm;
   };

};