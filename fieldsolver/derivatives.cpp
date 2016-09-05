/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute
*/

#include <cstdlib>

#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h. Uses RHO, RHOV[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method, and RHO_DT2, RHOV[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   cint i,
   cint j,
   cint k,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool& doMoments
) {
   std::array<Real, fsgrids::dperb::N_DPERB> * dPerB = dPerBGrid.get(i,j,k);
   std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
   
   std::array<Real, fsgrids::moments::N_MOMENTS> * leftMoments = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD> * leftPerB = NULL;
   std::array<Real, fsgrids::moments::N_MOMENTS> * centMoments = momentsGrid.get(i,j,k);
   std::array<Real, fsgrids::bfield::N_BFIELD> * centPerB = perBGrid.get(i,j,k);
   #ifdef DEBUG_SOLVERS
   if (centMoments->at(fsgrids::moments::RHO) <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (centMoments->at(fsgrids::moments::RHO) < 0 ? " Negative" : " Zero") << " density in spatial cell " << cellID
         << std::endl;
      abort();
   }
   #endif
   std::array<Real, fsgrids::moments::N_MOMENTS> * rghtMoments = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * rghtPerB = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * botRght = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topLeft = NULL;
   std::array<Real, fsgrids::bfield::N_BFIELD>  * topRght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i-1,j,k);
      rghtPerB = perBGrid.get(i+1,j,k);
      if (doMoments) {
         leftMoments = momentsGrid.get(i-1,j,k);
         rghtMoments = momentsGrid.get(i+1,j,k);
      }
      
      if (doMoments) {
         dMoments->at(fsgrids::dmoments::drhodx) = limiter(leftMoments->at(fsgrids::moments::RHO),centMoments->at(fsgrids::moments::RHO),rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dp11dx) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dx) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dx) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));

         dMoments->at(fsgrids::dmoments::dVxdx)  = limiter(leftMoments->at(fsgrids::moments::RHOVX), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVX), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVX), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVydx)  = limiter(leftMoments->at(fsgrids::moments::RHOVY), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVY), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVY), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVzdx)  = limiter(leftMoments->at(fsgrids::moments::RHOVZ), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVZ), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVZ), rghtMoments->at(fsgrids::moments::RHO));
         }
      dPerB->at(fsgrids::dperb::dPERBydx)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
      dPerB->at(fsgrids::dperb::dPERBzdx)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));
      if(Parameters::ohmHallTerm < 2) {
         dPerB->at(fsgrids::dperb::dPERBydxx) = 0.0;
         dPerB->at(fsgrids::dperb::dPERBzdxx) = 0.0;
      } else {
         dPerB->at(fsgrids::dperb::dPERBydxx) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
         dPerB->at(fsgrids::dperb::dPERBzdxx) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 0);
      }
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j-1,k);
      rghtPerB = perBGrid.get(i,j+1,k);
      if (doMoments) {
         leftMoments = momentsGrid.get(i,j-1,k);
         rghtMoments = momentsGrid.get(i,j+1,k);
      }
      
      if (doMoments) {
         dMoments->at(fsgrids::dmoments::drhody) = limiter(leftMoments->at(fsgrids::moments::RHO),centMoments->at(fsgrids::moments::RHO),rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dp11dy) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dy) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dy) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdy)  = limiter(leftMoments->at(fsgrids::moments::RHOVX), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVX), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVX), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVydy)  = limiter(leftMoments->at(fsgrids::moments::RHOVY), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVY), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVY), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVzdy)  = limiter(leftMoments->at(fsgrids::moments::RHOVZ), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVZ), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVZ), rghtMoments->at(fsgrids::moments::RHO));
      }
      dPerB->at(fsgrids::dperb::dPERBxdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
      dPerB->at(fsgrids::dperb::dPERBzdy)  = limiter(leftPerB->at(fsgrids::bfield::PERBZ),centPerB->at(fsgrids::bfield::PERBZ),rghtPerB->at(fsgrids::bfield::PERBZ));

      if(Parameters::ohmHallTerm < 2) {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = 0.0;
         dPerB->at(fsgrids::dperb::dPERBzdyy) = 0.0;
      } else {
         dPerB->at(fsgrids::dperb::dPERBxdyy) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
         dPerB->at(fsgrids::dperb::dPERBzdyy) = leftPerB->at(fsgrids::bfield::PERBZ) + rghtPerB->at(fsgrids::bfield::PERBZ) - 2.0*centPerB->at(fsgrids::bfield::PERBZ);
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      
      leftPerB = perBGrid.get(i,j,k-1);
      rghtPerB = perBGrid.get(i,j,k+1);
      if (doMoments) {
         leftMoments = momentsGrid.get(i,j,k-1);
         rghtMoments = momentsGrid.get(i,j,k+1);
      }
      
      if (doMoments) {
         dMoments->at(fsgrids::dmoments::drhodz) = limiter(leftMoments->at(fsgrids::moments::RHO),centMoments->at(fsgrids::moments::RHO),rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dp11dz) = limiter(leftMoments->at(fsgrids::moments::P_11),centMoments->at(fsgrids::moments::P_11),rghtMoments->at(fsgrids::moments::P_11));
         dMoments->at(fsgrids::dmoments::dp22dz) = limiter(leftMoments->at(fsgrids::moments::P_22),centMoments->at(fsgrids::moments::P_22),rghtMoments->at(fsgrids::moments::P_22));
         dMoments->at(fsgrids::dmoments::dp33dz) = limiter(leftMoments->at(fsgrids::moments::P_33),centMoments->at(fsgrids::moments::P_33),rghtMoments->at(fsgrids::moments::P_33));
         dMoments->at(fsgrids::dmoments::dVxdz)  = limiter(leftMoments->at(fsgrids::moments::RHOVX), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVX), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVX), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVydz)  = limiter(leftMoments->at(fsgrids::moments::RHOVY), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVY), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVY), rghtMoments->at(fsgrids::moments::RHO));
         dMoments->at(fsgrids::dmoments::dVzdz)  = limiter(leftMoments->at(fsgrids::moments::RHOVZ), leftMoments->at(fsgrids::moments::RHO),
                                                       centMoments->at(fsgrids::moments::RHOVZ), centMoments->at(fsgrids::moments::RHO),
                                                       rghtMoments->at(fsgrids::moments::RHOVZ), rghtMoments->at(fsgrids::moments::RHO));
      }
      dPerB->at(fsgrids::dperb::dPERBxdz)  = limiter(leftPerB->at(fsgrids::bfield::PERBX),centPerB->at(fsgrids::bfield::PERBX),rghtPerB->at(fsgrids::bfield::PERBX));
      dPerB->at(fsgrids::dperb::dPERBydz)  = limiter(leftPerB->at(fsgrids::bfield::PERBY),centPerB->at(fsgrids::bfield::PERBY),rghtPerB->at(fsgrids::bfield::PERBY));
      if(Parameters::ohmHallTerm < 2) {
         dPerB->at(fsgrids::dperb::dPERBxdzz) = 0.0;
         dPerB->at(fsgrids::dperb::dPERBydzz) = 0.0;
      } else {
         dPerB->at(fsgrids::dperb::dPERBxdzz) = leftPerB->at(fsgrids::bfield::PERBX) + rghtPerB->at(fsgrids::bfield::PERBX) - 2.0*centPerB->at(fsgrids::bfield::PERBX);
         dPerB->at(fsgrids::dperb::dPERBydzz) = leftPerB->at(fsgrids::bfield::PERBY) + rghtPerB->at(fsgrids::bfield::PERBY) - 2.0*centPerB->at(fsgrids::bfield::PERBY);
      }
   } else {
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 2);
      }
   }
   
   if (Parameters::ohmHallTerm < 2) {
      dPerB->at(fsgrids::dperb::dPERBxdyz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBydxz) = 0.0;
      dPerB->at(fsgrids::dperb::dPERBzdxy) = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j-1,k);
         botRght = perBGrid.get(i+1,j-1,k);
         topLeft = perBGrid.get(i-1,j+1,k);
         topRght = perBGrid.get(i+1,j+1,k);
         
         dPerB->at(fsgrids::dperb::dPERBzdxy) = FOURTH * (botLeft->at(fsgrids::bfield::PERBZ) + topRght->at(fsgrids::bfield::PERBZ) - botRght->at(fsgrids::bfield::PERBZ) - topLeft->at(fsgrids::bfield::PERBZ));
         
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 3);
         }
      }
      
      // Calculate xz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i-1,j,k-1);
         botRght = perBGrid.get(i+1,j,k-1);
         topLeft = perBGrid.get(i-1,j,k+1);
         topRght = perBGrid.get(i+1,j,k+1);
         
         dPerB->at(fsgrids::dperb::dPERBydxz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBY) + topRght->at(fsgrids::bfield::PERBY) - botRght->at(fsgrids::bfield::PERBY) - topLeft->at(fsgrids::bfield::PERBY));
         
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 4);
         }
      }
      
      // Calculate yz mixed derivatives:
      if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
         
         botLeft = perBGrid.get(i,j-1,k-1);
         botRght = perBGrid.get(i,j+1,k-1);
         topLeft = perBGrid.get(i,j-1,k+1);
         topRght = perBGrid.get(i,j+1,k+1);
         
         dPerB->at(fsgrids::dperb::dPERBxdyz) = FOURTH * (botLeft->at(fsgrids::bfield::PERBX) + topRght->at(fsgrids::bfield::PERBX) - botRght->at(fsgrids::bfield::PERBX) - topLeft->at(fsgrids::bfield::PERBX));
         
      } else {
         if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 5);
         } else {
            sysBoundaries.getSysBoundary(sysBoundaryFlag)->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, 5);
         }
      }
   }
}


/*! \brief High-level derivative calculation wrapper function.
 * 

 * B has to be updated because after the system boundary update in propagateMagneticFieldSimple there is no consistent state of B yet everywhere.
 * 
 * Then the derivatives are calculated.
 * 
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param doMoments If true, the derivatives of moments (rho, V, P) are computed.
 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase,
   const bool& doMoments
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate face derivatives");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   
   switch (RKCase) {
      case RK_ORDER1:
         // Means initialising the solver as well as RK_ORDER1
         // standard case exchange PERB*,RHO,RHOV,P with neighbours
         // The update of PERB[XYZ] is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBGrid.updateGhostCells();
         if(doMoments) {
            momentsGrid.updateGhostCells();
         }
         break;
      case RK_ORDER2_STEP1:
         // Exchange PERB*_DT2,RHO_DT2,RHOV*_DT2,P*DT2 with neighbours The
         // update of PERB[XYZ]_DT2 is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBDt2Grid.updateGhostCells();
         if(doMoments) {
            momentsDt2Grid.updateGhostCells();
         }
         break;
      case RK_ORDER2_STEP2:
         // Exchange PERB*,RHO,RHOV*,P* with neighbours The update of B
         // is needed after the system boundary update of
         // propagateMagneticFieldSimple.
         perBGrid.updateGhostCells();
         if(doMoments) {
            momentsGrid.updateGhostCells();
         }
      break;
    default:
      cerr << __FILE__ << ":" << __LINE__ << " Went through switch, this should not happen." << endl;
      abort();
   }
   
   phiprof::stop(timer);

   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);

   // Calculate derivatives
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               calculateDerivatives(i,j,k, perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase, doMoments);
            } else {
               calculateDerivatives(i,j,k, perBDt2Grid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase, doMoments);
            }
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate face derivatives",N_cells,"Spatial Cells");   
}

/*! \brief Low-level spatial derivatives calculation.
 * 
 * For the cell with ID cellID calculate the spatial derivatives of BVOL or apply the derivative boundary conditions defined in project.h.
 * 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   cint i,
   cint j,
   cint k,
   SysBoundary& sysBoundaries
) {
   std::array<Real, fsgrids::volfields::N_VOL> * array = volGrid.get(i,j,k);
   
   std::array<Real, fsgrids::volfields::N_VOL> * left = NULL;
   std::array<Real, fsgrids::volfields::N_VOL> * rght = NULL;
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i-1,j,k);
      rght = volGrid.get(i+1,j,k);
      
      array->at(fsgrids::volfields::dPERBYVOLdx) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
      array->at(fsgrids::volfields::dPERBZVOLdx) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 0);
      }
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i,j-1,k);
      rght = volGrid.get(i,j+1,k);
      
      array->at(fsgrids::volfields::dPERBXVOLdy) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBZVOLdy) = limiter(left->at(fsgrids::volfields::PERBZVOL),array->at(fsgrids::volfields::PERBZVOL),rght->at(fsgrids::volfields::PERBZVOL));
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 1);
      }
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
      left = volGrid.get(i,j,k-1);
      rght = volGrid.get(i,j,k+1);
      
      array->at(fsgrids::volfields::dPERBXVOLdz) = limiter(left->at(fsgrids::volfields::PERBXVOL),array->at(fsgrids::volfields::PERBXVOL),rght->at(fsgrids::volfields::PERBXVOL));
      array->at(fsgrids::volfields::dPERBYVOLdz) = limiter(left->at(fsgrids::volfields::PERBYVOL),array->at(fsgrids::volfields::PERBYVOL),rght->at(fsgrids::volfields::PERBYVOL));
   } else {
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
      } else {
         sysBoundaries.getSysBoundary(technicalGrid.get(i,j,k)->sysBoundaryFlag)->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, 2);
      }
   }
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 * 
 * \param sysBoundaries System boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries
) {
   int timer;
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate volume derivatives");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   volGrid.updateGhostCells();
   
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   
   // Calculate derivatives
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
         }
      }
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}
