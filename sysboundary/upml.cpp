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
/*!\file upml.cpp
 * \builds the arrays needed for the Uniaxial PML.
 */

#include "upml.h"
#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <cmath>
#include <omp.h>
#include "fsgrid.hpp"
#include "../definitions.h"
#include "../common.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"



bool assingFSGridCells(FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,
                       FsGrid< fsgrids::technical, 2> & technicalGrid,
                       dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                       int width,int  side, int comp){

   std::array<Real, fsgrids::pml::N_PML> *pmlValue;
   const int *globalDims = &pmlGrid.getGlobalSize()[0];
   const int *pmlDims=&pmlGrid.getLocalSize()[0];
   std::array<int32_t, 3> pos;
   Real xnum, xd;
   Real xxn, xn;
   Real alpha = P::pmlAlpha;
   typedef Parameters P;
   int start = P::pmlStart ;
   int enumOrder =5*comp;
   std:: vector<CellID> cells = mpiGrid.get_cells();
   bool doAssign;
   
   switch(side){
      case (-1):

         for (int kk = 0; kk < pmlDims[2] ; kk++){
            for (int jj = 0; jj < pmlDims[1]; jj++){
               for (int ii = 0; ii <pmlDims[0]; ii++){

                  // Get Global Arrays Indeces
                  pos=pmlGrid.getGlobalIndices(ii,jj,kk);
                  doAssign = width>0 &&  pos[comp]>=start && pos[comp]<width+start;
                  if (doAssign){
                     
                     // Get Local  Arrays
                     pmlValue = pmlGrid.get(ii, jj, kk);

                     xnum =width-pos[comp] +start;
                     xd = width;
                     xxn =xnum/xd;
                     xn =alpha*(xxn*xxn*xxn);
                     pmlValue->at( fsgrids::pml(enumOrder) )=1/(1+xn);
                     pmlValue->at(fsgrids::pml(enumOrder+1))=(1-xn)/(1+xn);

                     //printf("comp=-%d pos=%d ---->xn=%f, xxn=%f \n" ,comp, pos[comp],xn,xxn  );

                     xxn=(xnum-0.5)/xd;
                     xn=alpha*(xxn*xxn*xxn);
                     pmlValue->at(fsgrids::pml(enumOrder+2))=xn;
                     pmlValue->at(fsgrids::pml(enumOrder+3))=1/(1+xn);
                     pmlValue->at(fsgrids::pml(enumOrder+4))=(1-xn)/(1+xn);

                     //Set the technical grid flag
                     technicalGrid.get(ii,jj,kk)->pmlFlag= sysboundarytype::PMLCELL; 
                     
                     
                  }
               }
            }
         }
         //Assing DCCRG Cells
         for(uint i=0; i<cells.size(); i++) {
         
            creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
            creal dx = cellParams[CellParams::DX];
            creal dy = cellParams[CellParams::DY];
            creal dz = cellParams[CellParams::DZ];
            creal x = cellParams[CellParams::XCRD] + 0.5*dx;
            creal y = cellParams[CellParams::YCRD] + 0.5*dy;
            creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
            std::array<creal,3> p={x,y,z};

            doAssign = width>0 &&  p[comp]>start && p[comp]<width+start;
            if(doAssign) {
               mpiGrid[cells[i]]->pmlFlag = sysboundarytype::PMLCELL ;
            }
         }
         break;
      
      case(1):
        
         for (int kk = 0; kk < pmlDims[2] ; kk++){
            for (int jj = 0; jj < pmlDims[1]; jj++){
               for (int ii = 0; ii <pmlDims[0]; ii++){

                  // Get Global Arrays Indeces
                  pos=pmlGrid.getGlobalIndices(ii,jj,kk);
                  doAssign = (width>0 && pos[comp]>globalDims[0]-start- width && pos[comp] <=globalDims[0]-start);
                  if (doAssign){
                     
                     // Get Local  Arrays
                     pmlValue = pmlGrid.get(ii, jj, kk);

                     xnum =start+ width-(globalDims[0]- pos[comp]);
                     xd = width;
                     xxn =xnum/xd;
                     xn =alpha*(xxn*xxn*xxn);
                     pmlValue->at(fsgrids::pml(enumOrder))=1/(1+xn);
                     pmlValue->at(fsgrids::pml(enumOrder+1))=(1-xn)/(1+xn);

                     //printf("comp=+%d pos=%d ---->xn=%f, xxn=%f \n" ,comp, pos[comp],xn,xxn  );

                     xxn=(xnum-0.5)/xd;
                     xn=alpha*(xxn*xxn*xxn);
                     pmlValue->at(fsgrids::pml(enumOrder+2))=xn;
                     pmlValue->at(fsgrids::pml(enumOrder+3))=1/(1+xn);
                     pmlValue->at(fsgrids::pml(enumOrder+4))=(1-xn)/(1+xn);

                     //Set the technical grid flag
                     technicalGrid.get(ii,jj,kk)->pmlFlag = sysboundarytype::PMLCELL; 
                     
                     
                  }
               }
            }
         }
         
         //Assing DCCRG Cells
         for(uint i=0; i<cells.size(); i++) {
            
            creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
            creal dx = cellParams[CellParams::DX];
            creal dy = cellParams[CellParams::DY];
            creal dz = cellParams[CellParams::DZ];
            creal x = cellParams[CellParams::XCRD] + 0.5*dx;
            creal y = cellParams[CellParams::YCRD] + 0.5*dy;
            creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
            std::array<creal,3> p={x,y,z};
           
            doAssign =  (width>0 && p[comp]>globalDims[0]-start- width && p[comp] <=globalDims[0]-start); 
            if(doAssign) {
               mpiGrid[cells[i]]->pmlFlag = sysboundarytype::PMLCELL ;
            }
         }

         break;
   }

   return true;
}


bool buildPMLGrid(FsGrid<std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,FsGrid< fsgrids::technical, 2> & technicalGrid,dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid){

   typedef Parameters P;
   
   /*----Get PML Array and Domain Size-----*/
   std::array<Real, fsgrids::pml::N_PML> *pmlValue;
   const int *pmlDims=&pmlGrid.getLocalSize()[0];

    /*-----Initially set all arrays to one-----*/
    /*Iterate over domain and set the PML arrays to 1.0 thus not affecting the fieldsolver*/ 
   
    for (int kk = 0; kk < pmlDims[2]; kk++){
        for (int jj = 0; jj < pmlDims[1]; jj++){
            for (int ii = 0; ii < pmlDims[0]; ii++){
            
               pmlValue=pmlGrid.get(ii,jj,kk);
               pmlValue->at(fsgrids::pml::PFI1) = 1.0;
               pmlValue->at(fsgrids::pml::PFI2) = 1.0;
               pmlValue->at(fsgrids::pml::PFI3) = 1.0;
               pmlValue->at(fsgrids::pml::PGI2) = 1.0;
               pmlValue->at(fsgrids::pml::PGI3) = 1.0;
                
               pmlValue->at(fsgrids::pml::PFJ1) = 1.0;
               pmlValue->at(fsgrids::pml::PFJ2) = 1.0;
               pmlValue->at(fsgrids::pml::PFJ3) = 1.0;
               pmlValue->at(fsgrids::pml::PGJ2) = 1.0;
               pmlValue->at(fsgrids::pml::PGJ3) = 1.0;
                   
               pmlValue->at(fsgrids::pml::PFK1) = 1.0;
               pmlValue->at(fsgrids::pml::PFK2) = 1.0;
               pmlValue->at(fsgrids::pml::PFK3) = 1.0;
               pmlValue->at(fsgrids::pml::PGK2) = 1.0;
               pmlValue->at(fsgrids::pml::PGK3) = 1.0;
            }
        }
    }



   printf(" pml widths = %d %d %d %d %d %d\n",P::pmlWidthXm,P::pmlWidthXp,P::pmlWidthYm,P::pmlWidthYp,P::pmlWidthZm,P::pmlWidthZp);
   printf(" pml start = %d ,alpha= %f \n",P::pmlStart,P::pmlAlpha);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthXm,-1,0);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthYm,-1,1);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthZm,-1,2);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthXp, 1,0);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthYp, 1,1);
   assingFSGridCells(pmlGrid,technicalGrid,mpiGrid,P::pmlWidthZp, 1,2);
   pmlGrid.updateGhostCells();

   return true;
}
