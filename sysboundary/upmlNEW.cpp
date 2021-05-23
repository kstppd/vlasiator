#include"upmlNEW.h"
#include "../readparameters.h"
#include "../parameters.h"
#include "../common.h"
#include "../spatial_cell.hpp"


ABC::UPML::UPML(FsGrid< std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,FsGrid< fsgrids::technical, 2> & technicalGrid,dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid){
   getParameters();
   classifyCells(pmlGrid,technicalGrid,mpiGrid);
   buildConductivity(pmlGrid,1.0);
   calculateParameters(pmlGrid,1.0);

   std::cout<< "PML built"<<std::endl;
}

bool ABC::UPML::resetPML(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt){
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if (myRank == MASTER_RANK){
      logFile<<" (PML)  Reseting Parameters..."<<"dt="<<dt<<std::endl;
   }
   buildConductivity(pmlGrid,dt);
   calculateParameters(pmlGrid,dt);
   return true;
}


bool ABC::UPML::getParameters(){

   typedef Parameters P;
   this->widthXM=P::pmlWidthXm;
   this->widthXP=P::pmlWidthXp;
   this->widthYM=P::pmlWidthYm;
   this->widthYP=P::pmlWidthYp;
   this->widthZM=P::pmlWidthZm;
   this->widthZP=P::pmlWidthZp;
   this->start=P::pmlStart;
   this->alpha=P::pmlAlpha;
   this-> logcells=P::pmlCells;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   if (myRank == MASTER_RANK){
         std::cout << "--------PML layers--------" << std::endl;
         std::cout << "start-\t" << this->start << std::endl;
         std::cout << "alpha\t" << this->alpha << std::endl;
         std::cout << "x-\t" << this->widthXM << std::endl;
         std::cout << "x+\t" << this->widthXP << std::endl;
         std::cout << "y-\t" << this->widthYM << std::endl;
         std::cout << "y+\t" << this->widthYP << std::endl;
         std::cout << "z-\t" << this->widthZM << std::endl;
         std::cout << "z+\t" << this->widthZP << std::endl;
         std::cout << "ACC " << this->logcells << std::endl;
         std::cout << "---------------------------" << std::endl;
   }
      return true;
   }

bool ABC::UPML::buildConductivity(FsGrid< std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt){

   std::array<Real, fsgrids::pml::N_PML> *val;
   const int *globalDims = &pmlGrid.getGlobalSize()[0];
   const int *localDims=&pmlGrid.getLocalSize()[0];
   std::array<int32_t, 3> pos;

   //Initialize all PML arrays to zero
   for (int k=0; k < localDims[2];k++){
      for (int j=0; j < localDims[1]; j++){
         for (int i=0; i< localDims[0]; i++){
               val=pmlGrid.get(i,j,k);
               for (int c=0; c< fsgrids::pml::N_PML; c++){
                  val->at(c) = 0.0;
               }
         }
      }
   }

   //Now let's fill in the conductivity tensors based on the cfg file
   bool isPmlCellXM,isPmlCellXP,isPmlCellYM,isPmlCellYP,isPmlCellZM,isPmlCellZP;
   int index;
   
   for (int k=0; k < localDims[2];k++){
      for (int j=0; j < localDims[1]; j++){
         for (int i=0; i< localDims[0]; i++){
            pos=pmlGrid.getGlobalIndices(i,j,k);
            isPmlCellXM=widthXM>0 && pos[0]>=start && pos[0]<widthXM + start;
            isPmlCellXP=widthXP>0 && pos[0]>globalDims[0]-start-widthXP -1&& pos[0]<=globalDims[0]-start-1; 
            isPmlCellYM=widthYM>0 && pos[1]>=start && pos[1]<widthYM + start;
            isPmlCellYP=widthYP>0 && pos[1]>globalDims[1]-start-widthYP -1&& pos[1]<=globalDims[1]-start-1; 
            isPmlCellZM=widthZM>0 && pos[2]>=start && pos[2]<widthZM + start;
            isPmlCellZP=widthZP>0 && pos[2]>globalDims[2]-start-widthZP -1&& pos[2]<=globalDims[2]-start-1; 

            val=pmlGrid.get(i,j,k);

            if(isPmlCellXM){
               index= widthXM-pos[0]+2*start;
               val=pmlGrid.get(index,j,k);
               val->at(fsgrids::pml::sigx) =1.0* (alpha*1.0/dt)*(pos[0]-start)/widthXM;
            }
            if(isPmlCellXP){
               index=pos[0]-(globalDims[0]-start-widthXP-1);
               val->at(fsgrids::pml::sigx) = (alpha*1.0/dt)*index/widthXP;
            }

            if(isPmlCellYM){
               index= widthYM-pos[1]+2*start;
               val=pmlGrid.get(i,index,k);
               val->at(fsgrids::pml::sigy) = (alpha*1.0/dt)*(pos[1]-start)/widthYM;
            }
            if(isPmlCellYP){
               index=pos[1]-(globalDims[1]-start-widthYP-1);
               val->at(fsgrids::pml::sigy) = (alpha*1.0/dt)*index/widthYP;
            }
            
            if(isPmlCellZM){
               index= widthZM-pos[2]+2*start;
               val=pmlGrid.get(i,j,index);
               val->at(fsgrids::pml::sigz) = (alpha*1.0/dt)*(pos[2]-start)/widthZM;
            }
            if(isPmlCellZP){
               index=pos[2]-(globalDims[2]-start-widthZP-1);
               val->at(fsgrids::pml::sigz) = (alpha*1.0/dt)*index/widthZP;
            }

         }
      }
   }
   pmlGrid.updateGhostCells();

   return true;
}


bool ABC::UPML::calculateParameters(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,Real dt){

   std::array<Real, fsgrids::pml::N_PML> *val;
   const int *globalDims = &pmlGrid.getGlobalSize()[0];
   const int *localDims=&pmlGrid.getLocalSize()[0];
   std::array<int32_t, 3> pos;
   Real Sx,Sy,Sz;
   Real mHx,mHy,mHz,mDx,mDy,mDz;

   Real c = 1.0;
   Real mu = 1.0;
   Real epsilon = 1.0;
   for (int k=0; k < localDims[2];k++){
      for (int j=0; j < localDims[1]; j++){
         for (int i=0; i< localDims[0]; i++){
            

            val=pmlGrid.get(i,j,k);
            Sx= val->at(fsgrids::pml::sigx);
            Sy= val->at(fsgrids::pml::sigy);
            Sz= val->at(fsgrids::pml::sigz);


            val->at(fsgrids::pml::mHx0) = (1.0/dt) + (Sy+Sz)/(2.0) +(Sy*Sz)*dt/(4.0);;
            val->at(fsgrids::pml::mHx1) = (1.0/val->at(fsgrids::pml::mHx0)) * ((1.0/dt) - (Sy+Sz)/(2.*epsilon) -(Sy*Sz)*dt/(4.*epsilon*epsilon));
            val->at(fsgrids::pml::mHx2) = (1.0/val->at(fsgrids::pml::mHx0)) * (c/mu);
            val->at(fsgrids::pml::mHx3) = (1.0/val->at(fsgrids::pml::mHx0)) * (c*dt*Sx /epsilon*mu);
            val->at(fsgrids::pml::mHx4) = (1.0/val->at(fsgrids::pml::mHx0)) * (dt*Sy*Sz /epsilon*epsilon);

            val->at(fsgrids::pml::mHy0) = (1.0/dt)   + (Sx+Sz)/(2*epsilon) +(Sx*Sz)*dt/(4*epsilon*epsilon);
            val->at(fsgrids::pml::mHy1) = (1.0/val->at(fsgrids::pml::mHy0)) * ((1.0/dt) - (Sx+Sz)/(2.0*epsilon) -(Sx*Sz)*dt/(4.0*epsilon*epsilon));
            val->at(fsgrids::pml::mHy2) = (1.0/val->at(fsgrids::pml::mHy0)) * (c/mu);
            val->at(fsgrids::pml::mHy3) = (1.0/val->at(fsgrids::pml::mHy0)) * (c*dt*Sy /epsilon*mu);
            val->at(fsgrids::pml::mHy4) = (1.0/val->at(fsgrids::pml::mHy0)) * (dt*Sx*Sz /epsilon*epsilon);

            val->at(fsgrids::pml::mHz0) =  (1.0/dt)   + (Sx+Sy)/(2.0*epsilon) +(Sx*Sy)*dt/(4.0*epsilon*epsilon);
            val->at(fsgrids::pml::mHz1) =  (1.0/val->at(fsgrids::pml::mHz0)) * ((1/dt) - (Sx+Sy)/(2.0*epsilon) -(Sx*Sy)*dt/(4.0*epsilon*epsilon));
            val->at(fsgrids::pml::mHz2) =  (1.0/val->at(fsgrids::pml::mHz0)) * (c/mu);
            val->at(fsgrids::pml::mHz3) =  (1.0/val->at(fsgrids::pml::mHz0)) * (c*dt*Sz /epsilon*mu);
            val->at(fsgrids::pml::mHz4) =  (1.0/val->at(fsgrids::pml::mHz0)) * (dt*Sx*Sy /epsilon*epsilon);

            val->at(fsgrids::pml::mDx0) =  (1.0/dt)   + (Sy+Sz)/(2.0*epsilon) +(Sy*Sz)*dt/(4.0*epsilon*epsilon);
            val->at(fsgrids::pml::mDx1) =  (1.0/val->at(fsgrids::pml::mDx0)) * ((1/dt) - (Sy+Sz)/(2.0*epsilon) -(Sy*Sz)*dt/(4.0*epsilon*epsilon));
            val->at(fsgrids::pml::mDx2) =  (1./dt)*(1.0/val->at(fsgrids::pml::mDx0)) * (c);
            val->at(fsgrids::pml::mDx3) =  (1.0/val->at(fsgrids::pml::mDx0)) * (c*dt*Sx /epsilon);
            val->at(fsgrids::pml::mDx4) =  (1.0/val->at(fsgrids::pml::mDx0)) * (dt*Sy*Sz /epsilon*epsilon);

            val->at(fsgrids::pml::mDy0) =  (1.0/dt)   + (Sx+Sz)/(2.0*epsilon) +(Sx*Sz)*dt/(4.0*epsilon*epsilon);
            val->at(fsgrids::pml::mDy1) =  (1.0/val->at(fsgrids::pml::mDy0)) * ((1/dt) - (Sx+Sz)/(2.0*epsilon) -(Sx*Sz)*dt/(4.0*epsilon*epsilon));
            val->at(fsgrids::pml::mDy2) =  (1./dt)*(1.0/val->at(fsgrids::pml::mDy0)) * (c);
            val->at(fsgrids::pml::mDy3) =  (1.0/val->at(fsgrids::pml::mDy0)) * (c*dt*Sy /epsilon);
            val->at(fsgrids::pml::mDy4) =  (1.0/val->at(fsgrids::pml::mDy0)) * (dt*Sx*Sz /epsilon*epsilon) ;
 
            val->at(fsgrids::pml::mDz0) =  (1.0/dt)   + (Sx+Sy)/(2.0*epsilon) +(Sx*Sy)*dt/(4.0*epsilon*epsilon);
            val->at(fsgrids::pml::mDz1) =  (1.0/val->at(fsgrids::pml::mDz0)) * ((1/dt) - (Sx+Sy)/(2.0*epsilon) -(Sx*Sy)*dt/(4.0*epsilon*epsilon));
            val->at(fsgrids::pml::mDz2) =  ((1./dt)*1.0/val->at(fsgrids::pml::mDz0)) * (c);
            val->at(fsgrids::pml::mDz3) =  (1.0/val->at(fsgrids::pml::mDz0)) * (c*dt*Sz /epsilon);
            val->at(fsgrids::pml::mDz4) =  (1.0/val->at(fsgrids::pml::mDz0)) * (dt*Sx*Sy /epsilon*epsilon);
             

         }
      }
   }

   pmlGrid.updateGhostCells();

   return true;
}



bool ABC::UPML::classifyCells(FsGrid <std::array<Real, fsgrids::pml::N_PML>, 2> &pmlGrid,
                               FsGrid< fsgrids::technical, 2> & technicalGrid,
                               dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid)
{


   const int *globalDims = &pmlGrid.getGlobalSize()[0];
   const int *localDims=&pmlGrid.getLocalSize()[0];
   std::array<int32_t, 3> pos;
   std:: vector<CellID> cells = mpiGrid.get_cells();

   //Now let's fill in the conductivity tensors based on the cfg file
   bool isPmlCellXM,isPmlCellXP,isPmlCellYM,isPmlCellYP,isPmlCellZM,isPmlCellZP;
   bool doAssign;
   
   //Assign FsGrid cells
   for (int k=0; k < localDims[2];k++){
      for (int j=0; j < localDims[1]; j++){
         for (int i=0; i< localDims[0]; i++){
            pos=pmlGrid.getGlobalIndices(i,j,k);
            isPmlCellXM=widthXM>0 && pos[0]>=start && pos[0]<widthXM + start;
            isPmlCellXP=widthXP>0 && pos[0]>globalDims[0]-start-widthXP -1&& pos[0]<=globalDims[0]-start-1; 
            isPmlCellYM=widthYM>0 && pos[1]>=start && pos[1]<widthYM + start;
            isPmlCellYP=widthYP>0 && pos[1]>globalDims[1]-start-widthYP -1&& pos[1]<=globalDims[1]-start-1; 
            isPmlCellZM=widthZM>0 && pos[2]>=start && pos[2]<widthZM + start;
            isPmlCellZP=widthZP>0 && pos[2]>globalDims[2]-start-widthZP -1&& pos[2]<=globalDims[2]-start-1; 

            
            doAssign = isPmlCellXM || isPmlCellXP || isPmlCellYM || isPmlCellYP || isPmlCellZM || isPmlCellZP;
            technicalGrid.get(i, j, k)->pmlFlag = 0;
                // printf("-> %i %i %i %i %i %i \n ", isPmlCellXP, isPmlCellXM, isPmlCellYP, isPmlCellYM, isPmlCellZP, isPmlCellZM);

                if (doAssign)
            {
               // std::cerr<< "I am a PML FSGRID cell at "<< i<<" "<<j<<" "<<k<<std::endl;
               technicalGrid.get(i,j,k)->pmlFlag=1;
               // technicalGrid.get(i, j, k)->sysBoundaryFlag = sysboundarytype::PMLCELL;
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
      #warning "TODO: Fix casting below"
      std::array<int32_t, 3> pos{ (x-Parameters::xmin)/dx ,(y-Parameters::ymin)/dy ,(z-Parameters::zmin)/dz  };
      

            isPmlCellXM=widthXM>0 && pos[0]>=start && pos[0]<widthXM + start;
            isPmlCellXP=widthXP>0 && pos[0]>globalDims[0]-start-widthXP -1&& pos[0]<=globalDims[0]-start-1; 
            isPmlCellYM=widthYM>0 && pos[1]>=start && pos[1]<widthYM + start;
            isPmlCellYP=widthYP>0 && pos[1]>globalDims[1]-start-widthYP -1&& pos[1]<=globalDims[1]-start-1; 
            isPmlCellZM=widthZM>0 && pos[2]>=start && pos[2]<widthZM + start;
            isPmlCellZP=widthZP>0 && pos[2]>globalDims[2]-start-widthZP -1&& pos[2]<=globalDims[2]-start-1; 


      // printf("-> %i %i %i %i %i %i \n " , isPmlCellXP, isPmlCellXM,  isPmlCellYP, isPmlCellYM, isPmlCellZP, isPmlCellZM);
      doAssign = isPmlCellXM || isPmlCellXP || isPmlCellYM || isPmlCellYP || isPmlCellZM || isPmlCellZP;
      mpiGrid[cells[i]]->pmlFlag = 0;

      if(doAssign && this->logcells) {
         // std::cerr << "I am a PML MPIGRID cell at " << x << " " << y << " " << z << std::endl;
         mpiGrid[cells[i]]->pmlFlag = 1 ;
      }
   }

   pmlGrid.updateGhostCells();
   return true;
}