#include "upml.h"
#include "../parameters.h"
#include "../readparameters.h"
#include <cmath>

namespace PML=PerfectlyMatchedLayer;




bool PML::madFilter(FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
      FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid){


   const int stencilWidth = 3;          //  1D kernel stencil
   const int kernelOffset = 1;          // offset of 3 point stencil 1D kernel
   const Real kernelWeight= 1.0/27.0;   // 3D boxcar kernel weight
   const int passes=9; 
   FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> swapGrid = EGrid;  
   
   const int* LocalSize=&EGrid.getLocalSize()[0];

   for (int pass=0; pass<passes; pass++){
      for (int k=0; k<LocalSize[2]; k++){
         for (int j=0; j<LocalSize[1]; j++){
            for (int i=0; i< LocalSize[0];i++){

               bool isPML = technicalGrid.get(i,j,k)->pmlCell==UPMLCELLS::UPMLCELL;
               if (!isPML ||
                   technicalGrid.get(i, j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
                   (technicalGrid.get(i, j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i, j,k)->sysBoundaryLayer == 2)
                     ){
                  continue;
               }

               std::array<Real, fsgrids::efield::N_EFIELD> *cell;  
               std::array<Real,fsgrids::efield::N_EFIELD> *swap;
              
               // Set Cell to zero before passing filter
               swap = swapGrid.get(i,j,k);
               for (int e = 0; e < fsgrids::efield::N_EFIELD; ++e) {
                 swap->at(e)=0.0;
               }
               
               for (int a=0; a<stencilWidth; a++){
                 for (int b=0; b<stencilWidth; b++){
                   for (int c=0; c<stencilWidth; c++){

                     int xn=i+a-kernelOffset;
                     int yn=j+b-kernelOffset;
                     int zn=k+c-kernelOffset;

                     cell = EGrid.get(xn, yn, zn);
                     swap = swapGrid.get(i,j,k);

                     for (int e = 0; e < fsgrids::efield::N_EFIELD; ++e) {
                       swap->at(e)+=cell->at(e) * kernelWeight;

                    } 
                  }
                }
              }


            }
         }
      }

      //Swap buffers
      EGrid=swapGrid;
   }


   return true;
}




PML::UPML::UPML(){
    this->addParameters();
     this->getParameters();
 }




void PML::UPML::addParameters(){
      Readparameters::add("UPML.width", "Uniaxial Pefectly Matched Layer's width", 0);    
       Readparameters::add("UPML.xp", "UPML xdim +", 0);
       Readparameters::add("UPML.yp", "UPML ydim +", 0);
       Readparameters::add("UPML.zp", "UPML zdim +", 0);
       Readparameters::add("UPML.xm", "UPML xdim -", 0);
       Readparameters::add("UPML.ym", "UPML ydim -", 0);
       Readparameters::add("UPML.zm", "UPML zdim -", 0);
}


void PML::UPML::getParameters(){

   using namespace std;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   this->upmlWidth = Parameters::upmlCells;
   this->pml_Xp  =  Parameters::upmlXp;
   this->pml_Xm  =  Parameters::upmlXm;
   this->pml_Yp  =  Parameters::upmlYp;
   this->pml_Ym  =  Parameters::upmlYm;
   this->pml_Zp  =  Parameters::upmlZp;
   this->pml_Zm  =  Parameters::upmlZm;



   if (myRank == MASTER_RANK){

      std::cout<< "***** UPML Summary **********\n";
      std::cout<<"upml Width= "<< this->upmlWidth<<"\n";;
      std::cout<<"upml x+= "<< this->pml_Xp<<"\n";
      std::cout<<"upml x-= "<< this->pml_Xm<<"\n";
      std::cout<<"upml y+= "<< this->pml_Yp<<"\n";
      std::cout<<"upml y-= "<< this->pml_Ym<<"\n";
      std::cout<<"upml z+= "<< this->pml_Zp<<"\n";
      std::cout<<"upml z-= "<< this->pml_Zm<<"\n";
      std::cout<< "*******************"<<std::endl;
   }
}


void PML::UPML::update(FsGrid<std::array<Real, fsgrids::upml::N_UPML>, FS_STENCIL_WIDTH>& fsUpml,
               FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                const Real dt)
{


   // Real muz = physicalconstants::MU_0;
   // Real epsz = physicalconstants::EPS_0;
   Real muz = 1.0;
   Real epsz = 1.0;
   Real etaz = sqrt(muz / epsz);
   Real mur = 1.0;
   Real epsr = 1.0;
   Real eta = etaz * sqrt(mur / epsr);
   // this->upmlWidth = 10;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   


   const int* globalDims = &fsUpml.getGlobalSize()[0];
   const int* localDims = &fsUpml.getLocalSize()[0];

   // Zero Out PML 
   for (int k = 0; k < localDims[2]; k++) {
      for (int j = 0; j < localDims[1]; j++) {
         for (int i = 0; i < localDims[0]; i++) {
            for (int e = 0; e < fsgrids::upml::N_UPML; ++e) {
               std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(e) = 0.0;
            }
         }
      }
   }

   // Zero Out PML
   for (int k = 0; k < localDims[2]; k++) {
      for (int j = 0; j < localDims[1]; j++) {
         for (int i = 0; i < localDims[0]; i++) {
            technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::NORMALCELL;
         }
      }
   }

   //Set some constants for inner part of domain
   Real C1 = 1.0;
   Real C2 = dt;
   Real C3 = 1.0;
   Real C4 = 1.0 / 2.0 / epsr / epsr / epsz / epsz;
   Real C5 = 2.0 * epsr * epsz;
   Real C6 = 2.0 * epsr * epsz;

   Real D1 = 1.0;
   Real D2 = dt;
   Real D3 = 1.0;
   Real D4 = 1.0 / 2.0 / epsr / epsz / mur / muz;
   Real D5 = 2.0 * epsr * epsz;
   Real D6 = 2.0 * epsr * epsz;

   for (int k = 0; k < localDims[2]; k++) {
      for (int j = 0; j < localDims[1]; j++) {
         for (int i = 0; i < localDims[0]; i++) {

            std::array<Real, fsgrids::upml::N_UPML>* val;
            val = fsUpml.get(i, j, k);

            val->at(fsgrids::upml::C1EX) = C1;
            val->at(fsgrids::upml::C1EY) = C1;
            val->at(fsgrids::upml::C1EZ) = C1;
            val->at(fsgrids::upml::C1BX) = D1;
            val->at(fsgrids::upml::C1BY) = D1;
            val->at(fsgrids::upml::C1BZ) = D1;

            val->at(fsgrids::upml::C2EX) = C2;
            val->at(fsgrids::upml::C2EY) = C2;
            val->at(fsgrids::upml::C2EZ) = C2;
            val->at(fsgrids::upml::C2BX) = D2;
            val->at(fsgrids::upml::C2BY) = D2;
            val->at(fsgrids::upml::C2BZ) = D2;

            val->at(fsgrids::upml::C3EX) = C3;
            val->at(fsgrids::upml::C3EY) = C3;
            val->at(fsgrids::upml::C3EZ) = C3;
            val->at(fsgrids::upml::C3BX) = D3;
            val->at(fsgrids::upml::C3BY) = D3;
            val->at(fsgrids::upml::C3BZ) = D3;

            val->at(fsgrids::upml::C4EX) = C4;
            val->at(fsgrids::upml::C4EY) = C4;
            val->at(fsgrids::upml::C4EZ) = C4;
            val->at(fsgrids::upml::C4BX) = D4;
            val->at(fsgrids::upml::C4BY) = D4;
            val->at(fsgrids::upml::C4BZ) = D4;

            val->at(fsgrids::upml::C5EX) = C5;
            val->at(fsgrids::upml::C5EY) = C5;
            val->at(fsgrids::upml::C5EZ) = C5;
            val->at(fsgrids::upml::C5BX) = D5;
            val->at(fsgrids::upml::C5BY) = D5;
            val->at(fsgrids::upml::C5BZ) = D5;

            val->at(fsgrids::upml::C6EX) = C6;
            val->at(fsgrids::upml::C6EY) = C6;
            val->at(fsgrids::upml::C6EZ) = C6;
            val->at(fsgrids::upml::C6BX) = D6;
            val->at(fsgrids::upml::C6BY) = D6;
            val->at(fsgrids::upml::C6BZ) = D6;
         }
      }
   }

   
   Real ds= Parameters::dx_ini;
   int offset = Parameters::upmlOffset;
   Real orderbc=5.0;
   Real rmax=1e-16;
   Real delbc = this->upmlWidth * ds;
   Real sigmam = -log(rmax) * (orderbc + 1.0) / (2.0 * eta * delbc);
   sigmam*=Parameters::upmlFactor ;
   Real sigfactor = sigmam / (ds * pow(delbc,orderbc) * (orderbc + 1.0));
   Real kmax = 1.0;
   Real kfactor = (kmax - 1.0) / ds / (orderbc + 1.0) / pow(delbc, orderbc);


 //X- component 
 if (this->pml_Xm){

      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[0]>=this->upmlWidth +offset || pos[0] < offset  ) continue;

               Real x1 = (this->upmlWidth - pos[0] + 1) * ds;
               Real x2 = (this->upmlWidth - pos[0]) * ds;

               Real ksi = (this->upmlWidth - pos[0]+offset) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));

               // if (pos[1]==15 && pos[2]==15 ) std::cout<<pos[0]<<" "<<ksi<<" "<<sigmam<<" "<<sigma<<std::endl;
               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EX) = facp;
               val->at(fsgrids::upml::C6EX) = facm;
               val->at(fsgrids::upml::C1BZ) = facm/facp;
               val->at(fsgrids::upml::C2BZ) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BY) = facm / facp;
               val->at(fsgrids::upml::C4BY) = 1.0 / facp / mur / muz;
         
               //----

               x1 = (this->upmlWidth - pos[0] + 1.5) * ds;
               x2 = (this->upmlWidth - pos[0] + 0.5) * ds;

               ksi = (this->upmlWidth - pos[0] + .5 +offset) / this->upmlWidth;
               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);
      
               // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);
     
               // std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EZ) = facm/facp;
               val->at(fsgrids::upml::C2EZ) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EY) = facm / facp;
               val->at(fsgrids::upml::C4EY) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BX) = facp;
               val->at(fsgrids::upml::C6BX) = facm;


               // Assign cells
               technicalGrid.get(i,j,k)->pmlCell=UPMLCELLS::UPMLCELL;
            }
         }
      }
      


      //X PEC wall
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[0] != offset  ) continue;

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i,j,k);
               val->at(fsgrids::upml::C1EZ) = -1.0;
               val->at(fsgrids::upml::C2EZ) = -0.0;
               val->at(fsgrids::upml::C3EY) = -1.0;
               val->at(fsgrids::upml::C4EY) = -0.0;
            }
         }
      }
 }


   
   //Y- component
   if (this->pml_Ym){
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){
         

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[1]>=this->upmlWidth + offset|| pos[1] < offset) continue;

               Real y1 = (this->upmlWidth - pos[1] + 1) * ds;
               Real y2 = (this->upmlWidth - pos[1]) * ds;

               Real ksi = (this->upmlWidth - pos[1] + offset) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EY) = facp;
               val->at(fsgrids::upml::C6EY) = facm;
               val->at(fsgrids::upml::C1BX) = facm/facp;
               val->at(fsgrids::upml::C2BX) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BZ) = facm / facp;
               val->at(fsgrids::upml::C4BZ) = 1.0 / facp / mur / muz;
         
         
               y1 = (this->upmlWidth - pos[1] + 1.5) * ds;
               y2 = (this->upmlWidth - pos[1] + 0.5) * ds;
               ksi = (this->upmlWidth - pos[1] + .5 +offset) / this->upmlWidth;
               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);

               // std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EX) = facm/facp;
               val->at(fsgrids::upml::C2EX) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EZ) = facm / facp;
               val->at(fsgrids::upml::C4EZ) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BY) = facp;
               val->at(fsgrids::upml::C6BY) = facm;
               technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::UPMLCELL;
            }
         }
      }
      




   //Y PEC wal/*l*/
      /*for (int i=0; i < localDims[0];i++){*/
         /*for (int k=0; k < localDims[2];k++){*/
            /*for (int j=0; j < localDims[1]; j++){*/

               /*std::array<int32_t, 3> pos;*/
               /*pos = fsUpml.getGlobalIndices(i, j, k);*/

               /*//Skip cells that do not belong to the X- PML layers*/
               /*if (pos[1] != offset  ) continue;*/

               /*std::array<Real, fsgrids::upml::N_UPML> *val;*/
               /*val = fsUpml.get(i,j,k);*/
               /*val->at(fsgrids::upml::C1EX) = -1.0;*/
               /*val->at(fsgrids::upml::C2EX) = 0.0;*/
               /*val->at(fsgrids::upml::C3EZ) = -1.0;*/
               /*val->at(fsgrids::upml::C4EZ) = 0.0;*/
            /*}*/
         /*}*/
      /*}*/
   }



  //Z- component
   if (this->pml_Zm){

      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){
         

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[2]>=this->upmlWidth + offset|| pos[2] < offset) continue;

               Real z1 = (this->upmlWidth - pos[2] + 1) * ds;
               Real z2 = (this->upmlWidth - pos[2]) * ds;

               Real ksi = (this->upmlWidth - pos[2] + offset) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EZ) = facp;
               val->at(fsgrids::upml::C6EZ) = facm;
               val->at(fsgrids::upml::C1BY) = facm/facp;
               val->at(fsgrids::upml::C2BY) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BX) = facm / facp;
               val->at(fsgrids::upml::C4BX) = 1.0 / facp / mur / muz;
         
         
               z1 = (this->upmlWidth - pos[2] + 1.5) * ds;
               z2 = (this->upmlWidth - pos[2] + 0.5) * ds;
               ksi = (this->upmlWidth - pos[2] + .5 +offset) / this->upmlWidth;
               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);

               // std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EY) = facm/facp;
               val->at(fsgrids::upml::C2EY) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EX) = facm / facp;
               val->at(fsgrids::upml::C4EX) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BZ) = facp;
               val->at(fsgrids::upml::C6BZ) = facm;
               technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::UPMLCELL;
            }
         }
      }
      




   //Z- PEC wall
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[2] != offset  ) continue;

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i,j,k);
               val->at(fsgrids::upml::C1EY) = -1.0;
               val->at(fsgrids::upml::C2EY) = 0.0;
               val->at(fsgrids::upml::C3EX) = -1.0;
               val->at(fsgrids::upml::C4EX) = 0.0;
            }
         }
      }
   }








 //X+ component 
 if (this->pml_Xp){
      int Nx=Parameters::xcells_ini;
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[0]>Nx-offset-1 || pos[0] < Nx-this->upmlWidth-offset  ) continue;

               Real x1 = (this->upmlWidth - pos[0] + 1) * ds;
               Real x2 = (this->upmlWidth - pos[0]) * ds;

               Real ksi = (this->upmlWidth - (Nx-pos[0]-offset)+1) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);
               // if (pos[1]==15 && pos[2]==15 ) std::cout<<pos[0]<<" "<<ksi<<" "<<sigmam<<" "<<sigma<<std::endl;

               // Real ksi = (this->upmlWidth - pos[0] + offset) / this->upmlWidth;

               // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));

               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EX) = facp;
               val->at(fsgrids::upml::C6EX) = facm;
               val->at(fsgrids::upml::C1BZ) = facm/facp;
               val->at(fsgrids::upml::C2BZ) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BY) = facm / facp;
               val->at(fsgrids::upml::C4BY) = 1.0 / facp / mur / muz;
         
               //----

               x1 = (this->upmlWidth - pos[0] + 1.5) * ds;
               x2 = (this->upmlWidth - pos[0] + 0.5) * ds;

               ksi = (this->upmlWidth - (Nx - pos[0] + offset +0.5)) / this->upmlWidth;

               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);
      
               // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);
     
               // std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EZ) = facm/facp;
               val->at(fsgrids::upml::C2EZ) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EY) = facm / facp;
               val->at(fsgrids::upml::C4EY) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BX) = facp;
               val->at(fsgrids::upml::C6BX) = facm;
               technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::UPMLCELL;
            }
         }
      }
      


      //X PEC wall
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[0] != Nx- offset  ) continue;

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i,j,k);
               val->at(fsgrids::upml::C1EZ) = -1.0;
               val->at(fsgrids::upml::C2EZ) = -0.0;
               val->at(fsgrids::upml::C3EY) = -1.0;
               val->at(fsgrids::upml::C4EY) = -0.0;
            }
         }
      }
 }



   //Y+ component
   if (this->pml_Yp){
      int Ny =Parameters::ycells_ini;
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){
         

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[1]>Ny-offset -1 || pos[1] < Ny-this->upmlWidth-offset  ) continue;


               Real y1 = (this->upmlWidth - pos[1] + 1) * ds;
               Real y2 = (this->upmlWidth - pos[1]) * ds;

               Real ksi = (this->upmlWidth - (Ny - pos[1] - offset)+1) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);
               // if (pos[0]==15 && pos[2]==15 ) std::cout<<pos[1]<<" "<<ksi<<" "<<sigmam<<" "<<sigma<<std::endl;

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EY) = facp;
               val->at(fsgrids::upml::C6EY) = facm;
               val->at(fsgrids::upml::C1BX) = facm/facp;
               val->at(fsgrids::upml::C2BX) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BZ) = facm / facp;
               val->at(fsgrids::upml::C4BZ) = 1.0 / facp / mur / muz;
         
         
               y1 = (this->upmlWidth - pos[1] + 1.5) * ds;
               y2 = (this->upmlWidth - pos[1] + 0.5) * ds;
               ksi = (this->upmlWidth - (Ny - pos[1] + offset+.5)) / this->upmlWidth;
               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);

               // std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EX) = facm/facp;
               val->at(fsgrids::upml::C2EX) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EZ) = facm / facp;
               val->at(fsgrids::upml::C4EZ) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BY) = facp;
               val->at(fsgrids::upml::C6BY) = facm;
               technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::UPMLCELL;
            }
         }
      }
      




   ////Y+ PEC wall
      //for (int i=0; i < localDims[0];i++){
         //for (int k=0; k < localDims[2];k++){
            //for (int j=0; j < localDims[1]; j++){

               //std::array<int32_t, 3> pos;
               //pos = fsUpml.getGlobalIndices(i, j, k);

               ////Skip cells that do not belong to the X- PML layers
               //if (pos[1] != Ny-offset  ) continue;

               //std::array<Real, fsgrids::upml::N_UPML> *val;
               //val = fsUpml.get(i,j,k);
               //val->at(fsgrids::upml::C1EX) = -1.0;
               //val->at(fsgrids::upml::C2EX) = 0.0;
               //val->at(fsgrids::upml::C3EZ) = -1.0;
               //val->at(fsgrids::upml::C4EZ) = 0.0;
            //}
         //}
      //}
   }

  //Z+ component
   if (this->pml_Zp){
      int Nz=Parameters::zcells_ini;
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){
         

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[2]>Nz-offset-1 || pos[2] < Nz-this->upmlWidth-offset  ) continue;


               Real z1 = (this->upmlWidth - pos[2] + 1) * ds;
               Real z2 = (this->upmlWidth - pos[2]) * ds;

               Real ksi = (this->upmlWidth - (Nz - pos[2] - offset)+1) / this->upmlWidth;
               Real sigma = sigmam * pow(ksi, orderbc);
               Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               Real facm = (2 * epsr * epsz * ki - sigma * dt);
               Real facp = (2 * epsr * epsz * ki + sigma * dt);

               std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C5EZ) = facp;
               val->at(fsgrids::upml::C6EZ) = facm;
               val->at(fsgrids::upml::C1BY) = facm/facp;
               val->at(fsgrids::upml::C2BY) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3BX) = facm / facp;
               val->at(fsgrids::upml::C4BX) = 1.0 / facp / mur / muz;
         
         
               z1 = (this->upmlWidth - pos[2] + 1.5) * ds;
               z2 = (this->upmlWidth - pos[2] + 0.5) * ds;
               ksi = (this->upmlWidth - (Nz - pos[2] + offset+.5)) / this->upmlWidth;
               sigma = sigmam * pow(ksi, orderbc);
               ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

               // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
               facm = (2 * epsr * epsz * ki - sigma * dt);
               facp = (2 * epsr * epsz * ki + sigma * dt);

               // std::array<Real, fsgrids::upml::N_UPML>* val;
               val = fsUpml.get(i, j, k);
               val->at(fsgrids::upml::C1EY) = facm/facp;
               val->at(fsgrids::upml::C2EY) = 2.0 * epsr * epsz * dt / facp;
               val->at(fsgrids::upml::C3EX) = facm / facp;
               val->at(fsgrids::upml::C4EX) = 1.0 / facp / epsr / epsz;
               val->at(fsgrids::upml::C5BZ) = facp;
               val->at(fsgrids::upml::C6BZ) = facm;
               technicalGrid.get(i, j, k)->pmlCell = UPMLCELLS::UPMLCELL;
            }
         }
      }
      




   //Z+ PEC wall
      for (int i=0; i < localDims[0];i++){
         for (int k=0; k < localDims[2];k++){
            for (int j=0; j < localDims[1]; j++){

               std::array<int32_t, 3> pos;
               pos = fsUpml.getGlobalIndices(i, j, k);

               //Skip cells that do not belong to the X- PML layers
               if (pos[2] !=Nz- offset  ) continue;

               std::array<Real, fsgrids::upml::N_UPML> *val;
               val = fsUpml.get(i,j,k);
               val->at(fsgrids::upml::C1EY) = -1.0;
               val->at(fsgrids::upml::C2EY) = 0.0;
               val->at(fsgrids::upml::C3EX) = -1.0;
               val->at(fsgrids::upml::C4EX) = 0.0;
            }
         }
      }
   }




}
