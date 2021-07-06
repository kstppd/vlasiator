#include "upml.h"
#include "../parameters.h"
#include "../readparameters.h"
#include <cmath>

namespace PML=PerfectlyMatchedLayer;


void PML::UPML::addParameters(){
      Readparameters::add("UPML.width", "Uniaxial Pefectly Matched Layer's width", 4);
      Readparameters::add("UPML.xp", "UPML xdim +", 0);
      Readparameters::add("UPML.yp", "UPML ydim +", 0);
      Readparameters::add("UPML.zp", "UPML zdim +", 0);
      Readparameters::add("UPML.xm", "UPML xdim -", 0);
      Readparameters::add("UPML.ym", "UPML ydim -", 0);
      Readparameters::add("UPML.zm", "UPML zdim -", 0);


}

bool PML::UPML::getParameters(){

   using namespace std;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   
   // if (!Readparameters::get("UPML.width", this->upmlWidth)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.xp", this->pml_Xp)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.xm", this->pml_Xm)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.yp", this->pml_Yp)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.Ym", this->pml_Ym)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.zp", this->pml_Zp)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // }
   // if (!Readparameters::get("UPML.zm", this->pml_Zm)) {
   //    if (myRank == MASTER_RANK)
   //       cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
   //    exit(1);
   // // }


   // if (myRank ==MASTER_RANK){
   //    cout<<"-----UPML-------\n";
   //    cout<<"Width= "<<this->upmlWidth<<"\n";
   //    cout<<"Xp= "<<this->pml_Xp<<"\n";
   //    cout<<"Xm= "<<this->pml_Xm<<"\n";
   //    cout<<"Yp= "<<this->pml_Yp<<"\n";
   //    cout<<"Ym= "<<this->pml_Ym<<"\n";
   //    cout<<"Zp= "<<this->pml_Zp<<"\n";
   //    cout<<"Zm= "<<this->pml_Zm<<endl;;
   // }

   return true;

}



 PML::UPML::UPML(){
    this->addParameters();
    this->getParameters();
 }



void PML::UPML::update(FsGrid<std::array<Real, fsgrids::upml::N_UPML>, FS_STENCIL_WIDTH>& fsUpml,
                const FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
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
   this->upmlWidth = 10;
   this->pml_Xm = true;
   this->pml_Xp = false;
   this->pml_Ym = true;
   this->pml_Yp = false;
   this->pml_Zm = false;
   this->pml_Zp = false;


   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   


   const int* globalDims = &fsUpml.getGlobalSize()[0];
   const int* localDims = &fsUpml.getLocalSize()[0];

   // Zero Out PML 
   // #pragma omp parallel for collapse(3)
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

   // Real ds = 3;
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


   // if (myRank==0){
      // std::cout<<"Factor= "<<Parameters::upmlFactor<<"\n";
   //    std::cout<<"Sfactor= "<<sigfactor<<"\n";
   //    std::cout<<"Kfactor= "<<kfactor<<"\n";
   //    std::cout<<"Smax= "<<sigmam<<"\n";
   // }
 //X component
   for (int i=0; i < localDims[0];i++){
      for (int k=0; k < localDims[2];k++){
         for (int j=0; j < localDims[1]; j++){

            std::array<int32_t, 3> pos;
            pos = fsUpml.getGlobalIndices(i, j, k);

            //Skip cells that do not belong to the X- PML layers
            if (pos[0]>this->upmlWidth +offset || pos[0] < offset  ) continue;

            Real x1 = (this->upmlWidth - pos[0] + 1) * ds;
            Real x2 = (this->upmlWidth - pos[0]) * ds;

            Real ksi = (this->upmlWidth - pos[0]+offset) / this->upmlWidth;
            Real sigma = sigmam * pow(ksi, orderbc);
            Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

            // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
            // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));

            if (pos[1]==15 ) std::cout<<pos[0]<<" "<<ksi<<" "<<sigmam<<" "<<sigma<<std::endl;
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
      
         }
      }
   }
   

   for (int i=0; i < localDims[0];i++){
      for (int k=0; k < localDims[2];k++){
         for (int j=0; j < localDims[1]; j++){

            std::array<int32_t, 3> pos;
            pos = fsUpml.getGlobalIndices(i, j, k);

            //Skip cells that do not belong to the X- PML layers
            if (pos[0]>this->upmlWidth || pos[0] < offset) continue;

            Real x1 = (this->upmlWidth - pos[0] + 1.5) * ds;
            Real x2 = (this->upmlWidth - pos[0] + 0.5) * ds;

            Real ksi = (this->upmlWidth - pos[0] + .5 +offset) / this->upmlWidth;
            Real sigma = sigmam * pow(ksi, orderbc);
            Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

            // Real sigma = sigfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
            // Real ki = 1 + kfactor * (pow(x1, orderbc + 1) - pow(x2, orderbc + 1));
            Real facm = (2 * epsr * epsz * ki - sigma * dt);
            Real facp = (2 * epsr * epsz * ki + sigma * dt);

            std::array<Real, fsgrids::upml::N_UPML> *val;
            val = fsUpml.get(i, j, k);
            val->at(fsgrids::upml::C1EZ) = facm/facp;
            val->at(fsgrids::upml::C2EZ) = 2.0 * epsr * epsz * dt / facp;
            val->at(fsgrids::upml::C3EY) = facm / facp;
            val->at(fsgrids::upml::C4EY) = 1.0 / facp / epsr / epsz;
            val->at(fsgrids::upml::C5BX) = facp;
            val->at(fsgrids::upml::C6BX) = facm;
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






   //Y component
   for (int i=0; i < localDims[0];i++){
      for (int k=0; k < localDims[2];k++){
         for (int j=0; j < localDims[1]; j++){
      

            std::array<int32_t, 3> pos;
            pos = fsUpml.getGlobalIndices(i, j, k);

            //Skip cells that do not belong to the X- PML layers
            if (pos[1]>this->upmlWidth || pos[1] < offset) continue;

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
      
      
      
         }
      }
   }
   
   for (int i=0; i < localDims[0];i++){
      for (int k=0; k < localDims[2];k++){
         for (int j=0; j < localDims[1]; j++){

            std::array<int32_t, 3> pos;
            pos = fsUpml.getGlobalIndices(i, j, k);

            //Skip cells that do not belong to the X- PML layers
            if (pos[1]>this->upmlWidth || pos[1] < offset) continue;

            Real y1 = (this->upmlWidth - pos[1] + 1.5) * ds;
            Real y2 = (this->upmlWidth - pos[1] + 0.5) * ds;
            Real ksi = (this->upmlWidth - pos[1] + .5 +offset) / this->upmlWidth;
            Real sigma = sigmam * pow(ksi, orderbc);
            Real ki = 1 + (kmax - 1.0) * pow(ksi, orderbc);

            // Real sigma = sigfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
            // Real ki = 1 + kfactor * (pow(y1, orderbc + 1) - pow(y2, orderbc + 1));
            Real facm = (2 * epsr * epsz * ki - sigma * dt);
            Real facp = (2 * epsr * epsz * ki + sigma * dt);

            std::array<Real, fsgrids::upml::N_UPML>* val;
            val = fsUpml.get(i, j, k);
            val->at(fsgrids::upml::C1EX) = facm/facp;
            val->at(fsgrids::upml::C2EX) = 2.0 * epsr * epsz * dt / facp;
            val->at(fsgrids::upml::C3EZ) = facm / facp;
            val->at(fsgrids::upml::C4EZ) = 1.0 / facp / epsr / epsz;
            val->at(fsgrids::upml::C5BY) = facp;
            val->at(fsgrids::upml::C6BY) = facm;
         }
      }
   }




//Y PEC wall
   for (int i=0; i < localDims[0];i++){
      for (int k=0; k < localDims[2];k++){
         for (int j=0; j < localDims[1]; j++){

            std::array<int32_t, 3> pos;
            pos = fsUpml.getGlobalIndices(i, j, k);

            //Skip cells that do not belong to the X- PML layers
            if (pos[1] != offset  ) continue;

            std::array<Real, fsgrids::upml::N_UPML> *val;
            val = fsUpml.get(i,j,k);
            val->at(fsgrids::upml::C1EX) = -1.0;
            val->at(fsgrids::upml::C2EX) = 0.0;
            val->at(fsgrids::upml::C3EZ) = -1.0;
            val->at(fsgrids::upml::C4EZ) = 0.0;
         }
      }
   }




   // //IO
   // std::map<int, std::string> outputData;
   // this->outputVarsPML(outputData);
   // std::vector<std::map<int, std::string>> outputPML {outputData};
   // std::vector<FsGrid<std::array<Real, fsgrids::upml::N_UPML>, FS_STENCIL_WIDTH>*> vecpml {&fsUpml};

   // hdfIO<Real> IO;
   // IO.writeHDF(vecpml,99999, outputPML);
   // IO.Hello();

}
                                       



// void PML::UPML::outputVarsPML(std::map<int, std::string> &output){
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1EX,"C1EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2EX,"C2EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3EX,"C3EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4EX,"C4EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5EX,"C5EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6EX,"C6EX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1EY,"C1EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2EY,"C2EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3EY,"C3EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4EY,"C4EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5EY,"C5EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6EY,"C6EY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1EZ,"C1EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2EZ,"C2EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3EZ,"C3EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4EZ,"C4EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5EZ,"C5EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6EZ,"C6EZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1BX,"C1BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2BX,"C2BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3BX,"C3BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4BX,"C4BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5BX,"C5BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6BX,"C6BX"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1BY,"C1BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2BY,"C2BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3BY,"C3BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4BY,"C4BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5BY,"C5BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6BY,"C6BY"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C1BZ,"C1BZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C2BZ,"C2BZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C3BZ,"C3BZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C4BZ,"C4BZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C5BZ,"C5BZ"));
//    output.insert(std::pair<int, std::string>(fsgrids::upml::C6BZ,"C6BZ"));
   
// }
