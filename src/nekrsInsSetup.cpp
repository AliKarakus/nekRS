/*
The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"


ins_t *nekrsInsSetup(mesh_t *mesh, setupAide &options){

 ins_t *ins = new ins_t();

 string DINS, DHOLMES; 

 DINS.assign(getenv("NEKRS_INSTALL_DIR")); DINS +="/ins"; 
 // DINS.assign(getenv("NEKRS_INSTALL_DIR")); // move to new kernels... 
 DHOLMES.assign(getenv("NEKRS_LIBP_DIR"));
  
  ins->mesh = mesh;
  ins->options = options;
  ins->kernelInfo = new occa::properties();

  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);
  
  ins->TOMBO = 0; 
  if(options.compareArgs("TIME INTEGRATOR", "EXTBDF") || options.compareArgs("TIME INTEGRATOR", "TOMBO")) 
    ins->TOMBO = 1; 

  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure

  ins->frame = 0;
  ins->SNrk = 0;
  
  mesh->Nfields = 1; 

  ins->g0 =  1.0;

  if (ins->TOMBO) {
    ins->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));

    ins->extC = (dfloat*) calloc(3, sizeof(dfloat));
  }

  if (options.compareArgs("TIME INTEGRATOR", "TOMBO1")) {
    ins->Nstages = 1;
    ins->temporalOrder = 1;
    ins->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "TOMBO2")) { 
    ins->Nstages = 2;
    ins->temporalOrder = 2;
    ins->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "TOMBO3")) { 
    ins->Nstages = 3;
    ins->temporalOrder = 3;
    ins->g0 = 11.f/6.f;
  }

  ins->readRestartFile = 0; 
  options.getArgs("RESTART FROM FILE", ins->readRestartFile);
  
  ins->writeRestartFile = 0; 
  options.getArgs("WRITE RESTART FILE", ins->writeRestartFile);



  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  ins->Ntotal = Ntotal;
  ins->fieldOffset = Ntotal;
  ins->Nblock = (Nlocal+blockSize-1)/blockSize;

  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(ins->NVfields*ins->Nstages*Ntotal,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(              ins->Nstages*Ntotal,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

  //additional field storage
  ins->NU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->LU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->GP   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  ins->GU   = (dfloat*) calloc(ins->NVfields*Ntotal*4,sizeof(dfloat));
  
  ins->rkU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  ins->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  ins->rkNU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkLU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkGP = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  ins->Vort = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->Div  = (dfloat*) calloc(              Nlocal,sizeof(dfloat));

  ins->FU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  //extra storage for interpolated fields
  if(ins->elementType==HEXAHEDRA)
    ins->cU = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    ins->cU = ins->U;

  ins->Nsubsteps = 0;
  if (ins->TOMBO)
    options.getArgs("SUBCYCLING STEPS",ins->Nsubsteps);

  if(ins->Nsubsteps){
    ins->Ud    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->Ue    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->resU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->rhsUd = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

    if(ins->elementType==HEXAHEDRA)
      ins->cUd = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      ins->cUd = ins->U;

    // Prepare RK stages for Subcycling Part
    
    int Sorder = 4; // Defaulting to LSERK 4(5) 
    
    options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder==2){
      ins->SNrk     = 2; 
      dfloat rka[2] = {0.0,     1.0 };
      dfloat rkb[2] = {0.5,     0.5 };
      dfloat rkc[2] = {0.0,     1.0 };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder ==3){
      // Using Williamson 3rd order scheme converted to low storage since the better truncation 
      ins->SNrk     = 3; 
      dfloat rka[3] = {0.0,     -5.0/9.0,  -153.0/128.0};
      dfloat rkb[3] = {1.0/3.0, 15.0/16.0,    8.0/15.0 };
      dfloat rkc[3] = {0.0,      1.0/3.0,     3.0/4.0  };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else{
      ins->SNrk     = 5; 
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      // Asumes initialized in mesh, can be moved here
      for(int rk=0; rk<ins->SNrk; rk++){
        ins->Srka[rk] = mesh->rka[rk]; 
        ins->Srkb[rk] = mesh->rkb[rk]; 
        ins->Srkc[rk] = mesh->rkc[rk]; 
      }
    }
  }

  dfloat rho  = 1.0 ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", ins->ubar);
  options.getArgs("VBAR", ins->vbar);
  if (ins->dim==3)
    options.getArgs("WBAR", ins->wbar);
  options.getArgs("PBAR", ins->pbar);
  options.getArgs("VISCOSITY", ins->nu);

  //Reynolds number
  ins->Re = ins->ubar/ins->nu;

  occa::properties& kernelInfo = *ins->kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(ins->dim==3){
    if(ins->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo); 
  } 
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;
  occa::properties kernelInfoS  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_pbar"]= ins->pbar;
  kernelInfo["defines/" "p_ubar"]= ins->ubar;
  kernelInfo["defines/" "p_vbar"]= ins->vbar;
  kernelInfo["defines/" "p_wbar"]= ins->wbar;
  kernelInfo["defines/" "p_nu"]= ins->nu;

  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  ins->Nstages;
  if(ins->Nsubsteps)
    kernelInfo["defines/" "p_SUBCYCLING"]=  1;
  else
    kernelInfo["defines/" "p_SUBCYCLING"]=  0;

  if(ins->elementType==QUADRILATERALS && mesh->dim==3){
    kernelInfo["defines/" "p_fainv"] = (dfloat) 0.0;
    kernelInfo["defines/" "p_invRadiusSq"] = (dfloat) 1./(mesh->sphereRadius*mesh->sphereRadius);
  }

  if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
    ins->ARKswitch = 1;   
  else 
    ins->ARKswitch = 0;

  // Struct for BC implementation
  kernelInfo["includes"] += DINS + "/data/insBcData.h";

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();
  
  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

#if 0
  if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);
#endif
  

  if(options.compareArgs("INITIAL CONDITION", "BROWN-MINION") &&
     (ins->elementType == QUADRILATERALS && ins->dim==3)){
    printf("Setting up initial condition for BROWN-MINION test case...");
    insBrownMinionQuad3D(ins);
    ins->o_U.copyFrom(ins->U);
    ins->o_P.copyFrom(ins->P);
    printf("done\n");

    char fname[BUFSIZ];
    string outName;
    ins->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, ins->frame++);
    insPlotVTU(ins, fname);
  }else{
      // set kernel name suffix
  string fileName, kernelName;
   

    for (int r=0;r<2;r++){
      if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
  if (ins->dim==2) {
    fileName = DINS + "/okl/insSetFlowField2D.okl"; 
    kernelName = "insSetFlowField2D"; 
  }
  else{
    fileName = DINS + "/okl/insSetFlowField3D.okl"; 
    kernelName = "insSetFlowField3D"; 
      }
    ins->setFlowFieldKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      }
      MPI_Barrier(mesh->comm);
    }

    ins->startTime =0.0;
    options.getArgs("START TIME", ins->startTime);
    ins->setFlowFieldKernel(mesh->Nelements,
          ins->startTime,
          mesh->o_x,
          mesh->o_y,
          mesh->o_z,
          ins->fieldOffset,
          ins->o_U,
          ins->o_P);
    ins->o_U.copyTo(ins->U);
    
  }
  
  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){

    if(ins->elementType==TRIANGLES || ins->elementType == TETRAHEDRA){
      for(int f=0;f<mesh->Nfaces;++f){
  dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
  dfloat sJ   = mesh->sgeo[sid + SJID];
  dfloat invJ = mesh->sgeo[sid + IJID];

  dfloat hest = 2./(sJ*invJ);

  hmin = mymin(hmin, hest);
  hmax = mymax(hmax, hest);
      }
    }else{
      for(int f=0;f<mesh->Nfaces;++f){
  for(int n=0; n<mesh->Nfp; n++){
    dlong sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n);
    dfloat sJ   = mesh->sgeo[sid + SJID];
    dfloat invJ = mesh->sgeo[sid + IJID];

    dfloat hest = 2./(sJ*invJ);

    hmin = mymin(hmin, hest);
    hmax = mymax(hmax, hest);
  }
      }
    }

    // dfloat maxMagVecLoc = 0;

    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = ins->U[id+0*ins->fieldOffset];
      dfloat uyn = ins->U[id+1*ins->fieldOffset];
      dfloat uzn = 0.0;
      if (ins->dim==3) uzn = ins->U[id+2*ins->fieldOffset];


      //Squared maximum velocity
      dfloat numax;
      if (ins->dim==2)
  numax = uxn*uxn + uyn*uyn;
      else 
  numax = uxn*uxn + uyn*uyn + uzn*uzn;

      umax = mymax(umax, numax);
    }
  }

  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  options.getArgs("CFL", ins->cfl);
  dfloat dt     = ins->cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
  
  ins->dtAdaptStep = 0; 
  options.getArgs("TSTEPS FOR TIME STEP ADAPT", ins->dtAdaptStep);

  options.getArgs("DT", dt);
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dti), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  // save initial time-step estimate 
  ins->dt = ins->dti; 

  options.getArgs("FINAL TIME", ins->finalTime);
  options.getArgs("START TIME", ins->startTime);
  
  if (ins->TOMBO ){
    ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;

    if(ins->Nsubsteps){
      ins->dt         = ins->Nsubsteps*ins->dt;
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
      ins->sdt        = ins->dt/ins->Nsubsteps;
    } else{
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
    }
  }

  ins->dtMIN = 1E-2*ins->dt; //minumum allowed timestep

  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  ins->lambda = ins->g0 / (ins->dt * ins->nu);

  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);

  ins->outputForceStep = 0;
  options.getArgs("TSTEPS FOR FORCE OUTPUT", ins->outputForceStep);

  
  //make option objects for elliptc solvers
  ins->vOptions = options;
  ins->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  ins->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
  ins->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  ins->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  ins->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  ins->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  ins->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  ins->vOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("VELOCITY MULTIGRID CHEBYSHEV DEGREE"));

  ins->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  ins->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  ins->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));
  ins->vOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("VELOCITY PARALMOND CHEBYSHEV DEGREE"));

  ins->vOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("VELOCITY PARALMOND AGGREGATION STRATEGY"));
  
  ins->vOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->vOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  ///  std::cout << "vOptions: " << ins->vOptions << std::endl;
  
  ins->pOptions = options;
  ins->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  ins->pOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("PRESSURE SOLVER TOLERANCE"));
  ins->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  ins->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  ins->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  ins->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  ins->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  ins->pOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE"));

  ins->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  ins->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  ins->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));
  ins->pOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("PRESSURE PARALMOND CHEBYSHEV DEGREE"));

  ins->pOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("PRESSURE PARALMOND AGGREGATION STRATEGY"));
  
  ins->pOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->pOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  ///  std::cout << "pOptions: " << ins->pOptions << std::endl;  
  
  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int uBCType[7] = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int vBCType[7] = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int wBCType[7] = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[7] = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.
  // int pBCType[7] = {0,2,2,2,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances 
  ins->presTOL = 1E-8;
  ins->velTOL  = 1E-8;


  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  ins->uSolver = new elliptic_t();
  ins->uSolver->mesh = mesh;
  ins->uSolver->options = ins->vOptions;
  ins->uSolver->dim = ins->dim;
  ins->uSolver->elementType = ins->elementType;
  ins->uSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->uSolver->BCType,uBCType,7*sizeof(int));

  ellipticSolveSetup(ins->uSolver, ins->lambda, kernelInfoV); 

  ins->vSolver = new elliptic_t();
  ins->vSolver->mesh = mesh;
  ins->vSolver->options = ins->vOptions;
  ins->vSolver->dim = ins->dim;
  ins->vSolver->elementType = ins->elementType;
  ins->vSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->vSolver->BCType,vBCType,7*sizeof(int));
  
  ellipticSolveSetup(ins->vSolver, ins->lambda, kernelInfoV); //!!!!!

  
  if (ins->dim==3) {
    ins->wSolver = new elliptic_t();
    ins->wSolver->mesh = mesh;
    ins->wSolver->options = ins->vOptions;
    ins->wSolver->dim = ins->dim;
    ins->wSolver->elementType = ins->elementType;
    ins->wSolver->BCType = (int*) calloc(7,sizeof(int));
    memcpy(ins->wSolver->BCType,wBCType,7*sizeof(int));
    
    ellipticSolveSetup(ins->wSolver, ins->lambda, kernelInfoV);  //!!!!!
  }
  
  if (mesh->rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  ins->pSolver = new elliptic_t();
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,7*sizeof(int));
  
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoP); //!!!!

  // TW: this code needs to be re-evaluated from here ====>
  //make node-wise boundary flags
  dfloat largeNumber = 1<<20;
  ins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  ins->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) ins->VmapB[n+e*mesh->Np] = largeNumber;
  }

  for (int e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
  for (int n=0;n<mesh->Nfp;n++) {
    int fid = mesh->faceNodes[n+f*mesh->Nfp];
    ins->VmapB[fid+e*mesh->Np] = mymin(bc,ins->VmapB[fid+e*mesh->Np]);
    ins->PmapB[fid+e*mesh->Np] = mymax(bc,ins->PmapB[fid+e*mesh->Np]);
  }
      }
    }
  }

  // this is potentially not good for Neumann (since it will )
  ogsGatherScatter(ins->VmapB, ogsInt, ogsMin, mesh->ogs); 
  ogsGatherScatter(ins->PmapB, ogsInt, ogsMax, mesh->ogs); 
#if 1
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {

    if (ins->VmapB[n] == largeNumber) {
      ins->VmapB[n] = 0.;
    }
  }
#endif
 
  ins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->VmapB);
  ins->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->PmapB);


  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // if (ins->TOMBO){
    kernelInfo["defines/" "p_EXTBDF"]= 0;
    kernelInfo["defines/" "p_TOMBO"]= 1;
  // }

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};
    ins->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 
    ins->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 
    ins->o_prkA = ins->o_extbdfC;
    ins->o_prkB = ins->o_extbdfC;

  // dummy decleration for scratch space 
  ins->Wrk     = (dfloat*) calloc(1, sizeof(dfloat));
  ins->o_Wrk   = mesh->device.malloc(1*sizeof(dfloat), ins->Wrk);
 
  // MEMORY ALLOCATION
  ins->o_rhsU  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsP);
  //storage for helmholtz solves
  ins->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  ins->o_FU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->FU);
  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);

  // Set Lumped Mass Matrix
    ins->o_NC = ins->o_GP; // Use GP storage to store curl(curl(u)) history
    // build lumped mass matrix for NEK
    dfloat *lumpedMassMatrix     = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
    dfloat *copyLumpedMassMatrix = (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
    for(hlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        lumpedMassMatrix[e*mesh->Np+n]     = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
        copyLumpedMassMatrix[e*mesh->Np+n] = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
      }
    }

    ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, mesh->ogs);

    for(int n=0;n<mesh->Np*mesh->Nelements;++n)
      lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];

  ins->o_invLumpedMassMatrix = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), lumpedMassMatrix);

  // Need to be revised for Tet/Tri
  ins->o_InvM = ins->o_invLumpedMassMatrix; 

  

  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);
  //plotting fields
  ins->o_Vort = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Vort);
  ins->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), ins->Div);


  if(ins->elementType==HEXAHEDRA)
    ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);
  else 
    ins->o_cU = ins->o_U;

  if(mesh->totalHaloPairs){//halo setup
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*ins->NVfields*sizeof(dfloat);
    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes = ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);

    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);
#if 1
    ins->vSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vSendBuffer, ins->h_vSendBuffer);
    ins->vRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vRecvBuffer, ins->h_vRecvBuffer);

    ins->pSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pSendBuffer, ins->h_pSendBuffer);
    ins->pRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pRecvBuffer, ins->h_pRecvBuffer);

    ins->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, ins->o_gatherTmpPinned, ins->h_gatherTmpPinned);
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  ins->velocityHaloGatherTmp);
#else
    ins->vSendBuffer = (dfloat*) malloc(vHaloBytes);
    ins->vRecvBuffer = (dfloat*) malloc(vHaloBytes);
    ins->pSendBuffer = (dfloat*) malloc(pHaloBytes);
    ins->pRecvBuffer = (dfloat*) malloc(pHaloBytes);
#endif
  }

  ins->Nscalar = 0; 
  options.getArgs("NUMBER of SCALAR FIELD", ins->Nscalar);

  // Now set the scalar Solver, make sure time step size etc are set before.... 
  // if(ins->Nscalar)
    // insSetScalarSolver(ins, options, kernelInfoS); 

   if(ins->Nscalar)
    nekrsInsSetScalarSolver(ins, options, kernelInfoS); 

  // set kernel name suffix
  string suffix, fileName, kernelName;
  
  if(ins->elementType==QUADRILATERALS)
     suffix = "Quad2D";
  if(ins->elementType==HEXAHEDRA)
     suffix = "Hex3D";

  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      
      //      std::cout << "NEW " << " " <<  kernelInfo << std::endl;
      fileName = DINS + "/okl/insHaloExchange.okl";
      kernelName = "insVelocityHaloExtract";
      ins->velocityHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insVelocityHaloScatter";
      ins->velocityHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insPressureHaloExtract";
      ins->pressureHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insPressureHaloScatter";
      ins->pressureHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // ===========================================================================
      fileName = DINS + "/okl/insAdvection" + suffix + ".okl";
      kernelName = "insAdvectionCubatureVolume" + suffix;
      ins->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insAdvectionCubatureSurface" + suffix;
      ins->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insAdvectionVolume" + suffix;
      ins->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insAdvectionSurface" + suffix;
      ins->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insStrongAdvectionVolume" + suffix;
      ins->advectionStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insStrongAdvectionCubatureVolume" + suffix;
      ins->advectionStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // ===========================================================================
      if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION")){
        fileName = DINS + "/okl/insFilter.okl";
        kernelName = "insFilterRT" + suffix;
        ins->filterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // if(options.compareArgs("TIME INTEGRATOR","TOMBO")){
        
        // Create lap(U) kernel for Ins for (lap(u) - lambda*u) operation
        fileName = DINS + "/okl/insAx" + suffix + ".okl"; 
        kernelName = "insPressureAx" + suffix;
        ins->pressureAxKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
        
        // needed for velocity increament
        kernelName = "insVelocityAx" + suffix;
        ins->velocityAxKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

        // Curl operations...
        fileName = DINS + "/okl/insCurl" + suffix + ".okl"; 
        kernelName = "insCurl" + suffix;
        ins->curlKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
        
        if(ins->dim==2){
          kernelName = "insCurlB " + suffix;
          ins->curlBKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
        }

  // Curl operations...
        fileName = DINS + "/okl/insMassMatrix.okl"; 
        kernelName = "insMassMatrix" + suffix;
        ins->massMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

        kernelName = "insInvMassMatrix" + suffix;
        ins->invMassMatrixKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  

      // }


      fileName = DINS + "/okl/insGradient" +suffix + ".okl"; 
      kernelName = "insGradientVolume" + suffix;
      ins->gradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName, "insGradientSurface" + suffix;
      ins->gradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================
      fileName = DINS + "/okl/insDivergence" + suffix + ".okl";

      if(options.compareArgs("TIME INTEGRATOR","TOMBO")){

        kernelName = "insDivergenceVolumeTOMBO" + suffix;
        ins->divergenceVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        kernelName = "insDivergenceSurfaceTOMBO" + suffix;
        ins->divergenceSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      }else{
        // kernelName = "insDivergenceVolume" + suffix;
        // ins->divergenceVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

        // kernelName = "insDivergenceSurface" + suffix;
        // ins->divergenceSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // ===========================================================================
      
      fileName = DINS + "/okl/insPressureRhs" + suffix + ".okl";
      if(ins->TOMBO){
        kernelName = "insPressureRhsTOMBO" + suffix;
      }else{
        kernelName = "insPressureRhs" + suffix;
      }

      ins->pressureRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
        
      // if(!(ins->dim==3 && ins->elementType==QUADRILATERALS) ){
  fileName = DINS + "/okl/insPressureBC" + suffix + ".okl" ; 
  // kernelName = "insPressureIpdgBC" + suffix;
  // ins->pressureRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "insPressureBC" + suffix;
  ins->pressureRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


  if(ins->TOMBO){
    kernelName = "insPressureAddBCTOMBO" + suffix;
  }else{
    // kernelName = "insPressureAddBC" + suffix;
  }
  ins->pressureAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      // }

      fileName = DINS + "/okl/insPressureUpdate.okl";
      if(ins->TOMBO){
  kernelName = "insPressureUpdateTOMBO";
      }else{
  // kernelName = "insPressureUpdate";
      }
      ins->pressureUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      // AK. Note that currently only HEX and QUAD is implemented in weak form !!!!!!!
      fileName = DINS + "/okl/insDivergence" + suffix + ".okl"; 
      // if(ins->elementType==HEXAHEDRA || (ins->elementType==QUADRILATERALS && ins->dim==2)){
  kernelName = "insStrongDivergenceVolume" + suffix;
  ins->divergenceStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  //     }else{ // implement other divergence operators in weak form also!!!!
  // // kernelName = "insDivergenceVolume" + suffix;
  // // ins->divergenceStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
  //     }

     
      
      // ===========================================================================
      fileName = DINS + "/okl/insVorticity" + suffix + ".okl";
      kernelName = "insVorticity" + suffix;
      ins->vorticityKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      // ===========================================================================

      fileName = DHOLMES + "/okl/scaledAdd.okl";
      kernelName = "scaledAddwOffset";
      ins->scaledAddKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
    
      // ===========================================================================
      fileName = DINS + "/okl/insVelocityRhs" + suffix + ".okl";
      kernelName = "insVelocityRhsTOMBO" + suffix;

      ins->velocityRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // ===========================================================================
        fileName = DINS + "/okl/insVelocityBC" + suffix + ".okl"; 

        kernelName = "insVelocityBC" + suffix;
        ins->velocityRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
        kernelName = "insVelocityAddBC" + suffix;
        ins->velocityAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      fileName = DINS + "/okl/insVelocityUpdate.okl";
      kernelName = "insVelocityUpdate";
      ins->velocityUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);   

      fileName = DINS + "/okl/insError" + suffix + ".okl"; 
      kernelName = "insError" + suffix;
      ins->errorKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      kernelName = "insSetFlowFieldCub" + suffix;
      ins->setFlowFieldCubKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // Not implemented for Quad 3D yet !!!!!!!!!!
      if(ins->Nsubsteps){
  // Note that resU and resV can be replaced with already introduced buffer
  ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
  ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
  ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
  ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

  if(ins->elementType==HEXAHEDRA)
    ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
  else 
    ins->o_cUd = ins->o_Ud;
      
  fileName =  DINS + "/okl/insSubCycle" + suffix + ".okl"; 
  kernelName = "insSubCycleVolume" + suffix;
  ins->subCycleVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  kernelName =  "insSubCycleSurface" + suffix; 
  ins->subCycleSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  kernelName =  "insSubCycleCubatureVolume" + suffix;
  ins->subCycleCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  
    kernelName =  "insSubCycleStrongCubatureVolume" + suffix;
    ins->subCycleStrongCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
    kernelName = "insSubCycleStrongVolume" + suffix;
    ins->subCycleStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "insSubCycleCubatureSurface" + suffix;
  ins->subCycleCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  fileName = DINS + "/okl/insSubCycle.okl";
  kernelName = "insSubCycleRKUpdate";
  ins->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "insSubCycleExt";
  ins->subCycleExtKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
      

      fileName = DINS + "/okl/insHalo.okl"; 
      kernelName = "insHaloGet"; 
      ins->haloGetKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      kernelName = "insHaloPut";       
      ins->haloPutKernel =mesh->device.buildKernel(fileName, kernelName, kernelInfo);


      fileName = DHOLMES + "/okl/addScalar.okl";
      kernelName = "setScalar";
      ins->setScalarKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // CFL related Kernels
      fileName = DINS + "/okl/insCfl" + suffix + ".okl"; 
      kernelName = "insCfl" + suffix;
      ins->cflKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      fileName = DINS + "/okl/insMax.okl";
      kernelName = "insMax";
      ins->maxKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

    }
    MPI_Barrier(mesh->comm);
  }

  // Initialize Cfl Stuff
  ins->cflComputed = 0; 
   
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    insFilterSetup(ins); 

  if(!ins->options.compareArgs("EXACT","NONE")){ // check if there is an exact solution
    ins->Uer  = (dfloat*) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp, sizeof(dfloat));
    ins->Per  = (dfloat*) calloc(              mesh->Nelements*mesh->cubNp, sizeof(dfloat));

    ins->o_Uex = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->Uer); 
    ins->o_Pex = mesh->device.malloc(              mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->Per); 
  }



  return ins;
}

