#if !defined(nekrs_nekrs_hpp_)
#define nekrs_nekrs_hpp_

#define NEKRS_VERSION "19"
#define NEKRS_SUBVERSION "0"

#define EXIT(a)  { MPI_Finalize(); exit(a); } 

#include "libParanumal.hpp"

// std::to_string might be not accurate enough 
static string to_string_f(double a) {
  stringstream s;
  s << std::scientific << a;
  return s.str();
}

void runPlan4(ins_t *ins);
void restartRead(ins_t *ins, setupAide &options);
void report(ins_t *ins, dfloat time, int tstep);
ins_t *nekrsInsSetup(mesh_t *mesh, setupAide &options);
void nekrsInsSetScalarSolver(ins_t *ins, setupAide &options, occa::properties &kernelInfo);

libParanumal::setupAide parRead(std::string &setupFile, MPI_Comm comm); 

#endif
