module ThreeBodyTransCoefIsospinSingle
  use omp_lib
  use MPIFunction
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpace
  implicit none
#define PRECISION Single
#define single_precision
#include "ThreeBodyTransCoefIsospin.inc"
#undef single_precision
#undef PRECISION
end module ThreeBodyTransCoefIsospinSingle
