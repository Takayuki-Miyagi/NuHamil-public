module ThreeBodyTransCoefIsospinDouble
  use omp_lib
  use MPIFunction
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpace
  implicit none
#define double_precision
#define PRECISION Double
#include "ThreeBodyTransCoefIsospin.inc"
#undef PRECISION
#undef double_precision
end module ThreeBodyTransCoefIsospinDouble
