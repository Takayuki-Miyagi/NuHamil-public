module ThreeBodyTransCoefIsospinHalf
  use omp_lib
  use MPIFunction
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpace
  use half_precision_floating_points
  implicit none
#define PRECISION Half
#define half_precision
#include "ThreeBodyTransCoefIsospin.inc"
#undef half_precision
#undef PRECISION
end module ThreeBodyTransCoefIsospinHalf
