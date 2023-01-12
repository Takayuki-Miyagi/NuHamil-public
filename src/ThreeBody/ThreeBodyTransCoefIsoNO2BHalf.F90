module ThreeBodyTransCoefIsoNO2BHalf
  use omp_lib
  use LinAlgLib
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoNO2B
  use ThreeBodyTransCoefIsospinHalf
  use half_precision_floating_points
  implicit none
#define half_precision
#define PRECISION HALF
#include "ThreeBodyTransCoefIsoNO2B.inc"
#undef PRECISION
#undef half_precision
end module ThreeBodyTransCoefIsoNO2BHalf

