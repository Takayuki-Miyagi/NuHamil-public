module ThreeBodyTransCoefIsoNO2BSingle
  use omp_lib
  use LinAlgLib
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoNO2B
  use ThreeBodyTransCoefIsospinSingle
  implicit none
#define single_precision
#define PRECISION Single
#include "ThreeBodyTransCoefIsoNO2B.inc"
#undef PRECISION
#undef single_precision
end module ThreeBodyTransCoefIsoNO2BSingle

