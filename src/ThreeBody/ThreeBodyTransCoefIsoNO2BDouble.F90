module ThreeBodyTransCoefIsoNO2BDouble
  use omp_lib
  use LinAlgLib
  use ClassSys
  use StoreCouplings
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoNO2B
  use ThreeBodyTransCoefIsospinDouble
  implicit none
#define double_precision
#define PRECISION Double
#include "ThreeBodyTransCoefIsoNO2B.inc"
#undef PRECISION
#undef double_precision
end module ThreeBodyTransCoefIsoNO2BDouble

