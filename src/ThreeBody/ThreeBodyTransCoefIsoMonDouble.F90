module ThreeBodyTransCoefIsoMonDouble
  use omp_lib
  use ClassSys
  use LinAlgLib
  use StoreCouplings
  use SingleParticleState
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyTransCoefIsospinDouble
  implicit none
#define double_precision
#define PRECISION Double
#include "ThreeBodyTransCoefIsoMon.inc"
#undef PRECISION
#undef double_precision
end module ThreeBodyTransCoefIsoMonDouble
