module ThreeBodyTransCoefIsoMonSingle
  use omp_lib
  use ClassSys
  use LinAlgLib
  use StoreCouplings
  use SingleParticleState
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyTransCoefIsospinSingle
  implicit none
#define single_precision
#define PRECISION Single
#include "ThreeBodyTransCoefIsoMon.inc"
#undef PRECISION
#undef single_precision
end module ThreeBodyTransCoefIsoMonSingle
