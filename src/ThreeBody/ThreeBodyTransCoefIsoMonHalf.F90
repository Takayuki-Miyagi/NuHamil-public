module ThreeBodyTransCoefIsoMonHalf
  use omp_lib
  use ClassSys
  use LinAlgLib
  use StoreCouplings
  use SingleParticleState
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyTransCoefIsospinHalf
  use half_precision_floating_points
  implicit none
#define half_precision
#define PRECISION Half
#include "ThreeBodyTransCoefIsoMon.inc"
#undef PRECISION
#undef half_precision
end module ThreeBodyTransCoefIsoMonHalf
