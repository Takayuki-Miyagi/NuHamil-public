module ThreeBodyLabFileHalf
  use omp_lib
  use ClassSys
  use half_precision_floating_points
  use Profiler
  use NuHamilInput, only: InputParameters
  use SingleParticleState
  use ThreeBodyLabSpace
  use ThreeBodyLabOpsIso
  implicit none
#define PRECISION Half
#define half_precision
#include "ThreeBodyLabFile.inc"
#undef half_precision
#undef PRECISION
end module ThreeBodyLabFileHalf
