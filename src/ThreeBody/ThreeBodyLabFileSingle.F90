module ThreeBodyLabFileSingle
  use omp_lib
  use ClassSys
  use Profiler
  use NuHamilInput, only: InputParameters
  use SingleParticleState
  use ThreeBodyLabSpace
  use ThreeBodyLabOpsIso
  implicit none
#define PRECISION Single
#define single_precision
#include "ThreeBodyLabFile.inc"
#undef single_precision
#undef PRECISION
end module ThreeBodyLabFileSingle
