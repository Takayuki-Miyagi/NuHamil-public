module ThreeBodyLabFileDouble
  use omp_lib
  use ClassSys
  use Profiler
  use NuHamilInput, only: InputParameters
  use SingleParticleState
  use ThreeBodyLabSpace
  use ThreeBodyLabOpsIso
  implicit none
#define PRECISION Double
#define double_precision
#include "ThreeBodyLabFile.inc"
#undef single_precision
#undef PRECISION
end module ThreeBodyLabFileDouble
