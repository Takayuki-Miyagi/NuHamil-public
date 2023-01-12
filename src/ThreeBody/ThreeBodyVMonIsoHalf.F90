module ThreeBodyVMonIsoHalf
  use omp_lib
  use Profiler, only: timer
  use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
  use ClassSys
  use SingleParticleState
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyJacOpsIso
  use half_precision_floating_points
  use TransformationCoefficient
  implicit none
#define half_precision
#define PRECISION Half
#include "ThreeBodyVMonIso.inc"
#undef PRECISION
#undef half_precision
end module ThreeBodyVMonIsoHalf
