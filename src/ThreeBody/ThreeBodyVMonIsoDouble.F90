module ThreeBodyVMonIsoDouble
  use omp_lib
  use Profiler, only: timer
  use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
  use ClassSys
  use SingleParticleState
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyJacOpsIso
  use TransformationCoefficient
  implicit none
#define double_precision
#define PRECISION Double
#include "ThreeBodyVMonIso.inc"
#undef PRECISION
#undef double_precision
end module ThreeBodyVMonIsoDouble
