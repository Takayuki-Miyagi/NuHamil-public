module ThreeBodyVMonIsoSingle
  use omp_lib
  use Profiler, only: timer
  use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
  use ClassSys
  use SingleParticleState
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyJacOpsIso
  use TransformationCoefficient
  implicit none
#define single_precision
#define PRECISION Single
#include "ThreeBodyVMonIso.inc"
#undef PRECISION
#undef single_precision
end module ThreeBodyVMonIsoSingle
