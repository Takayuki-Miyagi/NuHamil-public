!
! Three-body matrix elements in laboratory frame
!
module ThreeBodyNO2BIsospinHalf
  use omp_lib
  use ClassSys
  use MyLibrary, only: gzip_open, gzip_writeline, gzip_close, gzip_readline
  use Profiler, only: timer
  use ThreeBodyJacOpsIso
  use ThreeBodyLabSpaceIsoNO2B
  use half_precision_floating_points
  use TransformationCoefficient
  implicit none
#define half_precision
#define PRECISION Half
#include "ThreeBodyNO2BIsospin.inc"
#undef PRECISION
#undef half_precision
end module ThreeBodyNO2BIsospinHalf
