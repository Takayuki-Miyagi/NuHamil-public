!
! Three-body matrix elements in laboratory frame
!
module ThreeBodyNO2BIsospinSingle
  use omp_lib
  use ClassSys
  use MyLibrary, only: gzip_open, gzip_writeline, gzip_close, gzip_readline
  use Profiler, only: timer
  use ThreeBodyJacOpsIso
  use ThreeBodyLabSpaceIsoNO2B
  use TransformationCoefficient
  implicit none
#define single_precision
#define PRECISION Single
#include "ThreeBodyNO2BIsospin.inc"
#undef PRECISION
#undef single_precision
end module ThreeBodyNO2BIsospinSingle
