!
! Three-body matrix elements in laboratory frame
!
module ThreeBodyTensorNO2BIsospinHalf
  use omp_lib
  use Profiler, only: timer
  use MyLibrary, only: triag, sjs, gzip_open, gzip_writeline, gzip_close, gzip_readline
  use ClassSys
  use OperatorDefinitions
  use ThreeBodyJacOpsIso
  use ThreeBodyLabSpaceIsoNO2B
  use half_precision_floating_points
  use TransformationCoefficient
  implicit none
#define half_precision
#define PRECISION Half
#include "ThreeBodyTensorNO2BIsospin.inc"
#undef PRECISION
#undef half_precision
end module ThreeBodyTensorNO2BIsospinHalf
