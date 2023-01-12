!
! Three-body matrix elements in laboratory frame
!
module ThreeBodyTensorNO2BIsospinSingle
  use omp_lib
  use Profiler, only: timer
  use MyLibrary, only: triag, sjs, gzip_open, gzip_writeline, gzip_close, gzip_readline
  use ClassSys
  use OperatorDefinitions
  use ThreeBodyJacOpsIso
  use ThreeBodyLabSpaceIsoNO2B
  use TransformationCoefficient
  implicit none
#define single_precision
#define PRECISION Single
#include "ThreeBodyTensorNO2BIsospin.inc"
#undef PRECISION
#undef single_precision
end module ThreeBodyTensorNO2BIsospinSingle
