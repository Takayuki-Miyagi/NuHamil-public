module ThreeBodyLabOpsIsoHalf
  use omp_lib
  use LinAlgLib
  use MPIFunction
  use ClassSys
  use MyLibrary, only: triag, sjs, snj, hat, gzip_open, gzip_writeline, gzip_close, gzip_readline
  use Profiler, only: timer
  use OperatorDefinitions
  use ThreeBodyJacobiSpace
  use ThreeBodyJacOpsIso
  use ThreeBodyLabSpace
  use half_precision_floating_points
  use TransformationCoefficient
  implicit none
#define PRECISION Half
#define half_precision
#include "ThreeBodyLabOpsIso.inc"
#undef half_precision
#undef PRECISION
end module ThreeBodyLabOpsIsoHalf
