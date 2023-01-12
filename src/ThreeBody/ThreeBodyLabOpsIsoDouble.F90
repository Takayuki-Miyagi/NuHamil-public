module ThreeBodyLabOpsIsoDouble
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
  use TransformationCoefficient
  implicit none
#define PRECISION Double
#define double_precision
#include "ThreeBodyLabOpsIso.inc"
#undef double_precision
#undef PRECISION
end module ThreeBodyLabOpsIsoDouble
