module OneBodyManager
  use ClassSys
  use Profiler, only: timer
  use NuHamilInput, only: InputParameters
  use MPIFunction, only: myrank
  use omp_lib
  implicit none
  public :: OneBodyManagers
  private
contains
  subroutine OneBodyManagers()
    use NuHamilInput, only: params
    use SingleParticleState
    use OperatorDefinitions
    use OneBodyLabOps
    type(Orbits) :: sps
    type(OneBodyLabOp) :: op
    type(str) :: f
    type(sys) :: s
    integer :: i

    call sps%init(params%emax, params%lmax)
    do i = 1, size(params%Operators)
      call op%init(sps, params%hw, params%Operators(i))

      call op%set()
      f = op%GetFileName(params%file_name_n)
      call op%WriteOperator(f)
      call op%fin()
    end do
    call sps%fin()
  end subroutine OneBodyManagers
end module OneBodyManager
