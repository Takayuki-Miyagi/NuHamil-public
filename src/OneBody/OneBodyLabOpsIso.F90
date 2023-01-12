module OneBodyLabOpsIso
  use omp_lib
  use ClassSys
  use Profiler, only: timer
  use LinAlgLib
  use SingleParticleState
  use OperatorDefinitions
  use OneBodyLabOps
  implicit none

  private
  public :: OneBodyLabOpIso

  type, extends(OperatorDef) :: OneBodyLabOpIso
    type(DMat) :: mat
    type(OrbitsIsospin), pointer :: sps
    real(8) :: hw
  contains
    procedure :: InitOneBodyLabOpIso
    procedure :: FinOneBodyLabOpIso
    procedure :: MakeIsospinSymmetric
    generic :: init => InitOneBodyLabOpIso
    generic :: fin => FinOneBodyLabOpIso
  end type OneBodyLabOpIso

contains
  subroutine FinOneBodyLabOpIso(this)
    class(OneBodyLabOpIso), intent(inout) :: this
    call this%mat%fin()
    this%sps => null()
  end subroutine FinOneBodyLabOpIso

  subroutine InitOneBodyLabOpIso(this, sps, jr, pr, tr, hw)
    class(OneBodyLabOpIso), intent(inout) :: this
    type(OrbitsIsospin), intent(in), target :: sps
    integer, intent(in) :: jr, pr, tr
    real(8), intent(in) :: hw
    this%sps => sps
    this%hw = hw
    call this%InitOpDef(.false., jr, pr, tr)
    call this%mat%zeros(sps%norbs, sps%norbs)
  end subroutine InitOneBodyLabOpIso

  subroutine MakeIsospinSymmetric(this, op)
    use MyLibrary, only: tjs
    class(OneBodyLabOpIso), intent(inout) :: this
    type(OneBodyLabOp), intent(in) :: op
    type(OrbitsIsospin), pointer :: isps
    type(Orbits), pointer :: sps
    integer :: i1, i2
    integer :: i1_pn, i2_pn
    type(SingleParticleOrbitIsospin), pointer :: o1, o2
    real(8) :: fact

    isps => this%sps
    sps => op%sps
    do i1 = 1, isps%norbs
      do i2 = 1, isps%norbs
        o1 => isps%GetOrbit(i1)
        o2 => isps%GetOrbit(i2)
        i1_pn = sps%nljz2idx(o1%n, o1%l, o1%j, 1)
        i2_pn = sps%nljz2idx(o2%n, o2%l, o2%j, 1)

        fact = tjs(1, 2*this%GetOpT(), 1, -1, 0, 1)
        this%mat%m(i1,i2) = op%mat%m(i1_pn,i2_pn) / fact
      end do
    end do
  end subroutine MakeIsospinSymmetric

end module OneBodyLabOpsIso
