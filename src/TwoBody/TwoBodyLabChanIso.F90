module TwoBodyLabChanIso
  use SingleParticleState
  implicit none

  public :: TwoBodyLabIsoChan

  private
  private :: InitTwoBodyLabIsoChan
  private :: FinTwoBodyLabIsoChan
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetNumberStates
  private :: GetIndex
  private :: GetPhase

  type :: TwoBodyLabIsoChan
    integer, private :: n_states = 0
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: t = -100
    integer, allocatable :: labels2n(:,:)
    integer, allocatable :: iphase(:,:)
    integer, allocatable :: n2label1(:)
    integer, allocatable :: n2label2(:)
  contains
    procedure :: InitTwoBodyLabIsoChan
    procedure :: FinTwoBodyLabIsoChan
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetPhase
    generic :: init => InitTwoBodyLabIsoChan
    generic :: fin => FinTwoBodyLabIsoChan
  end type TwoBodyLabIsoChan
contains

  function GetJ(this) result(J)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer :: J
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(T)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer :: T
    t = this%T
  end function GetT

  function GetNumberStates(this) result(n_states)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetIndex(this, a, b) result(idx)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: a, b
    integer :: idx
    idx = this%labels2n(a,b)
  end function GetIndex

  function GetPhase(this, a, b) result(phase)
    class(TwoBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: a, b
    real(8) :: phase
    phase = this%iphase(a,b)
  end function GetPhase

  subroutine InitTwoBodyLabIsoChan(this, sps, j, p, t, e2max, n_states)
    use MyLibrary, only: triag
    class(TwoBodyLabIsoChan), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: j, p, t, e2max
    integer, intent(inout) :: n_states
    type(SingleParticleOrbitIsospin), pointer :: o1, o2
    integer :: m, i1, i2, j1, l1, tz1, j2, l2, tz2
    integer :: ii1, ii2
    logical :: set

    set = .false.
    m = sps%norbs
    if(n_states > 0) then
      set = .true.
      this%n_states = n_states
      this%j = j
      this%p = p
      this%t = t
      allocate(this%n2label1( this%GetNumberStates() ))
      allocate(this%n2label2( this%GetNumberStates() ))
      allocate(this%labels2n(m,m))
      allocate(this%iphase(m,m))
      this%n2label1 = 0
      this%n2label2 = 0
      this%labels2n = 0
      this%iphase = 0
    end if

    n_states = 0
    do i1 = 1, m
      o1 => sps%orb(i1)
      j1 = o1%j
      l1 = o1%l
      do i2 = 1, m
        o2 => sps%orb(i2)
        j2 = o2%j
        l2 = o2%l
        if(o1%e + o2%e > e2max) cycle
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1 + l2) /= p) cycle
        if(i1 == i2 .and. mod(j+t, 2) == 1) cycle
        n_states = n_states + 1
        if(set) then
          ii1 = i1
          ii2 = i2
          this%n2label1(n_states) = ii1
          this%n2label2(n_states) = ii2
          this%labels2n(ii1,ii2) = n_states
          this%labels2n(ii2,ii1) = n_states
          this%iphase(ii1,ii2) = 1
          this%iphase(ii2,ii1) = (-1) ** ((j1 + j2) / 2 - j - t)
        end if
      end do
    end do
  end subroutine InitTwoBodyLabIsoChan

  subroutine FinTwoBodyLabIsoChan(this)
    class(TwoBodyLabIsoChan), intent(inout) :: this
    if(this%n_states == 0) return
    deallocate(this%n2label1)
    deallocate(this%n2label2)
    deallocate(this%labels2n)
    deallocate(this%iphase)
    this%n_states = 0
    this%j = -1
    this%p = -100
    this%t = -100
  end subroutine FinTwoBodyLabIsoChan

end module TwoBodyLabChanIso
