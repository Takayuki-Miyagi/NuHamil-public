module TwoBodyLabChanPN
  use SingleParticleState
  implicit none

  public :: TwoBodyLabPNChan

  private
  private :: InitTwoBodyLabPNChan
  private :: FinTwoBodyLabPNChan
  private :: GetJ
  private :: GetParity
  private :: GetZ
  private :: GetNumberStates
  private :: GetIndex
  private :: GetPhase

  type :: TwoBodyLabPNChan
    integer, private :: n_states = 0
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: z = -100
    integer, allocatable :: labels2n(:,:)
    integer, allocatable :: iphase(:,:)
    integer, allocatable :: n2label1(:)
    integer, allocatable :: n2label2(:)
  contains
    procedure :: InitTwoBodyLabPNChan
    procedure :: FinTwoBodyLabPNChan
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetPhase
    generic :: init => InitTwoBodyLabPNChan
    generic :: fin => FinTwoBodyLabPNChan
  end type TwoBodyLabPNChan
contains

  function GetJ(this) result(J)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer :: J
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetNumberStates(this) result(n_states)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetIndex(this, a, b) result(idx)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer, intent(in) :: a, b
    integer :: idx
    idx = this%labels2n(a,b)
  end function GetIndex

  function GetPhase(this, a, b) result(phase)
    class(TwoBodyLabPNChan), intent(in) :: this
    integer, intent(in) :: a, b
    real(8) :: phase
    phase = this%iphase(a,b)
  end function GetPhase

  subroutine InitTwoBodyLabPNChan(this, sps, j, p, z, e2max, n_states)
    use MyLibrary, only: triag
    class(TwoBodyLabPNChan), intent(inout) :: this
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: j, p, z, e2max
    integer, intent(inout) :: n_states
    type(SingleParticleOrbit), pointer :: o1, o2
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
      this%z = z
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
      tz1 = o1%z
      do i2 = 1, m
        o2 => sps%orb(i2)
        j2 = o2%j
        l2 = o2%l
        tz2 = o2%z
        if(z == 0 .and. (tz1 == 1 .and. tz2 == -1)) cycle
        if(o1%e + o2%e > e2max) cycle
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1 + l2) /= p) cycle
        if(tz1 + tz2 /= 2 * z) cycle
        if(abs(z) == 1 .and. i1 < i2) cycle
        if(i1 == i2 .and. mod(j, 2) == 1) cycle
        n_states = n_states + 1
        if(set) then
          ii1 = i1
          ii2 = i2
          this%n2label1(n_states) = ii1
          this%n2label2(n_states) = ii2
          this%labels2n(ii1,ii2) = n_states
          this%labels2n(ii2,ii1) = n_states
          this%iphase(ii1,ii2) = 1
          this%iphase(ii2,ii1) = -(-1) ** ((j1 + j2) / 2 - j)
        end if
      end do
    end do
  end subroutine InitTwoBodyLabPNChan

  subroutine FinTwoBodyLabPNChan(this)
    class(TwoBodyLabPNChan), intent(inout) :: this
    if(this%n_states == 0) return
    deallocate(this%n2label1)
    deallocate(this%n2label2)
    deallocate(this%labels2n)
    deallocate(this%iphase)
    this%n_states = 0
    this%j = -1
    this%p = -100
    this%z = -100
  end subroutine FinTwoBodyLabPNChan

end module TwoBodyLabChanPN
