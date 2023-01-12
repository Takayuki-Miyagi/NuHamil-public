module ThreeBodyLabChanIso
  use SingleParticleState
  implicit none

  public :: ThreeBodyLabIsoChan
  public :: ThreeBodyLabIsoSubChan

  private :: InitThreeBodyLabIsoSubChan
  private :: FinThreeBodyLabIsoSubChan
  private :: SetThreeBodyLabIsoSubChan
  private :: GetThreeBodyLabIsoSubIndex
  private :: InitThreeBodyLabIsoChan
  private :: FinThreeBodyLabIsoChan
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetABCJT
  private :: GetIndexFromInts
  private :: GetIndexFromArray
  private :: GetNumberStates
  private :: GetNumberABC

  type :: ThreeBodyLabIsoSubChan
    integer, private, allocatable :: SubIndex(:,:)
    integer, private :: jmin
    integer, private :: jmax
    integer, private :: tmin
    integer, private :: tmax
  contains
    procedure :: InitThreeBodyLabIsoSubChan
    procedure :: FinThreeBodyLabIsoSubChan
    procedure :: SetThreeBodyLabIsoSubChan
    procedure :: GetThreeBodyLabIsoSubIndex
    generic :: init => InitThreeBodyLabIsoSubChan
    generic :: fin => FinThreeBodyLabIsoSubChan
    generic :: SetIndex => SetThreeBodyLabIsoSubChan
    generic :: GetIndex => GetThreeBodyLabIsoSubIndex
  end type ThreeBodyLabIsoSubChan

  type :: ThreeBodyLabIsoChan
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: t = -1
    integer, private :: n_states = 0
    integer, private :: n_abc_states = 0
    type(OrbitsIsospin), pointer :: sps
    type(ThreeBodyLabIsoSubChan), allocatable, private :: idx(:)
    integer, private, allocatable :: ABC2Index(:,:,:)
    integer, private, allocatable :: IndexABCJT(:,:)
  contains
    procedure :: InitThreeBodyLabIsoChan
    procedure :: FinThreeBodyLabIsoChan
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetABCJT
    procedure :: GetIndexFromInts
    procedure :: GetIndexFromArray
    procedure :: GetNumberStates
    procedure :: GetNumberABC
    generic :: init => InitThreeBodyLabIsoChan
    generic :: fin => FinThreeBodyLabIsoChan
    generic :: GetIndex => GetIndexFromInts, GetIndexFromArray
  end type ThreeBodyLabIsoChan
  integer, private :: jmin_init = 10000
  integer, private :: jmax_init = -1
contains

  function GetJ(this) result(J)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer :: J
    J = this%j
  end function GetJ

  function GetParity(this) result(P)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer :: P
    P = this%p
  end function GetParity

  function GetT(this) result(T)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer :: T
    T = this%t
  end function GetT

  function GetABCJT(this, idx) result(ABCJT)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: idx
    integer :: ABCJT(5)
    ABCJT(:) = this%IndexABCJT(:,idx)
  end function GetABCJT

  function GetIndexFromInts(this, i1, i2, i3, j12, t12) result(idx)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: i1, i2, i3, j12, t12
    integer :: idx, iabc
    idx = 0
    iabc = this%ABC2Index(i1,i2,i3)
    if(iabc == 0) return
    if(i1==i2 .and. mod(j12+t12,2) == 0) return
    idx = this%idx(iabc)%GetIndex(j12,t12)
  end function GetIndexFromInts

  function GetIndexFromArray(this, abcJT ) result(idx)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: abcJT(5)
    integer :: idx, iabc
    idx = 0
    iabc = this%ABC2Index(abcJT(1),abcJT(2),abcJT(3))
    if(iabc == 0) return
    if(abcJT(1)==abcJT(2) .and. mod(abcJT(4)+abcJT(5),2) == 0) return
    idx = this%idx(iabc)%GetIndex(abcJT(4),abcJT(5))
  end function GetIndexFromArray

  function GetNumberStates(this) result(n_states)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetNumberABC(this) result(n_states)
    class(ThreeBodyLabIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_abc_states
  end function GetNumberABC

  function GetThreeBodyLabIsoSubIndex(this, j12, t12) result(idx)
    class(ThreeBodyLabIsoSubChan), intent(in) :: this
    integer, intent(in) :: j12, t12
    integer :: idx
    idx = this%SubIndex(j12,t12)
  end function GetThreeBodyLabIsoSubIndex

  subroutine FinThreeBodyLabIsoSubChan(this)
    class(ThreeBodyLabIsoSubChan), intent(inout) :: this
    this%jmin = jmin_init
    this%jmax = jmax_init
    this%tmin = jmin_init
    this%tmax = jmax_init
    if(.not. allocated(this%SubIndex) ) return
    deallocate(this%SubIndex)
  end subroutine FinThreeBodyLabIsoSubChan

  subroutine InitThreeBodyLabIsoSubChan(this, chan, i1, i2, i3)
    use MyLibrary, only: triag
    class(ThreeBodyLabIsoSubChan), intent(inout) :: this
    type(ThreeBodyLabIsoChan), intent(in) :: chan
    integer, intent(in) :: i1, i2, i3
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: j12, t12

    sps => chan%sps
    o1 => sps%GetOrbit(i1)
    o2 => sps%GetOrbit(i2)
    o3 => sps%GetOrbit(i3)
    this%jmin = jmin_init
    this%jmax = -1
    this%tmin = jmin_init
    this%tmax = -1

    do j12 = abs(o1%j-o2%j)/2, (o1%j+o2%j)/2
      if( triag(2*j12,o3%j,chan%GetJ() )) cycle
      do t12 = 0, 1
        if(triag(2*t12, 1, chan%GetT() )) cycle
        if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
        this%jmin = min(j12, this%jmin)
        this%jmax = max(j12, this%jmax)
        this%tmin = min(t12, this%tmin)
        this%tmax = max(t12, this%tmax)
      end do
    end do
    if(this%jmin > this%jmax) return
    if(this%tmin > this%tmax) return
    allocate( this%SubIndex(this%jmin:this%jmax, this%tmin:this%tmax) )
    this%SubIndex(:,:) = 0
  end subroutine InitThreeBodyLabIsoSubChan

  subroutine SetThreeBodyLabIsoSubChan(this, j12, t12, idx)
    class(ThreeBodyLabIsoSubChan), intent(inout) :: this
    integer, intent(in) :: j12, t12, idx
    this%SubIndex(j12,t12) = idx
  end subroutine SetThreeBodyLabIsoSubChan

  subroutine FinThreeBodyLabIsoChan(this)
    class(ThreeBodyLabIsoChan), intent(inout) :: this
    integer :: idx
    do idx = 1, this%GetNumberABC()
      call this%idx(idx)%fin()
    end do
    deallocate(this%ABC2Index)
    deallocate(this%IndexABCJT)
    this%j = -1
    this%p = -100
    this%t = -100
    this%n_states = 0
    this%n_abc_states = 0
  end subroutine FinThreeBodyLabIsoChan

  subroutine InitThreeBodyLabIsoChan(this, sps, j, p, t, e2max, e3max, set_mode, n_states, n_abc)
    class(ThreeBodyLabIsoChan), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: j, p, t, e2max, e3max
    logical, intent(in) :: set_mode
    integer, intent(inout) :: n_states, n_abc
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    type(ThreeBodyLabIsoSubChan) :: sub
    integer :: i1, i2, i3, j12, t12

    this%j = j
    this%p = p
    this%t = t
    this%sps => sps
    if(set_mode) then
      this%n_states = n_states
      this%n_abc_states = n_abc
      allocate(this%ABC2Index( sps%norbs, sps%norbs, sps%norbs) )
      allocate(this%IndexABCJT( 5,this%GetNumberStates() ))
      allocate(this%idx( this%GetNumberABC() ))
      this%ABC2Index(:,:,:) = 0
      this%IndexABCJT(:,:) = 0
    end if

    n_states = 0
    n_abc = 0
    do i1 = 1, sps%norbs
      o1 => sps%GetOrbit(i1)
      do i2 = 1, i1
        o2 => sps%GetOrbit(i2)
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, i2
          o3 => sps%GetOrbit(i3)
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          if( (-1) ** (o1%l + o2%l + o3%l) .ne. p) cycle

          call sub%init(this, i1, i2, i3)
          if(sub%jmin > sub%jmax) cycle
          if(sub%tmin > sub%tmax) cycle
          n_abc = n_abc + 1
          if(set_mode) then
            call this%idx(n_abc)%init(this, i1, i2, i3)
            do j12 = sub%jmin, sub%jmax
              do t12 = sub%tmin, sub%tmax
                if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
                n_states = n_states + 1
                this%ABC2Index( i1,i2,i3 ) = n_abc
                call this%idx(n_abc)%SetIndex(j12,t12,n_states)
                this%IndexABCJT(:,n_states) = [i1,i2,i3,j12,t12]
              end do
            end do
          end if

          if(.not. set_mode) then
            do j12 = sub%jmin, sub%jmax
              do t12 = sub%tmin, sub%tmax
                if(i1 == i2 .and. mod(j12+t12,2) == 0) cycle
                n_states = n_states + 1
              end do
            end do
          end if
          call sub%fin()
        end do
      end do
    end do
  end subroutine InitThreeBodyLabIsoChan

end module ThreeBodyLabChanIso
