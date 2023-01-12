module ThreeBodyLabSpaceIso
  use SingleParticleState
  use ThreeBodyLabChanIso
  implicit none

  public ThreeBodyLabIsoSpace

  private :: InitThreeBodyLabIsoSpace
  private :: FinThreeBodyLabIsoSpace
  private :: GetIndex
  private :: GetChannelFromJPT
  private :: GetChannelFromIndex
  private :: GetNumberChannels
  private :: GetEmax
  private :: GetE2max
  private :: GetE3max
  private :: GetFrequency

  type :: ThreeBodyLabIsoSpace
    type(ThreeBodyLabIsoChan), private, allocatable :: jpt(:)
    type(OrbitsIsospin), pointer :: sps
    integer, private, allocatable :: jpt2idx(:,:,:)
    integer, private :: emax = -1
    integer, private :: e2max = -1
    integer, private :: e3max = -1
    integer, private :: NChan = 0
    real(8), private :: hw = 0.d0
    logical :: is_Constructed=.false.
  contains
    procedure :: InitThreeBodyLabIsoSpace
    procedure :: FinThreeBodyLabIsoSpace
    procedure :: GetIndex
    procedure :: GetChannelFromJPT
    procedure :: GetChannelFromIndex
    procedure :: GetNumberChannels
    procedure :: GetEmax
    procedure :: GetE2max
    procedure :: GetE3max
    procedure :: GetFrequency
    generic :: init => InitThreeBodyLabIsoSpace
    generic :: fin => FinThreeBodyLabIsoSpace
    generic :: GetChannel => GetChannelFromJPT, GetChannelFromIndex
  end type ThreeBodyLabIsoSpace
contains

  function GetIndex(this, j, p, t) result(idx)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j,p,t)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(r)
    class(ThreeBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    type(ThreeBodyLabIsoChan), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%jpt(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this, j, p, t) result(r)
    class(ThreeBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: j, p, t
    type(ThreeBodyLabIsoChan), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,t) )
  end function GetChannelFromJPT

  function GetNumberChannels(this) result(NChan)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetEmax(this) result(emax)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    integer :: emax
    emax = this%emax
  end function GetEmax

  function GetE2max(this) result(e2max)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    integer :: e2max
    e2max = this%e2max
  end function GetE2max

  function GetE3max(this) result(e3max)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    integer :: e3max
    e3max = this%e3max
  end function GetE3max

  function GetFrequency(this) result(hw)
    class(ThreeBodyLabIsoSpace), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  subroutine FinThreeBodyLabIsoSpace(this)
    class(ThreeBodyLabIsoSpace), intent(inout) :: this
    integer :: ich
    if(.not. this%is_Constructed) return
    do ich = 1, this%GetNumberChannels()
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%sps => null()
    this%emax = -1
    this%e2max = -1
    this%e3max = -1
    this%NChan = 0
    this%hw = 0.d0
    this%is_Constructed=.false.
  end subroutine FinThreeBodyLabIsoSpace

  subroutine InitThreeBodyLabIsoSpace(this, sps, e2max, e3max, hw)
    use MyLibrary, only: triag
    class(ThreeBodyLabIsoSpace), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    real(8), intent(in) :: hw
    type(ThreeBodyLabIsoChan) :: sub
    integer :: ich, jtot, ptot, ttot, n_states, n_abc
    this%sps => sps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max
    this%hw = hw
    allocate(this%jpt2idx(2 * e3max + 3, -1:1, 1:3))
    this%jpt2idx(:,:,:) = 0

    ich = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2
          call sub%init(sps, jtot, ptot, ttot, e2max, e3max, .false., n_states, n_abc)
          if(n_states == 0) cycle
          if(n_abc == 0) cycle
          ich = ich + 1
          this%jpt2idx(jtot,ptot,ttot) = ich
        end do
      end do
    end do
    this%NChan = ich
    allocate(this%jpt( this%GetNumberChannels() ))

    do ttot = 1, 3, 2
      do jtot = 1, 2 * e3max + 3, 2
        do ptot = 1, -1, -2
          ich = this%GetIndex(jtot,ptot,ttot)
          if(ich == 0) cycle
          call sub%init(sps, jtot, ptot, ttot, e2max, e3max, .false., n_states, n_abc)
          call this%jpt(ich)%init(sps, jtot, ptot, ttot, e2max, e3max, .true., n_states, n_abc)
        end do
      end do
    end do
    this%is_Constructed=.true.
  end subroutine InitThreeBodyLabIsoSpace

end module ThreeBodyLabSpaceIso
