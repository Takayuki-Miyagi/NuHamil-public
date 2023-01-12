module TwoBodyLabSpaceIso
  use SingleParticleState
  use TwoBodyLabChanIso
  implicit none
  public :: TwoBodyLabIsoSpace

  private :: InitTwoBodyLabIsoSpace
  private :: FinTwoBodyLabIsoSpace
  private :: GetEmax
  private :: GetLmax
  private :: GetE2max
  private :: GetNumberChannels
  private :: GetFrequency
  private :: GetIndex
  private :: GetChannelFromIndex
  private :: GetChannelFromJPT

  type :: TwoBodyLabIsoSpace
    integer, private :: emax =-1
    integer, private :: e2max=-1
    integer, private :: lmax=-1
    integer, private :: NChan= 0
    real(8), private :: hw = -1.d0
    type(OrbitsIsospin), pointer :: sps
    type(TwoBodyLabIsoChan), allocatable :: jpt(:)
    integer, allocatable :: jpt2idx(:,:,:)
    integer :: mass_snt
    logical :: is_Constructed=.false.
  contains
    procedure :: InitTwoBodyLabIsoSpace
    procedure :: FinTwoBodyLabIsoSpace
    procedure :: GetEmax
    procedure :: GetE2max
    procedure :: GetLmax
    procedure :: GetNumberChannels
    procedure :: GetFrequency
    procedure :: GetIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromJPT
    generic :: init => InitTwoBodyLabIsoSpace
    generic :: fin => FinTwoBodyLabIsoSpace
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromJPT
  end type TwoBodyLabIsoSpace
contains

  function GetEmax(this) result(emax)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    integer :: emax
    emax = this%emax
  end function GetEmax

  function GetLmax(this) result(lmax)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    integer :: lmax
    lmax = this%lmax
  end function GetLmax

  function GetE2max(this) result(e2max)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    integer :: e2max
    e2max = this%e2max
  end function GetE2max

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, t) result(idx)
    class(TwoBodyLabIsoSpace), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j,p,t)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(chan)
    class(TwoBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    type(TwoBodyLabIsoChan), pointer :: chan
    chan => null()
    if(idx == 0) return
    chan => this%jpt(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this,j,p,t) result(chan)
    class(TwoBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: j, p, t
    type(TwoBodyLabIsoChan), pointer :: chan
    chan => this%GetChannel( this%GetIndex(j,p,t) )
  end function GetChannelFromJPT

  subroutine InitTwoBodyLabIsoSpace(this, hw, sps, e2max)
    use MyLibrary, only: triag
    class(TwoBodyLabIsoSpace), intent(inout) :: this
    real(8) :: hw
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: e2max
    type(TwoBodyLabIsoChan) :: dummy
    type(TwoBodyLabIsoChan), pointer :: r
    integer :: NChan, j, p, t, n_states

    this%hw = hw
    this%emax = sps%emax
    this%lmax = sps%lmax
    this%e2max = e2max
    this%sps => sps
    allocate(this%jpt2idx(0:this%GetE2max()+1, -1:1, 0:1) )
    this%jpt2idx(:,:,:) = 0

    NChan = 0
    do j = 0, this%GetE2max() + 1
      do p = 1, -1, -2
        do t = 0, 1
          n_states = 0
          call dummy%init(sps, j, p, t, this%GetE2max(), n_states)
          if(n_states == 0) cycle
          NChan = NChan + 1
          this%jpt2idx(j,p,t) = NChan
        end do
      end do
    end do

    this%NChan = NChan
    allocate( this%jpt(NChan) )

    do j = 0, this%GetE2max() + 1
      do p = 1, -1, -2
        do t = 0, 1
          NChan = this%jpt2idx(j,p,t)
          if(NChan == 0) cycle
          r => this%GetChannel(j,p,t)
          n_states = 0
          call r%init(sps, j, p, t, this%GetE2max(), n_states)
          call r%init(sps, j, p, t, this%GetE2max(), n_states)
        end do
      end do
    end do
    this%is_Constructed = .true.
  end subroutine InitTwoBodyLabIsoSpace

  subroutine FinTwoBodyLabIsoSpace(this)
    class(TwoBodyLabIsoSpace), intent(inout) :: this
    integer :: ch
    if(.not. this%is_Constructed) return
    do ch = 1, this%GetNumberChannels()
      call this%jpt(ch)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%is_Constructed = .false.
    this%emax =-1
    this%e2max=-1
    this%NChan= 0
    this%hw = -1.d0
    this%sps => null()
  end subroutine FinTwoBodyLabIsoSpace
end module TwoBodyLabSpaceIso
