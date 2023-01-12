module TwoBodyLabSpacePN
  use SingleParticleState
  use TwoBodyLabChanPN
  implicit none
  public :: TwoBodyLabPNSpace

  private :: InitTwoBodyLabPNSpace
  private :: FinTwoBodyLabPNSpace
  private :: GetEmax
  private :: GetLmax
  private :: GetE2max
  private :: GetNumberChannels
  private :: GetFrequency
  private :: GetIndex
  private :: GetChannelFromIndex
  private :: GetChannelFromJPZ

  type :: TwoBodyLabPNSpace
    integer, private :: emax =-1
    integer, private :: e2max=-1
    integer, private :: lmax=-1
    integer, private :: NChan= 0
    real(8), private :: hw = -1.d0
    type(Orbits), pointer :: sps
    type(TwoBodyLabPNChan), allocatable :: jpz(:)
    integer, allocatable :: jpz2idx(:,:,:)
    integer :: mass_snt
    logical :: is_Constructed=.false.
  contains
    procedure :: InitTwoBodyLabPNSpace
    procedure :: FinTwoBodyLabPNSpace
    procedure :: GetEmax
    procedure :: GetE2max
    procedure :: GetLmax
    procedure :: GetNumberChannels
    procedure :: GetFrequency
    procedure :: GetIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromJPZ
    generic :: init => InitTwoBodyLabPNSpace
    generic :: fin => FinTwoBodyLabPNSpace
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromJPZ
  end type TwoBodyLabPNSpace
contains

  function GetEmax(this) result(emax)
    class(TwoBodyLabPNSpace), intent(in) :: this
    integer :: emax
    emax = this%emax
  end function GetEmax

  function GetLmax(this) result(lmax)
    class(TwoBodyLabPNSpace), intent(in) :: this
    integer :: lmax
    lmax = this%lmax
  end function GetLmax

  function GetE2max(this) result(e2max)
    class(TwoBodyLabPNSpace), intent(in) :: this
    integer :: e2max
    e2max = this%e2max
  end function GetE2max

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyLabPNSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyLabPNSpace), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, z) result(idx)
    class(TwoBodyLabPNSpace), intent(in) :: this
    integer, intent(in) :: j, p, z
    integer :: idx
    idx = this%jpz2idx(j,p,z)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(chan)
    class(TwoBodyLabPNSpace), intent(in), target :: this
    integer, intent(in) :: idx
    type(TwoBodyLabPNChan), pointer :: chan
    chan => null()
    if(idx == 0) return
    chan => this%jpz(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPZ(this,j,p,z) result(chan)
    class(TwoBodyLabPNSpace), intent(in), target :: this
    integer, intent(in) :: j, p, z
    type(TwoBodyLabPNChan), pointer :: chan
    chan => this%GetChannel( this%GetIndex(j,p,z) )
  end function GetChannelFromJPZ

  subroutine InitTwoBodyLabPNSpace(this, hw, sps, e2max, mass_snt)
    use MyLibrary, only: triag
    class(TwoBodyLabPNSpace), intent(inout) :: this
    real(8) :: hw
    type(Orbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    type(TwoBodyLabPNChan) :: dummy
    type(TwoBodyLabPNChan), pointer :: r
    integer :: mass_snt
    integer :: NChan, j, p, z, n_states

    this%hw = hw
    this%emax = sps%emax
    this%lmax = sps%lmax
    this%e2max = e2max
    this%mass_snt = mass_snt
    this%sps => sps
    allocate(this%jpz2idx(0:this%GetE2max()+1, -1:1, -1:1) )
    this%jpz2idx(:,:,:) = 0

    NChan = 0
    do j = 0, this%GetE2max() + 1
      do p = 1, -1, -2
        do z = -1, 1
          n_states = 0
          call dummy%init(sps, j, p, z, this%GetE2max(), n_states)
          if(n_states == 0) cycle
          NChan = NChan + 1
          this%jpz2idx(j,p,z) = NChan
        end do
      end do
    end do

    this%NChan = NChan
    allocate( this%jpz(NChan) )

    do j = 0, this%GetE2max() + 1
      do p = 1, -1, -2
        do z = -1, 1
          NChan = this%jpz2idx(j,p,z)
          if(NChan == 0) cycle
          r => this%GetChannel(j,p,z)
          n_states = 0
          call r%init(sps, j, p, z, this%GetE2max(), n_states)
          call r%init(sps, j, p, z, this%GetE2max(), n_states)
        end do
      end do
    end do
    this%is_Constructed = .true.
  end subroutine InitTwoBodyLabPNSpace

  subroutine FinTwoBodyLabPNSpace(this)
    class(TwoBodyLabPNSpace), intent(inout) :: this
    integer :: ch
    if(.not. this%is_Constructed) return
    do ch = 1, this%GetNumberChannels()
      call this%jpz(ch)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2idx)
    this%is_Constructed = .false.
    this%emax =-1
    this%e2max=-1
    this%NChan= 0
    this%hw = -1.d0
    this%sps => null()
  end subroutine FinTwoBodyLabPNSpace
end module TwoBodyLabSpacePN
