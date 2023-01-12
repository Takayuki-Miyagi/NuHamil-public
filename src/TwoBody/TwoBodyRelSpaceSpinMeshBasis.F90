module TwoBodyRelSpaceSpinMeshBasis
  use TwoBodyRelChanSpinMeshBasis
  implicit none

  public :: TwoBodyRelSpaceSpinMBasis
  private :: InitTwoBodyRelSpaceSpinMBasis
  private :: FinTwoBodyRelSpaceSpinMBasis
  private :: GetTwoBodyRelChanSpinMBasisFromJPSTZ
  private :: GetTwoBodyRelChanSpinMBasisFromCh
  private :: SetMeshWeight
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetNMesh
  private :: GetJmax

  type :: TwoBodyRelSpaceSpinMBasis
    integer :: NMesh, Jmax, NChan
    integer, allocatable :: jpsz2idx(:,:,:,:)
    real(8) :: xmin = 0.d0, xmax = 0.d0
    logical :: initialized = .false.
    type(TwoBodyRelChanSpinMBasis), allocatable :: jpsz(:)
  contains
    procedure :: InitTwoBodyRelSpaceSpinMBasis
    procedure :: FinTwoBodyRelSpaceSpinMBasis
    procedure :: GetTwoBodyRelChanSpinMBasisFromJPSTZ
    procedure :: GetTwoBodyRelChanSpinMBasisFromCh
    procedure :: SetMeshWeight
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetNMesh
    procedure :: GetJmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceSpinMBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceSpinMBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanSpinMBasisFromJPSTZ, &
        & GetTwoBodyRelChanSpinMBasisFromCh
  end type TwoBodyRelSpaceSpinMBasis

contains

  subroutine FinTwoBodyRelSpaceSpinMBasis(this)
    class(TwoBodyRelSpaceSpinMBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpsz(ich)%fin()
    end do
    deallocate(this%jpsz)
    deallocate(this%jpsz2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceSpinMBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceSpinMBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, s, z) result(idx)
    class(TwoBodyRelSpaceSpinMBasis), intent(in) :: this
    integer, intent(in) :: j, p, s, z
    integer :: idx
    idx = this%jpsz2idx(j, p, s, z)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelSpaceSpinMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceSpinMBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceSpinMBasis(this, xmin, xmax, NMesh, Jmax, j_in, p_in, s_in, tz_in)
    class(TwoBodyRelSpaceSpinMBasis), intent(inout) :: this
    integer, intent(in) :: Nmesh, Jmax
    real(8), intent(in) :: xmin, xmax
    integer, intent(in), optional :: j_in, p_in, s_in, tz_in
    integer :: j, ipar, itz
    integer :: ich, cnt
    integer :: loop, s, l
    this%NMesh = NMesh
    this%Jmax = Jmax
    this%xmin = xmin
    this%xmax = xmax
    do loop = 1, 2
      ich = 0
      do itz = -1, 1
        do j = 0, this%jmax
          do ipar = 1, -1, -2

            if(present(j_in)) then
              if(j /= j_in) cycle
            end if
            if(present(p_in)) then
              if(ipar /= p_in) cycle
            end if
            if(present(tz_in)) then
              if(itz /= tz_in) cycle
            end if
            do s = 0, 1
              if(present(s_in)) then
                if(s /= s_in) cycle
              end if
              cnt = 0
              do l = iabs(j-s), j+s
                if((-1) ** l /= ipar) cycle
                if(l < 0) cycle
                if(iabs(itz) == 1 .and. (-1) **(l + s) == -1) cycle
                cnt = cnt + 1
              end do
#ifdef TwoBodyRelSpaceSpinMBasisDebug
              if(cnt /= 0 .and. loop == 2) then
                write(*,'(a,4i4)') 'InitTwoBodyRelSpaceSpinMBasis |Parity,S,J,Tz>', &
                    & ipar,s,j,itz
              end if
#endif
              if(cnt > 0) then
                ich = ich + 1
                if(loop == 2) then
                  this%jpsz2idx(j, ipar, s, itz) = ich
                  call this%jpsz(ich)%init(NMesh, j, ipar, s, itz)
                end if
              end if
            end do
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpsz(ich))
        allocate(this%jpsz2idx(0:jmax, -1:1, 0:1, -1:1))
        this%jpsz2idx(:,:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(s_in) .and. present(tz_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpaceSpin, spin, parity, and pn are not correct'
      stop
    end if
    this%initialized = .true.
  end subroutine InitTwoBodyRelSpaceSpinMBasis

  subroutine SetMeshWeight(this, x, w)
    class(TwoBodyRelSpaceSpinMBasis), intent(inout) :: this
    real(8), intent(in) :: x(:), w(:)
    integer :: ch
    do ch = 1, this%NChan
      call this%jpsz(ch)%SetMeshWeight(x,w)
    end do
  end subroutine SetMeshWeight

  function GetTwoBodyRelChanSpinMBasisFromJPSTZ(this, j, p, s, z) result(r)
    class(TwoBodyRelSpaceSpinMBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, s, z
    type(TwoBodyRelChanSpinMBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,s,z) )
  end function GetTwoBodyRelChanSpinMBasisFromJPSTZ

  function GetTwoBodyRelChanSpinMBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceSpinMBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanSpinMBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpsz(ch)
  end function GetTwoBodyRelChanSpinMBasisFromCh
end module TwoBodyRelSpaceSpinMeshBasis

!program test
!  use RelativeSpaceSpinHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceSpinMBasis) :: ho
!  type(TwoBodyRelativeChannelSpinMeshBasis), pointer :: ch
!  type(HarmonicOscillator), pointer :: m
!  integer :: ich, n
!
!  call ho%init(20,8)
!  do ich = 1, ho%NChan
!    ch => ho%getp(ich)
!    write(*,'(4i3)') ch%j, ch%p, ch%s, ch%z
!  end do
!  call ho%fin()
!
!  write(*,*)
!  call ho%init(20,8,j_in=1,p_in=1,s_in=1,tz_in=0)
!  do ich = 1, ho%NChan
!    ch => ho%getp(ich)
!    write(*,'(3i3)') ch%j, ch%p, ch%z
!    do n = 1, ch%n_state
!      m => ch%getp(n)
!      write(*,'(3i3)') m%n, m%l, m%s
!    end do
!  end do
!  call ho%fin()
!end program test
