module TwoBodyRelSpaceSpinIsoMeshBasis
  use TwoBodyRelChanSpinIsoMeshBasis
  implicit none

  public :: TwoBodyRelSpaceSpinIsoMBasis
  private :: InitTwoBodyRelSpaceSpinIsoMBasis
  private :: FinTwoBodyRelSpaceSpinIsoMBasis
  private :: GetTwoBodyRelChanSpinIsoMBasisFromJPST
  private :: GetTwoBodyRelChanSpinIsoMBasisFromCh
  private :: SetMeshWeight
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetNMesh
  private :: GetJmax

  type :: TwoBodyRelSpaceSpinIsoMBasis
    integer :: NMesh, Jmax, NChan
    integer, allocatable :: jpst2idx(:,:,:,:)
    real(8) :: xmin = 0.d0, xmax = 0.d0
    logical :: initialized = .false.
    type(TwoBodyRelChanSpinIsoMBasis), allocatable :: jpst(:)
  contains
    procedure :: InitTwoBodyRelSpaceSpinIsoMBasis
    procedure :: FinTwoBodyRelSpaceSpinIsoMBasis
    procedure :: GetTwoBodyRelChanSpinIsoMBasisFromJPST
    procedure :: GetTwoBodyRelChanSpinIsoMBasisFromCh
    procedure :: SetMeshWeight
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetNMesh
    procedure :: GetJmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceSpinIsoMBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceSpinIsoMBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanSpinIsoMBasisFromJPST, &
        & GetTwoBodyRelChanSpinIsoMBasisFromCh
  end type TwoBodyRelSpaceSpinIsoMBasis

contains

  subroutine FinTwoBodyRelSpaceSpinIsoMBasis(this)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpst(ich)%fin()
    end do
    deallocate(this%jpst)
    deallocate(this%jpst2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceSpinIsoMBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, s, t) result(idx)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(in) :: this
    integer, intent(in) :: j, p, s, t
    integer :: idx
    idx = this%jpst2idx(j, p, s, t)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceSpinIsoMBasis(this, xmin, xmax, NMesh, Jmax, j_in, p_in, s_in, t_in)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(inout) :: this
    integer, intent(in) :: Nmesh, Jmax
    real(8), intent(in) :: xmin, xmax
    integer, intent(in), optional :: j_in, p_in, s_in, t_in
    integer :: j, ipar, t
    integer :: ich, cnt
    integer :: loop, s, l
    this%NMesh = NMesh
    this%Jmax = Jmax
    this%xmin = xmin
    this%xmax = xmax
    do loop = 1, 2
      ich = 0
      do t = 0, 1
        do j = 0, this%jmax
          do ipar = 1, -1, -2

            if(present(j_in)) then
              if(j /= j_in) cycle
            end if
            if(present(p_in)) then
              if(ipar /= p_in) cycle
            end if
            if(present(t_in)) then
              if(t /= t_in) cycle
            end if
            do s = 0, 1
              if(present(s_in)) then
                if(s /= s_in) cycle
              end if
              cnt = 0
              do l = iabs(j-s), j+s
                if((-1) ** l /= ipar) cycle
                if(l < 0) cycle
                if((-1) **(l + s + t) == 1) cycle
                cnt = cnt + 1
#ifdef TwoBodyRelSpaceSpinIsoMBasisDebug
                if(cnt /= 0) then
                  write(*,'(a,4i4)') 'InitTwoBodyRelSpaceSpinIsoMBasis |L,S,J,Tz>', &
                      & l,s,j,itz
                end if
#endif
              end do
              if(cnt > 0) then
                ich = ich + 1
                if(loop == 2) then
                  this%jpst2idx(j, ipar, s, t) = ich
                  call this%jpst(ich)%init(xmin, xmax, NMesh, j, ipar, s, t)
                end if
              end if
            end do
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpst(ich))
        allocate(this%jpst2idx(0:jmax, -1:1, 0:1, -1:1))
        this%jpst2idx(:,:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(s_in) .and. present(t_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpaceSpinIsospin, spin, parity, and pn are not correct'
      stop
    end if
    this%initialized = .true.
  end subroutine InitTwoBodyRelSpaceSpinIsoMBasis

  subroutine SetMeshWeight(this, x, w)
    class(TwoBodyRelSpaceSpinIsoMBasis), intent(inout) :: this
    real(8), intent(in) :: x(:), w(:)
    integer :: ch
    do ch = 1, this%NChan
      call this%jpst(ch)%SetMeshWeight(x,w)
    end do
  end subroutine SetMeshWeight

  function GetTwoBodyRelChanSpinIsoMBasisFromJPST(this, j, p, s, t) result(r)
    class(TwoBodyRelSpaceSpinIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, s, t
    type(TwoBodyRelChanSpinIsoMBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,s,t) )
  end function GetTwoBodyRelChanSpinIsoMBasisFromJPST

  function GetTwoBodyRelChanSpinIsoMBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceSpinIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanSpinIsoMBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpst(ch)
  end function GetTwoBodyRelChanSpinIsoMBasisFromCh
end module TwoBodyRelSpaceSpinIsoMeshBasis

!program test
!  use RelativeSpaceSpinIsospinHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceSpinIsoMBasis) :: ho
!  type(TwoBodyRelativeChannelSpinIsospinMeshBasis), pointer :: ch
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
