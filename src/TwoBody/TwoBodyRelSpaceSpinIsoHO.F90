module TwoBodyRelSpaceSpinIsoHO
  use TwoBodyRelChanSpinIsoHO
  implicit none

  public :: TwoBodyRelSpaceSpinIsoHOBasis
  private :: InitTwoBodyRelSpaceSpinIsoHOBasis
  private :: FinTwoBodyRelSpaceSpinIsoHOBasis
  private :: GetTwoBodyRelChanSpinIsoHOBasisFromJPST
  private :: GetTwoBodyRelChanSpinIsoHOBasisFromCh
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetFrequency
  private :: GetJmax
  private :: GetNmax

  type :: TwoBodyRelSpaceSpinIsoHOBasis
    real(8) :: hw
    integer :: Nmax, Jmax, NChan
    integer, allocatable :: jpst2idx(:,:,:,:)
    logical :: initialized = .false.
    type(TwoBodyRelChanSpinIsoHOBasis), allocatable :: jpst(:)
  contains
    procedure :: InitTwoBodyRelSpaceSpinIsoHOBasis
    procedure :: FinTwoBodyRelSpaceSpinIsoHOBasis
    procedure :: GetTwoBodyRelChanSpinIsoHOBasisFromJPST
    procedure :: GetTwoBodyRelChanSpinIsoHOBasisFromCh
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetFrequency
    procedure :: GetJmax
    procedure :: GetNmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceSpinIsoHOBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceSpinIsoHOBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanSpinIsoHOBasisFromJPST, &
        & GetTwoBodyRelChanSpinIsoHOBasisFromCh
  end type TwoBodyRelSpaceSpinIsoHOBasis

contains

  subroutine FinTwoBodyRelSpaceSpinIsoHOBasis(this)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpst(ich)%fin()
    end do
    deallocate(this%jpst)
    deallocate(this%jpst2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceSpinIsoHOBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, s, t) result(idx)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(in) :: this
    integer, intent(in) :: j, p, s, t
    integer :: idx
    idx = this%jpst2idx(j, p, s, t)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax


  subroutine InitTwoBodyRelSpaceSpinIsoHOBasis(this, hw, Nmax, Jmax, j_in, p_in, s_in, t_in)
    class(TwoBodyRelSpaceSpinIsoHOBasis), intent(inout) :: this
    real(8) :: hw
    integer, intent(in) :: Nmax, Jmax
    integer, intent(in), optional :: j_in, p_in, s_in, t_in
    integer :: j, ipar, t
    integer :: ich, cnt
    integer :: loop, s, l, n
    this%hw = hw
    this%Nmax = Nmax
    this%Jmax = Jmax
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
                do n = 0, this%Nmax
                  if(2*n + l > this%Nmax) cycle
                  if((-1) ** (l+s+t) == 1) cycle
                  cnt = cnt + 1
                end do
#ifdef TwoBodyRelSpaceSpinIsoHOBasisDebug
                if(cnt /= 0) then
                  write(*,'(a,4i4)') 'InitTwoBodyRelSpaceSpinIsoHOBasis |L,S,J,Tz>', &
                      & l,s,j,itz
                end if
#endif
              end do
              if(cnt > 0) then
                ich = ich + 1
                if(loop == 2) then
                  this%jpst2idx(j, ipar, s, t) = ich
                  call this%jpst(ich)%init(this%Nmax, j, ipar, s, t)
                end if
              end if
            end do
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpst(ich))
        allocate(this%jpst2idx(0:jmax, -1:1, 0:1, 0:1))
        this%jpst2idx(:,:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(s_in) .and. present(t_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpaceSpinIsospin, spin, parity, and pn are not correct'
      stop
    end if
    this%initialized = .true.
  end subroutine InitTwoBodyRelSpaceSpinIsoHOBasis

  function GetTwoBodyRelChanSpinIsoHOBasisFromJPST(this, j, p, s, t) result(r)
    class(TwoBodyRelSpaceSpinIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, s, t
    type(TwoBodyRelChanSpinIsoHOBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,s,t) )
  end function GetTwoBodyRelChanSpinIsoHOBasisFromJPST

  function GetTwoBodyRelChanSpinIsoHOBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceSpinIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanSpinIsoHOBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpst(ch)
  end function GetTwoBodyRelChanSpinIsoHOBasisFromCh
end module TwoBodyRelSpaceSpinIsoHO

!program test
!  use RelativeSpaceSpinIsospinHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceSpinIsoHOBasis) :: ho
!  type(TwoBodyRelativeChannelSpinIsospinHOBasis), pointer :: ch
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
