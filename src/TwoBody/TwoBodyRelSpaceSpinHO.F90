module TwoBodyRelSpaceSpinHO
  use TwoBodyRelChanSpinHO
  implicit none

  public :: TwoBodyRelSpaceSpinHOBasis
  private :: InitTwoBodyRelSpaceSpinHOBasis
  private :: FinTwoBodyRelSpaceSpinHOBasis
  private :: GetTwoBodyRelChanSpinHOBasisFromJPSTZ
  private :: GetTwoBodyRelChanSpinHOBasisFromCh
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetFrequency
  private :: GetJmax
  private :: GetNmax

  type :: TwoBodyRelSpaceSpinHOBasis
    real(8), private :: hw
    integer, private :: Nmax, Jmax, NChan
    integer, private, allocatable :: jpsz2idx(:,:,:,:)
    logical, private :: initialized = .false.
    type(TwoBodyRelChanSpinHOBasis), allocatable :: jpsz(:)
  contains
    procedure :: InitTwoBodyRelSpaceSpinHOBasis
    procedure :: FinTwoBodyRelSpaceSpinHOBasis
    procedure :: GetTwoBodyRelChanSpinHOBasisFromJPSTZ
    procedure :: GetTwoBodyRelChanSpinHOBasisFromCh
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetFrequency
    procedure :: GetJmax
    procedure :: GetNmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceSpinHOBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceSpinHOBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanSpinHOBasisFromJPSTZ, &
        & GetTwoBodyRelChanSpinHOBasisFromCh
  end type TwoBodyRelSpaceSpinHOBasis

contains

  subroutine FinTwoBodyRelSpaceSpinHOBasis(this)
    class(TwoBodyRelSpaceSpinHOBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpsz(ich)%fin()
    end do
    deallocate(this%jpsz)
    deallocate(this%jpsz2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceSpinHOBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceSpinHOBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyRelSpaceSpinHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, s, z) result(idx)
    class(TwoBodyRelSpaceSpinHOBasis), intent(in) :: this
    integer, intent(in) :: j, p, z, s
    integer :: idx
    idx = this%jpsz2idx(j, p, s, z)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelSpaceSpinHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceSpinHOBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceSpinHOBasis(this, hw, Nmax, Jmax, j_in, p_in, s_in, tz_in)
    class(TwoBodyRelSpaceSpinHOBasis), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: Nmax, Jmax
    integer, intent(in), optional :: j_in, p_in, s_in, tz_in
    integer :: j, ipar, itz
    integer :: ich, cnt
    integer :: loop, s, l, n
    this%hw = hw
    this%Nmax = Nmax
    this%Jmax = Jmax
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
                do n = 0, this%Nmax
                  if(2*n + l > this%Nmax) cycle
                  if(iabs(itz) == 1 .and. (-1) **(l + s) == -1) cycle
                  cnt = cnt + 1
                end do
              end do
#ifdef TwoBodyRelSpaceSpinHOBasisDebug
              if(cnt /= 0) then
                write(*,'(a,5i4)') 'InitTwoBodyRelSpaceSpinHOBasis |parity,S,J,Tz>', &
                    & ipar,s,j,itz, cnt
              end if
#endif
              if(cnt > 0) then
                ich = ich + 1
                if(loop == 2) then
                  this%jpsz2idx(j, ipar, s, itz) = ich
                  call this%jpsz(ich)%init(this%Nmax, j, ipar, s, itz)
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
  end subroutine InitTwoBodyRelSpaceSpinHOBasis

  function GetTwoBodyRelChanSpinHOBasisFromJPSTZ(this, j, p, s, z) result(r)
    class(TwoBodyRelSpaceSpinHOBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, s, z
    type(TwoBodyRelChanSpinHOBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,s,z) )
  end function GetTwoBodyRelChanSpinHOBasisFromJPSTZ

  function GetTwoBodyRelChanSpinHOBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceSpinHOBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanSpinHOBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpsz(ch)
  end function GetTwoBodyRelChanSpinHOBasisFromCh
end module TwoBodyRelSpaceSpinHO

!program test
!  use RelativeSpaceSpinHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceSpinHOBasis) :: ho
!  type(TwoBodyRelativeChannelSpinHOBasis), pointer :: ch
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
