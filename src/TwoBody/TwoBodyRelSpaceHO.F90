module TwoBodyRelSpaceHO
  use TwoBodyRelChanHO
  implicit none

  public :: TwoBodyRelSpaceHOBasis
  private :: InitTwoBodyRelSpaceHOBasis
  private :: FinTwoBodyRelSpaceHOBasis
  private :: GetTwoBodyRelChanHOBasisFromJPTZ
  private :: GetTwoBodyRelChanHOBasisFromCh
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetFrequency
  private :: GetJmax
  private :: GetNmax

  type :: TwoBodyRelSpaceHOBasis
    real(8), private :: hw
    integer, private :: Nmax, Jmax, NChan
    integer, private, allocatable :: jpz2idx(:,:,:)
    logical, private :: initialized = .false.
    type(TwoBodyRelChanHOBasis), allocatable :: jpz(:)
  contains
    procedure :: InitTwoBodyRelSpaceHOBasis
    procedure :: FinTwoBodyRelSpaceHOBasis
    procedure :: GetTwoBodyRelChanHOBasisFromJPTZ
    procedure :: GetTwoBodyRelChanHOBasisFromCh
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetFrequency
    procedure :: GetJmax
    procedure :: GetNmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceHOBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceHOBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanHOBasisFromJPTZ, &
        & GetTwoBodyRelChanHOBasisFromCh
  end type TwoBodyRelSpaceHOBasis

contains

  subroutine FinTwoBodyRelSpaceHOBasis(this)
    class(TwoBodyRelSpaceHOBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz(ich)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceHOBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceHOBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyRelSpaceHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, z) result(idx)
    class(TwoBodyRelSpaceHOBasis), intent(in) :: this
    integer, intent(in) :: j, p, z
    integer :: idx
    idx = this%jpz2idx(j, p, z)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelSpaceHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceHOBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceHOBasis(this, hw, Nmax, Jmax, j_in, p_in, tz_in)
    class(TwoBodyRelSpaceHOBasis), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: Nmax, Jmax
    integer, intent(in), optional :: j_in, p_in, tz_in
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
            cnt = 0
            do s = 0, 1
              do l = iabs(j-s), j+s
                if((-1) ** l /= ipar) cycle
                if(l < 0) cycle
                do n = 0, this%Nmax
                  if(2*n + l > this%Nmax) cycle
                  if(iabs(itz) == 1 .and. (-1) **(l + s) == -1) cycle
                  cnt = cnt + 1
                end do
              end do
            end do
#ifdef TwoBodyRelativeSpaceDebug
            if(cnt /= 0) then
              write(*,'(a,3i4)') 'InitTwoBodyRelSpaceHOBasis |Parity,J,Tz>', &
                  & ipar,j,itz
            end if
#endif
            if(cnt > 0) then
              ich = ich + 1
              if(loop == 2) then
                this%jpz2idx(j, ipar, itz) = ich
                call this%jpz(ich)%init(this%Nmax, j, ipar, itz, hw)
              end if
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpz(ich))
        allocate(this%jpz2idx(0:jmax, -1:1, -1:1))
        this%jpz2idx(:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(tz_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpace, spin, parity, and pn are not correct'
      stop
    end if
    this%initialized = .true.
  end subroutine InitTwoBodyRelSpaceHOBasis

  function GetTwoBodyRelChanHOBasisFromJPTZ(this, j, p, z) result(r)
    class(TwoBodyRelSpaceHOBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, z
    type(TwoBodyRelChanHOBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,z) )
  end function GetTwoBodyRelChanHOBasisFromJPTZ

  function GetTwoBodyRelChanHOBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceHOBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanHOBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpz(ch)
  end function GetTwoBodyRelChanHOBasisFromCh

end module TwoBodyRelSpaceHO

!program test
!  use RelativeSpaceHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceHOBasis) :: ho
!  type(TwoBodyRelativeChannelHOBasis), pointer :: ch
!  type(HarmonicOscillator), pointer :: m
!  integer :: ich, n
!
!  call ho%init(20,8)
!  do ich = 1, ho%NChan
!    ch => ho%getp(ich)
!    write(*,'(3i3)') ch%j, ch%p, ch%z
!  end do
!  call ho%fin()
!
!  call ho%init(20,8,1,1,0)
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
