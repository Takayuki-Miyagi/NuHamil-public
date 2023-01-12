module TwoBodyRelSpaceIsoHO
  use TwoBodyRelChanIsoHO
  implicit none

  public :: TwoBodyRelSpaceIsoHOBasis
  private :: InitTwoBodyRelSpaceIsoHOBasis
  private :: FinTwoBodyRelSpaceIsoHOBasis
  private :: GetTwoBodyRelChanIsospinHOBasisFromJPT
  private :: GetTwoBodyRelChanIsospinHOBasisFromCh
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetFrequency
  private :: GetJmax
  private :: GetNmax

  type :: TwoBodyRelSpaceIsoHOBasis
    real(8), private :: hw
    integer, private :: Nmax, Jmax, NChan
    integer, private, allocatable :: jpt2idx(:,:,:)
    logical, private :: initialized = .false.
    type(TwoBodyRelChanIsoHOBasis), allocatable :: jpt(:)
  contains
    procedure :: InitTwoBodyRelSpaceIsoHOBasis
    procedure :: FinTwoBodyRelSpaceIsoHOBasis
    procedure :: GetTwoBodyRelChanIsospinHOBasisFromJPT
    procedure :: GetTwoBodyRelChanIsospinHOBasisFromCh
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetFrequency
    procedure :: GetJmax
    procedure :: GetNmax
    ! overriding
    generic :: init => InitTwoBodyRelSpaceIsoHOBasis ! constructor
    generic :: fin => FinTwoBodyRelSpaceIsoHOBasis   ! destructor
    generic :: GetChannel => GetTwoBodyRelChanIsospinHOBasisFromJPT, &
        & GetTwoBodyRelChanIsospinHOBasisFromCh
  end type TwoBodyRelSpaceIsoHOBasis

contains

  subroutine FinTwoBodyRelSpaceIsoHOBasis(this)
    class(TwoBodyRelSpaceIsoHOBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceIsoHOBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceIsoHOBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetFrequency(this) result(hw)
    class(TwoBodyRelSpaceIsoHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, t) result(idx)
    class(TwoBodyRelSpaceIsoHOBasis), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j, p, t)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelSpaceIsoHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceIsoHOBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceIsoHOBasis(this, hw, Nmax, Jmax, j_in, p_in, t_in)
    class(TwoBodyRelSpaceIsoHOBasis), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: Nmax, Jmax
    integer, intent(in), optional :: j_in, p_in, t_in
    integer :: j, ipar, it
    integer :: ich, cnt
    integer :: loop, s, l, n
    this%hw = hw
    this%Nmax = Nmax
    this%Jmax = Jmax
    do loop = 1, 2
      ich = 0
      do it = 0, 1
        do j = 0, this%jmax
          do ipar = 1, -1, -2

            if(present(j_in)) then
              if(j /= j_in) cycle
            end if
            if(present(p_in)) then
              if(ipar /= p_in) cycle
            end if
            if(present(t_in)) then
              if(it /= t_in) cycle
            end if
            cnt = 0
            do s = 0, 1
              do l = iabs(j-s), j+s
                if((-1) ** l /= ipar) cycle
                if(l < 0) cycle
                do n = 0, this%Nmax
                  if(2*n + l > this%Nmax) cycle
                  if((-1) **(l + s + it) == 1) cycle
                  cnt = cnt + 1
                end do
#ifdef TwoBodyRelativeSpaceDebug
                if(cnt /= 0) then
                  write(*,'(a,4i4)') 'InitTwoBodyRelSpaceIsoHOBasis |L,S,J,T>', &
                      & l,s,j,it
                end if
#endif
              end do
            end do
            if(cnt > 0) then
              ich = ich + 1
              if(loop == 2) then
                this%jpt2idx(j, ipar, it) = ich
                call this%jpt(ich)%init(this%Nmax, j, ipar, it, hw)
              end if
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpt(ich))
        allocate(this%jpt2idx(0:jmax, -1:1, 0:1))
        this%jpt2idx(:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(t_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpace, spin, parity, and pn are not correct'
      stop
    end if
    this%initialized = .true.
  end subroutine InitTwoBodyRelSpaceIsoHOBasis

  function GetTwoBodyRelChanIsospinHOBasisFromJPT(this, j, p, t) result(r)
    class(TwoBodyRelSpaceIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, t
    type(TwoBodyRelChanIsoHOBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j, p, t) )
  end function GetTwoBodyRelChanIsospinHOBasisFromJPT

  function GetTwoBodyRelChanIsospinHOBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanIsoHOBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpt(ch)
  end function GetTwoBodyRelChanIsospinHOBasisFromCh

end module TwoBodyRelSpaceIsoHO

!program test
!  use RelativeSpaceHarmonicOscillator
!  implicit none
!  type(TwoBodyRelativeSpaceHOBasis) :: ho
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
