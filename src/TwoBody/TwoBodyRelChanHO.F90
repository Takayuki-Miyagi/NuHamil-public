module TwoBodyRelChanHO
  use TwoBodyRelRadQNumbers
  implicit none

  public :: TwoBodyRelChanHOBasis
  private :: InitTwoBodyRelChanHOBasis
  private :: FinTwoBodyRelChanHOBasis
  private :: GetTwoBodyRelativeHORadial
  private :: GetTwoBodyRelativeHORadialFromNLS
  private :: GetJ
  private :: GetParity
  private :: GetZ
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNmax
  private :: GetFrequency

  type :: TwoBodyRelChanHOBasis
    integer, private :: j
    integer, private :: p
    integer, private :: z
    integer, private :: n_state
    real(8), private :: hw
    integer, private :: lmin, lmax, smin, smax, Nmax
    integer, private, allocatable :: nls2idx(:,:,:)
    type(HarmonicOscillator), allocatable :: HOn(:)
  contains
    procedure :: InitTwoBodyRelChanHOBasis
    procedure :: FinTwoBodyRelChanHOBasis
    procedure :: GetTwoBodyRelativeHORadial
    procedure :: GetTwoBodyRelativeHORadialFromNLS
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNmax
    procedure :: GetFrequency

    ! overriding
    generic :: init => InitTwoBodyRelChanHOBasis  ! constructor
    generic :: fin => FinTwoBodyRelChanHOBasis    ! destructor
    generic :: getp => GetTwoBodyRelativeHORadial,&
      & GetTwoBodyRelativeHORadialFromNLS         ! getter of HOn (pointer returned)
  end type TwoBodyRelChanHOBasis

contains

  subroutine FinTwoBodyRelChanHOBasis(this)
    ! two-body relative channel in harmonic oscilator basis
    ! destructor
    class(TwoBodyRelChanHOBasis), intent(inout) :: this
    if( this%n_state == 0 ) return
    deallocate(this%nls2idx)
    deallocate(this%HOn)
  end subroutine FinTwoBodyRelChanHOBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetFrequency(this) result(hw)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l,s) result(idx)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer, intent(in) :: n, l, s
    integer :: idx
    idx = this%nls2idx(n,l,s)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelChanHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetTwoBodyRelativeHORadialFromNLS(this, n, l, s) result(r)
    class(TwoBodyRelChanHOBasis), target, intent(in) :: this
    integer, intent(in) :: n, l, s
    type(HarmonicOscillator), pointer :: r
    r => this%getp( this%nls2idx(n, l, s) )
  end function GetTwoBodyRelativeHORadialFromNLS

  function GetTwoBodyRelativeHORadial(this, idx) result(r)
    class(TwoBodyRelChanHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(HarmonicOscillator), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%HOn(idx)
  end function GetTwoBodyRelativeHORadial

  subroutine InitTwoBodyRelChanHOBasis(this, Nmax, j, p, z, hw)
    ! two-body relative channel in harmonic oscilator basis
    ! constructor
    class(TwoBodyRelChanHOBasis), intent(inout) :: this
    integer, intent(in) :: Nmax, j, p, z
    real(8), intent(in) :: hw
    integer :: lmin, lmax
    integer :: smin, smax
    integer :: n12max, n12min
    integer :: n, l, cnt, s
    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%z = z
    this%hw = hw
    lmin = Nmax
    lmax = 0

    smin = 1
    smax = 0

    n12min = Nmax
    n12max = 0

    cnt = 0
    do s = 0, 1
      do l = iabs(j-s), j+s
        if((-1) ** l /= p) cycle
        if(l < 0) cycle
        do n = 0, Nmax
          if(2*n + l > Nmax) cycle
          if(iabs(z) == 1 .and. (-1) **(l + s) == -1) cycle
          cnt = cnt + 1
          if(smin > s) smin = s
          if(smax < s) smax = s

          if(lmin > l) lmin = l
          if(lmax < l) lmax = l

          if(n12min > n) n12min = n
          if(n12max < n) n12max = n
        end do
      end do
    end do
    this%n_state = cnt
    this%smin = smin
    this%smax = smax
    this%lmin = lmin
    this%lmax = lmax
    if(this%n_state == 0) return


    allocate(this%HOn(cnt))
    allocate(this%nls2idx(0:n12max, lmin:lmax, smin:smax))
    this%nls2idx(:,:,:) = 0

    cnt = 0
    do s = 0, 1
      do l = iabs(j-s), j+s
        if((-1) ** l /= p) cycle
        if(l < 0) cycle
        do n = 0, Nmax
          if(2*n + l > Nmax) cycle
          if(iabs(z) == 1 .and. (-1) **(l + s) == -1) cycle
          cnt = cnt + 1
          call this%HOn(cnt)%set(n,l,s)
          this%nls2idx(n,l,s) = cnt
#ifdef TwoBodyRelChanHOBasisDebug
          write(*,'(a,5i3)') 'InitTwoBodyRelChanHOBasis, |n,l,s,j,tz>:', &
              & n,l,s,j,z
#endif
        end do
      end do
    end do
  end subroutine InitTwoBodyRelChanHOBasis
end module TwoBodyRelChanHO
