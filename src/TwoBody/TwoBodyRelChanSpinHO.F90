module TwoBodyRelChanSpinHO
  use TwoBodyRelRadQNumbers
  implicit none

  public :: TwoBodyRelChanSpinHOBasis
  private :: InitTwoBodyRelChanSpinHOBasis
  private :: FinTwoBodyRelChanSpinHOBasis
  private :: GetTwoBodyRelativeHORadial
  private :: GetTwoBodyRelativeHORadialFromNL
  private :: GetJ
  private :: GetParity
  private :: GetZ
  private :: GetS
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNmax

  type :: TwoBodyRelChanSpinHOBasis
    integer, private :: j, p, s, z, n_state
    integer, private :: lmin, lmax, Nmax
    integer, private, allocatable :: nl2idx(:,:)
    type(HarmonicOscillator), allocatable :: HOn(:)
  contains
    procedure :: InitTwoBodyRelChanSpinHOBasis
    procedure :: FinTwoBodyRelChanSpinHOBasis
    procedure :: GetTwoBodyRelativeHORadial
    procedure :: GetTwoBodyRelativeHORadialFromNL
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetS
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNmax

    ! overriding
    generic :: init => InitTwoBodyRelChanSpinHOBasis ! constructor
    generic :: fin => FinTwoBodyRelChanSpinHOBasis   ! destructor
    generic :: getp => GetTwoBodyRelativeHORadial, &
      & GetTwoBodyRelativeHORadialFromNL             ! getter of HOn (pointer returned)
  end type TwoBodyRelChanSpinHOBasis

  contains

  subroutine FinTwoBodyRelChanSpinHOBasis(this)
    ! two-body relative channel in harmonic oscillator basis (spin symmetry is assumed)
    ! destructor
    class(TwoBodyRelChanSpinHOBasis), intent(inout) :: this
    deallocate(this%nl2idx)
    deallocate(this%HOn)
  end subroutine FinTwoBodyRelChanSpinHOBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetS(this) result(s)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: s
    s = this%s
  end function GetS

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l) result(idx)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer, intent(in) :: n, l
    integer :: idx
    idx = this%nl2idx(n,l)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelChanSpinHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetTwoBodyRelativeHORadialFromNL(this, n, l) result(r)
    class(TwoBodyRelChanSpinHOBasis), target, intent(in) :: this
    integer, intent(in) :: n, l
    type(HarmonicOscillator), pointer :: r
    r => this%getp( this%nl2idx(n,l) )
  end function GetTwoBodyRelativeHORadialFromNL

  function GetTwoBodyRelativeHORadial(this, idx) result(r)
    class(TwoBodyRelChanSpinHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(HarmonicOscillator), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%HOn(idx)
  end function GetTwoBodyRelativeHORadial

  subroutine InitTwoBodyRelChanSpinHOBasis(this, Nmax, j, p, s, z)
    ! two-body relative channel in harmonic oscillator basis (spin symmetry is assumed)
    ! constructor
    class(TwoBodyRelChanSpinHOBasis), intent(inout) :: this
    integer, intent(in) :: j, p, s, z, Nmax
    integer :: lmin, lmax, n, l, num
    integer :: n12max, n12min
    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%s = s
    this%z = z
    n12min = 0
    n12max = (Nmax - abs(j-s)) / 2
    lmin = j+s
    lmax = abs(j-s)
    num = 0
    do l = abs(j-s), j+s
      if((-1) ** l /= p) cycle
      if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
      lmin = min(lmin, l)
      lmax = max(lmax, l)
      do n = n12min, n12max
        if(2 * n + l > Nmax) cycle
        num = num + 1
      end do
    end do
    this%n_state = num
    this%lmin = lmin
    this%lmax = lmax
    allocate(this%HOn(num))
    allocate(this%nl2idx(n12min:n12max, lmin:lmax))
    this%nl2idx(:,:) = 0

    num = 0
    do l = lmin, lmax
      if((-1) ** l /= p) cycle
      if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
      do n = n12min, n12max
        if(2 * n + l > Nmax) cycle
        num = num + 1
        if(num > this%n_state) stop 'error'
        call this%HOn(num)%set(n,l,s)
        this%nl2idx(n, l) = num
#ifdef TwoBodyRelChanSpinHOBasisDebug
        write(*,'(a,5i3)') 'InitTwoBodyRelChanSpinHOBasis, |n,l,s,j,tz>: ', &
            & n,l,s,j,z
#endif
      end do
    end do
  end subroutine InitTwoBodyRelChanSpinHOBasis

end module TwoBodyRelChanSpinHO
