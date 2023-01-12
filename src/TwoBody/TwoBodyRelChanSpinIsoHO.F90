module TwoBodyRelChanSpinIsoHO
  use TwoBodyRelRadQNumbers
  implicit none

  public :: TwoBodyRelChanSpinIsoHOBasis
  private :: InitTwoBodyRelChanSpinIsoHOBasis
  private ::  FinTwoBodyRelChanSpinIsoHOBasis
  private :: GetTwoBodyRelativeHORadial
  private :: GetTwoBodyRelativeHORadialFromNL
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetS
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNmax

  type :: TwoBodyRelChanSpinIsoHOBasis
    integer, private :: j, p, s, t, n_state
    integer, private :: lmin, lmax, Nmax
    integer, private, allocatable :: nl2idx(:,:)
    type(HarmonicOscillator), allocatable :: HOn(:)
  contains
    procedure :: InitTwoBodyRelChanSpinIsoHOBasis
    procedure ::  FinTwoBodyRelChanSpinIsoHOBasis
    procedure :: GetTwoBodyRelativeHORadial
    procedure :: GetTwoBodyRelativeHORadialFromNL
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetS
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNmax

    ! overriding
    generic :: init => InitTwoBodyRelChanSpinIsoHOBasis ! constructor
    generic :: fin => FinTwoBodyRelChanSpinIsoHOBasis   ! destructor
    generic :: getp => GetTwoBodyRelativeHORadial                   ! getter of HOn (pointer returned)
  end type TwoBodyRelChanSpinIsoHOBasis
contains
  subroutine FinTwoBodyRelChanSpinIsoHOBasis(this)
    ! two-body relative channel in harmonic oscillator basis (spin symmetry is assumed)
    ! destructor
    class(TwoBodyRelChanSpinIsoHOBasis), intent(inout) :: this
    deallocate(this%nl2idx)
    deallocate(this%HOn)
  end subroutine FinTwoBodyRelChanSpinIsoHOBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(t)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetS(this) result(s)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: s
    s = this%s
  end function GetS

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l) result(idx)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer, intent(in) :: n, l
    integer :: idx
    idx = this%nl2idx(n,l)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetTwoBodyRelativeHORadialFromNL(this, n, l) result(r)
    class(TwoBodyRelChanSpinIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: n, l
    type(HarmonicOscillator), pointer :: r
    r => this%getp( this%nl2idx(n,l) )
  end function GetTwoBodyRelativeHORadialFromNL

  function GetTwoBodyRelativeHORadial(this, idx) result(r)
    class(TwoBodyRelChanSpinIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(HarmonicOscillator), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%HOn(idx)
  end function GetTwoBodyRelativeHORadial

  subroutine InitTwoBodyRelChanSpinIsoHOBasis(this, Nmax, j, p, s, t)
    ! two-body relative channel in harmonic oscillator basis (spin and isospin symmetries are assumed)
    ! constructor
    class(TwoBodyRelChanSpinIsoHOBasis), intent(inout) :: this
    integer, intent(in) :: j, p, s, t, Nmax
    integer :: lmin, lmax, n, l, num
    integer :: n12max, n12min
    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%s = s
    this%t = t
    lmin = j+s
    lmax = abs(j-s)
    n12min = 0
    n12max = (Nmax - abs(j-s)) / 2
    num = 0
    do l = abs(j-s), j+s
      if((-1) ** l /= p) cycle
      if((-1)**(l+s+t) == 1) cycle
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
    !this%nmax = n12max
    allocate(this%HOn(num))
    allocate(this%nl2idx(n12min:n12max, lmin:lmax))
    this%nl2idx(:,:) = 0

    num = 0
    do l = lmin, lmax
      if((-1) ** l /= p) cycle
      if((-1)**(l+s+t) == 1) cycle
      do n = n12min, n12max
        if(2 * n + l > Nmax) cycle
        num = num + 1
        if(num > this%n_state) stop 'error'
        call this%HOn(num)%set(n,l,s)
        this%nl2idx(n,l) = num
#ifdef TwoBodyRelChanSpinIsoHOBasisDebug
        write(*,'(a,5i3)') 'InitTwoBodyRelChanSpinIsoHOBasis, |n,l,s,j,t>: ', &
            & n,l,s,j,t
#endif
      end do
    end do
  end subroutine InitTwoBodyRelChanSpinIsoHOBasis
end module TwoBodyRelChanSpinIsoHO
