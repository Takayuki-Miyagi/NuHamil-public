module TwoBodyRelChanIsoHO
  use TwoBodyRelRadQNumbers
  implicit none

  public :: TwoBodyRelChanIsoHOBasis

  private :: InitTwoBodyRelChanIsoHOBasis
  private :: FinTwoBodyRelChanIsoHOBasis
  private :: GetTwoBodyRelativeHORadial
  private :: GetTwoBodyRelativeHORadialFromNLS
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNmax
  private :: GetFrequency

  type :: TwoBodyRelChanIsoHOBasis
    integer, private :: j, p, t, n_state
    integer, private :: lmin, lmax, smin, smax, Nmax
    real(8) :: hw = -1.d0
    integer, private, allocatable :: nls2idx(:,:,:)
    type(HarmonicOscillator), allocatable :: HOn(:)
  contains
    procedure :: InitTwoBodyRelChanIsoHOBasis
    procedure :: FinTwoBodyRelChanIsoHOBasis
    procedure :: GetTwoBodyRelativeHORadial
    procedure :: GetTwoBodyRelativeHORadialFromNLS
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNmax
    procedure :: GetFrequency

    ! overriding
    generic :: init => InitTwoBodyRelChanIsoHOBasis ! constructor
    generic :: fin => FinTwoBodyRelChanIsoHOBasis   ! destructor
    generic :: getp => GetTwoBodyRelativeHORadial, &
      & GetTwoBodyRelativeHORadialFromNLS           ! getter of HOn (pointer returned)
  end type TwoBodyRelChanIsoHOBasis

contains

  subroutine FinTwoBodyRelChanIsoHOBasis(this)
    ! two-body relative channel in harmonic oscilator basis
    ! destructor
    class(TwoBodyRelChanIsoHOBasis), intent(inout) :: this
    if( this%n_state == 0 ) return
    deallocate(this%nls2idx)
    deallocate(this%HOn)
  end subroutine FinTwoBodyRelChanIsoHOBasis

  function GetJ(this) result(j)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(t)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetFrequency(this) result(hw)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetNumberStates(this) result(n)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l,s) result(idx)
    class(TwoBodyRelCHanIsoHOBasis), intent(in) :: this
    integer, intent(in) :: n, l, s
    integer :: idx
    !write(*,'(7i4)') n, l, s, this%GetNmax(), this%GetJ(), this%GetParity(), this%GetT()
    idx = this%nls2idx(n,l,s)
  end function GetIndex

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelChanIsoHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetTwoBodyRelativeHORadialFromNLS(this, n, l, s) result(r)
    class(TwoBodyRelChanIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: n, l, s
    type(HarmonicOscillator), pointer :: r
    r => this%getp( this%nls2idx(n, l, s) )
  end function GetTwoBodyRelativeHORadialFromNLS

  function GetTwoBodyRelativeHORadial(this, idx) result(r)
    class(TwoBodyRelChanIsoHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(HarmonicOscillator), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%HOn(idx)
  end function GetTwoBodyRelativeHORadial

  subroutine InitTwoBodyRelChanIsoHOBasis(this, Nmax, j, p, t, hw)
    ! two-body relative channel in harmonic oscilator basis
    ! constructor
    class(TwoBodyRelChanIsoHOBasis), intent(inout) :: this
    integer, intent(in) :: Nmax, j, p, t
    real(8), intent(in) :: hw
    integer :: lmin, lmax
    integer :: smin, smax
    integer :: n12max, n12min
    integer :: n, l, cnt, s
    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%t= t
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
          if((-1) **(l + s + t) == 1) cycle
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
    allocate(this%nls2idx(n12min:n12max, lmin:lmax, smin:smax))
    this%nls2idx(:,:,:) = 0

    cnt = 0
    do s = 0, 1
      do l = iabs(j-s), j+s
        if((-1) ** l /= p) cycle
        if(l < 0) cycle
        do n = 0, Nmax
          if(2*n + l > Nmax) cycle
          if((-1) **(l + s + t) == 1) cycle
          cnt = cnt + 1
          call this%HOn(cnt)%set(n,l,s)
          this%nls2idx(n,l,s) = cnt
#ifdef TwoBodyRelChanIsoHOBasisDebug
          write(*,'(a,5i3)') 'InitTwoBodyRelChanIsoHOBasis, |n,l,s,j,t>: ', &
              & n,l,s,j,t
#endif
        end do
      end do
    end do
  end subroutine InitTwoBodyRelChanIsoHOBasis
end module TwoBodyRelChanIsoHO
