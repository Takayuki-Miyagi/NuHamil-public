module TwoBodyRelChanSpinIsoMeshBasis
  use TwoBodyRelRadQNumbers
  implicit none

  public :: TwoBodyRelChanSpinIsoMBasis
  private :: InitTwoBodyRelChanSpinIsoMBasis
  private :: FinTwoBodyRelChanSpinIsoMBasis
  private :: GetTwoBodyRelativeMeshRadialFromNL
  private :: GetTwoBodyRelativeMeshRadial
  private :: SetMeshWeight
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetS
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNMesh

  type :: TwoBodyRelChanSpinIsoMBasis
    ! Relative Meshentum basis for given J,P,S,Tz
    ! used only for calc. of two-body interaction
    integer, private :: j, p, s, t, n_state
    integer, private :: lmin, lmax, NMesh
    real(8), private :: pmin, pmax
    integer, private, allocatable :: nl2idx(:,:)
    type(Mesh), allocatable :: points(:)
  contains
    procedure :: InitTwoBodyRelChanSpinIsoMBasis
    procedure :: FinTwoBodyRelChanSpinIsoMBasis
    procedure :: GetTwoBodyRelativeMeshRadialFromNL
    procedure :: GetTwoBodyRelativeMeshRadial
    procedure :: SetMeshWeight
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetS
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNMesh

    ! overriding
    generic :: init => InitTwoBodyRelChanSpinIsoMBasis ! constructor
    generic :: fin => FinTwoBodyRelChanSpinIsoMBasis   ! destructor
    generic :: getp => GetTwoBodyRelativeMeshRadial, GetTwoBodyRelativeMeshRadialFromNL
  end type TwoBodyRelChanSpinIsoMBasis

  contains

  subroutine FinTwoBodyRelChanSpinIsoMBasis(this)
    ! two-body relative channel in harmonic oscillator basis (SpinIsospin symmetry is assumed)
    ! destructor
    class(TwoBodyRelChanSpinIsoMBasis), intent(inout) :: this
    deallocate(this%nl2idx)
    deallocate(this%points)
  end subroutine FinTwoBodyRelChanSpinIsoMBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(t)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetS(this) result(s)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: s
    s = this%s
  end function GetS

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l) result(idx)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer, intent(in) :: n, l
    integer :: idx
    idx = this%nl2idx(n,l)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelChanSpinIsoMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetTwoBodyRelativeMeshRadialFromNL(this, n, l) result(r)
    class(TwoBodyRelChanSpinIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: n, l
    type(Mesh), pointer :: r
    integer :: idx
    r => null()
    idx = this%nl2idx(n,l)
    if(idx == 0) return
    r => this%getp(idx)
  end function GetTwoBodyRelativeMeshRadialFromNL

  function GetTwoBodyRelativeMeshRadial(this, idx) result(r)
    class(TwoBodyRelChanSpinIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(Mesh), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%points(idx)
  end function GetTwoBodyRelativeMeshRadial

  subroutine SetMeshWeight(this, p, w)
    class(TwoBodyRelChanSpinIsoMBasis), intent(inout) :: this
    real(8), intent(in) :: p(:), w(:)
    integer :: l, i, num

    if(this%NMesh /= size(p)) then
      write(*,*) "Warning: NMesh is wrong"
      return
    end if

    if(this%NMesh /= size(w)) then
      write(*,*) "Warning: NMesh is wrong"
      return
    end if

    this%pmin = minval(p)
    this%pmax = maxval(p)
    do l = this%lmin, this%lmax
      if((-1) ** l /= this%p) cycle
      if((-1)**(l+this%s+this%t) == 1) cycle
      do i = 1, this%NMesh
        num = this%nl2idx(i,l)
        call this%points(num)%set(i,l,this%s,p(i),w(i))
      end do
    end do
  end subroutine SetMeshWeight

  subroutine InitTwoBodyRelChanSpinIsoMBasis(this, pmin, pmax, NMesh, j, p, s, t)
    ! two-body relative channel in Meshentum basis (SpinIsospin symmetry is assumed)
    ! constructor
    use MyLibrary, only: gauss_legendre
    class(TwoBodyRelChanSpinIsoMBasis), intent(inout) :: this
    real(8), intent(in) :: pmin, pmax
    integer, intent(in) :: NMesh
    integer, intent(in) :: j, p, s, t
    integer :: lmin, lmax, n, i, l, num
    real(8), allocatable :: pMesh(:), wmom(:)

    this%pmin = pmin
    this%pmax = pmax
    this%NMesh = NMesh
    call gauss_legendre(pmin,pmax,pMesh,wmom,Nmesh)
    this%j = j
    this%p = p
    this%s = s
    this%t = t
    lmin = j+s
    lmax = abs(j-s)
    n = 0
    do l = abs(j-s), j+s
      if((-1) ** l /= p) cycle
      if((-1)**(l+s+t) == 1) cycle
      lmin = min(lmin, l)
      lmax = max(lmax, l)
      do i = 1, nmesh
        n = n + 1
      end do
    end do
    this%n_state = n
    this%lmin = lmin
    this%lmax = lmax
    allocate(this%points(n))
    allocate(this%nl2idx(NMesh, lmin:lmax))
    this%nl2idx(:,:) = 0

    num = 0
    do l = lmin, lmax
      if((-1) ** l /= p) cycle
      if((-1)**(l+s+t) == 1) cycle
      do i = 1, NMesh
        num = num + 1
        call this%points(num)%set(i,l,s,pMesh(i),wmom(i))
        this%nl2idx(i,l) = num
#ifdef TwoBodyRelChanSpinIsoMBasisDebug
        write(*,'(a,6i4)') 'InitTwoBodyRelChanSpinIsoMBasis, num, |i,l,s,j,tz>: ', &
            & num,i,l,s,j,z
#endif
      end do
    end do
    deallocate(pMesh,wmom)
  end subroutine InitTwoBodyRelChanSpinIsoMBasis

end module TwoBodyRelChanSpinIsoMeshBasis
