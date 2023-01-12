module TwoBodyRelChanSpinMeshBasis
  use TwoBodyRelRadQNumbers, only: Mesh
  implicit none

  public :: TwoBodyRelChanSpinMBasis
  private :: InitTwoBodyRelChanSpinMBasis
  private :: FinTwoBodyRelChanSpinMBasis
  private :: GetTwoBodyRelativeMeshRadialFromNL
  private :: GetTwoBodyRelativeMeshRadial
  private :: SetMeshWeight
  private :: GetJ
  private :: GetParity
  private :: GetZ
  private :: GetS
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNMesh

  type :: TwoBodyRelChanSpinMBasis
    ! Relative Momentum basis for given J,P,S,Tz
    ! used only for calc. of two-body interaction
    integer, private :: j, p, s, z, n_state
    integer, private :: lmin, lmax, NMesh
    real(8), private :: pmin, pmax
    integer, private, allocatable :: nl2idx(:,:)
    type(Mesh), allocatable :: points(:)
  contains
    procedure :: InitTwoBodyRelChanSpinMBasis
    procedure :: FinTwoBodyRelChanSpinMBasis
    procedure :: GetTwoBodyRelativeMeshRadialFromNL
    procedure :: GetTwoBodyRelativeMeshRadial
    procedure :: SetMeshWeight
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetS
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNMesh

    ! overriding
    generic :: init => InitTwoBodyRelChanSpinMBasis ! constructor
    generic :: fin => FinTwoBodyRelChanSpinMBasis   ! destructor
    generic :: getp => GetTwoBodyRelativeMeshRadial, GetTwoBodyRelativeMeshRadialFromNL
  end type TwoBodyRelChanSpinMBasis

  contains

  subroutine FinTwoBodyRelChanSpinMBasis(this)
    ! two-body relative channel in harmonic oscillator basis (spin symmetry is assumed)
    ! destructor
    class(TwoBodyRelChanSpinMBasis), intent(inout) :: this
    deallocate(this%nl2idx)
    deallocate(this%points)
  end subroutine FinTwoBodyRelChanSpinMBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetS(this) result(s)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: s
    s = this%s
  end function GetS

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l) result(idx)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer, intent(in) :: n, l
    integer :: idx
    idx = this%nl2idx(n,l)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelChanSpinMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetTwoBodyRelativeMeshRadialFromNL(this, n, l) result(r)
    class(TwoBodyRelChanSpinMBasis), target, intent(in) :: this
    integer, intent(in) :: n, l
    type(Mesh), pointer :: r
    integer :: idx
    r => null()
    idx = this%nl2idx(n,l)
    r => this%getp(idx)
  end function GetTwoBodyRelativeMeshRadialFromNL

  function GetTwoBodyRelativeMeshRadial(this, idx) result(r)
    class(TwoBodyRelChanSpinMBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(Mesh), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%points(idx)
  end function GetTwoBodyRelativeMeshRadial

  subroutine SetMeshWeight(this, p, w)
    class(TwoBodyRelChanSpinMBasis), intent(inout) :: this
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
      if(iabs(this%z) == 1 .and. (-1)**(l+this%s) == -1) cycle
      do i = 1, this%NMesh
        num = this%nl2idx(i,l)
        call this%points(num)%set(i,l,this%s,p(i),w(i))
      end do
    end do
  end subroutine SetMeshWeight

  subroutine InitTwoBodyRelChanSpinMBasis(this, NMesh, j, p, s, z)
    ! two-body relative channel in Meshentum basis (spin symmetry is assumed)
    ! constructor
    class(TwoBodyRelChanSpinMBasis), intent(inout) :: this
    integer, intent(in) :: NMesh
    integer, intent(in) :: j, p, s, z
    integer :: lmin, lmax, n, i, l, num

    this%NMesh = NMesh
    this%j = j
    this%p = p
    this%s = s
    this%z = z
    lmin = j+s
    lmax = abs(j-s)
    n = 0
    do l = abs(j-s), j+s
      if((-1) ** l /= p) cycle
      if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
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
      if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
      do i = 1, NMesh
        num = num + 1
        call this%points(num)%set(i,l,s,0.d0,0.d0)
        this%nl2idx(i,l) = num
#ifdef TwoBodyRelChanSpinMBasisDebug
        write(*,'(a,6i4)') 'InitTwoBodyRelChanSpinMBasis, num, |i,l,s,j,tz>: ', &
            & num,i,l,s,j,z
#endif
      end do
    end do
  end subroutine InitTwoBodyRelChanSpinMBasis
end module TwoBodyRelChanSpinMeshBasis

!program test
!  use RelativeChannelSpinMeshentum
!  implicit none
!  type(TwoBodyRelativeChannelSpinMeshBasis) :: mom
!  type(Meshentum), pointer :: n
!  integer :: i
!  call Mesh%init(0.d0, 5.d0, 20, 1, 1, 1, 0)
!  do i = 1, Mesh%n_state
!    n => Mesh%getp(i)
!    write(*,'(3i3,2f12.6)') n%n, n%l, n%s, n%p, n%w
!  end do
!  call Mesh%fin()
!end program test
