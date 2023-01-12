module TwoBodyRelChanMeshBasis
  use TwoBodyRelRadQNumbers, only: Mesh
  implicit none

  private :: InitTwoBodyRelChanMBasis
  private :: FinTwoBodyRelChanMBasis
  private :: GetTwoBodyRelativeMeshRadialFromILS
  private :: GetTwoBodyRelativeMeshRadial
  private :: SetMeshWeight
  private :: GetJ
  private :: GetParity
  private :: GetZ
  private :: GetNumberStates
  private :: GetIndex
  private :: GetNMesh

  type :: TwoBodyRelChanMBasis
    integer, private :: j, p, z, n_state
    integer, private :: lmin, lmax, smin, smax, NMesh
    real(8), private :: pmin, pmax
    integer, private, allocatable :: ils2idx(:,:,:)
    type(Mesh), allocatable :: points(:)
  contains
    procedure :: InitTwoBodyRelChanMBasis
    procedure :: FinTwoBodyRelChanMBasis
    procedure :: GetTwoBodyRelativeMeshRadialFromILS
    procedure :: GetTwoBodyRelativeMeshRadial
    procedure :: SetMeshWeight
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNMesh
    generic :: init => InitTwoBodyRelChanMBasis
    generic :: fin => FinTwoBodyRelChanMBasis
    generic :: getp => GetTwoBodyRelativeMeshRadialFromILS, &
        & GetTwoBodyRelativeMeshRadial
    generic :: setp => SetMeshWeight
  end type TwoBodyRelChanMBasis
contains
  subroutine FinTwoBodyRelChanMBasis(this)
    class(TwoBodyRelChanMBasis), intent(inout) :: this
    deallocate(this%ils2idx)
    deallocate(this%points)
  end subroutine FinTwoBodyRelChanMBasis

  function GetJ(this) result(j)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer :: z
    Z = this%z
  end function GetZ

  function GetNumberStates(this) result(n)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer :: n
    n = this%n_state
  end function GetNumberStates

  function GetIndex(this,n,l,s) result(idx)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer, intent(in) :: n, l, s
    integer :: idx
    idx = this%ils2idx(n,l,s)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelChanMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetTwoBodyRelativeMeshRadialFromILS(this, i, l, s) result(r)
    class(TwoBodyRelChanMBasis), target, intent(in) :: this
    integer, intent(in) :: i, l, s
    type(Mesh), pointer :: r
    r => this%getp( this%ils2idx(i,l,s) )
  end function GetTwoBodyRelativeMeshRadialFromILS

  function GetTwoBodyRelativeMeshRadial(this, idx) result(r)
    class(TwoBodyRelChanMBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(Mesh), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%points(idx)
  end function GetTwoBodyRelativeMeshRadial

  subroutine SetMeshWeight(this, p, w)
    class(TwoBodyRelChanMBasis), intent(inout) :: this
    real(8), intent(in) :: p(:), w(:)
    integer :: l, i, s, num

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
    do s = this%smin, this%smax
      do l = this%lmin, this%lmax
        if((-1) ** l /= this%p) cycle
        if(iabs(this%z) == 1 .and. (-1)**(l+s) == -1) cycle
        do i = 1, this%NMesh
          num = this%ils2idx(i,l,s)
          call this%points(num)%set(i,l,s,p(i),w(i))
        end do
      end do
    end do
  end subroutine SetMeshWeight

  subroutine InitTwoBodyRelChanMBasis(this, NMesh, j, p, z)
    ! two-body relative channel in Meshentum basis (spin symmetry is assumed)
    ! constructor
    class(TwoBodyRelChanMBasis), intent(inout) :: this
    integer, intent(in) :: NMesh
    integer, intent(in) :: j, p, z
    integer :: lmin, lmax, n, i, l, s, num, smin, smax

    this%NMesh = NMesh
    this%j = j
    this%p = p
    this%z = z
    smin = 1
    smax = 0
    lmin = j+1
    lmax = abs(j-1)
    n = 0
    do s = 0, 1
      do l = abs(j-s), j+s
        if((-1) ** l /= p) cycle
        if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
        smin = min(smin, s)
        smax = max(smax, s)
        lmin = min(lmin, l)
        lmax = max(lmax, l)
        do i = 1, NMesh
          n = n + 1
        end do
      end do
    end do
    this%n_state = n
    this%lmin = lmin
    this%lmax = lmax
    this%smin = smin
    this%smax = smax
    allocate(this%points(n))
    allocate(this%ils2idx(NMesh, lmin:lmax, smin:smax))
    this%ils2idx(:,:,:) = 0

    num = 0
    do s = smin, smax
      do l = lmin, lmax
        if((-1) ** l /= p) cycle
        if(iabs(z) == 1 .and. (-1)**(l+s) == -1) cycle
        do i = 1, NMesh
          num = num + 1
          call this%points(num)%set(i,l,s,0.d0,0.d0)
          this%ils2idx(i,l,s) = num
#ifdef TwoBodyRelChanSpinMBasisDebug
          write(*,'(a,6i4)') 'InitTwoBodyRelChanSpinMBasis, num, |i,l,s,j,tz>: ', &
              & num,i,l,s,j,z
#endif
        end do
      end do
    end do
  end subroutine InitTwoBodyRelChanMBasis
end module TwoBodyRelChanMeshBasis

!program test
!  use TwoBodyRelChanMeshBasis
!  implicit none
!  type(TwoBodyRelChanMBasis) :: mom
!  type(Mesh), pointer :: n
!  integer :: i
!  call mom%init(5, 1, 1, 0)
!  call mom%setp([1.d0,2.d0,3.d0,4.d0,5.d0],[0.1d0,0.2d0,4.d0,5.d0,1.d0])
!  do i = 1, mom%n_state
!    n => mom%getp(i)
!    write(*,'(3i3,2f12.6)') n%n, n%l, n%s, n%p, n%w
!  end do
!  call mom%fin()
!end program test
