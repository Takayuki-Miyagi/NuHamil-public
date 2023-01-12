module TwoBodyRelCMChanMeshBasis
  use TwoBodyRelCMRadQNumbers
  implicit none
  public :: TwoBodyRelCMChanMBasis
  private

  ! | plSJ_rel, PL: J >
  type :: TwoBodyRelCMChanMBasis
    integer, private :: j, p, z, j_rel, l_cm, n_states, NMesh_rel, NMesh_cm
    integer, private :: l_rel_min, l_rel_max, s_min, s_max
    integer, private, allocatable :: ils2idx(:,:,:,:) ! i_rel, i_cm, l_rel, spin
    type(RelCMMesh), allocatable :: points(:)
  contains
    procedure :: InitTwoBodyRelCMChanMBasis
    procedure :: FinTwoBodyRelCMChanMBasis
    procedure :: GetJ
    procedure :: GetJRel
    procedure :: GetLCM
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetNumberStates
    procedure :: GetIndex
    procedure :: GetNMeshRel
    procedure :: GetNMeshCM
    procedure :: GetTwoBodyRelCMMesh
    procedure :: GetTwoBodyRelCMMeshFromQNumbers
    procedure :: SetMeshWeight
    generic :: init => InitTwoBodyRelCMChanMBasis
    generic :: fin => FinTwoBodyRelCMChanMBasis
    generic :: GetRelCMMesh => GetTwoBodyRelCMMesh, GetTwoBodyRelCMMeshFromQNumbers
  end type TwoBodyRelCMChanMBasis
contains

  function GetJ(this) result(j)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetJRel(this) result(j_rel)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: j_rel
    j_rel = this%j_rel
  end function GetJRel

  function GetLCM(this) result(l_cm)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: l_cm
    l_cm = this%l_cm
  end function GetLCM

  function GetIndex(this,i_rel,l_rel,spin,i_cm) result(idx)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer, intent(in) :: i_rel, l_rel, spin, i_cm
    integer :: idx
    idx = this%ils2idx(i_rel,i_cm,l_rel,spin)
  end function GetIndex

  function GetNumberStates(this) result(Nstates)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: Nstates
    Nstates = this%n_states
  end function GetNumberStates

  function GetNMeshRel(this) result(NMesh)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh_rel
  end function GetNMeshRel

  function GetNMeshCM(this) result(NMesh)
    class(TwoBodyRelCMChanMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh_cm
  end function GetNMeshCM

  function GetTwoBodyRelCMMesh(this, idx) result(r)
    class(TwoBodyRelCMChanMBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(RelCMMesh), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%points(idx)
  end function GetTwoBodyRelCMMesh

  function GetTwoBodyRelCMMeshFromQNumbers(this, i_rel, j_rel, i_cm, l_cm) result(r)
    class(TwoBodyRelCMChanMBasis), target, intent(in) :: this
    integer, intent(in) :: i_rel, j_rel, i_cm, l_cm
    type(RelCMMesh), pointer :: r
    r => this%GetRelCMMesh( this%GetIndex(i_rel, j_rel, i_cm, l_cm) )
  end function GetTwoBodyRelCMMeshFromQNumbers

  subroutine FinTwoBodyRelCMChanMBasis(this)
    class(TwoBodyRelCMChanMBasis), intent(inout) :: this
    deallocate(this%ils2idx)
    deallocate(this%points)
  end subroutine FinTwoBodyRelCMChanMBasis

  subroutine InitTwoBodyRelCMChanMBasis(this, NMesh_rel, NMesh_cm, j, p, z, j_rel, l_cm)
    class(TwoBodyRelCMChanMBasis), intent(inout) :: this
    integer, intent(in) :: NMesh_rel, NMesh_cm, j, p, z, j_rel, l_cm
    integer :: l_rel, spin
    integer :: s_min, s_max, l_rel_min, l_rel_max
    integer :: cnt, i_rel, i_cm

    this%NMesh_rel = NMesh_rel
    this%NMesh_cm = NMesh_cm
    this%j = j
    this%p = p
    this%z = z
    this%j_rel = j_rel
    this%l_cm = l_cm
    s_min = 9999
    s_max =-9999
    l_rel_min = 9999
    l_rel_max =-9999
    cnt = 0
    do spin = 0, 1
      do l_rel = abs(j_rel-spin), j_rel+spin
        if((-1)**(l_rel+l_cm) /= p) cycle
        if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
        s_min = min(s_min, spin)
        s_max = max(s_max, spin)
        l_rel_min = min(l_rel, l_rel_min)
        l_rel_max = max(l_rel, l_rel_max)
        do i_rel = 1, NMesh_rel
          do i_cm = 1, NMesh_cm
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%n_states = cnt
    this%l_rel_min = l_rel_min
    this%l_rel_max = l_rel_max
    this%s_min = s_min
    this%s_max = s_max
    allocate(this%points(cnt))
    allocate(this%ils2idx(NMesh_rel, NMesh_cm, l_rel_min:l_rel_max, s_min:s_max))
    this%ils2idx(:,:,:,:) = 0

    cnt = 0
    do spin = 0, 1
      do l_rel = abs(j_rel-spin), j_rel+spin
        if((-1)**(l_rel+l_cm) /= p) cycle
        if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
        do i_rel = 1, NMesh_rel
          do i_cm = 1, NMesh_cm
            cnt = cnt + 1
            call this%points(cnt)%set(i_rel,l_rel,spin,i_cm,0.d0,0.d0,0.d0,0.d0)
            this%ils2idx(i_rel,i_cm,l_rel,spin) = cnt
            !write(*,'(a,8i6)') 'num, |i_rel,l_rel,spin,j_rel,i_cm,l_cm,tz>: ', &
            !    & cnt,i_rel,l_rel,spin,j_rel,i_cm,l_cm,z
          end do
        end do
      end do
    end do
  end subroutine InitTwoBodyRelCMChanMBasis

  subroutine SetMeshWeight(this, p_rel, w_rel, p_cm, w_cm)
    class(TwoBodyRelCMChanMBasis), intent(inout) :: this
    real(8), intent(in) :: p_rel(:), w_rel(:), p_cm(:), w_cm(:)
    integer :: i_rel, l_rel, spin, i_cm, num


    if( this%NMesh_rel /= size(p_rel) .or.  this%NMesh_rel /= size(w_rel)) then
      write(*,*) "Warning: NMesh is wrong"
      return
    end if

    if( this%NMesh_cm /= size(p_cm) .or.  this%NMesh_cm /= size(w_cm)) then
      write(*,*) "Warning: NMesh is wrong"
      return
    end if
    do spin = 0, 1
      do l_rel = abs(this%GetJRel()-spin), this%GetJRel()+spin
        if((-1)**(l_rel+this%GetLCM()) /= this%GetParity()) cycle
        if(abs(this%GetZ())==1 .and. (-1)**(l_rel+spin)==-1) cycle

        do i_rel = 1, this%GetNMeshRel()
          do i_cm = 1, this%GetNMeshCM()
            num = this%GetIndex(i_rel,l_rel,spin,i_cm)
            call this%points(num)%set(i_rel,l_rel,spin,i_cm,p_rel(i_rel),w_rel(i_rel),p_cm(i_cm),w_cm(i_cm))
          end do
        end do

      end do
    end do
  end subroutine SetMeshWeight
end module TwoBodyRelCMChanMeshBasis
