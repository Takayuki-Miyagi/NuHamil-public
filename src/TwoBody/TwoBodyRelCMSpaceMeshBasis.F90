module TwoBodyRelCMSpaceMeshBasis
  use TwoBodyRelCMChanMeshBasis
  implicit none
  public :: TwoBodyRelCMSpaceMBasis
  public :: TestTwoBodyRelCMSpaceMBasis
  private

  type :: TwoBodyRelCMSpaceMBasis
    type(TwoBodyRelCMChanMBasis), allocatable :: jpz_jrel_lcm(:)
    integer, private :: NMesh_rel, NMesh_cm, Jmax, NChan, j_rel_max, Lcm_max
    integer, private, allocatable :: jpz_jrel_lcm2idx(:,:,:,:,:)
    real(8), private :: xmin_rel = 0.d0, xmax_rel = 0.d0
    real(8), private :: xmin_cm = 0.d0, xmax_cm = 0.d0
    logical, private :: initialized = .false.
    real(8), allocatable :: x_rel(:), w_rel(:), x_cm(:), w_cm(:)
  contains
    procedure :: InitTwoBodyRelCMSpaceMBasis
    procedure :: FinTwoBodyRelCMSpaceMBasis
    procedure :: GetJmax
    procedure :: GetNMeshRel
    procedure :: GetNMeshCM
    procedure :: GetNumberChannels
    procedure :: GetJRelMax
    procedure :: GetLcmMax
    procedure :: GetIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromQNumbers
    procedure :: SetMeshWeight
    procedure :: GetPrelMin
    procedure :: GetPrelMax
    procedure :: GetPcmMin
    procedure :: GetPcmMax
    generic :: init => InitTwoBodyRelCMSpaceMBasis
    generic :: fin => FinTwoBodyRelCMSpaceMBasis
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromQNumbers
  end type TwoBodyRelCMSpaceMBasis

contains

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  function GetNMeshRel(this) result(NMesh_rel)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: NMesh_rel
    NMesh_rel = this%NMesh_rel
  end function GetNMeshRel

  function GetNMeshCM(this) result(NMesh_cm)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: NMesh_cm
    NMesh_cm = this%NMesh_cm
  end function GetNMeshCM

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetJRelMax(this) result(j_rel_max)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: j_rel_max
    j_rel_max = this%j_rel_max
  end function GetJrelMax

  function GetLcmMax(this) result(Lcm_max)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer :: Lcm_max
    Lcm_max = this%Lcm_max
  end function GetLcmMax

  function GetPrelMin(this) result(p)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    real(8) :: p
    p = this%xmin_rel
  end function GetPrelMin

  function GetPrelMax(this) result(p)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    real(8) :: p
    p = this%xmax_rel
  end function GetPrelMax

  function GetPcmMin(this) result(p)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    real(8) :: p
    p = this%xmin_cm
  end function GetPcmMin

  function GetPcmMax(this) result(p)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    real(8) :: p
    p = this%xmax_cm
  end function GetPcmMax

  function GetIndex(this, j, p, z, jrel, lcm) result(idx)
    class(TwoBodyRelCMSpaceMBasis), intent(in) :: this
    integer, intent(in) :: j, p, z, jrel, lcm
    integer :: idx
    idx = this%jpz_jrel_lcm2idx(j,p,z,jrel,lcm)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(r)
    class(TwoBodyRelCMSpaceMBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(TwoBodyRelCMChanMBasis), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%jpz_jrel_lcm(idx)
  end function GetChannelFromIndex

  function GetChannelFromQNumbers(this, j, p, z, jrel, lcm) result(r)
    class(TwoBodyRelCMSpaceMBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, z, jrel, lcm
    type(TwoBodyRelCMChanMBasis), pointer :: r
    r => this%GetChannel(this%GetIndex(j,p,z,jrel,lcm))
  end function GetChannelFromQNumbers

  subroutine FinTwoBodyRelCMSpaceMBasis(this)
    class(TwoBodyRelCMSpaceMBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz_jrel_lcm(ich)%fin()
    end do
    deallocate(this%jpz_jrel_lcm)
    deallocate(this%jpz_jrel_lcm2idx)
    if(allocated(this%x_rel)) deallocate(this%x_rel)
    if(allocated(this%w_rel)) deallocate(this%w_rel)
    if(allocated(this%x_cm )) deallocate(this%x_cm)
    if(allocated(this%w_cm )) deallocate(this%w_cm)
    this%xmin_rel = 0.d0
    this%xmax_rel = 0.d0
    this%xmin_cm = 0.d0
    this%xmax_cm = 0.d0
    this%initialized = .false.
  end subroutine FinTwoBodyRelCMSpaceMBasis

  subroutine InitTwoBodyRelCMSpaceMBasis(this, xmin_rel, xmax_rel, NMesh_rel, &
        & xmin_cm, xmax_cm, NMesh_cm, Jmax, j_rel_max, Lcm_max)
    class(TwoBodyRelCMSpaceMBasis), intent(inout) :: this
    integer, intent(in) :: NMesh_rel, NMesh_cm, Jmax, j_rel_max
    integer, intent(in), optional :: Lcm_max
    real(8), intent(in) :: xmin_rel, xmax_rel, xmin_cm, xmax_cm
    integer :: j, p, z, j_rel, l_cm, spin, l_rel, cnt, nch, loop

    this%initialized = .true.
    this%xmin_rel = xmin_rel
    this%xmax_rel = xmax_rel
    this%NMesh_rel = NMesh_rel
    this%xmin_cm = xmin_cm
    this%xmax_cm = xmax_cm
    this%NMesh_cm = NMesh_cm
    this%Jmax = Jmax
    this%J_rel_max = j_rel_max
    this%Lcm_max = Jmax+j_rel_max
    if(present(Lcm_max)) this%Lcm_max = Lcm_max

    do loop = 1, 2
      nch = 0
      do z = -1, 1
        do j = 0, this%GetJmax()
          do p = 1, -1, -2

            do j_rel = 0, this%GetJRelMax()
              do l_cm = abs(j-j_rel), min(j+j_rel, this%GetLcmMax())

                cnt = 0
                do spin = 0, 1
                  do l_rel = abs(j_rel-spin), j_rel+spin
                    if((-1)**(l_rel+l_cm) /= p) cycle
                    if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
                    cnt = cnt + 1
                  end do
                end do
                if(cnt > 0) then
                  nch = nch + 1
                  if(loop==2) then
                    this%jpz_jrel_lcm2idx(j,p,z,j_rel,l_cm) = nch
                    call this%jpz_jrel_lcm(nch)%init(NMesh_rel, NMesh_cm, j, p, z, j_rel, l_cm)
                  end if
                end if

              end do
            end do

          end do
        end do
      end do
      if(loop==1) then
        this%NChan = nch
        allocate(this%jpz_jrel_lcm(nch))
        allocate(this%jpz_jrel_lcm2idx(0:Jmax,-1:1,-1:1,0:j_rel_max,0:this%Lcm_max))
        this%jpz_jrel_lcm2idx(:,:,:,:,:) = 0
      end if
    end do

  end subroutine InitTwoBodyRelCMSpaceMBasis

  subroutine SetMeshWeight(this, x_rel, w_rel, x_cm, w_cm)
    class(TwoBodyRelCMSpaceMBasis), intent(inout) :: this
    real(8), intent(in) :: x_rel(:), w_rel(:), x_cm(:), w_cm(:)
    integer :: ch
    allocate(this%x_rel, source=x_rel)
    allocate(this%w_rel, source=w_rel)
    allocate(this%x_cm, source=x_cm)
    allocate(this%w_cm, source=w_cm)
    do ch = 1, this%NChan
      call this%jpz_jrel_lcm(ch)%SetMeshWeight(x_rel,w_rel,x_cm,w_cm)
    end do
  end subroutine SetMeshWeight

  subroutine TestTwoBodyRelCMSpaceMBasis()
    ! | (Pcm, Lcm), (prel, lrel, s, jrel) J >
    type(TwoBodyRelCMSpaceMBasis) :: ms
    type(TwoBodyRelCMChanMBasis), pointer :: chan
    integer :: ich
    call ms%init(0.d0, 10.d0, 50, 0.d0, 8.d0, 50, 8, 8, 2)
    do ich = 1, ms%GetNumberChannels()
      chan => ms%GetChannel(ich)
      write(*,'(a,5i4,i8)') 'Index, Lcm, Jrel, J, Z, Num: ', &
          & ich, chan%GetLCM(), chan%GetJRel(), chan%GetJ(), chan%GetZ(), chan%GetNumberStates()
    end do

  end subroutine TestTwoBodyRelCMSpaceMBasis
end module TwoBodyRelCMSpaceMeshBasis
