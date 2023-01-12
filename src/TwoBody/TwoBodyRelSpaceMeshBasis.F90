module TwoBodyRelSpaceMeshBasis
  use TwoBodyRelChanMeshBasis
  implicit none

  private :: InitTwoBodyRelSpaceMBasis
  private :: FinTwoBodyRelSpaceMBasis
  private :: GetTwoBodyRelChanMBasisFromJPZ
  private :: GetTwoBodyRelChanMBasisFromCh
  private :: SetMeshWeight
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetNMesh
  private :: GetJmax

  type :: TwoBodyRelSpaceMBasis
    type(TwoBodyRelChanMBasis), allocatable :: jpz(:)
    integer, private :: NMesh, Jmax, NChan
    integer, private, allocatable :: jpz2idx(:,:,:)
    real(8), private :: xmin = 0.d0, xmax = 0.d0
    logical, private :: initialized = .false.
  contains
    procedure :: InitTwoBodyRelSpaceMBasis
    procedure :: FinTwoBodyRelSpaceMBasis
    procedure :: GetTwoBodyRelChanMBasisFromJPZ
    procedure :: GetTwoBodyRelChanMBasisFromCh
    procedure :: SetMeshWeight
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetNMesh
    procedure :: GetJmax
    generic :: init => InitTwoBodyRelSpaceMBasis
    generic :: fin => FinTwoBodyRelSpaceMBasis
    generic :: GetChannel => GetTwoBodyRelChanMBasisFromJPZ, GetTwoBodyRelChanMBasisFromCh
    generic :: setp => SetMeshWeight
  end type TwoBodyRelSpaceMBasis

contains

  subroutine FinTwoBodyRelSpaceMBasis(this)
    class(TwoBodyRelSpaceMBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz(ich)%fin()
    end do
    deallocate(this%jpz)
    deallocate(this%jpz2idx)
    this%xmin = 0.d0
    this%xmax = 0.d0
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceMBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceMBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, z) result(idx)
    class(TwoBodyRelSpaceMBasis), intent(in) :: this
    integer, intent(in) :: j, p, z
    integer :: idx
    idx = this%jpz2idx(j, p, z)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelSpaceMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceMBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceMBasis(this, xmin, xmax, NMesh, Jmax, &
        & j_in, p_in, z_in)
    class(TwoBodyRelSpaceMBasis), intent(inout) :: this
    integer, intent(in) :: Nmesh, Jmax
    real(8), intent(in) :: xmin, xmax
    integer, intent(in), optional :: j_in, p_in, z_in
    integer :: j, p, z, s, n, ich, cnt, loop, l

    this%initialized = .true.
    this%NMesh = NMesh
    this%Jmax = Jmax
    this%xmin = xmin
    this%xmax = xmax
    do loop = 1, 2
      ich = 0
      do z = -1, 1
        do j = 0, this%jmax
          do p = 1, -1, -2
            if(present(j_in)) then
              if(j /= j_in) cycle
            end if
            if(present(p_in)) then
              if(p /= p_in) cycle
            end if
            if(present(z_in)) then
              if(z /= z_in) cycle
            end if

            cnt = 0
            do s = 0, 1
              do l = iabs(j-s), j+s
                if((-1) ** l /= p) cycle
                if(l < 0) cycle
                do n = 0, this%NMesh
                  if(iabs(z) == 1 .and. (-1) **(l + s) == -1) cycle
                  cnt = cnt + 1
                end do
              end do
            end do
#ifdef TwoBodyRelativeSpaceDebug
            if(cnt /= 0) then
              write(*,'(a,4i4)') 'InitTwoBodyRelSpaceHOBasis |Parity,J,Tz>', &
                  & p,j,z
            end if
#endif
            if(cnt > 0) then
              ich = ich + 1
              if(loop == 2) then
                this%jpz2idx(j, p, z) = ich
                call this%jpz(ich)%init(this%NMesh, j, p, z)
              end if
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpz(ich))
        allocate(this%jpz2idx(0:jmax, -1:1, -1:1))
        this%jpz2idx(:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(z_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpace, spin, parity, and pn are not correct'
      stop
    end if

  end subroutine InitTwoBodyRelSpaceMBasis

  subroutine SetMeshWeight(this, x, w)
    class(TwoBodyRelSpaceMBasis), intent(inout) :: this
    real(8), intent(in) :: x(:), w(:)
    integer :: ch
    do ch = 1, this%NChan
      call this%jpz(ch)%SetMeshWeight(x,w)
    end do
  end subroutine SetMeshWeight

  function GetTwoBodyRelChanMBasisFromJPZ(this, j, p, z) result(r)
    class(TwoBodyRelSpaceMBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, z
    type(TwoBodyRelChanMBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,z) )
  end function GetTwoBodyRelChanMBasisFromJPZ

  function GetTwoBodyRelChanMBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceMBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanMBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpz(ch)
  end function GetTwoBodyRelChanMBasisFromCh

end module TwoBodyRelSpaceMeshBasis

!program main
!  use TwoBodyRelSpaceMeshBasis
!  type(TwoBodyRelSpaceMBasis) :: ms
!  type(TwoBodyRelChanMBasis), pointer :: ch
!  type(Mesh), pointer :: point
!  integer :: ich, i_st
!  integer :: j, p, z
!  call ms%init(0.d0, 10.d0, 5, 1)
!  do ich = 1, ms%NChan
!    ch => ms%getp(ich)
!    j = ch%j
!    p = ch%p
!    z = ch%z
!    write(*,"(3i3)") j, p, z
!    do i_st = 1, ch%n_state
!      point => ch%getp(i_st)
!      write(*,"(3i3,2f12.6)") point%GetI(), point%GetL(), point%GetS(), &
!        & point%GetP(), point%GetW()
!    end do
!  end do
!end program main
