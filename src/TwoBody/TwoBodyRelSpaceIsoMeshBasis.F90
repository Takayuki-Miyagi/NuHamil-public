module TwoBodyRelSpaceIsoMeshBasis
  use TwoBodyRelChanIsoMeshBasis
  implicit none

  private :: InitTwoBodyRelSpaceIsoMBasis
  private :: FinTwoBodyRelSpaceIsoMBasis
  private :: GetTwoBodyRelChanIsoMBasisFromJPZ
  private :: GetTwoBodyRelChanIsoMBasisFromCh
  private :: SetMeshWeight
  private :: GetNumberChannels
  private :: GetIndex
  private :: GetNMesh
  private :: GetJmax

  type :: TwoBodyRelSpaceIsoMBasis
    type(TwoBodyRelChanIsoMBasis), allocatable :: jpt(:)
    integer :: NMesh, Jmax, NChan
    integer, allocatable :: jpt2idx(:,:,:)
    real(8) :: xmin = 0.d0, xmax = 0.d0
    logical :: initialized = .false.
  contains
    procedure :: InitTwoBodyRelSpaceIsoMBasis
    procedure :: FinTwoBodyRelSpaceIsoMBasis
    procedure :: GetTwoBodyRelChanIsoMBasisFromJPZ
    procedure :: GetTwoBodyRelChanIsoMBasisFromCh
    procedure :: SetMeshWeight
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetNMesh
    procedure :: GetJmax
    generic :: init => InitTwoBodyRelSpaceIsoMBasis
    generic :: fin => FinTwoBodyRelSpaceIsoMBasis
    generic :: GetChannel => GetTwoBodyRelChanIsoMBasisFromJPZ, GetTwoBodyRelChanIsoMBasisFromCh
    generic :: setp => SetMeshWeight
  end type TwoBodyRelSpaceIsoMBasis

contains

  subroutine FinTwoBodyRelSpaceIsoMBasis(this)
    class(TwoBodyRelSpaceIsoMBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%xmin = 0.d0
    this%xmax = 0.d0
    this%initialized = .false.
  end subroutine FinTwoBodyRelSpaceIsoMBasis

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelSpaceIsoMBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, t) result(idx)
    class(TwoBodyRelSpaceIsoMBasis), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j, p, t)
  end function GetIndex

  function GetNMesh(this) result(NMesh)
    class(TwoBodyRelSpaceIsoMBasis), intent(in) :: this
    integer :: NMesh
    NMesh = this%NMesh
  end function GetNMesh

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelSpaceIsoMBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine InitTwoBodyRelSpaceIsoMBasis(this, xmin, xmax, NMesh, Jmax, &
        & j_in, p_in, t_in)
    class(TwoBodyRelSpaceIsoMBasis), intent(inout) :: this
    integer, intent(in) :: Nmesh, Jmax
    real(8), intent(in) :: xmin, xmax
    integer, intent(in), optional :: j_in, p_in, t_in
    integer :: j, p, t, s, n, ich, cnt, loop, l

    this%initialized = .true.
    this%NMesh = NMesh
    this%Jmax = Jmax
    this%xmin = xmin
    this%xmax = xmax
    do loop = 1, 2
      ich = 0
      do t = 0, 1
        do j = 0, this%jmax
          do p = 1, -1, -2
            if(present(j_in)) then
              if(j /= j_in) cycle
            end if
            if(present(p_in)) then
              if(p /= p_in) cycle
            end if
            if(present(t_in)) then
              if(t /= t_in) cycle
            end if
            cnt = 0
            do s = 0, 1
              do l = iabs(j-s), j+s
                if((-1) ** l /= p) cycle
                if(l < 0) cycle
                do n = 0, this%NMesh
                  if((-1) **(l + s + t) == 1) cycle
                  cnt = cnt + 1
                end do
#ifdef TwoBodyRelativeSpaceDebug
                if(cnt /= 0) then
                  write(*,'(a,4i4)') 'InitTwoBodyRelSpaceHOBasis |L,S,J,Tz>', &
                      & l,s,j,itz
                end if
#endif
              end do
            end do
            if(cnt > 0) then
              ich = ich + 1
              if(loop == 2) then
                this%jpt2idx(j, p, t) = ich
                call this%jpt(ich)%init(this%NMesh, j, p, t)
              end if
            end if
          end do
        end do
      end do
      if(loop == 1) then
        this%Nchan = ich
        allocate(this%jpt(ich))
        allocate(this%jpt2idx(0:jmax, -1:1, -1:1))
        this%jpt2idx(:,:,:) = 0
      end if
    end do
    if(present(j_in) .and. present(p_in) .and. present(t_in) .and. this%NChan == 0) then
      write(*,*) 'In InitTwoBodyRelSpace, spin, parity, and pn are not correct'
      stop
    end if

  end subroutine InitTwoBodyRelSpaceIsoMBasis

  subroutine SetMeshWeight(this, x, w)
    class(TwoBodyRelSpaceIsoMBasis), intent(inout) :: this
    real(8), intent(in) :: x(:), w(:)
    integer :: ch
    do ch = 1, this%NChan
      call this%jpt(ch)%SetMeshWeight(x,w)
    end do
  end subroutine SetMeshWeight

  function GetTwoBodyRelChanIsoMBasisFromJPZ(this, j, p, t) result(r)
    class(TwoBodyRelSpaceIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, t
    type(TwoBodyRelChanIsoMBasis), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,t) )
  end function GetTwoBodyRelChanIsoMBasisFromJPZ

  function GetTwoBodyRelChanIsoMBasisFromCh(this, ch) result(r)
    class(TwoBodyRelSpaceIsoMBasis), target, intent(in) :: this
    integer, intent(in) :: ch
    type(TwoBodyRelChanIsoMBasis), pointer :: r
    r => null()
    if(ch == 0) return
    r => this%jpt(ch)
  end function GetTwoBodyRelChanIsoMBasisFromCh

end module TwoBodyRelSpaceIsoMeshBasis
