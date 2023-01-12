module ThreeBodyJacobiSpaceIso
  use ThreeBodyJacChanIso
  implicit none
  public :: ThreeBodyJacIsoSpace

  private :: InitThreeBodyJacIsoSpace
  private :: FinThreeBodyJacIsoSpace
  private :: GetFrequency
  private :: GetNumberChannels
  private :: GetJmax
  private :: GetIndex
  private :: GetChannelFromIndex
  private :: GetChannelFromJPT
  private :: GetMemory

  Type :: ThreeBodyJacIsoSpace
    real(8), private :: hw = -1.d0
    integer, private :: NChan = 0
    integer, private :: Jmax = -1
    type(ThreeBodyJacIsoChan), allocatable :: jpt(:)
    integer, private, allocatable :: jpt2idx(:,:,:)
    logical :: is_Constructed=.false.
  contains
    procedure :: InitThreeBodyJacIsoSpace
    procedure :: FinThreeBodyJacIsoSpace
    procedure :: GetFrequency
    procedure :: GetNumberChannels
    procedure :: GetJmax
    procedure :: GetIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromJPT
    procedure :: GetMemory
    generic :: init => InitThreeBodyJacIsoSpace
    generic :: fin => FinThreeBodyJacIsoSpace
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromJPT
  end Type ThreeBodyJacIsoSpace

contains

  function GetMemory(this) result(r)
    class(ThreeBodyJacIsoSpace), intent(in) :: this
    real(8) :: r
    integer :: ich
    r = 0.d0
    do ich = 1, this%GetNumberChannels()
      r = r + this%jpt(ich)%GetMemory()
    end do
  end function GetMemory

  function GetChannelFromIndex(this,i) result(r)
    class(ThreeBodyJacIsoSpace), intent(in), target :: this
    integer, intent(in) :: i
    type(ThreeBodyJacIsoChan), pointer :: r
    r => null()
    if(i == 0) return
    r => this%jpt(i)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this,j,p,t) result(r)
    class(ThreeBodyJacIsoSpace), intent(in), target :: this
    integer, intent(in) :: j,p,t
    type(ThreeBodyJacIsoChan), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,t) )
  end function GetChannelFromJPT

  function GetFrequency(this) result(hw)
    class(ThreeBodyJacIsoSpace), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetIndex(this, j, p, t) result(idx)
    class(ThreeBodyJacIsoSpace), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j,p,t)
  end function GetIndex

  function GetNumberChannels(this) result(NChan)
    class(ThreeBodyJacIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetJmax(this) result(Jmax)
    class(ThreeBodyJacIsoSpace), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  subroutine FinThreeBodyJacIsoSpace(this)
    class(ThreeBodyJacIsoSpace), intent(inout) :: this
    integer :: ich
    if(.not. this%is_Constructed) return
    do ich = 1, this%NChan
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%hw = -1.d0
    this%NChan = 0
    this%Jmax = -1
    this%is_Constructed=.false.
  end subroutine FinThreeBodyJacIsoSpace

  subroutine InitThreeBodyJacIsoSpace(this, hw, Jmax, Nmax, ramp, Nmax_for_file, path_to_dir)
    use ClassSys, only: str
    class(ThreeBodyJacIsoSpace), intent(inout) :: this
    real(8) :: hw
    integer, intent(in) :: Jmax, Nmax
    type(str) :: ramp, path_to_dir
    integer, intent(in), optional :: Nmax_for_file
    integer :: NChan, ich, j, p, t
    integer :: Nmax_file
    integer :: runit = 22, wunit = 23
    type(str) :: f
    logical :: ex

    this%Jmax = Jmax
    this%hw = hw
    NChan = ((this%Jmax + 1) / 2) * 4
    this%NChan = NChan

    allocate(this%jpt(NChan))
    allocate(this%jpt2idx(1:Jmax, -1:1, 1:3))
    this%jpt2idx(:,:,:) = 0

    Nmax_file = Nmax
    if(present(Nmax_for_file)) then
      Nmax_file = Nmax_for_file
      if(Nmax_file < Nmax) then
        write(*,*) "Error, Three Body Jacobi basis: Nmax_file < Nmax"
        stop
      end if
    end if

    ich = 0
    do t = 1, 3, 2
      do j = 1, this%Jmax, 2
        do p = 1, -1, -2
          ich = ich + 1
          this%jpt2idx(j,p,t) = ich
          Nmax_file = GetRampNmax(j, ramp)
          f = this%jpt(ich)%GetFileName(j,p,t,Nmax_file,path_to_dir)
          inquire(file = f%val, exist = ex)
          if(.not. ex) then
            call this%jpt(ich)%init(hw,j,p,t,Nmax_file,path_to_dir) ! constructor

            ! -- write to file
            open(wunit,file=f%val,status='replace',form='unformatted',access='stream')
            call this%jpt(ich)%writef(wunit)
            close(wunit)
            call this%jpt(ich)%fin()

          end if
          ! -- read from file
          open(runit,file=f%val,status='old',form='unformatted',access='stream')
          call this%jpt(ich)%readf(hw,runit,Nmax)
          close(runit)

        end do
      end do
    end do
    this%is_Constructed=.true.
  end subroutine InitThreeBodyJacIsoSpace

  function GetRampNmax(j, ramp) result(Nmax)
    integer, intent(in) :: j
    integer :: Nmax
    type(str), intent(in) :: ramp
    type(sys) :: s
    type(str), allocatable :: ramps(:)
    integer :: i, j_read, Nmax_read

    Nmax = 0
    if(s%find(ramp, s%str("flat"))) then
      read(ramp%val(5:),*) Nmax
      return
    end if

    if(s%find(ramp, s%str("ramp"))) then
      call s%split(ramp, s%str("-"), ramps)
      read(ramps(1)%val(5:),*) Nmax_read
      Nmax = Nmax_read
      do i = 1, size(ramps)/2
        read(ramps(2*i)%val,*) j_read
        read(ramps(2*i+1)%val,*) Nmax_read
        if(j > j_read) Nmax = Nmax_read
      end do
      return
    end if

    write(*,*) "ramp error"
    stop
  end function GetRampNmax

end module ThreeBodyJacobiSpaceIso
