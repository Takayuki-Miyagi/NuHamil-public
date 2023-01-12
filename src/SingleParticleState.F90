module SingleParticleState
  implicit none

  public :: SingleParticleOrbitIsospin
  public :: SingleParticleOrbit
  public :: OrbitsIsospin
  public :: Orbits

  private :: SetSingleParticleOrbitIsospin
  private :: SetSingleParticleOrbit
  private :: InitOrbitsIsospin
  private :: FinOrbitsIsospin
  private :: InitOrbits
  private :: FinOrbits
  private :: iso2pn
  private :: pn2iso
  private :: GetLabelFromIndex
  private :: GetLabelFromIndexIsospin
  private :: GetIndexFromLabel
  private :: GetIndexFromLabelIsospin
  private :: GetOrbit

  type :: SingleParticleOrbitIsospin
    integer :: n = -1
    integer :: l = -1
    integer :: j = -1
    integer :: e = -1
    integer :: idx = -1
  contains
    procedure :: set => SetSingleParticleOrbitIsospin
  end type SingleParticleOrbitIsospin

  type :: SingleParticleOrbit
    integer :: n = -1
    integer :: l = -1
    integer :: j = -1
    integer :: z =  0
    integer :: e = -1
    integer :: idx = -1
    real(8) :: occ = -1.d0
    integer :: ph = -1 ! 0:hole,1:particle,2:partial occ
  contains
    procedure :: set => SetSingleParticleOrbit
    procedure :: SetOccupation
    procedure :: SetParticleHole

  end type SingleParticleOrbit

  type :: OrbitsIsospin
    integer, allocatable :: nlj2idx(:,:,:)
    type(SingleParticleOrbitIsospin), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
  contains
    procedure :: init => InitOrbitsIsospin
    procedure :: fin => FinOrbitsIsospin
    procedure :: iso2pn
    procedure :: GetLabelFromIndexIsospin
    procedure :: GetIndexFromLabelIsospin
    procedure :: GetOrbitIsospin
    generic :: GetOrbit => GetOrbitIsospin
  end type OrbitsIsospin

  type :: Orbits
    integer, allocatable :: nljz2idx(:,:,:,:)
    type(SingleParticleOrbit), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
    character(:), allocatable :: mode
  contains
    procedure :: init => InitOrbits
    procedure :: fin => FinOrbits
    procedure :: pn2iso
    procedure :: GetLabelFromIndex
    procedure :: GetIndexFromLabel
    procedure :: GetOrbit
  end type Orbits

  character(1), private :: OrbitalAngMom(21)=['s',&
      & 'p','d','f','g','h','i','k','l','m','n',&
      & 'o','q','r','t','u','v','w','x','y','z']

contains

  subroutine FinOrbitsIsospin(this)
    class(OrbitsIsospin), intent(inout) :: this
    if(.not. this%is_constructed) return
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In FinOrbitsIsospin'
#endif
    deallocate(this%orb)
    deallocate(this%nlj2idx)
    this%is_constructed = .false.
  end subroutine FinOrbitsIsospin

  subroutine InitOrbitsIsospin(this, emax, lmax)
    class(OrbitsIsospin), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, cnt

#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In InitOrbitsIsospin:'
#endif
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    if(this%lmax < this%emax) then
      write(*,"(a,i3)") "# Lmax truncation is used in the single-particle orbits: Lmax=", this%lmax
    end if
    allocate(this%nlj2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1))
    this%nlj2idx(:,:,:) = 0
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
          call this%orb(cnt)%set(n,l,j,cnt)
          this%nlj2idx(n,l,j) = cnt
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine InitOrbitsIsospin

  function GetOrbit(this, idx) result(r)
    class(Orbits), intent(in), target :: this
    integer, intent(in) :: idx
    type(SingleParticleOrbit), pointer :: r
    r => this%orb(idx)
  end function GetOrbit

  function GetOrbitIsospin(this, idx) result(r)
    class(OrbitsIsospin), intent(in), target :: this
    integer, intent(in) :: idx
    type(SingleParticleOrbitIsospin), pointer :: r
    r => this%orb(idx)
  end function GetOrbitIsospin

  subroutine FinOrbits(this)
    class(Orbits), intent(inout) :: this
    if(.not. this%is_constructed) return
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In FinOrbits'
#endif
    deallocate(this%orb)
    deallocate(this%nljz2idx)
    this%is_constructed = .false.
  end subroutine FinOrbits

  subroutine InitOrbits(this, emax, lmax, mode)
    class(Orbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    character(*), optional, intent(in) :: mode
    integer :: e, l, n, s, j, z, cnt

    this%mode = ""
    if(present(mode)) this%mode = mode
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In InitOrbits:'
#endif
    select case(this%mode)
    case("kshell","KSHELL","Kshell")
      call init_orbits_kshell_format(this, emax, lmax)
      return
    case default
      this%emax = emax
      this%lmax = emax
      if(present(lmax)) this%lmax = lmax
      if(this%lmax < this%emax) then
        write(*,"(a,i3)") "# Lmax truncation is used in the single-particle orbits: Lmax=", this%lmax
      end if
      allocate(this%nljz2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1, -1:1))
      this%nljz2idx(:,:,:,:) = 0
      cnt = 0
      do e = 0, this%emax
        do l = 0, min(e,this%lmax)
          if(mod(e - l, 2) == 1) cycle
          n = (e - l) / 2
          do s = -1, 1, 2
            j = 2*l + s
            if(j < 1) cycle
            do z = -1, 1, 2
              cnt = cnt + 1
            end do
          end do
        end do
      end do
      this%norbs = cnt
      allocate(this%orb(this%norbs))
      cnt = 0
      do e = 0, this%emax
        do l = 0, min(e,this%lmax)
          if(mod(e - l, 2) == 1) cycle
          n = (e - l) / 2
          do s = -1, 1, 2
            j = 2*l + s
            if(j < 1) cycle
            do z = -1, 1, 2
              cnt = cnt + 1
              call this%orb(cnt)%set(n,l,j,z,cnt)
              this%nljz2idx(n,l,j,z) = cnt
            end do
          end do
        end do
      end do
      this%is_constructed = .true.
    end select
  end subroutine InitOrbits

  subroutine init_orbits_kshell_format(this, emax, lmax)
    type(Orbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, z, cnt

    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nljz2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1, -1:1))
    this%nljz2idx(:,:,:,:) = 0
    cnt = 0
    do z = -1, 1, 2
      do e = 0, this%emax
        do l = 0, min(e,this%lmax)
          if(mod(e - l, 2) == 1) cycle
          n = (e - l) / 2
          do s = -1, 1, 2
            j = 2*l + s
            if(j < 1) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do z = -1, 1, 2
      do e = 0, this%emax
        do l = 0, min(e,this%lmax)
          if(mod(e - l, 2) == 1) cycle
          n = (e - l) / 2
          do s = -1, 1, 2
            j = 2*l + s
            if(j < 1) cycle
            cnt = cnt + 1
            call this%orb(cnt)%set(n,l,j,z,cnt)
            this%nljz2idx(n,l,j,z) = cnt
          end do
        end do
      end do
    end do
    this%is_constructed = .true.

  end subroutine init_orbits_kshell_format


  subroutine SetSingleParticleOrbitIsospin(this, n, l, j, idx)
    class(SingleParticleOrbitIsospin), intent(inout) :: this
    integer, intent(in) :: n, l, j, idx
    this%n = n
    this%l = l
    this%j = j
    this%e = 2*n + l
    this%idx = idx
#ifdef SingleParticleStateDebug
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
        & 'index=', idx, ', n=', n, ', l=', l, ', j=', j, ', e=', 2*n+l
#endif
  end subroutine SetSingleParticleOrbitIsospin

  subroutine SetSingleParticleOrbit(this, n, l, j, z, idx)
    class(SingleParticleOrbit), intent(inout) :: this
    integer, intent(in) :: n, l, j, z, idx
    this%n = n
    this%l = l
    this%j = j
    this%z = z
    this%e = 2*n + l
    this%idx = idx
#ifdef SingleParticleStateDebug
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
        & 'indx=',idx,', n=', n, ', l=', l, ', j=', j, ', z=', z, ', e=', 2*n+l
#endif
  end subroutine SetSingleParticleOrbit

  subroutine SetOccupation(this,occ)
    class(SingleParticleOrbit), intent(inout) :: this
    real(8), intent(in) :: occ
    this%occ = occ
  end subroutine SetOccupation

  subroutine SetParticleHole(this,ph_label)
    class(SingleParticleOrbit), intent(inout) :: this
    integer, intent(in) :: ph_label
    this%ph = ph_label
  end subroutine SetParticleHole

  function iso2pn(this, sps, idx, z) result(r)
    class(OrbitsIsospin), intent(in) :: this
    type(Orbits), intent(in) :: sps
    integer, intent(in) :: idx, z
    integer :: r

    r = 0
    if(z /= -1 .and. z /= 1) then
      write(*,'(a,i3)') "Error in iso2pn, tz has to be -1 or 1. tz = ", z
      stop
    end if
    if(this%orb(idx)%e > sps%emax) return
    if(this%orb(idx)%l > sps%lmax) return
    r=sps%nljz2idx(this%orb(idx)%n,this%orb(idx)%l,this%orb(idx)%j,z)
  end function iso2pn

  function pn2iso(this, sps, idx) result(r)
    class(Orbits), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: idx
    integer :: r

    r = 0
    if(this%orb(idx)%e > sps%emax) return
    if(this%orb(idx)%l > sps%lmax) return
    r=sps%nlj2idx(this%orb(idx)%n,this%orb(idx)%l,this%orb(idx)%j)
  end function pn2iso

  function GetLabelFromIndex(this, idx) result(r)
    use ClassSys
    class(Orbits), intent(in) :: this
    integer, intent(in) :: idx
    type(str) :: r
    type(sys) :: s
    integer :: n, l, j, z

    if(idx > this%norbs) then
      write(*,'(a)') "Warning in GetLabelFromIndex"
      r = ''
      return
    end if

    n = this%orb(idx)%n
    l = this%orb(idx)%l
    j = this%orb(idx)%j
    z = this%orb(idx)%z
    if(z == -1) then
      r = 'p' + s%str(n) + OrbitalAngMom(l+1) + s%str(j) + '/2'
      return
    end if

    if(z ==  1) then
      r = 'n' + s%str(n) + OrbitalAngMom(l+1) + s%str(j) + '/2'
      return
    end if
    write(*,'(a)') 'Error in GetLabelFromIndex'
  end function GetLabelFromIndex

  function GetLabelFromIndexIsospin(this,idx) result(r)
    use ClassSys
    class(OrbitsIsospin), intent(in) :: this
    integer, intent(in) :: idx
    type(str) :: r
    type(sys) :: s
    integer :: n, l, j

    if(idx > this%norbs) then
      write(*,'(a)') "Warning in GetLabelFromIndexIsospin"
      r = ''
      return
    end if

    n = this%orb(idx)%n
    l = this%orb(idx)%l
    j = this%orb(idx)%j
    r = s%str(n) + OrbitalAngMom(l+1) + s%str(j) + '/2'
  end function GetLabelFromIndexIsospin

  function GetIndexFromLabel(this, label) result(r)
    class(Orbits), intent(in) :: this
    character(*), intent(in) :: label
    integer :: r
    character(:), allocatable :: str, str_n, str_l, str_j
    integer :: n, l, j, z, i
    r = 0
    if(label == '') then
      write(*,'(a)') "Warning in GetIndexFromLabel"
      return
    end if
    if(label(1:1) /= 'p' .and. label(1:1) /= 'n') then
      write(*,'(a)') 'Error in GetIndexFromLabel: label format error'
      return
    end if

    if(label(1:1) == 'p') z = -1
    if(label(1:1) == 'n') z =  1

    str = label(2:)
    str_n = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_n = trim(str_n) // str(1:1)
      str = str(2:)
    end do
    str_l = str(1:1)

    str = str(2:)
    str_j = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_j = trim(str_j) // str(1:1)
      str = str(2:)
    end do

    do i = 1, size(OrbitalAngMom)
      if(str_l == OrbitalAngMom(i)) then
        l = i - 1
        exit
      end if
      if(i == size(OrbitalAngMom)) then
        write(*,'(a)') "Error in GetIndexFromLabel: orbital angular momentum limitation"
        l = 0
        exit
      end if
    end do

    read(str_n,*) n
    read(str_j,*) j
    r = this%nljz2idx(n,l,j,z)
  end function GetIndexFromLabel

  function GetIndexFromLabelIsospin(this, label) result(r)
    class(OrbitsIsospin), intent(in) :: this
    character(*), intent(in) :: label
    integer :: r
    character(:), allocatable :: str, str_n, str_l, str_j
    integer :: n, l, j, i
    r = 0

    if(label == '') then
      write(*,'(a)') "Warning in GetIndexFromLabelIsospin"
      return
    end if

    str = label(1:)
    str_n = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_n = trim(str_n) // str(1:1)
      str = str(2:)
    end do
    str_l = str(1:1)

    str = str(2:)
    str_j = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_j = trim(str_j) // str(1:1)
      str = str(2:)
    end do

    do i = 1, size(OrbitalAngMom)
      if(str_l == OrbitalAngMom(i)) then
        l = i - 1
        exit
      end if
      if(i == size(OrbitalAngMom)) then
        write(*,'(a)') "Error in GetIndexFromLabelIsospin: orbital angular momentum limitation"
        l = 0
        exit
      end if
    end do

    read(str_n,*) n
    read(str_j,*) j
    r = this%nlj2idx(n,l,j)
  end function GetIndexFromLabelIsospin
end module SingleParticleState

! main program for check
!program test
!  use SingleParticleState
!  type(Orbits) :: o
!  type(OrbitsIsospin) :: io
!  character(:), allocatable :: str
!  call o%init(4)
!  str = o%GetLabelFromIndex(1)
!  write(*,*) 1, o%GetIndexFromLabel(str)
!  str = o%GetLabelFromIndex(20)
!  write(*,*) 20, o%GetIndexFromLabel(str)
!  str = o%GetLabelFromIndex(15)
!  write(*,*) 15, o%GetIndexFromLabel(str)
!  call o%fin()
!
!
!  call o%init(4,2)
!  call o%fin()
!
!  call io%init(4)
!  str = io%GetLabelFromIndexIsospin(1)
!  write(*,*) 1, io%GetIndexFromLabelIsospin(str)
!  str = io%GetLabelFromIndexIsospin(20)
!  write(*,*) 20, io%GetIndexFromLabelIsospin(str)
!  str = io%GetLabelFromIndexIsospin(15)
!  write(*,*) 15, io%GetIndexFromLabelIsospin(str)
!  call io%fin()
!
!  call io%init(4,2)
!  call io%fin()
!end program test


