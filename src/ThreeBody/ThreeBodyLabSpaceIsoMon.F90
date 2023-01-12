!
! Three-body state with isospin formalism | ab (Jab Tab) c: JT>.
! Normally, I can't calculate e3max > 18 or so. v('_')v
! For three-body matrix elements, I want to calculate
! <ab (Jab Tab) c: JT | V | de (Jde Tde) f: JT>
! where, ja = jd, jb = je, jc = jf, Jab = Jde.
! This limitation allows to calculate up to e3max = 40 or so, hopefully. b(^_^)d
!
module ThreeBodyLabSpaceIsoMon
  !$ use omp_lib
  use ClassSys
  use SingleParticleState
  use ThreeBodyLabSpaceIsoNO2B, only: OneBodyChannels
  use MPIFunction, only: myrank
  implicit none

  public MonThreeBodyLabIsoChan
  public MonThreeBodyLabIsoSpace
  private

  type :: MonThreeBodyLabIsoChan
    type(OrbitsIsospin), pointer :: sps
    integer, private :: t  = -1
    integer, private :: j1 = -1
    integer, private :: p1 = -100
    integer, private :: j2 = -1
    integer, private :: p2 = -100
    integer, private :: j3 = -1
    integer, private :: p3 = -100
    integer, private :: n_states = 0
    integer, allocatable :: IndexABCTE(:,:)
    integer, allocatable :: TNs2Index(:,:,:,:) ! [tab,na,nb,nc] => idx
    integer, allocatable :: TNs2Phase(:,:,:,:) ! [tab,na,nb,nc] => idx
  contains
    procedure :: InitMonThreeBodyLabIsoChan
    procedure :: FinMonThreeBodyLabIsoChan
    procedure :: GetT
    procedure :: GetParity
    procedure :: GetJ1
    procedure :: GetP1
    procedure :: GetJ2
    procedure :: GetP2
    procedure :: GetJ3
    procedure :: GetP3
    procedure :: GetNumberStates
    procedure :: GetABCTE
    procedure :: GetIndexFromABCT
    generic :: init => InitMonThreeBodyLabIsoChan
    generic :: fin => FinMonThreeBodyLabIsoChan
    generic :: GetIndex => GetIndexFromABCT
  end type MonThreeBodyLabIsoChan

  type :: MonThreeBodyLabIsoSpace
    type(MonThreeBodyLabIsoChan), allocatable :: chan(:)
    type(OrbitsIsospin), pointer :: sps
    type(OneBodyChannels) :: one
    integer, allocatable, private :: idx_sorted(:)
    integer, private :: emax = -1
    integer, private :: e2max= -1
    integer, private :: e3max= -1
    integer, private :: NChan= 0
    integer, private :: NT = 0
    integer, private :: NMon = 0
    integer, allocatable :: TMon2Ch(:,:,:,:)
    logical, allocatable :: Flip(:,:,:,:)
    logical :: is_Constructed=.false.
  contains
    procedure :: InitMonThreeBodyLabIsoSpace
    procedure :: FinMonThreeBodyLabIsoSpace
    procedure :: GetNumberChannels
    procedure :: GetEmax
    procedure :: GetE2max
    procedure :: GetE3max
    procedure :: GetIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromQNs
    procedure :: GetSortedChannelIndex
    generic :: init => InitMonThreeBodyLabIsoSpace
    generic :: fin => FinMonThreeBodyLabIsoSpace
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromQNs
  end type MonThreeBodyLabIsoSpace

contains

  !
  ! MonThreeBodyLabIsoChan
  !
  function GetT(this) result(T)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetParity(this) result(prty)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: prty
    prty = this%GetP1() * this%GetP2() * this%GetP3()
  end function GetParity

  function GetJ1(this) result(J1)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: J1
    J1 = this%J1
  end function GetJ1

  function GetP1(this) result(P1)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: P1
    P1 = this%P1
  end function GetP1

  function GetJ2(this) result(J2)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: J2
    J2 = this%J2
  end function GetJ2

  function GetP2(this) result(P2)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: P2
    P2 = this%P2
  end function GetP2

  function GetJ3(this) result(J3)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: J3
    J3 = this%J3
  end function GetJ3

  function GetP3(this) result(P3)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: P3
    P3 = this%P3
  end function GetP3

  function GetNumberStates(this) result(n_states)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetABCTE(this,i) result(QNs)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer :: i
    integer :: QNs(5)
    QNs(:) = this%IndexABCTE(:,i)
  end function GetABCTE

  function GetIndexFromABCT(this,i1,i2,i3,t12) result(idx)
    class(MonThreeBodyLabIsoChan), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,t12
    integer :: idx
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    o1 => this%sps%GetOrbit(i1)
    o2 => this%sps%GetOrbit(i2)
    o3 => this%sps%GetOrbit(i3)
    idx = this%TNs2Index(t12, o1%n, o2%n, o3%n)
  end function GetIndexFromABCT

  !
  ! MonThreeBodyLabIsoSpace
  !
  function GetNumberChannels(this) result(NChan)
    class(MonThreeBodyLabIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetEmax(this) result(Emax)
    class(MonThreeBodyLabIsoSpace), intent(in) :: this
    integer :: Emax
    Emax = this%emax
  end function GetEmax

  function GetE2max(this) result(E2max)
    class(MonThreeBodyLabIsoSpace), intent(in) :: this
    integer :: E2max
    E2max = this%e2max
  end function GetE2max

  function GetE3max(this) result(E3max)
    class(MonThreeBodyLabIsoSpace), intent(in) :: this
    integer :: E3max
    E3max = this%e3max
  end function GetE3max

  function GetIndex(this,t,j1,p1,j2,p2,j3,p3) result(idx)
    class(MonThreeBodyLabIsoSpace), intent(in) :: this
    integer, intent(in) :: t,j1,p1,j2,p2,j3,p3
    integer :: idx, ch1, ch2, ch3
    idx = 0
    ch1 = this%one%GetIndex(j1,p1)
    ch2 = this%one%GetIndex(j2,p2)
    ch3 = this%one%GetIndex(j3,p3)
    if(ch1 * ch2 * ch3 == 0) return
    idx = this%TMon2Ch(t,ch1,ch2,ch3)
  end function GetIndex

  function GetChannelFromIndex(this,idx) result(r)
    class(MonThreeBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    type(MonThreeBodyLabIsoChan), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%Chan(idx)
  end function GetChannelFromIndex

  function GetChannelFromQNs(this,t,j1,p1,j2,p2,j3,p3) result(r)
    class(MonThreeBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: t,j1,p1,j2,p2,j3,p3
    type(MonThreeBodyLabIsoChan), pointer :: r
    r => this%GetChannel( this%GetIndex(t,j1,p1,j2,p2,j3,p3) )
  end function GetChannelFromQNs

  function GetSortedChannelIndex(this, idx) result(sorted_idx)
    class(MonThreeBodyLabIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    integer :: sorted_idx
    sorted_idx = 0
    if(idx == 0) return
    sorted_idx = this%idx_sorted(idx)
  end function GetSortedChannelIndex
  !
  !
  !

  subroutine FinMonThreeBodyLabIsoChan(this)
    class(MonThreeBodyLabIsoChan), intent(inout) :: this
    this%sps => null()
    this%t = -1
    this%j1 = -1
    this%p1 = -100
    this%j2 = -1
    this%p2 = -100
    this%j3 = -1
    this%p3 = -100
    this%n_states = 0
    deallocate(this%IndexABCTE)
    deallocate(this%TNs2Index)
  end subroutine FinMonThreeBodyLabIsoChan

  subroutine InitMonThreeBodyLabIsoChan(this,sps,t,j1,p1,j2,p2,j3,p3,e2max,e3max)
    class(MonThreeBodyLabIsoChan), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: t, j1, p1, j2, p2, j3, p3, e2max, e3max
    integer :: n_states
    this%sps => sps
    this%t = t
    this%j1 = j1
    this%p1 = p1
    this%j2 = j2
    this%p2 = p2
    this%j3 = j3
    this%p3 = p3
    call set_mon_red_channel_indices(this, sps, e2max, e3max, n_states, .false.)
    call set_mon_red_channel_indices(this, sps, e2max, e3max, n_states, .true.)
  end subroutine InitMonThreeBodyLabIsoChan

  subroutine set_mon_red_channel_indices(this, sps, e2max, e3max, n_states, set_mode)
    type(MonThreeBodyLabIsoChan), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer, intent(inout) :: n_states
    logical, intent(in) :: set_mode
    integer :: i1, i2, i3, t12
    integer :: l1, l2, l3, n1, n2, n3, e1, e2, e3
    integer :: n2max

    if(set_mode) then
      this%n_states = n_states
      allocate(this%IndexABCTE(5,n_states))
      allocate(this%TNs2Index(0:1, 0:sps%emax/2, 0:sps%emax/2, 0:sps%emax/2))
      allocate(this%TNs2Phase(0:1, 0:sps%emax/2, 0:sps%emax/2, 0:sps%emax/2))
      this%TNs2Index(:,:,:,:) = 0
      this%TNs2Phase(:,:,:,:) = 0
    end if
    l1 = (this%j1+1)/2
    l2 = (this%j2+1)/2
    l3 = (this%j3+1)/2
    if(this%p1 == (-1)**((this%j1-1)/2)) l1 = (this%j1-1)/2
    if(this%p2 == (-1)**((this%j2-1)/2)) l2 = (this%j2-1)/2
    if(this%p3 == (-1)**((this%j3-1)/2)) l3 = (this%j3-1)/2

    n_states = 0
    do n1 = 0, (sps%emax - l1)/2
      e1 = 2*n1 + l1
      i1 = sps%nlj2idx(n1,l1,this%j1)

      n2max = (sps%emax - l2)/2
      !if(this%j1==this%j2 .and. this%p1==this%p2) n2max = n1
      do n2 = 0, n2max
        e2 = 2*n2 + l2
        i2 = sps%nlj2idx(n2,l2,this%j2)
        if( e1 + e2 > e2max ) cycle

        do n3 = 0, (sps%emax - l3)/2
          e3 = 2*n3 + l3
          i3 = sps%nlj2idx(n3,l3,this%j3)
          if( e1 + e3 > e2max ) cycle
          if( e2 + e3 > e2max ) cycle
          if( e1 + e2 + e3 > e3max ) cycle

          do t12 = 0, 1
            n_states = n_states+1
            if(.not. set_mode) cycle
            this%IndexABCTE(:,n_states) = [i1,i2,i3,t12,e1+e2+e3]
            this%TNs2Index(t12,n1,n2,n3) = n_states
            this%TNs2Phase(t12,n1,n2,n3) = 1
            !if(this%j1==this%j2 .and. this%p1==this%p2) then
            !  this%TNs2Index(t12,n2,n1,n3) = n_states
            !  if(n1/=n2) this%TNs2Index(t12,n2,n1,n3) = (-1)**t12
            !end if
          end do
        end do
      end do
    end do
  end subroutine set_mon_red_channel_indices

  subroutine FinMonThreeBodyLabIsoSpace(this)
    class(MonThreeBodyLabIsoSpace), intent(inout) :: this
    integer :: ch
    do ch = 1, this%GetNumberChannels()
      call this%Chan(ch)%fin()
    end do
    deallocate(this%Chan)
    this%sps => null()
    call this%one%fin()
    this%emax = -1
    this%e2max= -1
    this%e3max= -1
    this%NChan= 0
    deallocate(this%TMon2Ch)
    this%is_Constructed = .false.
  end subroutine FinMonThreeBodyLabIsoSpace

  subroutine InitMonThreeBodyLabIsoSpace(this, sps, e2max, e3max)
    use Profiler
    class(MonThreeBodyLabIsoSpace), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: e2max, e3max
    integer :: ch, ch1, ch2, ch3, nch1
    integer :: Ttot, j1, p1, j2, p2, j3, p3, n_states, idx, ich
    integer, allocatable :: idx_loop(:,:)
    integer, allocatable :: ndim_array(:)
    type(sys) :: s
    real(8) :: mem, time

    time = omp_get_wtime()
    this%is_Constructed = .true.
    this%sps => sps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max
    call this%one%init(sps)
    nch1 = this%one%GetNumberChannels()
    allocate(this%TMon2Ch(1:3, nch1, nch1, nch1))
    allocate(this%flip(1:3, nch1, nch1, nch1))
    this%TMon2Ch(:,:,:,:) = 0
    ch = 0
    do Ttot = 1,3,2
      do ch1 = 1, nch1
        j1 = this%one%GetJ(ch1)
        p1 = this%one%GetP(ch1)
        do ch2 = 1, ch1
          j2 = this%one%GetJ(ch2)
          p2 = this%one%GetP(ch2)
          do ch3 = 1, nch1
            j3 = this%one%GetJ(ch3)
            p3 = this%one%GetP(ch3)

            n_states = count_number_states(sps, e2max, e3max, j1, p1, j2, p2, j3, p3, Ttot)
            if(n_states==0) cycle
            ch = ch + 1
            this%TMon2Ch(Ttot, ch1, ch2, ch3) = ch
            this%TMon2Ch(Ttot, ch2, ch1, ch3) = ch
            this%flip(Ttot, ch1, ch2, ch3) = .false.
            if(ch1/=ch2) this%flip(Ttot, ch2, ch1, ch3) = .true.

          end do
        end do
      end do
    end do
    this%NChan = ch
    allocate(this%Chan( this%GetNumberChannels() ))
    allocate(idx_loop(4,this%GetNumberChannels() ))
    do Ttot = 1,3,2
      do ch1 = 1, nch1
        do ch2 = 1, ch1
          do ch3 = 1, nch1
            ch = this%TMon2Ch(Ttot, ch1, ch2, ch3)
            if(ch == 0) cycle
            idx_loop(:,ch) = [Ttot, ch1, ch2, ch3]
          end do
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(ch, Ttot, ch1, ch2, ch3, j1, p1, j2, p2, j3, p3)
    do ch = 1, this%GetNumberChannels()
      Ttot = idx_loop(1,ch)
      ch1  = idx_loop(2,ch)
      ch2  = idx_loop(3,ch)
      ch3  = idx_loop(4,ch)

      j1 = this%one%GetJ(ch1)
      p1 = this%one%GetP(ch1)
      j2 = this%one%GetJ(ch2)
      p2 = this%one%GetP(ch2)
      j3 = this%one%GetJ(ch3)
      p3 = this%one%GetP(ch3)

      call this%Chan(ch)%init(sps,ttot,j1,p1,j2,p2,j3,p3,e2max,e3max)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(idx_loop)
    allocate(this%idx_sorted( this%GetNumberChannels() ))
    allocate(ndim_array( this%GetNumberChannels() ))
    do ich = 1, this%GetNumberChannels()
      ndim_array(ich) = this%Chan(ich)%GetNumberStates()
    end do
    do ich = 1, this%GetNumberChannels()
      idx = maxloc(ndim_array,dim=1)
      this%idx_sorted(ich) = idx
      ndim_array(idx) = -1
    end do
    deallocate(ndim_array)
    mem = dble(this%GetNumberChannels()) * dble(sps%emax/2+1)**3 * 4.d0 / 1024.d0**3
    if(myrank == 0) write(*,"(a, f10.4,a)") "# Memory used for MonThreeBodyIsoSpace is ", mem, " GB"
    call timer%add(s%str("InitMonThreeBodyIsoSpace"), omp_get_wtime()-time)
  end subroutine InitMonThreeBodyLabIsoSpace

  function count_number_states(sps, e2max, e3max, j1, p1, j2, p2, j3, p3, T) result(r)
    use MyLibrary, only: triag
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max, e3max, j1, p1, j2, p2, j3, p3, T
    integer :: r
    integer :: t12
    integer :: l1, l2, l3, n1, n2, n3, e1, e2, e3
    integer :: n2max

    l1 = (j1+1)/2
    l2 = (j2+1)/2
    l3 = (j3+1)/2
    if(p1 == (-1)**((j1-1)/2)) l1 = (j1-1)/2
    if(p2 == (-1)**((j2-1)/2)) l2 = (j2-1)/2
    if(p3 == (-1)**((j3-1)/2)) l3 = (j3-1)/2
    r = 0
    do n1 = 0, (sps%emax - l1)/2
      e1 = 2*n1 + l1

      n2max = (sps%emax - l2)/2
      !if(j1==j2 .and. p1==p2) n2max = n1
      do n2 = 0, n2max
        e2 = 2*n2 + l2
        if( e1 + e2 > e2max ) cycle

        do n3 = 0, (sps%emax - l3)/2
          e3 = 2*n3 + l3
          if( e1 + e3 > e2max ) cycle
          if( e2 + e3 > e2max ) cycle
          if( e1 + e2 + e3 > e3max ) cycle

          do t12 = 0, 1
            if(triag( 2*t12, 1, T )) cycle
            r = r + 1
          end do
        end do
      end do
    end do
  end function count_number_states

end module ThreeBodyLabSpaceIsoMon

