!
! Three-body state with isospin formalism | ab (Jab Tab) c: JT>.
! Normally, I can't calculate e3max > 18 or so. v('_')v
! All three-body matrix elements are not needed if we use the NO2B approximation.
! What I want is <ab (Jab Tab) c: T | V | de (Jde Tde) f: T>
! with jc = jf, lc = lf, Jab = Jde.
! This limitation allows to calculate up to e3max = 24 or higher. b(^_^)d
!
module ThreeBodyLabSpaceIsoNO2B
  !$ use omp_lib
  use ClassSys
  use MPIFunction, only: myrank
  use SingleParticleState
  implicit none

  private :: InitOneBodyChannels
  private :: FinOneBodyChannels
  private :: GetIndexOneBodyChannel
  private :: InitTwoBodyChannels
  private :: FinTwoBodyChannels
  private :: GetIndexTwoBodyChannel
  private :: InitThreeBodyLabIsoChanNO2B
  private :: FinThreeBodyLabIsoChanNO2B
  private :: GetABCT
  private :: GetJNO2B
  private :: GetParityNO2B
  private :: GetTNO2B
  private :: GetJ12NO2B
  private :: GetP12NO2B
  private :: GetJ3NO2B
  private :: GetP3NO2B
  private :: GetNumberStatesNO2B
  private :: InitThreeBodyLabIsoSpaceNO2B
  private :: FinThreeBodyLabIsoSpaceNO2B
  private :: GetNumberChannelsNO2B
  private :: GetEmaxNO2B
  private :: GetE2maxNO2B
  private :: GetE3maxNO2B
  private :: GetIndexNO2B
  private :: GetChannelNO2BFromQNs
  private :: GetChannelNO2BFromIndex
  private :: InitNO2BThreeBodyIsoChan
  private :: FinNO2BThreeBodyIsoChan
  private :: GetT
  private :: GetJ12
  private :: GetP12
  private :: GetJ3
  private :: GetP3
  private :: GetNumberStates
  private :: GetIndexChanFromQNs
  private :: GetPhase
  private :: InitNO2BThreeBodyIsoSpace
  private :: FinNO2BThreeBodyIsoSpace
  private :: GetNumberChannels
  private :: GetEmax
  private :: GetE2max
  private :: GetE3max
  private :: GetIndex
  private :: GetChannelFromQNs
  private :: GetChannelFromIndex
  private :: GetSortedChannelIndex
  private :: GetSortedChannelNO2BIndex

  type :: OneBodyChannels
    integer, private :: emax = -1
    integer, private :: lmax = -1
    integer, private :: NChan = 0
    integer, private, allocatable :: j(:)
    integer, private, allocatable :: p(:)
    integer, private, allocatable :: jp2ch(:,:)
  contains
    procedure :: InitOneBodyChannels
    procedure :: FinOneBodyChannels
    procedure :: GetIndexOneBodyChannel
    procedure :: GetNumberOneBodyChannels
    procedure :: GetJOneBodyChannel
    procedure :: GetPOneBodyChannel
    generic :: init => InitOneBodyChannels
    generic :: fin => FinOneBodyChannels
    generic :: GetIndex => GetIndexOneBodyChannel
    generic :: GetNumberChannels => GetNumberOneBodyChannels
    generic :: GetJ => GetJOneBodyChannel
    generic :: GetP => GetPOneBodyChannel
  end type OneBodyChannels

  type :: TwoBodyChannels
    integer, private :: emax = -1
    integer, private :: lmax = -1
    integer, private :: e2max = -1
    integer, private :: NChan = 0
    integer, private, allocatable :: j(:)
    integer, private, allocatable :: p(:)
    integer, private, allocatable :: jp2ch(:,:)
  contains
    procedure :: InitTwoBodyChannels
    procedure :: FinTwoBodyChannels
    procedure :: GetNumberTwoBodyChannels
    procedure :: GetIndexTwoBodyChannel
    procedure :: GetJTwoBodyChannel
    procedure :: GetPTwoBodyChannel
    generic :: init => InitTwoBodyChannels
    generic :: fin => FinTwoBodyChannels
    generic :: GetIndex => GetIndexTwoBodyChannel
    generic :: GetNumberChannels => GetNumberTwoBodyChannels
    generic :: GetJ => GetJTwoBodyChannel
    generic :: GetP => GetPTwoBodyChannel
  end type TwoBodyChannels

  type :: ThreeBodyLabIsoChanNO2B
    type(OrbitsIsospin), pointer :: sps
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: t = -1
    integer, private :: j12 = -1
    integer, private :: p12 = -100
    integer, private :: j3 = -1
    integer, private :: p3 = -100
    integer, private :: n_states = 0
    integer, private, allocatable :: IndexABCT(:,:)
  contains
    procedure :: InitThreeBodyLabIsoChanNO2B
    procedure :: FinThreeBodyLabIsoChanNO2B
    procedure :: GetABCTNO2B
    procedure :: GetJNO2B
    procedure :: GetParityNO2B
    procedure :: GetTNO2B
    procedure :: GetJ12NO2B
    procedure :: GetP12NO2B
    procedure :: GetJ3NO2B
    procedure :: GetP3NO2B
    procedure :: GetNumberStatesNO2B
    generic :: init => InitThreeBodyLabIsoChanNO2B
    generic :: fin => FinThreeBodyLabIsoChanNO2B
    generic :: GetJ => GetJNO2B
    generic :: GetParity => GetParityNO2B
    generic :: GetT => GetTNO2B
    generic :: GetNumberStates => GetNumberStatesNO2B
    generic :: GetJ12 => GetJ12NO2B
    generic :: GetP12 => GetP12NO2B
    generic :: GetJ3 => GetJ3NO2B
    generic :: GetP3 => GetP3NO2B
    generic :: GetABCT => GetABCTNO2B
  end type ThreeBodyLabIsoChanNO2B

  type :: ThreeBodyLabIsoSpaceNO2B
    type(ThreeBodyLabIsoChanNO2B), allocatable :: Chan(:)
    type(OrbitsIsospin), pointer :: sps
    type(OneBodyChannels) :: one
    type(TwoBodyChannels) :: two
    integer, private, allocatable :: JPTNO2B2Ch(:,:)! ch_jpt, ch_no2b
    integer, private, allocatable :: JPT2Ch(:,:,:)  ! J, P, T
    integer, private, allocatable :: NO2B2Ch(:,:)   ! ch_two, ch_one
    integer, private, allocatable :: idx_sorted(:)
    integer, private :: emax = -1
    integer, private :: e2max= -1
    integer, private :: e3max= -1
    integer, private :: NChan= 0
    integer, private :: NNO2B= 0
    integer, private :: NJPT = 0
    logical :: is_Constructed=.false.
  contains
    procedure :: InitThreeBodyLabIsoSpaceNO2B
    procedure :: FinThreeBodyLabIsoSpaceNO2B
    procedure :: GetNumberChannelsNO2B
    procedure :: GetEmaxNO2B
    procedure :: GetE2maxNO2B
    procedure :: GetE3maxNO2B
    procedure :: GetIndexNO2B
    procedure :: GetChannelNO2BFromQNs
    procedure :: GetChannelNO2BFromIndex
    procedure :: GetSortedChannelNO2BIndex
    generic :: init => InitThreeBodyLabIsoSpaceNO2B
    generic :: fin => FinThreeBodyLabIsoSpaceNO2B
    generic :: GetNumberChannels => GetNumberChannelsNO2B
    generic :: GetEmax => GetEmaxNO2B
    generic :: GetE2max => GetE2maxNO2B
    generic :: GetE3max => GetE3maxNO2B
    generic :: GetIndex => GetIndexNO2B
    generic :: GetChannel => GetChannelNO2BFromQNs, GetChannelNO2BFromIndex
    generic :: GetSortedChannelIndex => GetSortedChannelNO2BIndex
  end type ThreeBodyLabIsoSpaceNO2B

  type :: NO2BThreeBodyIsoChan
    type(OrbitsIsospin), pointer :: sps
    integer, private :: t = -1
    integer, private :: j12 = -1
    integer, private :: p12 = -100
    integer, private :: j3  = -1
    integer, private :: p3  = -100
    integer, private :: n_states = 0
    integer, allocatable :: IndexABCT(:,:)
    !integer, allocatable :: ABNcT2Index(:,:,:,:) ! a, b, nc, tab => n
    !integer, allocatable :: iphase(:,:,:,:)    ! phase for a <-> b
    integer, allocatable :: TNcAB2Index(:,:,:,:) ! a, b, nc, tab => n
    integer, allocatable :: iphase(:,:,:,:)    ! phase for a <-> b
  contains
    procedure :: InitNO2BThreeBodyIsoChan
    procedure :: FinNO2BThreeBodyIsoChan
    procedure :: GetT
    procedure :: GetJ12
    procedure :: GetP12
    procedure :: GetJ3
    procedure :: GetP3
    procedure :: GetNumberStates
    procedure :: GetIndexChanFromQNs
    procedure :: GetABCT
    procedure :: GetPhase
    generic :: init => InitNO2BThreeBodyIsoChan
    generic :: fin => FinNO2BThreeBodyIsoChan
    generic :: GetIndex => GetIndexChanFromQNs
  end type NO2BThreeBodyIsoChan

  type :: NO2BThreeBodyIsoSpace
    type(NO2BThreeBodyIsoChan), allocatable :: Chan(:)
    type(OrbitsIsospin), pointer :: sps
    type(OneBodyChannels) :: one
    type(TwoBodyChannels) :: two
    integer, private, allocatable :: idx_sorted(:)
    integer, private, allocatable :: NO2B2Ch(:,:)   ! ch_ab (jab, pab), ch_c
    integer, private, allocatable :: T2Ch(:)        ! T_total
    integer, private, allocatable :: TNO2B2Ch(:,:) ! T, NO2B
    integer, private :: emax = -1
    integer, private :: e2max= -1
    integer, private :: e3max= -1
    integer, private :: NChan = 0
    integer, private :: NNO2B = 0
    integer, private :: NT = 0
    logical :: is_Constructed=.false.
  contains
    procedure :: InitNO2BThreeBodyIsoSpace
    procedure :: FinNO2BThreeBodyIsoSpace
    procedure :: GetNumberChannels
    procedure :: GetEmax
    procedure :: GetE2max
    procedure :: GetE3max
    procedure :: GetIndex
    procedure :: GetChannelFromQNs
    procedure :: GetChannelFromIndex
    procedure :: GetSortedChannelIndex
    generic :: init => InitNO2BThreeBodyIsoSpace
    generic :: fin => FinNO2BThreeBodyIsoSpace
    generic :: GetChannel => GetChannelFromQNs, GetChannelFromIndex
  end type NO2BThreeBodyIsoSpace
contains

  function GetNumberChannels(this) result(NChan)
    class(NO2BThreeBodyIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetEmax(this) result(Emax)
    class(NO2BThreeBodyIsoSpace), intent(in) :: this
    integer :: Emax
    Emax = this%emax
  end function GetEmax

  function GetE2max(this) result(E2max)
    class(NO2BThreeBodyIsoSpace), intent(in) :: this
    integer :: E2max
    E2max = this%e2max
  end function GetE2max

  function GetE3max(this) result(E3max)
    class(NO2BThreeBodyIsoSpace), intent(in) :: this
    integer :: E3max
    E3max = this%e3max
  end function GetE3max

  function GetIndex(this,t,j12,p12,j3,p3) result(idx)
    class(NO2BThreeBodyIsoSpace), intent(in) :: this
    integer, intent(in) :: t, j12, p12, j3, p3
    integer :: idx, idx_t, idx_no2b
    idx_t = this%T2Ch(t)
    idx_no2b = this%NO2B2Ch( this%two%GetIndex(j12,p12), this%one%GetIndex(j3,p3) )
    idx = this%TNO2B2Ch(idx_t, idx_no2b)
  end function GetIndex

  function GetChannelFromIndex(this,idx) result(r)
    class(NO2BThreeBodyIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    type(NO2BThreeBodyIsoChan), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%Chan(idx)
  end function GetChannelFromIndex

  function GetSortedChannelIndex(this,idx) result(sorted_idx)
    class(NO2BThreeBodyIsoSpace), intent(in), target :: this
    integer, intent(in) :: idx
    integer :: sorted_idx
    sorted_idx = 0
    if(idx == 0) return
    sorted_idx = this%idx_sorted(idx)
  end function GetSortedChannelIndex

  function GetChannelFromQNs(this,t,j12,p12,j3,p3) result(r)
    class(NO2BThreeBodyIsoSpace), intent(in), target :: this
    integer, intent(in) :: t,j12,p12,j3,p3
    type(NO2BThreeBodyIsoChan), pointer :: r
    r => this%GetChannel( this%GetIndex(t,j12,p12,j3,p3) )
  end function GetChannelFromQNs

  function GetT(this) result(T)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: T
    T = this%T
  end function GetT

  function GetJ12(this) result(j12)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: j12
    j12 = this%j12
  end function GetJ12

  function GetP12(this) result(p12)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: p12
    p12 = this%p12
  end function GetP12

  function GetJ3(this) result(j3)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: j3
    j3 = this%j3
  end function GetJ3

  function GetP3(this) result(p3)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: p3
    p3 = this%p3
  end function GetP3

  function GetNumberStates(this) result(n_states)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetABCT(this, idx) result(ABCT)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer, intent(in) :: idx
    integer :: ABCT(4)
    ABCT = this%IndexABCT(:,idx)
  end function GetABCT

  function GetIndexChanFromQNs(this,a,b,c,t) result(idx)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer, intent(in) :: a, b, c, t
    integer :: idx
    type(SingleParticleOrbitIsospin), pointer :: oc
    oc => this%sps%GetOrbit(c)
    !idx = this%ABNcT2Index(a,b,oc%n,t)
    idx = this%TNcAB2Index(t,oc%n,a,b)
  end function GetIndexChanFromQNs

  function GetPhase(this,a,b,c,t) result(phase)
    class(NO2BThreeBodyIsoChan), intent(in) :: this
    integer, intent(in) :: a, b, c, t
    real(8) :: phase
    type(SingleParticleOrbitIsospin), pointer :: oc
    oc => this%sps%GetOrbit(c)
    !phase = dble( this%iphase(a,b,oc%n,t) )
    phase = dble( this%iphase(t,oc%n,a,b) )
  end function GetPhase

  function GetNumberChannelsNO2B(this) result(NChan)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannelsNO2B

  function GetEmaxNO2B(this) result(Emax)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in) :: this
    integer :: Emax
    Emax = this%Emax
  end function GetEmaxNO2B

  function GetE2maxNO2B(this) result(E2max)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in) :: this
    integer :: E2max
    E2max = this%E2max
  end function GetE2maxNO2B

  function GetE3maxNO2B(this) result(E3max)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in) :: this
    integer :: E3max
    E3max = this%E3max
  end function GetE3maxNO2B

  function GetIndexNO2B(this, j, p, t, j12, p12, j3, p3) result(idx)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in) :: this
    integer, intent(in) :: j, p, t, j12, p12, j3, p3
    integer :: idx_jpt, idx_no2b, idx
    idx_jpt = this%JPT2Ch(j,p,t)
    idx_no2b = this%NO2B2Ch(this%two%GetIndex(j12,p12), this%one%GetIndex(j3,p3))
    idx = this%JPTNO2B2Ch(idx_jpt, idx_no2b)
  end function GetIndexNO2B

  function GetChannelNO2BFromIndex(this, idx) result(r)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in), target :: this
    integer, intent(in) :: idx
    type(ThreeBodyLabIsoChanNO2B), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%Chan(idx)
  end function GetChannelNO2BFromIndex

  function GetSortedChannelNO2BIndex(this, idx) result(sorted_idx)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in), target :: this
    integer, intent(in) :: idx
    integer :: sorted_idx
    sorted_idx = 0
    if(idx == 0) return
    sorted_idx = this%idx_sorted(idx)
  end function GetSortedChannelNO2BIndex

  function GetChannelNO2BFromQNs(this, j, p, t, j12, p12, j3, p3) result(r)
    class(ThreeBodyLabIsoSpaceNO2B), intent(in), target :: this
    integer, intent(in) :: j, p, t, j12, p12, j3, p3
    type(ThreeBodyLabIsoChanNO2B), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,t,j12,p12,j3,p3) )
  end function GetChannelNO2BFromQNs

  function GetABCTNO2B(this, idx) result(ABCT)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer, intent(in) :: idx
    integer :: ABCT(4)
    ABCT = this%IndexABCT(:,idx)
  end function GetABCTNO2B

  function GetJNO2B(this) result(J)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: J
    J = this%J
  end function GetJNO2B

  function GetParityNO2B(this) result(P)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: P
    P = this%P
  end function GetParityNO2B

  function GetTNO2B(this) result(T)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: T
    T = this%T
  end function GetTNO2B

  function GetJ12NO2B(this) result(J12)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: J12
    J12 = this%J12
  end function GetJ12NO2B

  function GetP12NO2B(this) result(P12)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: P12
    P12 = this%P12
  end function GetP12NO2B

  function GetJ3NO2B(this) result(J)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: J
    J = this%J3
  end function GetJ3NO2B

  function GetP3NO2B(this) result(P)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: P
    P = this%P3
  end function GetP3NO2B

  function GetNumberStatesNO2B(this) result(n_states)
    class(ThreeBodyLabIsoChanNO2B), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStatesNO2B

  function GetIndexOneBodyChannel(this, j, p) result(idx)
    class(OneBodyChannels), intent(in) :: this
    integer, intent(in) :: j, p
    integer :: idx
    idx = this%jp2ch(j,p)
  end function GetIndexOneBodyChannel

  function GetJOneBodyChannel(this, idx) result(j)
    class(OneBodyChannels), intent(in) :: this
    integer, intent(in) :: idx
    integer :: j
    j = this%j(idx)
  end function GetJOneBodyChannel

  function GetPOneBodyChannel(this, idx) result(p)
    class(OneBodyChannels), intent(in) :: this
    integer, intent(in) :: idx
    integer :: p
    p = this%p(idx)
  end function GetPOneBodyChannel

  function GetNumberOneBodyChannels(this) result(NChan)
    class(OneBodyChannels), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberOneBodyChannels

  function GetNumberTwoBodyChannels(this) result(NChan)
    class(TwoBodyChannels), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberTwoBodyChannels

  function GetIndexTwoBodyChannel(this, j, p) result(idx)
    class(TwoBodyChannels), intent(in) :: this
    integer, intent(in) :: j, p
    integer :: idx
    idx = this%jp2ch(j,p)
  end function GetIndexTwoBodyChannel

  function GetJTwoBodyChannel(this, idx) result(j)
    class(TwoBodyChannels), intent(in) :: this
    integer, intent(in) :: idx
    integer :: j
    j = this%j(idx)
  end function GetJTwoBodyChannel

  function GetPTwoBodyChannel(this, idx) result(p)
    class(TwoBodyChannels), intent(in) :: this
    integer, intent(in) :: idx
    integer :: p
    p = this%p(idx)
  end function GetPTwoBodyChannel

  subroutine FinOneBodyChannels(this)
    class(OneBodyChannels), intent(inout) :: this
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%jp2ch)
    this%emax = -1
    this%NChan = 0
  end subroutine FinOneBodyChannels

  subroutine InitOneBodyChannels(this, sps)
    class(OneBodyChannels), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer :: la, ja, cnt
    this%emax = sps%emax
    this%lmax = sps%lmax
    allocate(this%jp2ch(2*sps%emax+1,-1:1))
    cnt = 0
    do la = 0, sps%lmax
      do ja = abs(2*la-1), 2*la+1, 2
        cnt = cnt + 1
      end do
    end do
    this%NChan = cnt
    allocate(this%j(this%NChan))
    allocate(this%p(this%NChan))

    cnt = 0
    do la = 0, sps%lmax
      do ja = abs(2*la-1), 2*la+1, 2
        cnt = cnt + 1
        this%jp2ch(ja,(-1)**la) = cnt
        this%j(cnt) = ja
        this%p(cnt) = (-1)**la
      end do
    end do
  end subroutine InitOneBodyChannels

  subroutine FinTwoBodyChannels(this)
    class(TwoBodyChannels), intent(inout) :: this
    deallocate(this%j)
    deallocate(this%p)
    deallocate(this%jp2ch)
    this%emax = -1
    this%e2max = -1
    this%NChan = 0
  end subroutine FinTwoBodyChannels

  subroutine InitTwoBodyChannels(this, sps, e2max)
    class(TwoBodyChannels), intent(inout) :: this
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: j, p, ich, jmax
    this%lmax = sps%lmax
    this%emax = sps%emax
    this%e2max = e2max
    jmax = min(2*sps%lmax, e2max)+1
    this%NChan = (jmax+1) * 2

    allocate(this%j(this%NChan))
    allocate(this%p(this%NChan))
    allocate(this%jp2ch(0:jmax,-1:1))
    this%jp2ch(:,:) = 0

    ich = 0
    do j = 0, jmax
      do p = 1, -1, -2
        ich = ich + 1
        this%j(ich) = j
        this%p(ich) = p
        this%jp2ch(j,p) = ich
      end do
    end do
  end subroutine InitTwoBodyChannels

  subroutine FinThreeBodyLabIsoChanNO2B(this)
    class(ThreeBodyLabIsoChanNO2B), intent(inout) :: this
    if(this%GetNumberStates() == 0) return
    deallocate(this%IndexABCT)
    this%sps => null()
    this%j = -1
    this%p = -100
    this%t = -1
    this%j12 = -1
    this%p12 = -100
    this%j3 = -1
    this%p3 = -100
    this%n_states = 0
  end subroutine FinThreeBodyLabIsoChanNO2B

  subroutine InitThreeBodyLabIsoChanNO2B(this, sps, j, p, t, j12, p12, j3, p3, e2max, e3max, n_states, set_mode)
    use MyLibrary, only: triag
    class(ThreeBodyLabIsoChanNO2B), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: j, p, t, j12, p12, j3, p3, e2max, e3max
    integer, intent(inout) :: n_states
    logical, intent(in) :: set_mode
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: i1, i2, i3, t12

    this%j = j
    this%p = p
    this%t = t
    this%j12 = j12
    this%p12 = p12
    this%j3 = j3
    this%p3 = p3
    this%sps => sps

    if( set_mode ) then
      this%n_states = n_states
      allocate( this%IndexABCT(4, this%GetNumberStates() ) )
    end if

    n_states = 0
    do i1 = 1, sps%norbs
      o1 => sps%GetOrbit(i1)
      do i2 = 1, i1
        o2 => sps%GetOrbit(i2)
        if(triag(o1%j,o2%j,2*j12)) cycle
        if((-1)**(o1%l+o2%l) /= p12) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, sps%norbs
          o3 => sps%GetOrbit(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do t12 = 0, 1
            if(triag(2*t12, 1, t)) cycle
            if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
            n_states = n_states + 1
            if(set_mode) this%IndexABCT(:,n_states) = [i1,i2,i3,t12]
          end do
        end do
      end do
    end do
  end subroutine InitThreeBodyLabIsoChanNO2B

  subroutine FinThreeBodyLabIsoSpaceNO2B(this)
    class(ThreeBodyLabIsoSpaceNO2B), intent(inout) :: this
    integer :: ch
    if(this%GetNumberChannels() == 0) return
    do ch = 1, this%GetNumberChannels()
      call this%Chan(ch)%fin()
    end do
    deallocate(this%Chan)
    call this%one%fin()
    call this%two%fin()
    this%sps => null()
    deallocate( this%JPTNO2B2Ch )
    deallocate( this%JPT2Ch )
    deallocate( this%NO2B2Ch )
    this%emax = -1
    this%e2max= -1
    this%e3max= -1
    this%NChan= 0
    this%NNO2B= 0
    this%NJPT = 0
    this%is_Constructed = .false.
  end subroutine FinThreeBodyLabIsoSpaceNO2B

  subroutine InitThreeBodyLabIsoSpaceNO2B(this, sps, e2max, e3max)
    use MyLibrary, only: triag
    use Profiler
    class(ThreeBodyLabIsoSpaceNO2B), intent(inout) :: this
    type(OrbitsIsospin), intent(in), target :: sps
    integer, intent(in) :: e2max, e3max
    integer :: ttot, jtot, ptot, ch12, j12, p12, ch3, j3, p3, n_states
    integer :: ich, ich_JPT, ich_NO2B, idx
    integer, allocatable :: jlist(:), plist(:), tlist(:)
    real(8) :: time
    type(ThreeBodyLabIsoChanNO2B) :: temp
    integer, allocatable :: ndim_array(:)
    type(sys) :: s

    time = omp_get_wtime()
    this%is_Constructed = .true.
    this%sps => sps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max
    call this%one%init(sps)
    call this%two%init(sps,e2max)

    allocate(jlist(4*(e3max+2)))
    allocate(plist(4*(e3max+2)))
    allocate(tlist(4*(e3max+2)))
    jlist(:) = 0; plist(:) = 0; tlist(:) = 0
    ich_JPT = 0
    do ttot = 1, 3, 2
      do jtot = 1, 2*e3max+3, 2
        do ptot = 1, -1, -2
           ich_JPT = ich_JPT + 1
           jlist(ich_JPT) = jtot
           plist(ich_JPT) = ptot
           tlist(ich_JPT) = ttot
        end do
      end do
    end do

    this%NJPT = (2*e3max+3) * 4
    this%NNO2B = this%two%NChan * this%one%NChan
    allocate(this%JPTNO2B2Ch(this%NJPT, this%NNO2B))
    allocate(this%JPT2Ch(2 * e3max + 3, -1:1, 1:3))
    allocate(this%NO2B2Ch(this%two%NChan, this%one%NChan ))
    this%JPTNO2B2Ch(:,:) = 0
    this%JPT2Ch(:,:,:) = 0
    this%NO2B2Ch(:,:) = 0
    ich = 0
    do ich_JPT = 1, size(jlist)
      jtot = jlist(ich_JPT)
      ptot = plist(ich_JPT)
      ttot = tlist(ich_JPT)
      this%JPT2Ch(jtot,ptot,ttot) = ich_JPT

      ich_NO2B = 0
      do ch12 = 1, this%two%NChan
        j12 = this%two%j(ch12)
        p12 = this%two%p(ch12)
        do ch3 = 1, this%one%NChan
          j3 = this%one%j(ch3)
          p3 = this%one%p(ch3)
          ich_NO2B = ich_NO2B + 1
          this%NO2B2Ch(ch12,ch3) = ich_NO2B

          if(triag(2*j12, j3, jtot)) cycle
          if(p12 * p3 /= ptot) cycle
          call temp%init(sps, jtot, ptot, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .false.)
          if(n_states == 0) cycle
          ich = ich + 1
          this%JPTNO2B2Ch(ich_JPT,ich_NO2B) = ich
        end do
      end do
    end do
    this%NChan = ich
    allocate(this%Chan( this%GetNumberChannels() ))
    allocate(this%idx_sorted( this%GetNumberChannels() ))

    ich = 0
    do ich_JPT = 1, size(jlist)
      jtot = jlist(ich_JPT)
      ptot = plist(ich_JPT)
      ttot = tlist(ich_JPT)
      this%JPT2Ch(jtot,ptot,ttot) = ich_JPT

      ich_NO2B = 0
      do ch12 = 1, this%two%NChan
        j12 = this%two%j(ch12)
        p12 = this%two%p(ch12)
        do ch3 = 1, this%one%NChan
          j3 = this%one%j(ch3)
          p3 = this%one%p(ch3)
          ich_NO2B = ich_NO2B + 1
          this%NO2B2Ch(ch12,ch3) = ich_NO2B

          if(triag(2*j12, j3, jtot)) cycle
          if(p12 * p3 /= ptot) cycle
          call temp%init(sps, jtot, ptot, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .false.)
          if(n_states == 0) cycle
          ich = ich + 1
          call this%Chan(ich)%init(sps, jtot, ptot, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .true.)
        end do
      end do
    end do

    deallocate(jlist, plist, tlist)
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
    !do ich = 1, this%GetNumberChannels()
    !  write(*,'(3i12)') ich, this%Chan(ich)%GetNumberStates(), this%Chan( this%idx_sorted(ich) )%GetNumberStates()
    !end do
    call timer%add(s%str("InitThreeBodyIsoSpaceNO2B"), omp_get_wtime()-time)
  end subroutine InitThreeBodyLabIsoSpaceNO2B

  subroutine FinNO2BThreeBodyIsoChan(this)
    class(NO2BThreeBodyIsoChan), intent(inout) :: this
    if(this%GetNumberStates() == 0) return
    deallocate(this%IndexABCT)
    deallocate(this%TNcAB2Index)
    !deallocate(this%ABNcT2Index)
    deallocate(this%iphase)
    this%sps => null()
    this%t = -1
    this%j12 = -1
    this%p12 = -100
    this%j3  = -1
    this%p3  = -100
    this%n_states = 0
  end subroutine FinNO2BThreeBodyIsoChan

  subroutine InitNO2BThreeBodyIsoChan(this, sps, t, j12, p12, j3, p3, e2max, e3max, n_states, set_mode)
    use MyLibrary, only: triag
    class(NO2BThreeBodyIsoChan), intent(inout) :: this
    type(OrbitsIsospin), target, intent(in) :: sps
    integer, intent(in) :: t, j12, p12, j3, p3, e2max, e3max
    integer, intent(inout) :: n_states
    logical, intent(in) :: set_mode
    type(SingleParticleOrbitIsospin), pointer :: o1,o2,o3
    integer :: i1,i2,i3,t12

    this%t = t
    this%j12 = j12
    this%p12 = p12
    this%j3 = j3
    this%p3 = p3
    this%sps => sps
    if( set_mode ) then
      this%n_states = n_states
      allocate( this%IndexABCT(4, this%GetNumberStates() ) )
      !allocate( this%ABNcT2Index( sps%norbs, sps%norbs, 0:sps%emax/2, 0:1) )
      !allocate( this%iphase(      sps%norbs, sps%norbs, 0:sps%emax/2, 0:1) )
      allocate( this%TNcAB2Index(  0:1, 0:sps%emax/2, sps%norbs, sps%norbs) )
      allocate( this%iphase(       0:1, 0:sps%emax/2, sps%norbs, sps%norbs) )
      !this%ABNcT2Index(:,:,:,:) = 0
      this%TNcAB2Index(:,:,:,:) = 0
      this%iphase(     :,:,:,:) = 0
    end if
    n_states = 0
    do i1 = 1, sps%norbs
      o1 => sps%GetOrbit(i1)
      do i2 = 1, i1
        o2 => sps%GetOrbit(i2)
        if(triag(o1%j,o2%j,2*j12)) cycle
        if((-1)**(o1%l+o2%l) /= p12) cycle
        if(o1%e + o2%e > e2max) cycle
        do i3 = 1, sps%norbs
          o3 => sps%GetOrbit(i3)
          if(j3 /= o3%j) cycle
          if(p3 /= (-1)**o3%l) cycle
          if(o1%e + o3%e > e2max) cycle
          if(o2%e + o3%e > e2max) cycle
          if(o1%e + o2%e + o3%e > e3max) cycle
          do t12 = 0, 1
            if(triag(2*t12, 1, t)) cycle
            if(i1 == i2 .and. mod(j12 + t12, 2) == 0) cycle
            n_states = n_states + 1
            if(set_mode) then
              this%IndexABCT(:,n_states) = [i1,i2,i3,t12]
              !this%ABNcT2Index(i1,i2,o3%n,t12) = n_states
              !this%ABNcT2Index(i2,i1,o3%n,t12) = n_states
              !this%iphase(i1,i2,o3%n,t12) = 1
              !this%iphase(i2,i1,o3%n,t12) = (-1)**((o1%j+o2%j)/2+j12+t12)
              this%TNcAB2Index(t12,o3%n,i1,i2) = n_states
              this%TNcAB2Index(t12,o3%n,i2,i1) = n_states
              this%iphase(t12,o3%n,i1,i2) = 1
              this%iphase(t12,o3%n,i2,i1) = (-1)**((o1%j+o2%j)/2+j12+t12)
            end if
          end do
        end do
      end do
    end do
  end subroutine InitNO2BThreeBodyIsoChan

  subroutine FinNO2BThreeBodyIsoSpace(this)
    class(NO2BThreeBodyIsoSpace), intent(inout) :: this
    integer :: ch
    if(this%GetNumberChannels() == 0) return
    do ch = 1, this%GetNumberChannels()
      call this%Chan(ch)%fin()
    end do
    deallocate(this%Chan)
    call this%one%fin()
    call this%two%fin()
    this%sps => null()
    deallocate( this%NO2B2Ch )
    deallocate( this%T2Ch )
    deallocate( this%TNO2B2Ch )
    this%emax = -1
    this%e2max= -1
    this%e3max= -1
    this%NChan = 0
    this%NNO2B = 0
    this%NT = 0
    this%is_Constructed=.false.
  end subroutine FinNO2BThreeBodyIsoSpace

  subroutine InitNO2BThreeBodyIsoSpace(this, sps, e2max, e3max)
    use Profiler
    class(NO2BThreeBodyIsoSpace), intent(inout) :: this
    type(OrbitsIsospin), intent(in), target :: sps
    integer, intent(in) :: e2max, e3max
    integer :: ch12, ch3, ich, ich_NO2B, ich_T, j12, p12, j3, p3, n_states, Ttot, idx
    type(NO2BThreeBodyIsoChan) :: temp
    real(8) :: time, mem
    integer, allocatable :: ndim_array(:)
    type(sys) :: s

    time = omp_get_wtime()
    this%is_Constructed = .true.
    this%sps => sps
    this%emax = sps%emax
    this%e2max = e2max
    this%e3max = e3max
    call this%one%init(sps)
    call this%two%init(sps,e2max)
    this%NT = 2
    this%NNO2B = this%one%NChan * this%two%NChan
    allocate(this%TNO2B2Ch(this%NT, this%NNO2B))
    allocate(this%T2Ch(1:3))
    allocate(this%NO2B2Ch(this%two%NChan, this%one%NChan))
    this%TNO2B2Ch(:,:) = 0
    this%T2Ch(:) = 0
    this%NO2B2Ch(:,:) = 0
    ich = 0
    ich_T = 0
    do Ttot = 1, 3, 2
      ich_T = ich_T + 1
      this%T2Ch(Ttot) = ich_T

      ich_NO2B = 0
      do ch12 = 1, this%two%NChan
        j12 = this%two%j(ch12)
        p12 = this%two%p(ch12)
        do ch3 = 1, this%one%NChan
          j3 = this%one%j(ch3)
          p3 = this%one%p(ch3)
          ich_NO2B = ich_NO2B + 1
          this%NO2B2Ch(ch12,ch3) = ich_NO2B
          call temp%init(sps, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .false.)
          if(n_states == 0) cycle
          ich = ich + 1
          this%TNO2B2Ch(ich_T,ich_NO2B) = ich

        end do
      end do
    end do
    this%NChan = ich
    allocate(this%Chan( this%GetNumberChannels() ))
    allocate(this%idx_sorted( this%GetNumberChannels() ))
    ich = 0
    ich_T = 0
    do Ttot = 1, 3, 2
      ich_T = ich_T + 1
      this%T2Ch(Ttot) = ich_T

      ich_NO2B = 0
      do ch12 = 1, this%two%NChan
        j12 = this%two%j(ch12)
        p12 = this%two%p(ch12)
        do ch3 = 1, this%one%NChan
          j3 = this%one%j(ch3)
          p3 = this%one%p(ch3)
          ich_NO2B = ich_NO2B + 1
          this%NO2B2Ch(ch12,ch3) = ich_NO2B
          call temp%init(sps, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .false.)
          if(n_states == 0) cycle
          ich = ich + 1
          call this%Chan(ich)%init(sps, ttot, j12, p12, j3, p3, e2max, e3max, n_states, .true.)

        end do
      end do
    end do

    allocate(ndim_array(this%GetNumberChannels()))
    do ich = 1, this%GetNumberChannels()
      ndim_array(ich) = this%Chan(ich)%GetNumberStates()
    end do

    do ich = 1, this%GetNumberChannels()
      idx = maxloc(ndim_array,dim=1)
      this%idx_sorted(ich) = idx
      ndim_array(idx) = -1
    end do
    deallocate(ndim_array)

    !do ich = 1, this%GetNumberChannels()
    !  write(*,'(3i12)') ich, this%Chan(ich)%GetNumberStates(), this%Chan( this%idx_sorted(ich) )%GetNumberStates()
    !end do
    mem = dble(this%GetNumberChannels()) * dble(sps%norbs) * dble(sps%norbs) * &
        & dble(sps%emax/2+1) * 16.d0 / 1024.d0**3
    if(myrank == 0) write(*,"(a, f10.4,a)") "# Memory used for NO2BThreeBodyIsoSpace is ", mem, " GB"
    call timer%add(s%str("InitNO2BThreeBodyIsoSpace"), omp_get_wtime()-time)
  end subroutine InitNO2BThreeBodyIsoSpace

end module ThreeBodyLabSpaceIsoNO2B

