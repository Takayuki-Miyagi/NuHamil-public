module TwoBodyTransCoefIso
  use omp_lib
  use ClassSys
  use SingleParticleState
  use TwoBodyLabSpaceIso
  use StoreCouplings
  implicit none

  public :: TransRel2LabSpaceIso
  !public :: TransRel2LabChannelIso
  !public :: TCoefsIso
  !public :: RelativeCMQNsIso
  private

  type :: RelativeCMQNs
    integer :: Nmax = -1
    integer :: Nrel = -1
    integer :: Lrel = -1
    integer :: s    = -1
    integer :: Jrel = -1
    integer :: Ncm  = -1
    integer :: Lcm  = -1
    integer :: idx  = -1
  contains
    procedure :: SetRelativeCMQNs
    generic :: set => SetRelativeCMQNs
  end type RelativeCMQNs

  type :: TwoBodyLabQNs
    type(SingleParticleOrbitIsospin), pointer :: o1, o2
    integer :: Nmax, idx
  contains
    procedure :: SetTwoBodyLabQNs
    generic :: set => SetTwoBodyLabQNs
  end type TwoBodyLabQNs

  type :: TCoefsIso
    real(8), allocatable :: mat(:,:)
    integer, private, allocatable :: RelCMIndex(:)
    integer, private, allocatable :: LabIndex(:)
    integer, private, allocatable :: RelCMQNs2Idx(:,:,:,:,:,:) ! ncm, lcm, nrel, lrel, s, jrel
    integer, private, allocatable :: LabIndex2Idx(:,:) ! (a,b)
    integer, private, allocatable :: LabIndex2Phase(:,:) ! (a,b)
    integer, private :: NumberRelCMStates = 0
    integer, private :: NumberLabStates = 0
  contains
    procedure :: InitTCoefs
    procedure :: FinTCoefs
    procedure :: GetNumberRelCMStatesTCoefs
    procedure :: GetNumberLabStatesTCoefs
    procedure :: GetRelCMIndexTCoefs
    procedure :: GetLabIndexTCoefs
    procedure :: GetMemoryTCoefs
    generic :: init => InitTCoefs
    generic :: fin => FinTCoefs
    generic :: GetNumberRelCMStates => GetNumberRelCMStatesTCoefs
    generic :: GetNumberLabStates => GetNumberLabStatesTCoefs
    generic :: GetRelCMIndex => GetRelCMIndexTCoefs
    generic :: GetLabIndex => GetLabIndexTCoefs
    generic :: GetMemory => GetMemoryTCoefs
  end type TCoefsIso

  type :: TransRel2LabChannelIso
    type(OrbitsIsospin), pointer :: sps
    type(TCoefsIso), allocatable :: N2max(:)
    type(RelativeCMQNs), allocatable :: RelCM(:)
    type(TwoBodyLabQNs), allocatable :: Lab(:)
    integer, private, pointer :: labels2n(:,:) => null()
    integer, private, pointer :: iphase(:,:) => null()
    integer, private, pointer :: Index2a(:) => null()
    integer, private, pointer :: Index2b(:) => null()
    integer, private :: NumberRelCMStates = 0
    integer, private :: NumberLabStates = 0
    integer, private :: Nmax = -1
    integer, private :: J = -1
    integer, private :: P = -100
    integer, private :: T = -100
  contains
    procedure :: InitTransRel2LabChannel
    procedure :: FinTransRel2LabChannel
    procedure :: GetNmaxChannel
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetRelCMFromIndex
    procedure :: GetRelCMFromNIndex
    procedure :: GetLabFromIndex
    procedure :: GetLabFromNIndex
    procedure :: GetNumberLabStates
    procedure :: GetNumberRelCMStates
    procedure :: GetMemoryChannel
    generic :: init => InitTransRel2LabChannel
    generic :: fin => FinTransRel2LabChannel
    generic :: GetRelCM => GetRelCMFromIndex, GetRelCMFromNIndex
    generic :: GetLab => GetLabFromIndex, GetLabFromNIndex
    generic :: GetNmax => GetNmaxChannel
    generic :: GetMemory => GetMemoryChannel
  end type TransRel2LabChannelIso

  type :: TransRel2LabSpaceIso
    type(OrbitsIsospin), pointer :: sps
    type(TransRel2LabChannelIso), allocatable :: jpt(:)
    integer, private, allocatable :: jpt2idx(:,:,:)
    integer :: NChan = 0
    integer :: Nmax = -1
    integer :: Jmax = -1
    logical :: is_Constructed=.False.
  contains
    procedure :: GetNmax
    procedure :: GetJmax
    procedure :: GetIndex
    procedure :: GetNumberChannels
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromJPT
    procedure :: GetMemory
    procedure :: GetMatrixElement
    procedure :: init => InitTransRel2LabSpace
    procedure :: fin => FinTransRel2LabSpace
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromJPT
    generic :: GetME => GetMatrixElement
  end type TransRel2LabSpaceIso

  type(SixJsStore), private :: sixjs
  type(NineJsStore), private :: ninejs
  type(TMbracketStore), private :: TMbks

contains

  function GetNumberRelCMStatesTcoefs(this) result(n_states)
    class(TCoefsIso), intent(in) :: this
    integer :: n_states
    n_states = this%NumberRelCMStates
  end function GetNumberRelCMStatesTCoefs

  function GetNumberLabStatesTcoefs(this) result(n_states)
    class(TCoefsIso), intent(in) :: this
    integer :: n_states
    n_states = this%NumberLabStates
  end function GetNumberLabStatesTCoefs

  function GetRelCMIndexTCoefs(this,i) result(idx)
    class(TCoefsIso), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%RelCMIndex(i)
  end function GetRelCMIndexTCoefs

  function GetLabIndexTCoefs(this,i) result(idx)
    class(TCoefsIso), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%LabIndex(i)
  end function GetLabIndexTCoefs

  !
  ! TransRel2LabChannel
  !
  function GetRelCMFromIndex(this,i) result(r)
    class(TransRel2LabChannelIso), intent(in), target :: this
    integer, intent(in) :: i
    type(RelativeCMQNs), pointer :: r
    r => null()
    if(i==0) return
    r => this%RelCM(i)
  end function GetRelCMFromIndex

  function GetRelCMFromNIndex(this,N,i) result(r)
    class(TransRel2LabChannelIso), intent(in), target :: this
    integer, intent(in) :: N,i
    type(RelativeCMQNs), pointer :: r
    r => this%RelCM( this%N2max(N)%GetRelCMIndex(i) )
  end function GetRelCMFromNIndex

  function GetLabFromIndex(this,i) result(r)
    class(TransRel2LabChannelIso), intent(in), target :: this
    integer, intent(in) :: i
    type(TwoBodyLabQNs), pointer :: r
    r => null()
    if(i==0) return
    r => this%Lab(i)
  end function GetLabFromIndex

  function GetLabFromNIndex(this,N,i) result(r)
    class(TransRel2LabChannelIso), intent(in), target :: this
    integer, intent(in) :: N,i
    type(TwoBodyLabQNs), pointer :: r
    r => this%Lab( this%N2max(N)%GetLabIndex(i) )
  end function GetLabFromNIndex

  function GetNumberRelCMStates(this) result(n_states)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: n_states
    n_states = this%NumberRelCMStates
  end function GetNumberRelCMStates

  function GetNumberLabStates(this) result(n_states)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: n_states
    n_states = this%NumberLabStates
  end function GetNumberLabStates

  function GetJ(this) result(J)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: J
    J = this%J
  end function GetJ

  function GetParity(this) result(P)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: P
    P = this%P
  end function GetParity

  function GetT(this) result(T)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: T
    T = this%T
  end function GetT

  function GetNmaxChannel(this) result(Nmax)
    class(TransRel2LabChannelIso), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmaxChannel
  !
  !
  !

  !
  ! TransRel2LabSpace
  !
  function GetNmax(this) result(Nmax)
    class(TransRel2LabSpaceIso), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TransRel2LabSpaceIso), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  function GetNumberChannels(this) result(NChan)
    class(TransRel2LabSpaceIso), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, t) result(idx)
    class(TransRel2LabSpaceIso), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j,p,t)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(r)
    class(TransRel2LabSpaceIso), intent(in), target :: this
    integer, intent(in) :: idx
    type(TransRel2LabChannelIso), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%jpt(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this, j, p, t) result(r)
    class(TransRel2LabSpaceIso), intent(in), target :: this
    integer, intent(in) :: j, p, t
    type(TransRel2LabChannelIso), pointer :: r
    r => this%jpt( this%GetIndex(j, p, t) )
  end function GetChannelFromJPT
  !
  !
  !

  subroutine FinTCoefs(this)
    class(TCoefsIso), intent(inout) :: this
    if(this%GetNumberLabStates() * this%GetNumberRelCMStates() == 0) return
    deallocate(this%mat)
    deallocate(this%RelCMIndex)
    deallocate(this%LabIndex)
    deallocate(this%LabIndex2Idx)
    deallocate(this%LabIndex2Phase)
    deallocate(this%RelCMQNs2Idx)
    this%NumberLabStates = 0
    this%NumberRelCMStates = 0
  end subroutine FinTCoefs

  subroutine InitTCoefs(this, chan, Nmax)
    use MyLibrary, only: triag
    class(TCoefsIso), intent(inout) :: this
    type(TransRel2LabChannelIso), intent(in) :: chan
    integer, intent(in) :: Nmax
    integer :: idx, cnt, i_relcm, i_lab
    type(SingleParticleOrbitIsospin), pointer :: oa, ob
    integer :: a, na, la, ja, b, nb, lb, jb
    integer :: nrel, lrel, ncm, lcm, s, jrel
    integer :: lam, lammin, lammax
    type(RelativeCMQNs), pointer :: relcmqn
    type(TwoBodyLabQNs), pointer :: labqn
    real(8) :: delab, ft, tc

    this%NumberLabStates = 0
    cnt = 0
    do idx = 1, chan%GetNumberLabStates()
      labqn => chan%GetLab(idx)
      if(labqn%Nmax /= Nmax ) cycle
      cnt = cnt + 1
    end do
    this%NumberLabStates = cnt

    this%NumberRelCMStates = 0
    cnt = 0
    do idx = 1, chan%GetNumberRelCMStates()
      relcmqn => chan%GetRelCM(idx)
      if(relcmqn%Nmax /= Nmax ) cycle
      cnt = cnt + 1
    end do
    this%NumberRelCMStates = cnt
    if(this%GetNumberLabStates() * this%GetNumberRelCMStates() == 0) return
    allocate(this%mat( this%GetNumberRelCMStates(), this%GetNumberLabStates() ))
    allocate(this%RelCMIndex( this%GetNumberRelCMStates() ))
    allocate(this%LabIndex( this%GetNumberLabStates() ))
    this%mat(:,:) = 0.d0
    this%RelCMIndex(:) = 0
    this%LabIndex(:) = 0

    cnt = 0
    do idx = 1, chan%GetNumberLabStates()
      labqn => chan%GetLab(idx)
      if(labqn%Nmax /= Nmax ) cycle
      cnt = cnt + 1
      labqn%idx = cnt
      this%LabIndex(cnt) = idx
    end do

    cnt = 0
    do idx = 1, chan%GetNumberRelCMStates()
      relcmqn => chan%GetRelCM(idx)
      if(relcmqn%Nmax /= Nmax ) cycle
      cnt = cnt + 1
      relcmqn%idx = cnt
      this%RelCMIndex(cnt) = idx
    end do

    allocate(this%RelCMQNs2Idx(0:Nmax/2, 0:Nmax, 0:Nmax/2, 0:Nmax, 0:1, 0:Nmax+1))
    allocate(this%LabIndex2Idx(chan%sps%norbs, chan%sps%norbs))
    allocate(this%LabIndex2Phase(chan%sps%norbs, chan%sps%norbs))
    this%RelCMQNs2Idx(:,:,:,:,:,:) = 0
    this%LabIndex2Idx(:,:) = 0
    this%LabIndex2Phase(:,:) = 0
    do i_lab = 1, this%GetNumberLabStates()
      labqn => chan%GetLab( this%GetLabIndex(i_lab) )
      a = chan%Index2a( this%GetLabIndex(i_lab) )
      b = chan%Index2b( this%GetLabIndex(i_lab) )
      oa => labqn%o1
      ob => labqn%o2
      this%LabIndex2Idx(a,b) = i_lab
      this%LabIndex2Idx(b,a) = i_lab
      this%LabIndex2Phase(a,b) = 1
      this%LabIndex2Phase(b,a) = (-1)**((oa%j + ob%j)/2 - chan%GetJ() - chan%GetT())
    end do

    do i_relcm = 1, this%GetNumberRelCMStates()
      relcmqn => chan%GetRelCM( this%RelCMIndex(i_relcm) )
      nrel = relcmqn%nrel
      lrel = relcmqn%lrel
      ncm  = relcmqn%ncm
      lcm  = relcmqn%lcm
      s    = relcmqn%s
      jrel = relcmqn%Jrel
      this%RelCMQNs2Idx(ncm, lcm, nrel, lrel, s, jrel) = i_relcm
    end do

    !$omp parallel
    !$omp do private(i_lab, labqn, a, b, oa, ob, na, la, ja, nb, lb, jb, delab, &
    !$omp &  i_relcm, relcmqn, nrel, lrel, ncm, lcm, s, jrel, ft, lammin, lammax, &
    !$omp &  tc, lam) schedule(dynamic)
    do i_lab = 1, this%GetNumberLabStates()
      labqn => chan%GetLab( this%GetLabIndex(i_lab) )
      a = chan%Index2a( this%GetLabIndex(i_lab) )
      b = chan%Index2b( this%GetLabIndex(i_lab) )
      oa => labqn%o1
      ob => labqn%o2
      na = oa%n; la = oa%l; ja = oa%j
      nb = ob%n; lb = ob%l; jb = ob%j
      delab = 1.d0
      if(a == b) delab = 1.d0 / sqrt(2.d0)

      do i_relcm = 1, this%GetNumberRelCMStates()
        relcmqn => chan%GetRelCM( this%RelCMIndex(i_relcm) )
        nrel = relcmqn%nrel
        lrel = relcmqn%lrel
        ncm  = relcmqn%ncm
        lcm  = relcmqn%lcm
        s    = relcmqn%s
        jrel = relcmqn%Jrel
        ft = (1.d0 - (-1.d0)**(la + lb + s - lcm + chan%GetT())) / sqrt(2.d0)
        if(ft == 0.d0) cycle
        lammin = max(abs(lcm - lrel), abs(chan%GetJ() - s), abs(la - lb))
        lammax = min(lcm + lrel, chan%GetJ() + s, la + lb)
        tc = 0.d0
        do lam = lammin, lammax
          tc = tc + (-1.d0) **(lcm + lrel + s + chan%GetJ() ) * &
              & sqrt(dble(2*lam+1) * dble(2*s+1) * dble(ja+1) * dble(jb+1)) * &
              & ninejs%get(2*la,2*lb,2*lam,1,1,2*s,ja,jb,2*chan%GetJ()) * &
              & sqrt(dble(2*lam+1) * dble(2*jrel+1)) * &
              & sixjs%get(2*lcm,2*lrel,2*lam,2*s,2*chan%GetJ(),2*jrel) * &
              & tmbks%get(ncm,lcm,nrel,lrel,na,la,nb,lb,lam) * &
              & ft * delab
        end do
        this%mat(i_relcm, i_lab) = tc
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine InitTCoefs

  function GetMemoryTCoefs(this) result(mem)
    class(TCoefsIso), intent(in) :: this
    real(8) :: mem
    mem = dble( this%GetNumberLabStates() ) * dble( this%GetNumberRelCMStates() ) * 8.d0 / 1024.d0**3
    mem = mem + dble(size(this%RelCMQNs2Idx)) * 4.d0 / 1024.d0**3
    mem = mem + dble(size(this%LabIndex2Idx)) * 8.d0 / 1024.d0**3
  end function GetMemoryTCoefs

  subroutine FinTransRel2LabChannel(this)
    class(TransRel2LabChannelIso), intent(inout) :: this
    integer :: N
    do N = 0, this%Nmax
      call this%N2max(N)%fin()
    end do
    deallocate(this%N2max)
    if(this%GetNumberRelCMStates() * this%GetNumberLabStates() == 0) return
    deallocate(this%RelCM)
    this%labels2n => null()
    this%iphase => null()
    this%Index2a => null()
    this%Index2b => null()
    this%Nmax = -1
    this%J = -1
    this%P = -100
    this%T = -100
    this%NumberRelCMStates = 0
    this%NumberLabStates = 0
  end subroutine FinTransRel2LabChannel

  subroutine InitTransRel2LabChannel(this, chan, sps, Nmax, Jmax)
    class(TransRel2LabChannelIso), intent(inout) :: this
    type(TwoBodyLabIsoChan), intent(in), target :: chan
    type(OrbitsIsospin), intent(in), target :: sps
    integer, intent(in) :: Nmax, Jmax
    integer :: N, n_states = 0
    integer :: a, b
    type(SingleParticleOrbitIsospin), pointer :: oa, ob
    this%sps => sps
    this%J = chan%GetJ()
    this%P = chan%GetParity()
    this%T = chan%GetT()
    this%Nmax = Nmax

    allocate(this%N2max(0:Nmax))
    call set_rel_cm( this, Jmax, this%NumberRelCMStates, .false. )
    this%NumberLabStates = chan%GetNumberStates()
    if(this%GetNumberRelCMStates() * this%GetNumberLabStates() == 0) return
    this%labels2n => chan%labels2n
    this%iphase => chan%iphase
    this%Index2a => chan%n2label1
    this%Index2b => chan%n2label2
    call set_rel_cm( this, Jmax, n_states, .true. )
    allocate(this%Lab( this%GetNumberLabStates() ))
    do n_states = 1, this%GetNumberLabStates()
      a = this%Index2a(n_states)
      b = this%Index2b(n_states)
      oa => sps%GetOrbit(a)
      ob => sps%GetOrbit(b)
      call this%Lab(n_states)%set(oa,ob)
    end do

    do N = 0, this%GetNmax()
      call this%N2max(N)%init( this, N )
    end do

  end subroutine InitTransRel2LabChannel

  function GetMemoryChannel(this) result(mem)
    class(TransRel2LabChannelIso), intent(in) :: this
    real(8) :: mem
    integer :: N
    mem = 0.d0
    do N = 0, this%GetNmax()
      mem = mem + this%N2max(N)%GetMemory()
    end do
  end function GetMemoryChannel

  subroutine set_rel_cm( RelCMChan, Jmax, n_states, set_mode )
    use MyLibrary, only: triag
    class(TransRel2LabChannelIso), intent(inout), target :: RelCMChan
    integer, intent(in) :: Jmax
    integer, intent(inout) :: n_states
    logical, intent(in) :: set_mode
    integer :: nrel, lrel, Nmax, s, Ncm, Lcm, jrel
    type(RelativeCMQNs), pointer :: tmp

    if( set_mode ) then
      allocate( RelCMChan%RelCM( RelCMChan%GetNumberRelCMStates()) )
    end if
    n_states = 0
    do Nmax = 0, RelCMChan%GetNmax()
      do nrel = 0, Nmax / 2
        do Ncm = 0, Nmax/2 - nrel

          do jrel = 0, Jmax
            do s = 0, 1
              do lrel = abs(jrel - s), jrel + s
                if(2*nrel + lrel > RelCMChan%GetNmax() ) cycle
                Lcm = Nmax - 2*nrel - 2*Ncm - lrel
                if(triag(Lcm, jrel, RelCMChan%GetJ()) ) cycle
                if(mod(lrel+s+RelCMChan%GetT(),2) == 0) cycle
                if((-1) ** (lrel + Lcm) /= RelCMChan%GetParity() ) cycle
                n_states = n_states + 1
                if( set_mode ) then
                  ! becuase of a bug of intel fortran
                  !call RelCMChan%RelCM( n_states )%set(nrel,lrel,s,jrel,Ncm,Lcm)
                  tmp => RelCMChan%GetRelCM( n_states )
                  tmp%nrel = nrel
                  tmp%lrel = lrel
                  tmp%s    = s
                  tmp%jrel = jrel
                  tmp%Ncm  = Ncm
                  tmp%Lcm  = Lcm
                  tmp%Nmax = 2*nrel + lrel + 2*Ncm + Lcm
                end if
              end do
            end do
          end do

        end do
      end do
    end do
  end subroutine set_rel_cm

  subroutine InitTransRel2LabSpace(this, ms, e2max, Jmax)
    class(TransRel2LabSpaceIso), intent(inout) :: this
    type(TwoBodyLabIsoSpace), target, intent(in) :: ms
    type(OrbitsIsospin), pointer :: sps
    integer, intent(in) :: e2max, Jmax
    integer :: ch
    type(sys) :: s
    real(8) :: time

    sps => ms%sps
    this%sps => ms%sps
    this%NChan = ms%GetNumberChannels()
    this%Nmax = e2max
    this%Jmax = min(Jmax,e2max+1)
    call sixjs%init(0,2*e2max,.false., 0,2*e2max,.false., 0,2, .false.)
    call ninejs%init(0,2*ms%GetEmax(),.false., 0,2*ms%GetEmax(),.false., 1,1,.true., 1,1,.true.)
    call TMbks%init(e2max, 1.d0)
    allocate(this%jpt2idx( 0:e2max+1, -1:1, -1:1 ))
    allocate(this%jpt( this%GetNumberChannels() ))
    this%jpt2idx(:,:,:) = 0

    time = omp_get_wtime()
    do ch = 1, this%GetNumberChannels()
      this%jpt2idx( ms%jpt(ch)%GetJ(), ms%jpt(ch)%GetParity(), ms%jpt(ch)%GetT() ) = ch
      call this%jpt(ch)%init( ms%jpt(ch), sps, this%GetNmax(), this%GetJmax() )
    end do
    call timer%add(s%str('Transformation coefficient'), omp_get_wtime() - time)
    call TMbks%fin()
    call ninejs%fin()
    call sixjs%fin()

    this%is_Constructed=.True.
  end subroutine InitTransRel2LabSpace

  subroutine FinTransRel2LabSpace(this)
    class(TransRel2LabSpaceIso), intent(inout) :: this
    integer :: ch
    if(.not. this%is_Constructed) return
    do ch = 1, this%GetNumberChannels()
      call this%jpt(ch)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%NChan = 0
    this%Nmax = -1
    this%Jmax = -1
    this%is_Constructed=.False.
  end subroutine FinTransRel2LabSpace

  function GetMemory(this) result(mem)
    class(TransRel2LabSpaceIso), intent(in) :: this
    real(8) :: mem
    integer :: ch
    mem = 0.d0
    do ch = 1, this%GetNumberChannels()
      mem = mem + this%jpt(ch)%GetMemory()
    end do
  end function GetMemory

  subroutine SetTwoBodyLabQNs(this, o1, o2)
    class(TwoBodyLabQNs), intent(inout) :: this
    type(SingleParticleOrbitIsospin), intent(in), target :: o1, o2
    this%o1 => o1
    this%o2 => o2
    this%Nmax = o1%e + o2%e
  end subroutine SetTwoBodyLabQNs

  subroutine SetRelativeCMQNs(this, nrel, lrel, s, jrel, Ncm, Lcm)
    class(RelativeCMQNs), intent(inout) :: this
    integer, intent(in) :: nrel, lrel, s, jrel, Ncm, Lcm
    this%nrel = nrel
    this%lrel = lrel
    this%s    = s
    this%jrel = jrel
    this%Ncm  = Ncm
    this%Lcm  = Lcm
    this%Nmax = 2*nrel + lrel + 2*Ncm + Lcm
  end subroutine SetRelativeCMQNs

  function GetMatrixElement(this, ncm, lcm, nrel, lrel, s, jrel, J, T, a, b) result(res)
    class(TransRel2LabSpaceIso), intent(in) :: this
    integer, intent(in) :: ncm, lcm, nrel, lrel, s, jrel, J, T, a, b
    real(8) :: res
    integer :: Prty
    type(TransRel2LabChannelIso), pointer :: channel
    type(TCoefsIso), pointer :: tc
    integer :: idx, phase
    integer :: idx_relcm, idx_lab
    integer :: Nmax

    res = 0.d0
    Prty = (-1)**(lrel+lcm)
    if(mod(lrel+s+T,2)==0) return
    if(a==b .and. mod(J+T,2)==0) return
    channel => this%GetChannel(J,Prty,T)
    Nmax = 2*ncm+lcm+2*nrel+lrel
    tc => channel%N2max(Nmax)
    idx_relcm = tc%RelCMQNs2Idx(ncm, lcm, nrel, lrel, s, jrel)
    idx_lab = tc%LabIndex2Idx(a,b)
    if(idx_relcm * idx_lab == 0) return
    phase = tc%LabIndex2Phase(a,b)
    res = tc%mat(idx_relcm, idx_lab) * dble(phase)
  end function GetMatrixElement

end module TwoBodyTransCoefIso

