#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION


  public :: MACRO_ADD_SUFX(ASJacobiCMQNs)
  public :: MACRO_ADD_SUFX(NASJacobiCMQNs)
  public :: MACRO_ADD_SUFX(ThreeBodyLabQNs)
  public :: MACRO_ADD_SUFX(TCoefs)
  public :: MACRO_ADD_SUFX(TransJac2LabChannelIsospin)
  public :: MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)
  private
  ! methods
  private :: FinTCoefs
  private :: InitTCoefs
  private :: set_lab_jacobi_cm_indices
  private :: set_as_tcoefs
  private :: set_nas_tcoefs
  private :: GetNASJacCMIndexTCoefs
  private :: GetASJacCMIndexTCoefs
  private :: GetLabIndexTCoefs
  private :: GetNumberNASJacCMStatesTCoefs
  private :: GetNumberASJacCMStatesTCoefs
  private :: GetNumberLabStatesTCoefs
  private :: GetMemoryTCoefs

  private :: InitTransJac2LabChannel
  private :: FinTransJac2LabChannel
  private :: set_channel_indices
  private :: GetNmax
  private :: GetNumberBlock
  private :: GetNumberASJacCMStates
  private :: GetNumberNASJacCMStates
  private :: GetNumberLabStates
  private :: GetASJacCMFromIndex
  private :: GetASJacCMFromEIndex
  private :: GetNASJacCMFromIndex
  private :: GetNASJacCMFromEIndex
  private :: GetLabFromIndex
  private :: GetLabFromEIndex
  private :: GetMemoryChannel

  private :: InitTransJac2LabSpace
  private :: FinTransJac2LabSpace
  private :: GetIndex
  private :: GetNmaxSpace
  private :: GetJmaxSpace
  private :: GetNumberBlockSpace
  private :: GetChannelFromJPT
  private :: GetChannelFromIndex
  private :: GetNumberChannels
  private :: GetMemory

  type :: MACRO_ADD_SUFX(ASJacobiCMQNs)
    integer :: Nmax, E, i, jrel, Ncm, Lcm, idx
  contains
    procedure :: SetASJacobiCMQNs
    generic :: set => SetASJacobiCMQNs
  end type MACRO_ADD_SUFX(ASJacobiCMQNs)

  type :: MACRO_ADD_SUFX(NASJacobiCMQNs)
    integer :: Nmax, n12, l12, s12, j12, t12
    integer :: n3, l3, j3, jrel, Ncm, Lcm, alpha, idx
  contains
    procedure :: SetNASJacobiCMQNs
    generic :: set => SetNASJacobiCMQNs
  end type MACRO_ADD_SUFX(NASJacobiCMQNs)

  type :: MACRO_ADD_SUFX(ThreeBodyLabQNs)
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3
    integer :: j12, t12
    integer :: Nmax, idx
  contains
    procedure :: SetThreeBodyLabQNs
    generic :: set => SetThreeBodyLabQNs
  end type MACRO_ADD_SUFX(ThreeBodyLabQNs)

  type :: MACRO_ADD_SUFX(TCoefs)
#if defined(half_precision)
    type(my_real16), allocatable :: mat(:,:)
#elif defined(single_precision)
    real(4), allocatable :: mat(:,:)
#elif defined(double_precision)
    real(8), allocatable :: mat(:,:)
#endif
    integer, allocatable :: NASJacCMIndex(:)
    integer, allocatable :: ASJacCMIndex(:)
    integer, allocatable :: LabIndex(:)
    integer :: NumberNASJacCMStates = 0
    integer :: NumberASJacCMStates = 0
    integer :: NumberLabStates = 0
  contains
    procedure :: InitTCoefs
    procedure :: FinTCoefs
    procedure :: GetNASJacCMIndexTCoefs
    procedure :: GetASJacCMIndexTCoefs
    procedure :: GetLabIndexTCoefs
    procedure :: GetNumberNASJacCMStatesTCoefs
    procedure :: GetNumberASJacCMStatesTCoefs
    procedure :: GetNumberLabStatesTCoefs
    procedure :: GetMemoryTCoefs
    generic :: init => InitTCoefs
    generic :: fin => FinTCoefs
    generic :: GetNASJacCMIndex => GetNASJacCMIndexTCoefs
    generic :: GetASJacCMIndex => GetASJacCMIndexTCoefs
    generic :: GetLabIndex => GetLabIndexTCoefs
    generic :: GetNumberNASJacCMStates => GetNumberNASJacCMStatesTCoefs
    generic :: GetNumberASJacCMStates => GetNumberASJacCMStatesTCoefs
    generic :: GetNumberLabStates => GetNumberLabStatesTCoefs
    generic :: GetMemory => GetMemoryTCoefs
  end type MACRO_ADD_SUFX(TCoefs)

  type :: MACRO_ADD_SUFX(TransJac2LabChannelIsospin)
    type(MACRO_ADD_SUFX(TCoefs)), allocatable :: N3max(:)
    type(ThreeBodyJacIsoSpace), pointer :: SpJac
    type(ThreeBodyLabIsoChan), pointer :: ChLab
    type(MACRO_ADD_SUFX(ASJacobiCMQNs)), allocatable :: ASJacCM(:)
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), allocatable :: NASJacCM(:)
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), allocatable :: Lab(:)
    integer, private :: Nmax = -1
    integer, private :: NBlock = 0
    integer, private :: NumberASJacCMStates = 0
    integer, private :: NumberNASJacCMStates = 0
    integer, private :: NumberLabStates = 0
  contains
    procedure :: InitTransJac2LabChannel
    procedure :: FinTransJac2LabChannel
    procedure :: GetNmax
    procedure :: GetNumberBlock
    procedure :: GetNumberASJacCMStates
    procedure :: GetNumberNASJacCMStates
    procedure :: GetNumberLabStates
    procedure :: GetASJacCMFromIndex
    procedure :: GetASJacCMFromEIndex
    procedure :: GetNASJacCMFromIndex
    procedure :: GetNASJacCMFromEIndex
    procedure :: GetLabFromIndex
    procedure :: GetLabFromEIndex
    procedure :: GetMemoryChannel
    generic :: init => InitTransJac2LabChannel
    generic :: fin => FinTransJac2LabChannel
    generic :: GetASJacCM => GetASJacCMFromIndex, GetASJacCMFromEIndex
    generic :: GetNASJacCM => GetNASJacCMFromIndex, GetNASJacCMFromEIndex
    generic :: GetLab => GetLabFromIndex, GetLabFromEIndex
    generic :: GetMemory => GetMemoryChannel
  end type MACRO_ADD_SUFX(TransJac2LabChannelIsospin)

  type :: MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), allocatable :: jpt(:)
    type(ThreeBodyLabIsoSpace), pointer :: SpLab
    integer, private :: NChan = 0
    integer, private :: NBlock = 0
    integer, private :: Nmax = -1
    integer, private :: Jmax = -1
    logical :: is_Constructed=.false.
  contains
    procedure :: InitTransJac2LabSpace
    procedure :: FinTransJac2LabSpace
    procedure :: GetIndex
    procedure :: GetNmaxSpace
    procedure :: GetJmaxSpace
    procedure :: GetNumberBlockSpace
    procedure :: GetChannelFromJPT
    procedure :: GetChannelFromIndex
    procedure :: GetNumberChannels
    procedure :: GetMemory
    generic :: init => InitTransJac2LabSpace
    generic :: fin => FinTransJac2LabSpace
    generic :: GetNmax => GetNmaxSpace
    generic :: GetJmax => GetJmaxSpace
    generic :: GetNumberBlock => GetNumberBlockSpace
    generic :: GetChannel => GetChannelFromJPT, GetChannelFromIndex
  end type MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)

#ifdef old_tcoefs
  type(SixJsStore), private :: Ls, LL
  type(NineJsStore), private :: ls2j_3
#else
  type(SixJsStore), private :: Lls, e3o3, e2o4
#endif
  type(NineJsStore), private :: ls2j_2
  type(TMbracketStore), private :: tmbk_1, tmbk_2
contains
  function GetMemory(this) result(memory)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    real(8) :: memory
    integer :: ch
    memory = 0.d0
    do ch = 1, this%GetNumberChannels()
      memory = memory + this%jpt(ch)%GetMemory()
    end do
  end function GetMemory

  function GetMemoryChannel(this) result(memory)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    real(8) :: memory
    integer :: N
    memory = 0.d0
    do N = 0, this%GetNmax()
      memory = memory + this%N3max(N)%GetMemory()
    end do
  end function GetMemoryChannel

  function GetMemoryTCoefs(this) result(memory)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    real(8) :: memory
#if defined(half_precision)
    memory = 2.d0 * dble(size(this%mat,1)) * dble(size(this%mat,2))
#elif defined(single_precision)
    memory = 4.d0 * dble(size(this%mat,1)) * dble(size(this%mat,2))
#elif defined(double_precision)
    memory = 8.d0 * dble(size(this%mat,1)) * dble(size(this%mat,2))
#endif
    memory = memory / 1024.d0**3
  end function GetMemoryTCoefs

  function GetNumberChannels(this) result(NChan)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetIndex(this, j, p, t) result(idx)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%SpLab%GetIndex(j,p,t)
  end function GetIndex

  function GetNmaxSpace(this) result(Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmaxSpace

  function GetJmaxSpace(this) result(Jmax)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmaxSpace

  function GetNumberBlockSpace(this) result(nblock)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: this
    integer :: nblock
    nblock = this%nblock
  end function GetNumberBlockSpace

  function GetChannelFromIndex(this, idx) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in), target :: this
    integer, intent(in) :: idx
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%jpt(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this, j, p, t) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in), target :: this
    integer, intent(in) :: j, p, t
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), pointer :: r
    r => this%GetChannel( this%GetIndex(j, p, t) )
  end function GetChannelFromJPT

  function GetNmax(this) result(Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetNumberBlock(this) result(NBlock)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    integer :: NBlock
    NBlock = this%NBlock
  end function GetNumberBlock

  function GetNumberASJacCMStates(this) result(n_states)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberASJacCMStates
  end function GetNumberASJacCMStates

  function GetNumberNASJacCMStates(this) result(n_states)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberNASJacCMStates
  end function GetNumberNASJacCMStates

  function GetNumberLabStates(this) result(n_states)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberLabStates
  end function GetNumberLabStates

  function GetASJacCMFromIndex(this, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: i
    type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: r
    r => null()
    if(i == 0) return
    r => this%ASJacCM(i)
  end function GetASJacCMFromIndex

  function GetASJacCMFromEIndex(this, E, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: E, i
    type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: r
    r => this%GetASJacCM( this%N3max(E)%GetASJacCMIndex(i) )
  end function GetASJacCMFromEIndex

  function GetNASJacCMFromIndex(this, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: i
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: r
    r => null()
    if(i == 0) return
    r => this%NASJacCM(i)
  end function GetNASJacCMFromIndex

  function GetNASJacCMFromEIndex(this, E, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: E, i
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: r
    r => this%GetNASJacCM( this%N3max(E)%GetNASJacCMIndex(i) )
  end function GetNASJacCMFromEIndex

  function GetLabFromIndex(this, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: i
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: r
    r => null()
    if(i == 0) return
    r => this%Lab(i)
  end function GetLabFromIndex

  function GetLabFromEIndex(this, E, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in), target :: this
    integer, intent(in) :: E, i
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: r
    r => this%GetLab( this%N3max(E)%GetLabIndex(i) )
  end function GetLabFromEIndex
  !
  ! TCoefs
  !
  function GetNASJacCMIndexTCoefs(this, i) result(idx)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%NASJacCMIndex(i)
  end function GetNASJacCMIndexTCoefs

  function GetASJacCMIndexTCoefs(this, i) result(idx)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%ASJacCMIndex(i)
  end function GetASJacCMIndexTCoefs

  function GetLabIndexTCoefs(this, i) result(idx)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%LabIndex(i)
  end function GetLabIndexTCoefs

  function GetNumberNASJacCMStatesTCoefs(this) result(n_states)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberNASJacCMStates
  end function GetNumberNASJacCMStatesTCoefs

  function GetNumberASJacCMStatesTCoefs(this) result(n_states)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberASJacCMStates
  end function GetNumberASJacCMStatesTCoefs

  function GetNumberLabStatesTCoefs(this) result(n_states)
    class(MACRO_ADD_SUFX(TCoefs)), intent(in) :: this
    integer :: n_states
    n_states = this%NumberLabStates
  end function GetNumberLabStatesTCoefs

  !
  ! TransJac2LabChannelIsospin
  !

  subroutine FinTCoefs(this)
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    if(this%GetNumberLabStates() * this%GetNumberASJacCMStates() * this%GetNumberNASJacCMStates() == 0) return
    deallocate(this%mat)
    deallocate(this%NASJacCMIndex)
    deallocate(this%ASJacCMIndex)
    deallocate(this%LabIndex)
    this%NumberASJacCMStates = 0
    this%NumberNASJacCMStates = 0
    this%NumberLabStates = 0
  end subroutine FinTCoefs

  subroutine InitTCoefs(this, chan, N)
    use MyLibrary, only: triag
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: chan
    integer, intent(in) :: N
    type(MACRO_ADD_SUFX(TCoefs)) :: temp
    integer :: n_nas, n_as, n_lab
    type(sys) :: s
    real(8) :: time

    call set_lab_jacobi_cm_indices(this, chan, N, n_nas, n_as, n_lab, .false.)
    if(n_nas * n_as * n_lab == 0) return

    call set_lab_jacobi_cm_indices(this, chan, N, n_nas, n_as, n_lab, .true.)
    call set_lab_jacobi_cm_indices(temp, chan, N, n_nas, n_as, n_lab, .true.)
    allocate(this%mat( this%GetNumberASJacCMStates(), this%GetNumberLabStates() ))
    allocate(temp%mat( this%GetNumberNASJacCMStates(),this%GetNumberLabStates() ))
#if defined(half_precision) || defined(single_precision)
    this%mat(:,:) = 0.0
    temp%mat(:,:) = 0.0
#elif defined(double_precision)
    this%mat(:,:) = 0.d0
    temp%mat(:,:) = 0.d0
#endif
    time = omp_get_wtime()
    call set_nas_tcoefs(temp, chan)
    call timer%Add(s%str("T coefficient (NAS|Lab)"), omp_get_wtime() - time)
    time = omp_get_wtime()
    call set_as_tcoefs(this, temp, chan)
    call timer%Add(s%str("T coefficient ( AS|Lab)"), omp_get_wtime() - time)
    call temp%fin()

  end subroutine InitTCoefs

#ifdef old_tcoefs
  subroutine set_nas_tcoefs(this, chan)
    use MyLibrary, only: hat
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: chan
    integer :: i_jaccm, i_lab
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: nas
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc
    integer :: na,la,ja,nb,lb,jb,nc,lc,jc,jab,tab,J
    integer :: n12,l12,s12,j12,t12,n3,l3,j3,jrel,ncm,lcm
    real(8) :: tc
    ! -- for tcoef --
    integer :: lab_min, lab_max, lab, lit_min, lit_max, lit, nit
    integer :: lam_min, lam_max, lam, ltot_min, ltot_max, ltot
    integer :: stot_min, stot_max, stot, lreltot_min, lreltot_max, lreltot
    real(8) :: snj_ab, rel_gmosh, cm_gmosh, sjs_lc, snj_c

    J = chan%ChLab%GetJ()
    !$omp parallel
    !$omp do private( i_lab, labqns, oa, ob, oc, jab, tab, na, la, ja, nb, lb, jb, nc, lc, jc, &
    !$omp &  i_jaccm, nas, n12, l12, s12, j12, t12, n3, l3, j3, jrel, ncm, lcm, tc, &
    !$omp &  lab_min, lab_max, lab, snj_ab, lit_min, lit_max, lit, nit, rel_gmosh, &
    !$omp &  lam_min, lam_max, lam, cm_gmosh, ltot_min, ltot_max, ltot, sjs_lc, &
    !$omp &  stot_min, stot_max, stot, snj_c, lreltot_min, lreltot_max, lreltot) schedule(dynamic)
    do i_lab = 1, this%GetNumberLabStates()
      labqns => chan%GetLab( this%GetLabIndex( i_lab ))
      oa => labqns%o1
      ob => labqns%o2
      oc => labqns%o3
      jab = labqns%j12
      tab = labqns%t12
      na = oa%n; la = oa%l; ja = oa%j
      nb = ob%n; lb = ob%l; jb = ob%j
      nc = oc%n; lc = oc%l; jc = oc%j

      do i_jaccm = 1, this%GetNumberNASJacCMStates()
        nas => chan%GetNASJacCM( this%GetNASJacCMIndex( i_jaccm ))
        t12 = nas%t12
        if(t12 /= tab) cycle
        n12 = nas%n12
        l12 = nas%l12
        s12 = nas%s12
        j12 = nas%j12
        n3  = nas%n3
        l3  = nas%l3
        j3  = nas%j3
        jrel= nas%jrel
        ncm = nas%ncm
        lcm = nas%lcm

        tc = 0.d0
        lab_min = max(iabs(jab-s12), iabs(la-lb))
        lab_max = min(     jab+s12,       la+lb)
        do lab = lab_min, lab_max
          snj_ab = ls2j_2%get(2*la, 1, ja, 2*lb, 1, jb, 2*lab, 2*s12, 2*jab)
          if(abs(snj_ab) < 1.d-16) cycle

          lit_min = max(iabs(lab-l12), iabs(jab-j12))
          lit_max = min(     lab+l12,       jab+j12, 2*na + la + 2*nb + lb - 2*n12 - l12)
          do lit = lit_min, lit_max
            nit = (2*na + la + 2*nb + lb - 2*n12 - l12 - lit) / 2

            rel_gmosh = tmbk_1%get(nit, lit, n12, l12, na, la, nb, lb, lab)
            if(abs(rel_gmosh) < 1.d-16) cycle

            lam_min = max(iabs(lit-lc), iabs(lcm-l3))
            lam_max = min(     lit+lc,       lcm+l3)
            do lam = lam_min, lam_max
              cm_gmosh = tmbk_2%get(ncm, lcm, n3, l3, nit, lit, nc, lc, lam)
              if(abs(cm_gmosh) < 1.d-16) cycle

              ltot_min = max(iabs(lam-l12), iabs(lab-lc))
              ltot_max = min(     lam+l12,       lab+lc)
              do ltot = ltot_min, ltot_max
                sjs_lc = LL%get(2*l12, 2*lit, 2*lab, 2*lc, 2*ltot, 2*lam)
                if(abs(sjs_lc) < 1.d-16) cycle

                stot_min = max(iabs(J-2*ltot), iabs(2*s12-1))
                stot_max = min(     J+2*ltot,       2*s12+1)
                do stot = stot_min, stot_max, 2
                  snj_c = ls2j_3%get(2*lab,2*s12,2*jab,2*lc,1,jc,2*ltot,stot,J )
                  if(abs(snj_c) < 1.d-16) cycle

                  lreltot_min = max(iabs(l12-l3), iabs(stot-jrel)/2, iabs(ltot-lcm))
                  lreltot_max = min(     l12+l3,      (stot+jrel)/2,      ltot+lcm)
                  do lreltot = lreltot_min, lreltot_max
                    tc = tc + (-1.d0) ** (lc + lam + lab + ltot + l12 + (stot + J)/2) * &
                        & hat(ja) * hat(jb) * hat(jc) * hat(2*jab) * hat(jrel) * &
                        & hat(2*j12) * hat(j3) * hat(2*s12) * dble(stot+1) * &
                        & dble(2*lab+1) * dble(2*ltot+1) * dble(2*lreltot+1) * dble(2*lam+1) * &
                        & rel_gmosh * cm_gmosh * snj_ab * snj_c * sjs_lc * &
                        & ls2j_3%get(2*l12,2*s12,2*j12,2*l3,1,j3,2*lreltot,stot,jrel) * &
                        & LL%get(2*lcm, 2*l3, 2*lam, 2*l12, 2*ltot, 2*lreltot) * &
                        & Ls%get(2*lcm, 2*lreltot, 2*ltot, stot, J, jrel)

                  end do
                end do
              end do
            end do
          end do
        end do
#if defined(half_precision) || defined(single_precision)
        this%mat(i_jaccm, i_lab) = real(tc)
#elif defined(double_precision)
        this%mat(i_jaccm, i_lab) = tc
#endif
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine set_nas_tcoefs

#else
  subroutine set_nas_tcoefs(this, chan)
    use MyLibrary, only: hat
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: chan
    integer :: i_jaccm, i_lab
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: nas
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc
    integer :: na,la,ja,nb,lb,jb,nc,lc,jc,jab,tab,J
    integer :: n12,l12,s12,j12,t12,n3,l3,j3,jrel,ncm,lcm
    real(8) :: tc
    ! -- for tcoef --
    integer :: lab_min, lab_max, lab, lit_min, lit_max, lit, nit
    integer :: LLam_min, LLam_max, LLam, lam_min, lam_max, lam
    real(8) :: flab, flit, fLLam, flam

    J = chan%ChLab%GetJ()
    !$omp parallel
    !$omp do private( i_lab, labqns, oa, ob, oc, jab, tab, na, la, ja, nb, lb, jb, nc, lc, jc, &
    !$omp &  i_jaccm, nas, n12, l12, s12, j12, t12, n3, l3, j3, jrel, ncm, lcm, tc, &
    !$omp &  lab_min, lab_max, lab, flab, lit_min, lit_max, lit, nit, flit, &
    !$omp &  LLam_min, LLam_max, LLam, fLLam, lam_min, lam_max, lam, flam ) schedule(dynamic)
    do i_lab = 1, this%GetNumberLabStates()
      labqns => chan%GetLab( this%GetLabIndex( i_lab ))
      oa => labqns%o1
      ob => labqns%o2
      oc => labqns%o3
      jab = labqns%j12
      tab = labqns%t12
      na = oa%n; la = oa%l; ja = oa%j
      nb = ob%n; lb = ob%l; jb = ob%j
      nc = oc%n; lc = oc%l; jc = oc%j

      do i_jaccm = 1, this%GetNumberNASJacCMStates()
        nas => chan%GetNASJacCM( this%GetNASJacCMIndex( i_jaccm ))
        t12 = nas%t12
        if(t12 /= tab) cycle
        n12 = nas%n12
        l12 = nas%l12
        s12 = nas%s12
        j12 = nas%j12
        n3  = nas%n3
        l3  = nas%l3
        j3  = nas%j3
        jrel= nas%jrel
        ncm = nas%ncm
        lcm = nas%lcm

        tc = 0.d0
        lab_min = max(iabs(jab-s12), iabs(la-lb))
        lab_max = min(     jab+s12,       la+lb)
        do lab = lab_min, lab_max
          flab = dble(2*lab+1) * ls2j_2%get(2*la, 1, ja, 2*lb, 1, jb, 2*lab, 2*s12, 2*jab)
          if(abs(flab) < 1.d-16) cycle

          lit_min = max(iabs(lab-l12), iabs(jab-j12))
          lit_max = min(     lab+l12,       jab+j12, 2*na+la+2*nb+lb-2*n12-l12, 2*Ncm+lcm+2*n3+l3-2*nc-lc)
          do lit = lit_min, lit_max
            nit = (2*na + la + 2*nb + lb - 2*n12 - l12 - lit) / 2
            flit = tmbk_1%get(nit, lit, n12, l12, na, la, nb, lb, lab) * &
                & Lls%get(2*lit, 2*l12, 2*lab, 2*s12, 2*jab, 2*j12)
            if(abs(flit) < 1.d-16) cycle

            LLam_min = max(iabs(2*j12-J), iabs(2*lit-jc), iabs(2*Lcm-j3))
            LLam_max = min(     2*j12+J,       2*lit+jc,       2*Lcm+j3)
            do LLam = LLam_min, LLam_max, 2
              fLLam = (-1.d0)**( (LLam-Jrel)/2 ) * &
                  & dble(LLam+1) * &
                  & e3o3%get(2*j12, 2*lit, 2*jab, jc, J, LLam) * &
                  & e2o4%get(2*Lcm, j3, LLam, 2*j12, J, jrel)
              if(abs(FLLam) < 1.d-16) cycle

              lam_min = max(iabs(lit-lc), iabs(1-LLam)/2, iabs(Lcm-l3))
              lam_max = min(     lit+lc,      (1+LLam)/2,     (Lcm+l3))
              do lam = lam_min, lam_max
                flam = dble(2*lam+1) * &
                    & e3o3%get(2*lit, 2*lc, 2*lam, 1, LLam, jc) * &
                    & e3o3%get(2*Lcm, 2*l3, 2*lam, 1, LLam, j3) * &
                    & tmbk_2%get(ncm, lcm, n3, l3, nit, lit, nc, lc, lam)
                if(abs(flam) < 1.d-16) cycle
                tc = tc + flam * fLLam * flab * flit

              end do
            end do
          end do
        end do
        tc = tc * hat(ja) * hat(jb) * hat(jc) * hat(2*jab) * hat(2*s12) * hat(2*j12) * hat(j3) * hat(jrel) * &
            & (-1.d0)**( s12+l12+1+j12+lc+l3+ (jc+J)/2 )

#if defined(half_precision) || defined(single_precision)
        this%mat(i_jaccm, i_lab) = real(tc)
#elif defined(double_precision)
        this%mat(i_jaccm, i_lab) = tc
#endif
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine set_nas_tcoefs

#endif


  subroutine set_as_tcoefs(this, temp, chan)
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TCoefs)), intent(in) :: temp
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: chan
    type(ThreeBodyJacIsoChan), pointer :: ChJac
    integer :: i_nas, i_as
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: nas
    type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: as
    integer :: m1, m2, m3
#if defined(half_precision) || defined(single_precision)
    real(4), allocatable :: mat(:,:)
#elif defined(double_precision)
    real(8), allocatable :: mat(:,:)
#endif
#if defined(half_precision)
    real(4), allocatable :: work1(:,:), work2(:,:)
#endif

    m1 = this%GetNumberASJacCMStates()
    m2 = this%GetNumberNASJacCMStates()
    m3 = this%GetNumberLabStates()
    allocate(mat(m1,m2))
#if defined(half_precision) || defined(single_precision)
    mat = 0.0
#elif defined(double_precision)
    mat = 0.d0
#endif

    !$omp parallel
    !$omp do private(i_as, as, i_nas, nas, ChJac)
    do i_as = 1, m1
      as => chan%GetASJacCM( this%GetASJacCMIndex(i_as) )
      do i_nas = 1, m2
        nas => chan%GetNASJacCM( this%GetNASJacCMIndex(i_nas) )
        if(as%Ncm /= nas%Ncm) cycle
        if(as%Lcm /= nas%Lcm) cycle
        if(as%Nmax /= nas%Nmax) cycle
        if(as%jrel /= nas%jrel) cycle
        ChJac => chan%SpJac%GetChannel(nas%jrel, (-1)**(nas%l12+nas%l3), chan%ChLab%GetT())
        if(.not. associated(ChJac) ) cycle
#if defined(half_precision) || defined(single_precision)
        mat(i_as, i_nas) = real( ChJac%cfp(as%E)%m( as%i, nas%alpha ) )
#elif defined(double_precision)
        mat(i_as, i_nas) = ChJac%cfp(as%E)%m( as%i, nas%alpha )
#endif
      end do
    end do
    !$omp end do
    !$omp end parallel
#if defined(half_precision)
    allocate(work1( m2, m3 ), work2( m1, m3) )
    work1 = temp%mat
    call sgemm('n','n',m1,m3,m2,1.0,mat,m1,work1,m2,0.0,work2,m1)
    this%mat = work2
    deallocate(work1, work2)
#elif defined(single_precision)
    call sgemm('n','n',m1,m3,m2,1.0,mat,m1,temp%mat,m2,0.0,this%mat,m1)
#elif defined(double_precision)
    call dgemm('n','n',m1,m3,m2,1.d0,mat,m1,temp%mat,m2,0.d0,this%mat,m1)
#endif
    deallocate(mat)
  end subroutine set_as_tcoefs

  subroutine set_lab_jacobi_cm_indices(this, chan, Nmax, n_nas, n_as, n_lab, set_mode)
    class(MACRO_ADD_SUFX(TCoefs)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(in) :: chan
    integer, intent(in) :: Nmax
    integer, intent(inout) :: n_nas, n_as, n_lab
    logical, intent(in) :: set_mode
    type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: asqns
    type(MACRO_ADD_SUFX(NASJacobiCMQNs)), pointer :: nasqns
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
    integer :: idx, cnt

    if(set_mode) then
      this%NumberNASJacCMStates = n_nas
      this%NumberASJacCMStates = n_as
      this%NumberLabStates = n_lab
      allocate(this%NASJacCMIndex( this%GetNumberNASJacCMStates() ))
      allocate(this%ASJacCMIndex( this%GetNumberASJacCMStates() ))
      allocate(this%LabIndex( this%GetNumberLabStates() ))
      this%NASJacCMIndex(:) = 0
      this%ASJacCMIndex(:) = 0
      this%LabIndex(:) = 0
    end if

    n_nas = 0
    n_as = 0
    n_lab = 0

    cnt = 0
    do idx = 1, chan%GetNumberLabStates()
      labqns => chan%GetLab(idx)
      if(labqns%Nmax /= Nmax ) cycle
      cnt = cnt + 1
      if(set_mode) then
        labqns%idx = cnt
        this%LabIndex(cnt) = idx
      end if
    end do
    n_lab = cnt

    cnt = 0
    do idx = 1, chan%GetNumberASJacCMStates()
      asqns => chan%GetASJacCM(idx)
      if(asqns%Nmax /= Nmax ) cycle
      cnt = cnt + 1
      if(set_mode) then
        asqns%idx = cnt
        this%ASJacCMIndex(cnt) = idx
      end if
    end do
    n_as = cnt

    cnt = 0
    do idx = 1, chan%GetNumberNASJacCMStates()
      nasqns => chan%GetNASJacCM(idx)
      if(nasqns%Nmax /= Nmax ) cycle
      cnt = cnt + 1
      if(set_mode) then
        nasqns%idx = cnt
        this%NASJacCMIndex(cnt) = idx
      end if
    end do
    n_nas = cnt
  end subroutine set_lab_jacobi_cm_indices

  subroutine FinTransJac2LabChannel(this)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(inout) :: this
    integer :: N

    do N = 0, this%Nmax
      call this%N3max(N)%fin()
    end do
    this%Nmax = -1
    deallocate(this%N3max)
    if( this%GetNumberASJacCMStates() == 0) return
    if( this%GetNumberNASJacCMStates() == 0) return
    if( this%GetNumberLabStates() == 0) return

    this%SpJac => null()
    this%ChLab => null()
    deallocate(this%ASJacCM)
    deallocate(this%NASJacCM)
    deallocate(this%Lab)
    this%NumberASJacCMStates = 0
    this%NumberNASJacCMStates = 0
    this%NumberLabStates = 0
  end subroutine FinTransJac2LabChannel

  subroutine InitTransJac2LabChannel(this, ChLab, SpJac, Nmax, NBlock)
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), intent(in), target :: SpJac
    type(ThreeBodyLabIsoChan), intent(in), target :: ChLab
    integer, intent(in) :: Nmax
    integer, intent(in) :: NBlock
    integer :: n_as, n_nas, n_lab, N
    this%ChLab => ChLab
    this%SpJac => SpJac
    this%Nmax = Nmax
    this%NBlock = NBlock
    call set_channel_indices(this, n_as, n_nas, n_lab, .false.)
    call set_channel_indices(this, n_as, n_nas, n_lab, .true.)
    allocate(this%N3max(0:this%GetNmax() ) )
    do N = 0, this%GetNmax()
      call this%N3max(N)%init( this, N )
    end do

  end subroutine InitTransJac2LabChannel

  subroutine set_channel_indices(this, n_as, n_nas, n_lab, set_mode)
    use MyLibrary, only: triag
    class(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), intent(inout) :: this
    integer, intent(inout) :: n_as, n_nas, n_lab
    logical, intent(in) :: set_mode
    type(OrbitsIsospin), pointer :: sps
    integer :: ABCJT(5)
    integer :: i1, i2, i3, j12, t12
    integer :: idx
    integer :: N123, n12, n3, ncm, l12, l3, lcm, s12, j3, jrel, E, i, prel, alpha
    type(ThreeBodyJacIsoChan), pointer :: ChJac

    if( set_mode ) then
      this%NumberASJacCMStates = n_as
      this%NumberNASJacCMStates = n_nas
      this%NumberLabStates = n_lab
      allocate( this%ASJacCM( this%GetNumberASJacCMStates() ))
      allocate( this%NASJacCM( this%GetNumberNASJacCMStates() ))
      allocate( this%Lab( this%GetNumberLabStates() ))
    end if

    n_as = 0
    n_nas = 0
    n_lab = 0

    n_lab = this%ChLab%GetNumberStates()
    if(set_mode) then
      sps => this%ChLab%sps
      do idx = 1, n_lab
        ABCJT = this%ChLab%GetABCJT(idx)
        i1 = ABCJT(1)
        i2 = ABCJT(2)
        i3 = ABCJT(3)
        j12= ABCJT(4)
        t12= ABCJT(5)
        call this%Lab(idx)%set(sps%GetOrbit(i1),sps%GetOrbit(i2),sps%GetOrbit(i3),j12,t12)
      end do
    end if

    do N123 = 0, this%GetNmax()
      do E = 0, N123
        do ncm = 0, (n123 - E)/2
          lcm = n123 - E - 2 * ncm
          if(E + 2*ncm + lcm > this%GetNmax() ) cycle

          do jrel = iabs(this%ChLab%GetJ() - 2*lcm), min(this%ChLab%GetJ() + 2*lcm, 2*this%GetNmax()+3, this%SpJac%GetJmax()), 2
            do prel = 1, -1, -2
              if(this%ChLab%GetParity() /= prel * (-1) ** lcm) cycle
              ChJac => this%SpJac%GetChannel(jrel, prel, this%ChLab%GetT())
              if(.not. associated(ChJac) ) cycle
              if( E > ChJac%GetNmax() ) cycle
              do i = 1, ChJac%cfp(E)%GetNumberAStates()
                n_as = n_as + 1
                if( .not. set_mode ) cycle
                call this%ASJacCM(n_as)%set(E,i,jrel,Ncm,Lcm)
              end do
            end do
          end do
        end do
      end do
    end do

    do N123 = 0, this%GetNmax()
      do n12 = 0, n123 / 2
        do n3 = 0, n123 / 2 - n12
          do ncm = 0, n123 / 2 - n12 - n3
            do l12 = 0, n123 - 2 * n12 - 2 * n3 - 2 * ncm
              do l3 = 0, n123 - 2 * n12 - 2 * n3 - 2 * ncm - l12
                lcm = n123 - 2 * n12 - 2 * n3 - 2 * ncm - l12 - l3
                if(2 * n12 + l12 + 2 * n3 + l3 > this%GetNmax()) cycle

                do s12 = 0, 1
                  do j12 = iabs(l12 - s12), l12 + s12
                    do j3 = iabs(2 * l3 - 1), 2 * l3 + 1, 2
                      do t12 = 0, 1
                        if(triag(2*t12, 1, this%ChLab%GetT())) cycle

                        do jrel = iabs(2 * j12 - j3), min(2 * j12 + j3, this%SpJac%GetJmax() ), 2

                          if(triag(2 * lcm, jrel, this%ChLab%GetJ() )) cycle
                          if((-1) ** (l12 + s12 + t12) /= -1) cycle
                          if((-1) ** (l12 + l3 + lcm) /= this%ChLab%GetParity() ) cycle
                          ChJac => this%SpJac%GetChannel(jrel, (-1)**(l12+l3), this%ChLab%GetT())
                          if(.not. associated(ChJac) ) cycle
                          if( 2*n12+l12+2*n3+l3 > ChJac%GetNmax() ) cycle
                          if( ChJac%GetNASIndex( n12,l12,s12,j12,t12,n3,l3,j3 ) == 0) cycle
                          alpha = ChJac%GetAlphaNAS( ChJac%GetNASIndex( n12, l12, s12, j12, t12, n3, l3, j3 ))
                          if(alpha == 0) cycle
                          n_nas = n_nas + 1
                          if( .not. set_mode ) cycle
                          call this%NASJacCM(n_nas)%set(n12,l12,s12,j12,t12,n3,l3,j3,jrel,Ncm,Lcm,alpha)
                        end do

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine set_channel_indices

  subroutine FinTransJac2LabSpace(this)
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(inout) :: this
    integer :: ch
    if(.not. this%is_Constructed) return
    do ch = 1, this%GetNumberChannels()
      call this%jpt(ch)%fin()
    end do
    deallocate(this%jpt)
    this%SpLab => null()
    this%NChan = 0
    this%NBlock = 0
    this%Nmax = -1
    this%is_Constructed=.false.
  end subroutine FinTransJac2LabSpace

  subroutine InitTransJac2LabSpace(this, SpJac, SpLab, e3max, Jmax, Nblock)
    use ThreeBodyJacobiSpace
    class(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), intent(in), target :: SpJac
    type(ThreeBodyLabIsoSpace), intent(in), target :: SpLab
    integer, intent(in) :: e3max, Jmax, Nblock
    integer :: ich
    if(myrank == 0) then
      write(*,'("# Calculating Transformation coefficient ...")')
    end if
    this%Jmax = Jmax
    this%Nmax = e3max
    call store_couplings_for_tcoefs(this%GetNmax())
    this%SpLab => SpLab
    this%NChan = SpLab%GetNumberChannels()
    this%Nblock = Nblock

    allocate(this%jpt(this%NChan))
    do ich = 1, this%GetNumberChannels()
      if(myrank==0) write(*,'(a,i4,a,i4)') "calculating T-Coefs channel: ", ich, "/", this%NChan
      call this%jpt(ich)%init( SpLab%GetChannel(ich), SpJac, this%GetNmax(), this%NBlock )
    end do

    call release_couplings_for_tcoefs()
    this%is_Constructed=.True.
  end subroutine InitTransJac2LabSpace

  subroutine SetASJacobiCMQNs(this, E, i, jrel, Ncm, Lcm)
    class(MACRO_ADD_SUFX(ASJacobiCMQNs)), intent(inout) :: this
    integer, intent(in) :: E, i, jrel, Ncm, Lcm
    this%E = E
    this%i = i
    this%jrel = jrel
    this%Ncm = Ncm
    this%Lcm = Lcm
    this%Nmax = this%E + 2*this%Ncm + this%Lcm
  end subroutine SetASJacobiCMQNs

  subroutine SetNASJacobiCMQNs(this, n12, l12, s12, j12, t12, n3, l3, j3, Jrel, Ncm, Lcm, alpha)
    class(MACRO_ADD_SUFX(NASJacobiCMQNs)), intent(inout) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3, Jrel, Ncm, Lcm, alpha
    this%n12 = n12
    this%l12 = l12
    this%s12 = s12
    this%j12 = j12
    this%t12 = t12
    this%n3  = n3
    this%l3  = l3
    this%j3  = j3
    this%jrel= jrel
    this%Ncm = Ncm
    this%Lcm = Lcm
    this%alpha = alpha
    this%Nmax = 2*n12 + l12 + 2*n3 + l3 + 2*Ncm + Lcm
  end subroutine SetNASJacobiCMQNs

  subroutine SetThreeBodyLabQNs(this, o1, o2, o3, j12, t12)
    class(MACRO_ADD_SUFX(ThreeBodyLabQNs)), intent(inout) :: this
    type(SingleParticleOrbitIsospin), intent(in), target :: o1, o2, o3
    integer, intent(in) :: j12, t12
    this%o1 => o1
    this%o2 => o2
    this%o3 => o3
    this%j12 = j12
    this%t12 = t12
    this%Nmax = o1%e + o2%e + o3%e
  end subroutine SetThreeBodyLabQNs

  subroutine store_couplings_for_tcoefs(Nmax)
    integer, intent(in) :: Nmax
#ifdef old_tcoefs
    call LL%init(0,2*Nmax,.false., 0,2*Nmax,.false., 0,2*Nmax,.false.)
    call Ls%init(0,2*Nmax,.false., 0,2*Nmax,.false., 1,3,.true.)
    call ls2j_3%init(0,2*Nmax,.false., 0,2,.false., 0,2*Nmax,.false., 1,1,.true.)
#else
    call Lls%init(0,2*Nmax,.false., 0,2*Nmax,.false., 0,2,.false.)
    call e3o3%init(0,2*Nmax+2,.false.,0,2*Nmax,.false., 1,2*Nmax+1,.true.)
    call e2o4%init(0,2*Nmax,.false.,1,2*Nmax+1,.true., 0,2*Nmax+2,.false.)
#endif
    call ls2j_2%init(0,2*Nmax,.false., 1,1,.true., 0,2*Nmax,.false., 1,1,.true.)
    call tmbk_1%init(Nmax, 1.d0)
    call tmbk_2%init(Nmax, 2.d0)
  end subroutine store_couplings_for_tcoefs

  subroutine release_couplings_for_tcoefs()
    call tmbk_2%fin()
    call tmbk_1%fin()
    call ls2j_2%fin()
#ifdef old_tcoefs
    call ls2j_3%fin()
    call LL%fin()
    call Ls%fin()
#else
    call Lls%fin()
    call e3o3%fin()
    call e2o4%fin()
#endif
  end subroutine release_couplings_for_tcoefs

