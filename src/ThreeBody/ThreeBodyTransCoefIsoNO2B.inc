#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION
!
! For NO2B approximation related Hamiltonian
!
  public :: MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)

  private
  type :: MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)
    type(MACRO_ADD_SUFX(TCoefs)), allocatable :: Chan(:)
    integer, private, allocatable :: JacobiChannel(:)
    integer, private, allocatable :: Ncm(:)
    integer, private, allocatable :: Lcm(:)
    integer, private, allocatable :: JacNLcm2index(:,:,:)
    type(ThreeBodyJacIsoSpace), pointer :: SpJac
    type(ThreeBodyLabIsoChanNO2B), pointer :: ChLab
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), allocatable :: Lab(:)
    integer, private :: Nmax = -1
    integer, private :: NumberChannels = 0
  contains
    procedure :: InitTransJac2LabChannelIsospinNO2B
    procedure :: FinTransJac2LabChannelIsospinNO2B
    procedure :: GetJacobiIndex
    procedure :: GetNcm
    procedure :: GetLcm
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetJ12
    procedure :: GetP12
    procedure :: GetJ3
    procedure :: GetP3
    procedure :: GetNmax
    procedure :: GetMemory
    procedure :: store_couplings
    procedure :: release_couplings
    procedure :: GetLabFromIndex
    generic :: init => InitTransJac2LabChannelIsospinNO2B
    generic :: fin => FinTransJac2LabChannelIsospinNO2B
    generic :: GetLab => GetLabFromIndex
  end type MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)

  type, private :: cfp_channel
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: cfp(:,:)
#elif defined(double_precision)
    real(8), allocatable :: cfp(:,:)
#endif
  end type cfp_channel
  type(cfp_channel), allocatable, private :: jpt(:)
  logical, private :: test_mem = .false.

  type(SixJsStore), private :: Lls, e3o3, e2o4
  type(NineJsStore), private :: ls2j_2
  type(TMbracketStore), private :: tmbk_1, tmbk_2
contains

  function GetMemory(this) result(memory)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    real(8) :: memory
    integer :: ch
    memory = 0.d0
    do ch = 1, this%GetNumberChannels()
      memory = memory + this%Chan(ch)%GetMemory()
    end do
  end function GetMemory

  function GetJacobiIndex(this, i) result(idx)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%JacobiChannel(i)
  end function GetJacobiIndex

  function GetNcm(this, i) result(Ncm)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer, intent(in) :: i
    integer :: Ncm
    Ncm = this%Ncm(i)
  end function GetNcm

  function GetLcm(this, i) result(Lcm)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer, intent(in) :: i
    integer :: Lcm
    Lcm = this%Lcm(i)
  end function GetLcm

  function GetNumberChannels(this) result(NChan)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: NChan
    NChan = this%NumberChannels
  end function GetNumberChannels

  function GetIndex(this, index_jacobi, Ncm, Lcm) result(idx)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer, intent(in) :: index_jacobi, Ncm, Lcm
    integer :: idx
    idx = this%JacNLcm2Index(index_jacobi, Ncm, Lcm)
  end function GetIndex

  function GetJ(this) result(J)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: j
    j = this%ChLab%GetJ()
  end function GetJ

  function GetParity(this) result(p)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: p
    p = this%ChLab%GetParity()
  end function GetParity

  function GetT(this) result(t)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: t
    t = this%ChLab%GetT()
  end function GetT

  function GetJ12(this) result(J12)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: j12
    j12 = this%ChLab%GetJ12()
  end function GetJ12

  function GetP12(this) result(p12)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: p12
    p12 = this%ChLab%GetP12()
  end function GetP12

  function GetJ3(this) result(J3)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: j3
    j3 = this%ChLab%GetJ3()
  end function GetJ3

  function GetP3(this) result(p3)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: p3
    p3 = this%ChLab%GetP3()
  end function GetP3

  function GetNmax(this) result(Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetLabFromIndex(this, i) result(r)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in), target :: this
    integer, intent(in) :: i
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: r
    r => null()
    if(i == 0) return
    r => this%Lab(i)
  end function GetLabFromIndex

  subroutine FinTransJac2LabChannelIsospinNO2B(this)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(inout) :: this
    integer :: ch
    do ch = 1, this%GetNumberChannels()
      call this%Chan(ch)%fin()
    end do
    deallocate(this%JacNLcm2index)
    deallocate(this%Lab)
    deallocate(this%Chan)
    deallocate(this%JacobiChannel)
    deallocate(this%Ncm)
    deallocate(this%Lcm)
  end subroutine FinTransJac2LabChannelIsospinNO2B

  subroutine InitTransJac2LabChannelIsospinNO2B(this, ChLab, SpJac, e3max, n_threads_tcoef_no2b)
    use MyLibrary, only: triag
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(inout), target :: this
    type(ThreeBodyLabIsoChanNO2B), intent(in), target :: ChLab
    type(ThreeBodyJacIsoSpace), intent(in), target :: SpJac
    integer, intent(in) :: e3max, n_threads_tcoef_no2b
    integer :: NChan, ch, n_lab, idx, i1, i2, i3, t12, j12, ABCT(4)
    type(OrbitsIsospin), pointer :: sps
    type(MACRO_ADD_SUFX(TCoefs)), pointer :: tcoef
    type(sys) :: s
    real(8) :: time


    time = omp_get_wtime()
    this%ChLab => ChLab
    this%SpJac => SpJac
    this%Nmax = e3max
    sps => ChLab%sps
    call set_transjac2lab_indices(this, NChan, .false.)
    call set_transjac2lab_indices(this, NChan, .true.)
    allocate(this%Lab(this%ChLab%GetNumberStates()))
    n_lab = this%ChLab%GetNumberStates()
    do idx = 1, n_lab
      ABCT = this%ChLab%GetABCT(idx)
      i1 = ABCT(1)
      i2 = ABCT(2)
      i3 = ABCT(3)
      t12= ABCT(4)
      j12= this%ChLab%GetJ12()
      call this%Lab(idx)%set(sps%GetOrbit(i1),sps%GetOrbit(i2),sps%GetOrbit(i3),j12,t12)
    end do
    call timer%add(s%str("T coefficient, set indices"), omp_get_wtime()-time)

    do ch = 1, this%GetNumberChannels()
      tcoef => this%Chan(ch)
      tcoef%NumberLabStates = ChLab%GetNumberStates()
      allocate(tcoef%NASJacCMIndex( tcoef%GetNumberNASJacCMStates() ))
      allocate(tcoef%ASJacCMIndex( tcoef%GetNumberASJacCMStates() ))
      allocate(tcoef%LabIndex( tcoef%GetNumberLabStates() ))
      allocate(tcoef%mat(tcoef%GetNumberLabStates(), tcoef%GetNumberASJacCMStates()))
      tcoef%NASJacCMIndex(:) = 0
      tcoef%ASJacCMIndex(:) = 0
      tcoef%LabIndex(:) = 0
#if defined(single_precision) || defined(half_precision)
      tcoef%mat(:,:) = 0.0
#elif defined(double_precision)
      tcoef%mat(:,:) = 0.d0
#endif
    end do

    time = omp_get_wtime()
    call set_transjac2lab_as_lab(this, n_threads_tcoef_no2b)
    call timer%add(s%str("T coefficient"), omp_get_wtime()-time)

    if(test_mem) write(*,"(a,f12.6,a)") "Memory for tcoef(AS):  ", this%GetMemory(), " GB"
    if(this%GetMemory() > 50.d0) &
        & write(*,"(a,f12.6,a)") "Warning, may be using too much memory for tcoef(AS):  ", this%GetMemory(), " GB"
  end subroutine InitTransJac2LabChannelIsospinNO2B

  subroutine set_transjac2lab_indices(this, NChan, set_mode)
    use MyLibrary, only: triag
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(inout), target :: this
    integer, intent(inout) :: NChan
    logical, intent(in) :: set_mode
    integer :: Ncm, Lcm, ichJac, n_jac_as_nmax, n_jac_nas_nmax, NmaxJac
    type(ThreeBodyJacIsoChan), pointer :: ChJac
    type(ThreeBodyLabIsoChanNO2B), pointer :: ChLab

    if(set_mode) then
      this%NumberChannels = NChan
      allocate(this%JacNLcm2index(this%SpJac%GetNumberChannels(), 0:this%GetNmax()/2, 0:this%GetNmax()))
      allocate(this%Chan( this%GetNumberChannels() ))
      allocate(this%JacobiChannel( this%GetNumberChannels() ))
      allocate(this%Ncm( this%GetNumberChannels() ))
      allocate(this%Lcm( this%GetNumberChannels() ))
      this%JacNLcm2index(:,:,:) = 0
    end if

    ChLab => this%ChLab
    NChan = 0
    do Ncm = 0, this%GetNmax()/2
      do Lcm = 0, this%GetNmax()
        if(2*Ncm + Lcm > this%GetNmax()) cycle
        do ichJac = 1, this%SpJac%GetNumberChannels()
          ChJac => this%SpJac%GetChannel(ichjac)
          if(triag( ChJac%GetJ(), 2*Lcm, ChLab%GetJ() )) cycle
          if(ChJac%GetParity() * (-1)**Lcm /= ChLab%GetParity() ) cycle
          if(ChJac%GetT() /= ChLab%GetT()) cycle
          n_jac_as_nmax = 0
          n_jac_nas_nmax= 0
          do NmaxJac = 0, this%GetNmax() - 2*Ncm - Lcm
            if(NmaxJac > ChJac%GetNmax() ) cycle
            n_jac_as_nmax = n_jac_as_nmax + ChJac%cfp(NmaxJac)%GetNumberAStates()
            n_jac_nas_nmax= n_jac_nas_nmax+ ChJac%cfp(NmaxJac)%GetNumberNAStates()
          end do
          if(n_jac_as_nmax * n_jac_nas_nmax == 0) cycle
          NChan = NChan + 1
          if(.not. set_mode) cycle
          this%JacNLcm2index(ichJac, Ncm, Lcm) = NChan
          this%JacobiChannel( NChan ) = ichJac
          this%Ncm( NChan ) = Ncm
          this%Lcm( NChan ) = Lcm
          this%Chan(NChan )%NumberASJacCMStates = n_jac_as_nmax
          this%Chan(NChan )%NumberNASJacCMStates= n_jac_nas_nmax
        end do
      end do
    end do
  end subroutine set_transjac2lab_indices

  subroutine set_transjac2lab_as_lab(this, nthreads)
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(inout) :: this
    integer, intent(in) :: nthreads
    integer :: J, Jab, idx_channel, ichJac, Ncm, Lcm, E_cm, Jr, idx_jacobi, E_jac, idx_lab, tab, E_lab
    integer :: m1, m2, m3
    type(ThreeBodyJacIsoChan), pointer :: ChJac
    TYPE(NONANTISYMMETRIZEDISOQNS3), POINTER :: NAS
    type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc
    real(8) :: tc
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable, save :: NASLab(:,:)
#elif defined(double_precision)
    real(8), allocatable, save :: NASLab(:,:)
#endif
#ifdef half_precision
    real(4), allocatable, save :: work(:,:)
    !$omp threadprivate(work)
#endif
    !$omp threadprivate(NASLab)
    real(8) :: time
    type(sys) :: s

    J = this%ChLab%GetJ()
    Jab = this%ChLab%GetJ12()
    !$omp parallel private(idx_channel, ichJac, Ncm, Lcm, E_cm, ChJac, Jr, idx_jacobi, &
    !$omp & nas, E_jac, idx_lab, labqns, tab, E_lab, oa, ob, oc, tc, m1, m2, m3) num_threads(nthreads)
    !$omp do schedule(dynamic)
    do idx_channel = 1, this%GetNumberChannels()
      allocate( NASLab( this%ChLab%GetNumberStates(), this%Chan(idx_channel)%GetNumberNASJacCMStates()) )
#if defined(single_precision) || defined(half_precision)
      NASLab(:,:) = 0.0
#elif defined(double_precision)
      NASLab(:,:) = 0.d0
#endif
      Ncm = this%GetNcm(idx_channel)
      Lcm = this%GetLcm(idx_channel)
      E_cm = 2*Ncm + Lcm
      ichJac = this%GetJacobiIndex(idx_channel)
      ChJac => this%SpJac%GetChannel(ichJac)
      Jr = ChJac%GetJ()
      !time = omp_get_wtime()
      do idx_jacobi = 1, this%Chan(idx_channel)%GetNumberNASJacCMStates()
        nas => ChJac%GetNAS(idx_jacobi)
        E_jac = ChJac%GetENAS(idx_jacobi) + E_cm
        do idx_lab = 1, this%ChLab%GetNumberStates()
          labqns => this%GetLab( idx_lab )
          tab = labqns%t12
          E_lab = labqns%Nmax
          if(E_jac /= E_lab) cycle
          oa => labqns%o1
          ob => labqns%o2
          oc => labqns%o3
          tc = func_tcoef(oa%n, oa%l, oa%j, ob%n, ob%l, ob%j, oc%n, oc%l, oc%j, jab, tab, &
              & Ncm, Lcm, nas%GetN12(), nas%GetL12(), nas%GetS12(), nas%GetJ12(), nas%GetT12(), &
              & nas%GetN3(), nas%GetL3(), nas%GetJ3(), Jr, J)
#if defined(single_precision) || defined(half_precision)
          NASLab(idx_lab,idx_jacobi) = real(tc)
#elif defined(double_precision)
          NASLab(idx_lab,idx_jacobi) = tc
#endif
        end do
      end do
      !call timer%add(s%str("T coefficient NAS"), omp_get_wtime()-time)

      !time = omp_get_wtime()
      m1 = size(NASLab,1)
      m2 = size(NASLab,2)
      m3 = this%Chan(idx_channel)%GetNumberASJacCMStates()
#if defined(half_precision)
      allocate(work( m1,m3 ))
      call sgemm('n','n',m1,m3,m2,1.0,NASLab,m1,jpt(ichJac)%cfp(:m2,:m3),m2,0.0,work,m1)
      this%Chan(idx_channel)%mat = work
      deallocate(work)
#elif defined(single_precision)
      call sgemm('n','n',m1,m3,m2,1.0,NASLab,m1,jpt(ichJac)%cfp(:m2,:m3),m2,0.0,this%Chan(idx_channel)%mat,m1)
#elif defined(double_precision)
      call dgemm('n','n',m1,m3,m2,1.d0,NASLab,m1,jpt(ichJac)%cfp(:m2,:m3),m2,0.d0,this%Chan(idx_channel)%mat,m1)
#endif
      !call timer%add(s%str("T coefficient NAS->AS"), omp_get_wtime()-time)
      deallocate( NASLab )
    end do
    !$omp end do
    !$omp end parallel
  contains
    pure function func_tcoef(na,la,ja,nb,lb,jb,nc,lc,jc,Jab,tab,&
          & Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jrel, J) result(r)
      use MyLibrary, only: hat
      integer, intent(in) :: na,la,ja,nb,lb,jb,nc,lc,jc,jab,tab
      integer, intent(in) :: Ncm, Lcm, n12,l12,s12,j12,t12,n3,l3,j3,jrel,J
      real(8) :: r
      ! --- for tcoef ---
      integer :: lab_min, lab_max, lab, lit_min, lit_max, lit, nit
      integer :: LLam_min, LLam_max, LLam, lam_min, lam_max, lam
      real(8) :: flab, flit, fLLam, flam

      r = 0.d0
      if(t12 /= tab) return
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
              r = r + flam * fLLam * flab * flit

            end do
          end do
        end do
      end do
      r = r * hat(ja) * hat(jb) * hat(jc) * hat(2*jab) * hat(2*s12) * hat(2*j12) * hat(j3) * hat(jrel) * &
          & (-1.d0)**( s12+l12+1+j12+lc+l3+ (jc+J)/2 )
    end function func_tcoef
  end subroutine set_transjac2lab_as_lab

  subroutine store_couplings(this, jac, Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    type(ThreeBodyJacIsoSpace), intent(in) :: jac
    integer, intent(in) :: Nmax
    integer :: ch, n_phys, n_orth, i, j, N
    type(DMat) :: cfp

    call Lls%init(0,2*Nmax,.false., 0,2*Nmax,.false., 0,2,.false.)
    call e3o3%init(0,2*Nmax+2,.false., 0,2*Nmax,.false., 1,2*Nmax+1,.true., &
      &         j12dmin_in=0, j12dmax_in=2*Nmax+2, j23dmin_in=1, j23dmax_in=2*Nmax+3)
    call e2o4%init(0,2*Nmax,.false.,1,2*Nmax+1,.true., 0,2*Nmax+2,.false., &
      &         j12dmin_in=1, j12dmax_in=2*Nmax+3, j23dmin_in=1, j23dmax_in=2*Nmax+3)
    call ls2j_2%init(0,2*Nmax,.false., 1,1,.true., 0,2*Nmax,.false., 1,1,.true.)
    call tmbk_1%init(Nmax, 1.d0)
    call tmbk_2%init(Nmax, 2.d0)
    allocate(jpt( jac%GetNumberChannels() ))
    do ch = 1, jac%GetNumberChannels()
      n_phys = 0
      n_orth = 0
      do N = 0, min(Nmax,jac%jpt(ch)%GetNmax())
        n_phys = n_phys + jac%jpt(ch)%cfp(N)%GetNumberNAStates()
        n_orth = n_orth + jac%jpt(ch)%cfp(N)%GetNumberAStates()
      end do
      if(n_phys * n_orth == 0) cycle
      allocate(jpt(ch)%cfp(n_phys, n_orth))
      cfp = jac%jpt(ch)%GetCFPMat()
      do i = 1, n_orth
        do j = 1, n_phys
          jpt(ch)%cfp(j,i) = real(cfp%m(j,i), kind(jpt(ch)%cfp))
        end do
      end do
    end do
  end subroutine store_couplings

  subroutine release_couplings(this)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in) :: this
    integer :: ch
    call tmbk_2%fin()
    call tmbk_1%fin()
    call ls2j_2%fin()
    call Lls%fin()
    call e3o3%fin()
    call e2o4%fin()
    do ch = 1, size(jpt)
      if(allocated(jpt(ch)%cfp)) deallocate(jpt(ch)%cfp)
    end do
    deallocate(jpt)
  end subroutine release_couplings
