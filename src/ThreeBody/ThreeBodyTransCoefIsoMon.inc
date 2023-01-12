#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION
!
! For HF Hamiltonian
!
  public ::MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)

  private
  type :: MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)
    type(MACRO_ADD_SUFX(TCoefs)), allocatable :: Chan(:)
    integer, private, allocatable :: Indices(:,:)
    integer, private, allocatable :: JacNLJJ2Index(:,:,:,:,:)
    type(ThreeBodyJacIsoSpace), pointer :: SpJac
    type(MonThreeBodyLabIsoChan), pointer :: ChLab
    integer, private :: Nmax = -1
    integer, private :: Jmax = -1
    integer, private :: J12max = -1
    integer, private :: NumberChannels = 0
  contains
    procedure :: InitTransJac2LabChannelIsospinMon
    procedure :: FinTransJac2LabChannelIsospinMon
    procedure :: GetJ12
    procedure :: GetJ12max
    procedure :: GetJ
    procedure :: GetJmax
    procedure :: GetJacobiIndex
    procedure :: GetNcm
    procedure :: GetLcm
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetJ1
    procedure :: GetP1
    procedure :: GetJ2
    procedure :: GetP2
    procedure :: GetJ3
    procedure :: GetP3
    procedure :: GetNmax
    procedure :: GetMemory
    procedure :: store_couplings
    procedure :: release_couplings
    generic :: init => InitTransJac2LabChannelIsospinMon
    generic :: fin => FinTransJac2LabChannelIsospinMon
  end type MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)

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
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    real(8) :: memory
    integer :: ch
    memory = 0.d0
    do ch = 1, this%GetNumberChannels()
      memory = memory + this%Chan(ch)%GetMemory()
    end do
  end function GetMemory

  function GetJacobiIndex(this, i) result(idx)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%Indices(1,i)
  end function GetJacobiIndex

  function GetNcm(this, i) result(Ncm)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: i
    integer :: Ncm
    Ncm = this%Indices(2,i)
  end function GetNcm

  function GetLcm(this, i) result(Lcm)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: i
    integer :: Lcm
    Lcm = this%Indices(3,i)
  end function GetLcm

  function GetJ12(this, i) result(J12)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: i
    integer :: J12
    J12 = this%Indices(4,i)
  end function GetJ12

  function GetJ(this, i) result(J)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: i
    integer :: J
    J = this%Indices(5,i)
  end function GetJ

  function GetJ12max(this) result(J12max)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: J12max
    J12max = this%J12max
  end function GetJ12max

  function GetNmax(this) result(Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  function GetNumberChannels(this) result(NChan)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: NChan
    NChan = this%NumberChannels
  end function GetNumberChannels

  function GetIndex(this, index_jacobi, Ncm, Lcm, J12, J) result(idx)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer, intent(in) :: index_jacobi, Ncm, Lcm, J12, J
    integer :: idx
    idx = this%JacNLJJ2Index(index_jacobi, Ncm, Lcm, J12, J)
  end function GetIndex

  function GetParity(this) result(p)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: p
    p = this%ChLab%GetParity()
  end function GetParity

  function GetT(this) result(t)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: t
    t = this%ChLab%GetT()
  end function GetT

  function GetJ1(this) result(J1)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: j1
    j1 = this%ChLab%GetJ1()
  end function GetJ1

  function GetP1(this) result(p1)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: p1
    p1 = this%ChLab%GetP1()
  end function GetP1

  function GetJ2(this) result(J2)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: j2
    j2 = this%ChLab%GetJ2()
  end function GetJ2

  function GetP2(this) result(p2)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: p2
    p2 = this%ChLab%GetP2()
  end function GetP2

  function GetJ3(this) result(J3)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: j3
    j3 = this%ChLab%GetJ3()
  end function GetJ3

  function GetP3(this) result(p3)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: p3
    p3 = this%ChLab%GetP3()
  end function GetP3

  subroutine FinTransJac2LabChannelIsospinMon(this)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(inout) :: this
    integer :: ch
    do ch = 1, this%GetNumberChannels()
      call this%Chan(ch)%fin()
    end do
    deallocate(this%JacNLJJ2Index)
    if(this%GetNumberChannels()==0) return
    deallocate(this%Chan)
    deallocate(this%Indices)
  end subroutine FinTransJac2LabChannelIsospinMon

  subroutine InitTransJac2LabChannelIsospinMon(this, ChLab, SpJac, e3max, n_threads_tcoef)
    use MyLibrary, only: triag
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(inout), target :: this
    type(MonThreeBodyLabIsoChan), intent(in), target :: ChLab
    type(ThreeBodyJacIsoSpace), intent(in), target :: SpJac
    integer, intent(in) :: e3max, n_threads_tcoef
    integer :: NChans(3), NChan, ch
    type(MACRO_ADD_SUFX(TCoefs)), pointer :: tcoef
    type(sys) :: s
    real(8) :: time

    time = omp_get_wtime()
    this%ChLab => ChLab
    this%SpJac => SpJac
    this%Nmax = e3max
    NChans = count_channels(this)
    NChan = NChans(1)
    this%NumberChannels = NChan
    this%J12max = NChans(2)
    this%Jmax = NChans(3)
    allocate(this%JacNLJJ2Index(this%SpJac%GetNumberChannels(), &
        & 0:this%Nmax/2, 0:this%Nmax, 0:this%J12max, this%Jmax))
    allocate(this%Chan(NChan))
    allocate(this%Indices(5,NChan))
    call set_transjac2lab_indices(this)
    call timer%add(s%str("T coefficient, set indices"), omp_get_wtime()-time)

    do ch = 1, this%GetNumberChannels()
      tcoef => this%Chan(ch)
      tcoef%NumberLabStates = ChLab%GetNumberStates()
      allocate(tcoef%NASJacCMIndex( tcoef%GetNumberNASJacCMStates() ))
      allocate(tcoef%ASJacCMIndex( tcoef%GetNumberASJacCMStates() ))
      allocate(tcoef%LabIndex( tcoef%GetNumberLabStates() ))
      allocate(tcoef%mat( tcoef%GetNumberLabStates(), tcoef%GetNumberNASJacCMStates()))
      tcoef%NASJacCMIndex(:) = 0
      tcoef%ASJacCMIndex(:) = 0
      tcoef%LabIndex(:) = 0
#if defined(single_precision) || defined(half_precision)
      tcoef%mat = 0.0
#elif defined(double_precision)
      tcoef%mat = 0.d0
#endif
    end do

    time = omp_get_wtime()
    call set_tcoef(this, n_threads_tcoef)
    call timer%add(s%str("T coefficient"), omp_get_wtime()-time)

    if(test_mem) write(*,"(a,f12.6,a)") "Memory for tcoef(AS):  ", this%GetMemory(), " GB"
    if(this%GetMemory() > 50.d0) &
        & write(*,"(a,f12.6,a)") "Warning, may be using too much memory for tcoef(AS):  ", this%GetMemory(), " GB"
  end subroutine InitTransJac2LabChannelIsospinMon

  subroutine set_tcoef(this, n_threads_tcoef)
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(inout) :: this
    integer, intent(in) :: n_threads_tcoef
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
    real(8) :: tcoef
    integer :: idx_channel, JJ, JJab, ichJac, Ncm, Lcm, idx_jacobi, idx_lab
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc
    type(ThreeBodyJacIsoChan), pointer :: ChJac
    integer :: ABCTE(5), m1, m2, m3, Jr, E_cm, E_jac
    integer :: qns(11)
    type(sys) :: s
    !real(8) :: time

    sps => this%ChLab%sps
    !$omp parallel num_threads(n_threads_tcoef)
    !$omp do private(idx_channel, ichJac, Ncm, Lcm, JJab, JJ, E_cm, ChJac, Jr, idx_jacobi, nas, E_jac, &
    !$omp &  idx_lab, ABCTE, oa, ob, oc, tcoef, m1, m2, m3) schedule(dynamic)
    do idx_channel = 1, this%GetNumberChannels()
      allocate(NASLab( this%Chan(idx_channel)%GetNumberLabStates(), this%Chan(idx_channel)%GetNumberNASJacCMStates()))
#if defined(single_precision) || defined(half_precision)
      NASLab(:,:) = 0.0
#elif defined(double_precision)
      NASLab(:,:) = 0.d0
#endif
      ichJac = this%Indices(1,idx_channel)
      Ncm    = this%Indices(2,idx_channel)
      Lcm    = this%Indices(3,idx_channel)
      JJab   = this%Indices(4,idx_channel)
      JJ     = this%Indices(5,idx_channel)
      E_cm = 2*Ncm + Lcm
      ChJac => this%SpJac%GetChannel(ichJac)
      Jr = ChJac%GetJ()
      !time = omp_get_wtime()

      do idx_jacobi = 1, this%Chan(idx_channel)%GetNumberNASJacCMStates()
        nas => ChJac%GetNAS(idx_jacobi)
        E_jac = ChJac%GetENAS(idx_jacobi) + E_cm
        do idx_lab = 1, this%Chan(idx_channel)%GetNumberLabStates()
          ABCTE = this%ChLab%GetABCTE(idx_lab)
          oa => sps%GetOrbit(ABCTE(1))
          ob => sps%GetOrbit(ABCTE(2))
          oc => sps%GetOrbit(ABCTE(3))
          if(E_jac /= ABCTE(5)) cycle
          tcoef = func_tcoef(oa%n, oa%l, oa%j, ob%n, ob%l, ob%j, oc%n, oc%l, oc%j, JJab, ABCTE(4), &
              & Ncm, Lcm, nas%GetN12(), nas%GetL12(), nas%GetS12(), nas%GetJ12(), nas%GetT12(), &
              & nas%GetN3(), nas%GetL3(), nas%GetJ3(), Jr, JJ)

#if defined(single_precision) || defined(half_precision)
          NASLab(idx_lab,idx_jacobi) = real(tcoef)
#else
          NASLab(idx_lab,idx_jacobi) = tcoef
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
#ifdef debug
    function func_tcoef(na, la, ja, nb, lb, jb, nc, lc, jc, jab, tab, &
          & Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jrel, J) result(r)
#else
    pure function func_tcoef(na, la, ja, nb, lb, jb, nc, lc, jc, jab, tab, &
          & Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jrel, J) result(r)
#endif
      integer, intent(in) :: na, la, ja, nb, lb, jb, nc, lc, jc, jab, tab
      integer, intent(in) :: Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jrel
      integer, intent(in) :: J
      real(8) :: r
      real(8) :: flab, flit, fllam, flam
      ! loop variables
      integer :: lab, lab_min, lab_max
      integer :: lit, lit_min, lit_max, nit
      integer :: LLam, LLam_min, LLam_max
      integer :: lam, lam_min, lam_max

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
      r = r * sqrt(dble((ja+1)*(jb+1)*(jc+1)*(2*jab+1)*(2*s12+1)*(2*j12+1)*(j3+1)*(jrel+1))) * &
          & (-1.d0)**(s12+l12+1+j12+lc+l3+(jc+J)/2)
    end function func_tcoef
  end subroutine set_tcoef

  function count_channels(this) result(r)
    use MyLibrary, only: triag
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
    integer :: r(3)
    integer :: Nchan, J12max, Jmax
    integer :: j1, j2, j3, j12, j
    integer :: Ncm, Lcm, ichJac
    integer :: n_jac_as_nmax, n_jac_nas_nmax, NmaxJac
    type(ThreeBodyJacIsoChan), pointer :: ChJac
    j1 = this%GetJ1()
    j2 = this%GetJ2()
    j3 = this%GetJ3()

    NChan = 0
    J12max = -1
    Jmax = -1
    do j12 = abs(j1-j2)/2, (j1+j2)/2
      do j = abs(2*j12-j3), 2*j12+j3, 2
        do Ncm = 0, this%GetNmax()/2
          do Lcm = 0, this%GetNmax()
            if(2*Ncm + Lcm > this%GetNmax()) cycle

            do ichJac = 1, this%SpJac%GetNumberChannels()
              ChJac => this%SpJac%GetChannel(ichJac)
              if(triag( ChJac%GetJ(), 2*Lcm, J )) cycle
              if(ChJac%GetParity() * (-1)**Lcm /= this%GetParity() ) cycle
              if(ChJac%GetT() /= this%GetT() ) cycle
              n_jac_as_nmax = 0
              n_jac_nas_nmax= 0
              do NmaxJac = 0, this%GetNmax() - 2*Ncm - Lcm
                if(NmaxJac > ChJac%GetNmax() ) cycle
                n_jac_as_nmax = n_jac_as_nmax + ChJac%cfp(NmaxJac)%GetNumberAStates()
                n_jac_nas_nmax= n_jac_nas_nmax+ ChJac%cfp(NmaxJac)%GetNumberNAStates()
              end do
              if(n_jac_as_nmax * n_jac_nas_nmax == 0) cycle
              NChan = NChan + 1
              J12max = max(j12, J12max)
              Jmax = max(J, Jmax)
            end do
          end do
        end do

      end do
    end do
    r = [NChan, J12max, Jmax]
  end function count_channels

  subroutine set_transjac2lab_indices(this)
    use MyLibrary, only: triag
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(inout) :: this
    integer :: Nchan
    integer :: j1, j2, j3, j12, j
    integer :: Ncm, Lcm, ichJac
    integer :: n_jac_as_nmax, n_jac_nas_nmax, NmaxJac
    type(ThreeBodyJacIsoChan), pointer :: ChJac

    j1 = this%GetJ1()
    j2 = this%GetJ2()
    j3 = this%GetJ3()
    NChan = 0
    do j12 = abs(j1-j2)/2, (j1+j2)/2
      do j = abs(2*j12-j3), 2*j12+j3, 2
        do Ncm = 0, this%GetNmax()/2
          do Lcm = 0, this%GetNmax()
            if(2*Ncm + Lcm > this%GetNmax()) cycle

            do ichJac = 1, this%SpJac%GetNumberChannels()
              ChJac => this%SpJac%GetChannel(ichJac)
              if(triag( ChJac%GetJ(), 2*Lcm, J )) cycle
              if(ChJac%GetParity() * (-1)**Lcm /= this%GetParity() ) cycle
              if(ChJac%GetT() /= this%GetT() ) cycle
              n_jac_as_nmax = 0
              n_jac_nas_nmax= 0
              do NmaxJac = 0, this%GetNmax() - 2*Ncm - Lcm
                if(NmaxJac > ChJac%GetNmax() ) cycle
                n_jac_as_nmax = n_jac_as_nmax + ChJac%cfp(NmaxJac)%GetNumberAStates()
                n_jac_nas_nmax= n_jac_nas_nmax+ ChJac%cfp(NmaxJac)%GetNumberNAStates()
              end do
              if(n_jac_as_nmax * n_jac_nas_nmax == 0) cycle
              NChan = NChan + 1
              this%JacNLJJ2Index(ichJac, Ncm, Lcm, J12, J) = NChan
              this%Indices(:,Nchan) = [ichJac,Ncm,Lcm,J12,J]
              this%Chan(NChan )%NumberASJacCMStates = n_jac_as_nmax
              this%Chan(NChan )%NumberNASJacCMStates= n_jac_nas_nmax
            end do
          end do
        end do

      end do
    end do
  end subroutine set_transjac2lab_indices

  subroutine store_couplings(this, jac, Nmax)
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
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
    class(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in) :: this
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
