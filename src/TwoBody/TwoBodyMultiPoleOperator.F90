!
! Here the some functions are empty. I will update it once I publish a paper.
!
module TwoBodyMultiPoleOperator
  use omp_lib
  use ClassSys, only: sys, str
  use Profiler, only: timer
  use LinAlgLib
  use StoreCouplings
  use TwoBodyRelCMRadQNumbers
  use TwoBodyRelCMChanMeshBasis
  use TwoBodyRelCMChanHO
  use TwoBodyRelCMSpaceMeshBasis
  use TwoBodyRelCMSpaceHO
  use TwoBodyRelCMOps
  implicit none
  public :: TwoBodyMultiPoleOp
  public :: test_integrals
  private
  integer, parameter :: tz_phase = -1 ! particle physics convention -> nuclear physics convention
  type :: TwoBodyMultiPoleOp
    type(TwoBodyRelCMSpaceMBasis), pointer :: ms_mom => null()
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms_ho => null()
    type(str) :: OpName, external_field
    integer :: rankJ, rankP, rankZ, ch_order
    logical :: Qdependent_mesh = .true.
    real(8) :: c1, c2, c3, c4, c6, cD ! LECs
    real(8) :: Q ! momentum transfer
    integer :: NZMesh=40, NQMesh=40
  contains
    procedure :: InitTwoBodyMultiPoleOp
    procedure :: FinTwoBodyMultiPoleOp
    procedure :: show_options
    procedure :: SetHOMatrix
    generic :: init => InitTwoBodyMultiPoleOp
    generic :: fin => FinTwoBodyMultiPoleOp
  end type TwoBodyMultiPoleOp

  type :: Regulator
    logical :: reg_non_local = .false., reg_local = .false.
    integer :: reg_pow
    real(8) :: lambda
  end type Regulator
  type(Regulator) :: freg

  type :: Ystore_q
    real(8), allocatable :: y(:,:,:)
  contains
    procedure :: init_ystore_q
    procedure :: release_ystore_q
  end type Ystore_q

  type :: Ystore_lam
    type(Ystore_q), allocatable :: lam(:)
  contains
    procedure :: init_ystore_lam
    procedure :: release_ystore_lam
  end type Ystore_lam

  type :: cm_momenta
    real(8), allocatable :: P(:), W(:)
    real(8) :: pp, ww
  end type cm_momenta

  type :: angular_integral_vals
    real(8), allocatable :: v(:,:) ! nrel_bra, nrel_ket
  end type angular_integral_vals

  type :: angular_integral_lam2
    type(angular_integral_vals), allocatable :: lam2(:)
  end type angular_integral_lam2

  type :: angular_integral_lam1
    type(angular_integral_lam2), allocatable :: lam1(:)
  end type angular_integral_lam1

  type :: angular_integral_lrels
    type(angular_integral_lam1), allocatable :: lrels(:,:)
  end type angular_integral_lrels

  type :: angular_integrals
    type(angular_integral_lrels), allocatable :: args(:,:,:) ! pow of q1 (or q2), q_rel, Qex
  end type angular_integrals

  type :: PWD_functions
    real(8) :: c1, c3, c4, c6, cD
    integer :: NZMesh = 40
    integer :: NQMesh = 40
    integer :: NMesh_rel = 0
    integer :: NMesh_cm = 0
    integer :: NMeshMom = 0
    integer :: rankJ
    real(8) :: Q
    type(Ystore_lam), allocatable :: Ystore_cm(:,:)
    real(8), allocatable :: rel_mom_array(:)  ! MeV
    real(8), allocatable :: rel_wmom_array(:) ! MeV
    type(cm_momenta), allocatable :: cm_mom_array(:) ! MeV
    integer, allocatable :: idx_to_irel(:)
    integer, allocatable :: idx_to_icm(:)
    integer, allocatable :: icm_irel_to_idx(:,:)
    logical :: stored = .false.
    type(angular_integrals) :: f_ope_q1 ! ~ 1 / (q1^2 + m^2_pi)
    type(angular_integrals) :: f_ope_q2 ! ~ 1 / (q2^2 + m^2_pi)
    type(angular_integrals) :: f_pif    ! pion-in-flight ~ 1 / [(q1^2 + m^2_pi) (q2^2 + m^2_pi)]
    type(angular_integral_lrels) :: fq_unit
    real(8) :: time_pwd
    logical :: verbose
  contains
    procedure :: InitPWDFiniteQ
    procedure :: FinPWDFiniteQ
    generic :: init => InitPWDFiniteQ
    generic :: fin => FinPWDFiniteQ
  end type PWD_functions
  type(PWD_functions), save, target :: PWD

  !
  ! model-space definition for temporal LS coupled op expression
  !
  type :: LSCoupledRadial
    integer :: Lcm, lrel, Ltot, spin
    integer :: n_states
    integer, allocatable :: Ncm(:), nrel(:), Ncm_nrel2idx(:,:)
  end type LSCoupledRadial
  type :: LSCoupledChannel
    integer :: j, p, z
    type(LSCoupledRadial), allocatable :: Lcm_lrel_Ltot_spin(:)
    integer, allocatable :: Lcm_lrel_Ltot_spin2idx(:,:,:,:)
    integer :: NChan
    logical :: zero = .true.
  end type LSCoupledChannel
  type :: LSCoupledModelSpace
    type(LSCoupledChannel), allocatable :: jpz(:)
    integer, allocatable :: jpz2idx(:,:,:)
    integer :: NChan
  end type LSCoupledModelSpace
  type(LSCoupledModelSpace), target :: ms_ho_ls

  !
  ! temporally LS coupled expression
  !
  type :: OpLSChannel
    type(DMat), allocatable :: MatCh(:,:)
    type(LSCoupledChannel), pointer :: chbra, chket
    logical :: zero = .true.
  end type OpLSChannel
  type :: OpLScoupled
    type(OpLSChannel), allocatable :: OpCh(:,:)
    type(LSCoupledModelSpace), pointer :: ms
  end type OpLScoupled
  type(OpLScoupled), target :: opls
  real(8), allocatable :: coef_y_exp(:,:) ! sqrt[(2K)! / (2k1)! (2(K-k1))!]
  real(8), private :: kcmmin, kcmmax, krelmin, krelmax


contains

  subroutine InitTwoBodyMultiPoleOp(this, op, relcm_mom, c1, c3, c4, c6, cD, NZMesh, NQMesh, verbose, &
        & ch_order, regulator, reg_pow, lambda)
    use MyLibrary, only: gamma_function, hc
    class(TwoBodyMultiPoleOp), intent(inout) :: this
    type(TwoBodyRelCMOp), intent(in), target :: op
    type(TwoBodyRelCMSpaceMBasis), intent(in), target :: relcm_mom
    real(8), intent(in) :: c1, c3, c4, c6, cD
    integer, intent(in), optional :: NZMesh, NQMesh, ch_order, reg_pow
    logical, intent(in), optional :: verbose
    type(str), intent(in), optional :: regulator
    real(8), intent(in), optional :: lambda
    logical :: verb
    integer :: K, k1, k2
    type(sys) :: s
    this%OpName = op%GetOpName()
    this%Q = op%GetQ()
    this%c1 = c1
    this%c3 = c3
    this%c4 = c4
    this%c6 = c6
    this%cD = cD
    this%ms_ho => op%ms
    this%ms_mom => relcm_mom
    this%rankJ = op%GetOpJ()
    this%rankP = op%GetOpP()
    this%rankZ = op%GetOpZ()
    this%ch_order = 0
    freg%reg_pow = 2
    freg%lambda = 500.d0
    freg%reg_non_local = .false.
    freg%reg_local = .false.
    if(s%find(this%OpName,    s%str('M5_2B'))) this%external_field='axial charge'
    if(s%find(this%OpName,    s%str('L5_2B'))) this%external_field='axial vector'
    if(s%find(this%OpName,  s%str('Tel5_2B'))) this%external_field='axial vector'
    if(s%find(this%OpName, s%str('Tmag5_2B'))) this%external_field='axial vector'
    if(s%find(this%OpName,     s%str('M_2B'))) this%external_field='charge'
    if(s%find(this%OpName,     s%str('L_2B'))) this%external_field='vector'
    if(s%find(this%OpName,   s%str('Tel_2B'))) this%external_field='vector'
    if(s%find(this%OpName,  s%str('Tmag_2B'))) this%external_field='vector'
    if(present(NZMesh)) this%NZMesh = NZMesh
    if(present(NQMesh)) this%NQMesh = NQMesh
    verb = .false.
    if(present(verbose)) verb = verbose
    if(present(ch_order)) this%ch_order = ch_order
    if(present(regulator)) then
      if(regulator%val == "NonLocal") freg%reg_non_local = .true.
      if(regulator%val == "Local") freg%reg_local = .true.
    end if
    if(present(reg_pow)) freg%reg_pow = reg_pow
    if(present(lambda)) freg%lambda = lambda
    if(verb) call this%show_options()
    call PWD%init(relcm_mom, op%ms%GetNmax(), this%NZMesh, this%NQMesh, this%rankJ, this%Q, &
        & this%c1, this%c3, this%c4, this%c6, this%cD, this%external_field, this%Qdependent_mesh, verb)
    allocate(coef_y_exp(0:4,0:4))
    do K = 0, 4
      do k1 = 0, K
        k2 = K-k1
        coef_y_exp(K,k1) = sqrt( gamma_function(dble(2*K+1)) / (gamma_function(dble(2*k1+1))*gamma_function(dble(2*k2+1))))
      end do
    end do
  end subroutine InitTwoBodyMultiPoleOp

  subroutine show_options(this)
    class(TwoBodyMultiPoleOp), intent(in) :: this
    write(*,"(a,f8.2,a)") "# Calculation options for 2b multipole operator @ Q=", this%Q, " MeV"
    write(*,"(2a)") "# Operator: ", trim(this%OpName%val)
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a)") "# c1= ", this%c1, " GeV-1, c2= ", &
      & this%c2, " GeV-1, c3=", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    write(*,"(a,f8.4,a,f8.4)") "# c6= ", this%c6, ", cd= ", this%cd
    write(*,"(a)") "#"
    write(*,"(a)") "# Truncations: "
    write(*,"(a,i4,a,i4)") "# NMesh (rel) = ", this%ms_mom%GetNMeshRel(), ", NMesh (cm) = ", this%ms_mom%GetNMeshCM()
    write(*,"(a,f8.4,a,f8.4,a)") "# pmin = ", this%ms_mom%GetPrelMin(), " fm-1, pmax = ", this%ms_mom%GetPrelMax(), " fm-1"
    write(*,"(a,f8.4,a,f8.4,a)") "# Pmin = ", this%ms_mom%GetPcmMin(), " fm-1, Pmax = ", this%ms_mom%GetPcmMax(), " fm-1"
    write(*,"(a,i4,a,i4)") "# NMesh for theta = ", this%NZMesh, ", NMesh for q = ", this%NQMesh
    write(*,"(a,i4,a,i4,a,i4)") "# Nmax (tot) = ", this%ms_ho%GetNMax(), ", Nmax (rel) = ", this%ms_ho%GetNmaxRel(), &
        & ", Nmax (cm) = ", this%ms_ho%GetNmaxCM()
    write(*,"(a,i4,a,i4,a,i4)") "# Jmax (tot) = ", this%ms_mom%GetJmax(), ", Jmax (rel) = ", this%ms_mom%GetJrelMax(), &
        & ", Lcm max = ", this%ms_mom%GetLcmMax()
    if(freg%reg_local) write(*,"(a,f8.4,a,i4)") "# Local regulator: Lambda= ", freg%lambda, " MeV, n= ", freg%reg_pow
    if(freg%reg_non_local) write(*,"(a,f8.4,a,i4)") "# Non-local regulator: Lambda= ", freg%lambda, " MeV, n= ", freg%reg_pow
    write(*,"(a)") "#"
  end subroutine show_options

  subroutine FinTwoBodyMultiPoleOp(this)
    class(TwoBodyMultiPoleOp), intent(inout) :: this
    deallocate(coef_y_exp)
    call PWD%fin()
  end subroutine FinTwoBodyMultiPoleOp

  function set_momentum_matrix(this, func, chbra, chket, Jbra, Jket, zbra, zket) result(res)
    use MyLibrary, only: hc
    interface
      function func(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res_func)
        use LinAlgLib, only: DMat
        integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
        integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
        type(DMat) :: res_func
      end function func
    end interface
    class(TwoBodyMultiPoleOp), intent(in) :: this
    type(LSCoupledRadial) :: chbra, chket
    integer, intent(in) :: Jbra, Jket, Zbra, Zket
    integer :: lcm_bra, lrel_bra, s_bra, Lt_bra
    integer :: lcm_ket, lrel_ket, s_ket, Lt_ket
    real(8) :: phase
    type(DMat) :: res

    lcm_bra = chbra%Lcm; lrel_bra = chbra%Lrel; Lt_bra = chbra%Ltot; s_bra = chbra%spin
    lcm_ket = chket%Lcm; lrel_ket = chket%Lrel; Lt_ket = chket%Ltot; s_ket = chket%spin
    phase = 1.d0
    if(mod(abs(lrel_bra+lcm_bra-lrel_ket-lcm_ket+this%rankJ),2)==0) &
        & phase = (-1.d0)**((lrel_bra+lcm_bra-lrel_ket-lcm_ket+this%rankJ)/2)
    if(mod(abs(lrel_bra+lcm_bra-lrel_ket-lcm_ket+this%rankJ),2)==1) &
        & phase = (-1.d0)**((lrel_bra+lcm_bra-lrel_ket-lcm_ket+this%rankJ-1)/2) ! in unit of i
    phase = (-1.d0)**((lrel_ket+lcm_ket-lrel_bra-lcm_bra)/2) * (-1.d0)**(pwd%rankJ/2)
    call res%ini(chbra%n_states, chket%n_states)
    res = func(this%rankJ, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res * phase * hc**6
  end function set_momentum_matrix

  subroutine SetHOMatrix(this, op, term)
    use MyLibrary, only: triag
    abstract interface
      function func(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res_func)
        use LinAlgLib, only: DMat
        integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
        integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
        type(DMat) :: res_func
      end function func
    end interface
    class(TwoBodyMultiPoleOp), intent(inout) :: this
    type(TwoBodyRelCMOp), intent(inout) :: op
    type(str), intent(in) :: term
    type(LSCoupledRadial) :: chbra_ls, chket_ls
    type(LSCoupledRadial), pointer :: bra_ho, ket_ho
    integer :: ichbra, ichket, ichbra_ls, ichket_ls
    integer :: Jbra, Jket, Zbra, Zket, ichket_ls_max, iloop, cnt
    integer, allocatable :: iloops(:,:)
    type(DMat) :: mat
    type(sys) :: s
    type(str) :: opname
    procedure(func), pointer :: f_ptr
    real(8) :: time, time_set_mat, time_mom2ho, t1

    opname = this%OpName
    f_ptr => null()
    if(this%ch_order < 1) then
      write(*,*) "There is no contribution at LO!"
      return
    end if
    if(this%external_field%val=='vector') then
      if(this%ch_order > 2) write(*,*) "N3LO vector current is not implemented!"
    end if
    if(this%external_field%val=='axial vector') then
      if(this%ch_order < 2) then
        write(*,*) "There is no axial vector contribution at NLO!"
        return
       end if
      if(this%ch_order > 2) write(*,*) "N3LO vector current is not implemented!"
    end if
    if(this%external_field%val=='charge') then
      if(this%ch_order > 2) write(*,*) "N3LO charge is not implemented!"
      return
    end if
    if(this%external_field%val=='axial charge') then
      write(*,*) "Any axial charge operator is not implemented!"
      return
    end if
    if(s%find(opname, s%str('L5_2B')) .and. s%find(term, s%str('c1'))) f_ptr => L_axial_vector_n2lo_c1
    if(s%find(opname, s%str('L5_2B')) .and. s%find(term, s%str('c3'))) f_ptr => L_axial_vector_n2lo_c3
    if(s%find(opname, s%str('Tel5_2B')) .and. s%find(term, s%str('c3'))) f_ptr => Tel_axial_vector_n2lo_c3
    if(s%find(opname, s%str('Tmag5_2B')) .and. s%find(term, s%str('c3'))) f_ptr => Tmag_axial_vector_n2lo_c3
    if(s%find(opname, s%str('L5_2B')) .and. s%find(term, s%str('c4'))) f_ptr => L_axial_vector_n2lo_c4
    if(s%find(opname, s%str('Tel5_2B')) .and. s%find(term, s%str('c4'))) f_ptr => Tel_axial_vector_n2lo_c4
    if(s%find(opname, s%str('Tmag5_2B')) .and. s%find(term, s%str('c4'))) f_ptr => Tmag_axial_vector_n2lo_c4
    if(s%find(opname, s%str('L5_2B')) .and. s%find(term, s%str('cD'))) f_ptr => L_axial_vector_n2lo_cD
    if(s%find(opname, s%str('Tel5_2B')) .and. s%find(term, s%str('cD'))) f_ptr => Tel_axial_vector_n2lo_cD
    if(s%find(opname, s%str('Tmag5_2B')) .and. s%find(term, s%str('cD'))) f_ptr => Tmag_axial_vector_n2lo_cD
    if(s%find(opname, s%str('Tel5_2B')) .and. s%find(term, s%str('c6'))) f_ptr => Tel_axial_vector_n2lo_c6
    if(s%find(opname, s%str('Tmag5_2B')) .and. s%find(term, s%str('c6'))) f_ptr => Tmag_axial_vector_n2lo_c6
    if(s%find(opname, s%str('L_2B'))) f_ptr => L_vector_nlo
    if(s%find(opname, s%str('Tel_2B'))) f_ptr => Tel_vector_nlo
    if(s%find(opname, s%str('Tmag_2B'))) f_ptr => Tmag_vector_nlo
    if(s%find(opname, s%str('PIFTmag_2B'))) f_ptr => Tmag_vector_nlo_pion_in_flight
    if(s%find(opname, s%str('SGTmag_2B'))) f_ptr => Tmag_vector_nlo_seagull
    if(s%find(opname, s%str('PIFTel_2B'))) f_ptr => Tel_vector_nlo_pion_in_flight
    if(s%find(opname, s%str('SGTel_2B'))) f_ptr => Tel_vector_nlo_seagull

    do ichbra = 1, op%ms%GetNumberChannels()
      do ichket = 1, ichbra
        if(op%MatCh(ichbra,ichket)%is_zero) cycle
        op%MatCh(ichbra,ichket)%m(:,:) = 0.d0
      end do
    end do
    if(term%val == "c1" .and. abs(pwd%c1)<1.d-8) return
    if(term%val == "c3" .and. abs(pwd%c3)<1.d-8) return
    if(term%val == "c4" .and. abs(pwd%c4)<1.d-8) return
    if(term%val == "c6" .and. abs(pwd%c6)<1.d-8) return
    if(term%val == "cD" .and. abs(pwd%cD)<1.d-8) return
    if(.not. associated(f_ptr)) then
      write(*,*) "The term ", trim(term%val), ' of ', trim(OpName%val), ' vanishes; I do nothing here.'
      return
    end if

    call init_lscoupled_space(op%ms)
    call init_op_lscoupled(this)
    do iloop = 1, 2
      cnt = 0
      do ichbra = 1, opls%ms%NChan
        do ichket = 1, ichbra
          if(opls%OpCh(ichbra,ichket)%zero) cycle
          do ichbra_ls = 1, opls%ms%jpz(ichbra)%NChan
            ichket_ls_max = opls%ms%jpz(ichket)%NChan
            if(ichbra==ichket) ichket_ls_max = ichbra_ls
            do ichket_ls = 1, ichket_ls_max
              cnt = cnt + 1
              if(iloop==2) iloops(:,cnt) = [ichbra,ichket,ichbra_ls,ichket_ls]
            end do
          end do
        end do
      end do
      if(iloop==1) allocate(iloops(4,cnt))
    end do

    if(pwd%verbose) write(*,'(a,i8,a)') 'Calculations will be done for ', size(iloops,2), ' channels'
    time = omp_get_wtime()
    time_set_mat = 0.d0
    time_mom2ho = 0.d0

    !$omp parallel
    !$omp do private(iloop, ichbra, ichket, ichbra_ls, ichket_ls, Jbra, Jket, Zbra, Zket, &
    !$omp & bra_ho, ket_ho, chbra_ls, chket_ls, mat) schedule(dynamic)
    do iloop = 1, size(iloops,2)
      ichbra = iloops(1,iloop)
      ichket = iloops(2,iloop)
      ichbra_ls = iloops(3,iloop)
      ichket_ls = iloops(4,iloop)
      Jbra = opls%ms%jpz(ichbra)%j
      Jket = opls%ms%jpz(ichket)%j
      Zbra = opls%ms%jpz(ichbra)%Z
      Zket = opls%ms%jpz(ichket)%Z
      bra_ho => ms_ho_ls%jpz(ichbra)%Lcm_lrel_Ltot_spin(ichbra_ls)
      ket_ho => ms_ho_ls%jpz(ichket)%Lcm_lrel_Ltot_spin(ichket_ls)
      call init_lscoupled_channel_mom(chbra_ls, bra_ho%Lcm, bra_ho%lrel, bra_ho%Ltot, bra_ho%spin)
      call init_lscoupled_channel_mom(chket_ls, ket_ho%Lcm, ket_ho%lrel, ket_ho%Ltot, ket_ho%spin)
      t1 = omp_get_wtime()
      mat = set_momentum_matrix(this, f_ptr, chbra_ls, chket_ls, Jbra, Jket, Zbra, Zket)
      time_set_mat = time_set_mat + omp_get_wtime()-t1
      t1 = omp_get_wtime()
      call transform_to_ho_ls_channel(opls%OpCh(ichbra,ichket)%MatCh(ichbra_ls,ichket_ls), bra_ho, ket_ho, &
          & mat, chbra_ls, chket_ls, Zbra, Zket, op%ms%GetFrequency())
      !call transform_to_ho_ls_channel_old(opls%OpCh(ichbra,ichket)%MatCh(ichbra_ls,ichket_ls), bra_ho, ket_ho, &
      !    & mat, chbra_ls, chket_ls, Zbra, Zket, op%ms%GetFrequency())
      time_mom2ho = time_mom2ho + omp_get_wtime()-t1
      opls%OpCh(ichket,ichbra)%MatCh(ichket_ls,ichbra_ls) = &
          & opls%OpCh(ichbra,ichket)%MatCh(ichbra_ls,ichket_ls)%T() * (-1.d0)**(Jbra-Jket)

      call release_lscoupled_channel_mom(chket_ls)
      call release_lscoupled_channel_mom(chbra_ls)
    end do
    !$omp end do
    !$omp end parallel

    call timer%add(s%str('Mom. 2b multipole operator'), time_set_mat)
    call timer%add(s%str('Mom -> HO 2b multipole operator'), time_mom2ho)
    !call timer%add(s%str('ls coupled op in 2b multipole operator'), omp_get_wtime()-time)
    deallocate(iloops)

    time = omp_get_wtime()
    !$omp parallel
    !$omp do private(ichbra, ichket) schedule(dynamic)
    do ichbra = 1, op%ms%GetNumberChannels()
      do ichket = 1, op%ms%GetNumberChannels()
        if(op%MatCh(ichbra,ichket)%is_zero) cycle
        op%MatCh(ichbra,ichket)%m(:,:) = 0.d0
        call ls_2_jj(op%MatCh(ichbra,ichket))
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%add(s%str('ls -> jj 2b multipole operator'), omp_get_wtime()-time)

    call release_op_lscoupled()
    call release_lscoupled_space()
  end subroutine SetHOMatrix

  subroutine ls_2_jj(opjj)
    use MyLibrary, only: sjs
    type(TwoBodyRelCMOpChan), intent(inout) :: opjj
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    integer :: ibra, ncm_bra, lcm_bra, nrel_bra, lrel_bra, spin_bra, jrel_bra, Lt_bra, Jbra, Pbra, Zbra
    integer :: iket, ncm_ket, lcm_ket, nrel_ket, lrel_ket, spin_ket, jrel_ket, Lt_ket, Jket, Pket, Zket
    integer :: ls_chan_bra, ls_chan_ket, ls_ltot_bra, ls_ltot_ket, ls_idx_bra, ls_idx_ket
    real(8) :: me

    chbra => opjj%rel_ch_bra
    chket => opjj%rel_ch_ket
    Jbra = chbra%GetJ(); Pbra = chbra%GetParity(); Zbra = chbra%GetZ(); lcm_bra = chbra%GetLcm(); jrel_bra = chbra%GetJRel()
    Jket = chket%GetJ(); Pket = chket%GetParity(); Zket = chket%GetZ(); lcm_ket = chket%GetLcm(); jrel_ket = chket%GetJRel()
    ls_chan_bra = opls%ms%jpz2idx(Jbra,Pbra,Zbra)
    ls_chan_ket = opls%ms%jpz2idx(Jket,Pket,Zket)
    do ibra = 1, chbra%GetNumberStates()
      hobra => chbra%GetHOs(ibra)
      ncm_bra = hobra%GetNCM()
      nrel_bra = hobra%GetNRel()
      lrel_bra = hobra%GetLRel()
      spin_bra = hobra%GetSpin()
      do iket = 1, chket%GetNumberStates()
        hoket => chket%GetHOs(iket)
        ncm_ket = hoket%GetNCM()
        nrel_ket = hoket%GetNRel()
        lrel_ket = hoket%GetLRel()
        spin_ket = hoket%GetSpin()

        me = 0.d0
        do Lt_bra = max(abs(Jbra-spin_bra), abs(lcm_bra-lrel_bra)), min(Jbra+spin_bra, lcm_bra+lrel_bra)
          do Lt_ket = max(abs(Jket-spin_ket), abs(lcm_ket-lrel_ket)), min(Jket+spin_ket, lcm_ket+lrel_ket)
            ls_ltot_bra = opls%ms%jpz(ls_chan_bra)%Lcm_lrel_Ltot_spin2idx(lcm_bra,lrel_bra,Lt_bra,spin_bra)
            ls_ltot_ket = opls%ms%jpz(ls_chan_ket)%Lcm_lrel_Ltot_spin2idx(lcm_ket,lrel_ket,Lt_ket,spin_ket)
            if(ls_ltot_bra * ls_ltot_ket==0) cycle
            ls_idx_bra = opls%ms%jpz(ls_chan_bra)%Lcm_lrel_Ltot_spin(ls_ltot_bra)%Ncm_nrel2idx(ncm_bra,nrel_bra)
            ls_idx_ket = opls%ms%jpz(ls_chan_ket)%Lcm_lrel_Ltot_spin(ls_ltot_ket)%Ncm_nrel2idx(ncm_ket,nrel_ket)
            if(ls_idx_bra * ls_idx_ket==0) cycle

            me = me + (-1.d0)**(lcm_bra+lrel_bra+spin_bra+Jbra+ lcm_ket+lrel_ket+spin_ket+Jket) * &
                & sjs(2*lcm_bra, 2*lrel_bra, 2*Lt_bra, 2*spin_bra, 2*Jbra, 2*jrel_bra) * &
                & sjs(2*lcm_ket, 2*lrel_ket, 2*Lt_ket, 2*spin_ket, 2*Jket, 2*jrel_ket) * &
                & sqrt(dble((2*Lt_bra+1)*(2*jrel_bra+1)*(2*Lt_ket+1)*(2*jrel_ket+1))) * &
                & opls%OpCh(ls_chan_bra,ls_chan_ket)%MatCh(ls_ltot_bra,ls_ltot_ket)%m(ls_idx_bra,ls_idx_ket)

          end do
        end do
        opjj%m(ibra,iket) = me

      end do
    end do
  end subroutine ls_2_jj

  subroutine transform_to_ho_ls_channel(mat_ho, bra_ho, ket_ho, mat_mom, bra_mom, ket_mom, Zbra, Zket, hw)
    use MyLibrary, only: hc, m_proton, m_neutron, ho_radial_wf_norm, gauss_legendre
    use NdSpline
    type(Dmat), intent(inout) :: mat_ho
    type(LSCoupledRadial), intent(in) :: bra_ho, ket_ho
    type(DMat), intent(in) :: mat_mom
    type(LSCoupledRadial), intent(in) :: bra_mom, ket_mom
    integer, intent(in) :: Zbra, Zket
    real(8), intent(in) :: hw
    integer :: ibra_ho, iket_ho, ibra_mom, iket_mom
    integer :: nrel_bra, lrel_bra, nrel_ket, lrel_ket, ncm_bra, lcm_bra, ncm_ket, lcm_ket
    integer :: irel_bra, irel_ket, icm_bra, icm_ket, i
    real(8) :: me, mrel_bra, mrel_ket, mcm_bra, mcm_ket, arel_bra, arel_ket, acm_bra, acm_ket
    real(8) :: time, norm
    type(sys) :: s
    real(8), allocatable :: radial_rel1(:,:), radial_rel2(:,:), radial_cm1(:,:), radial_cm2(:,:)
    type(spline) :: sp
    integer, parameter :: k = 4
    integer, parameter :: NMesh = 80
    logical, parameter :: extrap_near_edge = .true.
    real(8) :: integral
    real(8), allocatable :: tmp1d(:), kcm_mesh(:), wcm_mesh(:), krel_mesh(:), wrel_mesh(:)
    real(8), allocatable :: Bcm(:,:), Brel(:,:), Mat_mid(:,:,:,:), Mat(:,:,:), tmp3d(:,:,:)
    real(8), allocatable :: work1(:,:,:), work2(:,:,:)

    time = omp_get_wtime()
    mrel_bra = (m_proton*m_neutron)/(m_proton+m_neutron)
    mcm_bra = m_proton + m_neutron
    mrel_ket = (m_proton*m_neutron)/(m_proton+m_neutron)
    mcm_ket = m_proton + m_neutron
    arel_bra = hc**2 / (mrel_bra*hw)
    arel_ket = hc**2 / (mrel_ket*hw)
    acm_bra =  hc**2 / (mcm_bra *hw)
    acm_ket =  hc**2 / (mcm_ket *hw)

    lcm_bra = bra_ho%lcm
    lrel_bra = bra_ho%lrel
    lcm_ket = ket_ho%lcm
    lrel_ket = ket_ho%lrel

    norm = 0.d0
    do ibra_mom = 1, pwd%NMeshMom
      do iket_mom = 1, pwd%NMeshMom
        norm = norm + mat_mom%m(ibra_mom, iket_mom)**2
      end do
    end do
    if(sqrt(norm) < 1.d-24) return

    allocate(radial_cm2(pwd%NMesh_cm, 0:maxval(ket_ho%Ncm)))
    allocate(Mat_mid(pwd%NMesh_rel, pwd%NMesh_cm, pwd%NMesh_rel, 0:maxval(ket_ho%Ncm)))
    allocate(tmp1d(pwd%NMesh_cm), kcm_mesh(NMesh), wcm_mesh(NMesh))
    ! This part is expensive. Optimize here
    do icm_bra = 1, pwd%NMesh_cm

      call gauss_legendre(abs(pwd%cm_mom_array(icm_bra)%pp-pwd%Q), pwd%cm_mom_array(icm_bra)%pp+pwd%Q, &
        & kcm_mesh, wcm_mesh, NMesh)
      tmp1d(:) = 0.d0
      call sp%init([k], [pwd%NMesh_cm], pwd%cm_mom_array(icm_bra)%p, tmp1d, extrapolation=extrap_near_edge)
      Bcm = sp%Get_Bmatrix(1,kcm_mesh)
      call sp%fin()
      radial_cm2(:,:) = 0.d0
      do ncm_ket = 0, maxval(ket_ho%Ncm)
        do i = 1, NMesh
          integral = ho_radial_wf_norm(ncm_ket, lcm_ket, acm_ket, kcm_mesh(i)/hc) * (-1.d0)**ncm_ket * &
            & kcm_mesh(i) * wcm_mesh(i) / hc**2
          do icm_ket = 1, pwd%NMesh_cm
            radial_cm2(icm_ket, ncm_ket) = radial_cm2(icm_ket, ncm_ket) + integral * Bcm(icm_ket, i)
          end do
        end do
      end do

      do irel_bra = 1, pwd%NMesh_rel
        do irel_ket = 1, pwd%NMesh_rel
          do icm_ket = 1, pwd%NMesh_cm
            ibra_mom = pwd%icm_irel_to_idx(icm_bra,irel_bra)
            iket_mom = pwd%icm_irel_to_idx(icm_ket,irel_ket)
            tmp1d(icm_ket) = mat_mom%m(ibra_mom, iket_mom)
          end do
          call sp%init([k], [pwd%NMesh_cm], pwd%cm_mom_array(icm_bra)%p, tmp1d, extrapolation=extrap_near_edge)
          tmp1d = sp%get_coefs()
          call sp%fin()
          Mat_mid(irel_bra, icm_bra, irel_ket, :) = matmul(tmp1d, radial_cm2)
        end do
      end do
    end do
    deallocate(tmp1d, kcm_mesh, wcm_mesh)

    allocate(radial_rel1(pwd%NMesh_rel, 0:maxval(bra_ho%Nrel)))
    allocate(radial_rel2(pwd%NMesh_rel, 0:maxval(ket_ho%Nrel)))
    allocate(radial_cm1(pwd%NMesh_cm, 0:maxval(bra_ho%Ncm)))
    allocate(tmp1d(pwd%NMesh_rel*pwd%NMesh_cm*pwd%NMesh_rel))
    allocate(tmp3d(pwd%NMesh_rel,pwd%NMesh_cm,pwd%NMesh_rel))
    allocate(work1(0:maxval(bra_ho%Nrel), pwd%NMesh_cm, pwd%NMesh_rel))
    allocate(work2(0:maxval(bra_ho%Nrel), 0:maxval(bra_ho%Ncm),pwd%NMesh_rel))
    allocate(Mat(0:maxval(bra_ho%Nrel), 0:maxval(bra_ho%Ncm), 0:maxval(ket_ho%Nrel)))
    call gauss_legendre(kcmmin, kcmmax, kcm_mesh, wcm_mesh, NMesh)
    call gauss_legendre(krelmin, krelmax, krel_mesh, wrel_mesh, NMesh)
    tmp1d(:) = 0.d0
    call sp%init([k,k,k], [pwd%NMesh_rel, pwd%NMesh_cm, pwd%NMesh_rel], &
      & [pwd%rel_mom_array, pwd%cm_mom_array(:)%pp, pwd%rel_mom_array], tmp1d, extrapolation=extrap_near_edge)
    Bcm = sp%Get_Bmatrix(1,kcm_mesh)
    Brel = sp%Get_Bmatrix(2,krel_mesh)
    call sp%fin()

    radial_cm1(:,:) = 0.d0
    do ncm_bra = 0, maxval(bra_ho%Ncm)
      do i = 1, NMesh
        integral = ho_radial_wf_norm(ncm_bra, lcm_bra, acm_bra, kcm_mesh(i)/hc) * (-1.d0)**ncm_bra * &
          & kcm_mesh(i) * wcm_mesh(i) /hc**2
        do icm_bra = 1, pwd%NMesh_cm
          radial_cm1(icm_bra, ncm_bra) = radial_cm1(icm_bra, ncm_bra) + integral * Bcm(icm_bra,i)
        end do
      end do
    end do

    radial_rel1(:,:) = 0.d0
    do nrel_bra = 0, maxval(bra_ho%Nrel)
      do i = 1, NMesh
        integral = ho_radial_wf_norm(nrel_bra, lrel_bra, arel_bra, krel_mesh(i)/hc) * (-1.d0)**nrel_bra * &
            & krel_mesh(i) * wrel_mesh(i) /hc**2
        if(freg%reg_non_local) integral = integral * regulator_function(krel_mesh(i))
        do irel_bra = 1, pwd%NMesh_rel
          radial_rel1(irel_bra, nrel_bra) = radial_rel1(irel_bra, nrel_bra) + integral * Brel(irel_bra,i)
        end do
      end do
    end do

    radial_rel2(:,:) = 0.d0
    do nrel_ket = 0, maxval(ket_ho%Nrel)
      do i = 1, NMesh
        integral = ho_radial_wf_norm(nrel_ket, lrel_ket, arel_ket, krel_mesh(i)/hc) * (-1.d0)**nrel_ket * &
            & krel_mesh(i) * wrel_mesh(i) /hc**2
        if(freg%reg_non_local) integral = integral * regulator_function(krel_mesh(i))
        do irel_ket = 1, pwd%NMesh_rel
          radial_rel2(irel_ket, nrel_ket) = radial_rel2(irel_ket, nrel_ket) + integral * Brel(irel_ket,i)
        end do
      end do
    end do

    do ncm_ket = 0, maxval(ket_ho%Ncm)
      tmp1d = reshape(Mat_mid(:,:,:,ncm_ket), shape(tmp1d))
      call sp%init([k,k,k], [pwd%NMesh_rel,pwd%NMesh_cm,pwd%NMesh_rel], &
        & [pwd%rel_mom_array, pwd%cm_mom_array(:)%pp, pwd%rel_mom_array], tmp1d, extrapolation=extrap_near_edge)
      tmp3d = reshape(sp%get_coefs(), shape(tmp3d))
      call sp%fin()

      do nrel_bra = 0, maxval(bra_ho%Nrel)
        do icm_bra = 1, pwd%NMesh_cm
          do irel_ket = 1, pwd%NMesh_rel
            work1(nrel_bra,icm_bra,irel_ket) = dot_product(radial_rel1(:,nrel_bra),tmp3d(:,icm_bra,irel_ket))
          end do
        end do
      end do

      do nrel_bra = 0, maxval(bra_ho%Nrel)
        do ncm_bra = 0, maxval(bra_ho%Ncm)
          do irel_ket = 1, pwd%NMesh_rel
            work2(nrel_bra,ncm_bra,irel_ket) = dot_product(radial_cm1(:,ncm_bra),work1(nrel_bra,:,irel_ket))
          end do
        end do
      end do

      do nrel_bra = 0, maxval(bra_ho%Nrel)
        do ncm_bra = 0, maxval(bra_ho%Ncm)
          do nrel_ket = 0, maxval(ket_ho%Nrel)
            Mat(nrel_bra,ncm_bra,nrel_ket) = dot_product(radial_rel2(:,nrel_ket),work2(nrel_bra,ncm_bra,:))
          end do
        end do
      end do

      do ncm_bra = 0, maxval(bra_ho%Ncm)
        do nrel_ket = 0, maxval(ket_ho%Nrel)
          do nrel_bra = 0, maxval(bra_ho%Nrel)
            ibra_ho = bra_ho%Ncm_nrel2idx(Ncm_bra,Nrel_bra)
            iket_ho = ket_ho%Ncm_nrel2idx(Ncm_ket,Nrel_ket)
            if(ibra_ho * iket_ho < 1) cycle
            mat_ho%m(ibra_ho,iket_ho) = Mat(nrel_bra,ncm_bra,nrel_ket)
          end do
        end do
      end do
    end do

    deallocate(tmp1d)
    deallocate(tmp3d)
    deallocate(work1)
    deallocate(work2)
    deallocate(Mat)
    deallocate(radial_rel1)
    deallocate(radial_rel2)
    deallocate(radial_cm1)
    deallocate(radial_cm2)
  end subroutine transform_to_ho_ls_channel

  subroutine transform_to_ho_ls_channel_old(mat_ho, bra_ho, ket_ho, mat_mom, bra_mom, ket_mom, Zbra, Zket, hw)
    use MyLibrary, only: hc, m_proton, m_neutron, ho_radial_wf_norm
    type(Dmat), intent(inout) :: mat_ho
    type(LSCoupledRadial), intent(in) :: bra_ho, ket_ho
    type(DMat), intent(in) :: mat_mom
    type(LSCoupledRadial), intent(in) :: bra_mom, ket_mom
    integer, intent(in) :: Zbra, Zket
    real(8), intent(in) :: hw
    integer :: ibra_ho, iket_ho, ibra_mom, iket_mom
    integer :: nrel_bra, lrel_bra, nrel_ket, lrel_ket, ncm_bra, lcm_bra, ncm_ket, lcm_ket
    integer :: irel_bra, irel_ket, icm_bra, icm_ket
    real(8) :: me, mrel_bra, mrel_ket, mcm_bra, mcm_ket, arel_bra, arel_ket, acm_bra, acm_ket
    real(8) :: time
    type(sys) :: s
    real(8), allocatable :: radial_rel1(:,:), radial_rel2(:,:), radial_cm1(:,:), radial_cm2(:,:,:)

    time = omp_get_wtime()
    mrel_bra = (m_proton*m_neutron)/(m_proton+m_neutron)
    mcm_bra = m_proton + m_neutron
    mrel_ket = (m_proton*m_neutron)/(m_proton+m_neutron)
    mcm_ket = m_proton + m_neutron
    arel_bra = hc**2 / (mrel_bra*hw)
    arel_ket = hc**2 / (mrel_ket*hw)
    acm_bra =  hc**2 / (mcm_bra *hw)
    acm_ket =  hc**2 / (mcm_ket *hw)

    lcm_bra = bra_ho%lcm
    lrel_bra = bra_ho%lrel
    lcm_ket = ket_ho%lcm
    lrel_ket = ket_ho%lrel

    allocate(radial_rel1(pwd%NMesh_rel, 0:maxval(bra_ho%Nrel)))
    allocate(radial_rel2(pwd%NMesh_rel, 0:maxval(ket_ho%Nrel)))
    allocate(radial_cm1(pwd%NMesh_cm, 0:maxval(bra_ho%Ncm)))
    allocate(radial_cm2(pwd%NMesh_cm, pwd%NMesh_cm, 0:maxval(ket_ho%Ncm)))
    do nrel_bra = 0, maxval(bra_ho%Nrel)
      do irel_bra = 1, pwd%NMesh_rel
        radial_rel1(irel_bra, nrel_bra) = &
            & ho_radial_wf_norm(nrel_bra, lrel_bra, arel_bra, pwd%rel_mom_array(irel_bra)/hc) * (-1.d0)**nrel_bra * &
            & pwd%rel_mom_array(irel_bra) * pwd%rel_wmom_array(irel_bra) /hc**2
      end do
    end do
    do nrel_ket = 0, maxval(ket_ho%Nrel)
      do irel_ket = 1, pwd%NMesh_rel
        radial_rel2(irel_ket, nrel_ket) = &
            & ho_radial_wf_norm(nrel_ket, lrel_ket, arel_ket, pwd%rel_mom_array(irel_ket)/hc) * (-1.d0)**nrel_ket * &
            & pwd%rel_mom_array(irel_ket) * pwd%rel_wmom_array(irel_ket) /hc**2
      end do
    end do
    do ncm_bra = 0, maxval(bra_ho%Ncm)
      do icm_bra = 1, pwd%NMesh_cm
        radial_cm1(icm_bra, ncm_bra) = &
            & ho_radial_wf_norm(ncm_bra, lcm_bra, acm_bra, pwd%cm_mom_array(icm_bra)%pp/hc) * (-1.d0)**ncm_bra * &
            & pwd%cm_mom_array(icm_bra)%pp * pwd%cm_mom_array(icm_bra)%ww /hc**2
      end do
    end do
    do ncm_ket = 0, maxval(ket_ho%Ncm)
      do icm_bra = 1, pwd%NMesh_cm
        do icm_ket = 1, pwd%NMesh_cm
          radial_cm2(icm_bra, icm_ket, ncm_ket) = &
              & ho_radial_wf_norm(ncm_ket, lcm_ket, acm_ket, pwd%cm_mom_array(icm_bra)%p(icm_ket)/hc) * (-1.d0)**ncm_ket * &
              & pwd%cm_mom_array(icm_bra)%p(icm_ket) * pwd%cm_mom_array(icm_bra)%w(icm_ket) /hc**2
        end do
      end do
    end do

    do ibra_ho = 1, bra_ho%n_states
      nrel_bra = bra_ho%Nrel(ibra_ho)
      ncm_bra = bra_ho%Ncm(ibra_ho)
      do iket_ho = 1, ket_ho%n_states
        nrel_ket = ket_ho%Nrel(iket_ho)
        ncm_ket = ket_ho%Ncm(iket_ho)

        me = 0.d0
        do ibra_mom = 1, bra_mom%n_states
          icm_bra = bra_mom%Ncm(ibra_mom)
          irel_bra = bra_mom%Nrel(ibra_mom)
          do iket_mom = 1, ket_mom%n_states
            icm_ket = ket_mom%Ncm(iket_mom)
            irel_ket = ket_mom%Nrel(iket_mom)

            me = me + mat_mom%m(ibra_mom, iket_mom) * &
                & radial_rel1(irel_bra, nrel_bra) * radial_rel2(irel_ket, nrel_ket) * &
                & radial_cm1(icm_bra, ncm_bra) * radial_cm2(icm_bra, icm_ket, ncm_ket)
          end do
        end do
        mat_ho%m(ibra_ho,iket_ho) = me

      end do
    end do
    deallocate(radial_rel1)
    deallocate(radial_rel2)
    deallocate(radial_cm1)
    deallocate(radial_cm2)
    !call timer%add(s%str('mom -> ho 2b multipole operator'), omp_get_wtime()-time)
  end subroutine transform_to_ho_ls_channel_old

  subroutine InitPWDFiniteQ(this, ms, Nmax, NZMesh, NQMesh, rankJ, Q, c1, c3, c4, c6, cD, external_field, Qdependent_mesh, verbose)
    use MyLibrary, only: legendre_polynomial, gauss_legendre, hc
    class(PWD_functions), intent(inout) :: this
    type(TwoBodyRelCMSpaceMBasis), intent(in) :: ms
    integer, intent(in) :: Nmax, NZMesh, NQMesh, rankJ
    real(8), intent(in) :: Q, c1, c3, c4, c6, cD
    type(str), intent(in) :: external_field
    logical, intent(in), optional :: Qdependent_mesh
    logical, intent(in), optional :: verbose
    integer :: i, j, k, l1, l2
    real(8), allocatable :: q_tmp(:), w_tmp(:)
    integer :: LrelMax, LcmMax, AlphaMax, BetaMax

    this%NMesh_cm = ms%GetNMeshCM()
    this%NMesh_rel = ms%GetNMeshRel()
    this%NMeshMom = this%NMesh_cm * this%NMesh_rel
    this%c1 = c1
    this%c3 = c3
    this%c4 = c4
    this%c6 = c6
    this%cD = cD
    this%Q = Q
    this%rankJ = rankJ
    this%verbose = .false.
    if(present(verbose)) this%verbose = verbose
    LrelMax = min(Nmax, ms%GetJrelMax() + 1)
    LcmMax = min(ms%GetLcmMax(), Nmax)
    AlphaMax = 4 ! Depending on the case, this might not be enough
    BetaMax = AlphaMax + rankJ + 1
    allocate(this%idx_to_icm(this%NMeshMom))
    allocate(this%idx_to_irel(this%NMeshMom))
    allocate(this%icm_irel_to_idx(this%NMesh_cm, this%NMesh_rel))
    this%icm_irel_to_idx(:,:) = 0
    k = 0
    do i = 1, size(ms%x_cm)
      do j = 1, size(ms%x_rel)
        k = k + 1
        this%idx_to_icm(k) = i
        this%idx_to_irel(k) = j
        this%icm_irel_to_idx(i,j) = k
      end do
    end do

    ! This should be better, but it seems compiler dependent
    !allocate(this%rel_mom_array, source=ms%x_rel*hc)
    !allocate(this%rel_wmom_array, source=ms%w_rel*hc)
    allocate(this%rel_mom_array(size(ms%x_rel)))
    allocate(this%rel_wmom_array(size(ms%w_rel)))
    this%rel_mom_array  = ms%x_rel*hc
    this%rel_wmom_array  = ms%w_rel*hc

    allocate(this%cm_mom_array(size(ms%x_cm)))
    do i = 1, size(ms%x_cm)
      this%cm_mom_array(i)%pp = ms%x_cm(i)*hc
      this%cm_mom_array(i)%ww = ms%w_cm(i)*hc
      allocate(this%cm_mom_array(i)%P(size(ms%x_cm)))
      allocate(this%cm_mom_array(i)%W(size(ms%x_cm)))
      if(present(Qdependent_mesh)) then
        if(Qdependent_mesh) then
          call gauss_legendre(abs(ms%x_cm(i)*hc-Q), (ms%x_cm(i)*hc+Q), q_tmp, w_tmp, size(ms%x_cm))
        else
          q_tmp(:) = ms%x_cm(:)*hc
          w_tmp(:) = ms%w_cm(:)*hc
        end if
      else
        q_tmp(:) = ms%x_cm(:)*hc
        w_tmp(:) = ms%w_cm(:)*hc
      end if
      this%cm_mom_array(i)%P(:) = q_tmp
      this%cm_mom_array(i)%W(:) = w_tmp
      deallocate(q_tmp, w_tmp)
    end do
    !kcmmin = ms%GetPcmMin() * hc
    !kcmmax = ms%GetPcmMax() * hc
    !krelmin = ms%GetPrelMin() * hc
    !krelmax = ms%GetPrelMax() * hc
    kcmmin = minval(this%cm_mom_array(:)%pp)
    kcmmax = maxval(this%cm_mom_array(:)%pp)
    krelmin = minval(this%rel_mom_array)
    krelmax = maxval(this%rel_mom_array)

    allocate(this%Ystore_cm(0:LcmMax,0:LcmMax))
    do l1 = 0, LcmMax
      do l2 = 0, LcmMax
        !call this%Ystore_cm(l2,l1)%init_ystore_lam(l1,l2,Q)
        call this%Ystore_cm(l1,l2)%init_ystore_lam(l1,l2,Q)
      end do
    end do
    if(external_field%val=="axial vector") call store_axial_vector_integrals(Q, LrelMax, AlphaMax, NZMesh, NQMesh, this%NMesh_rel)
    if(external_field%val=="vector") call store_vector_integrals(Q, LrelMax, AlphaMax, NZMesh, NQMesh, this%NMesh_rel)
    this%stored=.true.
    this%time_pwd = 0.d0
  end subroutine InitPWDFiniteQ

  subroutine FinPWDFiniteQ(this)
    class(PWD_functions), intent(inout) :: this
    integer :: l1, l2
    call release_integrals()
    do l1 = 0, ubound(this%Ystore_cm,1)
      do l2 = 0, ubound(this%Ystore_cm,2)
        call this%Ystore_cm(l1,l2)%release_ystore_lam()
      end do
    end do
    deallocate(this%Ystore_cm)
    do l1 = 1, size(this%cm_mom_array)
      deallocate(this%cm_mom_array(l1)%p)
      deallocate(this%cm_mom_array(l1)%w)
    end do
    deallocate(this%cm_mom_array)
    deallocate(this%rel_mom_array, this%rel_wmom_array)
    deallocate(this%idx_to_icm)
    deallocate(this%idx_to_irel)
    deallocate(this%icm_irel_to_idx)
    this%stored=.false.
    this%time_pwd = 0.d0
  end subroutine FinPWDFiniteQ

  function PWD_orbital_mat(fq, lcm_bra, lrel_bra, Lt_bra, lcm_ket, lrel_ket, Lt_ket, alpha, beta, kappa) result(mat)
    use MyLibrary, only: hc, pi, dcg, sjs, snj
    type(angular_integral_lam1), intent(in) :: fq
    integer, intent(in) :: lcm_bra, lrel_bra, Lt_bra
    integer, intent(in) :: lcm_ket, lrel_ket, Lt_ket
    integer, intent(in) :: alpha, beta, kappa
    real(8) :: mat(pwd%NMeshMom,pwd%NMeshMom)
    real(8) :: tmp(pwd%NMeshMom,pwd%NMeshMom)
    integer :: lam1, lam2, Llam
    integer :: ibra, iket, icm_bra, irel_bra, icm_ket, irel_ket
    real(8) :: f_q, mom_fact, fact
    real(8) :: time

    time = omp_get_wtime()
    mat(:,:) = 0.d0
    do Llam = abs(lcm_bra-lcm_ket), lcm_bra+lcm_ket
      do lam1 = max(abs(lrel_bra-lrel_ket), abs(kappa-Llam)), min(lrel_bra+lrel_ket, kappa+Llam)
        do lam2 = max(abs(lam1-alpha), abs(Llam-beta)), min(lam1+alpha, Llam+beta)
          if(mod(lam1+alpha+lam2, 2)==1) cycle
          if(mod(Llam+beta+lam2, 2)==1) cycle
          fact = (-1.d0)**lam2 * &
              & dcg(2*lam1,0,2*alpha,0,2*lam2,0) * dcg(2*Llam,0,2*beta,0,2*lam2,0) * &
              & sjs(2*lam1, 2*alpha, 2*lam2, 2*beta, 2*Llam, 2*kappa) * &
              & snj(2*lcm_bra, 2*lrel_bra, 2*Lt_bra, 2*lcm_ket, 2*lrel_ket, 2*Lt_ket, 2*Llam, 2*lam1, 2*kappa) * &
              & sqrt(dble((2*Llam+1)*(2*lam1+1)))
          if(abs(fact) < 1.d-10) cycle

          do ibra = 1, pwd%NMeshMom
            icm_bra = pwd%idx_to_icm(ibra)
            irel_bra = pwd%idx_to_irel(ibra)
            do iket = 1, pwd%NMeshMom
              icm_ket = pwd%idx_to_icm(iket)
              irel_ket = pwd%idx_to_irel(iket)
              mom_fact = 1.d0 / (PWD%cm_mom_array(icm_bra)%pp * PWD%cm_mom_array(icm_bra)%P(icm_ket) * &
                  & PWD%rel_mom_array(irel_bra) * PWD%rel_mom_array(irel_ket) * PWD%Q)
              f_q = fq%lam1(lam1)%lam2(lam2)%v(irel_bra,irel_ket)
              !tmp(ibra,iket) = f_q * PWD%Ystore_cm(lcm_ket,lcm_bra)%lam(Llam)%y(icm_ket,icm_bra,1) * mom_fact
              tmp(ibra,iket) = f_q * PWD%Ystore_cm(lcm_bra,lcm_ket)%lam(Llam)%y(icm_bra,icm_ket,1) * mom_fact
            end do
          end do

          mat(:,:) = mat(:,:) + tmp(:,:) * fact
        end do
      end do
    end do
    fact = (2.d0*pi)**3 * (-1.d0)**(lcm_bra+lrel_bra) * &
        & sqrt(dble((2*alpha+1)*(2*beta+1)*(2*Lt_bra+1)*(2*Lt_ket+1)*(2*kappa+1)))
    mat(:,:) = mat(:,:) * fact
    pwd%time_pwd = pwd%time_pwd + omp_get_wtime() - time
  end function PWD_orbital_mat

  subroutine init_ystore_lam(this,l1,l2,Q)
    class(Ystore_lam), intent(inout) :: this
    integer, intent(in) :: l1, l2
    real(8), intent(in), optional :: Q
    integer :: lam, lam_min, lam_max
    lam_min = abs(l1-l2)
    lam_max = l1 + l2
    allocate(this%lam(lam_min:lam_max))
    do lam = lam_min, lam_max
      call this%lam(lam)%init_ystore_q(l1,l2,lam,Q)
    end do
  end subroutine init_ystore_lam

  subroutine release_ystore_lam(this)
    class(Ystore_lam), intent(inout) :: this
    integer :: lam
    do lam = lbound(this%lam,1), ubound(this%lam,1)
      call this%lam(lam)%release_ystore_q()
    end do
    deallocate(this%lam)
  end subroutine release_ystore_lam

  subroutine init_ystore_q(this,l1,l2,lam,Q)
    class(Ystore_q), intent(inout) :: this
    integer, intent(in) :: l1, l2, lam
    real(8), intent(in) :: Q
    integer :: i, j

    allocate(this%y(size(pwd%cm_mom_array), size(pwd%cm_mom_array),1))
    do i = 1, size(pwd%cm_mom_array)
      do j = 1, size(pwd%cm_mom_array)
        !this%y(j,i,1) = Yfunc(pwd%cm_mom_array(i)%pp, l1, pwd%cm_mom_array(i)%P(j), l2, Q, lam)
        this%y(i,j,1) = Yfunc(pwd%cm_mom_array(i)%pp, l1, pwd%cm_mom_array(i)%P(j), l2, Q, lam)
      end do
    end do
  end subroutine init_ystore_q

  subroutine release_ystore_q(this)
    class(Ystore_q), intent(inout) :: this
    if(allocated(this%y)) deallocate(this%y)
  end subroutine release_ystore_q

  function Yfunc(p1, l1, p2, l2, q, L) result(res)
    use MyLibrary, only: dcg, spherical_harmonics0
    integer, intent(in) :: l1, l2, L
    real(8), intent(in) :: p1, p2, q
    integer :: m
    real(8) :: cost1, cost2
    real(8) :: res

    cost1 = (p1**2 - p2**2 + q**2) / (2.d0*p1*q)
    cost2 = (p1**2 - p2**2 - q**2) / (2.d0*p2*q)

    res = 0.d0
    do m = -min(l1, l2), min(l1, l2)
      res = res + dcg(2*l1, -2*m, 2*l2, 2*m, 2*L, 0) * &
          & spherical_harmonics0(l1,-m,cost1) * spherical_harmonics0(l2,m,cost2)
    end do
  end function Yfunc

  subroutine store_axial_vector_integrals(Q, LrelMax, AlphaMax, NZMesh, NQMesh, NMesh_rel)
    use MyLibrary, only: legendre_polynomial, gauss_legendre
    real(8), intent(in) :: Q
    integer, intent(in) :: LrelMax, AlphaMax, NZMesh, NQMesh, NMesh_rel
    integer :: K, k1, lam1_min, lam1_max, lam2_min, lam2_max
    integer :: lrel_bra, lrel_ket, lam1, lam2, nz, nq, ibra, iket
    real(8) :: integral, integral1, integral2, prel_bra, prel_ket, fz, fz1, fz2
    real(8), allocatable :: ZMesh(:), WZMesh(:)
    real(8), allocatable, save :: QMesh(:), WQMesh(:)
    !$omp threadprivate(QMesh, WQMesh)
    integer :: iloop
    integer, allocatable :: loops(:,:)
    type(sys) :: s
    real(8) :: time

    time = omp_get_wtime()
    if(pwd%verbose) write(*,'(a,f10.6,a)') 'Estimated memory for angular integrals in 2b multipole: ', &
        & 8.d0 * 28.d0 * dble(LrelMax+1)**2 * dble(2*LrelMax+1) * &
        & dble(2*LrelMax+AlphaMax+1) * dble(NMesh_rel)**2 / (1024.d0)**3, ' GB'
    call gauss_legendre(-1.d0, 1.d0, ZMesh, WZMesh, NZMesh)
    allocate(PWD%fq_unit%lrels(0:LrelMax,0:LrelMax))
    allocate(loops(2,NMesh_rel**2))
    iloop = 0
    do ibra = 1, NMesh_rel
      do iket = 1, NMesh_rel
        iloop = iloop+1
        loops(:,iloop) = [ibra,iket]
      end do
    end do

    do lrel_bra = 0, LrelMax
      do lrel_ket = 0, LrelMax
        lam1_min = abs(lrel_bra - lrel_ket)
        lam1_max = lrel_bra+lrel_ket
        allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
        do lam1 = lam1_min, lam1_max
          lam2_min = 0
          lam2_max = lam1+AlphaMax
          allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))

          do lam2 = lam2_min, lam2_max
            allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))

            !$omp parallel
            !$omp do private(iloop, ibra, iket, prel_bra, prel_ket, integral, nq, fz, nz)
            do iloop = 1, NMesh_rel**2
              ibra = loops(1,iloop)
              iket = loops(2,iloop)
              prel_bra = pwd%rel_mom_array(ibra)
              prel_ket = pwd%rel_mom_array(iket)
              call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, QMesh, WQMesh, NQMesh)

              integral = 0.d0
              do nq = 1, NQMesh
                fz = 0.d0
                if(lam2==0) fz=2.d0
                !do nz = 1, NZMesh
                !  fz = fz + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * 1.d0
                !end do
                integral = integral + fz * WQMesh(nq) * QMesh(nq) * &
                    & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
              end do
              PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral
            end do
            !$omp end do
            !$omp end parallel

          end do
        end do
      end do
    end do

    allocate(PWD%f_ope_q1%args(0:2,0:2,0:2))
    allocate(PWD%f_ope_q2%args(0:2,0:2,0:2))
    do K = 0, 2
      do k1 = 0, 2
        allocate(PWD%f_ope_q1%args(K,k1,0)%lrels(0:LrelMax,0:LrelMax))
        allocate(PWD%f_ope_q2%args(K,k1,0)%lrels(0:LrelMax,0:LrelMax))
        allocate(PWD%f_ope_q1%args(K,k1,1)%lrels(0:LrelMax,0:LrelMax))
        allocate(PWD%f_ope_q2%args(K,k1,1)%lrels(0:LrelMax,0:LrelMax))
        allocate(PWD%f_ope_q1%args(K,k1,2)%lrels(0:LrelMax,0:LrelMax))
        allocate(PWD%f_ope_q2%args(K,k1,2)%lrels(0:LrelMax,0:LrelMax))
        do lrel_bra = 0, LrelMax
          do lrel_ket = 0, LrelMax
            lam1_min = abs(lrel_bra - lrel_ket)
            lam1_max = lrel_bra+lrel_ket
            allocate(PWD%f_ope_q1%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            allocate(PWD%f_ope_q2%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            allocate(PWD%f_ope_q1%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            allocate(PWD%f_ope_q2%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            allocate(PWD%f_ope_q1%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            allocate(PWD%f_ope_q2%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            do lam1 = lam1_min, lam1_max
              lam2_min = 0
              lam2_max = lam1+AlphaMax
              allocate(PWD%f_ope_q1%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
              allocate(PWD%f_ope_q2%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
              allocate(PWD%f_ope_q1%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
              allocate(PWD%f_ope_q2%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
              allocate(PWD%f_ope_q1%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
              allocate(PWD%f_ope_q2%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))

              do lam2 = lam2_min, lam2_max
                allocate(PWD%f_ope_q1%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
                allocate(PWD%f_ope_q2%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
                allocate(PWD%f_ope_q1%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
                allocate(PWD%f_ope_q2%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
                allocate(PWD%f_ope_q1%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
                allocate(PWD%f_ope_q2%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))

                !$omp parallel
                !$omp do private(iloop, ibra, iket, prel_bra, prel_ket, integral1, integral2, nq, fz1, fz2, nz)
                do iloop = 1, NMesh_rel**2
                  ibra = loops(1,iloop)
                  iket = loops(2,iloop)
                  prel_bra = pwd%rel_mom_array(ibra)
                  prel_ket = pwd%rel_mom_array(iket)
                  call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, QMesh, WQMesh, NQMesh)

                  integral1 = 0.d0
                  integral2 = 0.d0
                  do nq = 1, NQMesh
                    fz1 = 0.d0
                    fz2 = 0.d0
                    do nz = 1, NZMesh
                      fz1 = fz1 + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * func_q1(Q, QMesh(nq), ZMesh(nz), [K,k1,0])
                      fz2 = fz2 + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * func_q2(Q, QMesh(nq), ZMesh(nz), [K,k1,0])
                    end do
                    integral1 = integral1 + fz1 * WQMesh(nq) * QMesh(nq) * &
                        & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
                    integral2 = integral2 + fz2 * WQMesh(nq) * QMesh(nq) * &
                        & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
                  end do
                  PWD%f_ope_q1%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral1
                  PWD%f_ope_q2%args(K,k1,0)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral2
                  PWD%f_ope_q1%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral1*Q
                  PWD%f_ope_q2%args(K,k1,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral2*Q
                  PWD%f_ope_q1%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral1*Q*Q
                  PWD%f_ope_q2%args(K,k1,2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral2*Q*Q
                end do
                !$omp end do
                !$omp end parallel

              end do
            end do
          end do
        end do
      end do
    end do

    deallocate(QMesh, WQMesh, ZMesh, WZMesh)
    call timer%add(s%str('Storing integrals in 2b multipole operator'), omp_get_wtime()-time)

  contains
    function func_q1(Qex, q_rel, z, args) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Qex, q_rel, z
      integer, intent(in) :: args(3)
      real(8) :: res, q1
      q1 = sqrt((0.5*Qex)**2 + q_rel**2 + Qex * q_rel * z)
      res = q1**args(1) * q_rel**args(2) * Qex**args(3) / (q1**2 + m_pi**2)
    end function func_q1

    function func_q2(Qex, q_rel, z, args) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Qex, q_rel, z
      integer, intent(in) :: args(3)
      real(8) :: res, q2
      q2 = sqrt((0.5*Qex)**2 + q_rel**2 - Qex * q_rel * z)
      res = q2**args(1) * q_rel**args(2) * Qex**args(3) / (q2**2 + m_pi**2)
    end function func_q2
  end subroutine store_axial_vector_integrals

  subroutine store_vector_integrals(Q, LrelMax, AlphaMax, NZMesh, NQMesh, NMesh_rel)
    use MyLibrary, only: legendre_polynomial, gauss_legendre
    real(8), intent(in) :: Q
    integer, intent(in) :: LrelMax, AlphaMax, NZMesh, NQMesh, NMesh_rel
    integer :: k1, k2, lam1_min, lam1_max, lam2_min, lam2_max
    integer :: lrel_bra, lrel_ket, lam1, lam2, nz, nq, ibra, iket
    real(8) :: integral, integral1, integral2, prel_bra, prel_ket, fz, fz1, fz2
    real(8), allocatable :: ZMesh(:), WZMesh(:)
    real(8), allocatable, save :: QMesh(:), WQMesh(:)
    !$omp threadprivate(QMesh, WQMesh)
    integer :: iloop
    integer, allocatable :: loops(:,:)
    type(sys) :: s
    real(8) :: time

    time = omp_get_wtime()
    if(pwd%verbose) write(*,'(a,f10.6,a)') 'Estimated memory for angular integrals in 2b multipole: ', &
        & 8.d0 * 28.d0 * dble(LrelMax+1)**2 * dble(2*LrelMax+1) * &
        & dble(2*LrelMax+AlphaMax+1) * dble(NMesh_rel)**2 / (1024.d0)**3, ' GB'
    call gauss_legendre(-1.d0, 1.d0, ZMesh, WZMesh, NZMesh)
    allocate(PWD%fq_unit%lrels(0:LrelMax,0:LrelMax))
    allocate(loops(2,NMesh_rel**2))
    iloop = 0
    do ibra = 1, NMesh_rel
      do iket = 1, NMesh_rel
        iloop = iloop+1
        loops(:,iloop) = [ibra,iket]
      end do
    end do

    do lrel_bra = 0, LrelMax
      do lrel_ket = 0, LrelMax
        lam1_min = abs(lrel_bra - lrel_ket)
        lam1_max = lrel_bra+lrel_ket
        allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
        do lam1 = lam1_min, lam1_max
          lam2_min = 0
          lam2_max = lam1+AlphaMax
          allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))

          do lam2 = lam2_min, lam2_max
            allocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))

            !$omp parallel
            !$omp do private(iloop, ibra, iket, prel_bra, prel_ket, integral, nq, fz, nz)
            do iloop = 1, NMesh_rel**2
              ibra = loops(1,iloop)
              iket = loops(2,iloop)
              prel_bra = pwd%rel_mom_array(ibra)
              prel_ket = pwd%rel_mom_array(iket)
              call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, QMesh, WQMesh, NQMesh)

              integral = 0.d0
              do nq = 1, NQMesh
                fz = 0.d0
                if(lam2==0) fz=2.d0
                integral = integral + fz * WQMesh(nq) * QMesh(nq) * &
                    & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
              end do
              PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral
            end do
            !$omp end do
            !$omp end parallel

          end do
        end do
      end do
    end do

    allocate(PWD%f_ope_q1%args(0:1,0:1,1))
    allocate(PWD%f_ope_q2%args(0:1,0:1,1))
    do k1 = 0, 1
      k2 = 1 - k1
      allocate(PWD%f_ope_q1%args(k1,k2,1)%lrels(0:LrelMax,0:LrelMax))
      allocate(PWD%f_ope_q2%args(k1,k2,1)%lrels(0:LrelMax,0:LrelMax))
      do lrel_bra = 0, LrelMax
        do lrel_ket = 0, LrelMax
          lam1_min = abs(lrel_bra - lrel_ket)
          lam1_max = lrel_bra+lrel_ket
          allocate(PWD%f_ope_q1%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
          allocate(PWD%f_ope_q2%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
          do lam1 = lam1_min, lam1_max
            lam2_min = 0
            lam2_max = lam1+AlphaMax
            allocate(PWD%f_ope_q1%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))
            allocate(PWD%f_ope_q2%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))

            do lam2 = lam2_min, lam2_max
              allocate(PWD%f_ope_q1%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))
              allocate(PWD%f_ope_q2%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))

              !$omp parallel
              !$omp do private(iloop, ibra, iket, prel_bra, prel_ket, integral1, integral2, nq, fz1, fz2, nz)
              do iloop = 1, NMesh_rel**2
                ibra = loops(1,iloop)
                iket = loops(2,iloop)
                prel_bra = pwd%rel_mom_array(ibra)
                prel_ket = pwd%rel_mom_array(iket)
                call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, QMesh, WQMesh, NQMesh)

                integral1 = 0.d0
                integral2 = 0.d0
                do nq = 1, NQMesh
                  fz1 = 0.d0
                  fz2 = 0.d0
                  do nz = 1, NZMesh
                    fz1 = fz1 + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * func_q1(Q, QMesh(nq), ZMesh(nz), [k1,k2])
                    fz2 = fz2 + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * func_q2(Q, QMesh(nq), ZMesh(nz), [k1,k2])
                  end do
                  integral1 = integral1 + fz1 * WQMesh(nq) * QMesh(nq) * &
                      & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
                  integral2 = integral2 + fz2 * WQMesh(nq) * QMesh(nq) * &
                      & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
                end do
                PWD%f_ope_q1%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral1
                PWD%f_ope_q2%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral2
              end do
              !$omp end do
              !$omp end parallel

            end do
          end do
        end do
      end do
    end do

    allocate(PWD%f_pif%args(0:2,1:3,1))
    do k1 = 0, 2
      do k2 = 1, 3
        !if(k1+k2 /= 3) cycle
        allocate(PWD%f_pif%args(k1,k2,1)%lrels(0:LrelMax,0:LrelMax))
        do lrel_bra = 0, LrelMax
          do lrel_ket = 0, LrelMax
            lam1_min = abs(lrel_bra - lrel_ket)
            lam1_max = lrel_bra+lrel_ket
            allocate(PWD%f_pif%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1_min:lam1_max))
            do lam1 = lam1_min, lam1_max
              lam2_min = 0
              lam2_max = lam1+AlphaMax
              allocate(PWD%f_pif%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2_min:lam2_max))

              do lam2 = lam2_min, lam2_max
                allocate(PWD%f_pif%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(NMesh_rel,NMesh_rel))

                !$omp parallel
                !$omp do private(iloop, ibra, iket, prel_bra, prel_ket, integral1, nq, fz1, nz)
                do iloop = 1, NMesh_rel**2
                  ibra = loops(1,iloop)
                  iket = loops(2,iloop)
                  prel_bra = pwd%rel_mom_array(ibra)
                  prel_ket = pwd%rel_mom_array(iket)
                  call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, QMesh, WQMesh, NQMesh)

                  integral1 = 0.d0
                  do nq = 1, NQMesh
                    fz1 = 0.d0
                    do nz = 1, NZMesh
                      fz1 = fz1 + WZmesh(nz) * legendre_polynomial(lam2,ZMesh(nz)) * &
                          & func_pif(Q, QMesh(nq), ZMesh(nz), [k1,k2])
                    end do
                    integral1 = integral1 + fz1 * WQMesh(nq) * QMesh(nq) * &
                        & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, QMesh(nq), lam1)
                  end do
                  PWD%f_pif%args(k1,k2,1)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v(ibra,iket) = integral1
                end do
                !$omp end do
                !$omp end parallel

              end do
            end do
          end do
        end do
      end do
    end do

    deallocate(QMesh, WQMesh, ZMesh, WZMesh)
    call timer%add(s%str('Storing integrals in 2b multipole operator'), omp_get_wtime()-time)

  contains
    function func_q1(Qex, q_rel, z, args) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Qex, q_rel, z
      integer, intent(in) :: args(2)
      real(8) :: res, q1
      q1 = sqrt((0.5*Qex)**2 + q_rel**2 + Qex * q_rel * z)
      res = Qex**args(1) * q_rel**args(2) / (q1**2 + m_pi**2)
    end function func_q1

    function func_q2(Qex, q_rel, z, args) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Qex, q_rel, z
      integer, intent(in) :: args(2)
      real(8) :: res, q2
      q2 = sqrt((0.5*Qex)**2 + q_rel**2 - Qex * q_rel * z)
      res = Qex**args(1) * q_rel**args(2) / (q2**2 + m_pi**2)
    end function func_q2

    function func_pif(Qex, q_rel, z, args) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Qex, q_rel, z
      integer, intent(in) :: args(2)
      real(8) :: res, q2, q1
      q1 = sqrt((0.5*Qex)**2 + q_rel**2 + Qex * q_rel * z)
      q2 = sqrt((0.5*Qex)**2 + q_rel**2 - Qex * q_rel * z)
      res = Qex**args(1) * q_rel**args(2) / ((q1**2 + m_pi**2) * (q2**2 + m_pi**2))
    end function func_pif
  end subroutine store_vector_integrals

  subroutine release_integrals()
    ! TODO
    integer :: k1, k2, K
    integer :: lrel_bra, lrel_ket
    integer :: lam1, lam2
    do lrel_bra = lbound(pwd%fq_unit%lrels,1), ubound(pwd%fq_unit%lrels,1)
      do lrel_ket = lbound(pwd%fq_unit%lrels,2), ubound(pwd%fq_unit%lrels,2)
        do lam1 = lbound(pwd%fq_unit%lrels(lrel_bra,lrel_ket)%lam1, 1), &
              & ubound(pwd%fq_unit%lrels(lrel_bra,lrel_ket)%lam1, 1)
          do lam2 = lbound(pwd%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2, 1), &
                & ubound(pwd%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2, 1)
            deallocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v)
          end do
          deallocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2)
        end do
        deallocate(PWD%fq_unit%lrels(lrel_bra,lrel_ket)%lam1)
      end do
    end do
    deallocate(PWD%fq_unit%lrels)

    do K = lbound(pwd%f_ope_q1%args, 1), ubound(pwd%f_ope_q1%args, 1)
      do k1 = lbound(pwd%f_ope_q1%args, 2), ubound(pwd%f_ope_q1%args, 2)
        do k2 = lbound(pwd%f_ope_q1%args, 3), ubound(pwd%f_ope_q1%args, 3)
          if(.not. allocated(PWD%f_ope_q1%args(K,k1,k2)%lrels)) cycle

          do lrel_bra = lbound(pwd%f_ope_q1%args(K,k1,k2)%lrels,1), ubound(pwd%f_ope_q1%args(K,k1,k2)%lrels,1)
            do lrel_ket = lbound(pwd%f_ope_q1%args(K,k1,k2)%lrels,2), ubound(pwd%f_ope_q1%args(K,k1,k2)%lrels,2)

              do lam1 = lbound(pwd%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1,1), &
                    & ubound(pwd%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1,1)
                do lam2 = lbound(pwd%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2,1), &
                      & ubound(pwd%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2,1)

                  deallocate(PWD%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v)
                  deallocate(PWD%f_ope_q2%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v)
                end do
                deallocate(PWD%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2)
                deallocate(PWD%f_ope_q2%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2)
              end do
              deallocate(PWD%f_ope_q1%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1)
              deallocate(PWD%f_ope_q2%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1)
            end do
          end do
          deallocate(PWD%f_ope_q1%args(K,k1,k2)%lrels)
          deallocate(PWD%f_ope_q2%args(K,k1,k2)%lrels)
        end do
      end do
    end do
    deallocate(PWD%f_ope_q1%args)
    deallocate(PWD%f_ope_q2%args)

    if(allocated(pwd%f_pif%args)) then
      do K = lbound(pwd%f_pif%args, 1), ubound(pwd%f_pif%args, 1)
        do k1 = lbound(pwd%f_pif%args, 2), ubound(pwd%f_pif%args, 2)
          do k2 = lbound(pwd%f_pif%args, 3), ubound(pwd%f_pif%args, 3)
            if(.not. allocated(PWD%f_pif%args(K,k1,k2)%lrels)) cycle

            do lrel_bra = lbound(pwd%f_pif%args(K,k1,k2)%lrels,1), ubound(pwd%f_pif%args(K,k1,k2)%lrels,1)
              do lrel_ket = lbound(pwd%f_pif%args(K,k1,k2)%lrels,2), ubound(pwd%f_pif%args(K,k1,k2)%lrels,2)

                do lam1 = lbound(pwd%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1,1), &
                      & ubound(pwd%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1,1)
                  do lam2 = lbound(pwd%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2,1), &
                        & ubound(pwd%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2,1)

                    deallocate(PWD%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2(lam2)%v)
                  end do
                  deallocate(PWD%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1(lam1)%lam2)
                end do
                deallocate(PWD%f_pif%args(K,k1,k2)%lrels(lrel_bra,lrel_ket)%lam1)
              end do
            end do
            deallocate(PWD%f_pif%args(K,k1,k2)%lrels)
          end do
        end do
      end do
      deallocate(PWD%f_pif%args)
    end if
  end subroutine release_integrals

  subroutine init_lscoupled_space(ms)
    type(TwoBodyRelCMSpaceHOBasis), intent(in) :: ms
    integer :: j, p, z, ich
    integer :: spin, Ltot, lrel, Lcm, n_state_ls
    integer :: ncm, nrel, n_state_rad
    ms_ho_ls%NChan = (1+ms%GetJmax()) * 6
    allocate(ms_ho_ls%jpz(ms_ho_ls%NChan))
    allocate(ms_ho_ls%jpz2idx(0:ms%GetJmax(), -1:1, -1:1))
    ms_ho_ls%jpz2idx(:,:,:) = 0
    ich = 0
    do z = -1, 1
      do j = 0, ms%GetJmax()
        do p = 1, -1, -2
          ich = ich + 1
          ms_ho_ls%jpz(ich)%j = j
          ms_ho_ls%jpz(ich)%p = p
          ms_ho_ls%jpz(ich)%z = z
          ms_ho_ls%jpz2idx(j,p,z) = ich
          allocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin2idx(0:ms%GetLcmMax(), 0:ms%GetJRelMax()+1, j-1:j+1, 0:1))
          ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin2idx(:,:,:,:) = 0

          n_state_ls = 0
          do spin = 0, 1
            do Ltot = abs(j-spin), j+spin
              do lrel = 0, ms%GetJRelMax()+1
                do Lcm = abs(Ltot-lrel), min(Ltot+lrel, ms%GetLcmMax())
                  if(abs(z)==1 .and. (-1)**(lrel+spin)==-1) cycle
                  if((-1)**(lrel+Lcm) /= p) cycle
                  if(Lcm+lrel > ms%GetNmax()) cycle
                  n_state_ls = n_state_ls + 1
                  ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin2idx(Lcm, lrel, Ltot, spin) = n_state_ls

                end do
              end do
            end do
          end do
          if(n_state_ls>0) ms_ho_ls%jpz(ich)%zero = .false.
          ms_ho_ls%jpz(ich)%NChan = n_state_ls
          allocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls))
        end do
      end do
    end do

    do ich = 1, ms_ho_ls%NChan
      j = ms_ho_ls%jpz(ich)%j
      p = ms_ho_ls%jpz(ich)%p
      z = ms_ho_ls%jpz(ich)%z

      do spin = 0, 1
        do Ltot = abs(j-spin), j+spin
          do lrel = 0, ms%GetJRelMax()+1
            do Lcm = abs(Ltot-lrel), min(Ltot+lrel, ms%GetLcmMax())
              if(abs(z)==1 .and. (-1)**(lrel+spin)==-1) cycle
              if((-1)**(lrel+Lcm) /= p) cycle
              if(Lcm+lrel > ms%GetNmax()) cycle
              n_state_ls = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin2idx(Lcm, lrel, Ltot, spin)
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Lcm = Lcm
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%lrel = lrel
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ltot = Ltot
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%spin = spin

              allocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm_nrel2idx(&
                  & 0:(ms%GetNmax()-Lcm)/2, 0:(ms%GetNmax()-lrel)/2))
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm_nrel2idx(:,:) = 0
              n_state_rad = 0
              do ncm = 0, (ms%GetNmax()-Lcm)/2
                do nrel = 0, (ms%GetNmax()-lrel)/2
                  if(2*(ncm+nrel)+(Lcm+lrel) > ms%GetNmax()) cycle
                  n_state_rad = n_state_rad + 1
                  ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm_nrel2idx(ncm,nrel) = n_state_rad
                end do
              end do
              ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%n_states = n_state_rad
              allocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm(n_state_rad))
              allocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Nrel(n_state_rad))
            end do
          end do
        end do
      end do
    end do

    do ich = 1, ms_ho_ls%NChan
      j = ms_ho_ls%jpz(ich)%j
      p = ms_ho_ls%jpz(ich)%p
      z = ms_ho_ls%jpz(ich)%z

      do n_state_ls = 1, ms_ho_ls%jpz(ich)%NChan
        Lcm = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Lcm
        lrel = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%lrel
        Ltot = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ltot
        spin = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%spin
        do ncm = 0, (ms%GetNmax()-Lcm)/2
          do nrel = 0, (ms%GetNmax()-lrel)/2
            if(2*(ncm+nrel)+(Lcm+lrel) > ms%GetNmax()) cycle
            n_state_rad = ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm_nrel2idx(ncm,nrel)
            ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm(n_state_rad) = Ncm
            ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%nrel(n_state_rad) = nrel
          end do
        end do
      end do
    end do
  end subroutine init_lscoupled_space

  subroutine release_lscoupled_space()
    integer :: ich, n_state_ls
    do ich = 1, ms_ho_ls%NChan
      do n_state_ls = 1, ms_ho_ls%jpz(ich)%NChan
        deallocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm_nrel2idx)
        deallocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%Ncm)
        deallocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin(n_state_ls)%nrel)
      end do
      deallocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin)
      deallocate(ms_ho_ls%jpz(ich)%Lcm_lrel_Ltot_spin2idx)
    end do
    deallocate(ms_ho_ls%jpz)
    deallocate(ms_ho_ls%jpz2idx)
  end subroutine release_lscoupled_space

  subroutine init_lscoupled_channel_mom(this, lcm, lrel, Ltot, spin)
    type(LSCoupledRadial), intent(inout) :: this
    integer, intent(in) :: lcm, lrel, Ltot, spin

    this%Lcm = Lcm
    this%lrel = lrel
    this%Ltot = Ltot
    this%spin = spin
    this%n_states = size(PWD%idx_to_icm)

    !allocate(this%Ncm_nrel2idx, source=PWD%icm_irel_to_idx)
    !allocate(this%Ncm, source=PWD%idx_to_icm)
    !allocate(this%Nrel, source=PWD%idx_to_irel)
    allocate(this%Ncm_nrel2idx(size(pwd%icm_irel_to_idx,1),size(pwd%icm_irel_to_idx,2)))
    allocate(this%Ncm(size(pwd%idx_to_icm)))
    allocate(this%Nrel(size(pwd%idx_to_irel)))
    this%Ncm_nrel2idx(:,:) = PWD%icm_irel_to_idx(:,:)
    this%Ncm(:) = PWD%idx_to_icm(:)
    this%Nrel(:) = PWD%idx_to_irel(:)
  end subroutine init_lscoupled_channel_mom

  subroutine release_lscoupled_channel_mom(this)
    type(LSCoupledRadial), intent(inout) :: this
    deallocate(this%Ncm_nrel2idx)
    deallocate(this%Ncm)
    deallocate(this%nrel)
  end subroutine release_lscoupled_channel_mom

  subroutine init_op_lscoupled(this)
    use MyLibrary, only: triag
    type(TwoBodyMultiPoleOp), intent(in) :: this
    type(OpLSChannel), pointer :: op_ch
    type(LSCoupledChannel), pointer :: chbra, chket
    type(LSCoupledRadial), pointer :: chbra_rad, chket_rad
    integer :: ichbra, ichket
    integer :: ichbra_ls, ichket_ls
    real(8) :: mem

    mem = 0.d0
    opls%ms => ms_ho_ls
    allocate(opls%OpCh(ms_ho_ls%NChan, ms_ho_ls%NChan))
    do ichbra = 1, ms_ho_ls%NChan
      do ichket = 1, ms_ho_ls%NChan

        chbra => ms_ho_ls%jpz(ichbra)
        chket => ms_ho_ls%jpz(ichket)
        if(chbra%zero) cycle
        if(chket%zero) cycle
        if(triag(this%rankJ, chbra%j, chket%j)) cycle
        if(chbra%p * chket%p * this%rankP==-1) cycle
        if(abs(chbra%z-chket%z) /= this%rankZ) cycle

        op_ch => opls%OpCh(ichbra, ichket)
        op_ch%zero = .false.
        allocate(op_ch%MatCh(chbra%NChan, chket%NChan))
        do ichbra_ls = 1, chbra%NChan
          do ichket_ls = 1, chket%NChan
            chbra_rad => chbra%Lcm_lrel_Ltot_spin(ichbra_ls)
            chket_rad => chket%Lcm_lrel_Ltot_spin(ichket_ls)

            call op_ch%MatCh(ichbra_ls, ichket_ls)%zeros(chbra_rad%n_states, chket_rad%n_states)
            mem = mem + chbra_rad%n_states * chket_rad%n_states * 8.d0 / (1024.d0**3)

          end do
        end do

      end do
    end do
    if(pwd%verbose) write(*,'(a,f10.6,a)') 'Estimated memory for ls coupled op in 2b multipole: ', mem, ' GB'
  end subroutine init_op_lscoupled

  subroutine release_op_lscoupled()
    type(OpLSChannel), pointer :: op_ch
    integer :: ichbra, ichket

    do ichbra = 1, opls%ms%NChan
      do ichket = 1, opls%ms%NChan

        op_ch => opls%OpCh(ichbra, ichket)
        if(allocated(op_ch%MatCh)) deallocate(op_ch%MatCh)
      end do
    end do
    deallocate(opls%OpCh)
  end subroutine release_op_lscoupled

  function L_axial_vector_n2lo_c1(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_iso, tau2_iso, pi, m_pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    real(8) :: fact, fk1, fn, ff, tau1, tau2
    real(8), allocatable :: f_k1(:,:), f_N(:,:)
    type(DMat) :: res
    integer :: lam, k1, k2, N

    tau1 = tau1_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    tau2 = tau2_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    fact = (-2.d0 * PWD%c1 * 1.d-3 * g_A * m_pi**2 / f_pi**2) / (PWD%Q**2 + m_pi**2) / (4.d0 * pi) / (2.d0*pi)**3
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function L_axial_vector_n2lo_c1

  function L_axial_vector_n2lo_c3(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_iso, tau2_iso, pi, g_A, f_pi, m_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact
    real(8) :: tau1, tau2

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    fact = -(PWD%c3 * 1.d-3 * g_A / f_pi**2) / sqrt(4.d0 * pi) / (2.d0*pi)**3
    fact = fact * (1.d0 - PWD%Q**2/(PWD%Q**2 + m_pi**2))
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_axial_vector_n2lo_c3(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble((kappa+1)*(2*kappa+3))) * (-1.d0)
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_c3(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa*(2*kappa-1)))
    end if
    res = res * fact
  end function L_axial_vector_n2lo_c3

  function Tel_axial_vector_n2lo_c3(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: tau1_iso, tau2_iso, pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: ff, fkn, fk1, fx, fact
    real(8), allocatable :: f_kn(:,:), f_k1(:,:), f_x(:,:)
    integer :: lam, K, N, k1, k2, X
    real(8) :: tau1, tau2

    tau1 = tau1_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    tau2 = tau2_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau1*tau2) < 1.d-16) return

    fact = -(PWD%c3 * 1.d-3 * g_A / f_pi**2) / sqrt(4.d0 * pi) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_axial_vector_n2lo_c3(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa*(2*kappa+3)))
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_c3(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble((kappa+1)*(2*kappa-1)))
    end if
    res = res * fact
  end function Tel_axial_vector_n2lo_c3

  function Tmag_axial_vector_n2lo_c3(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: tau1_iso, tau2_iso, pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: ff, fkn, fk1, fx, fact
    real(8), allocatable :: f_kn(:,:), f_k1(:,:), f_x(:,:)
    integer :: lam, K, N, k1, k2, X
    real(8) :: tau1, tau2

    tau1 = tau1_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    tau2 = tau2_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau1*tau2) < 1.d-16) return

    fact = -(PWD%c3 * 1.d-3 * g_A / f_pi**2) / sqrt(4.d0 * pi) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * dble(2*kappa+1)
    res = func_axial_vector_n2lo_c3(kappa, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res * fact
  end function Tmag_axial_vector_n2lo_c3

  function L_axial_vector_n2lo_c4(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi, m_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: ff, fkn, fk1, fx, fact
    real(8), allocatable :: f_kn(:,:), f_k1(:,:), f_x(:,:)
    integer :: lam, K, N, k1, k2, X
    real(8) :: tau_x_tau

    tau_x_tau = tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-16) return

    fact = -(PWD%c4 * 1.d-3 * g_A / f_pi**2) * sqrt(3.d0) / sqrt(8.d0 * pi) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * tau_x_tau * (1.d0 - PWD%Q**2/(PWD%Q**2 + m_pi**2))
    res = func_axial_vector_n2lo_c4(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa+3))) * (-1.d0)
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_c4(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa-1)))
    end if
    res = res * fact
  end function L_axial_vector_n2lo_c4

  function Tel_axial_vector_n2lo_c4(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: ff, fkn, fk1, fx, fact
    real(8), allocatable :: f_kn(:,:), f_k1(:,:), f_x(:,:)
    integer :: lam, K, N, k1, k2, X
    real(8) :: tau_x_tau

    tau_x_tau = tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-16) return

    fact = -(PWD%c4 * 1.d-3 * g_A / f_pi**2) * sqrt(3.d0) / sqrt(8.d0 * pi) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * tau_x_tau
    res = func_axial_vector_n2lo_c4(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa+3)))
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_c4(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa-1)))
    end if
    res = res * fact
  end function Tel_axial_vector_n2lo_c4

  function Tmag_axial_vector_n2lo_c4(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: ff, fkn, fk1, fx, fact
    real(8), allocatable :: f_kn(:,:), f_k1(:,:), f_x(:,:)
    integer :: lam, K, N, k1, k2, X
    real(8) :: tau_x_tau

    tau_x_tau = tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-16) return

    fact = -(PWD%c4 * 1.d-3 * g_A / f_pi**2) * sqrt(3.d0) / sqrt(8.d0 * pi) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * tau_x_tau * dble(2*kappa+1)
    res = func_axial_vector_n2lo_c4(kappa, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket)
    res = res * fact
  end function Tmag_axial_vector_n2lo_c4

  function L_axial_vector_n2lo_cD(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi, m_pi, Lambda_chi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    fact = (PWD%cD / (Lambda_chi * f_pi**2)) / (8.d0 * sqrt(pi)) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * (1.d0 - PWD%Q**2/(PWD%Q**2 + m_pi**2))
    res = func_axial_vector_n2lo_cD(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa+1)) * (-1.d0)
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_cD(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa))
    end if
    res = res * fact
  end function L_axial_vector_n2lo_cD

  function Tel_axial_vector_n2lo_cD(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi, Lambda_chi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    fact = (PWD%cD / (Lambda_chi * f_pi**2)) / (8.d0 * sqrt(pi)) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_axial_vector_n2lo_cD(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa))
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_cD(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) * sqrt(dble(kappa+1))
    end if
    res = res * fact
  end function Tel_axial_vector_n2lo_cD

  function Tmag_axial_vector_n2lo_cD(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi, Lambda_chi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    fact = (PWD%cD / (Lambda_chi * f_pi**2)) / (8.d0 * sqrt(pi)) / (2.d0*pi)**3
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)*(2*kappa+1)))
    res = func_axial_vector_n2lo_cD(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res * fact
  end function Tmag_axial_vector_n2lo_cD

  function Tel_axial_vector_n2lo_c6(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_tau2_tensor_iso, pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact
    real(8) :: tau_x_tau

    tau_x_tau = tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-16) return
    fact =-3.d0 * pwd%c6 * 1.d-3 * g_A / (f_pi**2 * sqrt(8.d0 * pi)) * (-1.d0)**kappa
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_axial_vector_n2lo_c6(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa+3)))
    if(kappa > 0) then
      res = res + func_axial_vector_n2lo_c6(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
          & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa-1)))
    end if
    res = res * fact
  end function Tel_axial_vector_n2lo_c6

  function Tmag_axial_vector_n2lo_c6(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_tau2_tensor_iso, pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact
    real(8) :: tau_x_tau

    tau_x_tau = tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-16) return
    fact = 3.d0 * pwd%c6 * 1.d-3 * g_A / (f_pi**2 * sqrt(8.d0 * pi)) * (-1.d0)**kappa
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * dble(2*kappa+1)
    res = func_axial_vector_n2lo_c6(kappa, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket)
    res = res * fact
  end function Tmag_axial_vector_n2lo_c6

  function func_axial_vector_n2lo_c3(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_iso, tau2_iso, gamma_function
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    integer :: K, N, R, k1, k2, X
    real(8) :: fKN, fR, fk1k2, fX, ts1, ts2
    type(angular_integral_lam1), pointer :: fp

    ts1 = tau1_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket) * tau1_iso(s_bra, s_ket)
    ts2 = tau2_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket) * tau2_iso(s_bra, s_ket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function func_axial_vector_n2lo_c3

  function func_axial_vector_n2lo_c4(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) result(res)
    use MyLibrary, only: dcg, sjs, snj, gamma_function, tau1_tau2_tensor_iso
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket
    type(DMat) :: res
    integer :: K, N, R, Y, k1, k2, Z
    real(8) :: fK, fNR, fY, fk1k2, fZ
    type(angular_integral_lam1), pointer :: fp

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function func_axial_vector_n2lo_c4

  function func_axial_vector_n2lo_c6(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) result(res)
    use MyLibrary, only: dcg, sjs, snj, gamma_function, tau1_tau2_tensor_iso
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket
    type(DMat) :: res
    integer :: K, X, R, Y, k1, k2, Z
    real(8) :: fK, fRX, fY, fk1k2, fZ
    type(angular_integral_lam1), pointer :: fp

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function func_axial_vector_n2lo_c6

  function func_axial_vector_n2lo_cD(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: dcg, sjs, snj, gamma_function, tau1_iso, tau2_iso
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, ts1, ts2
    type(angular_integral_lam1), pointer :: fp

    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    ts1 = tau1_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket) * tau1_iso(s_bra, s_ket)
    ts2 = tau2_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket) * tau2_iso(s_bra, s_ket)
    if(abs(ts1 + ts2) < 1.d-16) return

    fact = snj(2*Lt_bra, 2*s_bra, 2*Jbra, 2*Lt_ket, 2*s_ket, 2*Jket, 2*kappa, 2, 2*lambda)
    fp => pwd%fq_unit%lrels(lrel_bra,lrel_ket)
    res%m(:,:) = pwd_orbital_mat(fp, lcm_bra, lrel_bra, Lt_bra, lcm_ket, lrel_ket, Lt_ket, 0, kappa, kappa) * &
        & (fact * (ts1 + ts2))
  end function func_axial_vector_n2lo_cD

  function L_vector_nlo(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res, res_seagull, res_pif
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    res_seagull = L_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res_pif = L_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res_seagull + res_pif
  end function L_vector_nlo

  function Tel_vector_nlo(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res, res_seagull, res_pif
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    res_seagull = Tel_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res_pif = Tel_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res_seagull + res_pif
  end function Tel_vector_nlo

  function Tmag_vector_nlo(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res, res_seagull, res_pif
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    res_seagull = Tmag_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res_pif = Tmag_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket)
    res = res_seagull + res_pif
  end function Tmag_vector_nlo

  function L_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact = g_A**2 / (4.d0 * f_pi**2 * sqrt(4.d0*pi)) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * (-1.d0)**kappa
    res = func_vector_seagull(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa+3))) * (-1.d0)
    if(kappa > 0) then
      res = res + func_vector_seagull(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa-1)))
    end if
    res = res * fact
  end function L_vector_nlo_seagull

  function L_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact = 3.d0 * g_A**2 / (2.d0 * f_pi**2 * (sqrt(4.d0*pi))) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_vector_pion_in_flight(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa+3))) * (-1.d0)
    if(kappa > 0) then
      res = res + func_vector_pion_in_flight(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa-1)))
    end if
    res = res * fact
  end function L_vector_nlo_pion_in_flight

  function Tel_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact = g_A**2 / (4.d0 * f_pi**2 * sqrt(4.d0*pi)) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * (-1.d0)**kappa
    res = func_vector_seagull(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa*(2*kappa+3)))
    if(kappa > 0) then
      res = res + func_vector_seagull(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa-1)))
    end if
    res = res * fact
  end function Tel_vector_nlo_seagull

  function Tel_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact =-3.d0 * g_A**2 / (2.d0 * f_pi**2 * (sqrt(4.d0*pi))) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_vector_pion_in_flight(kappa+1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble(kappa * (2*kappa+3)))
    if(kappa > 0) then
      res = res + func_vector_pion_in_flight(kappa-1, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) * sqrt(dble((kappa+1)*(2*kappa-1)))
    end if
    res = res * fact
  end function Tel_vector_nlo_pion_in_flight

  function Tmag_vector_nlo_seagull(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact = -g_A**2 / (4.d0 * f_pi**2 * sqrt(4.d0*pi)) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * dble(2*kappa+1) * sqrt(dble((2*Jbra+1)*(2*Jket+1))) * (-1.d0)**kappa
    res = func_vector_seagull(kappa, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket)
    res = res * fact
  end function Tmag_vector_nlo_seagull

  function Tmag_vector_nlo_pion_in_flight(kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, zbra, &
        & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket, zket) result(res)
    use MyLibrary, only: pi, g_A, f_pi
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra, zbra, kappa
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket, zket
    type(DMat) :: res
    real(8) :: fact, tau_x_tau
    tau_x_tau = -tau_x_tau_op(lrel_bra, s_bra, zbra, lrel_ket, s_ket, zket)
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
    if(abs(tau_x_tau) < 1.d-8) return
    fact = 3.d0 * g_A**2 / (2.d0 * f_pi**2 * (sqrt(4.d0*pi))) / (2.d0*pi)**3 * tau_x_tau
    fact = fact * dble(2*kappa+1) * sqrt(dble((2*Jbra+1)*(2*Jket+1)))
    res = func_vector_pion_in_flight(kappa, kappa, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket)
    res = res * fact
  end function Tmag_vector_nlo_pion_in_flight

  function func_vector_seagull(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_tau2_tensor_iso
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket
    type(DMat) :: res
    integer :: k1, k2, K, N, X
    real(8) :: fkN, fx
    type(angular_integral_lam1), pointer :: fp1, fp2
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function func_vector_seagull

  function func_vector_pion_in_flight(kappa, lambda, lcm_bra, lrel_bra, Lt_bra, s_bra, Jbra, &
      & lcm_ket, lrel_ket, Lt_ket, s_ket, Jket) result(res)
    use MyLibrary, only: dcg, sjs, snj, tau1_tau2_tensor_iso
    integer, intent(in) :: kappa, lambda
    integer, intent(in) :: lcm_bra, lrel_bra, s_bra, Lt_bra, Jbra
    integer, intent(in) :: lcm_ket, lrel_ket, s_ket, Lt_ket, Jket
    type(DMat) :: res
    integer :: k1, k2, k3, k4, K, k13, k24, N, X, Y
    real(8) :: fks, fxy
    type(angular_integral_lam1), pointer :: fp
    call res%zeros(pwd%NMeshMom, pwd%NMeshMom)
  end function func_vector_pion_in_flight

  function tau1_op(lbra, sbra, zbra, lket, sket, zket) result(res)
    use MyLibrary, only: tau_1, tau_x, tau_y, tau_z
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    integer :: ibra, iket
    integer, allocatable :: z1bras(:), z2bras(:), z1kets(:), z2kets(:)
    integer :: z1bra, z2bra, z1ket, z2ket
    real(8) :: norm, phbra, phket, res, me
    call set_pn_combination(zbra, z1bras, z2bras)
    call set_pn_combination(zket, z1kets, z2kets)
    res = 0.d0
    norm = 1.d0 / sqrt(dble(size(z1bras) * size(z1kets)) )
    do ibra = 1, size(z1bras)
      phbra = 1.d0
      z1bra = z1bras(ibra); z2bra = z2bras(ibra)
      if(z1bra == 1 .and. z2bra ==-1) phbra = (-1.d0)**(lbra+sbra)
      do iket = 1, size(z1kets)
        phket = 1.d0
        z1ket = z1kets(iket); z2ket = z2kets(iket)
        if(z1ket == 1 .and. z2ket ==-1) phket = (-1.d0)**(lket+sket)
        me = 0.d0
        if(zbra-zket==0) then
          me = tau_z(z1bra,z1ket,phase=tz_phase) * tau_1(z2bra,z2ket)
        else if(zbra-zket==1) then
          me = (tau_x(z1bra,z1ket) + tau_y(z1bra,z1ket,phase=tz_phase) ) * tau_1(z2bra,z2ket) / sqrt(2.d0)
        else if(zbra-zket==-1) then
          me = (tau_x(z1bra,z1ket) - tau_y(z1bra,z1ket,phase=tz_phase) ) * tau_1(z2bra,z2ket) / sqrt(2.d0)
        else
          write(*,*) "Error:", __LINE__, __FILE__
          stop
        end if
        res = res + me * phbra * phket
      end do
    end do
    res = res * norm
    deallocate(z1bras, z2bras, z1kets, z2kets)
  end function tau1_op

  function tau2_op(lbra, sbra, zbra, lket, sket, zket) result(res)
    use MyLibrary, only: tau_1, tau_x, tau_y, tau_z
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    integer :: ibra, iket
    integer, allocatable :: z1bras(:), z2bras(:), z1kets(:), z2kets(:)
    integer :: z1bra, z2bra, z1ket, z2ket
    real(8) :: norm, phbra, phket, res, me
    call set_pn_combination(zbra, z1bras, z2bras)
    call set_pn_combination(zket, z1kets, z2kets)
    res = 0.d0
    norm = 1.d0 / sqrt(dble(size(z1bras) * size(z1kets)) )
    do ibra = 1, size(z1bras)
      phbra = 1.d0
      z1bra = z1bras(ibra); z2bra = z2bras(ibra)
      if(z1bra == 1 .and. z2bra ==-1) phbra = (-1.d0)**(lbra+sbra)
      do iket = 1, size(z1kets)
        phket = 1.d0
        z1ket = z1kets(iket); z2ket = z2kets(iket)
        if(z1ket == 1 .and. z2ket ==-1) phket = (-1.d0)**(lket+sket)
        me = 0.d0
        if(zbra-zket==0) then
          me = tau_z(z2bra,z2ket,phase=tz_phase) * tau_1(z1bra,z1ket)
        else if(zbra-zket==1) then
          me = (tau_x(z2bra,z2ket) + tau_y(z2bra,z2ket,phase=tz_phase) ) * tau_1(z1bra,z1ket) / sqrt(2.d0)
        else if(zbra-zket==-1) then
          me = (tau_x(z2bra,z2ket) - tau_y(z2bra,z2ket,phase=tz_phase) ) * tau_1(z1bra,z1ket) / sqrt(2.d0)
        else
          write(*,*) "Error:", __LINE__, __FILE__
          stop
        end if
        res = res + me * phbra * phket
      end do
    end do
    res = res * norm
    deallocate(z1bras, z2bras, z1kets, z2kets)
  end function tau2_op

  function tau_x_tau_op(lbra, sbra, zbra, lket, sket, zket) result(res)
    use MyLibrary, only: tau_1, tau_x, tau_y, tau_z
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    integer :: ibra, iket
    integer, allocatable :: z1bras(:), z2bras(:), z1kets(:), z2kets(:)
    integer :: z1bra, z2bra, z1ket, z2ket
    real(8) :: norm, phbra, phket, res, me
    call set_pn_combination(zbra, z1bras, z2bras)
    call set_pn_combination(zket, z1kets, z2kets)
    res = 0.d0
    norm = 1.d0 / sqrt(dble(size(z1bras) * size(z1kets)) )
    do ibra = 1, size(z1bras)
      phbra = 1.d0
      z1bra = z1bras(ibra); z2bra = z2bras(ibra)
      if(z1bra == 1 .and. z2bra ==-1) phbra = (-1.d0)**(lbra+sbra)
      do iket = 1, size(z1kets)
        phket = 1.d0
        z1ket = z1kets(iket); z2ket = z2kets(iket)
        if(z1ket == 1 .and. z2ket ==-1) phket = (-1.d0)**(lket+sket)
        me = 0.d0
        if(zbra-zket==0) then
          me = (tau_x(z1bra,z1ket)*tau_y(z2bra,z2ket,phase=tz_phase) - &
              & tau_y(z1bra,z1ket,phase=tz_phase)*tau_x(z2bra,z2ket)) * (-1.d0) ! tau_x_tau = i (tau1 x tau2)
        else if(zbra-zket==1) then
          me = ((tau_x(z1bra,z1ket) - tau_y(z1bra,z1ket,phase=tz_phase)) * tau_z(z2bra,z2ket,phase=tz_phase) - &
              & tau_z(z1bra,z1ket,phase=tz_phase) * (tau_x(z2bra,z2ket) - tau_y(z2bra,z2ket,phase=tz_phase))) / sqrt(2.d0)
        else if(zbra-zket==-1) then
          me = -((tau_x(z1bra,z1ket) + tau_y(z1bra,z1ket,phase=tz_phase)) * tau_z(z2bra,z2ket,phase=tz_phase) - &
              & tau_z(z1bra,z1ket,phase=tz_phase) * (tau_x(z2bra,z2ket) + tau_y(z2bra,z2ket,phase=tz_phase))) / sqrt(2.d0)
        else
          write(*,*) "Error:", __LINE__, __FILE__
          stop
        end if
        res = res + me * phbra * phket
      end do
    end do
    res = res * norm
    deallocate(z1bras, z2bras, z1kets, z2kets)
  end function tau_x_tau_op

  subroutine set_pn_combination(Tz, res1, res2)
    integer, intent(in) :: Tz
    integer, allocatable, intent(out) :: res1(:), res2(:)
    if(Tz==-1) then
      allocate(res1(1), res2(1))
      res1 = [-1]; res2 = [-1]
    elseif(Tz==1) then
      allocate(res1(1), res2(1))
      res1 = [1]; res2 = [1]
    elseif(Tz==0) then
      allocate(res1(2), res2(2))
      res1 = [-1,1]; res2 = [1,-1]
    else
      write(*,*) __LINE__, __FILE__
      stop
    end if
  end subroutine set_pn_combination

  function regulator_function(p) result(f)
    real(8), intent(in) :: p
    real(8) :: f
    f = exp(-(p**2 / freg%lambda**2)**freg%reg_pow)
  end function regulator_function


  !
  ! test for the integral
  !
  subroutine test_integrals()
    use MyLibrary, only: gauss_legendre, spherical_harmonics, pi, dcg, triag
    integer :: kappa, mu
    real(8) :: Q, p_rel_bra, p_cm_bra, p_rel_ket, p_cm_ket
    real(8), allocatable :: zmesh(:), zwmesh(:), phi_mesh(:), phi_wmesh(:)
    integer :: NMesh, loop, idx
    integer :: Lt_bra, Lt_ket
    integer :: mbra, mket
    integer :: l_rel_bra, l_cm_bra, l_rel_ket, l_cm_ket
    integer :: m_rel_bra, m_cm_bra, m_rel_ket, m_cm_ket
    integer :: alpha, beta
    integer, allocatable :: l_loops(:,:)
    complex(8) :: rr, integral
    real(8) :: r, cgs
    integer :: k1, k2, k3, k4, k13, k24
    integer :: m1, m2, m3, m4
    NMesh=12
    p_rel_bra = 500.d0
    p_cm_bra = 400.d0
    p_rel_ket = 200.d0
    p_cm_ket = 300.d0
    alpha = 1
    beta = 1
    k1 = 1; k2 = 1; k3 = 1; k4 = 1
    k13 = 1; k24 = 2
    do loop = 1, 2
      idx = 0
      do l_rel_bra = 0, 2
        do l_rel_ket = 0, 2
          do l_cm_bra = 0, 2
            do l_cm_ket = 0, 2
              do Lt_bra = 0, 2
                do Lt_ket = 0, 2
                  do kappa = 0, 2, 2
                    if(triag(l_rel_bra, l_cm_bra, Lt_bra)) cycle
                    if(triag(l_rel_ket, l_cm_ket, Lt_ket)) cycle
                    if(triag(Lt_bra, Lt_ket, kappa)) cycle
                    !if(triag(alpha, kappa, kappa)) cycle
                    idx = idx + 1
                    if(loop==1) cycle
                    l_loops(:,idx) = [l_rel_bra, l_rel_ket, l_cm_bra, l_cm_ket, Lt_bra, Lt_ket, kappa]
                    write(*,*) l_loops(:,idx)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      if(loop==1) allocate(l_loops(7,idx))
    end do
    write(*,*) size(l_loops,2)

    Q = 200.d0
    call gauss_legendre(-1.d0, 1.d0, zmesh, zwmesh, NMesh)
    call gauss_legendre(0.d0, 2.d0*pi, phi_mesh, phi_wmesh, NMesh)

    mbra = 0
    do loop = 1, size(l_loops,2)
      l_rel_bra = l_loops(1,loop)
      l_rel_ket = l_loops(2,loop)
      l_cm_bra = l_loops(3,loop)
      l_cm_ket = l_loops(4,loop)
      Lt_bra = l_loops(5,loop)
      Lt_ket = l_loops(6,loop)
      kappa = l_loops(7,loop)
      r = PWDintegral(p_cm_bra, p_rel_bra, l_cm_bra, l_rel_bra, Lt_bra, &
          & p_cm_ket, p_rel_ket, l_cm_ket, l_rel_ket, Lt_ket, alpha, beta, kappa) * Q * p_cm_bra * p_cm_ket
      !r = PWDintegral_gen(p_cm_bra, p_rel_bra, l_cm_bra, l_rel_bra, Lt_bra, &
      !    & p_cm_ket, p_rel_ket, l_cm_ket, l_rel_ket, Lt_ket, k1, k2, k3, k4, k13, k24, kappa) * Q * p_cm_bra * p_cm_ket
      if(abs(r)<1.d-8) cycle
      rr = 0.d0
      mbra = Lt_bra
      mket = Lt_ket
      mu = mbra-mket
      do m_rel_bra = -l_rel_bra, l_rel_bra
        do m_rel_ket = -l_rel_ket, l_rel_ket
          do m_cm_bra = -l_cm_bra, l_cm_bra
            do m_cm_ket = -l_cm_ket, l_cm_ket
              if(m_rel_bra+m_cm_bra /= mbra) cycle
              if(m_rel_ket+m_cm_ket /= mket) cycle

              cgs = sqrt(dble(2*Lt_bra+1)) / dcg(2*Lt_ket, 2*mket, 2*kappa, 2*mu, 2*Lt_bra, 2*mbra) * &
                  & dcg(2*l_cm_bra, 2*m_cm_bra, 2*l_rel_bra, 2*m_rel_bra, 2*Lt_bra, 2*mbra) * &
                  & dcg(2*l_cm_ket, 2*m_cm_ket, 2*l_rel_ket, 2*m_rel_ket, 2*Lt_ket, 2*mket)
              if(abs(cgs)<1.d-8) cycle
              integral = explicit_integral(l_cm_bra, m_cm_bra, l_rel_bra, m_rel_bra, &
                  & l_cm_ket, m_cm_ket, l_rel_ket, m_rel_ket, kappa, mu)
              !write(*,*) integral
              rr = rr + cgs * integral

            end do
          end do
        end do
      end do
      write(*,*) loop, l_cm_bra, l_rel_bra, Lt_bra, l_cm_ket, l_rel_ket, Lt_ket, kappa, real(rr), r
    end do

  contains
    function get_r_theta_phi(x, y, z) result(rtp)
      use MyLibrary, only: pi
      real(8), intent(in) :: x, y, z
      real(8) :: rtp(3)
      rtp(1) = sqrt(x**2 + y**2 + z**2)
      rtp(2) = acos(z/rtp(1))
      rtp(3) = acos(x/sqrt(x**2+y**2))
      if(y<0.d0) rtp(3) = - rtp(3)
    end function get_r_theta_phi

    function explicit_integral(l_cm_bra, m_cm_bra, l_rel_bra, m_rel_bra, l_cm_ket, m_cm_ket, l_rel_ket, m_rel_ket, kappa, mu) &
          & result(res)
      integer, intent(in) :: l_cm_bra, m_cm_bra, l_rel_bra, m_rel_bra, l_cm_ket, m_cm_ket, l_rel_ket, m_rel_ket, kappa, mu
      integer :: ibra_cm_theta, ibra_cm_phi
      integer :: i_cm_phi
      integer :: ibra_rel_theta, ibra_rel_phi
      integer :: iket_rel_theta, iket_rel_phi
      real(8) :: cm_theta_bra, cm_phi_bra
      real(8) :: cm_theta_ket, cm_phi_ket
      real(8) :: cm_theta, cm_phi
      real(8) :: rel_theta_bra, rel_phi_bra
      real(8) :: rel_theta_ket, rel_phi_ket
      real(8) :: c
      complex(8) :: res, f1
      real(8) :: rot1(3,3), rot2(3,3), p_cm_bra_vec(3), p_cm_ket_vec(3), rtp_Qcm(3), QcmVec(3), QrelVec(3)
      real(8) :: p_rel_bra_vec(3), p_rel_ket_vec(3), rtp_qrel(3)
      integer :: m_alpha, m_beta
      res = 0.d0
      c = (p_cm_bra**2 + p_cm_ket**2 - Q**2)/(2.d0*p_cm_bra*p_cm_ket)
      cm_theta = acos(c)
      do i_cm_phi = 1, NMesh
        cm_phi = phi_mesh(i_cm_phi)

        do ibra_cm_theta = 1, NMesh
          do ibra_cm_phi = 1, NMesh
            cm_theta_bra = acos(zmesh(ibra_cm_theta))
            cm_phi_bra = phi_mesh(ibra_cm_phi)

            rot1(1,:) = [cos(-cm_theta_bra),0.d0, -sin(-cm_theta_bra)]
            rot1(2,:) = [0.d0, 1.d0, 0.d0]
            rot1(3,:) = [sin(-cm_theta_bra),0.d0, cos(-cm_theta_bra)]

            rot2(1,:) = [cos(-cm_phi_bra), sin(-cm_phi_bra), 0.d0]
            rot2(2,:) = [-sin(-cm_phi_bra), cos(-cm_phi_bra), 0.d0]
            rot2(3,:) = [0.d0, 0.d0, 1.d0]


            p_cm_bra_vec = [0.d0, 0.d0, p_cm_bra]
            p_cm_ket_vec = [p_cm_ket*sin(cm_theta)*cos(cm_phi), p_cm_ket*sin(cm_theta)*sin(cm_phi), p_cm_ket*cos(cm_theta)]

            p_cm_bra_vec = matmul(rot2, matmul(rot1, p_cm_bra_vec))
            p_cm_ket_vec = matmul(rot2, matmul(rot1, p_cm_ket_vec))
            QcmVec = p_cm_bra_vec - p_cm_ket_vec
            rtp_Qcm = get_r_theta_phi(p_cm_ket_vec(1), p_cm_ket_vec(2), p_cm_ket_vec(3))
            cm_theta_ket = rtp_Qcm(2)
            cm_phi_ket = rtp_Qcm(3)
            rtp_Qcm = get_r_theta_phi(QcmVec(1), QcmVec(2), QcmVec(3))

            do ibra_rel_theta = 1, NMesh
              do ibra_rel_phi = 1, NMesh
                rel_theta_bra = acos(zmesh(ibra_rel_theta))
                rel_phi_bra = phi_mesh(ibra_rel_phi)
                do iket_rel_theta = 1, NMesh
                  do iket_rel_phi = 1, NMesh
                    rel_theta_ket = acos(zmesh(iket_rel_theta))
                    rel_phi_ket = phi_mesh(iket_rel_phi)

                    p_rel_bra_vec = [p_rel_bra*sin(rel_theta_bra)*cos(rel_phi_bra), &
                        & p_rel_bra*sin(rel_theta_bra)*sin(rel_phi_bra), p_rel_bra*cos(rel_theta_bra)]
                    p_rel_ket_vec = [p_rel_ket*sin(rel_theta_ket)*cos(rel_phi_ket), &
                        & p_rel_ket*sin(rel_theta_ket)*sin(rel_phi_ket), p_rel_ket*cos(rel_theta_ket)]
                    QrelVec = p_rel_bra_vec - p_rel_ket_vec
                    rtp_qrel = get_r_theta_phi(QrelVec(1), QrelVec(2), QrelVec(3))

                    f1 = 0.d0
                    ! special case
                    do m_alpha = -alpha, alpha
                      m_beta = mu-m_alpha
                      if(abs(m_beta)>beta) cycle
                      f1 = f1 + spherical_harmonics(alpha, m_alpha, cos(rtp_qrel(2)), rtp_qrel(3)) * &
                          & spherical_harmonics(beta, m_beta, cos(rtp_Qcm(2)), rtp_Qcm(3)) * &
                          & dcg(2*alpha, 2*m_alpha, 2*beta, 2*m_beta, 2*kappa, 2*mu)
                    end do

                    ! general case
                    !do m1 = -k1, k1
                    !  do m2 = -k2, k2
                    !    do m3 = -k3, k3
                    !      do m4 = -k4, k4
                    !        if(abs(m1+m3)>k13) cycle
                    !        if(abs(m2+m4)>k24) cycle
                    !        if(m1+m2+m3+m4/=mu) cycle
                    !        if(abs(m1+m2+m3+m4)>kappa) cycle
                    !        f1 = f1 + &
                    !            & spherical_harmonics(k1, m1, cos(cm_theta_bra), cm_phi_bra) * &
                    !            & spherical_harmonics(k2, m2, cos(rel_theta_bra), rel_phi_bra) * &
                    !            & spherical_harmonics(k3, m3, cos(cm_theta_ket), cm_phi_ket) * &
                    !            & spherical_harmonics(k4, m4, cos(rel_theta_ket), rel_phi_ket) * &
                    !            & dcg(2*k1, 2*m1, 2*k3, 2*m3, 2*k13, 2*(m1+m3)) * &
                    !            & dcg(2*k2, 2*m2, 2*k4, 2*m4, 2*k24, 2*(m2+m4)) * &
                    !            & dcg(2*k13, 2*(m1+m3), 2*k24, 2*(m2+m4), 2*kappa, 2*mu)
                    !      end do
                    !    end do
                    !  end do
                    !end do
                    if(abs(f1)<1.d-8) cycle
                    res = res + f1 * func2(Q, rtp_qrel(1), dot_product(QcmVec, QrelVec)/(Q*rtp_qrel(1))) * &
                        & conjg(spherical_harmonics(l_cm_bra, m_cm_bra, cos(cm_theta_bra), cm_phi_bra)) * &
                        & conjg(spherical_harmonics(l_rel_bra, m_rel_bra, cos(rel_theta_bra), rel_phi_bra)) * &
                        & spherical_harmonics(l_cm_ket, m_cm_ket, cos(cm_theta_ket), cm_phi_ket) * &
                        & spherical_harmonics(l_rel_ket, m_rel_ket, cos(rel_theta_ket), rel_phi_ket) * &
                        & phi_wmesh(i_cm_phi) * zwmesh(ibra_cm_theta) * phi_wmesh(ibra_cm_phi) * &
                        & zwmesh(ibra_rel_theta) * phi_wmesh(ibra_rel_phi) * zwmesh(iket_rel_theta) * phi_wmesh(iket_rel_phi)
                  end do
                end do
              end do
            end do

          end do
        end do

      end do
    end function explicit_integral

    !function func1(Q, qrel, c) result(res)
    !  real(8), intent(in) :: Q, qrel, c
    !  real(8) :: res
    !  res = 1.d0
    !end function func1

    function func2(Q, qrel, c) result(res)
      use MyLibrary, only: m_pi
      real(8), intent(in) :: Q, qrel, c
      real(8) :: q1
      real(8) :: res
      q1 = sqrt(Q**2+qrel**2+2.d0*Q*qrel*c)
      res = q1**2 / (q1**2 + m_pi**2)
    end function func2

    function PWDintegral(pcm_bra, prel_bra, lcm_bra, lrel_bra, Lt_bra, &
          & pcm_ket, prel_ket, lcm_ket, lrel_ket, Lt_ket, alpha, beta, kappa) result(res)
      use MyLibrary, only: pi, dcg, sjs, snj, legendre_polynomial, hc
      real(8), intent(in) :: pcm_bra, prel_bra, pcm_ket, prel_ket
      integer, intent(in) :: lcm_bra, lrel_bra, Lt_bra
      integer, intent(in) :: lcm_ket, lrel_ket, Lt_ket
      integer, intent(in) :: alpha, beta, kappa
      real(8) :: res
      real(8), allocatable :: c_mesh(:), c_wmesh(:), q_mesh(:), q_wmesh(:)
      integer :: lam1, lam2, Llam, nq, nz
      real(8) :: f_q, f_z, mom_fact

      call gauss_legendre(-1.d0, 1.d0, c_mesh, c_wmesh, NMesh)
      call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, q_mesh, q_wmesh, NMesh)
      res = 0.d0
      do Llam = abs(lcm_bra-lcm_ket), lcm_bra+lcm_ket
        do lam1 = max(abs(lrel_bra-lrel_ket), abs(kappa-Llam)), min(lrel_bra+lrel_ket, kappa+Llam)
          do lam2 = max(abs(lam1-alpha), abs(Llam-beta)), min(lam1+alpha, Llam+beta)
            if(mod(lam1+alpha+lam2, 2)==1) cycle
            if(mod(Llam+beta+lam2, 2)==1) cycle

            f_q = 0.d0
            do nq = 1, NMesh
              f_z = 0.d0
              do nz = 1, NMesh
                f_z = f_z + c_wmesh(nz) * legendre_polynomial(lam2,c_mesh(nz)) * func2(Q,q_mesh(nq),c_mesh(nz))
              end do
              f_q = f_q + f_z * q_wmesh(nq) * q_mesh(nq) * &
                  & Yfunc(prel_bra, lrel_bra, prel_ket, lrel_ket, q_mesh(nq), lam1)
              if(freg%reg_local) f_q = f_q * regulator_function(q_mesh(nq)*hc)
            end do
            res = res + f_q * (-1.d0)**lam2 * &
                & Yfunc(pcm_bra, lcm_bra, pcm_ket, lcm_ket, Q, Llam) * &
                & dcg(2*lam1, 0, 2*alpha, 0, 2*lam2, 0) * dcg(2*Llam, 0, 2*beta, 0, 2*lam2, 0) * &
                & sjs(2*lam1, 2*alpha, 2*lam2, 2*beta, 2*Llam, 2*kappa) * &
                & snj(2*lcm_bra, 2*lrel_bra, 2*Lt_bra, 2*lcm_ket, 2*lrel_ket, 2*Lt_ket, 2*Llam, 2*lam1, 2*kappa) * &
                & sqrt(dble((2*Llam+1)*(2*lam1+1)))
          end do
        end do
      end do
      res = res * (2.d0*pi)**3 * (-1.d0)**(lcm_bra+lrel_bra) * &
          & sqrt(dble((2*alpha+1)*(2*beta+1)*(2*Lt_bra+1)*(2*Lt_ket+1)*(2*kappa+1)))
      mom_fact = pcm_bra * pcm_ket * prel_bra * prel_ket * Q
      res = res / mom_fact
    end function PWDintegral

    function PWDintegral_gen(pcm_bra, prel_bra, lcm_bra, lrel_bra, Lt_bra, &
          & pcm_ket, prel_ket, lcm_ket, lrel_ket, Lt_ket, k1, k2, k3, k4, k13, k24, kappa) result(res)
      use MyLibrary, only: pi, dcg, sjs, snj, legendre_polynomial, hc
      real(8), intent(in) :: pcm_bra, prel_bra, pcm_ket, prel_ket
      integer, intent(in) :: lcm_bra, lrel_bra, Lt_bra
      integer, intent(in) :: lcm_ket, lrel_ket, Lt_ket
      integer, intent(in) :: k1, k2, k3, k4, k13, k24, kappa
      real(8) :: res
      real(8), allocatable :: c_mesh(:), c_wmesh(:), q_mesh(:), q_wmesh(:)
      integer :: lam, X, Llam, nq, nz
      integer :: l1, l2, l3, l4
      real(8) :: f_q, f_z, mom_fact

      call gauss_legendre(-1.d0, 1.d0, c_mesh, c_wmesh, NMesh)
      call gauss_legendre(abs(prel_bra-prel_ket), prel_bra+prel_ket, q_mesh, q_wmesh, NMesh)
      res = 0.d0
      do l1 = abs(lcm_bra-k1), lcm_bra+k1
        do l2 = abs(lrel_bra-k2), lrel_bra+k2
          do l3 = abs(lcm_ket-k3), lcm_ket+k3
            do l4 = abs(lrel_ket-k4), lrel_ket+k4

              do Llam = abs(lcm_bra-lcm_ket), lcm_bra+lcm_ket
                do lam = max(abs(lrel_bra-lrel_ket), abs(kappa-Llam)), min(lrel_bra+lrel_ket, kappa+Llam)
                  do X = max(abs(l1-l3), abs(l2-l4), abs(Llam-k13), abs(lam-k24)), min(l1+l3, l2+l4, Llam+k13, lam+k24)

                    f_q = 0.d0
                    do nq = 1, NMesh
                      f_z = 0.d0
                      do nz = 1, NMesh
                        f_z = f_z + c_wmesh(nz) * legendre_polynomial(X,c_mesh(nz)) * func2(Q,q_mesh(nq),c_mesh(nz))
                      end do
                      f_q = f_q + f_z * q_wmesh(nq) * q_mesh(nq) * &
                          & Yfunc(prel_bra, l2, prel_ket, l4, q_mesh(nq), X)
                      if(freg%reg_local) f_q = f_q * regulator_function(q_mesh(nq)*hc)
                    end do
                    res = res + f_q * (-1.d0)**lam * &
                        & Yfunc(pcm_bra, l1, pcm_ket, l3, Q, X) * &
                        & dble((2*Llam+1)*(2*lam+1)) * &
                        & snj(2*lcm_bra, 2*lrel_bra, 2*Lt_bra, 2*lcm_ket, 2*lrel_ket, 2*Lt_ket, 2*Llam, 2*lam, 2*kappa) * &
                        & snj(2*lcm_bra, 2*k1, 2*l1, 2*lcm_ket, 2*k3, 2*l3, 2*Llam, 2*k13, 2*X) * &
                        & snj(2*lrel_bra, 2*k2, 2*l2, 2*lrel_ket, 2*k4, 2*l4, 2*lam, 2*k24, 2*X) * &
                        & sjs(2*Llam, 2*k13, 2*X, 2*k24, 2*lam, 2*kappa) * &
                        & dcg(2*lcm_bra, 0, 2*k1, 0, 2*l1, 0) * dcg(2*lrel_bra, 0, 2*k2, 0, 2*l2, 0) * &
                        & dcg(2*lcm_ket, 0, 2*k3, 0, 2*l3, 0) * dcg(2*lrel_ket, 0, 2*k4, 0, 2*l4, 0)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      res = res * 2.d0*pi**2 * (-1.d0)**(lcm_bra+lrel_bra+k13) * &
          & sqrt(dble((2*lcm_bra+1)*(2*lrel_bra+1)*(2*Lt_bra+1)*(2*lcm_ket+1)*(2*lrel_ket+1)*(2*Lt_ket+1)*&
          & (2*k1+1)*(2*k2+1)*(2*k3+1)*(2*k4+1)*(2*k13+1)*(2*k24+1)*(2*kappa+1)))
      mom_fact = pcm_bra * pcm_ket * prel_bra * prel_ket * Q
      res = res / mom_fact
    end function PWDintegral_gen
  end subroutine test_integrals

end module TwoBodyMultiPoleOperator
