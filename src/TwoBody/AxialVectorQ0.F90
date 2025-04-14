!
! A clear test is needed, but it looks like I can reproduce the figure 2 in Eur. Phys. J. A 54, 76 (2018).
!
! For Tz = +/- case with pn formalism:
! Following a spherical tensor definition, we should compute A(+/-) as (-/+) [Ax +/- i Ay]/sqrt(2).
! However, this definition is conflict with GT operator. Also, it is not Hermitian, for example, {-[tx + ity] /sqrt(2)}^t != [tx - ity] / sqrt(2)
! To avoid this, A(+/-) is defined as
!                        A(+/-) = [Ax +/- iAy] / sqrt(2)
! Then, the GT operator correction is given as GT(+/-) = - A(+/-) sqrt(2) / gA
! (Note that LO 1B GT operator is GT(+/-) = sigma [tx +/- ity]/2 = - A(+/-) sqrt(2) / gA with A(1B,+/-) = -gA sigma t(+/-) / sqrt(2) )
!
module AxialVectorQ0
  !$ use omp_lib
  use StoreCouplings
  use PartialWaveDecomposition, only: PWD, &
      & identity_helicity_rep, sigma_dot_sigma, sigma1q_sigma2q, &
      & sigma_plus_q, sigma_plus_QQ, sigma_minus_q, sigma_minus_QQ, sigma_cross_q, sigma_cross_QQ
  implicit none
  public :: AxialCurrentQ0
  private

  type, extends(PWD) :: AxialCurrentQ0
    integer :: ch_order = 0
    character(:), allocatable :: OpName
    character(:), allocatable :: NNint
    real(8) :: c1, c3, c4, cD
  contains
    procedure :: InitAxialCurrentQ0
    procedure :: FinAxialCurrentQ0
    procedure :: show_options
    procedure :: calc_matrix_element
    generic :: init => InitAxialCurrentQ0
    generic :: fin => FinAxialCurrentQ0
  end type AxialCurrentQ0

  logical, private :: SFRegulator = .false.
  real(8), private :: LambdaSFR
  type :: MomFunctions
    real(8), allocatable :: op(:,:,:,:)
    real(8), allocatable :: QMesh(:)
    real(8), allocatable :: pi_prop(:)
    real(8), allocatable :: loop_L(:)
    real(8), allocatable :: loop_H(:)
    real(8), allocatable :: loop_A(:)
    real(8), allocatable :: s_var(:)
  end type MomFunctions

  type(CGsStore) :: cgs
  real(8) :: y1_fact, y1_fact_inv
  logical, parameter :: relativistic_correction = .true.
  !logical, parameter :: relativistic_correction = .false.
contains

  subroutine show_options(this)
    class(AxialCurrentQ0), intent(in) :: this
    write(*,"(a)") "# Calculation options for axial current operator @ Q=0"
    write(*,"(2a)") "# Operator: ", trim(this%OpName)
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a)") "# c3= ", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    write(*,"(a,f8.4)") "# cD= ", this%cD
    if(relativistic_correction) write(*,"(a,f8.4)") "# Relativistic correction: on"
    if(.not. relativistic_correction) write(*,"(a,f8.4)") "# Relativistic correction: off"
    call this%show_pwd_options()
  end subroutine show_options

  subroutine FinAxialCurrentQ0(this)
    class(AxialCurrentQ0), intent(inout) :: this
    this%ch_order = 0
    this%OpName = "None"
    this%c1 = 0.d0
    this%c3 = 0.d0
    this%c4 = 0.d0
    this%cD = 0.d0
    call cgs%fin()
    call this%fin_pwd()
  end subroutine FinAxialCurrentQ0

  subroutine InitAxialCurrentQ0(this, OpName, LECs, ch_order, Jmax, &
        & Regulator, RegulatorPower, RegulatorLambda, LamSFR, mode)
    use MyLibrary, only: pi
    class(AxialCurrentQ0), intent(inout) :: this
    character(*), intent(in) :: OpName
    real(8), intent(in) :: LECs(:)
    real(8), intent(in), optional :: RegulatorLambda, LamSFR
    integer, intent(in), optional :: ch_order, Jmax, RegulatorPower
    character(*), intent(in), optional :: Regulator
    character(*), intent(in), optional :: mode
    this%OpName = OpName
    call this%init_pwd(Jmax=Jmax, Rank=1, Regulator=Regulator, &
        & RegulatorPower=RegulatorPower, RegulatorLambda=RegulatorLambda)
    y1_fact = sqrt(3.d0 / (4.d0*pi))
    y1_fact_inv = 1.d0 / y1_fact
    call cgs%init(0, 2, .false., 0, 2, .false.)
    if( present(ch_order) ) this%ch_order = ch_order
    if( present(LamSFR)) then
      SFRegulator = .true.
      LambdaSFR = LamSFR
    end if

    this%c1 = LECs(1)
    this%c3 = LECs(2)
    this%c4 = LECs(3)
    this%cD = LECs(4)
  end subroutine InitAxialCurrentQ0

  function calc_matrix_element(this, pbra, lbra, sbra, jbra, zbra, pket, lket, sket, jket, zket, pn_formalism) result(r)
    use MyLibrary, only: pi, tau_1, tau_m, tau1_cross_tau2
    class(AxialCurrentQ0), intent(in) :: this
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: lbra, sbra, jbra, zbra, lket, sket, jket, zket
    logical, intent(in) :: pn_formalism
    real(8) :: r
    type(MomFunctions) :: fq

    call set_mom_functions(this, fq, pbra, pket)
    r = 0.d0
    if(pn_formalism) then
      fq%op(:,:,:,:) = 0.d0
      call set_helicity_rep_pn(this, fq, pbra, pket, lbra, sbra, zbra, lket, sket, zket)
      r = r + this%do_pwd(fq%op, pbra, lbra, sbra, jbra, pket, lket, sket, jket)
    else
      fq%op(:,:,:,:) = 0.d0
      call set_helicity_rep_isospin(this, fq, pbra, pket, zbra, zket)
      r = this%do_pwd(fq%op, pbra, lbra, sbra, jbra, pket, lket, sket, jket)
    end if

    call release_mom_functions(fq)
  end function calc_matrix_element

  subroutine set_helicity_rep_pn(this, fq, pbra, pket, lbra, sbra, zbra, lket, sket, zket)
    use MyLibrary
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    integer :: ibra, iket, i, m, mu_Tz
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:), v3(:), v4(:)
    real(8) :: tau1, tau2, tau_x_tau

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))
    allocate(v4(this%GetNMesh()))

    mu_Tz = zket-zbra
    tau1      = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_m, 1, mu_Tz, phase=-1) 
    tau2      = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau2_m, 1, mu_Tz, phase=-1) 
    tau_x_tau = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_cross_tau2, 1, mu_Tz, phase=-1) * (-1.d0)

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0; v4(:) = 0.d0
    call qdep_term_q( this, fq, v1 )
    call qdep_term_q_x_sigma( this, fq, v2 )
    call qdep_term_sigma( this, fq, v3 )
    !if(relativistic_correction) call qdep_term_p_s_dot_q( this, fq, v4 ) ! this term is not tested

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do m = -1, 1
          do i = 1, this%GetNMesh()
            fq%op(i,m,ibra,iket) = v1(i) * &
                & (s1_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau1 + &
                &  s2_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau2 )
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v2(i) * tau_x_tau * &
                & (s1_dot_q_q_x_s2(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i) - &
                &  s2_dot_q_q_x_s1(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i))
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v3(i) * &
                & (s1m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau1 + &
                &  s2m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau2)
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v4(i) * tau_x_tau * &
                & p_s_dot_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)
          end do
        end do
      end do
    end do
    deallocate(v1, v2, v3, v4)

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
  end subroutine set_helicity_rep_pn

  subroutine set_helicity_rep_isospin(this, fq, pbra, pket, tbra, tket)
    use MyLibrary, only: pi, tau1_iso, tau2_iso, tau1_tau2_tensor_iso
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: tbra, tket
    integer :: ibra, iket, i, m
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:), v3(:), v4(:)
    real(8) :: tau1, tau2, tau_x_tau

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))
    allocate(v4(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0; v4(:) = 0.d0
    call qdep_term_q( this, fq, v1 )
    call qdep_term_q_x_sigma( this, fq, v2 )
    call qdep_term_sigma( this, fq, v3 )
    !if(relativistic_correction) call qdep_term_p_s_dot_q( this, fq, v4 ) ! this term is not tested
    tau1 = tau1_iso(tbra,tket)
    tau2 = tau2_iso(tbra,tket)
    tau_x_tau = tau1_tau2_tensor_iso(tbra,tket,1) * sqrt(2.d0) ! tau_x_tau = i(t1 x t2) = sqrt(2) [t1 t2]_1

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do m = -1, 1
          do i = 1, this%GetNMesh()
            fq%op(i,m,ibra,iket) = v1(i) * &
                & (s1_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau1 + &
                &  s2_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau2 )
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v2(i) * tau_x_tau * &
                & (s1_dot_q_q_x_s2(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i) - &
                &  s2_dot_q_q_x_s1(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i))
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v3(i) * &
                & (s1m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau1 + &
                &  s2m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau2)
            fq%op(i,m,ibra,iket) = fq%op(i,m,ibra,iket) + &
                & v4(i) * tau_x_tau * &
                & p_s_dot_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)
          end do
        end do
      end do
    end do
    deallocate(v1, v2, v3, v4)

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
  end subroutine set_helicity_rep_isospin

  subroutine qdep_term_q(this, fq, v)
    use MyLibrary, only: g_A, f_pi
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer :: i

    if(this%ch_order < 2) return
    do i = 1, this%GetNMesh()
      v(i) = v(i) + g_A / f_pi**2 * fq%pi_prop(i) * (this%c3*1.d-3)
    end do
    if(this%ch_order < 3) return
  end subroutine qdep_term_q

  subroutine qdep_term_q_x_sigma(this, fq, v)
    use MyLibrary, only: g_A, f_pi, m_nucleon
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer :: i

    if(this%ch_order < 2) return
    do i = 1, this%GetNMesh()
       if(relativistic_correction) v(i) = v(i) + 0.5d0 * g_A / f_pi**2 * fq%pi_prop(i) * (this%c4*1.d-3 + 0.25d0 / m_nucleon) ! w/ relativistic correction
       if(.not. relativistic_correction) v(i) = v(i) + 0.5d0 * g_A / f_pi**2 * fq%pi_prop(i) * this%c4*1.d-3
    end do
    if(this%ch_order < 3) return
  end subroutine qdep_term_q_x_sigma

  subroutine qdep_term_sigma(this, fq, v)
    use MyLibrary, only: f_pi, lambda_chi, g_A, m_nucleon
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer :: i

    if(this%ch_order < 2) return
    do i = 1, this%GetNMesh()
      v(i) = v(i) - 0.25d0 * this%cD / (lambda_chi*f_pi**2)
    end do
    if(this%ch_order < 3) return
  end subroutine qdep_term_sigma

  subroutine qdep_term_p_s_dot_q(this, fq, v)
    use MyLibrary, only: g_A, f_pi, m_nucleon
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer :: i

    if(this%ch_order < 2) return
    do i = 1, this%GetNMesh()
      v(i) = v(i) - g_A / (8.d0 * m_nucleon * f_pi**2) * fq%pi_prop(i)
    end do
    if(this%ch_order < 3) return
  end subroutine qdep_term_p_s_dot_q

  function s1m(lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i) result(res)
    !
    ! s1_m
    !
    use PartialWaveDecomposition, only: helicity_expectation_value_1, helicity_expectation_value_sigma
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8) :: res
    res = helicity_expectation_value_sigma(lam1_bra, lam1_ket, i, m) * &
        & helicity_expectation_value_1(-lam2_bra,-lam2_ket, i)
  end function s1m

  function s2m(lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i) result(res)
    !
    ! s2_m
    !
    use PartialWaveDecomposition, only: helicity_expectation_value_1, helicity_expectation_value_sigma
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8) :: res
    res = helicity_expectation_value_1(lam1_bra, lam1_ket, i) * &
        & helicity_expectation_value_sigma(-lam2_bra,-lam2_ket, i, m)
  end function s2m

  function s1_dot_q_q(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, m, i) result(res)
    !
    ! (s1 . q) q_m
    !
    use PartialWaveDecomposition, only: identity_helicity_rep, get_y1
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8), intent(in) :: pbra, pket
    real(8) :: res
    res = y1_fact_inv * (&
        & lam1_bra*pbra*get_y1(i, m) * pbra - &
        & lam1_bra*pbra*y1_theta0(m) * pket - &
        & lam1_ket*pket*get_y1(i, m) * pbra + &
        & lam1_ket*pket*y1_theta0(m) * pket) * &
        & identity_helicity_rep(lam1_bra, lam2_bra, lam1_ket, lam2_ket, i)
  end function s1_dot_q_q

  function s2_dot_q_q(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, m, i) result(res)
    !
    ! (s2 . q) q_m
    !
    use PartialWaveDecomposition, only: get_y1, identity_helicity_rep
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8), intent(in) :: pbra, pket
    real(8) :: res
    res =-y1_fact_inv * ( &
        & lam2_bra*pbra*get_y1(i, m) * pbra - &
        & lam2_bra*pbra*y1_theta0(m) * pket - &
        & lam2_ket*pket*get_y1(i, m) * pbra + &
        & lam2_ket*pket*y1_theta0(m) * pket) * &
        & identity_helicity_rep(lam1_bra, lam2_bra, lam1_ket, lam2_ket, i)
  end function s2_dot_q_q

  function s1_dot_q_q_x_s2(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, m, i) result(res)
    !
    ! -i(s1 . q) (q x s2)_m
    !
    use PartialWaveDecomposition, only: get_y1, helicity_expectation_value_1, helicity_expectation_value_sigma
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8), intent(in) :: pbra, pket
    real(8) :: res
    integer :: m1, m2
    res = 0.d0
    do m1 = -1, 1
      m2 = m-m1
      if(abs(m2)>1) cycle
      res = res + y1_fact_inv * (&
          & lam1_bra*pbra*get_y1(i, m1) * pbra - &
          & lam1_bra*pbra*y1_theta0(m1) * pket - &
          & lam1_ket*pket*get_y1(i, m1) * pbra + &
          & lam1_ket*pket*y1_theta0(m1) * pket) * &
          & helicity_expectation_value_1(lam1_bra, lam1_ket, i) * &
          & helicity_expectation_value_sigma(-lam2_bra,-lam2_ket, i, m2) * &
          & cgs%get(2, 2*m1, 2, 2*m2, 2, 2*m) * (-sqrt(2.d0))
    end do
  end function s1_dot_q_q_x_s2

  function s2_dot_q_q_x_s1(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, m, i) result(res)
    !
    ! -i(s2 . q) (q x s1)_m
    !
    use PartialWaveDecomposition, only: helicity_expectation_value_1, helicity_expectation_value_sigma, get_y1
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8), intent(in) :: pbra, pket
    real(8) :: res
    integer :: m1, m2
    res = 0.d0
    do m1 = -1, 1
      m2 = m-m1
      if(abs(m2)>1) cycle
      res = res - y1_fact_inv * (&
          & lam2_bra*pbra*get_y1(i, m1) * pbra - &
          & lam2_bra*pbra*y1_theta0(m1) * pket - &
          & lam2_ket*pket*get_y1(i, m1) * pbra + &
          & lam2_ket*pket*y1_theta0(m1) * pket) * &
          & helicity_expectation_value_sigma(lam1_bra, lam1_ket, i, m2) * &
          & helicity_expectation_value_1(-lam2_bra,-lam2_ket, i) * &
          & cgs%get(2, 2*m1, 2, 2*m2, 2, 2*m) * (-sqrt(2.d0))
    end do
  end function s2_dot_q_q_x_s1

  function p_s_dot_q(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, m, i) result(res)
    !
    ! -[(p1 + p1') (s2 . q) + (p2 + p2') (s1 . q)] = -[(P/2 + p + p') (s2 . q) + (P/2 - p - p') (s1 . q)]
    ! (s1-s2) . q (p' + p) + (s1 + s2) . q P
    ! Since the isospin part is t1 x t2, i.e., T' != T, S' != S has to be satisfied (operator is parity even).
    ! From < S'=0 | (s1 + s2) | S=1 > = 0, P dependent term vanishes
    !
    ! (s1-s2) . q (p' + p)_m
    !
    use PartialWaveDecomposition, only: get_y1, sigma_minus_q
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, m, i
    real(8), intent(in) :: pbra, pket
    real(8) :: res
    res = sigma_minus_q(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i)
    res = res * (pbra*get_y1(i,m) + pket*y1_theta0(m))
    res = res * y1_fact_inv
  end function p_s_dot_q

  function y1_theta0(m) result(res)
    integer, intent(in) :: m
    real(8) :: res
    res = 0.d0
    if(m == 0) res = y1_fact
  end function y1_theta0

  subroutine set_mom_functions(this, fq, pbra, pket)
    use MyLibrary, only: m_pi
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    real(8) :: q, s, l
    integer :: i

    allocate(fq%op(this%GetNMesh(), -1:1, this%GetNumHeli(), this%GetNumHeli()))
    allocate(fq%QMesh(this%GetNMesh()))
    allocate(fq%pi_prop(this%GetNMesh()))
    allocate(fq%loop_L(this%GetNMesh()))
    allocate(fq%loop_H(this%GetNMesh()))
    allocate(fq%loop_A(this%GetNMesh()))
    allocate(fq%s_var(this%GetNMesh()))

    fq%op(:,:,:,:) = 0.d0
    do i = 1, this%GetNMesh()
      q = sqrt(pbra**2 + pket**2 - 2.d0*pbra*pket*this%ZMesh(i))
      fq%QMesh(i) = q
      fq%pi_prop(i) = 1.d0 / (q**2 + m_pi**2)
      fq%s_var(i) = sqrt(4.d0*m_pi*m_pi + q*q)
    end do

    if(SFRegulator) then
      do i = 1, this%GetNMesh()
        q = fq%QMesh(i)
        s = fq%s_var(i)
        if( LambdaSFR < 2*m_pi ) cycle
        l = sqrt(LambdaSFR**2 - 4.d0*m_pi**2)
        fq%loop_L(i) = 0.5d0*(s/q) * log( (LambdaSFR**2*s**2 + q**2*l**2+ 2.d0*q*s*l*LambdaSFR) / &
            &(4.d0*m_pi**2*(LambdaSFR**2 + q**2) ) )
        fq%loop_H(i) = 4.d0 * (m_pi**2/s**2) * fq%loop_L(i)
        fq%loop_A(i) = 0.5d0 / q * atan( q*(LambdaSFR-2.d0*m_pi) / (q**2+2.d0*LambdaSFR*m_pi) )
      end do
      return
    end if

    do i = 1, this%GetNMesh()
      q = fq%QMesh(i)
      s = fq%s_var(i)
      fq%loop_L(i) = 0.5d0 * (s/q) * log( (s+q) / (s-q) )
      fq%loop_H(i) = 4.d0 * (m_pi**2/s**2) * fq%loop_L(i)
      fq%loop_A(i) = 0.5d0/q * atan( 0.5d0 * q / m_pi )
    end do
  end subroutine set_mom_functions

  subroutine release_mom_functions(fq)
    type(MomFunctions), intent(inout) :: fq
    deallocate(fq%op)
    deallocate(fq%QMesh)
    deallocate(fq%pi_prop)
    deallocate(fq%loop_L)
    deallocate(fq%loop_H)
    deallocate(fq%loop_A)
    deallocate(fq%s_var)
  end subroutine release_mom_functions

!  !
!  ! test functions
!  !
!  function c3_term_iso(this, pbra, lbra, sbra, jbra, tbra, pket, lket, sket, jket, tket) result(res)
!    use MyLibrary, only: tau1_iso, tau2_iso, dcg, pi, g_A, f_pi
!    type(AxialCurrentQ0), intent(in) :: this
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, tbra, lket, sket, jket, tket
!    real(8) :: res, tau1, tau2, falpha
!    integer :: alpha
!
!    tau1 = tau1_iso(tbra,tket)
!    tau2 = tau2_iso(tbra,tket)
!    res = 0.d0
!    falpha = 0.d0
!    do alpha = 0, 2, 2
!      falpha = falpha + dcg(2, 0, 2, 0, 2*alpha, 0) * &
!          & (general_qfunction(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, 1, 0, 1, 1, ffunc) * tau1 + &
!          &  general_qfunction(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, 0, 1, 1, 1, ffunc) * tau2)
!    end do
!    res = falpha * (-sqrt(4*pi / 3)) * g_A * this%c3*1.d-3 / f_pi**2
!    res = res / (2.d0 * pi)**3
!  contains
!    function ffunc(q) result(res)
!      use MyLibrary, only: m_pi
!      real(8), intent(in) :: q
!      real(8) :: res
!      res = q**2 / (q**2 + m_pi**2)
!    end function ffunc
!  end function c3_term_iso
!
!  function c4_term_iso(this, pbra, lbra, sbra, jbra, tbra, pket, lket, sket, jket, tket) result(res)
!    use MyLibrary, only: tau1_tau2_tensor_iso, dcg, sjs, pi, g_A, f_pi
!    type(AxialCurrentQ0), intent(in) :: this
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, tbra, lket, sket, jket, tket
!    real(8) :: res, tau_x_tau, falpha
!    integer :: alpha
!
!    tau_x_tau = tau1_tau2_tensor_iso(tbra,tket,1) * sqrt(2.d0) ! tau_x_tau = i(t1 x t2) = sqrt(2) [t1 t2]_1
!    res = 0.d0
!    falpha = 0.d0
!    do alpha = 0, 2, 2
!      falpha = falpha + dcg(2, 0, 2, 0, 2*alpha, 0) * sjs(2, 2, 2*alpha, 2, 2, 2) * tau_x_tau * &
!          & general_qfunction(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, 1, 1, 1, 1, ffunc)
!    end do
!    res = falpha * sqrt(24.d0*pi) * g_A * this%c4*1.d-3 / f_pi**2
!    res = res / (2.d0 * pi)**3
!  contains
!    function ffunc(q) result(res)
!      use MyLibrary, only: m_pi
!      real(8), intent(in) :: q
!      real(8) :: res
!      res = q**2 / (q**2 + m_pi**2)
!    end function ffunc
!  end function c4_term_iso
!
!  function cD_term_iso(this, pbra, lbra, sbra, jbra, tbra, pket, lket, sket, jket, tket) result(res)
!    use MyLibrary, only: tau1_iso, tau2_iso, f_pi, lambda_chi, snj, pi
!    type(AxialCurrentQ0), intent(in) :: this
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, tbra, lket, sket, jket, tket
!    real(8) :: res, tau1, tau2
!
!    tau1 = tau1_iso(tbra,tket)
!    tau2 = tau2_iso(tbra,tket)
!    res = 0.d0
!    if(lbra /= lket) return
!    if(lket /= 0) return
!    res = 3.d0 * sqrt(dble((2*jbra+1)*(2*jket+1)*(2*lket+1)*(2*sbra+1)*(2*sket+1))) * &
!        & snj(2*lbra, 2*sbra, 2*jbra, 2*lket, 2*sket, 2*jket, 0, 2, 2) * &
!        & (snj(1, 1, 2*sbra, 1, 1, 2*sket, 2, 0, 2)*2.d0*sqrt(3.d0)*tau1 + &
!        &  snj(1, 1, 2*sbra, 1, 1, 2*sket, 0, 2, 2)*2.d0*sqrt(3.d0)*tau2)
!    res = res * (-0.25d0) * this%cD / (Lambda_chi * f_pi**2) * 4.d0*pi
!    res = res / (2.d0 * pi)**3
!  end function cD_term_iso
!
!  function Pterm_iso(this, pbra, lbra, sbra, jbra, tbra, pket, lket, sket, jket, tket) result(res)
!    use MyLibrary, only: tau1_tau2_tensor_iso, dcg, pi, g_A, f_pi, m_nucleon
!    type(AxialCurrentQ0), intent(in) :: this
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, tbra, lket, sket, jket, tket
!    real(8) :: res, tau_x_tau, falpha, tmp
!    integer :: alpha
!
!    tau_x_tau = tau1_tau2_tensor_iso(tbra,tket,1) * sqrt(2.d0) ! tau_x_tau = i(t1 x t2) = sqrt(2) [t1 t2]_1
!    res = 0.d0
!    falpha = 0.d0
!    do alpha = 0, min(2, lbra+lket)
!      tmp = (-4.d0 * pi) / sqrt(3.d0) * pbra**2 * dcg(2, 0, 2, 0, 2*alpha, 0) * &
!           general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, 0, alpha, 1, 0, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '1: ', tmp
!
!      tmp =  (4.d0 * pi) / sqrt(3.d0) * pbra * pket * sqrt(dble(2*alpha+1)) / 3.d0 * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 1, 1, alpha, 1, 0, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '2: ', tmp
!
!      tmp = (4.d0 * pi) / sqrt(3.d0) * pbra**2 * dcg(2, 0, 2, 0, 2*alpha, 0) * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, 0, alpha, 0, 1, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '3: ', tmp
!
!      tmp = (-4.d0 * pi) / sqrt(3.d0) * pbra * pket * sqrt(dble(2*alpha+1)) / 3.d0 * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 1, 1, alpha, 0, 1, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '4: ', tmp
!
!      tmp = (-4.d0 * pi) / sqrt(3.d0) * pbra * pket * sqrt(dble(2*alpha+1)) / 3.d0 * (-1.d0)**alpha * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 1, 1, alpha, 1, 0, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '5: ', tmp
!
!      tmp = (4.d0 * pi) / sqrt(3.d0) * pket**2 * dcg(2, 0, 2, 0, 2*alpha, 0) * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 0, alpha, alpha, 1, 0, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '6: ', tmp
!
!      tmp = (4.d0 * pi) / sqrt(3.d0) * pbra * pket * sqrt(dble(2*alpha+1)) / 3.d0 * (-1.d0)**alpha * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 1, 1, alpha, 0, 1, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '7: ', tmp
!
!      tmp = (-4.d0 * pi) / sqrt(3.d0) * pket**2 * dcg(2, 0, 2, 0, 2*alpha, 0) * &
!          & general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, 0, alpha, alpha, 0, 1, 1, 1, ffunc)
!      falpha = falpha + tmp
!      write(*,*) '8: ', tmp
!    end do
!    res = -falpha * tau_x_tau * g_A / (8.d0 * m_nucleon * f_pi**2)
!    res = res / (2.d0 * pi)**3
!  contains
!    function ffunc(q) result(res)
!      use MyLibrary, only: m_pi
!      real(8), intent(in) :: q
!      real(8) :: res
!      res = 1.d0 / (q**2 + m_pi**2)
!    end function ffunc
!  end function Pterm_iso
!
!  function general_qfunction(pbra, lbra, sbra, jbra, pket, lket, sket, jket, alpha, is1, is2, beta, K, func) result(res)
!    !
!    ! <p'l'S'J' || f(q) [[Y_alpha(q)] [s1 s2]_beta]_K || plSJ>
!    !
!    use MyLibrary, only: tau1_iso, tau2_iso, tau1_tau2_tensor_iso, sjs, dcg, pi, gamma_function, snj
!    interface
!      function func(q) result(r)
!        real(8), intent(in) :: q
!        real(8) :: r
!      end function func
!    end interface
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, alpha, beta, is1, is2, K
!    real(8) :: res
!    real(8) :: spin_part, orb_part, falpha, fL
!    integer :: alpha1, alpha2, L
!
!    res = 0.d0
!    spin_part = 0.d0
!    if(is1 == 0 .and. is2 == 0 .and. sbra==sket) spin_part = sqrt(dble(2*sket+1))
!    if(is1 == 1 .and. is2 == 0) spin_part = tau1_iso(sbra,sket)
!    if(is1 == 0 .and. is2 == 1) spin_part = tau2_iso(sbra,sket)
!    if(is1 == 1 .and. is2 == 1) spin_part = tau1_tau2_tensor_iso(sbra,sket,beta)
!    if(abs(spin_part) < 1.d-8) return
!
!    orb_part = 0.d0
!    do alpha1 = 0, alpha
!      alpha2 = alpha - alpha1
!      falpha = sqrt(4.d0 * pi * gamma_function(dble(2*alpha+2)) / &
!          & (gamma_function(dble(2*alpha1+2)) * gamma_function(dble(2*alpha2+2)))) * &
!          & pbra ** alpha1 * pket ** alpha2 * &
!          & sqrt(dble((2*alpha+1)*(2*alpha1+1)*(2*alpha2+1))) * &
!          & (-1.d0)**(lbra + lket + alpha)
!
!      do L = max(abs(lbra-alpha1), abs(lket-alpha2)), min(lbra+alpha1, lket+alpha2)
!        fL = sjs(2*alpha1, 2*lbra, 2*L, 2*lket, 2*alpha2, 2*alpha) * &
!            & dcg(2*L, 0, 2*alpha1, 0, 2*lbra, 0) * &
!            & dcg(2*L, 0, 2*alpha2, 0, 2*lket, 0) * decomposition(pbra, pket, L)
!        orb_part = orb_part + falpha * fL
!      end do
!    end do
!
!    res = sqrt(dble((2*jbra+1)*(2*jket+1)*(2*K+1))) * &
!        & snj(2*lbra, 2*sbra, 2*jbra, 2*lket, 2*sket, 2*jket, 2*alpha, 2*beta, 2*K) * &
!        & orb_part * spin_part * (-1.d0)**((lbra - lket)/2)
!  contains
!    function decomposition(pbra, pket, L) result(res)
!      use MyLibrary, only: gauss_legendre, legendre_polynomial
!      real(8), intent(in) :: pbra, pket
!      integer, intent(in) :: L
!      integer :: i, nm
!      real(8) :: res, q
!      real(8), allocatable :: z(:), w(:)
!      nm = 50
!      call gauss_legendre(-1.d0, 1.d0, z, w, nm)
!      res = 0.d0
!      do i = 1, nm
!        q = sqrt(pbra**2 + pket**2 - 2.d0 * z(i) * pbra * pket)
!        res = res + legendre_polynomial(L,z(i)) * func(q) * q**(-alpha) * w(i)
!      end do
!      res = res * dble(2*L+1) * 0.5d0
!      deallocate(z, w)
!    end function decomposition
!  end function general_qfunction
!
!  function general_function(pbra, lbra, sbra, jbra, pket, lket, sket, jket, &
!        & alpha1, alpha2, alpha, is1, is2, beta, K, func) result(res)
!    !
!    ! <p'l'S'J' || f(q) [[Y_alpha1(p') Y_alpha2(p)]_alpha [s1 s2]_beta]_K || plSJ>
!    !
!    use MyLibrary, only: tau1_iso, tau2_iso, tau1_tau2_tensor_iso, sjs, dcg, pi, gamma_function, snj
!    interface
!      function func(q) result(r)
!        real(8), intent(in) :: q
!        real(8) :: r
!      end function func
!    end interface
!    real(8), intent(in) :: pbra, pket
!    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, alpha1, alpha2, alpha, is1, is2, beta, K
!    real(8) :: res
!    real(8) :: spin_part, orb_part, falpha, fL
!    integer :: L
!
!    res = 0.d0
!    spin_part = 0.d0
!    if(is1 == 0 .and. is2 == 0 .and. sbra==sket) spin_part = sqrt(dble(2*sket+1))
!    if(is1 == 1 .and. is2 == 0) spin_part = tau1_iso(sbra,sket)
!    if(is1 == 0 .and. is2 == 1) spin_part = tau2_iso(sbra,sket)
!    if(is1 == 1 .and. is2 == 1) spin_part = tau1_tau2_tensor_iso(sbra,sket,beta)
!    if(abs(spin_part) < 1.d-8) return
!
!    orb_part = 0.d0
!    falpha = sqrt(dble((2*alpha+1)*(2*alpha1+1)*(2*alpha2+1))) * &
!        & (-1.d0)**(lbra + lket + alpha2 + alpha)
!    do L = max(abs(lbra-alpha1), abs(lket-alpha2)), min(lbra+alpha1, lket+alpha2)
!      fL = (-1.d0)**L * sjs(2*alpha1, 2*lbra, 2*L, 2*lket, 2*alpha2, 2*alpha) * &
!          & dcg(2*L, 0, 2*alpha1, 0, 2*lbra, 0) * &
!          & dcg(2*L, 0, 2*alpha2, 0, 2*lket, 0) * decomposition(pbra, pket, L)
!      orb_part = orb_part + falpha * fL
!    end do
!
!    res = sqrt(dble((2*jbra+1)*(2*jket+1)*(2*K+1))) * &
!        & snj(2*lbra, 2*sbra, 2*jbra, 2*lket, 2*sket, 2*jket, 2*alpha, 2*beta, 2*K) * &
!        & orb_part * spin_part * (-1.d0)**((lbra - lket)/2)
!  contains
!    function decomposition(pbra, pket, L) result(res)
!      use MyLibrary, only: gauss_legendre, legendre_polynomial
!      real(8), intent(in) :: pbra, pket
!      integer, intent(in) :: L
!      integer :: i, nm
!      real(8) :: res, q
!      real(8), allocatable :: z(:), w(:)
!      nm = 50
!      call gauss_legendre(-1.d0, 1.d0, z, w, nm)
!      res = 0.d0
!      do i = 1, nm
!        q = sqrt(pbra**2 + pket**2 - 2.d0 * z(i) * pbra * pket)
!        res = res + legendre_polynomial(L,z(i)) * func(q) * w(i)
!      end do
!      res = res * dble(2*L+1) * 0.5d0
!      deallocate(z, w)
!    end function decomposition
!  end function general_function
end module AxialVectorQ0
