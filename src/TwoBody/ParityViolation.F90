! Important note: after the partial-wave decomposition all the matrix elements are given in unit of i
! Reference:
! [1] J. de Vries, E. Epelbaum, L. Girlanda, A. Gnech, E. Mereghetti, and M. Viviani, Front. Phys. 8, (2020).
!
module ParityViolation
  !$ use omp_lib
  use PartialWaveDecomposition, only: PWD, &
      & identity_helicity_rep, sigma_dot_sigma, sigma1q_sigma2q, &
      & sigma_plus_q, sigma_plus_QQ, sigma_minus_q, sigma_minus_QQ, sigma_cross_q, sigma_cross_QQ
  implicit none

  ! physics constants
  real(8), parameter, private :: chi_V = 3.7d0
  real(8), parameter, private :: chi_S =-0.12d0

  ! pv_couplings in chiral EFT
  ! idx,  coupling name,      note
  !   1,           hpi1, LO  (OPE),     Isovector
  !   2,            Ct1, NLO (Contact), Isoscalar
  !   3,            Ct2, NLO (Contact), Isoscalar
  !   4,            Ct3, NLO (Contact), Isovector
  !   5,            Ct4, NLO (Contact), Isovector
  !   6,            Ct5, NLO (Contact), Isotensor
  !   7,            hV0, N2LO,          Isoscalar
  !   8,            hV1, N2LO           Isovector
  !   9,            hV2, N2LO           Isotensor
  !  10,            hA1, N2LO           Isovector
  !  11,            hA2, N2LO           Isotensor

  ! pv_couplings in meson-exchange theory (DDH pot.)
  ! idx,  coupling name,      note
  !   1,           hpi1,   Pion, Isovector
  !   2,          hrho0,    Rho, Isoscalar
  !   3,          hrho1,    rho, Isovector
  !   4,         hrho1',    rho, Isovector
  !   5,          hrho2,    rho, Isotensor
  !   6,        homega0,  omega, Isoscalar
  !   7,        homega1,  omega, Isovector
  type, extends(PWD) :: PViolation
    integer :: ch_order = 0
    character(:), allocatable :: OpName
    character(:), allocatable :: NNint
    real(8) :: c1, c2, c3, c4
    real(8) :: pv_couplings(11) = 0.d0
    logical :: meson_exchange = .false.
  contains
    procedure :: InitPViolation
    procedure :: FinPViolation
    procedure :: show_pv_options
    procedure :: calc_matrix_element
    !procedure :: set_helicity_rep
    generic :: init => InitPViolation
    generic :: fin => FinPViolation
  end type PViolation

  logical, private :: SFRegulator = .false.
  real(8), private :: LambdaSFR
  type :: MomFunctions
    real(8), allocatable :: op(:,:,:,:)
    real(8), allocatable :: QMesh(:)
    real(8), allocatable :: pi_prop(:)
    real(8), allocatable :: rho_prop(:)
    real(8), allocatable :: omega_prop(:)
    real(8), allocatable :: loop_L(:)
    real(8), allocatable :: loop_H(:)
    real(8), allocatable :: loop_A(:)
    real(8), allocatable :: s_var(:)
  end type MomFunctions
contains

  subroutine show_pv_options(this)
    class(PViolation), intent(in) :: this
    write(*,"(a)") "# Calculation options for parity violation operator"
    write(*,"(2a)") "# Operator: ", trim(this%OpName)
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a)") "# c1= ", this%c1, " GeV-1, c2= ", &
      & this%c2, " GeV-1, c3=", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    if(this%meson_exchange) then
      write(*,'(a, es12.4)') "#    hpi1 = ", this%pv_couplings(1)
      write(*,'(a, es12.4)') "#   hrho0 = ", this%pv_couplings(2)
      write(*,'(a, es12.4)') "#   hrho1 = ", this%pv_couplings(3)
      write(*,'(a, es12.4)') "#  hrho1' = ", this%pv_couplings(4)
      write(*,'(a, es12.4)') "#   hrho2 = ", this%pv_couplings(5)
      write(*,'(a, es12.4)') "# homega0 = ", this%pv_couplings(6)
      write(*,'(a, es12.4)') "# homega1 = ", this%pv_couplings(7)
    else
      write(*,'(a, es12.4)') "# hpi1 = ", this%pv_couplings( 1)
      write(*,'(a, es12.4)') "#  Ct1 = ", this%pv_couplings( 2)
      write(*,'(a, es12.4)') "#  Ct2 = ", this%pv_couplings( 3)
      write(*,'(a, es12.4)') "#  Ct3 = ", this%pv_couplings( 4)
      write(*,'(a, es12.4)') "#  Ct4 = ", this%pv_couplings( 5)
      write(*,'(a, es12.4)') "#  Ct5 = ", this%pv_couplings( 6)
      write(*,'(a, es12.4)') "#  hV0 = ", this%pv_couplings( 7)
      write(*,'(a, es12.4)') "#  hV1 = ", this%pv_couplings( 8)
      write(*,'(a, es12.4)') "#  hV2 = ", this%pv_couplings( 9)
      write(*,'(a, es12.4)') "#  hA1 = ", this%pv_couplings(10)
      write(*,'(a, es12.4)') "#  hA2 = ", this%pv_couplings(11)
    end if
    write(*,*)
    call this%show_pwd_options()
  end subroutine show_pv_options

  subroutine FinPViolation(this)
    class(PViolation), intent(inout) :: this
    this%ch_order = 0
    this%OpName = "None"
    this%c1 = 0.d0
    this%c2 = 0.d0
    this%c3 = 0.d0
    this%c4 = 0.d0
    call this%fin_pwd()
  end subroutine FinPViolation

  subroutine InitPViolation(this, OpName, NNint, pv_couplings, Rank, ch_order, Jmax, &
        & Regulator, RegulatorPower, RegulatorLambda, meson_exchange, LamSFR)
    class(PViolation), intent(inout) :: this
    character(*), intent(in) :: OpName, NNint
    real(8), intent(in) :: pv_couplings(11)
    real(8), intent(in), optional :: RegulatorLambda, LamSFR
    integer, intent(in), optional :: ch_order, Jmax, RegulatorPower, Rank
    character(*), intent(in), optional :: Regulator
    logical, intent(in), optional :: meson_exchange
    this%OpName = OpName
    this%NNint = NNint
    this%pv_couplings = pv_couplings

    call this%init_pwd(Jmax=Jmax, Rank=Rank, Regulator=Regulator, &
        & RegulatorPower=RegulatorPower, RegulatorLambda=RegulatorLambda)
    if( present(ch_order) ) this%ch_order = ch_order
    if( present(meson_exchange)) this%meson_exchange = meson_exchange
    if( present(LamSFR)) then
      SFRegulator = .true.
      LambdaSFR = LamSFR
    end if

    select case(this%NNint)
    case("N3LO_EM500")
      this%c1 = -0.81d0
      this%c2 =  2.8d0
      this%c3 = -3.2d0
      this%c4 =  5.4d0
    case("LO_EMN500", "NLO_EMN500")
      this%c1 = 0.d0
      this%c2 = 0.d0
      this%c3 = 0.d0
      this%c4 = 0.d0
    case("N2LO_EMN500")
      this%c1 = -0.74d0
      this%c2 =  0.d0
      this%c3 = -3.61d0
      this%c4 =  2.44d0
    case("N3LO_EMN500")
      this%c1 = -1.07d0
      this%c2 =  3.2d0
      this%c3 = -5.32d0
      this%c4 =  3.56d0
    case("N4LO_EMN500")
      this%c1 = -1.10d0
      this%c2 =  3.57d0
      this%c3 = -5.54d0
      this%c4 =  4.17d0
    case default
      write(*,*) "Unknown NNint:", __FILE__, __LINE__
    end select
  end subroutine InitPViolation

  function calc_matrix_element(this, pbra, lbra, sbra, jbra, zbra, pket, lket, sket, jket, zket, pn_formalism) result(r)
    class(PViolation), intent(in) :: this
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
    use MyLibrary, only: pi
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    integer :: ibra, iket, i
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:), v3(:), v4(:)

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))
    allocate(v4(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0; v4(:) = 0.d0
    call sigma_plus_q_term_pn(  this, fq, v1, lbra, sbra, zbra, lket, sket, zket)
    call sigma_cross_q_term_pn( this, fq, v2, lbra, sbra, zbra, lket, sket, zket)
    call sigma_minus_QQ_term_pn(this, fq, v3, lbra, sbra, zbra, lket, sket, zket)
    call sigma_plus_QQ_term_pn( this, fq, v4, lbra, sbra, zbra, lket, sket, zket)

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = v1(i) * sigma_plus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              & v2(i) * sigma_cross_q( lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              & v3(i) * sigma_minus_QQ(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)*0.5d0 + &
              & v4(i) * sigma_plus_QQ( lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)*0.5d0
        end do
      end do
    end do

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
    deallocate(v1, v2, v3, v4)
  end subroutine set_helicity_rep_pn

  subroutine set_helicity_rep_isospin(this, fq, pbra, pket, tbra, tket)
    use MyLibrary, only: pi
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: tbra, tket
    integer :: ibra, iket, i
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:), v3(:), v4(:)

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))
    allocate(v4(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0; v4(:) = 0.d0
    call sigma_plus_q_term_isospin(  this, fq, v1, tbra, tket)
    call sigma_cross_q_term_isospin( this, fq, v2, tbra, tket)
    call sigma_minus_QQ_term_isospin(this, fq, v3, tbra, tket)
    call sigma_plus_QQ_term_isospin( this, fq, v4, tbra, tket)

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = v1(i) * sigma_plus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              & v2(i) * sigma_cross_q( lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              & v3(i) * sigma_minus_QQ(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)*0.5d0 + &
              & v4(i) * sigma_plus_QQ( lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)*0.5d0
        end do
      end do
    end do

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
    deallocate(v1, v2, v3, v4)
  end subroutine set_helicity_rep_isospin

  subroutine sigma_plus_q_term_pn(this, fq, v, lbra, sbra, zbra, lket, sket, zket)
    use MyLibrary
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_cross_tau2, 1, 0, phase=-1) * (-1.d0)
    t_dot   = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_dot_tau2, 1, 0, phase=-1)
    t_plus  = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_plus_tau2, 1, 0, phase=-1)
    t_tensor= asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_tau2_tensor, 2, 0, phase=-1) * sqrt(6.d0)
    t_1     = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau_identity, 0, 0, phase=-1) 
    t_minus = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_minus_tau2, 1, 0, phase=-1) 
    call qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_q_term_pn

  subroutine sigma_plus_q_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0) ! -1 is from particle -> nuclear physics convention
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_q_term_isospin

  subroutine qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: pi, m_pi, m_nucleon, lambda_chi, g_pi, g_rho, g_A=>g_A_GT, f_pi
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        ! Note there is the sign difference in pion contribution, see [1] C. H. Hyun and B. Desplanques, Phys. Lett. B 552, 41 (2003).
        v(i) = v(i) - g_pi * this%pv_couplings(1) / (sqrt(8.d0) * m_nucleon) * t_cross * fq%pi_prop(i)
        v(i) = v(i) + g_rho * this%pv_couplings(4) / (2.d0 * m_nucleon) * t_cross * fq%rho_prop(i)
      end do
    else
      do i = 1, this%GetNMesh()
        v(i) = v(i) - g_A * this%pv_couplings(1) / (sqrt(8.d0) * f_pi) * t_cross * fq%pi_prop(i)
      end do
      if(this%ch_order < 1) return

      ! Contact
      do i = 1, this%GetNMesh()
        v(i) = v(i) - this%pv_couplings(4) / (Lambda_chi * Lambda_chi * f_pi) * t_cross
      end do

      ! Two-Pion Exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) + g_A * this%pv_couplings(1) / (sqrt(8.d0) * f_pi * Lambda_chi**2) * &
            & t_cross * fq%loop_L(i) - &
            & g_A**3 * this%pv_couplings(1) /(sqrt(8.d0) * f_pi * Lambda_chi**2) * &
            & t_cross * ( fq%loop_H(i) - 3.d0*fq%loop_L(i) )
      end do

      if(this%ch_order < 2) return
      ! Two-Pion Exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) + g_A**3 * pi * this%pv_couplings(8) * 0.25d0 / (f_pi**2 * Lambda_chi**2) * t_cross * &
            & (1.d0 - 2.d0 * m_pi**2 / fq%s_var(i)**2) * fq%s_var(i)**2 * fq%loop_A(i)
      end do
    end if
  end subroutine qdep_sigma_plus_q

  subroutine sigma_cross_q_term_pn(this, fq, v, lbra, sbra, zbra, lket, sket, zket)
    use MyLibrary
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_cross_tau2, 1, 0, phase=-1) * (-1.d0)
    t_dot   = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_dot_tau2, 1, 0, phase=-1)
    t_plus  = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_plus_tau2, 1, 0, phase=-1)
    t_tensor= asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_tau2_tensor, 2, 0, phase=-1) * sqrt(6.d0)
    t_1     = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau_identity, 0, 0, phase=-1) 
    t_minus = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_minus_tau2, 1, 0, phase=-1) 
    call qdep_sigma_cross_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_cross_q_term_pn

  subroutine sigma_cross_q_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0) ! -1 is from particle -> nuclear physics convention
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_cross_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_cross_q_term_isospin

  subroutine qdep_sigma_cross_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: pi, m_pi, m_nucleon, lambda_chi, g_omega, g_rho, g_A=>g_A_GT, f_pi
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        v(i) = v(i) - g_rho / m_nucleon * (1.d0 + chi_v) * fq%rho_prop(i) * &
            & (t_dot * this%pv_couplings(2) + t_plus * this%pv_couplings(3) * 0.5d0 + &
            & t_tensor * this%pv_couplings(5) * 0.5d0 / sqrt(6.d0))
        v(i) = v(i) - g_omega / m_nucleon * (1.d0 + chi_s) * fq%omega_prop(i) * &
            & (t_1 * this%pv_couplings(6) + t_plus * this%pv_couplings(7) * 0.5d0)
      end do
    else
      if(this%ch_order < 1) return
      ! Contact
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 1.d0 / (Lambda_chi * Lambda_chi * f_pi) * &
            & (this%pv_couplings(2) * t_1 + &
            &  this%pv_couplings(3) * t_dot + &
            &  this%pv_couplings(5) * t_plus + &
            &  this%pv_couplings(6) * t_tensor)
      end do

      ! Two-Pion Exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) - g_A**3 * this%pv_couplings(1) * sqrt(2.d0) / (f_pi * Lambda_chi**2) * t_plus * fq%loop_L(i)
      end do

      if(this%ch_order < 2) return
      ! Three-Pion Exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) - (this%c4*1.d-3) * this%pv_couplings(1) * g_A * pi / (sqrt(2.d0) * f_pi * Lambda_chi**2) * &
            & t_plus * fq%s_var(i)**2 * fq%loop_A(i) - &
            & g_A**2 * pi * 0.5d0 / (f_pi**2 * Lambda_chi**2) * (&
            & 0.75d0*g_A*this%pv_couplings(7) + 0.5d0*g_A*this%pv_couplings(7)*t_dot + &
            & (0.25d0*g_A*this%pv_couplings(8) - this%pv_couplings(10)) * t_plus - &
            & (this%pv_couplings(11)+g_A*this%pv_couplings(9)/3.d0) * t_tensor) * fq%s_var(i)**2 * fq%loop_A(i)
      end do
    end if
  end subroutine qdep_sigma_cross_q

  subroutine sigma_minus_QQ_term_pn(this, fq, v, lbra, sbra, zbra, lket, sket, zket)
    use MyLibrary
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_cross_tau2, 1, 0, phase=-1) * (-1.d0)
    t_dot   = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_dot_tau2, 1, 0, phase=-1)
    t_plus  = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_plus_tau2, 1, 0, phase=-1)
    t_tensor= asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_tau2_tensor, 2, 0, phase=-1) * sqrt(6.d0)
    t_1     = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau_identity, 0, 0, phase=-1) 
    t_minus = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_minus_tau2, 1, 0, phase=-1) 
    call qdep_sigma_minus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_minus_QQ_term_pn

  subroutine sigma_minus_QQ_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0) ! -1 is from particle -> nuclear physics convention
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_minus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_minus_QQ_term_isospin

  subroutine qdep_sigma_minus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: pi, m_nucleon, g_omega, g_rho
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        v(i) = v(i) - 2.d0 * g_rho / m_nucleon * fq%rho_prop(i) * &
            & (t_dot * this%pv_couplings(2) + t_plus * this%pv_couplings(3) * 0.5d0 + &
            & t_tensor * this%pv_couplings(5) * 0.5d0 / sqrt(6.d0))
        v(i) = v(i) - 2.d0 * g_omega / m_nucleon * fq%omega_prop(i) * &
            & (t_1 * this%pv_couplings(6) + t_plus * this%pv_couplings(7) * 0.5d0)
      end do
    else
      return
    end if
  end subroutine qdep_sigma_minus_QQ

  subroutine sigma_plus_QQ_term_pn(this, fq, v, lbra, sbra, zbra, lket, sket, zket)
    use MyLibrary
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_cross_tau2, 1, 0, phase=-1) * (-1.d0)
    t_dot   = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_dot_tau2, 1, 0, phase=-1)
    t_plus  = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_plus_tau2, 1, 0, phase=-1)
    t_tensor= asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_tau2_tensor, 2, 0, phase=-1) * sqrt(6.d0)
    t_1     = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau_identity, 0, 0, phase=-1) 
    t_minus = asym_isospin_func_pn(lbra, sbra, zbra, lket, sket, zket, tau1_minus_tau2, 1, 0, phase=-1) 
    call qdep_sigma_plus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_QQ_term_pn

  subroutine sigma_plus_QQ_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0) ! -1 is from particle -> nuclear physics convention
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_plus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_QQ_term_isospin

  subroutine qdep_sigma_plus_QQ(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: m_nucleon, g_omega, g_rho
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i
    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        v(i) = v(i) + g_rho / m_nucleon * fq%rho_prop(i) * t_minus * this%pv_couplings(3)
        v(i) = v(i) - g_omega / m_nucleon * fq%omega_prop(i) * t_minus * this%pv_couplings(7)
      end do
    else
      return
    end if
  end subroutine qdep_sigma_plus_QQ

  subroutine set_mom_functions(this, fq, pbra, pket)
    use MyLibrary, only: m_pi, m_rho, m_omega
    type(PViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    real(8) :: q, s, l
    integer :: i

    allocate(fq%op(this%GetNMesh(), -this%GetRank():this%GetRank(), this%GetNumHeli(), this%GetNumHeli()))
    allocate(fq%QMesh(this%GetNMesh()))
    allocate(fq%pi_prop(this%GetNMesh()))
    allocate(fq%rho_prop(this%GetNMesh()))
    allocate(fq%omega_prop(this%GetNMesh()))
    allocate(fq%loop_L(this%GetNMesh()))
    allocate(fq%loop_H(this%GetNMesh()))
    allocate(fq%loop_A(this%GetNMesh()))
    allocate(fq%s_var(this%GetNMesh()))

    fq%op(:,:,:,:) = 0.d0
    do i = 1, this%GetNMesh()
      q = sqrt(pbra**2 + pket**2 - 2.d0*pbra*pket*this%ZMesh(i))
      fq%QMesh(i) = q
      fq%pi_prop(i) = 1.d0 / (q**2 + m_pi**2)
      fq%rho_prop(i) = 1.d0 / (q**2 + m_rho**2)
      fq%omega_prop(i) = 1.d0 / (q**2 + m_omega**2)
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
    deallocate(fq%rho_prop)
    deallocate(fq%omega_prop)
    deallocate(fq%loop_L)
    deallocate(fq%loop_H)
    deallocate(fq%loop_A)
    deallocate(fq%s_var)
  end subroutine release_mom_functions
end module ParityViolation

!program test
!  use ParityViolation
!  type(PViolation) :: check
!  integer :: lbra = 0, sbra = 1, jbra = 1
!  integer :: lket = 1, sket = 0, jket = 1
!  real(8) :: pv_couplings(11)
!  real(8) :: pbra = 10.d0, pket = 10.d0
!  real(8) :: me
!  pv_couplings = [1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
!  !call check%init("0vbb-Fermi", "", pv_couplings, meson_exchange=.true.)
!  call check%init("0vbb-Fermi", "", pv_couplings, ch_order=0)
!  call check%show_pv_options()
!  !call check%set_helicity_rep(pbra, pket, -1, 1)
!  !me = check%do_pwd(pbra, lbra, sbra, jbra, pket, lket, sket, jket)
!  me = check%calc_matrix_element(pbra, lbra, sbra, jbra, 0, pket, lket, sket, jket, 0, .true.)
!  write(*,'(6i4,es18.8)') lbra, sbra, jbra, lket, sket, jket, me*197.32705d0**3
!end program test
