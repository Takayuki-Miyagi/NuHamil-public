!
! Reference:
! [1] J. de Vries, E. Epelbaum, L. Girlanda, A. Gnech, E. Mereghetti, and M. Viviani, Front. Phys. 8, (2020).
!
module ParityTimeReViolation
  !$ use omp_lib
  use PartialWaveDecomposition, only: PWD, &
      & identity_helicity_rep, sigma_dot_sigma, sigma1q_sigma2q, &
      & sigma_plus_q, sigma_plus_QQ, sigma_minus_q, sigma_minus_QQ, sigma_cross_q, sigma_cross_QQ
  implicit none

  ! pvtv_couplings in chiral EFT
  ! idx,  coupling name,      note
  !   1,           gpi0, LO  (1PE),     Isoscalar
  !   2,           gpi1, LO  (1PE),     Isovector
  !   3,           gpi2, LO  (1PE),     Isotensor
  !   4,          delta, NLO (3PE),     Isovector
  !   5,            Ct1, N2LO(Contact), Isoscalar
  !   6,            Ct2, N2LO(Contact), Isoscalar
  !   7,            Ct3, N2LO(Contact), Isovector
  !   8,            Ct4, N2LO(Contact), Isovector
  !   9,            Ct5, N2LO(Contact), Isotensor

  ! pvtv_couplings in meson-exchange theory
  ! idx,  coupling name,      note
  !   1,           gpi0,   Pion, Isoscalar
  !   2,           gpi1,   Pion, Isovector
  !   3,           gpi2,   Pion, Isotensor
  !   4,          grho0,    rho, Isoscalar
  !   5,          grho1,    rho, Isovector
  !   6,          grho2,    rho, Isotensor
  !   7,          geta0,    eta, Isoscalar
  !   8,          geta1,    eta, Isovector
  !   9,        gomega0,  omega, Isoscalar
  !  10,        gomega1,  omega, Isovector
  type, extends(PWD) :: PTViolation
    integer :: ch_order = 0
    character(:), allocatable :: OpName
    character(:), allocatable :: NNint
    real(8) :: c1, c2, c3, c4
    real(8) :: pvtv_couplings(10) = 0.d0
    logical :: meson_exchange = .false.
  contains
    procedure :: InitPTViolation
    procedure :: FinPTViolation
    procedure :: show_pvtv_options
    procedure :: calc_matrix_element
    !procedure :: set_helicity_rep
    generic :: init => InitPTViolation
    generic :: fin => FinPTViolation
  end type PTViolation

  logical, private :: SFRegulator = .false.
  real(8), private :: LambdaSFR
  type :: MomFunctions
    real(8), allocatable :: op(:,:,:,:)
    real(8), allocatable :: QMesh(:)
    real(8), allocatable :: pi_prop(:)
    real(8), allocatable :: eta_prop(:)
    real(8), allocatable :: rho_prop(:)
    real(8), allocatable :: omega_prop(:)
    real(8), allocatable :: loop_L(:)
    real(8), allocatable :: loop_H(:)
    real(8), allocatable :: loop_A(:)
    real(8), allocatable :: s_var(:)
  end type MomFunctions
contains

  subroutine show_pvtv_options(this)
    class(PTViolation), intent(in) :: this
    write(*,"(a)") "# Calculation options for parity and time-reversal violation operator"
    write(*,"(2a)") "# Operator: ", trim(this%OpName)
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a)") "# c1= ", this%c1, " GeV-1, c2= ", &
      & this%c2, " GeV-1, c3=", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    if(this%meson_exchange) then
      write(*,'(a, f8.4)') "#    gpi0 = ", this%pvtv_couplings( 1)
      write(*,'(a, f8.4)') "#    gpi1 = ", this%pvtv_couplings( 2)
      write(*,'(a, f8.4)') "#    gpi2 = ", this%pvtv_couplings( 3)
      write(*,'(a, f8.4)') "#   grho0 = ", this%pvtv_couplings( 4)
      write(*,'(a, f8.4)') "#   grho1 = ", this%pvtv_couplings( 5)
      write(*,'(a, f8.4)') "#   grho2 = ", this%pvtv_couplings( 6)
      write(*,'(a, f8.4)') "#   geta0 = ", this%pvtv_couplings( 7)
      write(*,'(a, f8.4)') "#   geta1 = ", this%pvtv_couplings( 8)
      write(*,'(a, f8.4)') "# gomega0 = ", this%pvtv_couplings( 9)
      write(*,'(a, f8.4)') "# gomega1 = ", this%pvtv_couplings(10)
    else
      write(*,'(a, f8.4)') "#    gpi0 = ", this%pvtv_couplings(1)
      write(*,'(a, f8.4)') "#    gpi1 = ", this%pvtv_couplings(2)
      write(*,'(a, f8.4)') "#    gpi2 = ", this%pvtv_couplings(3)
      write(*,'(a, f8.4)') "#   delta = ", this%pvtv_couplings(4)
      write(*,'(a, f8.4)') "#     Ct1 = ", this%pvtv_couplings(5)
      write(*,'(a, f8.4)') "#     Ct2 = ", this%pvtv_couplings(6)
      write(*,'(a, f8.4)') "#     Ct3 = ", this%pvtv_couplings(7)
      write(*,'(a, f8.4)') "#     Ct4 = ", this%pvtv_couplings(8)
      write(*,'(a, f8.4)') "#     Ct5 = ", this%pvtv_couplings(9)
    end if
    write(*,*)
    call this%show_pwd_options()
  end subroutine show_pvtv_options

  subroutine FinPTViolation(this)
    class(PTViolation), intent(inout) :: this
    this%ch_order = 0
    this%OpName = "None"
    this%c1 = 0.d0
    this%c2 = 0.d0
    this%c3 = 0.d0
    this%c4 = 0.d0
    call this%fin_pwd()
  end subroutine FinPTViolation

  subroutine InitPTViolation(this, OpName, NNint, pvtv_couplings, Rank, ch_order, Jmax, &
        & Regulator, RegulatorPower, RegulatorLambda, meson_exchange, LamSFR)
    class(PTViolation), intent(inout) :: this
    character(*), intent(in) :: OpName, NNint
    real(8), intent(in) :: pvtv_couplings(10)
    real(8), intent(in), optional :: RegulatorLambda, LamSFR
    integer, intent(in), optional :: ch_order, Jmax, RegulatorPower, Rank
    character(*), intent(in), optional :: Regulator
    logical, intent(in), optional :: meson_exchange
    this%OpName = OpName
    this%NNint = NNint
    this%pvtv_couplings = pvtv_couplings

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
  end subroutine InitPTViolation

  function calc_matrix_element(this, pbra, lbra, sbra, jbra, zbra, pket, lket, sket, jket, zket, pn_formalism) result(r)
    class(PTViolation), intent(in) :: this
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: lbra, sbra, jbra, zbra, lket, sket, jket, zket
    logical, intent(in) :: pn_formalism
    real(8) :: r, phbra, phket, norm
    integer :: ibra, iket
    integer, allocatable :: z1bras(:), z2bras(:), z1kets(:), z2kets(:)
    type(MomFunctions) :: fq

    if(zbra==-1) then
      allocate(z1bras(1), z2bras(1))
      z1bras = [-1]
      z2bras = [-1]
    elseif(zbra==1) then
      allocate(z1bras(1), z2bras(1))
      z1bras = [1]
      z2bras = [1]
    elseif(zbra==0) then
      allocate(z1bras(2), z2bras(2))
      z1bras = [-1,1]
      z2bras = [1,-1]
    end if

    if(zket==-1) then
      allocate(z1kets(1), z2kets(1))
      z1kets = [-1]
      z2kets = [-1]
    elseif(zket==1) then
      allocate(z1kets(1), z2kets(1))
      z1kets = [1]
      z2kets = [1]
    elseif(zket==0) then
      allocate(z1kets(2), z2kets(2))
      z1kets = [-1,1]
      z2kets = [1,-1]
    end if

    call set_mom_functions(this, fq, pbra, pket)
    r = 0.d0
    if(pn_formalism) then
      norm = 1.d0 / sqrt(dble(size(z1bras) * size(z1kets)) )
      do ibra = 1, size(z1bras)
        phbra = 1.d0
        if(z1bras(ibra)== 1 .and. z2bras(ibra)==-1) phbra = (-1.d0)**(lbra+sbra)
        do iket = 1, size(z1kets)
          phket = 1.d0
          if(z1kets(iket)== 1 .and. z2kets(iket)==-1) phket = (-1.d0)**(lket+sket)
          fq%op(:,:,:,:) = 0.d0
          call set_helicity_rep_pn(this, fq, pbra, pket, z1bras(ibra), z2bras(ibra), z1kets(iket), z2kets(iket))
          r = r + this%do_pwd(fq%op, pbra, lbra, sbra, jbra, pket, lket, sket, jket) * phbra * phket
        end do
      end do
      r = r * norm
    else
      fq%op(:,:,:,:) = 0.d0
      call set_helicity_rep_isospin(this, fq, pbra, pket, zbra, zket)
      r = this%do_pwd(fq%op, pbra, lbra, sbra, jbra, pket, lket, sket, jket)
    end if

    call release_mom_functions(fq)
    deallocate(z1bras, z2bras, z1kets, z2kets)
  end function calc_matrix_element

  subroutine set_helicity_rep_pn(this, fq, pbra, pket, z1bra, z2bra, z1ket, z2ket)
    use MyLibrary, only: pi
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket
    integer :: ibra, iket, i
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:)

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0
    call sigma_plus_q_term_pn( this, fq, v1, z1bra, z2bra, z1ket, z2ket)
    call sigma_minus_q_term_pn(this, fq, v2, z1bra, z2bra, z1ket, z2ket)

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = v1(i) * sigma_plus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              &                 v2(i) * sigma_minus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)
        end do
      end do
    end do

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
    deallocate(v1, v2)
  end subroutine set_helicity_rep_pn

  subroutine set_helicity_rep_isospin(this, fq, pbra, pket, tbra, tket)
    use MyLibrary, only: pi
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: tbra, tket
    integer :: ibra, iket, i
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:)

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0
    call sigma_plus_q_term_isospin( this, fq, v1, tbra, tket)
    call sigma_minus_q_term_isospin(this, fq, v2, tbra, tket)

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = v1(i) * sigma_plus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) + &
              &                 v2(i) * sigma_minus_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i)
        end do
      end do
    end do

    fq%op(:,:,:,:) = fq%op(:,:,:,:) / (2.d0*pi)**3
    deallocate(v1, v2)
  end subroutine set_helicity_rep_isospin

  subroutine sigma_plus_q_term_pn(this, fq, v, z1bra, z2bra, z1ket, z2ket)
    use MyLibrary, only: tau1_cross_tau2, tau1_dot_tau2, tau1_plus_tau2, tau1_tau2_tensor, tau_1, tau1_minus_tau2, &
        & tau_x, tau_y, tau_z
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_cross_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1) * (-1.d0)
    t_dot   = tau1_dot_tau2(z1bra,z2bra,z1ket,z2ket,0,0,phase=-1)
    t_plus  = tau1_plus_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1)
    t_tensor= 2.d0 * tau_z(z1bra,z1ket,phase=-1) * tau_z(z2bra,z2ket,phase=-1) - &
        & tau_x(z1bra,z1ket) * tau_x(z2bra,z2ket) + tau_y(z1bra,z1ket,phase=-1) * tau_y(z2bra,z2ket,phase=-1)
    t_1     = tau_1(z1bra,z1ket) * tau_1(z2bra,z2ket)
    t_minus = tau1_minus_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1)
    call qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_q_term_pn

  subroutine sigma_plus_q_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0)
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_plus_q_term_isospin

  subroutine qdep_sigma_plus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: pi, m_pi, m_nucleon, lambda_chi, g_pi, g_rho, g_eta, g_omega, g_A=>g_A_GT, f_pi
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        v(i) = v(i) + g_pi * this%pvtv_couplings(2) / (4.d0 * m_nucleon) * t_minus * fq%pi_prop(i)
        v(i) = v(i) + g_rho * this%pvtv_couplings(5) / (4.d0 * m_nucleon) * t_minus * fq%rho_prop(i)
        v(i) = v(i) - g_eta * this%pvtv_couplings(8) / (4.d0 * m_nucleon) * t_minus * fq%eta_prop(i)
        v(i) = v(i) - g_omega * this%pvtv_couplings(10) / (4.d0 * m_nucleon) * t_minus * fq%omega_prop(i)
      end do
    else
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 0.25d0 * g_A * this%pvtv_couplings(2) / f_pi * t_minus * fq%pi_prop(i)
      end do
      if(this%ch_order < 1) return

      ! three-pion exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) - 1.25d0 * g_A**3 * this%pvtv_couplings(4) * m_nucleon * pi / (f_pi * lambda_chi**2) * &
            & t_minus * fq%pi_prop(i) * ((1.d0 - 2.d0*m_pi**2/fq%s_var(i)**2)*fq%s_var(i)**2*fq%loop_A(i) + m_pi)
      end do

      if(this%ch_order < 2) return

      ! Contact
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 1.d0 / (Lambda_chi**2 * f_pi) * 0.5d0 * t_minus * (this%pvtv_couplings(7) - this%pvtv_couplings(8))
      end do

      ! three-pion  exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 2.5d0 * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c1*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_minus * 4.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i)
        v(i) = v(i) - (5.d0/6.d0) * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c2*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_minus * (2.d0 * fq%loop_L(i) + 6.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i))
        v(i) = v(i) - 1.25d0 * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c3*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_minus * (3.d0 * fq%loop_L(i) + 5.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i))
      end do
    end if
  end subroutine qdep_sigma_plus_q

  subroutine sigma_minus_q_term_pn(this, fq, v, z1bra, z2bra, z1ket, z2ket)
    use MyLibrary, only: tau1_cross_tau2, tau1_dot_tau2, tau1_plus_tau2, tau1_tau2_tensor, tau_1, tau1_minus_tau2, &
        & tau_x, tau_y, tau_z
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_cross_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1) * (-1.d0)
    t_dot   = tau1_dot_tau2(z1bra,z2bra,z1ket,z2ket,0,0,phase=-1)
    t_plus  = tau1_plus_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1)
    t_tensor= 2.d0 * tau_z(z1bra,z1ket,phase=-1) * tau_z(z2bra,z2ket,phase=-1) - &
        & tau_x(z1bra,z1ket) * tau_x(z2bra,z2ket) + tau_y(z1bra,z1ket,phase=-1) * tau_y(z2bra,z2ket,phase=-1)
    t_1     = tau_1(z1bra,z1ket) * tau_1(z2bra,z2ket)
    t_minus = tau1_minus_tau2(z1bra,z2bra,z1ket,z2ket,1,0,phase=-1)
    call qdep_sigma_minus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_minus_q_term_pn

  subroutine sigma_minus_q_term_isospin(this, fq, v, tbra, tket)
    use MyLibrary, only: tau1_tau2_tensor_iso, tau1_plus_tau2_iso, tau1_minus_tau2_iso, tau_identity_iso
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer, intent(in) :: tbra, tket
    real(8) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    t_cross = tau1_tau2_tensor_iso(tbra,tket,1) * (-sqrt(2.d0)) * (-1.d0)
    t_dot   = tau1_tau2_tensor_iso(tbra,tket,0) * (-sqrt(3.d0))
    t_tensor= tau1_tau2_tensor_iso(tbra,tket,2) * sqrt(6.d0)
    t_1     = tau_identity_iso(tbra,tket)
    t_plus  = tau1_plus_tau2_iso(tbra,tket) *  (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    t_minus = tau1_minus_tau2_iso(tbra,tket) * (-1.d0) ! -1 is from conversion from particle- to nuclear-physics convention
    call qdep_sigma_minus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
  end subroutine sigma_minus_q_term_isospin

  subroutine qdep_sigma_minus_q(this, fq, v, t_1, t_dot, t_plus, t_minus, t_cross, t_tensor)
    use MyLibrary, only: pi, m_pi, m_nucleon, lambda_chi, g_pi, g_eta, g_omega, g_rho, g_A=>g_A_GT, f_pi
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    real(8), intent(in) :: t_cross, t_dot, t_plus, t_tensor, t_1, t_minus
    integer :: i

    if(this%meson_exchange) then
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 0.5d0 * g_pi/m_nucleon * (this%pvtv_couplings(1)*t_dot + 0.5d0 * this%pvtv_couplings(2) * t_plus + &
            & this%pvtv_couplings(3) * t_tensor) * fq%pi_prop(i)
        v(i) = v(i) - 0.5d0 * g_rho/m_nucleon * (this%pvtv_couplings(4)*t_dot + 0.5d0 * this%pvtv_couplings(5) * t_plus + &
            & this%pvtv_couplings(6) * t_tensor) * fq%rho_prop(i)
        v(i) = v(i) + 0.5d0 * g_eta/m_nucleon * (this%pvtv_couplings(7)*t_dot + 0.5d0 * this%pvtv_couplings(8) * t_plus) * &
            & fq%eta_prop(i)
        v(i) = v(i) - 0.5d0 * g_omega/m_nucleon * (this%pvtv_couplings(9)*t_dot + 0.5d0 * this%pvtv_couplings(10) * t_plus) * &
            & fq%omega_prop(i)
      end do
    else
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 0.5d0 * g_A * this%pvtv_couplings(1) / f_pi * t_dot * fq%pi_prop(i)
        v(i) = v(i) + 0.25d0 * g_A * this%pvtv_couplings(2) / f_pi * t_plus * fq%pi_prop(i)
        v(i) = v(i) + (1.d0/6.d0) * g_A * this%pvtv_couplings(3) / f_pi * t_tensor * fq%pi_prop(i)
      end do


      if(this%ch_order < 1) return

      ! three pion
      do i = 1, this%GetNMesh()
        v(i) = v(i) - 1.25d0 * g_A**3 * this%pvtv_couplings(4) * m_nucleon * pi / (f_pi * lambda_chi**2) * &
            & t_plus * fq%pi_prop(i) * ((1.d0 - 2.d0*m_pi**2/fq%s_var(i)**2) * fq%s_var(i)**2 * fq%loop_A(i) + m_pi)
      end do

      if(this%ch_order < 2) return
      ! Contact
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 1.d0 / (Lambda_chi**2 * f_pi) * &
            & (this%pvtv_couplings(5)*t_1 + this%pvtv_couplings(6)*t_dot + &
            & (this%pvtv_couplings(7)+this%pvtv_couplings(8))*t_plus*0.5d0 + this%pvtv_couplings(9)*t_tensor)
      end do

      ! Two-Pion Exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) + g_A * this%pvtv_couplings(1) / (f_pi * lambda_chi**2) * t_dot * fq%loop_L(i)
        v(i) = v(i) + g_A**3 * this%pvtv_couplings(1) / (f_pi * lambda_chi**2) * t_dot * (fq%loop_H(i) - 3.d0*fq%loop_L(i))
        v(i) = v(i) - (1.d0/3.d0) * g_A * this%pvtv_couplings(3) / &
            & (f_pi * lambda_chi**2) * t_tensor * fq%loop_L(i)
        v(i) = v(i) - (1.d0/3.d0) * g_A**3 * this%pvtv_couplings(3) / &
            & (f_pi * lambda_chi**2) * t_tensor * (fq%loop_H(i) - 3.d0*fq%loop_L(i))
      end do

      ! three-pion  exchange
      do i = 1, this%GetNMesh()
        v(i) = v(i) + 2.5d0 * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c1*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_plus * 4.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i)
        v(i) = v(i) - (5.d0/6.d0) * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c2*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_plus * (2.d0 * fq%loop_L(i) + 6.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i))
        v(i) = v(i) - 1.25d0 * g_A * this%pvtv_couplings(4) * m_nucleon * (this%c3*1.d-3) / (f_pi * lambda_chi**2) * &
            & t_plus * (3.d0 * fq%loop_L(i) + 5.d0 * m_pi**2 * fq%pi_prop(i) * fq%loop_L(i))
      end do
    end if
  end subroutine qdep_sigma_minus_q

  function rc_helicity_term_plus(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) result(r)
    !
    ! (q . s1) [(q x Q) . s2] + (q . s2) [(q x Q) . s1]
    !
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, i
    real(8), intent(in) :: pbra, pket
    real(8) :: r
    r = rc_helicity_term1(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) + &
        & rc_helicity_term2(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i)
  end function rc_helicity_term_plus

  function rc_helicity_term_minus(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) result(r)
    !
    ! (q . s1) [(q x Q) . s2] + (q . s2) [(q x Q) . s1]
    !
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, i
    real(8), intent(in) :: pbra, pket
    real(8) :: r
    r = rc_helicity_term1(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) - &
        & rc_helicity_term2(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i)
  end function rc_helicity_term_minus

  function rc_helicity_term1(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) result(r)
    !
    ! (q . s1) [(q x Q) . s2]
    ! q = p'-p, Q = (p'+p)/4
    !
    use PartialWaveDecomposition, only: helicity_expectation_value_1, helicity_expectation_value_y, get_zmesh
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, i
    real(8), intent(in) :: pbra, pket
    real(8) :: r
    real(8) :: s, z
    z = get_zmesh(i)
    s = sqrt(1.d0 - z**2)
    r = 0.5d0 * (pbra*lam1_bra - pket*lam1_ket) * pbra * pket * (-s) * &
        & helicity_expectation_value_1( lam1_bra, lam1_ket, i) * &
        & helicity_expectation_value_y(-lam2_bra,-lam2_ket, i)
  end function rc_helicity_term1

  function rc_helicity_term2(lam1_bra, lam2_bra, lam1_ket, lam2_ket, pbra, pket, i) result(r)
    !
    ! (q . s2) [(q x Q) . s1]
    ! q = p'-p, Q = (p'+p)/4
    !
    use PartialWaveDecomposition, only: helicity_expectation_value_1, helicity_expectation_value_y, get_zmesh
    integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, i
    real(8), intent(in) :: pbra, pket
    real(8) :: r
    real(8) :: s, z
    z = get_zmesh(i)
    s = sqrt(1.d0 - z**2)
    r = 0.5d0 * (-pbra*lam2_bra + pket*lam2_ket) * pbra * pket * (-s) * &
        & helicity_expectation_value_y( lam1_bra, lam1_ket, i) * &
        & helicity_expectation_value_1(-lam2_bra,-lam2_ket, i)
  end function rc_helicity_term2

  subroutine set_mom_functions(this, fq, pbra, pket)
    use MyLibrary, only: m_pi, m_eta, m_rho, m_omega
    type(PTViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    real(8) :: q, s, l
    integer :: i

    allocate(fq%op(this%GetNMesh(), -this%GetRank():this%GetRank(), this%GetNumHeli(), this%GetNumHeli()))
    allocate(fq%QMesh(this%GetNMesh()))
    allocate(fq%pi_prop(this%GetNMesh()))
    allocate(fq%eta_prop(this%GetNMesh()))
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
      fq%eta_prop(i) = 1.d0 / (q**2 + m_eta**2)
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
    deallocate(fq%eta_prop)
    deallocate(fq%rho_prop)
    deallocate(fq%omega_prop)
    deallocate(fq%loop_L)
    deallocate(fq%loop_H)
    deallocate(fq%loop_A)
    deallocate(fq%s_var)
  end subroutine release_mom_functions
end module ParityTimeReViolation
