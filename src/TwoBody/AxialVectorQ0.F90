!
! This is for an independent test of the two-body multipole operators.
! A clear test is needed, but it looks like I can reproduce the figure 2 in Eur. Phys. J. A 54, 76 (2018).
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
contains

  subroutine show_options(this)
    class(AxialCurrentQ0), intent(in) :: this
    write(*,"(a)") "# Calculation options for axial current operator @ Q=0"
    write(*,"(2a)") "# Operator: ", trim(this%OpName)
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a)") "# c3= ", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    write(*,"(a,f8.4)") "# cD= ", this%cD
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
        & Regulator, RegulatorPower, RegulatorLambda, LamSFR)
    use MyLibrary, only: pi
    class(AxialCurrentQ0), intent(inout) :: this
    character(*), intent(in) :: OpName
    real(8), intent(in) :: LECs(:)
    real(8), intent(in), optional :: RegulatorLambda, LamSFR
    integer, intent(in), optional :: ch_order, Jmax, RegulatorPower
    character(*), intent(in), optional :: Regulator
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
    class(AxialCurrentQ0), intent(in) :: this
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
        if(z1bras(ibra)==1 .and. z2bras(ibra)==-1) phbra = (-1.d0)**(lbra+sbra)
        do iket = 1, size(z1kets)
          phket = 1.d0
          if(z1kets(iket)==1 .and. z2kets(iket)==-1) phket = (-1.d0)**(lket+sket)
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
    use MyLibrary, only: pi, tau_1, tau_m, tau1_cross_tau2
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket
    integer :: ibra, iket, i, m
    integer :: lams_bra(2), lams_ket(2)
    real(8), allocatable :: v1(:), v2(:), v3(:)
    real(8) :: tau1, tau2, tau_x_tau

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0
    call qdep_term_q( this, fq, v1 )
    call qdep_term_q_x_sigma( this, fq, v2 )
    call qdep_term_sigma( this, fq, v3 )
    tau1 = tau_m(z1bra,z1ket,(z1ket-z1bra)/2,phase=-1) * tau_1(z2bra,z2ket)
    tau2 = tau_m(z2bra,z2ket,(z2ket-z2bra)/2,phase=-1) * tau_1(z1bra,z1ket)
    tau_x_tau = tau1_cross_tau2(z1bra,z2bra,z1ket,z2ket,1,(z1ket+z2ket-z1bra-z2bra)/2,phase=-1) * (-1.d0) ! -1 is from i^2

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do m = -1, 1
          do i = 1, this%GetNMesh()
            fq%op(i,m,ibra,iket) = &
                & v1(i) * &
                & (s1_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau1 + &
                &  s2_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau2 ) + &
                & v2(i) * tau_x_tau * &
                & (s1_dot_q_q_x_s2(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i) - &
                &  s2_dot_q_q_x_s1(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)) + &
                & v3(i) * &
                & (s1m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau1 + &
                &  s2m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau2)
          end do
        end do
      end do
    end do
    deallocate(v1, v2, v3)

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
    real(8), allocatable :: v1(:), v2(:), v3(:)
    real(8) :: tau1, tau2, tau_x_tau

    allocate(v1(this%GetNMesh()))
    allocate(v2(this%GetNMesh()))
    allocate(v3(this%GetNMesh()))

    v1(:) = 0.d0; v2(:) = 0.d0; v3(:) = 0.d0
    call qdep_term_q( this, fq, v1 )
    call qdep_term_q_x_sigma( this, fq, v2 )
    call qdep_term_sigma( this, fq, v3 )
    tau1 = tau1_iso(tbra,tket)
    tau2 = tau2_iso(tbra,tket)
    tau_x_tau = tau1_tau2_tensor_iso(tbra,tket,1) * sqrt(2.d0)

    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do m = -1, 1
          do i = 1, this%GetNMesh()
            fq%op(i,m,ibra,iket) = v1(i) * &
                & (s1_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau1 + &
                &  s2_dot_q_q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)*tau2 ) + &
                & v2(i) * tau_x_tau * &
                & (s1_dot_q_q_x_s2(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i) - &
                &  s2_dot_q_q_x_s1(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, m, i)) + &
                & v3(i) * &
                & (s1m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau1 + &
                &  s2m(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), m, i)*tau2)
          end do
        end do
      end do
    end do
    deallocate(v1, v2, v3)

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
    use MyLibrary, only: g_A, f_pi
    type(AxialCurrentQ0), intent(in) :: this
    type(MomFunctions), intent(in) :: fq
    real(8), intent(inout) :: v(:)
    integer :: i

    if(this%ch_order < 2) return
    do i = 1, this%GetNMesh()
      v(i) = v(i) + 0.5d0 * g_A / f_pi**2 * fq%pi_prop(i) * (this%c4*1.d-3)
    end do
    if(this%ch_order < 3) return
  end subroutine qdep_term_q_x_sigma

  subroutine qdep_term_sigma(this, fq, v)
    use MyLibrary, only: f_pi, lambda_chi
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
    ! (s1 . q) (q x s2)_m, i is not taken into account here
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
    ! (s2 . q) (q x s1)_m, i is not taken into account here
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
end module AxialVectorQ0
