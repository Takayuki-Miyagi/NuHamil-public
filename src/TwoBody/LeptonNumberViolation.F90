!
! module for neutrinoless double beta decay matrix element
! See the following for references::
! [1] F. Šimkovic, A. Faessler, H. Müther, V. Rodin, and M. Stauf, Phys. Rev. C 79, 055501 (2009).
! [2] J. Engel and J. Menéndez, Reports Prog. Phys. 80, 046301 (2017).
! [3] V. Cirigliano, W. Dekens, E. Mereghetti, and A. Walker-Loud, Phys. Rev. C 97, 1 (2018).
! [4] R. Wirth, J. M. Yao, and H. Hergert, arXiv:2105.05415 (2021).
!
! Fermi:
!   O_{F}(q) = \frac{1}{2\pi^{2}} \frac{h_{F}(q)}{q(q+E_{c})} \tau_{1,+}\tau_{2,+}
! GT:
!   O_{GT}(q) = \frac{1}{2\pi^{2}} \frac{h_{GT}(q)}{q(q+E_{c})} (\sigma_{1}\cdot \sigma_{2})\tau_{1,+}\tau_{2,+}
! Tensor:
!   O_{T}(q) = -\frac{1}{2\pi^{2}} \frac{h_{T}(q)}{q(q+E_{c})} [3(\sigma_{1}\cdot \hat{\mathbf q})(\sigma_{2}\cdot\hat{\mathbf q})-(\sigma_{1}\cdot \sigma_{2})]\tau_{1,+}\tau_{2,+}
!
! Neutrino potentials:
!    h_{F} &= \frac{g_{V}^{2}(q)}{g_{V}^{2}} \\
!    h_{GT} &= \frac{1}{g_{A}^{2}}\left[g_{A}^{2}(q) - \frac{g_{A}(q)g_{P}(q) q^{2}}{3m_{N}} + \frac{g^{2}_{P}(q) q^{4}}{12m_{N}}
!       + \frac{g^{2}_{M}(q) q^{2}}{6m_{N}} \right] \\
!    h_{T} &= \frac{1}{g_{A}^{2}}\left[ \frac{g_{A}(q)g_{P}(q) q^{2}}{3m_{N}} - \frac{g^{2}_{P}(q) q^{4}}{12m_{N}}
!       + \frac{g^{2}_{M}(q) q^{2}}{12m_{N}} \right]
! with
!    g_{V}(q) = \frac{g_{V}}{(1 + q^{2}/\Lambda_{V}^{2})^{2}}
!    g_{A}(q) = \frac{g_{A}}{(1 + q^{2}/\Lambda_{A}^{2})^{2}}
!    g_{P}(q) = -\frac{2m_{N}g_{A}(q)}{q^{2}+m_{\pi}^{2}}
!    g_{M}(q) = \frac{\mu_{p}-\mu{n}}{\mu_{N}} g_{V}(q)
! Note that we should have the following at the leading order [see JHEP 82 2017 for the details]:
!    g_{V}(q) = g_{V} = 1, g_{A}(q) = g_{A}, g_{M}(q) = 1 + \kappa = 4.7
!    g_{P}(q) = -g_{A} \frac{2m_{N}}{q^{2}+m_{\pi}^{2}}
!
! For the fitting of the contact term
!    Following [4] we can use the transition amplitude: A = < pp | V_{long} + V_{short} | nn > using the scattering wave function.
!    Don't forget to multiply (2\pi)^{3} for all terms and g_{A}^{2} for GT and T pieces.
!
! Tested with A. Belly
! To get the widely used dimensionless 0vbb NME, you need to multiply R = 1.2 A^(1/3) after all the calculation steps.
!
module LeptonNumberViolation
  !$ use omp_lib
  use StoreCouplings, only: CGsStore
  use PartialWaveDecomposition, only: PWD, &
      & identity_helicity_rep, sigma_dot_sigma, sigma1q_sigma2q
  implicit none

  public :: LNViolation
  private

  ! physics constants
  real(8), parameter, private :: lambda_V = 850.d0  ! MeV, scale used for the momentum dependence of g_V
  real(8), parameter, private :: lambda_A = 1086.d0 ! MeV, scale used the the momentum dependence of g_A
  type, extends(PWD) :: LNViolation
    real(8) :: Ec = 0.d0
    integer :: ch_order = 0
    character(:), allocatable :: OpName
    character(:), allocatable :: NNint
    real(8) :: c1, c2, c3, c4, gvv
  contains
    procedure :: InitLNViolation
    procedure :: FinLNViolation
    procedure :: show_lnv_options
    procedure :: calc_matrix_element
    !procedure :: set_helicity_rep
    generic :: init => InitLNViolation
    generic :: fin => FinLNViolation
  end type LNViolation

  type :: MomFunctions
    real(8), allocatable :: op(:,:,:,:)
    real(8), allocatable :: QMesh(:)
  end type MomFunctions
  character(32), private, parameter :: name_0vbb_f ="0vbbFermi"
  character(32), private, parameter :: name_0vbb_gt="0vbbGamowTeller"
  character(32), private, parameter :: name_0vbb_t ="0vbbTensor"
  character(32), private, parameter :: name_0vbb_ct="0vbbContact"
  character(32), private, parameter :: name_0veeC_f="0veeCapGamowFermi"
  character(32), private, parameter :: name_0veeC_gt0="0veeCapGamowTeller0"
  character(32), private, parameter :: name_0veeC_gt1="0veeCapGamowTeller1"
  character(32), private, parameter :: name_0veeC_gt2="0veeCapGamowTeller2"
  character(32), private, parameter :: name_0veeC_t="0veeCapTensor0"
  type(CGsStore), private :: cg_store

contains

  subroutine show_lnv_options(this)
    class(LNViolation), intent(in) :: this
    write(*,"(a)") "# Calculation options for Lepton-number violation operator"
    write(*,"(2a)") "# Operator: ", trim(this%OpName)
    write(*,"(a,f10.4,a)") "# Closure Energy: ", this%Ec, " MeV"
    write(*,"(a,i3)") "# Order of chiral expansion: ", this%ch_order
    write(*,"(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a)") "# c1= ", this%c1, " GeV-1, c2= ", &
        & this%c2, " GeV-1, c3=", this%c3, " GeV-1, c4= ", this%c4, " GeV-1"
    if(this%OpName==name_0vbb_ct) write(*,"(a,es12.4,a)") "# Overall contact strength: ", this%gvv, " MeV-3 fm-1"
    write(*,*)
    call this%show_pwd_options()
  end subroutine show_lnv_options

  subroutine FinLNViolation(this)
    class(LNViolation), intent(inout) :: this
    this%Ec = 0.d0
    this%ch_order = 0
    this%OpName = "None"
    this%c1 = 0.d0
    this%c2 = 0.d0
    this%c3 = 0.d0
    this%c4 = 0.d0
    this%gvv = 0.d0
    call cg_store%fin()
    call this%fin_pwd()
  end subroutine FinLNViolation

  subroutine InitLNViolation(this, OpName, NNint, Rank, Ec, ch_order, Jmax, Regulator, RegulatorPower, RegulatorLambda)
    use MyLibrary, only: pi, hc, m_nucleon, g_A, f_pi
    class(LNViolation), intent(inout) :: this
    character(*), intent(in) :: OpName, NNint
    real(8), intent(in), optional :: Ec, RegulatorLambda
    integer, intent(in), optional :: ch_order, Jmax, RegulatorPower, Rank
    real(8) :: gvv
    character(*), intent(in), optional :: Regulator
    this%OpName = OpName
    this%NNint = NNint

    call this%init_pwd(Jmax=Jmax, Rank=Rank, Regulator=Regulator, &
        & RegulatorPower=RegulatorPower, RegulatorLambda=RegulatorLambda)
    if( present(Ec) ) this%Ec = Ec
    if( present(ch_order) ) this%ch_order = ch_order
    call cg_store%init(2, 2, .false., 2, 2, .false.)

    gvv = 1.d0
    select case(this%NNint)
    case("N3LO_EM500")
      this%c1 = -0.81d0
      this%c2 =  2.8d0
      this%c3 = -3.2d0
      this%c4 =  5.4d0
      !if(this%ch_order==0) gvv = 0.597d0 ! see [4]
      !if(this%ch_order==2) gvv = 0.606d0 ! see [4]
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
      !if(this%ch_order==0) gvv = 0.755d0 ! see [4]
      !if(this%ch_order==2) gvv = 0.765d0 ! see [4]
    case("N4LO_EMN500")
      this%c1 = -1.10d0
      this%c2 =  3.57d0
      this%c3 = -5.54d0
      this%c4 =  4.17d0
    case default
      write(*,*) "Unknown NNint:", __FILE__, __LINE__
    end select
    gvv = gvv * 0.5d0 ! <p|tau|n>=sqrt(2)
    this%gvv = (m_nucleon*g_A / (4.d0*f_pi**2))**2 / hc
    this%gvv = this%gvv * 1.d0/(2.d0*pi)**3 ! Fourier trans convention
    this%gvv = this%gvv * 2.d0 * gvv * (-1.d0)
    this%gvv = this%gvv / (4.d0 * sqrt(2.d0))  ! I stole this factor from Antoine's implementation; 1/(2 sqrt(2)) should relate with the cm delta function convention
  end subroutine InitLNViolation

  function calc_matrix_element(this, pbra, lbra, sbra, jbra, zbra, pket, lket, sket, jket, zket, pn_formalism) result(r)
    class(LNViolation), intent(in) :: this
    real(8), intent(in) :: pbra, pket ! relative momenta of outgoing and incoming states
    integer, intent(in) :: lbra, sbra, jbra, zbra, lket, sket, jket, zket
    logical, intent(in) :: pn_formalism
    real(8) :: r, phbra, phket, norm
    integer :: ibra, iket, i
    integer, allocatable :: z1bras(:), z2bras(:), z1kets(:), z2kets(:)
    type(MomFunctions) :: fq

    allocate(fq%op(this%GetNMesh(), -this%GetRank():this%GetRank(), this%GetNumHeli(), this%GetNumHeli()))
    allocate(fq%QMesh(this%GetNMesh()))
    fq%op(:,:,:,:) = 0.d0
    do i = 1, this%GetNMesh()
      fq%QMesh(i) = sqrt(pbra**2 + pket**2 - 2.d0*pbra*pket*this%ZMesh(i))
    end do

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
      z1bras = [1,-1]
      z2bras = [-1,1]
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
      z1kets = [1,-1]
      z2kets = [-1,1]
    end if

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

    deallocate(fq%op)
    deallocate(fq%QMesh)
    deallocate(z1bras, z2bras, z1kets, z2kets)
  end function calc_matrix_element

  subroutine set_helicity_rep_pn(this, fq, pbra, pket, z1bra, z2bra, z1ket, z2ket)
    !
    !   returns the 0v double-beta decay matrix element in the unit of MeV-3 fm-1
    !
    class(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket ! relative momenta of outgoing and incoming states
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket
    integer :: i
    real(8) :: t_factor

    if(z1bra /= z2bra) return
    if(z1ket /= z2ket) return
    if(z1bra == z1ket) return
    t_factor = 2.d0
    select case(this%OpName)
    case(name_0vbb_f, name_0veeC_f)
      call set_helicity_rep_0vbb_fermi(this, fq, pbra, pket)
    case(name_0vbb_gt, name_0veeC_gt0, name_0veeC_gt1, name_0veeC_gt2)
      call set_helicity_rep_0vbb_gt(this, fq, pbra, pket)
    case(name_0vbb_t, name_0veeC_t)
      call set_helicity_rep_0vbb_tensor(this, fq, pbra, pket)
    case(name_0vbb_ct)
      call set_helicity_rep_0vbb_contact(this, fq, pbra, pket)
    case default
      write(*,*) "Unknown!", trim(this%OpName)
      stop
    end select
    fq%op(:,:,:,:) = fq%op(:,:,:,:) * t_factor
  end subroutine set_helicity_rep_pn

  subroutine set_helicity_rep_isospin(this, fq, pbra, pket, tbra, tket)
    !
    !   returns the 0v double-beta decay matrix element in the unit of MeV-3 fm-1
    !
    class(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket ! relative momenta of outgoing and incoming states
    integer, intent(in) :: tbra, tket ! -1 or 1 [pp and nn]
    integer :: i
    real(8) :: t_factor

    if(tbra /= 1 .or. tket /= 1) return
    t_factor = 2.d0 * sqrt(5.d0)
    select case(this%OpName)
    case(name_0vbb_f, name_0veeC_f)
      call set_helicity_rep_0vbb_fermi(this, fq, pbra, pket)
    case(name_0vbb_gt, name_0veeC_gt0, name_0veeC_gt1, name_0veeC_gt2)
      call set_helicity_rep_0vbb_gt(this, fq, pbra, pket)
    case(name_0vbb_t, name_0veeC_t)
      call set_helicity_rep_0vbb_tensor(this, fq, pbra, pket)
    case(name_0vbb_ct)
      call set_helicity_rep_0vbb_contact(this, fq, pbra, pket)
    case default
      write(*,*) "Unknown!", trim(this%OpName)
      stop
    end select
    fq%op(:,:,:,:) = fq%op(:,:,:,:) * t_factor
  end subroutine set_helicity_rep_isospin

  subroutine set_helicity_rep_0vbb_fermi(this, fq, pbra, pket)
    use MyLibrary, only: pi, hc
    type(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer :: lams_bra(2), lams_ket(2)
    integer :: i, ibra, iket
    real(8), allocatable :: v(:)
    allocate(v(this%GetNMesh()))
    call q_dependence_fermi()
    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = v(i) * identity_helicity_rep(lams_bra(1),lams_bra(2),lams_ket(1),lams_ket(2),i) / &
              & (2.d0 * pi**2 * hc)
        end do
      end do
    end do
    deallocate(v)
  contains
    subroutine q_dependence_fermi()
      !
      !  h_f(q) / q ( q + Ec)
      !
      real(8) :: q
      integer :: i

      do i = 1, this%GetNMesh()
        q = fq%QMesh(i)
        v(i) = neutrino_pot_fermi(q) / (q * (q+this%Ec))
      end do
    end subroutine q_dependence_fermi

    function neutrino_pot_fermi(q) result(hf)
      use MyLibrary, only: g_V
      real(8), intent(in) :: q
      real(8) :: hf
      hf = 1.d0
      if(this%ch_order<2)  return

      hf = gv_func(q)**2 / g_V**2
      if(this%ch_order>3) write(*,*) "Higher order is not implemented yet"
    end function neutrino_pot_fermi
  end subroutine set_helicity_rep_0vbb_fermi

  subroutine set_helicity_rep_0vbb_gt(this, fq, pbra, pket)
    use MyLibrary, only: pi, hc
    class(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer :: lams_bra(2), lams_ket(2)
    integer :: i, ibra, iket, mu
    real(8), allocatable :: v(:)
    allocate(v(this%GetNMesh()))
    call q_dependence_gt()
    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do mu = -this%GetRank(), this%GetRank()
          do i = 1, this%GetNMesh()
            fq%op(i,mu,ibra,iket) = v(i) * sigma1_sigma2(lams_bra(1),lams_bra(2),lams_ket(1),lams_ket(2),mu,i) / &
                & (2.d0 * pi**2 * hc)
          end do
        end do
      end do
    end do
    deallocate(v)
  contains

    function sigma1_sigma2(lam1_bra, lam2_bra, lam1_ket, lam2_ket, kappa, i) result(r)
      !
      ! [s1 s2]^rank_kappa
      !
      !use MyLibrary, only: dcg
      use PartialWaveDecomposition, only: helicity_expectation_value_sigma
      integer, intent(in) :: lam1_bra, lam2_bra, lam1_ket, lam2_ket, kappa, i
      real(8) :: r
      integer :: m1, m2
      r = 0.d0
      do m1 = -1, 1
        m2 = kappa - m1
        if(abs(m2) > 1) cycle
        r = r + cg_store%get(2, 2*m1, 2, 2*m2, 2*this%GetRank(), 2*kappa) * &
            & helicity_expectation_value_sigma( lam1_bra, lam1_ket, i, m1) * &
            & helicity_expectation_value_sigma(-lam2_bra,-lam2_ket, i, m2)
      end do
      !if(this%Rank==0) r = r * (-1.d0) * sqrt(3.d0) ! s1 . s2 = -sqrt(3) [s1 s2]^0
      if(this%OpName==name_0vbb_gt) r = r * (-1.d0) * sqrt(3.d0) ! s1 . s2 = -sqrt(3) [s1 s2]^0
    end function sigma1_sigma2

    subroutine q_dependence_gt()
      real(8) :: q
      integer :: i
      do i = 1, this%GetNMesh()
        q = fq%QMesh(i)
        v(i) = neutrino_pot_GT(q) / ( q * (q + this%Ec) )
      end do
    end subroutine q_dependence_gt

    function neutrino_pot_GT(q) result(hgt)
      use MyLibrary, only: g_A, m_nucleon
      real(8), intent(in) :: q
      real(8) :: hgt
      hgt = 1.d0
      if(this%ch_order<2) return
      hgt = ga_func(q)**2 + ga_func(q) * gp_func(q) * q**2 / (3.d0 * m_nucleon) + &
          & gp_func(q)**2 * q**4 / (12.d0 * m_nucleon**2) + gm_func(q)**2 * q**2 / (6.d0 * m_nucleon**2)
      hgt = hgt / g_A**2
      if(this%ch_order>3) write(*,*) "Higher order is not implemented yet"
    end function neutrino_pot_GT
  end subroutine set_helicity_rep_0vbb_gt

  subroutine set_helicity_rep_0vbb_tensor(this, fq, pbra, pket)
    use MyLibrary, only: pi, hc
    type(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer :: lams_bra(2), lams_ket(2)
    integer :: i, ibra, iket
    real(8) :: t_factor, z, f1, f2
    real(8), allocatable :: v(:)
    allocate(v(this%GetNMesh()))
    call q_dependence_tensor()
    f1 = 1.d0; f2 = 0.d0
    if(this%ch_order==2) then
      f1 = 3.d0; f2 = 1.d0
    end if
    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          t_factor = f1 * sigma1q_sigma2q(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), pbra, pket, i) / fq%QMesh(i)**2- &
              & f2 * sigma_dot_sigma(lams_bra(1), lams_bra(2), lams_ket(1), lams_ket(2), i)
          fq%op(i,0,ibra,iket) = v(i) * t_factor / (2.d0 * pi**2 * hc)
        end do
      end do
    end do
    deallocate(v)
  contains
    subroutine q_dependence_tensor()
      real(8) :: q
      integer :: i
      do i = 1, this%GetNMesh()
        q = fq%QMesh(i)
        v(i) = neutrino_pot_T(q) / (q * (q+this%Ec))
      end do
    end subroutine q_dependence_tensor

    function neutrino_pot_T(q) result(ht)
      use MyLibrary, only: g_A, m_pi, m_nucleon
      real(8), intent(in) :: q
      real(8) :: ht
      ht = - q**2 * (2.d0 * m_pi**2 + q**2) / (q**2 + m_pi**2)**2
      if(this%ch_order<2) return
      ht = - gp_func(q) * ga_func(q) * q**2 / (3.d0*m_nucleon) - &
          & gp_func(q)**2 * q**4 / (12.d0*m_nucleon**2) + &
          & gm_func(q)**2 * q**2 / (12.d0*m_nucleon**2)
      ht = ht / g_A**2 * (-1.d0) ! this -1 is from fourier trans r -> p space
      if(this%ch_order>3) write(*,*) "Higher order is not implemented yet"
    end function neutrino_pot_T
  end subroutine set_helicity_rep_0vbb_tensor

  subroutine set_helicity_rep_0vbb_contact(this, fq, pbra, pket)
    type(LNViolation), intent(in) :: this
    type(MomFunctions), intent(inout) :: fq
    real(8), intent(in) :: pbra, pket
    integer :: lams_bra(2), lams_ket(2)
    integer :: i, ibra, iket
    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        do i = 1, this%GetNMesh()
          fq%op(i,0,ibra,iket) = identity_helicity_rep(lams_bra(1),lams_bra(2),lams_ket(1),lams_ket(2),i) * this%gvv
        end do
      end do
    end do
  end subroutine set_helicity_rep_0vbb_contact

  function gv_func( q ) result(gv)
    use MyLibrary, only: g_V
    real(8), intent(in) :: q
    real(8) :: gv
    gv = g_V / (1.d0 + q**2 / Lambda_V**2)**2
  end function gv_func

  function ga_func( q ) result(ga)
    use MyLibrary, only: g_A
    real(8), intent(in) :: q
    real(8) :: ga
    ga = g_A / (1.d0 + q**2 / Lambda_A**2)**2
  end function ga_func

  function gm_func( q ) result(gm)
    use MyLibrary, only: kappa=>gv
    real(8), intent(in) :: q
    real(8) :: gm
    gm = kappa * gv_func(q)
  end function gm_func

  function gp_func( q ) result(gp)
    use MyLibrary, only: m_nucleon, m_pi
    real(8), intent(in) :: q
    real(8) :: gp
    gp = -(2.d0 * m_nucleon * ga_func(q) ) / ( q**2 + m_pi**2 )
  end function gp_func

end module LeptonNumberViolation
