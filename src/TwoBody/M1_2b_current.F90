!
!
! This is just an ad hoc code. Only chiral NLO (LO) two-body current
! In the future, I want to replace it with more general one
!
module M1_2b_current
  use, intrinsic :: iso_c_binding
  use MyLibrary, only: pi, hc, m_proton, m_neutron, m_pi, f_pi, g_A, &
      & ho_radial_wf, tjs, sjs, snj, dcg, triag, tau1_tau2_tensor_iso, gauss_legendre
  implicit none
  public :: init_M1_2b_current
  public :: fin_M1_2b_current
  public :: me_M1_2b_current
  public :: me_M1_2b_current_sachs
  public :: me_M1_2b_current_intr
  private

  integer, private :: NMesh = 200
  real(8), private :: rmax = 50.d0
  real(8), private, allocatable :: r_mesh(:), w_mesh(:)
  real(8), private, allocatable :: radial_wf_rel(:,:,:)
  real(8), private :: hw
  real(8), private :: mu_at_rrel = -9999.d0
  real(8), private :: mu_at_Rcm = -9999.d0
  integer, private :: Nmax, Jmax
  logical, private :: initialized = .false.

contains
  subroutine init_M1_2b_current(hbar_omega, N_max, J_max, N_Mesh, r_max, op_at_rrel, op_at_Rcm)
    real(8), intent(in) :: hbar_omega
    integer, intent(in) :: N_max, J_max
    integer, intent(in), optional :: N_Mesh
    real(8), intent(in), optional :: r_max
    real(8), intent(in), optional :: op_at_rrel
    real(8), intent(in), optional :: op_at_Rcm
    real(8) :: a, v
    integer :: i, l, n
    if(initialized) return
    hw = hbar_omega
    Nmax = N_max
    Jmax = J_max
    if(present(N_Mesh)) NMesh = N_Mesh
    if(present(r_max)) rmax = r_max
    if(present(op_at_rrel)) mu_at_rrel = op_at_rrel
    if(present(op_at_Rcm)) mu_at_Rcm = op_at_Rcm
    call gauss_legendre(0.d0, rmax, r_mesh, w_mesh, NMesh)
    allocate(radial_wf_rel(NMesh,0:Nmax/2,0:Jmax+1))
    radial_wf_rel(:,:,:) = 0.d0
    a = 0.25d0 * (m_proton + m_neutron) * hw / hc**2 ! fm-2
    do n = 0, Nmax/2
      do l = 0, Jmax+1
        if(2*n + l > Nmax) cycle
        do i = 1, NMesh
          radial_wf_rel(i,n,l) = ho_radial_wf(n,l,a,r_mesh(i))
        end do
      end do
    end do
    initialized=.true.
  end subroutine init_M1_2b_current

  subroutine fin_M1_2b_current()
    deallocate(radial_wf_rel)
    deallocate(r_mesh)
    deallocate(w_mesh)
    initialized = .false.
  end subroutine fin_M1_2b_current

  function me_M1_2b_current(nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra, &
        & nket, lket, sket, jrket, ncmket, lcmket, jket, zket, pn) result(me)
    integer, intent(in) :: nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra
    integer, intent(in) :: nket, lket, sket, jrket, ncmket, lcmket, jket, zket
    logical, intent(in) :: pn
    real(8) :: me

    me = me_M1_2b_current_sachs(nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra, &
        & nket, lket, sket, jrket, ncmket, lcmket, jket, zket, pn)
    if(ncmbra/=ncmket) return
    if(lcmbra/=lcmket) return
    me = me + (-1.d0)**(Lcmket+jrket+Jbra+1) * sjs(2*jrbra, 2*Jbra, 2*Lcmket, 2*Jket, 2*jrket, 2) * &
        & sqrt(dble((2*Jbra+1)*(2*Jket+1))) * &
        & me_M1_2b_current_intr(nbra, lbra, sbra, jrbra, zbra, nket, lket, sket, jrket, zket, pn)
  end function me_M1_2b_current

  function me_M1_2b_current_intr(nbra, lbra, sbra, jbra, zbra, nket, lket, sket, jket, zket, pn) result(me)
    integer, intent(in) :: nbra, lbra, sbra, jbra, zbra
    integer, intent(in) :: nket, lket, sket, jket, zket
    logical, intent(in) :: pn
    integer :: k, i
    real(8) :: a, me, x, prefact, v
    real(8) :: tau_part, angular_part1, angular_part2
    me = 0.d0
    if(pn) tau_part = -tau1_tau2_vector_0_pn(lbra, sbra, zbra, lket, sket, zket)
    if(.not. pn) tau_part = -tau1_tau2_tensor_iso(zbra, zket, 1) * (-sqrt(2.d0))
    if(abs(tau_part) < 1.d-8) return
    prefact = -g_A**2 * m_pi / (8.d0 * pi * (2.d0*f_pi)**2)
    prefact = prefact * (-1.d0) * (2.d0 * m_proton) ! (-1) is from i^2, MeV-1 -> mu_N
    prefact = prefact * sqrt(3.d0 / (4.d0 * pi))

    angular_part1 = 0.d0 ! [(s1 x s2) . r] r
    do k = 0, 2, 2
      if(triag(lbra,lket,k)) cycle
      angular_part1 = angular_part1 + sqrt(8.d0 * pi / 3.d0) * dcg(2, 0, 2, 0, 2*K, 0) * &
          & Y_sigma_lsj(lbra, sbra, jbra, lket, sket, jket, K, 1, 1)
    end do
    angular_part2 = sqrt(2.d0) * sigma_cross_sigma_lsj(lbra, sbra, jbra, lket, sket, jket, 1) ! s1 x s2 term

    ! r space -> HO
    if(mu_at_rrel < 0.d0) then
      do i = 1, NMesh
      x = m_pi * r_mesh(i) / hc
      v = prefact * ((1.d0 + 1.d0/x) * angular_part1 + angular_part2) * exp(-x) * tau_part
      me = me + &
        & w_mesh(i) * r_mesh(i)**2 * radial_wf_rel(i,nbra,lbra) * v * radial_wf_rel(i,nket,lket)
      end do
    else
      a = 0.25d0 * (m_proton + m_neutron) * hw / hc**2 ! fm-2
      x = m_pi * mu_at_rrel / hc
      v = prefact * ((1.d0 + 1.d0/x) * angular_part1 + angular_part2) * exp(-x) * tau_part
      me = mu_at_rrel**2 * ho_radial_wf(nbra,lbra,a,mu_at_rrel) * ho_radial_wf(nket,lket,a,mu_at_rrel) * v
    end if
  end function me_M1_2b_current_intr

  function me_M1_2b_current_sachs(nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra, &
        & nket, lket, sket, jrket, ncmket, lcmket, jket, zket, pn) result(me)
    integer, intent(in) :: nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra
    integer, intent(in) :: nket, lket, sket, jrket, ncmket, lcmket, jket, zket
    logical, intent(in) :: pn
    integer :: i
    real(8) :: me, v, a
    real(8) :: tau_part, b_cm, cm_part, ang_part1, ang_part2, x, prefact
    integer :: alpha, beta
    me = 0.d0
    if(pn) tau_part = -tau1_tau2_vector_0_pn(lbra, sbra, zbra, lket, sket, zket)
    if(.not. pn) tau_part = -tau1_tau2_tensor_iso(zbra, zket, 1) * (-sqrt(2.d0))
    if(abs(tau_part) < 1.d-8) return
    b_cm = sqrt( hc**2 / ( (m_proton+m_neutron)*hw ))
    cm_part = red_r_l(ncmbra, lcmbra, ncmket, lcmket) * b_cm
    if(abs(cm_part) < 1.d-10) return

    prefact = g_A**2 * m_pi**2 / (8.d0*pi*(2.d0*f_pi)**2)
    prefact = prefact * 2.d0 * m_proton / hc * (-1.d0) ! (-1) is from i^2
    prefact = prefact * sqrt(3.d0 / (4.d0 * pi))
    ang_part1 = 0.d0
    do alpha = 0, 2, 2
      do beta = abs(1-alpha), 1+alpha
        ang_part1 = ang_part1 + sqrt(8.d0 * pi / 3.d0) * dcg(2, 0, 2, 0, 2*alpha, 0) * dcg(2, 0, 2*alpha, 0, 2*beta, 0) * &
            & R_Y_sigma_lsLj(ncmbra, lcmbra, lbra, sbra, jrbra, jbra, &
            & ncmket, lcmket, lket, sket, jrket, jket, beta, alpha, 1, 1) * b_cm
      end do
    end do
    ang_part2 = sqrt(8.d0 * pi) * R_Y_sigma_lsLj(ncmbra, lcmbra, lbra, sbra, jrbra, jbra, &
        & ncmket, lcmket, lket, sket, jrket, jket, 1, 0, 1, 1) * b_cm

    !
    ! r space -> HO
    if(mu_at_rrel<0.d0) then
      me = 0.d0
      do i = 1, NMesh
        x = m_pi * r_mesh(i) / hc
        v = (1.d0 + 3.d0 / x + 3.d0 / x**2) * ang_part1 + (1.d0 / x + 1.d0 / x**2) * ang_part2
        v = v * prefact * exp(-x)
        me = me + r_mesh(i)**2 * w_mesh(i) * radial_wf_rel(i,nbra,lbra) * v * radial_wf_rel(i,nket,lket)
      end do
    else
      a = 0.25d0 * (m_proton + m_neutron) * hw / hc**2 ! fm-2
      x = m_pi * mu_at_rrel / hc
      v = (1.d0 + 3.d0 / x + 3.d0 / x**2) * ang_part1 + (1.d0 / x + 1.d0 / x**2) * ang_part2
      v = v * prefact * exp(-x)
      me = mu_at_rrel**2 * ho_radial_wf(nbra,lbra,a,mu_at_rrel) * ho_radial_wf(nket,lket,a,mu_at_rrel) * v
    end if
    me = me * tau_part
  end function me_M1_2b_current_sachs

  function R_Y_sigma_lsLj(ncmbra, lcmbra, lbra, sbra, jrbra, jbra, &
        & ncmket, lcmket, lket, sket, jrket, jket, rellrank, srank, relrank, rank) result(r)
    !
    ! <Ncm'Lcm' l's'j': J' || [R [Y^k [s1 s2]^srank]^relrank]^rank || Ncm Lcm lsj : J >
    !
    integer, intent(in) :: ncmbra, lcmbra, lbra, sbra, jrbra, jbra
    integer, intent(in) :: ncmket, lcmket, lket, sket, jrket, jket
    integer, intent(in) :: rellrank, srank, relrank, rank
    real(8) :: r
    r = j_to_ls(lcmbra, jrbra, jbra, lcmket, jrket, jket, 1, relrank, rank) * &
        & red_r_l(ncmbra, lcmbra, ncmket, lcmket) * &
        & Y_sigma_lsj(lbra, sbra, jrbra, lket, sket, jrket, rellrank, srank, relrank)
  end function R_Y_sigma_lsLj

  function Y_sigma_lsj(lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank) result(r)
    !
    ! < l's'j' || [Y^lrank [s1 s2]^srank]^jrank || lsj >
    !
    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank
    real(8) :: r
    r = j_to_ls(lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank) * &
        & reduced_Yk(lbra, lket, lrank) * sigma_sigma_tensor(sbra, sket, srank)
  end function Y_sigma_lsj

  function sigma_cross_sigma_lsj(lbra, sbra, jbra, lket, sket, jket, srank) result(r)
    !
    ! < l's'j' || [s1 s2]^srank || lsj >
    !
    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, srank
    real(8) :: r
    r = 0.d0
    if(lbra /= lket) return
    r = j_to_ls(lbra, sbra, jbra, lket, sket, jket, 0, srank, srank) * &
        & sqrt(dble(2*lket+1)) * sigma_sigma_tensor(sbra, sket, srank)
  end function sigma_cross_sigma_lsj

  function sigma_sigma_tensor(sbra, sket, srank) result(r)
    !
    ! < s' || [s1 s2]^srank || s >
    !
    integer, intent(in) :: sbra, sket, srank
    real(8) :: r
    r = 6.d0 * sqrt(dble((2*srank+1)*(2*sbra+1)*(2*sket+1))) * &
        & snj(1, 1, 2*sbra, 1, 1, 2*sket, 2, 2, 2*srank)
  end function sigma_sigma_tensor

  function reduced_Yk(lbra, lket, lrank) result(r)
    !
    ! < l' || Y^lrank || l >
    !
    integer, intent(in) :: lbra, lket, lrank
    real(8) :: r
    r = ((-1.d0)**lbra) * sqrt( dble((2*lrank+1)*(2*lbra+1)*(2*lket+1)) / (4.d0 * pi) ) * &
        & tjs(2*lbra, 2*lrank, 2*lket, 0, 0, 0)
  end function reduced_Yk

  function j_to_ls(lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank) result(r)
    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank
    real(8) :: r
    r = sqrt(dble((2*jrank+1)*(2*jbra+1)*(2*jket+1))) * &
        & snj(2*lbra, 2*sbra, 2*jbra, 2*lket, 2*sket, 2*jket, 2*lrank, 2*srank, 2*jrank)
  end function j_to_ls

  real(8) function red_r_l(n1, l1, n2, l2) result(rl)
    integer, intent(in) :: n1, l1, n2, l2
    real(8) :: a
    if(mu_at_Rcm<0.d0) then
      if (n1 == n2 .and. l1 == l2-1) then
        rl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
      elseif (n1 == n2-1 .and. l1 == l2+1) then
        rl = -dsqrt(dble(l2 + 1)*dble(n2))
      elseif (n1 == n2+1 .and. l1 == l2-1) then
        rl = dsqrt(dble(l2)*(dble(n2  + 1)))
      elseif (n1 == n2 .and. l1==l2+1) then
        rl = dsqrt(dble(l2+1)*(dble(n2 +l2)+1.5d0))
      else
        rl = 0.d0
      end if
    else
      if(abs(l1-l2) == 1) then
        a = (m_proton + m_neutron) * hw / hc**2 ! fm-2
        rl = (-1.d0)**l1 * sqrt(dble(2*l1+1)*dble(2*l2+1)) * tjs(2*l1, 2, 2*l2, 0, 0, 0) * &
          & mu_at_Rcm**3 * ho_radial_wf(n1,l1,a,mu_at_Rcm) * ho_radial_wf(n2,l2,a,mu_at_Rcm) * sqrt(a)
      else
        rl = 0.d0
      end if
    end if
  end function red_r_l

  function tau1_tau2_vector_0_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
    ! tau_{1,x}tau_{2,y} - tau_{1,y}tau_{2,x}
    ! i is not taken into account here
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: r
    r = 0.d0
    if(zbra /= zket) return
    if(abs(zket) == 1) return
    r = (-1.d0)**(lbra+sbra) - (-1.d0)**(lket+sket)
  end function tau1_tau2_vector_0_pn
end module M1_2b_current
