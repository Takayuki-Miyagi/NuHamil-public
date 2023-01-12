!
! This is just an ad hoc code. Only chiral NLO (LO) two-body current
! In the future, I want to replace it with more general one
! Here some functions are empty. I will update it once I publish a paper.
!
module M1_2b_current
  use, intrinsic :: iso_c_binding
  use MyLibrary, only: pi, hc, m_proton, m_neutron, m_pi, f_pi, g_A
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
  integer, private :: Nmax, Jmax
  logical, private :: initialized = .false.

  interface
    function gauss_legendre_allocate(n) &
          & bind(c,name='gsl_integration_glfixed_table_alloc')
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: n
      type(c_ptr) :: gauss_legendre_allocate
    end function gauss_legendre_allocate
    function gauss_legendre_ith_point_weight(a,b,i,xi,wi,t) &
          & bind(c,name='gsl_integration_glfixed_point')
      import c_int, c_double, c_ptr
      real(c_double), value, intent(in) :: a, b
      integer(c_int), value, intent(in) :: i
      real(c_double) :: xi, wi
      type(c_ptr), value, intent(in) :: t
      integer(c_int) :: gauss_legendre_ith_point_weight
    end function gauss_legendre_ith_point_weight
    subroutine gauss_legendre_release(t) &
          & bind(c,name='gsl_integration_glfixed_table_free')
      import c_ptr
      type(c_ptr), value :: t
    end subroutine gauss_legendre_release

    function coupling_3j(j1,j2,j3,m1,m2,m3) bind(c,name='gsl_sf_coupling_3j')
      import c_int, c_double
      real(c_double) :: coupling_3j
      integer(c_int), value, intent(in) :: j1,j2,j3,m1,m2,m3
    end function coupling_3j

    ! 6-j symbol
    function coupling_6j(j1,j2,j3,j4,j5,j6) bind(c,name='gsl_sf_coupling_6j')
      import c_int, c_double
      real(c_double) :: coupling_6j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6
    end function coupling_6j

    ! 9-j symbol
    function coupling_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) bind(c,name='gsl_sf_coupling_9j')
      import c_int, c_double
      real(c_double) :: coupling_9j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    end function coupling_9j

    ! log of Gamma function ln \Gamma(x)
    function ln_gamma(x) bind(c,name='gsl_sf_lngamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: ln_gamma
    end function ln_gamma

    ! associated Laguerre polynomial L^{(a)}_{n}(x)
    function laguerre(n,a,x) bind(c,name='gsl_sf_laguerre_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: laguerre
    end function laguerre

  end interface

contains
  subroutine init_M1_2b_current(hbar_omega, N_max, J_max, N_Mesh, r_max)
    real(8), intent(in) :: hbar_omega
    integer, intent(in) :: N_max, J_max
    integer, intent(in), optional :: N_Mesh
    real(8), intent(in), optional :: r_max
    real(8) :: a, v
    integer :: i, l, n
    if(initialized) return
    hw = hbar_omega
    Nmax = N_max
    Jmax = J_max
    if(present(N_Mesh)) NMesh = N_Mesh
    if(present(r_max)) rmax = r_max
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

    me = 0.d0
    me = me_M1_2b_current_sachs(nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra, &
        & nket, lket, sket, jrket, ncmket, lcmket, jket, zket, pn)
    if(ncmbra==ncmket .and. lcmbra==lcmket) then
      me = me + (-1)**(Lcmket+jrket+Jbra+1) * sjs(2*jrbra, 2*Jbra, 2*Lcmket, 2*Jket, 2*jrket, 2) * &
          & me_M1_2b_current_intr(nbra, lbra, sbra, jrbra, zbra, nket, lket, sket, jrket, zket, pn)
    end if
  end function me_M1_2b_current

  function me_M1_2b_current_intr(nbra, lbra, sbra, jbra, zbra, nket, lket, sket, jket, zket, pn) result(me)
    integer, intent(in) :: nbra, lbra, sbra, jbra, zbra
    integer, intent(in) :: nket, lket, sket, jket, zket
    logical, intent(in) :: pn
    real(8) :: me
    me = 0.d0
  end function me_M1_2b_current_intr

  function me_M1_2b_current_sachs(nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra, &
        & nket, lket, sket, jrket, ncmket, lcmket, jket, zket, pn) result(me)
    integer, intent(in) :: nbra, lbra, sbra, jrbra, ncmbra, lcmbra, jbra, zbra
    integer, intent(in) :: nket, lket, sket, jrket, ncmket, lcmket, jket, zket
    logical, intent(in) :: pn
    real(8) :: me
    me = 0.d0
  end function me_M1_2b_current_sachs

  subroutine gauss_legendre(x1,x2,x,w,n)
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out), allocatable :: x(:), w(:)
    real(8) :: xi, wi
    integer :: info, i
    type(c_ptr) :: t

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    t = gauss_legendre_allocate(n)
    do i = 1, n
      info = gauss_legendre_ith_point_weight(x1,x2,i-1,xi,wi,t)
      x(i) = xi
      w(i) = wi
    end do
    call gauss_legendre_release(t)
  end subroutine gauss_legendre

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: s
    s = coupling_3j(j1,j2,j3,m1,m2,-m3) * sqrt(dble(j3+1)) * (-1.d0) ** ((j1-j2+m3)/2)
  end function dcg

  function tjs(j1, j2, j3, m1, m2, m3) result(r)
    real(8) :: r
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    r = coupling_3j(j1,j2,j3,m1,m2,m3)
  end function tjs

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    real(8) :: s
    s = coupling_6j(j1,j2,j3,l1,l2,l3)
  end function sjs

  function snj(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
    !
    !  9j coefficient
    !  d9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
    !    {(j11)/2 (j12)/2 (j13)/2}
    !  = {(j21)/2 (j22)/2 (j23)/2}
    !    {(j31)/2 (j32)/2 (j33)/2}
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    real(8) :: s
    s = coupling_9j(j11,j12,j13,j21,j22,j23,j31,j32,j33)
  end function snj

  function ho_radial_wf(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu**3 * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer,intent(in) :: n,l
    real(8),intent(in) :: anu,r
    real(8) :: s
    real(8) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu**3 )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf

  function tau1_minus_tau2_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
    ! tau_{1,z} - tau_{2,z}
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
    real(8) :: r
    r = 0.d0
    if(zbra /= zket) return
    if(abs(zket) == 1) return
    r = -1.d0 + (-1.d0)**(lbra+sbra+lket+sket)
  end function tau1_minus_tau2_pn

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

  function tau1_tau2_tensor_iso(tbra, tket, Trank) result(r)
    !  < Tbra || [tau1 x tau2]^{Trank} || Tket >
    use MyLibrary, only: snj
    integer, intent(in) :: tbra, tket, Trank
    real(8) :: r
    r = 6.d0 * sqrt(dble((2*Tbra+1)*(2*Tket+1)*(2*Trank+1))) * &
        & snj(1, 1, 2*Tbra, 1, 1, 2*Tket, 2, 2, 2*Trank)
  end function tau1_tau2_tensor_iso

  function sigma1(Sbra,Sket) result(r)
    ! < Sbra || sigma_{1} || Sket >
    integer, intent(in) :: Sbra, Sket
    real(8) :: r
    r = sqrt(6.d0) * (-1.d0)**Sket * &
        & sqrt(dble( (2*Sbra+1) * (2*Sket+1) ) ) * &
        & sjs(1, 1, 2, 2*Sket, 2*Sbra, 1)
  end function sigma1

  function sigma2(Sbra,Sket) result(r)
    ! < Sbra || sigma_{2} || Sket >
    integer, intent(in) :: Sbra, Sket
    real(8) :: r
    r = sqrt(6.d0) * (-1.d0)**Sbra * &
        & sqrt(dble( (2*Sbra+1) * (2*Sket+1) ) ) * &
        & sjs(1, 1, 2, 2*Sket, 2*Sbra, 1)
  end function sigma2

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

  pure real(8) function red_r_l(n1, l1, n2, l2) result(rl)
    integer, intent(in) :: n1, l1, n2, l2
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
  end function red_r_l

  pure logical function triag(i,j,k)
      implicit none
      integer,intent(in)::i,j,k
      triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag
end module M1_2b_current
