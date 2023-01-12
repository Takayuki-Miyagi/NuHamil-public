module MyLibrary
  use, intrinsic :: iso_c_binding
  use ClassSys, only: sys
  implicit none

  real(8), parameter, public :: pi = acos(-1.d0)         ! \pi
  real(8), parameter, public :: hc = 197.32705d0         ! \hbar c [MeV fm] or [eV nm]
  real(8), parameter, public :: m_e = 510.9989461        ! electron mass [keV]
  real(8), parameter, public :: alpha = 137.035999d0     ! electric fine structure constant
  real(8), parameter, public :: m_proton = 938.27231d0   ! proton mass [MeV] 
  real(8), parameter, public :: m_neutron = 939.56563d0  ! neutron mass [MeV]
  real(8), parameter, public :: m_nucleon = 938.91897d0  ! MeV, nucleon mass
  real(8), parameter, public :: m_red_pn = 469.45926d0   ! reduced mass of pn
  real(8), parameter, public :: g_V = 1.0d0              ! vector coupling
  real(8), parameter, public :: g_A = 1.27d0             ! axial vector coupling
  real(8), parameter, public :: g_A_GT = 1.29d0          ! axial vector coupling, corrected Goldberger–Treiman relation: g_{pi NN} = g_A M_N / f_pi
  real(8), parameter, public :: f_pi = 92.2d0            ! pion decay constatnt [MeV]
  real(8), parameter, public :: m_pi = 138.038d0         ! MeV, pion mass
  real(8), parameter, public :: m_eta = 547.862d0        ! MeV, eta mass
  real(8), parameter, public :: m_rho = 775.45d0         ! MeV, rho mass
  real(8), parameter, public :: m_omega = 782.66d0       ! MeV, omega mass
  real(8), parameter, public :: g_pi  = 13.07d0          ! coupling NNpi
  real(8), parameter, public :: g_eta  = 2.24d0          ! coupling NNeta
  real(8), parameter, public :: g_rho = 2.75d0           ! coupling NNrho
  real(8), parameter, public :: g_omega = 8.25d0         ! coupling NNomega
  real(8), parameter, public :: lambda_chi = 700.d0      ! ChEFT break down scale [MeV]
  !real(8), parameter, public :: lambda_chi = 1158.61937d0 ! ChEFT break down scale [MeV] 4pi * f_pi
  real(8), parameter, public :: gs = 0.880               ! nucleon's magnetic moemnt g-factor isoscalar
  real(8), parameter, public :: gv = 4.706               ! nucleon's magnetic moment g-factor isovector  mu_{p/n} = (gs +/- gv)/2

  ! cache for Talmi-Moshinsky bracket
  integer, private, parameter :: n_trinomial = 100
  real(8), private, allocatable  :: dtrinomial(:,:,:)

  private :: dtrinomial_func

  !
  ! C interfaces
  !
  interface
    ! 3-j symbol
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

    ! factorial n!
    function factorial(n) bind(c,name='gsl_sf_fact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: factorial
    end function factorial

    ! double factorial n!!
    function double_factorial(n) bind(c,name='gsl_sf_doublefact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: double_factorial
    end function double_factorial

    ! Gamma function \Gamma(x)
    function gamma_function(x) bind(c,name='gsl_sf_gamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: gamma_function
    end function gamma_function

    ! log of Gamma function ln \Gamma(x)
    function ln_gamma(x) bind(c,name='gsl_sf_lngamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: ln_gamma
    end function ln_gamma

    ! spherical bessel function j_l(x)
    function spherical_bessel_c(l,x) bind(c,name='gsl_sf_bessel_jl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: spherical_bessel_c
    end function spherical_bessel_c

    ! Legendre polynomial P_l(x)
    function legendre_polynomial(l,x) bind(c,name='gsl_sf_legendre_Pl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: legendre_polynomial
    end function legendre_polynomial

    ! Associated Legendre polynomial sutable for pherical harmonics
    ! \sqrt{\frac{2l+1}{4\pi}} \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(x)
    function assoc_legendre_spharm(l,m,x) bind(c,name='gsl_sf_legendre_sphPlm')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l,m
      real(c_double), value, intent(in) :: x
      real(c_double) :: assoc_legendre_spharm
    end function assoc_legendre_spharm

    ! associated Laguerre polynomial L^{(a)}_{n}(x)
    function laguerre(n,a,x) bind(c,name='gsl_sf_laguerre_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: laguerre
    end function laguerre

    ! Gegenbauer polynomial C^{(lambda)}_{n}(x)
    function Gegenbauer_polynomial(n,lambda,x) bind(c,name='gsl_sf_gegenpoly_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: lambda, x
      real(c_double) :: gegenbauer_polynomial
    end function Gegenbauer_polynomial

    ! hypergeometric function 2F1(a,b,c,x)
    function hyperg_2F1(a, b, c, x) bind(c,name='gsl_sf_hyperg_2F1')
      import c_double
      real(c_double), value, intent(in) :: a, b, c, x
      real(c_double) :: hyperg_2F1
    end function hyperg_2F1

    ! Pochhammer symbol (a)_x
    function poch(a, x) bind(c,name='gsl_sf_poch')
      import c_double
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: poch
    end function poch

    ! Gauss-Legendre quadrature
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

    ! fixed-point quadrature
    function integration_fixed_allocate(T, n, a, b, alpha, beta) &
          & bind(c,name='gsl_integration_fixed_alloc')
      import c_int, c_double, c_ptr
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, b, alpha, beta
      type(c_ptr), value :: T
      type(c_ptr) :: integration_fixed_allocate
    end function integration_fixed_allocate
    subroutine integration_fixed_free( workspace ) bind(c, name="gsl_integration_fixed_free")
      import c_ptr
      type(c_ptr), value :: workspace
    end subroutine integration_fixed_free
    function integration_fixed_n( workspace ) bind(c, name="gsl_integration_fixed_n")
      import c_ptr, c_int
      type(c_ptr), value :: workspace
      integer(c_int) :: integration_fixed_n
    end function integration_fixed_n
    function integration_fixed_nodes( workspace ) bind(c, name="gsl_integration_fixed_nodes")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_nodes
    end function integration_fixed_nodes
    function integration_fixed_weights( workspace ) bind(c, name="gsl_integration_fixed_weights")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_weights
    end function integration_fixed_weights

    ! open, read, write, and close gzip file (additional interface gzip_open below)
    ! When you pass string to C, you need add NULL (achar(0)) in the end of strings.
    ! It is done in gzip_open.
    function gz_open(filename, mode) bind(c, name='gzopen')
      import c_char, c_ptr
      character(c_char) :: filename(*), mode(*)
      type(c_ptr) :: gz_open
    end function gz_open
    function gzip_read(f, buf, len) bind(c, name='gzgets')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_read
    end function gzip_read
    function gzip_write(f, buf, len) bind(c, name='gzwrite')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_write
    end function gzip_write
    function gzip_close( f ) bind(c, name='gzclose')
      import c_ptr
      type(c_ptr), value :: f
      type(c_ptr) :: gzip_close
    end function gzip_close
  end interface
  type(c_ptr) :: integration_fixed_legendre; bind(c, name="gsl_integration_fixed_legendre") :: integration_fixed_legendre
  type(c_ptr) :: integration_fixed_chebyshev; bind(c, name="gsl_integration_fixed_chebyshev") :: integration_fixed_chebyshev
  type(c_ptr) :: integration_fixed_gegenbauer; bind(c, name="gsl_integration_fixed_gegenbauer") :: integration_fixed_gegenbauer
  type(c_ptr) :: integration_fixed_jacobi; bind(c, name="gsl_integration_fixed_jacobi") :: integration_fixed_jacobi
  type(c_ptr) :: integration_fixed_laguerre; bind(c, name="gsl_integration_fixed_laguerre") :: integration_fixed_laguerre
  type(c_ptr) :: integration_fixed_hermite; bind(c, name="gsl_integration_fixed_hermite") :: integration_fixed_hermite
  type(c_ptr) :: integration_fixed_exponential; bind(c, name="gsl_integration_fixed_exponential") :: integration_fixed_exponential
  type(c_ptr) :: integration_fixed_rational; bind(c, name="gsl_integration_fixed_rational") :: integration_fixed_rational
  type(c_ptr) :: integration_fixed_chebyshev2; bind(c, name="gsl_integration_fixed_chebyshev2") :: integration_fixed_chebyshev2
  !
  ! end C interfaces
  !
contains

  subroutine skip_comment(nfile, comment)
    implicit none
    integer,intent(in)::nfile
    character(*), intent(in) :: comment
    type(sys) :: s
    character(20) :: line
    read(nfile,'(a)') line
    do while  (s%find(s%str(line), s%str(comment)))
      read(nfile,'(a)') line
    end do
    backspace(nfile)
  end subroutine skip_comment

  pure real(8) function delta(i1, i2)
    integer, intent(in) :: i1, i2
    delta = 0.d0
    if(i1 == i2) delta = 1.d0
  end function delta

  pure function hat(j) result(r)
    real(8) :: r
    integer, intent(in) :: j
    r = dsqrt(dble(j + 1))
  end function hat

  ! 3-j symbol for Wigner-Eckert theorem
  real(8) function geometry_part(jbra, jop, jket, mbra, mop, mket) result(r)
    integer, intent(in) :: jbra, mbra, jop, mop, jket, mket
    r = (-1.d0) ** ((jbra - mbra)/2) * &
      & tjs(jbra, jop, jket, -mbra, mop, mket)
  end function geometry_part

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

  pure real(8) function red_nab_l(n1, l1, n2, l2) result(nl)
    integer, intent(in) :: n1, l1, n2, l2
    if(n1 == n2 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*(dble(n2 + l2) + 1.5d0))
    elseif(n1 == n2-1 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif(n1 == n2 .and. l1 == l2-1) then
      nl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif(n1 == n2+1 .and. l1==l2-1) then
      nl = -dsqrt(dble(l2)*dble(n2 + 1))
    else
      nl = 0.d0
    end if
  end function red_nab_l

  pure logical function triag(i,j,k)
      implicit none
      integer,intent(in)::i,j,k
      triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: s
    s = coupling_3j(j1,j2,j3,m1,m2,-m3) * hat(j3) * (-1.d0) ** ((j1-j2+m3)/2)
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

  !!!
  ! For HO transformation bracket
  !!!
  function gmosh(nl, ll, nr, lr, n1, l1, n2, l2, lm, d) result(r)
    !
    ! reference G.P.Kamuntavicius, R.K.Kalinauskas, B.R.Barrett, S.Mickevicius, and D. Germanas, Nucl. Phys. A 695, 191 (2001).
    !
    real(8) :: r
    integer, intent(in) :: nl, ll, nr, lr, n1, l1, n2, l2, lm
    real(8), intent(in) :: d
    integer :: ee, er, e1, e2, m, ed, eb, ec, ea, ld, lb, lc, la
    real(8) :: t, s

    r = 0.d0
    if(.not. allocated(dtrinomial) ) then
      write(*,'(a)') "you need to call init_dtrinomial first!"
      return
    end if
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if(ee + er /= e1 + e2) return
    if(triag(ll, lr, lm)) return
    if(triag(l1, l2, lm)) return
    t = dsqrt((d ** (e1 - er)) / ((1.d0 + d) ** (e1 + e2)))
    m = min(er, e2)
    s = 1.d0
    do ed = 0, m
      eb = er - ed
      ec = e2 - ed
      ea = e1 - er + ed

      do ld = ed, 0, -2
        do lb = eb, 0, -2
          if(triag(ld,lb,lr)) cycle
          do lc = ec, 0, -2
            if(triag(ld,lc,l2)) cycle
            do la = ea, 0, -2
              if(triag(la,lb,l1)) cycle
              if(triag(la,ll,lc)) cycle

              r = r + s * t * &
                  & snj(2*la, 2*lb, 2*l1, 2*lc, 2*ld, 2*l2, 2*ll, 2*lr, 2*lm) * &
                  & g(e1, l1, ea, la, eb, lb) * g(e2, l2, ec, lc, ed, ld) * &
                  & g(ee, ll, ea, la, ec, lc) * g(er, lr, eb, lb, ed, ld)

            end do
          end do
        end do
      end do
      s = s * (-d)
    end do
    r = r * (-1.d0) ** (n1 + n2 + nr + nl)
  contains
    function g(e1, l1, ea, la, eb, lb) result(r)
      real(8) :: r
      integer, intent(in) :: e1, l1, ea, la, eb, lb

      r = dcg(2*la, 0, 2*lb, 0, 2*l1, 0) * dsqrt((2*la + 1) * (2*lb + 1) * &
          & dtrinomial(e1 - l1, ea - la, eb - lb) * &
          & dtrinomial(e1 + l1 + 1, ea + la + 1, eb + lb + 1))

    end function g
  end function gmosh

  subroutine init_dtrinomial()
    ! cache for Talmi-Moshinksy bracket
    integer :: i, j, k, n, m, info
    real(8) :: d
    allocate(dtrinomial(0:n_trinomial, 0:n_trinomial, 0:n_trinomial))
    !$omp parallel do private( i, j, k ) schedule (dynamic)
    do k = 0, n_trinomial
      do j = 0, n_trinomial
        do i = 0, n_trinomial
          dtrinomial(i, j, k) = dtrinomial_func(i, j, k)
        end do
      end do
    end do
  end subroutine init_dtrinomial

  subroutine fin_dtrinomial()
    deallocate(dtrinomial)
  end subroutine fin_dtrinomial


  function dtrinomial_func(i, j, k) result(s)
    !
    !  trinomial coefficient: i!! / j!! / k!!
    !
    integer, intent(in) :: i, j, k
    real(8) :: s
    integer :: m
    s = 1.d0
    m = max(i, j, k)
    if(m == 0) return
    if(m > n_trinomial) then

      write(*,'(a)') 'in trinomial_func, index is too large'
      return

    end if
    s = double_factorial(i) / (double_factorial(j) * double_factorial(k))
  end function dtrinomial_func
  !!!
  ! end  HO transformation bracket
  !!!

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

  function ho_radial_wf_norm(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l+1} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
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
    prefact = sqrt( 2.d0 * nu )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l+1) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf_norm

  function hydrogen_radial_wf(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* )**3 (n-l-1)!/2n(n+l)! ) x^{l} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n,l
    real(8),intent(in) :: a_star,r
    real(8) :: s
    real(8) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 4.d0/ a_star**3 ) / dble(n**2)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf

  function hydrogen_radial_wf_norm(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* ) (n-l-1)!/2n(n+l)! ) x^{l+1} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n,l
    real(8),intent(in) :: a_star,r
    real(8) :: s
    real(8) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 1.d0/a_star ) / dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l+1)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf_norm

  function Laguerre_radial_wf_norm(n,l,b,r) result(s)
    ! From Eq.(13) in A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,r
    real(8) :: s
    real(8) :: prefact, exp_component, x, x2
    x = r / b
    x2 = 2.d0 * x
    prefact = sqrt( 2.d0 / b )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - &
        & 0.5d0*ln_gamma(dble(n+2*l+3)) - x  + dble(l+1) * log(x2)
    s = prefact * exp(exp_component) * laguerre(n,dble(2*l+2),x2)
  end function Laguerre_radial_wf_norm

  function hydrogen_radial_wf_mom_norm(n,l,a_star,p) result(s)
    ! from wikipedia https://en.wikipedia.org/wiki/Hydrogen_atom
    ! Do not forget phase (-i)^{l} when you calclate the matrix elements!!
    ! x = p * (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n, l
    real(8),intent(in) :: a_star, p
    real(8) :: s
    real(8) :: prefact, exp_component, x, z
    x = p * dble(n) * a_star
    z = (x**2-1.d0) / (x**2+1.d0)
    prefact = sqrt( 2.d0/pi * a_star) * dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) + dble(2*l+2)*log(2.d0) + &
        & ln_gamma(dble(l+1)) + dble(l+1)*log(x) - dble(l+2)*log(x**2+1)
    s = prefact * exp(exp_component) * Gegenbauer_polynomial(n-l-1,dble(l+1),z)
  end function hydrogen_radial_wf_mom_norm

  function spherical_harmonics(l,m,cth,phi) result(r)
    !
    ! Y_{l}^{m}(\theta,\phi) = (-1)^{(m-|m|}/2} \sqrt{\frac{2l+1}{4\pi}}
    !            \times \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(\cos(\theta)) e^{i m \phi}
    ! cth = \cos(\theta)
    !
    integer, intent(in) :: l, m
    real(8), intent(in) :: cth, phi
    complex(8) :: ex, ei=(0.d0,1.d0)
    complex(8) :: r
    ex = cos(dble(m)*phi) + ei * sin(dble(m)*phi)
    r = (-1.d0)**( (m-abs(m))/2 ) * assoc_legendre_spharm(l,abs(m),cth) * ex ! Condon-Shortley phase
  end function spherical_harmonics

  function spherical_harmonics0(l,m,cth) result(r)
    !
    ! Y_{l}^{m}(\theta,\phi) = (-1)^{(m-|m|}/2} \sqrt{\frac{2l+1}{4\pi}}
    !            \times \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(\cos(\theta))
    ! cth = \cos(\theta)
    !
    integer, intent(in) :: l, m
    real(8), intent(in) :: cth
    real(8) :: r
    r = (-1.d0)**( (m-abs(m))/2 ) * assoc_legendre_spharm(l,abs(m),cth)
  end function spherical_harmonics0

  function spherical_harmonics_unnorm(l,m,cth,phi) result(r)
    !
    ! Y_{l}^{m}(\theta,\phi) = (-1)^{(m+|m|}/2} \sqrt{\frac{2l+1}{4\pi}}
    !            \times \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(\cos(\theta)) e^{i m \phi}
    ! cth = \cos(\theta)
    !
    integer, intent(in) :: l, m
    real(8), intent(in) :: cth, phi
    complex(8) :: r
    r = spherical_harmonics(l,m,cth,phi) * sqrt( 4.d0*pi / (dble(2*l+1)) )
  end function spherical_harmonics_unnorm

  subroutine gauss_legendre(x1,x2,x,w,n)
    ! input:
    ! x1   : lower limit of the integration interval
    ! x2   : upper limit ---------- "" -------------
    ! n    : the desired number of mesh points
    ! output :
    ! x     : gauss-legendre mesh points on the interval (x1,x2)
    ! w     : the corresponding weights
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

  subroutine gauss_legendre_(x1,x2,x,w,n)
    ! input:
    ! x1   : lower limit of the integration interval
    ! x2   : upper limit ---------- "" -------------
    ! n    : the desired number of mesh points
    ! output :
    ! x     : gauss-legendre mesh points on the interval (x1,x2)
    ! w     : the corresponding weights
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out), allocatable :: x(:), w(:)
    call fixed_point_quadrature("legendre", n, x, w, a_in=x1, b_in=x2)
  end subroutine gauss_legendre_

  subroutine fixed_point_quadrature(quad_name, n, x, w, a_in, b_in, alpha_in, beta_in, weight_renorm)
    character(*), intent(in) :: quad_name
    integer, intent(in) :: n
    real(8), intent(in), optional :: a_in, b_in, alpha_in, beta_in
    logical, intent(in), optional :: weight_renorm
    real(8), allocatable :: x(:), w(:)
    real(8) :: a=0.d0, b=0.d0, alpha=0.d0, beta=0.d0
    integer :: i
    type(c_ptr) :: workspace, nodes, weights
    real(c_double), pointer :: x_(:), w_(:)

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    if(present(a_in)) a = a_in
    if(present(b_in)) b = b_in
    if(present(alpha_in)) alpha = alpha_in
    if(present(beta_in)) beta = beta_in

    select case(quad_name)
    case("legendre","Legendre")
      workspace = integration_fixed_allocate(integration_fixed_legendre, n, a, b, alpha, beta)
    case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev, n, a, b, alpha, beta)
    case("gegenbauer","Gegenbauer")
      workspace = integration_fixed_allocate(integration_fixed_gegenbauer, n, a, b, alpha, beta)
    case("jacobi","Jacobi")
      workspace = integration_fixed_allocate(integration_fixed_jacobi, n, a, b, alpha, beta)
    case("laguerre","Laguerre")
      workspace = integration_fixed_allocate(integration_fixed_laguerre, n, a, b, alpha, beta)
    case("hermite","Hermite")
      workspace = integration_fixed_allocate(integration_fixed_hermite, n, a, b, alpha, beta)
    case("exponential","Exponential")
      workspace = integration_fixed_allocate(integration_fixed_exponential, n, a, b, alpha, beta)
    case("rational","Rational")
      workspace = integration_fixed_allocate(integration_fixed_rational, n, a, b, alpha, beta)
    case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev2, n, a, b, alpha, beta)
    case default
      write(*,"(a)") "Unknown quadrature name"
      stop
    end select

    nodes = integration_fixed_nodes( workspace )
    weights = integration_fixed_weights( workspace )
    call c_f_pointer(nodes, x_, [n])
    call c_f_pointer(weights, w_, [n])
    x(:) = x_(:)
    w(:) = w_(:)
    call integration_fixed_free(workspace)
    if(.not. present(weight_renorm)) return
    if(.not. weight_renorm) return
    do i = 1, n
      select case(quad_name)
      case("legendre","Legendre")
        w(i) = w(i) * 1.d0
      case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
        w(i) = w(i) * sqrt((b-x(i)) * (x(i)-a))
      case("gegenbauer","Gegenbauer")
        w(i) = w(i) / ( (b-x(i)) * (x(i)-a) )**alpha
      case("jacobi","Jacobi")
        w(i) = w(i) / ( (b-x(i))**alpha * (x(i)-a)**beta )
      case("laguerre","Laguerre")
        w(i) = w(i) * exp( b * (x(i)-a) ) / (x(i)-a)**alpha
      case("hermite","Hermite")
        w(i) = w(i) * exp( b * (x(i)-a)**2 ) / abs(x(i)-a)**alpha
      case("exponential","Exponential")
        w(i) = w(i) / abs( x(i)-(a+b)*0.5d0 )**alpha
      case("rational","Rational")
        w(i) = w(i) / ( (x(i)-a)**alpha * (x(i)+b)*beta )
      case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
        w(i) = w(i) / sqrt( (b-x(i)) * (x(i)-a) )
      case default
        write(*,"(a)") "Unknown quadrature name"
        stop
      end select
    end do
  end subroutine fixed_point_quadrature

  function gzip_open( filename, mode ) result(p)
    character(*), intent(in) :: filename, mode
    type(c_ptr) :: p
    p = gz_open(trim(filename)//achar(0), trim(mode)//achar(0))
  end function gzip_open

  function gzip_writeline( f, buf, len) result(p)
    type(c_ptr) :: f, p
    character(*), intent(in) :: buf
    integer, intent(in) :: len
    p = gzip_write( f, trim(buf)//achar(10), len+1)
  end function gzip_writeline

  function gzip_readline( f, buf, len) result(p)
    ! note
    ! len_trim returns length removed space (32 in ascii code)
    ! -2 means removing null (0) and line feed (10)
    type(c_ptr) :: f, p
    character(*), intent(inout) :: buf
    integer, intent(in) :: len
    buf=""
    p = gzip_read( f, buf, len)
    buf = buf(1:len_trim(buf) - 2)
  end function gzip_readline

  function spherical_bessel(l,x) result(r)
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: r, a = 0.d0
    r = 0.d0
    ! threashhold that jl(x) ~ 10^-200 for small x, avoiding the under flow error
    a = exp(-200.d0 / dble(l) * log(10.d0) + dble(2*l+1)/dble(l) * log(dble(2*l+1)) - &
      & dble(2*l+1)/dble(l) - log(dble(l)) + 1 - log(2.d0) )
    if(x < a) return
    r = spherical_bessel_c(l,x)
  end function spherical_bessel

  !
  ! isospin function pn formalism
  !
  function tau_1(zbra,zket) result(r)
    ! <zbra | zket >
    integer, intent(in) :: zbra, zket
    real(8) :: r
    r = 0.d0
    if(zbra /= zket) r = 0.d0
    if(zbra == zket) r = 1.d0
  end function tau_1

  function tau_x(zbra,zket) result(r)
    ! < zbra| tau_{x} | zket >
    integer, intent(in) :: zbra, zket
    real(8) :: r
    r = 0.d0
    if(zbra==zket) then
      r = 0.d0
    else if(zbra/=zket) then
      r = 1.d0
    end if
  end function tau_x

  function tau_y(zbra,zket,phase) result(r)
    ! < zbra| tau_{y} | zket >, i is not taken into account
    integer, intent(in) :: zbra, zket
    integer, intent(in), optional :: phase
    real(8) :: r, ph
    r = 0.d0
    ph = 1.d0
    if(present(phase)) ph = dble(phase)
    if(zbra==zket) then
      r = 0.d0
    else if(zbra==1 .and. zket==-1) then
      r = -1.d0 * ph
    else if(zbra==-1 .and. zket==1) then
      r = 1.d0 * ph
    end if
  end function tau_y

  function tau_z(zbra,zket,phase) result(r)
    ! < zbra| tau_{z} | zket >
    integer, intent(in) :: zbra, zket
    integer, intent(in), optional :: phase
    real(8) :: r, ph
    r = 0.d0
    ph = 1.d0
    if(present(phase)) ph = dble(phase)
    if(zbra/=zket) then
      r = 0.d0
    else if(zbra== 1) then
      r = 1.d0 * ph
    else if(zbra==-1) then
      r =-1.d0 * ph
    end if
  end function tau_z

  function tau_m(zbra,zket,m,phase) result(r)
    ! < zbra| tau_{-1,0,1} | zket >
    integer, intent(in) :: zbra, zket, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = 0.d0
    if(m==1) then
      r = -(tau_x(zbra,zket) - tau_y(zbra,zket,phase)) / sqrt(2.d0)
    else if(m== 0) then
      r = tau_z(zbra,zket,phase)
    else if(m==-1) then
      r = ( tau_x(zbra,zket) + tau_y(zbra,zket,phase)) / sqrt(2.d0)
    end if
  end function tau_m

  function tau1_tau2_tensor(z1bra,z2bra,z1ket,z2ket,rank,m,phase) result(r)
    !
    ! <z1bra, z2bra | [tau_1 tau_2]^rank_m | z1ket, z2ket >
    !
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    integer :: m1, m2
    real(8) :: r
    r = 0.d0
    do m1 = -1, 1
      m2 = m-m1
      if(abs(m2)>1) cycle
      r = r + dcg(2, 2*m1, 2, 2*m2, 2*rank, 2*m) * tau_m(z1bra,z1ket,m1,phase) * tau_m(z2bra,z2ket,m2,phase)
    end do
  end function tau1_tau2_tensor

  function tau1_dot_tau2(z1bra,z2bra,z1ket,z2ket,rank,m,phase) result(r)
    !
    ! <z1bra, z2bra | tau_1x tau_2x + tau_1y tau_2y + tau_1z tau_2z | z1ket, z2ket >
    ! m has to be zero
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = -sqrt(3.d0) * tau1_tau2_tensor(z1bra,z2bra,z1ket,z2ket,0,0,phase=phase)
  end function tau1_dot_tau2

  function tau1_cross_tau2(z1bra,z2bra,z1ket,z2ket,rank,m,phase) result(r)
    !
    ! <z1bra, z2bra | (tau1 x tau2)_m | z1ket, z2ket > = -i sqrt(2) <z1bra, z2bra | [tau_1 tau_2]^1_m | z1ket, z2ket >
    ! note: i is not taken into account here
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = -sqrt(2.d0) * tau1_tau2_tensor(z1bra,z2bra,z1ket,z2ket,1,m,phase=phase)
  end function tau1_cross_tau2

  function tau1_m(z1bra, z2bra, z1ket, z2ket, rank, m, phase) result(r)
    !
    ! <z1bra, z2bra | tau1_m | z1ket, z2ket >
    !
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = tau_m(z1bra,z1ket,m,phase)*tau_1(z2bra,z2ket)
  end function tau1_m

  function tau2_m(z1bra, z2bra, z1ket, z2ket, rank, m, phase) result(r)
    !
    ! <z1bra, z2bra | tau1_m | z1ket, z2ket >
    !
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = tau_1(z1bra,z1ket)*tau_m(z2bra,z2ket,m,phase)
  end function tau2_m

  function tau1_minus_tau2(z1bra,z2bra,z1ket,z2ket,rank,m,phase) result(r)
    !
    ! <z1bra, z2bra | (tau1 - tau2)_m | z1ket, z2ket >
    !
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = tau1_m(z1bra,z2bra,z1ket,z2ket,1,m,phase) - tau2_m(z1bra,z2bra,z1ket,z2ket,1,m,phase)
  end function tau1_minus_tau2

  function tau1_plus_tau2(z1bra,z2bra,z1ket,z2ket,rank,m,phase) result(r)
    !
    ! <z1bra, z2bra | (tau1 + tau2)_m | z1ket, z2ket >
    !
    integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
    integer, intent(in), optional :: phase
    real(8) :: r
    r = tau1_m(z1bra,z2bra,z1ket,z2ket,1,m,phase) + tau2_m(z1bra,z2bra,z1ket,z2ket,1,m,phase)
  end function tau1_plus_tau2

  function asym_isospin_func_pn(lbra,sbra,zbra,lket,sket,zket,func,rank,m,phase) result(r)
    interface
      function func(z1bra, z2bra, z1ket, z2ket, rank, m, phase) result(r)
        integer, intent(in) :: z1bra, z2bra, z1ket, z2ket, rank, m
        integer, intent(in), optional :: phase
        real(8) :: r
      end function func
    end interface
    integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket, rank, m
    integer, intent(in), optional :: phase
    integer, allocatable :: z1bra(:), z2bra(:), z1ket(:), z2ket(:)
    real(8) :: r
    integer :: ibra, iket
    real(8) :: norm, ph_bra, ph_ket
    call pn_combinations()
    norm = 1.d0 / sqrt(dble(size(z1bra)*size(z1ket)))
    r = 0.d0
    do ibra = 1, size(z1bra)
      ph_bra = 1.d0
      if(z1bra(ibra)== 1 .and. z2bra(ibra)==-1) ph_bra = (-1.d0)**(lbra+sbra)
      do iket = 1, size(z1ket)
        ph_ket = 1.d0
        if(z1ket(iket)== 1 .and. z2ket(iket)==-1) ph_ket = (-1.d0)**(lket+sket)
        r = r + func(z1bra(ibra),z2bra(ibra),z1ket(iket),z2ket(iket),rank,m,phase) * ph_bra*ph_ket
      end do
    end do
    r = r*norm
    deallocate(z1bra,z2bra,z1ket,z2ket)
  contains
    subroutine pn_combinations()
      if(zbra==-1) then
        allocate(z1bra(1), z2bra(1))
        z1bra = [-1]
        z2bra = [-1]
      else if(zbra==1) then
        allocate(z1bra(1), z2bra(1))
        z1bra = [1]
        z2bra = [1]
      else if(zbra==0) then
        allocate(z1bra(2), z2bra(2))
        z1bra = [-1,1]
        z2bra = [1,-1]
      else
        write(*,*) "Error:", __LINE__, __FILE__
        stop
      end if

      if(zket==-1) then
        allocate(z1ket(1), z2ket(1))
        z1ket = [-1]
        z2ket = [-1]
      else if(zket==1) then
        allocate(z1ket(1), z2ket(1))
        z1ket = [1]
        z2ket = [1]
      else if(zket==0) then
        allocate(z1ket(2), z2ket(2))
        z1ket = [-1,1]
        z2ket = [1,-1]
      else
        write(*,*) "Error:", __LINE__, __FILE__
        stop
      end if
    end subroutine pn_combinations

  end function asym_isospin_func_pn

  !
  ! isospin function isospin formalism
  !
  function tau_identity_iso(tbra,tket) result(r)
    ! <tbra || 1 ||tket >
    integer, intent(in) :: tbra, tket
    real(8) :: r
    if(tbra /= tket) r = 0.d0
    if(tbra == tket) r = sqrt(dble(2*tket+1))
  end function tau_identity_iso

  function tau1_tau2_tensor_iso(tbra,tket,rank) result(r)
    !
    ! <1/2 1/2: tbra || [tau1 tau2]^rank || 1/2 1/2: tket >
    !
    integer, intent(in) :: tbra, tket, rank
    real(8) :: r
    r = 6.d0 * sqrt( dble( (2*rank+1) * (2*tbra+1) * (2*tket+1)) ) * snj(1,1,2*tbra,1,1,2*tket,2,2,2*rank)
  end function tau1_tau2_tensor_iso

  function tau1_iso(tbra,tket) result(r)
    !
    ! <1/2 1/2: tbra || tau1 || 1/2 1/2: tket >
    !
    integer, intent(in) :: tbra, tket
    real(8) :: r
    r = sqrt(dble( 6 * (2*tbra+1) * (2*tket+1) )) * (-1.d0)**tket * sjs(1, 1, 2, 2*tket, 2*tbra, 1)
  end function tau1_iso

  function tau2_iso(tbra,tket) result(r)
    !
    ! <1/2 1/2: tbra || tau2 || 1/2 1/2: tket >
    !
    integer, intent(in) :: tbra, tket
    real(8) :: r
    r = sqrt(dble( 6 * (2*tbra+1) * (2*tket+1) )) * (-1.d0)**tbra * sjs(1, 1, 2, 2*tket, 2*tbra, 1)
  end function tau2_iso

  function tau1_minus_tau2_iso(tbra,tket) result(r)
    !
    ! <1/2 1/2: tbra || (tau1 - tau2) || 1/2 1/2: tket >
    !
    integer, intent(in) :: tbra, tket
    real(8) :: r
    r = tau1_iso(tbra,tket) - tau2_iso(tbra,tket)
  end function tau1_minus_tau2_iso

  function tau1_plus_tau2_iso(tbra,tket) result(r)
    !
    ! <1/2 1/2: tbra || (tau1 + tau2) || 1/2 1/2: tket >
    !
    integer, intent(in) :: tbra, tket
    real(8) :: r
    r = tau1_iso(tbra,tket) + tau2_iso(tbra,tket)
  end function tau1_plus_tau2_iso

end module MyLibrary

