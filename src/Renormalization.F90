module Renormalization
  use omp_lib
  use LinAlgLib
  use ClassSys
  use Profiler, only: timer
  implicit none

  public :: LeeSuzukiSolver
  public :: SRGSolver

  private :: Inils, finls, LeeSuzuki
  private :: inisrg, finsrg, SRGHFlow, SRGHUFlow, SRGOmegaFlow

  type :: LeeSuzukiSolver
    ! H: transformed Hamiltonian
    type(DMat) :: S, H, expS
  contains
    procedure :: init => inils
    procedure :: fin => finls
    procedure :: LeeSuzuki
  end type LeeSuzukiSolver

  type :: SRGSolver
    ! H: evolved hamiltonian
    type(DMat) :: U, H, Omega
    character(20) :: mode = 'Hflow' ! default
    real(8) :: atol = 1.d-8
    real(8) :: rtol = 1.d-8
  contains
    procedure :: init => inisrg
    procedure :: fin => finsrg
    procedure :: SRGFlow
    procedure :: SRGHFlow
    procedure :: SRGHUFlow
    procedure :: SRGOmegaFlow
  end type SRGSolver

  type(sys), private :: sy
contains
  subroutine inisrg(this, Hin, flow_mode, atol, rtol)
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin
    character(*), intent(in) :: flow_mode
    real(8), intent(in), optional :: atol, rtol
    integer(4) :: n
    n = size(Hin%m, 1)
    call this%omega%ini(n,n)
    call this%U%eye(n)
    call this%H%ini(n,n)
    this%H = Hin
    this%mode = flow_mode
    if(present(atol)) this%atol = atol
    if(present(rtol)) this%rtol = rtol
  end subroutine inisrg

  subroutine finsrg(this)
    class(SRGSolver), intent(inout) :: this
    call this%omega%fin()
    call this%u%fin()
    call this%H%fin()
  end subroutine finsrg

  subroutine SRGFlow(this, Hin, Gen, alpha, rhs_zero_index)
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin, Gen
    real(8), intent(inout) :: alpha
    integer(4), intent(in), optional :: rhs_zero_index(:)

    select case(this%mode)
    case("HFlow", "Hflow")
      call SRGHflow(this, Hin, Gen, alpha)
    case("ModifiedHFlow", "ModifiedHflow")
      call ModifiedSRGHflow(this, Hin, Gen, alpha, rhs_zero_index)
    case("OmegaFlow", "Omegaflow")
      call SRGOmegaflow(this, Hin, Gen, alpha)
    case("HUFlow", "HUflow")
      call SRGHUflow(this, Hin, Gen, alpha)
    case default
      call SRGHflow(this, Hin, Gen, alpha)
    end select
  end subroutine SRGFlow

  subroutine SRGHflow(this, Hin, Gen, alpha)
    use dvode_f90_m, only: set_opts, get_stats, release_arrays, vode_opts, dvode_f90
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin, Gen
    type(DMat) :: eta, rhs
    type(EigenSolSymD) :: sbare, seff
    real(8), intent(inout) :: alpha
    real(8), allocatable :: vec_evolution(:)
    real(8) :: start_point, prod
    integer :: itask, istate, n_ode, i, j
    integer(4) :: n
    real(8) :: rstats(22)
    integer :: istats(31)
    type(vode_opts) :: options
    real(8) :: ti

    ti = omp_get_wtime()

    n = size(Hin%m, 1)
    call eta%ini(n,n); call rhs%ini(n,n)
    n_ode = n * (n + 1) / 2
    allocate(vec_evolution(n_ode))
    this%H = Hin
    start_point = 0.d0
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        vec_evolution(i * (i - 1) / 2 + j) = hin%m(i,j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    itask = 1
    istate = 1
    options = set_opts(method_flag=10, abserr=this%atol, relerr=this%rtol, mxstep=50000)
    !options = set_normal_opts(dense_j = .true., abserr_vector = atol, relerr = relerr)
    !options = set_opts(sparse_j=.true., abserr=this%atol, relerr=this%rtol, mxstep=100000, nzswag=20000)
    call dvode_f90(h_equation, n_ode, vec_evolution, start_point, alpha, itask, istate, options)
    call get_stats(rstats, istats)
    call release_arrays()
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        this%h%m(i,j) = vec_evolution(i * (i - 1) / 2 + j)
        this%h%m(j,i) = vec_evolution(i * (i - 1) / 2 + j)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(vec_evolution)
    call eta%fin(); call rhs%fin()

    call sbare%init(hin)
    call sbare%DiagSym(hin)
    call seff%init(this%h)
    call seff%DiagSym(this%h)
    do i = 1, n
      prod = dot_product(sbare%vec%m(:,i),seff%vec%m(:,i))
      if(prod < 0.d0) then
        ! assuming that the transformation does not much chagne the eigen vectors.
        seff%vec%m(:,i) = seff%vec%m(:,i) * (-1.d0)
      end if
    end do
    this%u = sbare%vec * seff%vec%T()

    call sbare%fin()
    call seff%fin()

    call timer%add(sy%str('SRG H flow'), omp_get_wtime() - ti)
  contains
    subroutine h_equation(neq, t, v, dv)
      integer, intent(in) :: neq
      real(8), intent(in) :: t
      real(8), intent(in) :: v(neq)
      real(8), intent(out) :: dv(neq)

      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          this%h%m(i,j) = v(i * (i - 1) / 2 + j)
          this%h%m(j,i) = v(i * (i - 1) / 2 + j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      !eta = (Gen * this%h) - (this%h * Gen)
      !rhs = (eta * this%h) - (this%h * eta)
      call dgemm('n','n', n, n, n, 1.d0, Gen%m, n, this%h%m, n, 0.d0, eta%m, n)
      call dgemm('n','n', n, n, n,-1.d0, this%h%m, n, Gen%m, n, 1.d0, eta%m, n)
      call dgemm('n','n', n, n, n, 1.d0, eta%m, n, this%h%m, n, 0.d0, rhs%m, n)
      call dgemm('n','n', n, n, n,-1.d0, this%h%m, n, eta%m, n, 1.d0, rhs%m, n)

      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          dv(i * (i - 1) / 2 + j) = rhs%m(i,j)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine h_equation
  end subroutine SRGHflow

  subroutine ModifiedSRGHflow(this, Hin, Gen, alpha, rhs_zero_index)
    use dvode_f90_m, only: set_opts, get_stats, release_arrays, vode_opts, dvode_f90
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin, Gen
    integer(4), intent(in), optional :: rhs_zero_index(:)
    type(DMat) :: eta, rhs, correction
    type(EigenSolSymD) :: sbare, seff
    real(8), intent(inout) :: alpha
    real(8), allocatable :: vec_evolution(:)
    real(8) :: start_point, prod
    integer :: itask, istate, n_ode, i, j
    integer(4) :: n
    real(8) :: rstats(22)
    integer :: istats(31)
    type(vode_opts) :: options
    real(8) :: ti

    ti = omp_get_wtime()

    n = size(Hin%m, 1)
    call eta%ini(n,n); call rhs%ini(n,n); call correction%zeros(n,n)
    n_ode = n * (n + 1) / 2
    allocate(vec_evolution(n_ode))
    this%H = Hin
    start_point = 0.d0
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        vec_evolution(i * (i - 1) / 2 + j) = hin%m(i,j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    itask = 1
    istate = 1
    options = set_opts(method_flag=10, abserr=this%atol, relerr=this%rtol, mxstep=50000)
    !options = set_normal_opts(dense_j = .true., abserr_vector = atol, relerr = relerr)
    !options = set_opts(sparse_j=.true., abserr=this%atol, relerr=this%rtol, mxstep=100000, nzswag=20000)
    call dvode_f90(h_equation, n_ode, vec_evolution, start_point, alpha, itask, istate, options)
    call get_stats(rstats, istats)
    call release_arrays()
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        this%h%m(i,j) = vec_evolution(i * (i - 1) / 2 + j)
        this%h%m(j,i) = vec_evolution(i * (i - 1) / 2 + j)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(vec_evolution)
    call eta%fin(); call rhs%fin(); call correction%fin()

    call sbare%init(hin)
    call sbare%DiagSym(hin)
    call seff%init(this%h)
    call seff%DiagSym(this%h)
    do i = 1, n
      prod = dot_product(sbare%vec%m(:,i),seff%vec%m(:,i))
      if(prod < 0.d0) then
        ! assuming that the transformation does not much chagne the eigen vectors.
        seff%vec%m(:,i) = seff%vec%m(:,i) * (-1.d0)
      end if
    end do
    this%u = sbare%vec * seff%vec%T()

    call sbare%fin()
    call seff%fin()

    call timer%add(sy%str('SRG H flow'), omp_get_wtime() - ti)
  contains
    subroutine h_equation(neq, t, v, dv)
      integer, intent(in) :: neq
      real(8), intent(in) :: t
      real(8), intent(in) :: v(neq)
      real(8), intent(out) :: dv(neq)

      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          this%h%m(i,j) = v(i * (i - 1) / 2 + j)
          this%h%m(j,i) = v(i * (i - 1) / 2 + j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      eta = (Gen * this%h) - (this%h * Gen)
      if( present(rhs_zero_index) ) then
        correction%m(:,:) = 0.d0
        do i = 1, size(rhs_zero_index)
          correction%m(rhs_zero_index(i),:) = this%h%m(rhs_zero_index(i),:) - gen%m(rhs_zero_index(i),:)
          correction%m(:,rhs_zero_index(i)) = this%h%m(:,rhs_zero_index(i)) - gen%m(:,rhs_zero_index(i))
        end do
      end if
      rhs = (eta * (this%h-correction)) - ((this%h-correction) * eta)
      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          dv(i * (i - 1) / 2 + j) = rhs%m(i,j)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine h_equation
  end subroutine ModifiedSRGHflow

  subroutine SRGHUflow(this, Hin, Gen, alpha)
    use dvode_f90_m, only: set_opts, get_stats, release_arrays, vode_opts, dvode_f90
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin, Gen
    type(DMat) :: eta, rhs
    real(8), intent(inout) :: alpha
    real(8), allocatable :: vec_evolution(:), atol(:)
    real(8) :: start_point, ti
    integer :: itask, istate, n_ode, i, j
    integer(4) :: n
    real(8) :: rstats(22)
    integer :: istats(31)
    type(vode_opts) :: options

    ti = omp_get_wtime()

    n = size(Hin%m, 1)
    call eta%ini(n,n); call rhs%ini(n,n)
    n_ode = n * (n + 1) / 2 + n ** 2
    allocate(vec_evolution(n_ode), atol(n_ode))
    atol = this%atol
    itask = 1
    istate = 1
    this%H = Hin
    start_point = 0.d0

    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        vec_evolution(i * (i - 1) / 2 + j) = hin%m(i,j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, n
        vec_evolution(n * (n+1)/2 + n * (i - 1)  + j) = this%u%m(i,j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    !options = set_normal_opts(dense_j = .true., abserr_vector = atol, relerr = relerr)
    options = set_opts(method_flag=10, abserr_vector = atol, &
        & relerr = this%rtol, mxstep=50000)
    call dvode_f90(hu_equation, n_ode, vec_evolution, start_point, alpha, itask, istate, options)
    call get_stats(rstats, istats)
    call release_arrays
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        this%h%m(i,j) = vec_evolution(i * (i - 1) / 2 + j)
        this%h%m(j,i) = vec_evolution(i * (i - 1) / 2 + j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, n
        this%u%m(i,j) = vec_evolution(n*(n+1)/2 + n*(i-1) + j)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(vec_evolution, atol)
    call eta%fin(); call rhs%fin()

    call timer%add(sy%str('SRG H and U flow'), omp_get_wtime() - ti)
  contains
    subroutine hu_equation(neq, t, v, dv)
      integer, intent(in) :: neq
      real(8), intent(in) :: t
      real(8), intent(in) :: v(neq)
      real(8), intent(out) :: dv(neq)
      type(DMat) :: u, du

      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          this%h%m(i,j) = v(i * (i - 1) / 2 + j)
          this%h%m(j,i) = v(i * (i - 1) / 2 + j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      u = this%u
      do i = 1, n
        do j = 1, n
          u%m(i,j) = v(n*(n+1)/2 + n*(i-1) + j)
        end do
      end do

      !eta = (Gen * this%h) - (this%h * Gen)
      !rhs = (eta * this%h) - (this%h * eta)
      call dgemm('n','n', n, n, n, 1.d0, Gen%m, n, this%h%m, n, 0.d0, eta%m, n)
      call dgemm('n','n', n, n, n,-1.d0, this%h%m, n, Gen%m, n, 1.d0, eta%m, n)
      call dgemm('n','n', n, n, n, 1.d0, eta%m, n, this%h%m, n, 0.d0, rhs%m, n)
      call dgemm('n','n', n, n, n,-1.d0, this%h%m, n, eta%m, n, 1.d0, rhs%m, n)
      du = (u * eta) * (-1.d0)
      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          dv(i * (i - 1) / 2 + j) = rhs%m(i,j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      do i = 1, n
        do j = 1, n
          dv(n*(n+1)/2 + n*(i-1)+j) = du%m(i,j)
        end do
      end do
      call u%fin(); call du%fin()
    end subroutine hu_equation
  end subroutine SRGHUflow

  subroutine SRGOmegaflow(this, Hin, Gen, alpha)
    use dvode_f90_m, only: set_opts, get_stats, release_arrays, vode_opts, dvode_f90
    class(SRGSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin, Gen
    type(DMat) :: eta, rhs, ad
    real(8), intent(inout) :: alpha
    real(8), allocatable :: vec_evolution(:), atol(:)
    real(8) :: start_point, ti
    integer :: itask, istate, n_ode, i, j
    integer(4) :: n
    real(8) :: rstats(22)
    integer :: istats(31)
    type(vode_opts) :: options
    integer, parameter :: nb = 14
    real(8) :: bn(0:nb)

    ti = omp_get_wtime()

    n = size(Hin%m, 1)
    call eta%ini(n,n); call rhs%ini(n,n)
    n_ode = n * (n + 1) / 2
    allocate(vec_evolution(n_ode), atol(n_ode))
    atol = this%atol
    itask = 1
    istate = 1
    this%H = Hin
    start_point = 0.d0
    vec_evolution = 0.d0

    call init_bernoulli()
    !options = set_normal_opts(dense_j = .true., abserr_vector = atol, relerr = relerr)
    options = set_opts(method_flag=10, abserr_vector = atol, &
        & relerr = this%rtol, mxstep=50000)
    call dvode_f90(omg_equation, n_ode, vec_evolution, start_point, alpha, itask, istate, options)
    call get_stats(rstats, istats)
    call release_arrays
    !$omp parallel
    !$omp do private(i, j)
    do i = 1, n
      do j = 1, i
        this%omega%m(i,j) =   vec_evolution(i * (i - 1) / 2 + j)
        this%omega%m(j,i) = - vec_evolution(i * (i - 1) / 2 + j)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(vec_evolution, atol)
    call eta%fin(); call rhs%fin()

    call timer%add(sy%str('SRG Omega flow'), omp_get_wtime() - ti)
  contains
    subroutine omg_equation(neq, t, v, dv)
      integer, intent(in) :: neq
      real(8), intent(in) :: t
      real(8), intent(in) :: v(neq)
      real(8), intent(out) :: dv(neq)
      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          this%omega%m(i,j) =  v(i * (i - 1) / 2 + j)
          this%omega%m(j,i) = -v(i * (i - 1) / 2 + j)
        end do
      end do
      !$omp end do
      !$omp end parallel

      ad = exp((-1.d0) * this%omega)
      this%h = ad%T() * Hin * ad
      eta = (Gen * this%h) - (this%h * Gen)
      ad = eta
      rhs = eta
      do i = 1, nb
        ad = ((this%omega * ad) - (ad * this%omega)) / dble(i)
        rhs = rhs + bn(i) * ad
        if(mod(i, 2) == 0 .and. maxval(abs(bn(i) * ad%m)) < 1.d-8) exit
        if(i == nb) then
          write(*,'(a)') "Warning, Magnus expansion may not converge!"
        end if
      end do

      !$omp parallel
      !$omp do private(i, j)
      do i = 1, n
        do j = 1, i
          dv(i * (i - 1) / 2 + j) = rhs%m(i,j)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine omg_equation

    subroutine init_bernoulli
      integer :: k, n
      bn(:) = 0.d0
      bn(0) = 1.d0
      bn(1) = -0.5d0
      do n = 2, nb - 1
        if(mod(n, 2) == 0) then
          do k = 0, n - 1
            bn(n) = bn(n) + binomial(n+1, k) * bn(k)
          end do
          bn(n) = - bn(n) / dble(n + 1)
        end if
      end do
    end subroutine init_bernoulli

    real(8) function binomial(n, k) result(b)
      integer, intent(in) :: n, k
      integer :: i
      real(8) :: nume, denm1, denm2

      nume = 0.d0; denm1 = 0.d0; denm2 = 0.d0
      if(n > 0) then
        do i = 1, n
          nume = nume + log(dble(i))
        end do
      end if

      if(k > 0) then
        do i = 1, k
          denm1 = denm1 + log(dble(i))
        end do
      end if

      if(n-k > 0) then
        do i = 1, n-k
          denm2 = denm2 + log(dble(i))
        end do
      end if

      b = exp(nume - denm1 - denm2)
    end function binomial

  end subroutine SRGOmegaflow

  subroutine inils(this, Hin)
    class(LeeSuzukiSolver), intent(inout) :: this
    type(DMat), intent(in) :: Hin
    integer(4) :: n
    n = size(Hin%m, 1)
    call this%S%ini(n,n); call this%H%ini(n,n)
  end subroutine inils

  subroutine finls(this)
    class(LeeSuzukiSolver), intent(inout) :: this
    call this%S%fin(); call this%H%fin()
  end subroutine finls

  subroutine LeeSuzuki(this, ndim_ex, npdim, nqdim, h)
    class(LeeSuzukiSolver), intent(inout) :: this
    integer(4), intent(in) :: ndim_ex, npdim, nqdim
    type(DMat), intent(in) :: h
    type(EigenSolSymD) :: sol
    type(DMat) :: hwork, omg
    integer(4) :: n, info
    real(8) :: ti

    ti = omp_get_wtime()

    n = ndim_ex + npdim + nqdim
    if(npdim < 1 .or. nqdim < 1) return
    n = npdim + nqdim
    call hwork%ini(n,n)
    hwork%m(:,:) = h%m(ndim_ex+1:, ndim_ex+1:)
    call sol%init(hwork)
    call sol%DiagSym(hwork)
    omg = Omega_LeeSuzuki(sol%vec, n, npdim, nqdim, info)
    call hwork%fin()
    if(info .ne. 0) return
    call sol%fin()
    this%S = Omega2S(ndim_ex, npdim, nqdim, omg)
    call omg%fin()
    this%expS = exp(this%S)
    this%H = this%expS%T() * h * this%expS
    call timer%add(sy%str('Lee-Suzuki mthod'), omp_get_wtime() - ti)
  contains
    type(DMat) function Omega_LeeSuzuki(wave, n, npdim, nqdim, info) result(omg)
      type(DMat), intent(in) :: wave
      integer(4), intent(in) :: n, npdim, nqdim
      integer(4), intent(out) :: info
      real(8), allocatable :: povlap(:)
      integer, allocatable :: lgt(:)
      real(8) :: ovlap
      type(DMat) :: a, ainv, ww_qp
      integer(4) :: i, j, k
      info = 0
      call omg%ini(nqdim, npdim)
      allocate(povlap(n), lgt(npdim))
      call a%Ini(npdim, npdim); call ww_qp%ini(nqdim, npdim)
      do i = 1, n
        ovlap = dot_product(wave%m(:npdim, i), wave%m(:npdim, i))
        povlap(i) = ovlap
      end do

      do i = 1, npdim
        k = 1
        do j = 1, n
          if(povlap(k) < povlap(j)) k = j
        end do
        lgt(i) = k
        povlap(k) = 0.d0
      end do

      do i = 1, npdim
        k = lgt(i)
        do j = 1, nqdim
          ww_qp%m(j,i) = wave%m(j+npdim, k)
        end do
      end do

      do i = 1, npdim
        do j = 1, npdim
          k = lgt(j)
          a%m(i,j) = wave%m(i,k)
        end do
      end do

      if(abs(a%Det()) < 1.d-8) then
        info = 1
        write(*,'(a, es18.6)') "in Omega_LeeSuzuki, linear dependence is not preserved: Det = ", a%Det()
        deallocate(povlap, lgt)
        call a%fin(); call ww_qp%fin()
        return
      end if
      ainv = a%Inv()
      omg = ww_qp * ainv
      deallocate(povlap, lgt)
    end function Omega_LeeSuzuki

    type(DMat) function Omega2S(nex, np, nq, omg) result(S)
      integer(4), intent(in) :: nex, np, nq
      type(DMat), intent(in) :: omg
      type(DMat) :: a, theta
      type(DVec) :: frac
      type(EigenSolSymD) :: sol
      integer(4) :: i, j

      a = omg%T() * omg
      call S%ini(nex + np + nq, nex + np + nq)
      call sol%Init(a); call sol%DiagSym(a)
      call frac%ini(np)
      do i = 1, np
        if(abs(sol%eig%v(i)) < 1.d-8) then
          frac%v(i) = 1.d0
        else
          frac%v(i) = atan(sqrt(sol%eig%v(i))) / sqrt(sol%eig%v(i))
        end if
      end do
      call theta%DiagMat(frac)
      a = omg * sol%vec * theta * sol%vec%T()
      call sol%fin(); call frac%fin(); call theta%fin()
      do j = 1, np
        do i = 1, nq
          S%m(nex + np + i, nex + j) =   a%m(i, j)
          S%m(nex + j, nex + np + i) = - a%m(i, j)
        end do
      end do
    end function Omega2S
  end subroutine LeeSuzuki
end module Renormalization
