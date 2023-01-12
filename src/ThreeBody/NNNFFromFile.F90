module NNNFFromFile
  ! This module is for three nucleon interaction in harmonics oscillator basis
  ! from the hdf5 files by K.Hebeler
  ! The main object of this code is ChEFT3NIntHO.
  ! To compile this file, you would need
  !    -I/(path to hdf5 library)/hdf5.mod -lhdf5_fortran
  !
  ! You can find the example how to use this module at the end of this file.
  !
  !  o          o
  !   \        /
  !    \      /   momentum q in original Jacobi basis
  !     \    /    pi2 = sqrt{3/2} q for the transformation to HO basis
  !      \  /
  !       \/
  !        \      momentum p in original Jacobi basis
  !         \     pi1 = sqrt{2} p for the transformation to HO basis
  !          o
  ! In this module, pi1 and pi2 are replaced as p and q, respectively.
  ! T.Miyagi (TRIUMF) 21 AUG. 2019
  use, intrinsic :: iso_fortran_env
  !$ use omp_lib
  use hdf5
  use MyLibrary, only: hc
  implicit none

  integer, parameter, private :: dp = real64
  integer, parameter, private :: sp = real32

  type, private :: PWChan
    integer :: l12, s12, j12, t12, l3, j3
    integer, allocatable :: n12(:), n3(:), ns2i(:,:)
    integer :: Nho
  end type PWChan

  type, private :: Momentum
    real(dp) :: p, w
  end type Momentum

  type, private :: MomPQ
    integer :: ip, iq
    real(dp) :: p, q, pw, qw
  end type MomPQ

  type, private :: ModelSpace
    ! Momenta
    type(Momentum), allocatable :: PMom(:)
    type(Momentum), allocatable :: QMom(:)
    type(MomPQ), allocatable :: pq(:)
    integer, allocatable :: ipiq2ipq(:,:)
    integer :: Np, Nq, Npq
    ! HO radial quantum number
    integer :: Nmax
    real(dp) :: hw
    ! partial wave channels
    type(PWChan), private, allocatable :: alpha(:)
    integer, private :: Nalpha, Nalpha_file
    integer, allocatable :: LSJTlj2idx(:,:,:,:,:,:)
    integer :: jtot, ptot, ttot, J12max=-1
  contains
    procedure :: prt => PrintModelSpace
  end type ModelSpace

  type, private :: MomMat
    real(dp), allocatable :: v(:,:)
    type(PWChan), pointer :: ch_bra, ch_ket
    type(MomPQ), pointer :: pq(:)
    logical :: zero=.false.
  contains
    procedure :: PrintMomMat
  end type MomMat

  type, private :: ChEFT3NIntMom
    type(MomMat), allocatable :: ch(:,:)
    type(ModelSpace), pointer :: ms
  contains
    procedure :: InitChEFT3NIntMom
    procedure :: FinChEFT3NIntMom
    procedure :: PrintChEFT3NIntMom
    generic :: init => InitChEFT3NIntMom
    generic :: fin => FinChEFT3NIntMom
    generic :: prt => PrintChEFT3NIntMom
  end type ChEFT3NIntMom
  type(ChEFT3NIntMom) :: vmom


  type, private :: HOMat
    real(dp), allocatable :: v(:,:)
  end type HOMat

  type :: ChEFT3NIntHO
    type(HOMat), allocatable :: ch(:,:)
    type(ModelSpace) :: ms
    real(dp) :: lec=1._dp
  contains
    procedure :: InitChEFT3NIntHO
    procedure :: FinChEFT3NIntHO
    procedure :: GetChEFT3NIntHOFromQN
    procedure :: GetChEFT3NIntHOFromMat
    procedure :: SetLEC
    generic :: init => InitChEFT3NIntHO
    generic :: fin => FinChEFT3NIntHO
    generic :: get => GetChEFT3NIntHOFromQN, GetChEFT3NIntHOFromMat
  end type ChEFT3NIntHO

  logical, private :: interpolate=.true.
  ! w/ interpolation
  real(dp), private, allocatable :: ho_radial_wf(:,:,:) ! mesh, n, l
  real(dp), private, allocatable :: kmesh(:), wmesh(:)
  real(dp), private :: kmin = 0._dp, kmax = 8._dp
  integer, private :: NMesh = 60

  ! w/o interpolation
  real(dp), private, allocatable :: ho_radial_wfp(:,:,:) ! mesh, n, l
  real(dp), private, allocatable :: ho_radial_wfq(:,:,:) ! mesh, n, l

  real(dp), private :: LambdaCut = 500._dp
  integer, private :: RegulatorPow = 2
  !real(dp), private :: lambda_chi = 700._dp ! MeV
  !real(dp), private :: lambda_chif ! fm

  ! methods
  private :: FinChEFT3NIntHO
  private :: InitChEFT3NIntHO
  private :: GetChEFT3NIntHOFromQN
  private :: GetChEFT3NIntHOFromMat
  private :: SetLEC
  private :: read_file
  private :: set_ho_radial_quantum_numbers
  private :: FinChEFT3NIntMom
  private :: InitChEFT3NIntMom
  private :: read_numbers
  private :: read_momentum_mesh
  private :: read_pw_chan
  private :: read_matrix_elements
  private :: PrintModelSpace
  private :: multiply_regulator_channel
  private :: non_local_regulator
  private :: TransformMomToHO
  private :: TransformMomToHO_Test
  private :: trans_mom_to_ho
  private :: GetOverlapMomHO
  private :: GetOverlap_no_interp
  !private :: InterpolateMomentumMesh
  private :: store_ho_radial
  private :: store_ho_radial_pq
  private :: PrintChEFT3NIntMom
  private :: PrintMomMat


contains

  subroutine FinChEFT3NIntHO(this)
    class(ChEFT3NIntHO), intent(inout) :: this
    integer :: chbra, chket, ch
    do chbra = 1, this%ms%Nalpha
      do chket = 1, this%ms%Nalpha
        deallocate(this%ch(chbra,chket)%v)
      end do
    end do
    deallocate(this%ch)
    deallocate(this%ms%PMom)
    deallocate(this%ms%QMom)
    deallocate(this%ms%pq)
    deallocate(this%ms%ipiq2ipq)
    deallocate(this%ms%LSJTlj2idx)
    do ch = 1, this%ms%Nalpha
      deallocate(this%ms%alpha(ch)%n12)
      deallocate(this%ms%alpha(ch)%n3)
      deallocate(this%ms%alpha(ch)%ns2i)
    end do
    deallocate(this%ms%alpha)

    if(interpolate) then
      deallocate(kmesh, wmesh, ho_radial_wf)
    end if

    if(.not. interpolate) then
      deallocate(ho_radial_wfp, ho_radial_wfq)
    end if
  end subroutine FinChEFT3NIntHO

  function InitChEFT3NIntHO(this, J, P, T, Nmax, hw, &
      &      fn, RegulatorPower, Lambda, MeshNum, intrp, J12max) result(times)
    ! constructor
    ! input:
    ! J: total angular momentum of three-body Jacobi state (twice of physical J)
    ! P: total parity of three-body Jacobi state
    ! T: total isospin of three-body Jacobi state (twice of physical T)
    ! hw: frequency of harmonic oscillator basis (in unit of MeV)
    ! fn: full file name we have to read (hdf5 file from K. Hebeler)
    ! RegulatorPower: power of regulator function, f(p,q) = exp(- ( (p^2 + q^2) / (2Lambda^2) )^{regulator power}), default is 2
    ! Lambda: cut off of regulator function in unit of MeV, default is 500 MeV
    ! MeshNum: Number of Gauss-Legendre mesh you have to interpolate
    ! intrp: switch of interpolation. if true (false), interpolation will (won't) be done. default is true
    ! J12max: maximum value of total angular momentum in two-body state. By default, the calc. tries to use all information in the file.
    !
    ! output:
    ! times: three dimensional real values.
    !        times(1) is time used to read the file
    !        times(2) is time used to interpolate the matrix element
    !        times(3) is time used to transform from momentum to harmonic oscillator
    !
    use MyLibrary, only: gauss_legendre
    class(ChEFT3NIntHO), intent(inout) :: this
    real(dp) :: times(3)
    integer, intent(in) :: J, P, T, Nmax
    real(dp), intent(in) :: hw
    character(*), intent(in) :: fn
    integer, intent(in), optional :: RegulatorPower, MeshNum
    real(dp), intent(in), optional :: Lambda
    logical, intent(in), optional :: intrp
    integer, intent(in), optional :: J12max
    integer :: ch, chbra, chket, nbra, nket
    !$ real(dp) :: ti

    times(:) = 0._dp
    this%ms%Jtot = J
    this%ms%Ptot = P
    this%ms%Ttot = T
    this%ms%Nmax = Nmax
    this%ms%hw = hw
    if(present(MeshNum)) NMesh = MeshNum
    if(present(Lambda)) LambdaCut = Lambda
    if(present(RegulatorPower)) RegulatorPow = RegulatorPower
    if(present(J12max)) this%ms%J12max = J12max
    if(present(intrp)) interpolate = intrp
    allocate(this%ms%LSJTlj2idx(0:Nmax,0:1,0:Nmax+1,0:1,0:Nmax,Nmax+1))
    this%ms%LSJTlj2idx(:,:,:,:,:,:) = 0

    !$ ti = omp_get_wtime()
    call read_file(this%ms, fn)
    !$ times(1) = omp_get_wtime() - ti
    do ch = 1, this%ms%Nalpha
      if(this%ms%alpha(ch)%l12 > Nmax) cycle
      if(this%ms%alpha(ch)%j12 > Nmax+1) cycle
      if(this%ms%alpha(ch)%l3 > Nmax) cycle
      if((this%ms%alpha(ch)%j3+1)/2 > Nmax+1) cycle
      this%ms%LSJTlj2idx(this%ms%alpha(ch)%l12,&
        & this%ms%alpha(ch)%s12, &
        & this%ms%alpha(ch)%j12, &
        & this%ms%alpha(ch)%t12, &
        & this%ms%alpha(ch)%l3, &
        & (this%ms%alpha(ch)%j3+1)/2) = ch
    end do


    ! ho wave function
    if(interpolate) then
      allocate(kmesh(Nmesh), wmesh(Nmesh))
      kmin = max(minval(this%ms%PMom(:)%p), minval(this%ms%QMom(:)%p))
      kmax = min(maxval(this%ms%PMom(:)%p), maxval(this%ms%QMom(:)%p))
      call gauss_legendre(kmin, kmax, kmesh, wmesh, NMesh)
      allocate(ho_radial_wf(NMesh, 0:this%ms%Nmax/2, 0:this%ms%Nmax))
      ho_radial_wf(:,:,:) = 0._dp
      call store_ho_radial(NMesh, kmesh, this%ms%Nmax, this%ms%hw)
    end if

    if(.not. interpolate) then
      allocate(ho_radial_wfp(this%ms%Np, 0:this%ms%Nmax/2, 0:this%ms%Nmax))
      allocate(ho_radial_wfq(this%ms%Nq, 0:this%ms%Nmax/2, 0:this%ms%Nmax))
      ho_radial_wfp(:,:,:) = 0._dp
      ho_radial_wfq(:,:,:) = 0._dp
      call store_ho_radial_pq(this%ms)
    end if

    allocate(this%ch(this%ms%Nalpha,this%ms%Nalpha))
    do chbra = 1, this%ms%Nalpha
      nbra = this%ms%alpha(chbra)%nho
      do chket = 1, this%ms%Nalpha
        nket = this%ms%alpha(chket)%nho
        allocate(this%ch(chbra,chket)%v(nbra,nket))
        this%ch(chbra,chket)%v(:,:) = 0._dp
      end do
    end do
    !write(*,*) trim(fn)
    !open(15,file="test.txt")
    !call vmom%prt(15)
    !close(15)
    !stop

    !! Mom => HO
    call TransformMomToHO(this, times(2), times(3))
    !call TransformMomToHO_old(this, times(2), times(3))
    !call TransformMomToHO_Test(this)
    call vmom%fin()
  end function InitChEFT3NIntHO

  function GetChEFT3NIntHOFromQN(this, n12, l12, s12, j12, t12, n3, l3, j3, &
      &                                n45, l45, s45, j45, t45, n6, l6, j6) result(me)
    ! use this method to get the matrix element in the harmonic oscillator basis
    ! (-1)^{l/2} comes from fourier transformation phase
    ! (-1)^{l3+l6}/sqrt(3)^3 comes from our Jacobi momentum definition
    class(ChEFT3NIntHO), intent(in) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3 ! bra
    integer, intent(in) :: n45, l45, s45, j45, t45, n6, l6, j6 ! ket
    real(dp) :: me
    integer :: chbra, chket, nbra, nket
    me = 0._dp
    chbra = this%ms%LSJTlj2idx(l12,s12,j12,t12,l3,(j3+1)/2)
    chket = this%ms%LSJTlj2idx(l45,s45,j45,t45,l6,(j6+1)/2)
    if(chbra * chket == 0) return
    nbra = this%ms%alpha(chbra)%ns2i(n12,n3)
    nket = this%ms%alpha(chket)%ns2i(n45,n6)
    if(nbra * nket == 0) return
    me = this%ch(chbra,chket)%v(nbra,nket) &
      & * (-1._dp) ** ((l45+l6-l12-l3)/2) &
      & * (1._dp/sqrt(3._dp))**3 * this%LEC &
      & * (-1._dp) ** (l3+l6) * hc
#ifdef NNNFFromFileDebug
    write(*,*)
    write(*,"(a, 8i3)") "(n', L', S', J', N', l', j'): ", n12, l12, s12, j12, t12, n3, l3, j3
    write(*,"(a, 8i3)") "(n,  L,  S,  J,  N,  l,  j ): ", n45, l45, s45, j45, t45, n6, l6, j6
    write(*,"(a, 3i3, f12.6)") "(Jtot, Ptot, Ttot, MatEl): ", this%ms%Jtot, this%ms%Ptot, this%ms%Ttot, me
#endif
  end function GetChEFT3NIntHOFromQN

  function GetChEFT3NIntHOFromMat(this, chbra, chket, nbra, nket) result(me)
    ! use this method to get the matrix element in the harmonic oscillator basis
    class(ChEFT3NIntHO), intent(in) :: this
    integer, intent(in) :: chbra, nbra, chket, nket
    integer :: n12, l12, s12, j12, t12, n3, l3, j3
    integer :: n45, l45, s45, j45, t45, n6, l6, j6
    real(dp) :: me
    me = 0._dp
    n12 = this%ms%alpha(chbra)%n12(nbra)
    l12 = this%ms%alpha(chbra)%l12
    s12 = this%ms%alpha(chbra)%s12
    j12 = this%ms%alpha(chbra)%j12
    t12 = this%ms%alpha(chbra)%t12
    n3  = this%ms%alpha(chbra)%n3( nbra)
    l3  = this%ms%alpha(chbra)%l3
    j3  = this%ms%alpha(chbra)%j3

    n45 = this%ms%alpha(chket)%n12(nket)
    l45 = this%ms%alpha(chket)%l12
    s45 = this%ms%alpha(chket)%s12
    j45 = this%ms%alpha(chket)%j12
    t45 = this%ms%alpha(chket)%t12
    n6  = this%ms%alpha(chket)%n3( nket)
    l6  = this%ms%alpha(chket)%l3
    j6  = this%ms%alpha(chket)%j3
    me = this%get(n12,l12,s12,j12,t12,n3,l3,j3,n45,l45,s45,j45,t45,n6,l6,j6)
  end function GetChEFT3NIntHOFromMat

  subroutine SetLEC(this, lec)
    class(ChEFT3NIntHO), intent(inout) :: this
    real(dp), intent(in) :: lec
    this%LEC = lec
  end subroutine SetLEC

  subroutine read_file(ms, fn)
    type(ModelSpace), intent(inout) :: ms
    character(*), intent(in) :: fn
    integer :: error
    integer :: Np, Nq, Nalpha, ch
    integer :: ip, iq, ipq

    call h5open_f(error)

    call read_numbers(fn, Np, Nq, Nalpha)
    ms%Np = Np
    ms%Nq = Nq
    ms%Nalpha_file = Nalpha
    allocate(ms%PMom(Np))
    allocate(ms%QMom(Nq))

    call read_momentum_mesh(fn, "p mesh", sqrt(2._dp), ms%Np, ms%PMom)
    call read_momentum_mesh(fn, "q mesh", sqrt(1.5_dp), ms%Nq, ms%QMom)
    !call read_momentum_mesh(fn, "p mesh", 1.d0, ms%Np, ms%PMom)
    !call read_momentum_mesh(fn, "q mesh", 1.d0, ms%Nq, ms%QMom)
    ms%Npq = Np * Nq
    allocate(ms%pq(ms%Npq))
    allocate(ms%ipiq2ipq(Np,Nq))
    ipq = 0
    do ip = 1, Np
      do iq = 1, Nq
        ipq = ipq + 1
        ms%pq(ipq)%ip = ip
        ms%pq(ipq)%iq = iq
        ms%pq(ipq)%p = ms%PMom(ip)%p
        ms%pq(ipq)%pw = ms%PMom(ip)%w
        ms%pq(ipq)%q = ms%QMom(iq)%p
        ms%pq(ipq)%qw = ms%QMom(iq)%w
        ms%ipiq2ipq(ip,iq) = ipq
      end do
    end do

    call read_pw_chan(fn, "pw channels", ms)
    do ch = 1, ms%Nalpha
      call set_ho_radial_quantum_numbers(ms%alpha(ch), ms%Nmax)
    end do
#ifdef NNNFFromFileDebug
    call ms%prt()
#endif

    call vmom%init(ms)
    call read_matrix_elements(fn, "matrix elements", vmom)

    call h5close_f(error)
  end subroutine read_file

  subroutine set_ho_radial_quantum_numbers(alpha, Nmax)
    type(PWChan), intent(inout) :: alpha
    integer, intent(in) :: Nmax
    integer :: l12, l3
    integer :: cnt, N, n12, n3
    l12 = alpha%l12
    l3  = alpha%l3
    cnt = 0
    do N = 0, (Nmax-l12-l3)/2
      do n12 = 0, N
        n3 = N - n12
        cnt = cnt + 1
      end do
    end do

    alpha%Nho = cnt
    allocate(alpha%n12(cnt), alpha%n3(cnt))
    allocate(alpha%ns2i(0:Nmax/2, 0:Nmax/2))
    alpha%ns2i(:,:) = 0
    cnt = 0
    do N = 0, (Nmax-l12-l3)/2
      do n12 = 0, N
        n3 = N - n12
        cnt = cnt + 1
        alpha%n12(cnt) = n12
        alpha%n3( cnt) = n3
        alpha%ns2i(n12, n3) = cnt
      end do
    end do
  end subroutine set_ho_radial_quantum_numbers

  subroutine FinChEFT3NIntMom(this)
    class(ChEFT3NIntMom), intent(inout) :: this
    integer :: chbra, chket
    do chbra = 1, this%ms%Nalpha
      do chket = 1, this%ms%Nalpha
        deallocate(this%ch(chbra,chket)%v)
        this%ch(chbra,chket)%pq => null()
        this%ch(chbra,chket)%ch_bra => null()
        this%ch(chbra,chket)%ch_ket => null()
      end do
    end do
    deallocate(this%ch)
    this%ms => null()
  end subroutine FinChEFT3NIntMom

  subroutine InitChEFT3NIntMom(this, ms)
    class(ChEFT3NIntMom), intent(inout) :: this
    type(ModelSpace), intent(in), target :: ms
    integer :: chbra, chket
    this%ms => ms
    allocate(this%ch(ms%Nalpha, ms%Nalpha))
    do chbra = 1, ms%Nalpha
      do chket = 1, ms%Nalpha
        this%ch(chbra, chket)%pq => ms%pq
        this%ch(chbra, chket)%ch_bra => ms%alpha(chbra)
        this%ch(chbra, chket)%ch_ket => ms%alpha(chket)
        allocate(this%ch(chbra,chket)%v(ms%Npq, ms%Npq))
        this%ch(chbra,chket)%v(:,:) = 0._dp
        this%ch(chbra,chket)%zero=.false.
      end do
    end do
  end subroutine InitChEFT3NIntMom

  subroutine read_numbers(fn, Np, Nq, Nalpha)
    character(*), intent(in) :: fn
    integer, intent(out) :: Np, Nq, Nalpha
    integer(hid_t) :: hfile
    integer(hid_t) :: hd_Np
    integer(hid_t) :: hd_Nq
    integer(hid_t) :: hd_Nalpha
    integer(hsize_t), dimension(1) :: ndim
    integer :: error

    call h5fopen_f(fn, h5f_acc_rdonly_f, hfile, error)

    call h5dopen_f(hfile, "Np",     hd_Np,     error)
    call h5dopen_f(hfile, "Nq",     hd_Nq,     error)
    call h5dopen_f(hfile, "Nalpha", hd_Nalpha, error)

    ndim(1) = 1
    call h5dread_f(hd_Np,     h5t_native_integer, Np,     ndim, error)
    call h5dread_f(hd_Nq,     h5t_native_integer, Nq,     ndim, error)
    call h5dread_f(hd_Nalpha, h5t_native_integer, Nalpha, ndim, error)

    call h5dclose_f(hd_Np,     error)
    call h5dclose_f(hd_Nq,     error)
    call h5dclose_f(hd_Nalpha, error)
    call h5fclose_f(hfile, error)
  end subroutine read_numbers

  subroutine read_momentum_mesh(fn, dname, factor, Np, PMom)
    character(*), intent(in) :: fn, dname
    integer, intent(in) :: Np
    real(dp), intent(in) :: factor
    type(Momentum), intent(inout) :: PMom(:)
    integer(hsize_t), dimension(1) :: ndim
    integer(hid_t) :: hfile
    integer(hid_t) :: hd
    integer(hid_t) :: dtype
    integer(size_t) :: type_sizei, type_sized, type_size, offset
    integer :: error, idx
    real(dp), allocatable :: p(:), w(:)

    ndim(1) = Np
    allocate(p(Np), w(Np))
    call h5fopen_f(fn, h5f_acc_rdonly_f, hfile, error)
    call h5dopen_f(hfile, dname, hd, error)

    call h5tget_size_f(h5t_native_integer, type_sizei, error)
    call h5tget_size_f(h5t_native_double,  type_sized, error)
    type_size = type_sizei + 2*type_sized

    call h5tcreate_f(h5t_compound_f, type_sized, dtype, error)
    offset = 0
    call h5tinsert_f(dtype, "mesh point", offset, h5t_native_double, error)
    call h5dread_f(hd, dtype, p, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sized, dtype, error)
    call h5tinsert_f(dtype, "mesh weight", offset, h5t_native_double, error)
    call h5dread_f(hd, dtype, w, ndim, error)

    do idx = 1, Np
      Pmom(idx)%p = p(idx) * factor
      Pmom(idx)%w = w(idx) * factor
    end do
    call h5dclose_f(hd, error)
    call h5fclose_f(hfile, error)
    deallocate(p, w)
  end subroutine read_momentum_mesh

  subroutine read_pw_chan(fn, dname, ms)
    character(*), intent(in) :: fn, dname
    type(ModelSpace), intent(inout) :: ms
    integer(hsize_t), dimension(1) :: ndim
    integer(hid_t) :: hfile
    integer(hid_t) :: hd
    integer(hid_t) :: dtype
    integer(size_t) :: type_sizei, offset
    integer :: error, idx, cnt
    integer, allocatable :: l12(:), s12(:), j12(:), t12(:), l3(:), j3(:)

    ms%Nalpha = ms%Nalpha_file
    ndim(1) = ms%Nalpha_file
    allocate(l12(ms%Nalpha_file), s12(ms%Nalpha_file), j12(ms%Nalpha_file), t12(ms%Nalpha_file))
    allocate(l3(ms%Nalpha_file), j3(ms%Nalpha_file))
    l12 = 0; s12 = 0; j12 = 0; t12 = 0; l3 = 0; j3 = 0

    call h5fopen_f(fn, h5f_acc_rdonly_f, hfile, error)
    call h5dopen_f(hfile, dname, hd, error)

    call h5tget_size_f(h5t_native_integer, type_sizei, error)

    offset = 0
    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "L_12", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, l12, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "S_12", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, s12, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "J_12", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, j12, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "T_12", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, t12, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "l_3", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, l3, ndim, error)

    call h5tcreate_f(h5t_compound_f, type_sizei, dtype, error)
    call h5tinsert_f(dtype, "2*j_3", offset, h5t_native_integer, error)
    call h5dread_f(hd, dtype, j3, ndim, error)

    if(ms%J12max/=-1) then
      cnt = 0
      do idx = 1, ms%Nalpha_file
        if(j12(idx) > ms%J12max) cycle
        cnt = cnt + 1
      end do
      ms%Nalpha = cnt
    end if
    allocate(ms%alpha(ms%Nalpha))

    cnt = 0
    do idx = 1, ms%Nalpha_file
      if(ms%J12max>-1 .and. j12(idx) > ms%J12max) cycle
      cnt = cnt + 1
      ms%alpha(cnt)%l12 = l12(idx)
      ms%alpha(cnt)%s12 = s12(idx)
      ms%alpha(cnt)%j12 = j12(idx)
      ms%alpha(cnt)%t12 = t12(idx)
      ms%alpha(cnt)%l3  = l3( idx)
      ms%alpha(cnt)%j3  = j3( idx)
    end do

    call h5dclose_f(hd, error)
    call h5fclose_f(hfile, error)
    deallocate(l12, s12, j12, t12, l3, j3)
  end subroutine read_pw_chan

  subroutine read_matrix_elements(fn, dname, vmom)
    ! orignally, file is written like <alpha' q' p' | v | alpha q p> format.
    ! But, fortran is column major <p' q' alpha' | v | p q alpha>
    ! Caution: loop for partial wave channels is not general. Since the order of index in
    ! the original files are ascending order with respect to J12, it should be OK.
    !
    character(*), intent(in) :: fn, dname
    type(ChEFT3NIntMom), intent(inout) :: vmom
    integer(hsize_t), dimension(6) :: ndim
    integer(hid_t) :: hfile
    integer(hid_t) :: hd
    integer :: chbra, chket, ip_bra, iq_bra, ip_ket, iq_ket, ipq_bra, ipq_ket
    integer :: error
    real(sp), allocatable :: v_in(:,:,:,:,:,:)
    real(dp) :: v

    ndim(1) = vmom%ms%Np
    ndim(2) = vmom%ms%Nq
    ndim(3) = vmom%ms%Nalpha_file
    ndim(4) = vmom%ms%Np
    ndim(5) = vmom%ms%Nq
    ndim(6) = vmom%ms%Nalpha_file

    allocate(v_in(ndim(1), ndim(2), ndim(3), ndim(4), ndim(5), ndim(6)))
    v_in = 0.0
    call h5fopen_f(fn, h5f_acc_rdonly_f, hfile, error)
    call h5dopen_f(hfile, dname, hd, error)
    call h5dread_f(hd, h5t_native_real, v_in, ndim, error)
    !$omp parallel
    !$omp do private(chbra, iq_bra, ip_bra, ipq_bra,&
    !$omp &          chket, iq_ket, ip_ket, ipq_ket, v)
    do chbra = 1, vmom%ms%Nalpha
      do iq_bra = 1, vmom%ms%Nq
        do ip_bra = 1, vmom%ms%Np
          ipq_bra = vmom%ms%ipiq2ipq(ip_bra,iq_bra)

          do chket = 1, vmom%ms%Nalpha
            do iq_ket = 1, vmom%ms%Nq
              do ip_ket = 1, vmom%ms%Np
                ipq_ket = vmom%ms%ipiq2ipq(ip_ket,iq_ket)

                v = dble(v_in(ip_ket, iq_ket, chket, ip_bra, iq_bra, chbra))
                vmom%ch(chbra,chket)%v(ipq_bra, ipq_ket) = v
              end do
            end do
          end do

        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(chbra, chket)
    do chbra = 1, vmom%ms%Nalpha
      do chket = 1, vmom%ms%Nalpha
        if(maxval(abs(vmom%ch(chbra,chket)%v) ) < 1.d-16) then
          vmom%ch(chbra,chket)%zero = .true.
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

    call h5dclose_f(hd, error)
    call h5fclose_f(hfile, error)
    deallocate(v_in)
  end subroutine read_matrix_elements

  subroutine PrintModelSpace(ms)
    class(ModelSpace), intent(in) :: ms
    real(dp), allocatable :: p(:), q(:), pw(:), qw(:)
    integer :: idx

    allocate(p(ms%Np), pw(ms%Np))
    allocate(q(ms%Nq), qw(ms%Nq))
    do idx = 1, ms%np
      p(idx) = ms%Pmom(idx)%p
      pw(idx) = ms%Pmom(idx)%w
    end do

    do idx = 1, ms%nq
      q(idx) = ms%Qmom(idx)%p
      qw(idx) = ms%Qmom(idx)%w
    end do

    write(*,"(a,i2,a,i2,a,i2,a)") "Three-Body Model Space (Mom Space): J=", &
      & ms%Jtot, "/2, P=", ms%Ptot, ", ", ms%Ttot, "/2"
    write(*,"(a)") "Momentum mesh (p): "
    write(*,"(10f12.6)") p(:)
    write(*,"(a)") "Momentum mesh (q): "
    write(*,"(10f12.6)") q(:)

    write(*,"(a)") "angular-momentum channels: "
    write(*,"(7x,a,3x,a,3x,a,3x,a,3x,a,4x,a,2x,a,3x,a,3x,a,3x,a)") &
        &"idx", "l12", "s12", "j12", "t12", "l3", "2*j3", "2*J", "Par", "2*T"
    do idx = 1, ms%Nalpha
      write(*,"(4x,10i6)") idx, ms%alpha(idx)%l12, ms%alpha(idx)%s12, &
          & ms%alpha(idx)%j12, ms%alpha(idx)%t12, ms%alpha(idx)%l3, &
          & ms%alpha(idx)%j3, ms%Jtot, ms%Ptot, ms%Ttot
    end do

    deallocate(p, pw, q, qw)
  end subroutine PrintModelSpace

  subroutine multiply_regulator()
    integer :: chbra, chket
    do chbra = 1, vmom%ms%Nalpha
      do chket = 1, vmom%ms%Nalpha
        if(vmom%ch(chbra,chket)%zero) cycle
        call multiply_regulator_channel(vmom%ch(chbra,chket))
      end do
    end do
  end subroutine multiply_regulator

  subroutine multiply_regulator_channel(this)
    type(MomMat), intent(inout) :: this
    integer :: bra, ket
    real(dp) :: p_bra, q_bra, p_ket, q_ket, f_bra, f_ket
    do bra = 1, size(this%v, 1)
      p_bra = this%pq(bra)%p
      q_bra = this%pq(bra)%q
      f_bra = non_local_regulator(p_bra, q_bra)
      do ket = 1, size(this%v, 2)
        p_ket = this%pq(ket)%p
        q_ket = this%pq(ket)%q
        f_ket = non_local_regulator(p_ket, q_ket)

        this%v(bra,ket) = this%v(bra,ket) * f_bra * f_ket

      end do
    end do
  end subroutine multiply_regulator_channel

  function non_local_regulator(p, q) result(f)
    real(dp), intent(in) :: p, q ! in units of fm-1
    real(dp) :: f
    if(RegulatorPow==-1) then
      f = 1.d0
      return
    end if
    f = exp( - ( (0.5_dp * (p**2 + q**2) * hc**2) / LambdaCut**2 )**RegulatorPow )
  end function non_local_regulator


!  subroutine TransformMomToHO_old(this, time_intrp, time_trans)
!    ! transformation from momentum to harmonic oscillator basis
!    ! <n'N' alpha'| v | nN alpha> = \sum <n'N' | p'q'> <p'q'alpha'| v | p q alpha > <pq | nN>
!    ! n: radial quantum number of nucleon 1 and 2
!    ! N: radial quantum number of nucleon 3, relative to the center of mass of nucleon 1 and 2
!    class(ChEFT3NIntHO), intent(inout) :: this
!    real(dp), intent(out) :: time_intrp, time_trans
!    real(dp), allocatable :: ovlp_bra(:,:), ovlp_ket(:,:), MatMom(:,:), work(:,:)
!    integer :: chbra, chket, l12, l3, l45, l6
!    integer :: nb, nk, m
!    real(dp), allocatable :: p(:), q(:)
!    type(ModelSpace), pointer :: ms
!    integer :: i
!    !$ real(dp) :: t1
!    !! regulator
!    !call multiply_regulator()
!    time_intrp = 0._dp
!    time_trans = 0._dp
!    ms => vmom%ms
!    if(interpolate) then
!
!      m = NMesh**2
!      allocate(MatMom(m,m))
!      allocate(p(ms%Np), q(ms%Nq))
!      do i = 1, ms%np
!        p(i) = ms%Pmom(i)%p
!      end do
!      do i = 1, ms%nq
!        q(i) = ms%Qmom(i)%p
!      end do
!
!      do chbra = 1, ms%Nalpha
!        l12 = ms%alpha(chbra)%l12
!        l3  = ms%alpha(chbra)%l3
!        nb  = ms%alpha(chbra)%nho
!        if(nb < 1) cycle
!        allocate(ovlp_bra(nb,m))
!        !$ t1 = omp_get_wtime()
!        call GetOverlapMomHO(ovlp_bra, ms%alpha(chbra)%n12, ms%alpha(chbra)%n3, l12, l3)
!        !$ time_trans = time_trans + (omp_get_wtime() - t1)
!
!        do chket = 1, chbra
!          l45 = ms%alpha(chket)%l12
!          l6  = ms%alpha(chket)%l3
!          nk  = ms%alpha(chket)%nho
!          if(nk < 1) cycle
!          if(vmom%ch(chbra,chket)%zero) cycle
!          allocate(ovlp_ket(nk,m))
!          !$ t1 = omp_get_wtime()
!          call GetOverlapMomHO(ovlp_ket, ms%alpha(chket)%n12, ms%alpha(chket)%n3, l45, l6)
!          !$ time_trans = time_trans + (omp_get_wtime() - t1)
!
!          !$ t1 = omp_get_wtime()
!          call InterpolateMomentumMesh(MatMom, vmom%Ch(chbra,chket)%v, p, q, p, q, ms%ipiq2ipq, kmesh)
!          !$ time_intrp = time_intrp + (omp_get_wtime() - t1)
!
!          !$ t1 = omp_get_wtime()
!          allocate(work(m,nk))
!          call dgemm('n','t',m,nk,m,1._dp,MatMom,m,ovlp_ket,nk,0._dp,work,m)
!          call dgemm('n','n',nb,nk,m,1._dp,ovlp_bra,nb,work,m,0._dp,this%Ch(chbra,chket)%v,nb)
!          deallocate(work,ovlp_ket)
!          this%Ch(chket,chbra)%v = transpose(this%Ch(chbra,chket)%v)
!          !$ time_trans = time_trans + (omp_get_wtime() - t1)
!
!        end do
!        deallocate(ovlp_bra)
!      end do
!      deallocate(MatMom, p, q)
!      return
!    end if
!
!    m = ms%Npq
!    allocate(MatMom(m,m))
!    do chbra = 1, ms%Nalpha
!      l12 = ms%alpha(chbra)%l12
!      l3  = ms%alpha(chbra)%l3
!      nb  = ms%alpha(chbra)%nho
!      if(nb < 1) cycle
!      allocate(ovlp_bra(nb,m))
!      !$ t1 = omp_get_wtime()
!      call GetOverlap_no_interp(ovlp_bra, ms, ms%alpha(chbra)%n12, ms%alpha(chbra)%n3, l12, l3)
!      !$ time_trans = time_trans + (omp_get_wtime() - t1)
!      do chket = 1, chbra
!        l45 = ms%alpha(chket)%l12
!        l6  = ms%alpha(chket)%l3
!        nk  = ms%alpha(chket)%nho
!        if(nk < 1) cycle
!        allocate(ovlp_ket(nk,m))
!        !$ t1 = omp_get_wtime()
!        call GetOverlap_no_interp(ovlp_ket, ms, ms%alpha(chket)%n12, ms%alpha(chket)%n3, l45, l6)
!        !$ time_trans = time_trans + (omp_get_wtime() - t1)
!        MatMom = vmom%ch(chbra,chket)%v
!
!        !$ t1 = omp_get_wtime()
!        allocate(work(m,nk))
!        call dgemm('n','t',m,nk,m,1._dp,MatMom,m,ovlp_ket,nk,0._dp,work,m)
!        call dgemm('n','n',nb,nk,m,1._dp,ovlp_bra,nb,work,m,0._dp,this%Ch(chbra,chket)%v,nb)
!        deallocate(work,ovlp_ket)
!        this%Ch(chket,chbra)%v = transpose(this%Ch(chbra,chket)%v)
!        !$ time_trans = time_trans + (omp_get_wtime() - t1)
!
!      end do
!      deallocate(ovlp_bra)
!    end do
!    deallocate(MatMom)
!
!  end subroutine TransformMomToHO_old


  subroutine TransformMomToHO_Test(this)
    class(ChEFT3NIntHO), intent(inout) :: this
    real(dp), allocatable :: MatMom(:,:)
    integer :: chbra, chket, l12, l3, l45, l6
    integer :: nb, nk, nbra, nket, m
    integer :: n12, n3, n45, n6
    type(ModelSpace), pointer :: ms
    real(dp) :: me

    if(interpolate) return
    ms => vmom%ms
    m = ms%Npq
    allocate(MatMom(m,m))
    do chbra = 1, ms%Nalpha
      l12 = ms%alpha(chbra)%l12
      l3  = ms%alpha(chbra)%l3
      nb  = ms%alpha(chbra)%nho
      if(nb < 1) cycle
      do chket = 1, chbra
        l45 = ms%alpha(chket)%l12
        l6  = ms%alpha(chket)%l3
        nk  = ms%alpha(chket)%nho
        if(nk < 1) cycle

        MatMom = vmom%ch(chbra,chket)%v
        do nbra = 1, nb
          n12 = ms%alpha(chbra)%n12(nbra)
          n3  = ms%alpha(chbra)%n3( nbra)
          do nket = 1, nk
            n45 = ms%alpha(chket)%n12(nket)
            n6  = ms%alpha(chket)%n3( nket)
            me = trans_mom_to_ho(MatMom, ms, n12, l12, n3, l3, n45, l45, n6, l6)
            this%Ch(chbra,chket)%v(nbra,nket) = me
            this%Ch(chket,chbra)%v(nket,nbra) = me
          end do
        end do

      end do
    end do
    deallocate(MatMom)

  end subroutine TransformMomToHO_Test

  function trans_mom_to_ho(v, ms, n12, l12, n3, l3, n45, l45, n6, l6) result(me)
    real(dp) :: me
    real(dp), intent(in) :: v(:,:)
    type(ModelSpace), intent(in) :: ms
    integer, intent(in) :: n12, l12, n3, l3, n45, l45, n6, l6
    integer :: i, j, k, l, bra, ket

    me = 0._dp
    !$omp parallel
    !$omp do private(bra, i, j, ket, k, l) reduction(+:me)
    do bra = 1, ms%Npq
      i = ms%pq(bra)%ip
      j = ms%pq(bra)%iq
      do ket = 1, ms%Npq
        k = ms%pq(ket)%ip
        l = ms%pq(ket)%iq
        me = me + ms%PMom(i)%w * ms%PMom(i)%p * &
          &       ms%QMom(j)%w * ms%QMom(j)%p * &
          &       ms%PMom(k)%w * ms%PMom(k)%p * &
          &       ms%QMom(l)%w * ms%QMom(l)%p * &
          & ho_radial_wfp(i,n12,l12) * ho_radial_wfq(j,n3,l3) * &
          & ho_radial_wfp(k,n45,l45) * ho_radial_wfq(l,n6,l6) * &
          & v(bra,ket)
      end do
    end do
    !$omp end do
    !$omp end parallel

  end function trans_mom_to_ho


  subroutine GetOverlapMomHO(m,n12,n3,l12,l3)
    ! Here, overlap <nN|pq> is calculated.
    ! <nN|pq> = R_{nl}(p) R_{NL}(q) p^2 q^2 dp dq
    ! ho_radial_wf is normalized so \int (ho_radial_wf)^2 dp = 1
    real(dp), intent(inout) :: m(:,:)
    integer, intent(in) :: n12(:), n3(:), l12, l3
    integer :: iho, imom12, imom3, imom
    integer :: n12_i, n3_i
    real(dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

    !$omp parallel
    !$omp do private(iho, n12_i, n3_i, imom, imom12, imom3)
    do iho = 1, size(n12)
      n12_i = n12(iho)
      n3_i  = n3(iho)
      imom = 0
      do imom12 = 1, size(kmesh)
        do imom3 = 1, size(kmesh)
          imom = imom + 1
          m(iho,imom) = wmesh(imom12)*kmesh(imom12)* wmesh(imom3)*kmesh(imom3) * &
              & ho_radial_wf(imom12,n12_i,l12) * ho_radial_wf(imom3,n3_i,l3) * &
              & non_local_regulator(kmesh(imom12), kmesh(imom3))
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    return

    ! check normalization
    allocate(tmp1(size(n12), size(kmesh)**2))
    allocate(tmp2(size(kmesh)**2, size(kmesh)**2))
    do iho = 1, size(n12)
      n12_i = n12(iho)
      n3_i  = n3(iho)
      imom = 0
      do imom12 = 1, size(kmesh)
        do imom3 = 1, size(kmesh)
          imom = imom + 1
          tmp1(iho,imom) = ho_radial_wf(imom12,n12_i,l12) * ho_radial_wf(imom3,n3_i,l3)
        end do
      end do
    end do
    tmp2(:,:) = 0._dp
    imom = 0
    do imom12 = 1, size(kmesh)
      do imom3 = 1, size(kmesh)
        imom = imom + 1
        tmp2(imom,imom) = wmesh(imom12) * wmesh(imom3)
      end do
    end do
    allocate(tmp3(size(n12), size(kmesh)**2))
    call dgemm('n','n',size(n12),size(kmesh)**2,size(kmesh)**2,1._dp,tmp1,size(n12),tmp2,size(kmesh)**2,0._dp,tmp3,size(n12))
    deallocate(tmp2)
    allocate(tmp2(size(n12), size(n12)))
    call dgemm('n','t',size(n12),size(n12),size(kmesh)**2,1._dp,tmp3,size(n12),tmp1,size(n12),0._dp,tmp2,size(n12))
    write(*,*)
    do iho = 1, size(n12)
      write(*,"(10f12.6)") tmp2(iho,:)
    end do
    deallocate(tmp1, tmp2, tmp3)

  end subroutine GetOverlapMomHO

  subroutine GetOverlap_no_interp(m,ms,n12,n3,l12,l3)
    ! Here, overlap <nN|pq> is calculated.
    ! <nN|pq> = R_{nl}(p) R_{NL}(q) p^2 q^2 dp dq
    ! ho_radial_wf is normalized so \int (ho_radial_wf)^2 dp = 1
    real(dp), intent(inout) :: m(:,:)
    type(ModelSpace), intent(in) :: ms
    integer, intent(in) :: n12(:), n3(:), l12, l3
    integer :: iho, imom12, imom3, imom
    integer :: n12_i, n3_i
    real(dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

    do imom = 1, ms%Npq
      imom12 = ms%pq(imom)%ip
      imom3  = ms%pq(imom)%iq
      do iho = 1, size(n12)
        n12_i = n12(iho)
        n3_i  = n3(iho)
        m(iho,imom) = ms%PMom(imom12)%w * ms%PMom(imom12)%p * &
          &           ms%QMom(imom3 )%w * ms%QMom(imom3 )%p * &
          &           ho_radial_wfp(imom12,n12_i,l12) * ho_radial_wfq(imom3,n3_i,l3)
      end do
    end do
    return

    ! check normalization
    allocate(tmp1(size(n12), ms%Npq))
    allocate(tmp2(ms%Npq, ms%Npq))
    do imom = 1, ms%Npq
      imom12 = ms%pq(imom)%ip
      imom3  = ms%pq(imom)%iq
      do iho = 1, size(n12)
        n12_i = n12(iho)
        n3_i  = n3(iho)
        tmp1(iho,imom) = ho_radial_wfp(imom12,n12_i,l12) * ho_radial_wfq(imom3,n3_i,l3)
      end do
    end do

    tmp2(:,:) = 0._dp
    do imom = 1, ms%Npq
      imom12 = ms%pq(imom)%ip
      imom3  = ms%pq(imom)%iq
      tmp2(imom,imom) = ms%PMom(imom12)%w *  ms%QMom(imom3 )%w
    end do

    allocate(tmp3(size(n12), ms%Npq))
    call dgemm('n','n',size(n12),ms%Npq,ms%Npq,1._dp,tmp1,size(n12),tmp2,ms%Npq,0._dp,tmp3,size(n12))
    deallocate(tmp2)
    allocate(tmp2(size(n12), size(n12)))
    call dgemm('n','t',size(n12),size(n12),ms%Npq,1._dp,tmp3,size(n12),tmp1,size(n12),0._dp,tmp2,size(n12))
    write(*,*)
    do iho = 1, size(n12)
      write(*,"(10f12.6)") tmp2(iho,:)
    end do
    deallocate(tmp1, tmp2, tmp3)
  end subroutine GetOverlap_no_interp

!  subroutine InterpolateMomentumMesh(f2d, f2din, x1, x2, x3, x4, x1x2, xout)
!    ! Here, interpolating the three-body interaction matrix elements
!    use NdSpline
!    real(dp), intent(inout) :: f2d(:,:)
!    real(dp), intent(in) :: f2din(:,:)
!    real(dp), intent(in) :: x1(:), x2(:), x3(:), x4(:), xout(:)
!    integer, intent(in) :: x1x2(:,:)
!    real(dp), allocatable :: f4d(:,:,:,:), f1d(:)
!    integer, parameter :: k = 4
!    type(spline) :: sp
!    integer :: i1, i2, i3, i4, nn
!    integer :: n1, n2, n3, n4
!    integer :: i12, i34
!
!    n1 = size(x1)
!    n2 = size(x2)
!    n3 = size(x3)
!    n4 = size(x4)
!    nn  = size(xout)
!    allocate(f4d(n1,n2,n3,n4), f1d(n1*n2*n3*n4))
!
!    do i1 = 1, n1
!      do i2 = 1, n2
!        i12 = x1x2(i1,i2)
!        do i3 = 1, n3
!          do i4 = 1, n4
!            i34 = x1x2(i3,i4)
!            f4d(i1,i2,i3,i4) = f2din(i12,i34)
!          end do
!        end do
!      end do
!    end do
!    f1d = reshape(f4d, shape(f1d))
!    call sp%init([k,k,k,k], [n1,n2,n3,n4], [x1,x2,x3,x4], f1d)
!    deallocate(f1d, f4d)
!
!    call sp%interpolate([nn,nn,nn,nn], [xout,xout,xout,xout])
!
!    allocate(f4d(nn,nn,nn,nn))
!    f4d = reshape(sp%fs, shape(f4d))
!
!    !$omp parallel
!    !$omp do private(i12,i1,i2,i34,i3,i4)
!    do i12 = 1, nn**2
!      i1 = (i12-1)/nn + 1
!      i2 = mod(i12-1,nn)+1
!      do i34 = 1, nn**2
!        i3 = (i34-1)/nn + 1
!        i4 = mod(i34-1,nn)+1
!        f2d(i12,i34) = f4d(i1,i2,i3,i4)
!      end do
!    end do
!    !$omp end do
!    !$omp end parallel
!    deallocate(f4d)
!    call sp%fin()
!  end subroutine InterpolateMomentumMesh

  subroutine store_ho_radial(nmesh, p, nmax, hw)
    use MyLibrary, only: m_red_pn, ho_radial_wf_norm
    integer, intent(in) :: nmax, nmesh
    real(dp), intent(in) :: hw, p(:)
    integer :: n, l, i
    real(dp) :: a, pp, wave

    a = hc**2 / (2.d0 * m_red_pn * hw)
    do n = 0, nmax/2
      do l = 0, nmax
        do i = 1, nmesh
          pp = p(i)
          ho_radial_wf(i,n,l) = (-1._dp) ** n * ho_radial_wf_norm(n, l, a, pp)
        end do
      end do
    end do
  end subroutine store_ho_radial

  subroutine store_ho_radial_pq(ms)
    use MyLibrary, only: m_red_pn, ho_radial_wf_norm
    type(ModelSpace), intent(in) :: ms
    integer :: n, l, i
    real(dp) :: a, pp, wave

    a = hc**2 / (2.d0 * m_red_pn * ms%hw)
    do n = 0, ms%Nmax/2
      do l = 0, ms%Nmax
        do i = 1, ms%Np
          pp = ms%Pmom(i)%p
          ho_radial_wfp(i,n,l) = (-1._dp) ** n * ho_radial_wf_norm(n, l, a, pp)
        end do
      end do
    end do

    do n = 0, ms%Nmax/2
      do l = 0, ms%Nmax
        do i = 1, ms%Nq
          pp = ms%Qmom(i)%p
          ho_radial_wfq(i,n,l) = (-1._dp) ** n * ho_radial_wf_norm(n, l, a, pp)
        end do
      end do
    end do
  end subroutine store_ho_radial_pq

  ! for the test
  subroutine PrintChEFT3NIntMom(this, iunit)
    class(ChEFT3NIntMom), intent(in) :: this
    integer, intent(in) :: iunit
    integer :: chbra, chket
    do chbra = 1, this%ms%Nalpha
      do chket = 1, this%ms%Nalpha
        write(iunit,"(a,2i4)") "Channel: ", chbra, chket
        call this%ch(chbra,chket)%PrintMomMat(iunit)
      end do
    end do
  end subroutine PrintChEFT3NIntMom

  subroutine PrintMomMat(this, iunit)
    class(MomMat), intent(in) :: this
    integer, intent(in) :: iunit
    integer :: bra, ket
    do bra = 1, size(this%pq)
      do ket = 1, size(this%pq)
        write(iunit,"(4i4,es18.8)") this%pq(bra)%ip, this%pq(bra)%iq, &
          & this%pq(ket)%ip, this%pq(ket)%iq, this%v(bra,ket)
      end do
    end do
  end subroutine PrintMomMat

  subroutine TransformMomToHO(this, time_intrp, time_trans)
    ! transformation from momentum to harmonic oscillator basis
    ! <n'N' alpha'| v | nN alpha> = \sum <n'N' | p'q'> <p'q'alpha'| v | p q alpha > <pq | nN>
    ! n: radial quantum number of nucleon 1 and 2
    ! N: radial quantum number of nucleon 3, relative to the center of mass of nucleon 1 and 2
    use NdSpline
    class(ChEFT3NIntHO), intent(inout) :: this
    real(dp), intent(out) :: time_intrp, time_trans
    real(dp), allocatable :: ovlp_bra(:,:), ovlp_ket(:,:), work(:,:), f1d(:)
    integer :: chbra, chket, l12, l3, l45, l6
    integer :: nb, nk, m
    real(dp), allocatable :: Bp(:,:), Bq(:,:), coefs(:,:)
    type(ModelSpace), pointer :: ms
    integer, parameter :: k = 4
    type(spline) :: sp
    integer :: i, j
    !$ real(dp) :: t1
    time_intrp = 0._dp
    time_trans = 0._dp
    ms => vmom%ms
    if(interpolate) then
      do chbra = 1, ms%Nalpha
        l12 = ms%alpha(chbra)%l12
        l3  = ms%alpha(chbra)%l3
        nb  = ms%alpha(chbra)%nho
        if(nb < 1) cycle

        do chket = 1, chbra
          l45 = ms%alpha(chket)%l12
          l6  = ms%alpha(chket)%l3
          nk  = ms%alpha(chket)%nho
          if(nk < 1) cycle
          if(vmom%ch(chbra,chket)%zero) cycle
          !$ t1 = omp_get_wtime()
          allocate(f1d((ms%Np*ms%Nq)**2), coefs(ms%Np*ms%Nq, ms%Np*ms%Nq))
          f1d = reshape(vmom%ch(chbra,chket)%v, shape(f1d))
          call sp%init([k,k,k,k], [ms%Np,ms%Nq,ms%Np,ms%Nq], [ms%Pmom(:)%p, ms%Qmom(:)%p, ms%Pmom(:)%p, ms%Qmom(:)%p], f1d)
          deallocate(f1d)
          coefs = reshape(sp%get_coefs(), shape(coefs))
          Bp = sp%get_Bmatrix(1,kmesh)
          Bq = sp%get_Bmatrix(2,kmesh)
          call sp%fin()
          allocate(ovlp_bra(nb,ms%Npq))
          allocate(ovlp_ket(nk,ms%Npq))
          call GetOverlapMomRegHO(ovlp_bra, transpose(Bp), transpose(Bq), ms%alpha(chbra)%n12, ms%alpha(chbra)%n3, l12, l3)
          call GetOverlapMomRegHO(ovlp_ket, transpose(Bp), transpose(Bq), ms%alpha(chket)%n12, ms%alpha(chket)%n3, l45, l6)
          !$ time_intrp = time_intrp + (omp_get_wtime() - t1)

          !$ t1 = omp_get_wtime()
          allocate(work(vmom%ms%Npq,nk))
          call dgemm('n','t',ms%Npq,nk,ms%Npq,1._dp,coefs,ms%Npq,ovlp_ket,nk,0._dp,work,ms%Npq)
          call dgemm('n','n',nb,nk,ms%Npq,1._dp,ovlp_bra,nb,work,ms%Npq,0._dp,this%Ch(chbra,chket)%v,nb)
          deallocate(coefs,work,ovlp_bra,ovlp_ket)
          this%Ch(chket,chbra)%v = transpose(this%Ch(chbra,chket)%v)
          !$ time_trans = time_trans + (omp_get_wtime() - t1)

        end do
      end do
      return
    end if
  end subroutine TransformMomToHO

  subroutine GetOverlapMomRegHO(m,Sp,Sq,n12,n3,l12,l3)
    real(dp), intent(inout) :: m(:,:)
    real(dp), intent(in) :: Sp(:,:), Sq(:,:)
    integer, intent(in) :: n12(:), n3(:), l12, l3
    integer, allocatable :: mom_idx(:,:)
    real(dp), allocatable :: work1(:,:), work2(:,:)
    integer :: imom, imom12, imom3, imom_old, imom12_old, imom3_old
    integer :: iho, n12_i, n3_i

    allocate(work1(size(n12), size(kmesh)**2))
    allocate(work2(size(kmesh)**2, vmom%ms%Npq))
    allocate(mom_idx(2,size(kmesh)**2))

    imom = 0
    do imom12 = 1, size(kmesh)
      do imom3 = 1, size(kmesh)
        imom = imom + 1
        mom_idx(:,imom) = [imom12,imom3]
      end do
    end do

    !$omp parallel
    !$omp do private(imom, imom12, imom3, iho, n12_i, n3_i, imom_old, imom12_old, imom3_old)
    do imom = 1, size(kmesh)**2
      imom12 = mom_idx(1,imom)
      imom3  = mom_idx(2,imom)

      do iho = 1, size(n12)
        n12_i = n12(iho)
        n3_i = n3(iho)
        work1(iho,imom) = wmesh(imom12)*kmesh(imom12)* wmesh(imom3)*kmesh(imom3) * &
            & ho_radial_wf(imom12,n12_i,l12) * ho_radial_wf(imom3,n3_i,l3) * &
            & non_local_regulator(kmesh(imom12), kmesh(imom3))
      end do

      do imom_old = 1, vmom%ms%Npq
        imom12_old = vmom%ms%pq(imom_old)%ip
        imom3_old  = vmom%ms%pq(imom_old)%iq
        work2(imom,imom_old) = Sp(imom12,imom12_old) * Sq(imom3,imom3_old)
      end do
    end do
    !$omp end do
    !$omp end parallel
    call dgemm('n','n',size(work1,1),size(work2,2),size(work2,1),1._dp,work1,size(work1,1),work2,size(work2,1),0._dp,m,size(m,1))
    deallocate(work1, work2, mom_idx)
  end subroutine GetOverlapMomRegHO

end module NNNFFromFile

!
! Example (comment out below when you combine with your own code)
!
!program test
!  ! To see details, compile with -DNNNFFromFileDebug
!  ! test
!  ! < n12 l12 s12 j12 t12 n3 l3 j3: Jtot Ttot | V | n45 l45 s45 j45 t45 n6 l6 j6: Jtot Ttot >
!  use omp_lib
!  use NNNFFromFile
!
!  implicit none
!  type(ChEFT3NIntHO) :: vho
!  character(:), allocatable :: filename
!  integer :: Jtot, Ptot, Ttot, Nmax
!  integer :: n12, l12, s12, j12, t12, n3, l3, j3
!  integer :: n45, l45, s45, j45, t45, n6, l6, j6
!  real(8) :: hw, me, times(3)
!  character(256) :: char_J, char_P, char_T
!  Nmax = 4
!  Jtot = 1
!  Ptot = 1
!  TTot = 1
!  char_J = ""
!  char_P = ""
!  char_T = ""
!  write(char_J,*) Jtot
!  write(char_P,*) Ptot
!  write(char_T,*) Ttot
!  hw = 25.d0
!  filename = "3NF_V_J3_" // trim(adjustl(char_J)) // "_PAR_" &
!    & // trim(adjustl(char_P)) // "_T3_" // trim(adjustl(char_T)) // "_N2LO_c1.h5"
!  times = vho%init(1,1,1,Nmax,25.d0,filename,MeshNum=10, intrp=.false.) ! w/o interpolation
!  !times = vho%init(Jtot,Ptot,Ttot,Nmax,hw,filename,MeshNum=10)
!  call vho%SetLEC(1.d0)
!
!  l12 = 0; s12 = 0; j12 = 0; t12 = 1
!  l3 = 0; j3 = 1
!
!  l45 = 0; s45 = 0; j45 = 0; t45 = 1
!  l6 = 0; j6 = 1
!
!  do n12 = 0, Nmax/2
!    do n3 = 0, (Nmax-2*n12)/2
!      do n45 = 0, Nmax/2
!        do n6 = 0, (Nmax-2*n45)/2
!          me = vho%get( n12,l12,s12,j12,t12,n3,l3,j3, n45,l45,s45,j45,t45,n6,l6,j6 )
!        end do
!      end do
!    end do
!  end do
!
!  write(*,"(a,f12.4,a,f12.4,a,f12.4,a)") "# reading file:", &
!    & times(1), " sec, interpolation:", times(2), " sec, mom->HO:", times(3), " sec"
!
!end program test
!
