!
! This module is for the conversion from partial-wave basis to cartesian basis.
! The output can be used as an input of the nuclear matter calculations on the lattice
!
module NNForceCartesian
  use omp_lib
  use ClassSys
  use Profiler, only: timer
  use LinAlgLib
  use StoreCouplings, only: CGsStore

  implicit none
  private :: interpolate_submatrix
  private :: interpolate_channel
  private :: interpolate
  private :: find_momenta_interpolated
  private :: SetNNCartesian
  private :: WriteFile
  private :: GetMeanSquaredError
  private :: FinNNCartesian
  private :: InitNNCartesian
  private :: FinNNCartesianChan
  private :: InitNNCartesianChan
  private :: GetChannelIndex
  private :: GetNumberChannels
  private :: FinSpinPNChannels
  private :: InitSpinPNChannels
  private :: FinSpinPNChannel
  private :: InitSpinPNChannel
  private :: GetGridIndex
  private :: GetNumberGrids
  private :: InitMomentumGrids
  private :: FinMomentumGrids
  private :: store_arrays
  private :: release_arrays
  private :: fix_phase

  type, private :: MomentumGrid
    integer :: ix, iy, iz    ! index of x, y, z ( Cartesian)
    real(8) :: px, py, pz    ! p = (px, py, pz )   Cartesian
    real(8) :: p, theta, phi ! p = (p, theta, phi) Spherical
  end type MomentumGrid

  type, private :: MomentumGrids
    type(MomentumGrid), allocatable :: grids(:) ! 1d momentum grids
    integer, allocatable :: grid_to_idx(:,:,:)  ! (ix, iy, iz) -> index
    integer :: n_abscissa, Nmax
    real(8), allocatable :: abscissa(:)
    real(8) :: Lbox = -1.d0
  contains
    procedure :: InitMomentumGrids
    procedure :: InitMomentumGridsFromAbscissa
    procedure :: FinMomentumGrids
    procedure :: GetNumberGrids
    procedure :: GetGridIndex
    generic :: init => InitMomentumGrids, InitMomentumGridsFromAbscissa
    generic :: fin => FinMomentumGrids
  end type MomentumGrids

  type, private :: SpinPNChannel
    integer :: Sz1, Sz2, Tz ! sz1, sz2: 3rd components of spin, Tz=-1, 0, 1, corresponding pp, pn, nn
    type(MomentumGrids), pointer :: mom
  contains
    procedure :: InitSpinPNChannel
    procedure :: FinSpinPNChannel
    generic :: init => InitSpinPNChannel
    generic :: fin => FinSpinPNChannel
  end type SpinPNChannel

  type, private :: SpinPNChannels
    type(SpinPNChannel), allocatable :: channels(:)
    integer, allocatable :: spin_pn_idx(:,:,:) ! (sz1, sz2, Tz) -> index
    type(MomentumGrids) :: mom
    logical :: rspace=.false.
  contains
    procedure :: InitSpinPNChannels
    procedure :: FinSpinPNChannels
    procedure :: GetNumberChannels
    procedure :: GetChannelIndex
    generic :: init => InitSpinPNChannels
    generic :: fin => FinSpinPNChannels
  end type SpinPNChannels

  type, extends(CMat), private :: NNCartesianChan
    type(SpinPNChannel), pointer :: ch_bra, ch_ket
    logical :: is_zero=.true.
  contains
    procedure :: InitNNCartesianChan
    procedure :: FinNNCartesianChan
    generic :: init => InitNNCartesianChan
    generic :: release => FinNNCartesianChan
  end type NNCartesianChan

  type :: NNCartesian
    type(NNCartesianChan), allocatable :: MatCh(:,:)
    type(SpinPNChannels) :: channels
  contains
    procedure :: InitNNCartesian
    procedure :: FinNNCartesian
    procedure :: SetNNCartesian
    procedure :: WriteFile
    procedure :: GetMeanSquaredError
    procedure :: GetMERelIndex
    procedure :: GetMESpIndex
    generic :: init => InitNNCartesian
    generic :: fin => FinNNCartesian
    generic :: GetME => GetMERelIndex, GetMESpIndex
  end type NNCartesian

  type(CGsStore), private :: cg_spin, cg_lsj
  complex(8), allocatable, private :: sphs(:,:,:)
contains

  subroutine FinMomentumGrids(this)
    class(MomentumGrids), intent(inout) :: this
    deallocate(this%grids)
    deallocate(this%grid_to_idx)
  end subroutine FinMomentumGrids

  subroutine InitMomentumGridsFromAbscissa(this, abscissa)
    class(MomentumGrids), intent(inout) :: this
    real(8), intent(in) :: abscissa(:)
    real(8) :: px, py, pz, p, theta, phi
    integer :: nx, ny, nz
    integer :: ix, iy, iz, idx
    nx = size(abscissa)
    ny = size(abscissa)
    nz = size(abscissa)
    this%Nmax = (size(abscissa)-1)/2
    allocate(this%grids( nx*ny*nz ))
    allocate(this%grid_to_idx( nx,ny,nz ))
    idx = 0
    do ix = 1, nx
      do iy = 1, ny
        do iz = 1, nz
          idx = idx + 1
          px = abscissa(ix)
          py = abscissa(iy)
          pz = abscissa(iz)
          this%grids(idx)%px = px
          this%grids(idx)%py = py
          this%grids(idx)%pz = pz
          this%grids(idx)%ix = ix
          this%grids(idx)%iy = iy
          this%grids(idx)%iz = iz
          p = sqrt(px**2 + py**2 + pz**2)
          if( px==0.d0 .and. py==0.d0 .and. pz==0.d0 ) then
            theta = 0.d0
          else
            theta = acos( pz / p )
          end if
          if( px==0.d0 .and. py==0.d0 ) then
            phi = 0.d0
          else
            phi = sign(acos( px / sqrt(px**2+py**2) ), py)
          end if
          this%grids(idx)%p = p
          this%grids(idx)%theta = theta
          this%grids(idx)%phi = phi
          !write(*,'(3f12.6,a,3f12.6)') px, py, pz, " -> ", p*sin(theta)*cos(phi), p*sin(theta)*sin(phi), p*cos(theta)
          this%grid_to_idx(ix,iy,iz) = idx
        end do
      end do
    end do
  end subroutine InitMomentumGridsFromAbscissa

  subroutine InitMomentumGrids(this, Lbox, Nmax)
    use MyLibrary, only: pi
    class(MomentumGrids), intent(inout) :: this
    real(8), intent(in) :: Lbox
    integer, intent(in) :: Nmax
    real(8), allocatable :: abscissa(:)
    real(8) :: px, py, pz, p, theta, phi
    integer :: nx, ny, nz
    integer :: ix, iy, iz, idx
    this%Lbox = Lbox
    this%Nmax = Nmax
    nx = 4*this%Nmax+1
    ny = 4*this%Nmax+1
    nz = 4*this%Nmax+1
    allocate(abscissa(-nx/2:nx/2))
    allocate(this%abscissa(-nx/2:nx/2))
    do ix = -nx/2, nx/2
      abscissa(ix) = pi * dble(ix) / Lbox ! not 2 pi n / L because of that the relative momentum is (p'-p)/2
    end do
    this%abscissa(:) = abscissa(:)
    allocate(this%grids( nx*ny*nz ))
    allocate(this%grid_to_idx( -nx:nx,-ny:ny,-nz:nz ))
    idx = 0
    do ix = -nx/2, nx/2
      do iy = -ny/2, ny/2
        do iz = -nz/2, nz/2
          idx = idx + 1
          px = abscissa(ix)
          py = abscissa(iy)
          pz = abscissa(iz)
          this%grids(idx)%px = px
          this%grids(idx)%py = py
          this%grids(idx)%pz = pz
          this%grids(idx)%ix = ix
          this%grids(idx)%iy = iy
          this%grids(idx)%iz = iz
          p = sqrt(px**2 + py**2 + pz**2)
          if( px==0.d0 .and. py==0.d0 .and. pz==0.d0 ) then
            theta = 0.d0
          else
            theta = acos( pz / p )
          end if
          if( px==0.d0 .and. py==0.d0 ) then
            phi = 0.d0
          else
            phi = sign(acos( px / sqrt(px**2+py**2) ), py)
          end if
          this%grids(idx)%p = p
          this%grids(idx)%theta = theta
          this%grids(idx)%phi = phi
          this%grid_to_idx(ix,iy,iz) = idx
        end do
      end do
    end do
    deallocate(abscissa)
  end subroutine InitMomentumGrids

  function GetNumberGrids(this) result(n)
    class(MomentumGrids), intent(in) :: this
    integer :: n
    n = size(this%grids)
  end function GetNumberGrids

  function GetGridIndex(this,ix,iy,iz) result(n)
    class(MomentumGrids), intent(in) :: this
    integer, intent(in) :: ix, iy, iz
    integer :: n
    n = this%grid_to_idx(ix,iy,iz)
  end function GetGridIndex

  subroutine InitSpinPNChannel(this, Sz1, Sz2, Tz, grids)
    class(SpinPNChannel), intent(inout) :: this
    integer, intent(in) :: Sz1, Sz2, Tz
    type(MomentumGrids), intent(in), target :: grids
    this%Sz1 = Sz1
    this%Sz2 = Sz2
    this%Tz = Tz
    this%mom => grids
  end subroutine InitSpinPNChannel

  subroutine FinSpinPNChannel(this)
    class(SpinPNChannel), intent(inout) :: this
    this%mom => null()
  end subroutine FinSpinPNChannel

  subroutine InitSpinPNChannels(this, abscissa, Lbox, Nmax, rspace)
    class(SpinPNChannels), intent(inout) :: this
    real(8), intent(in), optional :: abscissa(:), Lbox
    integer, intent(in), optional :: Nmax
    logical, intent(in), optional :: rspace
    integer :: Sz1, Sz2, Tz, idx

    if( present( abscissa ) ) then
      call this%mom%init( abscissa )
    elseif( present( Lbox ) .and. present( Nmax )) then
      call this%mom%init( Lbox, Nmax )
    else
      write(*,*) "Error, see at ", __LINE__, " in ", __FILE__
    end if
    if( present( rspace )) this%rspace = rspace
    allocate( this%spin_pn_idx(-1:1,-1:1,-1:1))
    this%spin_pn_idx(:,:,:) = 0
    idx = 0
    do Sz1 = -1, 1, 2
      do Sz2 = -1, 1, 2
        do Tz = -1, 1
          idx = idx + 1
          this%spin_pn_idx(Sz1,Sz2,Tz) = idx
        end do
      end do
    end do
    allocate( this%channels( maxval(this%spin_pn_idx) ) )
    do Sz1 = -1, 1, 2
      do Sz2 = -1, 1, 2
        do Tz = -1, 1
          call this%channels( this%spin_pn_idx(Sz1,Sz2,Tz) )%init(Sz1,Sz2,Tz,this%mom)
        end do
      end do
    end do
  end subroutine InitSpinPNChannels

  subroutine FinSpinPNChannels(this)
    class(SpinPNChannels), intent(inout) :: this
    integer :: ich

    do ich = 1, this%GetNumberChannels()
      call this%channels(ich)%fin()
    end do
    deallocate( this%channels )
    deallocate( this%spin_pn_idx )
    call this%mom%fin()
    this%rspace = .false.
  end subroutine FinSpinPNChannels

  function GetNumberChannels(this) result(n)
    class(SpinPNChannels), intent(in) :: this
    integer :: n
    n = size(this%channels)
  end function GetNumberChannels

  function GetChannelIndex(this,Sz1,Sz2,Tz) result(n)
    class(SpinPNChannels), intent(in) :: this
    integer, intent(in) :: Sz1, Sz2, Tz
    integer :: n
    n = this%spin_pn_idx(Sz1,Sz2,Tz)
  end function GetChannelIndex

  subroutine FinNNCartesianChan( this )
    class(NNCartesianChan), intent(inout) :: this
    call this%fin()
    this%ch_bra => null()
    this%ch_ket => null()
  end subroutine FinNNCartesianChan

  subroutine InitNNCartesianChan( this, chbra, chket )
    class(NNCartesianChan), intent(inout) :: this
    type(SpinPNChannel), intent(in), target :: chbra, chket
    this%is_zero = .false.
    this%ch_bra => chbra
    this%ch_ket => chket
    call this%zeros( chbra%mom%GetNumberGrids(), chket%mom%GetNumberGrids() )
  end subroutine InitNNCartesianChan

  subroutine InitNNCartesian( this, abscissa, Lbox, Nmax, rspace )
    class(NNCartesian), intent(inout), target :: this
    real(8), intent(in), optional :: abscissa(:), Lbox
    integer, intent(in), optional :: Nmax
    logical, intent(in), optional :: rspace
    integer :: ichbra, ichket
    type(SpinPNChannel), pointer :: chbra, chket

    call this%channels%init( abscissa, Lbox, Nmax, rspace )
    allocate( this%MatCh( this%channels%GetNumberChannels(), this%channels%GetNumberChannels() ))
    do ichbra = 1, this%channels%GetNumberChannels()
      chbra => this%channels%channels(ichbra)
      do ichket = 1, this%channels%GetNumberChannels()
        chket => this%channels%channels(ichket)
        if( chbra%Tz /= chket%Tz ) cycle
        call this%MatCh(ichbra,ichket)%init( this%channels%channels(ichbra), this%channels%channels(ichket) )
      end do
    end do
  end subroutine InitNNCartesian

  subroutine FinNNCartesian( this )
    class(NNCartesian), intent(inout) :: this
    integer :: ichbra, ichket
    do ichbra = 1, this%channels%GetNumberChannels()
      do ichket = 1, this%channels%GetNumberChannels()
        call this%MatCh(ichbra,ichket)%release()
      end do
    end do
    deallocate( this%MatCh )
    call this%channels%fin()
  end subroutine FinNNCartesian

  function GetMeanSquaredError(this, that) result(r)
    class(NNCartesian), intent(in), target :: this
    type(NNCartesian), intent(in) :: that
    integer :: ichbra, ichket, bra, ket
    real(8) :: r
    type(MomentumGrids), pointer :: g

    g => this%channels%mom
    r = 0.d0
    do ichbra = 1, this%channels%GetNumberChannels()
      do ichket = 1, this%channels%GetNumberChannels()
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        do bra = 1, g%GetNumberGrids()
          do ket = 1, g%GetNumberGrids()
            r = r + abs(this%MatCh(ichbra,ichket)%m(bra,ket)-that%MatCh(ichbra,ichket)%m(bra,ket))**2
          end do
        end do
      end do
    end do
  end function GetMeanSquaredError

  function GetMERelIndex( this, ix, iy, iz, si1, si2, zi, jx, jy, jz, sj1, sj2, zj) result(r)
    class(NNCartesian), intent(in), target :: this
    integer, intent(in) :: ix, iy, iz, si1, si2, zi
    integer, intent(in) :: jx, jy, jz, sj1, sj2, zj
    integer :: ichbra, ichket, bra, ket
    type(MomentumGrids), pointer :: g
    complex(8) :: r
    r = 0.d0
    g => this%channels%mom
    if( abs(ix) > 2*g%Nmax ) return
    if( abs(iy) > 2*g%Nmax ) return
    if( abs(iz) > 2*g%Nmax ) return
    if( abs(jx) > 2*g%Nmax ) return
    if( abs(jy) > 2*g%Nmax ) return
    if( abs(jz) > 2*g%Nmax ) return
    ichbra = this%channels%GetChannelIndex(si1, si2, zi)
    ichket = this%channels%GetChannelIndex(sj1, sj2, zj)
    bra = g%GetGridIndex(ix, iy, iz)
    ket = g%GetGridIndex(jx, jy, jz)
    r = this%MatCh(ichbra,ichket)%m(bra,ket)
  end function GetMERelIndex

  function GetMESpIndex( this, i, j, k, l ) result(r)
    !
    ! The relative momentum is P_i = (p1_i - p2_i)/2, using the single-particle momenta:
    ! p1_i = (2 pi n1_i + theta_i ) / L
    ! p2_i = (2 pi n2_i + theta_i ) / L
    ! P_i =  ( n1_i - n2_i ) / L
    ! i = x, y, z
    ! n1_i, n2_i are integer
    ! Input:
    !   i, j, k, l: 1d array ex) n1_x = i(1), n1_y = i(2), n1_z = i(3), sz1 = i(4), tz1 = i(5)
    ! Output:
    !   r: matrix element
    !
    class(NNCartesian), intent(in), target :: this
    integer, intent(in) :: i(5), j(5), k(5), l(5)
    integer :: ii(5), jj(5), kk(5), ll(5)
    complex(8) :: r
    real(8) :: phase

    r = 0.d0
    phase = 1.d0
    ii = i; jj = j; kk = k; ll = l
    if( i(5) == 1 .and. j(5) == -1 ) then
      ii = j; jj = i; phase = phase * (-1.d0)
    end if

    if( k(5) == 1 .and. l(5) == -1 ) then
      kk = l; ll = k; phase = phase * (-1.d0)
    end if
    r = this%GetME( ii(1)-jj(1), ii(2)-jj(2), ii(3)-jj(3), ii(4), jj(4), (ii(5)+jj(5))/2, &
        &           kk(1)-ll(1), kk(2)-ll(2), kk(3)-ll(3), kk(4), ll(4), (kk(5)+ll(5))/2 ) * phase
  end function GetMESpIndex

  subroutine WriteFile( this, filename )
    use MyLibrary, only: hc, pi
    class(NNCartesian), intent(in) , target:: this
    character(*), intent(in), optional :: filename
    character(:), allocatable :: fn
    integer :: wunit = 20
    character(255) :: header
    integer :: ichbra, ichket, bra, ket, ket_max
    type(SpinPNChannel), pointer :: chbra, chket
    type(MomentumGrids), pointer :: g

    g => this%channels%mom
    fn = "nn_cartesian.dat"
    if( present(filename) ) fn = filename
    open(wunit, file=fn)
    header = "# NN Force on 3d momentum grid, generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // " " // trim(VERSION)
#endif
    write(wunit,"(a)") trim(header)
    write(wunit,'(a,f8.4,a,i3)') "Box size: ", g%Lbox, " fm, Nmax: ", g%Nmax
    if(.not. this%channels%rspace) then
      write(wunit,'(a)')  "   idx,   ix,   iy,   iz,               px,               py,               pz"
    end if
    if(this%channels%rspace) then
      write(wunit,'(a)')  "   idx,   ix,   iy,   iz,                x,                y,                z"
    end if
    do ket = 1, g%GetNumberGrids()
      write(wunit,'(4i6,3es18.8)') ket, g%grids(ket)%ix, g%grids(ket)%iy, g%grids(ket)%iz, &
          & g%grids(ket)%px, g%grids(ket)%py, g%grids(ket)%pz
    end do
    if( g%Lbox > 0.d0 ) write(wunit,'(a)') "  idx', s'z1, s'z2,  idx,  sz1,  sz2,   Tz,    ME real (MeV),    ME imag (MeV)"
    if( g%Lbox <= 0.d0 ) write(wunit,'(a)') "  idx', s'z1, s'z2,  idx,  sz1,  sz2,   Tz,  ME real (MeV-2),  ME imag (MeV-2)"
    if( .not. this%channels%rspace ) then
      do ichbra = 1, this%channels%GetNumberChannels()
        chbra => this%channels%channels(ichbra)
        do ichket = 1, ichbra
          chket => this%channels%channels(ichket)
          if( this%MatCh(ichbra,ichket)%is_zero ) cycle

          do bra = 1, g%GetNumberGrids()
            ket_max = g%GetNumberGrids()
            if( ichbra==ichket ) ket_max = bra
            do ket = 1, ket_max

              write(wunit,'(7i6,2es18.8)') &
                  & bra, chbra%Sz1, chbra%Sz2, &
                  & ket, chket%Sz1, chket%Sz2, &
                  & chket%Tz, this%MatCh(ichbra,ichket)%m(bra,ket)

            end do
          end do

        end do
      end do
    end if

    if( this%channels%rspace ) then
      do ichbra = 1, this%channels%GetNumberChannels()
        chbra => this%channels%channels(ichbra)
        do ichket = 1, ichbra
          chket => this%channels%channels(ichket)
          if( this%MatCh(ichbra,ichket)%is_zero ) cycle

          do ket = 1, g%GetNumberGrids()

            write(wunit,'(7i6,2es18.8)') &
                & ket, chbra%Sz1, chbra%Sz2, &
                & ket, chket%Sz1, chket%Sz2, &
                & chket%Tz, this%MatCh(ichbra,ichket)%m(ket,ket)

          end do

        end do
      end do
    end if
    close(wunit)
  end subroutine WriteFile

  subroutine SetNNCartesian( this, params, interap )
    use MyLibrary, only: gauss_legendre
    use NuHamilInput, only: InputParameters
    use TwoBodyRelativeSpace
    use MyLibrary, only: gauss_legendre, pi, spherical_harmonics
    use NNForce, only: NNForceMom
    class(NNCartesian), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: interap
    type(InputParameters) :: input
    type(TwoBodyRelSpaceSpinMBasis) :: mom, mom_inter, rspace
    type(NNForceMom) :: vmom, vmom_inter, vr
    integer :: npmesh, nqmesh
    integer, allocatable :: grididx_to_pmom(:)
    real(8), allocatable :: p(:), w(:), pmom(:), wmom(:), r(:)
    logical :: test_mode = .false.
    real(8) :: time, rmax=10.d0
    integer :: i, n_r  = 200
    type(sys) :: s

    if( present(interap )) test_mode = .not. interap
    input = params
    input%renorm_space2 = 'mom'
    write(*,'(2a)',advance='no') 'NNint = ', trim(input%NNint%val)
    write(*,'(a,i3)',advance='no') ', Jmax = ', input%J2max_NNint
    write(*,*)

    ! Prepare the interaction for 3d grid
    call find_momenta_interpolated( pmom, grididx_to_pmom, this%channels%mom )
    call mom_inter%init(minval(pmom),maxval(pmom),size(pmom),input%J2max_NNint)
    allocate( wmom( size(pmom) ))
    wmom = 0.d0
    call mom_inter%SetMeshWeight(pmom,wmom)
    deallocate( pmom, wmom )
    call vmom_inter%init(mom_inter)

    if( test_mode ) then
      ! test mode, no interpolation, i.e., renormalization producere cannot be used
      call vmom_inter%SetNNForceMom(input)
      call fix_phase( vmom_inter )
    end if

    if( .not. test_mode ) then
      ! Here get the original interaction on Gauss-Legendre mesh
      call mom%init(0.d0,input%pmax2,input%NMesh2,input%J2max_NNint)
      allocate(pmom(input%NMesh2), wmom(input%NMesh2))
      if(input%renorm%val == 'vlowk') then
        npmesh = int(dble(input%NMesh2) * input%lambda / input%pmax2)
        nqmesh = input%NMesh2 - npmesh
        call gauss_legendre(0.d0, input%lambda, p, w, npmesh)
        pmom(:npmesh) = p
        wmom(:npmesh) = w
        call gauss_legendre(input%lambda, input%pmax2, p, w, nqmesh)
        pmom(npmesh+1:) = p
        wmom(npmesh+1:) = w
        deallocate(p,w)
      else
        call gauss_legendre(0.d0, input%pmax2, p, w, input%NMesh2)
        pmom = p
        wmom = w
        deallocate(p,w)
      end if
      call mom%SetMeshWeight(pmom,wmom)
      deallocate(pmom,wmom)
      call vmom%init(mom)
      call vmom%SetNNForceMom(input)
      call fix_phase( vmom )

      if( this%channels%rspace ) then
        allocate( r(n_r), wmom(n_r) )
        do i = 1, n_r
          r(i) = rmax * dble(i) / dble(n_r)
          wmom(i) = 0.d0
        end do
        call rspace%init(minval(r),maxval(r),size(r),input%J2max_NNint)
        call rspace%SetMeshWeight(r,wmom)
        deallocate( r, wmom )
        call vr%init(rspace)
        call transform_to_rspace( vmom, vr )
        call interpolate( vmom_inter, vr )
        call vr%fin()
        call rspace%fin()
      else
        call interpolate( vmom_inter, vmom )
      end if
      call vmom%fin()
      call mom%fin()
    end if

    ! Here partial-wave -> cartesian
    time = omp_get_wtime()
    call to_3d_spin_basis( this, vmom_inter, grididx_to_pmom )
    call timer%add(s%str("NN convert to on 3d grid"), omp_get_wtime()-time)
    call vmom_inter%fin()
    call mom_inter%fin()
  end subroutine SetNNCartesian

  subroutine to_3d_spin_basis( this, v, gidx_to_pidx )
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceMom
    type(NNCartesian), intent(inout), target :: this
    type(NNForceMom), intent(in) :: v
    integer, intent(in) :: gidx_to_pidx(:)
    integer :: ichbra, ichket
    type(MomentumGrids), pointer :: g
    type(SpinPNChannel), pointer :: chbra, chket

    g => this%channels%mom
    call store_arrays( v%ms%Jmax, g )
    do ichbra = 1, this%channels%GetNumberChannels()
      chbra => this%channels%channels(ichbra)
      do ichket = 1, this%channels%GetNumberChannels()
        chket => this%channels%channels(ichket)
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        call to_3d_spin_basis_channel( this%MatCH(ichbra,ichket), v, chbra, chket, gidx_to_pidx )
      end do
    end do
    call release_arrays()
  end subroutine to_3d_spin_basis

  subroutine to_3d_spin_basis_channel(this, v, chbra, chket, gidx_to_pidx )
    use MyLibrary, only: pi, hc
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceMom
    class(NNCartesianChan), intent(inout) :: this
    type(NNForceMom), intent(in) :: v
    type(SpinPNChannel), intent(in), target :: chbra, chket
    integer, intent(in) :: gidx_to_pidx(:)
    type(MomentumGrids), pointer :: g
    integer :: bra, ket, ip_bra, ip_ket, sz1, sz2, sz3, sz4, Tz, ms1, ms2
    integer :: ich_lsj, j, s, lbra, lket, M, m1, m2, ibra, iket
    type(TwoBodyRelChanSpinMBasis), pointer :: ch_lsj
    real(8) :: fact, spin_part
    complex(8) :: r, orb_part

    g => chket%mom
    fact = 1.d0 / hc**3
    if( g%Lbox > 0.d0 ) fact = (2.d0*pi)**3 / g%Lbox**3
    sz1 = chbra%sz1
    sz2 = chbra%sz2
    sz3 = chket%sz1
    sz4 = chket%sz2
    Tz  = chket%Tz
    ms1 = (sz1+sz2)/2
    ms2 = (sz3+sz4)/2

    !$omp parallel
    !$omp do private(bra, ip_bra, ket, ip_ket, r, ich_lsj, ch_lsj, j, s, lbra, ibra, lket, iket, &
    !$omp &  orb_part, M, m1, m2, spin_part)
    do bra = 1, g%GetNumberGrids()
      ip_bra = gidx_to_pidx(bra)
      do ket = 1, g%GetNumberGrids()
        ip_ket = gidx_to_pidx(ket)

        r = 0.d0
        do ich_lsj = 1, v%ms%GetNumberChannels()
          ch_lsj => v%ms%GetChannel( ich_lsj )
          if( ch_lsj%GetZ() /= chket%Tz ) cycle
          j = ch_lsj%GetJ()
          s = ch_lsj%GetS()
          if( 2*s < abs( sz1+sz2 )) cycle
          if( 2*s < abs( sz3+sz4 )) cycle
          do lbra = abs(j-s), j+s
            if((-1)**lbra /= ch_lsj%GetParity() ) cycle
            ibra = ch_lsj%GetIndex(ip_bra,lbra)
            do lket = abs(j-s), j+s
              if((-1)**lket /= ch_lsj%GetParity() ) cycle
              iket = ch_lsj%GetIndex(ip_ket,lket)

              orb_part = 0.d0
              do M = -J, J
                m1 = M-ms1
                m2 = M-ms2
                if( abs(m1) > lbra ) cycle
                if( abs(m2) > lket ) cycle
                orb_part = orb_part + (cg_lsj%get(2*lbra,2*m1,2*S,2*ms1,2*J,2*M) * cg_lsj%get(2*lket,2*m2,2*S,2*ms2,2*J,2*M) ) * &
                    & ( sphs(m1, lbra, bra) * conjg( sphs(m2, lket, ket)) )
              end do
              spin_part = cg_spin%get(1,sz1,1,sz2,2*S,2*ms1) * cg_spin%get(1,sz3,1,sz4,2*S,2*ms2)
              r = r + v%MatCh( ich_lsj )%m(ibra,iket) * &
                  & spin_part * orb_part * &
                  & ( 1.d0 + abs(Tz) * (-1.d0)**( lket + s ) )

            end do
          end do
        end do
        this%m(bra,ket) = r * fact

      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine to_3d_spin_basis_channel

  !function angular_part( theta1, phi1, theta2, phi2, s1, s2, s3, s4, l1, l2, S, J ) result(r)
  function angular_part( idx1, idx2, s1, s2, s3, s4, l1, l2, S, J ) result(r)
    use MyLibrary, only: dcg, spherical_harmonics, spherical_harmonics_unnorm
    !real(8), intent(in) :: theta1, phi1, theta2, phi2
    integer, intent(in) :: idx1, idx2
    integer, intent(in) :: s1, s2, s3, s4, l1, l2, S, J
    real(8) :: spin_part
    complex(8) :: r, orb_part
    integer :: M, m1, m2, ms1, ms2
    ms1 = (s1+s2)/2
    ms2 = (s3+s4)/2
    orb_part = 0.d0
    do M = -J, J
      m1 = M-ms1
      m2 = M-ms2
      if( abs(m1) > l1 ) cycle
      if( abs(m2) > l2 ) cycle
      !orb_part = orb_part + ( dcg( 2*l1, 2*m1, 2*S, 2*ms1, 2*J, 2*M )*dcg( 2*l2, 2*m2, 2*S, 2*ms2, 2*J, 2*M ) ) * &
      !    & ( spherical_harmonics(l1, m1, cos(theta1), phi1)*conjg(spherical_harmonics(l2, m2, cos(theta2), phi2)) )
      orb_part = orb_part + (cg_lsj%get(2*l1,2*m1,2*S,2*ms1,2*J,2*M) * cg_lsj%get(2*l2,2*m2,2*S,2*ms2,2*J,2*M) ) * &
          & ( sphs(m1, l1, idx1) * conjg( sphs(m2, l2, idx2)) )
    end do
    !spin_part = dcg( 1, s1, 1, s2, 2*S, 2*ms1 ) * dcg( 1, s3, 1, s4, 2*S, 2*ms2 )
    spin_part = cg_spin%get(1,s1,1,s2,2*S,2*ms1) * cg_spin%get(1,s3,1,s4,2*S,2*ms2)
    r = spin_part * orb_part * (-1.d0)**((l2-l1)/2) ! Phase convention such that the r-space spherical harmonics is phaseless..
  end function angular_part

  subroutine store_arrays(Jmax, g)
    use MyLibrary, only: spherical_harmonics
    integer, intent(in) :: Jmax
    type(MomentumGrids), intent(in) :: g
    integer :: idx, l, m
    real(8) :: theta, phi
    call cg_spin%init(1, 1, .true., 1, 1, .true.)
    call cg_lsj%init( 0,2*Jmax+2, .false., 0,2, .false. )
    allocate( sphs(-Jmax-1:Jmax+1, 0:Jmax+1, g%GetNumberGrids() ))
    sphs(:,:,:) = 0.d0
    do idx = 1, g%GetNumberGrids()
      theta = g%grids(idx)%theta
      phi = g%grids(idx)%phi
      do l = 0, Jmax+1
        do m = -l, l
          sphs(m,l,idx) = spherical_harmonics(l, m, cos(theta), phi)
        end do
      end do
    end do
  end subroutine store_arrays

  subroutine release_arrays()
    deallocate( sphs )
    call cg_spin%fin()
    call cg_lsj%fin()
  end subroutine release_arrays

  subroutine find_momenta_interpolated( p, grid_to_p, grids )
    real(8), allocatable :: p(:)
    integer, allocatable :: grid_to_p(:)
    type(MomentumGrids), intent(in) :: grids
    integer :: idx, idx_1, cnt
    integer, allocatable :: i_tmp(:)
    real(8) :: px, py, pz, pp
    real(8), allocatable :: p_tmp(:)

    allocate( p_tmp( grids%GetNumberGrids() ))
    do idx = 1, grids%GetNumberGrids()
      px = grids%grids(idx)%px
      py = grids%grids(idx)%py
      pz = grids%grids(idx)%pz
      p_tmp(idx) = sqrt(px**2+py**2+pz**2)
    end do

    allocate( i_tmp( grids%GetNumberGrids() ))
    cnt = 0
    i_tmp(:) = 0
    do idx = 1, grids%GetNumberGrids()
      pp = minval(p_tmp)
      if(pp == 1.d100) cycle
      cnt = cnt+1
      do idx_1 = 1, grids%GetNumberGrids()
        if( p_tmp(idx_1) == 1.d100 ) cycle
        if( pp == p_tmp(idx_1) ) then
          i_tmp(idx_1) = cnt
          p_tmp(idx_1) = 1.d100
        end if
      end do
    end do
    allocate( p(cnt) )
    allocate( grid_to_p(grids%GetNumberGrids() ) )
    grid_to_p = i_tmp
    do idx = 1, cnt
      do idx_1 = 1, grids%GetNumberGrids()
        if( i_tmp(idx_1) == idx ) then
          px = grids%grids(idx_1)%px
          py = grids%grids(idx_1)%py
          pz = grids%grids(idx_1)%pz
          pp = sqrt(px**2+py**2+pz**2)
          p(idx) = pp
        end if
      end do
    end do

    deallocate(p_tmp, i_tmp)
  end subroutine find_momenta_interpolated

  subroutine interpolate( v, v_original )
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceMom
    type(NNForceMom), intent(inout) :: v
    type(NNForceMom), intent(in) :: v_original
    integer :: ich, j, p, s, tz, idx
    do ich = 1, v%ms%GetNumberChannels()
      j = v%ms%jpsz(ich)%GetJ()
      p = v%ms%jpsz(ich)%GetParity()
      s = v%ms%jpsz(ich)%GetS()
      tz= v%ms%jpsz(ich)%GetZ()
      idx = v_original%ms%GetIndex(j,p,s,tz)
      call interpolate_channel( v%MatCh(ich), v%ms%GetChannel(ich), v_original%MatCh(idx), v_original%ms%GetChannel(idx) )
      !call v_original%MatCh(idx)%prt("original")
      !call v%MatCh(idx)%prt("interpolation")
    end do
  end subroutine interpolate

  subroutine interpolate_channel( v, chan, v_original, chan_original )
    use TwoBodyRelativeSpace
    type(DMat), intent(inout) :: v
    type(DMat), intent(in) :: v_original
    type(TwoBodyRelChanSpinMBasis), intent(in) :: chan, chan_original
    type(DMat) :: v_sub_in, v_sub_out
    type(Mesh), pointer :: bra_m, ket_m
    integer :: NMesh_in, NMesh_out
    integer :: bra, ket, ibra, iket, lbra, lket
    logical :: coup

    NMesh_in = chan_original%GetNMesh()
    NMesh_out = chan%GetNMesh()
    coup = .false.
    if( chan_original%GetNumberStates() > NMesh_in ) coup=.true.
    call v_sub_in%ini( NMesh_in, NMesh_in )
    call v_sub_out%ini( NMesh_out, NMesh_out )
    ibra = 0
    do bra = 1, chan_original%GetNumberStates()
      bra_m => chan_original%GetP(bra)
      lbra = bra_m%GetL()
      if( coup .and. lbra> chan_original%GetJ() ) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, chan_original%GetNumberStates()
        ket_m => chan_original%GetP(ket)
        lket = ket_m%GetL()
        if( coup .and. lket> chan_original%GetJ() ) cycle
        iket = iket + 1
        v_sub_in%m(ibra,iket) = v_original%m(bra,ket)
        v_sub_in%m(iket,ibra) = v_sub_in%m(ibra,iket)
      end do
    end do
    call interpolate_submatrix( v_sub_out, chan, v_sub_in, chan_original )
    v%m(:NMesh_out, :NMesh_out) = v_sub_out%m(:,:)
    call v_sub_in%fin()
    call v_sub_out%fin()
    if(.not. coup) return

    ! L'=J+1, L=J+1
    call v_sub_in%ini( NMesh_in, NMesh_in )
    call v_sub_out%ini( NMesh_out, NMesh_out )
    ibra = 0
    do bra = 1, chan_original%GetNumberStates()
      bra_m => chan_original%GetP(bra)
      lbra = bra_m%GetL()
      if( coup .and. lbra< chan_original%GetJ() ) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, chan_original%GetNumberStates()
        ket_m => chan_original%GetP(ket)
        lket = ket_m%GetL()
        if( coup .and. lket< chan_original%GetJ() ) cycle
        iket = iket + 1
        v_sub_in%m(ibra,iket) = v_original%m(bra,ket)
        v_sub_in%m(iket,ibra) = v_sub_in%m(ibra,iket)
      end do
    end do
    call interpolate_submatrix( v_sub_out, chan, v_sub_in, chan_original )
    v%m(NMesh_out+1:, NMesh_out+1:) = v_sub_out%m(:,:)

    ! L'=J+1, L=J-1
    ibra = 0
    do bra = 1, chan_original%GetNumberStates()
      bra_m => chan_original%GetP(bra)
      lbra = bra_m%GetL()
      if( coup .and. lbra< chan_original%GetJ() ) cycle
      ibra = ibra + 1
      iket = 0
      do ket = 1, chan_original%GetNumberStates()
        ket_m => chan_original%GetP(ket)
        lket = ket_m%GetL()
        if( coup .and. lket> chan_original%GetJ() ) cycle
        iket = iket + 1
        v_sub_in%m(ibra,iket) = v_original%m(bra,ket)
      end do
    end do
    call interpolate_submatrix( v_sub_out, chan, v_sub_in, chan_original )
    v%m(NMesh_out+1:, :NMesh_out) = v_sub_out%m(:,:)
    v%m(:NMesh_out, NMesh_out+1:) = transpose(v_sub_out%m(:,:))
    call v_sub_in%fin()
    call v_sub_out%fin()
  end subroutine interpolate_channel

  subroutine interpolate_submatrix( v, chan, v_original, chan_original )
    use TwoBodyRelativeSpace
    use NdSpline
    type(DMat), intent(inout) :: v
    type(DMat), intent(in) :: v_original
    type(TwoBodyRelChanSpinMBasis), intent(in) :: chan, chan_original
    type(spline) :: sp
    integer, parameter :: k = 4
    integer :: i, n, nn
    real(8), allocatable :: f2din(:,:), f1d(:), pmom(:), p(:)

    n = size(v_original%m,1)
    allocate(f2din(n,n), f1d(n*n))
    allocate(pmom(n))
    do i = 1, n
      pmom(i) = chan_original%points(i)%GetP()
    end do
    f2din(:,:) = v_original%m(:,:)
    f1d = reshape(f2din, shape(f1d))
    call sp%init([k,k], [n,n], [pmom,pmom], f1d, extrapolation=.true.)
    deallocate(f1d)

    nn = size(v%m,1)
    allocate(p(nn))
    do i = 1, nn
      p(i) = chan%points(i)%GetP()
    end do
    call sp%interpolate([nn,nn], [p,p])
    v%m(:,:) = reshape(sp%fs, shape(v%m))

    deallocate(p,pmom)
    call sp%fin()
  end subroutine interpolate_submatrix

  subroutine test_one_pion_exchange(this)
    use MyLibrary, only: hc, pi
    type(NNCartesian), intent(inout), target :: this
    integer :: ichbra, ichket, bra, ket
    type(MomentumGrids), pointer :: g
    type(SpinPNChannel), pointer :: chbra, chket
    integer :: s1_bra, s2_bra, z1_bra, z2_bra
    integer :: s1_ket, s2_ket, z1_ket, z2_ket
    real(8) :: px_bra, py_bra, pz_bra
    real(8) :: px_ket, py_ket, pz_ket
    real(8) :: fact

    g => this%channels%mom
    fact = 1.d0 / hc**3
    if( g%Lbox > 0.d0 ) fact = (2.d0*pi)**3 / g%Lbox**3
    do ichbra = 1, this%channels%GetNumberChannels()
      chbra => this%channels%channels(ichbra)
      s1_bra = chbra%sz1
      s2_bra = chbra%sz2
      do ichket = 1, this%channels%GetNumberChannels()
        chket => this%channels%channels(ichket)
        s1_ket = chket%sz1
        s2_ket = chket%sz2
        if( chbra%Tz /= chket%Tz ) cycle

        if( chket%Tz == -1 ) then
          z1_bra =-1; z2_bra =-1; z1_ket =-1; z2_ket =-1
        elseif( chket%Tz == 1) then
          z1_bra = 1; z2_bra = 1; z1_ket = 1; z2_ket = 1
        else
          z1_bra =-1; z2_bra = 1; z1_ket =-1; z2_ket = 1
        end if

        do bra = 1, g%GetNumberGrids()
          px_bra = g%grids(bra)%px * hc
          py_bra = g%grids(bra)%py * hc
          pz_bra = g%grids(bra)%pz * hc
          do ket = 1, g%GetNumberGrids()
            px_ket = g%grids(ket)%px * hc
            py_ket = g%grids(ket)%py * hc
            pz_ket = g%grids(ket)%pz * hc

            this%MatCh(ichbra,ichket)%m(bra,ket) = &
                & ( func_ope( px_bra-px_ket, py_bra-py_ket, pz_bra-pz_ket, &
                & s1_bra, s2_bra, s1_ket, s2_ket, z1_bra, z2_bra, z1_ket, z2_ket ) - &
                & func_ope( px_bra+px_ket, py_bra+py_ket, pz_bra+pz_ket, &
                & s1_bra, s2_bra, s2_ket, s1_ket, z1_bra, z2_bra, z2_ket, z1_ket ) ) * fact

          end do
        end do
      end do
    end do
  end subroutine test_one_pion_exchange

  function func_ope( qx, qy, qz, s1, s2, s3, s4, t1, t2, t3, t4 ) result(r)
    use MyLibrary, only: g_A, f_pi, m_pi, hc, pi
    real(8), intent(in) :: qx, qy, qz
    integer, intent(in) :: s1, s2, s3, s4, t1, t2, t3, t4
    complex(8) :: r
    complex(8), parameter :: ei=(0.d0, 1.d0)
    complex(8) :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2), qdot(2,2)
    complex(8) :: chi1(2), chi2(2), chi3(2), chi4(2), qdot1, qdot2
    complex(8) :: tchi1(2), tchi2(2), tchi3(2), tchi4(2), tau_dot, sig_dot
    complex(8) :: sig_eye, tau_eye
    sigma_x(:,:) = 0.d0
    sigma_y(:,:) = 0.d0
    sigma_z(:,:) = 0.d0

    sigma_x(1,2) = 1.d0
    sigma_x(2,1) = 1.d0
    sigma_z(1,1) = 1.d0
    sigma_z(2,2) = -1.d0
    sigma_y(1,2) = -ei
    sigma_y(2,1) = ei

    tchi1 = 0.d0; tchi2=0.d0
    tchi3 = 0.d0; tchi4=0.d0
    tchi1( (1-t1)/2+1 ) = 1.d0
    tchi2( (1-t2)/2+1 ) = 1.d0
    tchi3( (1-t3)/2+1 ) = 1.d0
    tchi4( (1-t4)/2+1 ) = 1.d0

    chi1 = 0.d0; chi2=0.d0
    chi3 = 0.d0; chi4=0.d0
    chi1( (1-s1)/2+1 ) = 1.d0
    chi2( (1-s2)/2+1 ) = 1.d0
    chi3( (1-s3)/2+1 ) = 1.d0
    chi4( (1-s4)/2+1 ) = 1.d0

    tau_eye = dot_product(tchi1, tchi3) * dot_product(tchi2, tchi4)
    sig_eye = dot_product( chi1,  chi3) * dot_product( chi2,  chi4)

    tau_dot = dot_product(tchi1, matmul( sigma_x, tchi3) ) * dot_product(tchi2, matmul( sigma_x, tchi4) ) + &
        &     dot_product(tchi1, matmul( sigma_y, tchi3) ) * dot_product(tchi2, matmul( sigma_y, tchi4) ) + &
        &     dot_product(tchi1, matmul( sigma_z, tchi3) ) * dot_product(tchi2, matmul( sigma_z, tchi4) )
    sig_dot = dot_product(chi1, matmul( sigma_x, chi3) ) * dot_product(chi2, matmul( sigma_x, chi4) ) + &
        &     dot_product(chi1, matmul( sigma_y, chi3) ) * dot_product(chi2, matmul( sigma_y, chi4) ) + &
        &     dot_product(chi1, matmul( sigma_z, chi3) ) * dot_product(chi2, matmul( sigma_z, chi4) )

    qdot = qx*sigma_x + qy*sigma_y + qz*sigma_z
    qdot1 = dot_product(chi1, matmul( qdot, chi3) )
    qdot2 = dot_product(chi2, matmul( qdot, chi4) )

    r = - qdot1*qdot2 / ( qx**2+qy**2+qz**2 + m_pi**2 ) * &
        & g_A**2/(4.0D0*f_pi**2) * tau_dot * hc**3 / (2*pi)**3
    !r = ( qx**2+qy**2+qz**2) * sig_dot * tau_dot * &
    !    & g_A**2/(4.0D0*f_pi**2) * hc**3 / (2*pi)**3
    !r = qdot1*qdot2 * tau_eye * &
    !    & g_A**2/(4.0D0*f_pi**2) * hc**3 / (2*pi)**3
  end function func_ope

  subroutine transform_to_rspace( vmom, vr )
    use MyLibrary, only: pi, spherical_bessel
    use NNForce, only: NNForceMom
    use TwoBodyRelativeSpace
    type(NNForceMom), intent(in), target :: vmom
    type(NNForceMom), intent(inout) :: vr
    type(TwoBodyRelChanSpinMBasis), pointer :: ch_lsj_p, ch_lsj_r
    type(DMat) :: ovlap
    type(Mesh), pointer :: r, p
    integer :: ch, bra, ket
    real(8) :: prefact

    prefact = sqrt( 2.d0/pi )
    do ch = 1, vmom%ms%GetNumberChannels()
      ch_lsj_p => vmom%ms%GetChannel(ch)
      ch_lsj_r => vr%ms%GetChannel(ch)
      call ovlap%zeros( ch_lsj_p%GetNumberStates(), ch_lsj_r%GetNumberStates() )

      do bra = 1, ch_lsj_p%GetNumberStates()
        p => ch_lsj_p%GetP(bra)
        do ket = 1, ch_lsj_r%GetNumberStates()
          r => ch_lsj_r%GetP(ket)
          if( p%GetL() /= r%GetL() ) cycle
          ovlap%m(bra,ket) = prefact * p%GetP()**2 * p%GetW() * &
              & spherical_bessel( p%GetL(), p%GetP() * r%GetP() ) ! I do not include the phase i^(3l) here.
        end do
      end do
      vr%MatCh(ch) = ovlap%t() * vmom%MatCh(ch) * ovlap
      call ovlap%fin()
    end do

  end subroutine transform_to_rspace

  subroutine print_pw_MEs( v )
    use MyLibrary, only: hc, m_proton, m_neutron
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceMom
    type(NNForceMom), intent(inout), target :: v
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: ch_lsj
    type(Mesh), pointer :: bra_m, ket_m
    integer :: ich, bra, ket, wunit=15
    real(8) :: factor

    factor = 0.5d0 * (m_proton + m_neutron) / hc**2 ! factor to scattering unit
    open(wunit,file="pw_me.txt")
    ms => v%ms
    do ich = 1, ms%GetNumberChannels()
      ch_lsj => ms%GetChannel(ich)
      if( ch_lsj%GetJ() > 1 ) cycle
      do bra = 1, ch_lsj%GetNumberStates()
        bra_m => ch_lsj%GetP(bra)
        do ket = 1, ch_lsj%GetNumberStates()
          ket_m => ch_lsj%GetP(ket)
          write(wunit,'(5i3,3es16.5)') ch_lsj%GetS(), bra_m%GetL(), ket_m%GetL(), ch_lsj%GetJ(), ch_lsj%GetZ(), &
              & bra_m%GetP(), ket_m%GetP(), v%MatCh(ich)%m(bra,ket)*factor
        end do
      end do
    end do
    close(wunit)
  end subroutine print_pw_MEs

  subroutine fix_phase( v )
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceMom
    type(NNForceMom), intent(inout), target :: v
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: ch_lsj
    type(Mesh), pointer :: bra_m, ket_m
    integer :: ich, bra, ket
    ms => v%ms
    do ich = 1, ms%GetNumberChannels()
      ch_lsj => ms%GetChannel(ich)
      do bra = 1, ch_lsj%GetNumberStates()
        bra_m => ch_lsj%GetP(bra)
        do ket = 1, ch_lsj%GetNumberStates()
          ket_m => ch_lsj%GetP(ket)
          v%MatCh(ich)%m(bra,ket) = v%MatCh(ich)%m(bra,ket) * (-1.d0)**( (bra_m%GetL() - ket_m%GetL())/2 )
        end do
      end do
    end do
  end subroutine fix_phase
end module NNForceCartesian
