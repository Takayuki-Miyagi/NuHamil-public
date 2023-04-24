module NNForce
  use omp_lib
  use ClassSys
  use MPIFunction,only: myrank, ierr, nprocs
  use Profiler, only: timer
  use LinAlgLib
  use TwoBodyRelativeSpace
  implicit none

  public :: NNForceHO
  public :: NNForceMom
  public :: phase_shift_analysis

  private :: FinNNForceHO
  private :: InitNNForceHO
  private :: SetNNForceHO
  private :: add_coulomb_force
  private :: trans_mom_to_ho
  private :: set_ho_radial
  private :: WriteNNForceHORelative
  private :: ReadNNForceHORelative
  private :: GetFileNameBareNNForceRelative
  private :: check_nnforce_ho_space
  private :: renorm_ho_space
  private :: HOSRGChannel
  private :: LeeSuzukiChannel

  private :: FinNNForceMom
  private :: InitNNForceMom
  private :: SetNNForceMom
  private :: PhaseShiftAnalysis
  private :: PhaseShiftAnalysisChannel
  private :: calc_phase_shift_uncouple
  private :: calc_phase_shift_couple
  private :: calc_NNinteraction
  private :: calc_propagator
  private :: interpolation_vsub
  private :: check_nnforce_mom_space
  private :: renorm_mom_space
  private :: MomSRGChannel
  private :: MomVlowkChannel

  private :: WeightNormalize
  private :: WeightUnNormalize

  type :: NNForceMom
    type(DMat), allocatable :: MatCh(:)
    type(DMat), allocatable :: UT(:)
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    logical :: weighted = .false.
    logical :: is_init = .false.
  contains
    procedure :: InitNNForceMom
    procedure :: FinNNForceMom
    procedure :: SetNNForceMom
    procedure :: SetBareNNForceMom
    procedure :: WeightNormalize
    procedure :: WeightUnNormalize
    procedure :: PhaseShiftAnalysis
    ! overriding
    generic :: init => InitNNForceMom
    generic :: fin => FinNNForceMom
  end type NNForceMom

  type :: NNForceHO
    type(DMat), allocatable :: MatCh(:)
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: ms
    logical :: is_init = .false.
    type(str) :: NNint
  contains
    procedure :: InitNNForceHO
    procedure :: FinNNForceHO
    procedure :: SetNNForceHO
    procedure :: WriteNNForceHORelative
    procedure :: ReadNNForceHORelative
    ! overriding
    generic :: init => InitNNForceHO
    generic :: fin => FinNNForceHO
    generic :: writef => WriteNNForceHORelative
    generic :: readf => ReadNNForceHORelative
  end type NNForceHO

  real(8), private, allocatable :: radial(:,:,:,:)
  type, private :: v_inter
    real(8), allocatable :: v(:,:,:,:,:)
  end type v_inter
  type(sys), private :: sy
contains

  subroutine FinNNForceHO(this)
    class(NNForceHO), intent(inout) :: this
    integer :: ich
    do ich = 1, this%ms%GetNumberChannels()
      call this%MatCh(ich)%fin()
    end do
    deallocate(this%MatCh)
    this%ms => null()
    this%is_init = .false.
  end subroutine FinNNForceHO

  subroutine InitNNForceHO(this, ms)
    class(NNForceHO), intent(inout) :: this
    type(TwoBodyRelSpaceSpinHOBasis), target :: ms
    type(TwoBodyRelChanSpinHOBasis), pointer :: ms_ch
    integer :: ch

    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    allocate(this%MatCh(ms%GetNumberChannels() ))
    do ch = 1, ms%GetNumberChannels()
      ms_ch => ms%GetChannel(ch)
      call this%MatCh(ch)%zeros( ms_ch%GetNumberStates(), ms_ch%GetNumberStates() )
    end do
    this%is_init = .true.
  end subroutine InitNNForceHO

  subroutine SetNNForceHO(this, U, params, verbose, fn)
    use MPIFunction
    use MyLibrary, only: gauss_legendre
    use NuHamilInput, only: InputParameters
    class(NNForceHO), intent(inout) :: this
    type(NNForceHO), intent(inout) :: U
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: verbose
    type(str), intent(in), optional :: fn
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelSpaceSpinHOBasis) :: twol
    type(NNForceHO) :: vnn_large
    type(TwoBodyRelSpaceSpinMBasis) :: mom
    type(NNForceMom) :: vmom
    type(sys) :: s
    integer :: ch, npmesh, nqmesh
    real(8) :: ti
    real(8), allocatable :: p(:), w(:), pmom(:), wmom(:)
    type(str) :: f

    if(.not. this%is_init) then
      write(*,*) "Initialize NNForceHO before calling SetNNForceHO"
      return
    end if

    if(.not. U%is_init) then
      write(*,*) "Initialize NNForceHO before calling SetNNForceHO"
      return
    end if
    two => this%ms
    f = GetFileNameBareNNForceRelative(params%NNInt, params%N2max, two%GetFrequency(), params%path_to_NNrel)
    ! -------------------------------------------
    if( (.not. s%isfile(f)) .or. params%renorm_space2%val=="mom") then
      if(.not. s%isfile(f) .and. myrank==nprocs-1) write(*,'(2a)') 'Not found ', trim(f%val)
      if(myrank==nprocs-1) write(*,'(a)') 'Calculating NN interaction in relative coordinate:'
      if(myrank==nprocs-1) write(*,'(2a)',advance='no') 'NNint = ', trim(params%NNint%val)
      if(myrank==nprocs-1) write(*,'(a,i4)',advance='no') ', Nmax = ', params%N2max
      if(myrank==nprocs-1) write(*,'(a,i3)',advance='no') ', Jmax = ', params%J2max_NNint
      if(myrank==nprocs-1) write(*,'(a,f6.2,a)') ', hw = ', two%GetFrequency(), ' MeV'
      if(params%N2max /= two%GetNmax() .or. params%J2max_NNint /= two%GetJmax()) then
        write(*,'(a,2i4,a,2i3)') "There is a inconsistency in Nmax:", params%N2max, two%GetNmax(), &
            &" Jmax:", params%J2max_NNint, two%GetJmax()
        stop
      end if
      call twol%init(two%GetFrequency(), params%N2max, params%J2max_NNint)
      call vnn_large%init(twol)
      call mom%init(0.d0,params%pmax2,params%NMesh2,params%J2max_NNint)
      allocate(pmom(params%NMesh2), wmom(params%NMesh2))
      if(params%renorm%val == 'vlowk') then
        npmesh = int(dble(params%NMesh2) * params%lambda / params%pmax2)
        nqmesh = params%NMesh2 - npmesh
        call gauss_legendre(0.d0, params%lambda, p, w, npmesh)
        pmom(:npmesh) = p
        wmom(:npmesh) = w
        call gauss_legendre(params%lambda, params%pmax2, p, w, nqmesh)
        pmom(npmesh+1:) = p
        wmom(npmesh+1:) = w
        deallocate(p,w)
      else
        call gauss_legendre(0.d0, params%pmax2, p, w, params%NMesh2)
        pmom = p
        wmom = w
        deallocate(p,w)
      end if
      call mom%SetMeshWeight(pmom,wmom)
      deallocate(pmom,wmom)

      ti = omp_get_wtime()
      call vmom%init(mom)
      call vmom%SetNNForceMom(params, fn)

      allocate(radial(params%NMesh2, 0:twol%GetNmax()/2, 0:min(twol%GetJmax()+1, twol%GetNmax()), -1:1 ))
      call set_ho_radial(mom, twol%GetNmax(), min(twol%GetJmax()+1, twol%GetNmax()), twol%GetFrequency(), params%pn_same_mass )
      call trans_mom_to_ho(vnn_large,vmom)
      deallocate(radial)

      call mom%fin()
      call vmom%fin()

      if( params%renorm_space2%val == 'ho') then
        if(myrank == nprocs-1) call vnn_large%writef(f)
        do ch = 1, two%GetNumberChannels()
          this%MatCh(ch)%m(:,:) = vnn_large%MatCh(ch)%m(:,:)
        end do
      end if

      if( params%renorm_space2%val == 'mom') then
        do ch = 1, two%GetNumberChannels()
          this%MatCh(ch)%m(:,:) = vnn_large%MatCh(ch)%m(:,:)
        end do
        if(params%coul) call add_coulomb_force(this)
        return
      end if
      call vnn_large%fin()
      call twol%fin()
    else
      call this%readf(f)
    end if
    ! -------------------------------------------
    do ch = 1, two%GetNumberChannels()
      call U%MatCh(ch)%eye( two%jpsz(ch)%GetNumberStates() )
    end do

    this%NNInt = params%NNInt
#ifdef NNForceDebug
    call check_nnforce_ho_space(this, U, two%GetFrequency() )
#endif
    if(params%renorm%val /= 'bare' .and. params%renorm_space2%val == 'ho') then
      if(params%coul .and. params%evolve_coulomb) then
        call add_coulomb_force(this)
        if(present(verbose)) write(*,*) 'Coulomb is also evolved.'
      end if
      call renorm_ho_space(this, U, params%renorm, params%lambda, params%particle_rank, &
          & params%N_ls, params%pn_same_mass, params%srg_generator, params%Nmax_srg_edge)
      if(present(fn)) call write_nn_ho(this, fn)
#ifdef NNForceDebug
      call check_nnforce_ho_space(this, U, two%GetFrequency())
#endif
    end if
    if( params%coul ) then
      if( params%renorm%val=='bare' .or. params%renorm_space2%val=='mom' .or. (.not. params%evolve_coulomb) ) then
        call add_coulomb_force(this)
        if(present(verbose)) write(*,*) 'Coulomb is added.'
      else
        if(present(verbose)) write(*,*) 'Coulomb is supposed to be evolved.'
      end if
    end if
  end subroutine SetNNForceHO

  subroutine phase_shift_analysis(params)
    use MPIFunction
    use MyLibrary, only: gauss_legendre
    use NuHamilInput, only: InputParameters
    type(InputParameters), intent(in) :: params
    type(TwoBodyRelSpaceSpinMBasis) :: mom
    type(NNForceMom) :: vmom
    integer :: npmesh, nqmesh
    real(8), allocatable :: p(:), w(:), pmom(:), wmom(:)

    write(*,'(a)') 'Calculating NN interaction in relative coordinate:'
    write(*,'(2a)') 'NNint = ', trim(params%NNint%val)
    call mom%init(0.d0,params%pmax2,params%NMesh2,params%J2max_NNint)
    allocate(pmom(params%NMesh2), wmom(params%NMesh2))
    if(params%renorm%val == 'vlowk') then
      npmesh = int(dble(params%NMesh2) * params%lambda / params%pmax2)
      nqmesh = params%NMesh2 - npmesh
      call gauss_legendre(0.d0, params%lambda, p, w, npmesh)
      pmom(:npmesh) = p
      wmom(:npmesh) = w
      call gauss_legendre(params%lambda, params%pmax2, p, w, nqmesh)
      pmom(npmesh+1:) = p
      wmom(npmesh+1:) = w
      deallocate(p,w)
    else
      call gauss_legendre(0.d0, params%pmax2, p, w, params%NMesh2)
      pmom = p
      wmom = w
      deallocate(p,w)
    end if
    call mom%SetMeshWeight(pmom,wmom)
    deallocate(pmom,wmom)

    call vmom%init(mom)
    call vmom%SetNNForceMom(params)
    call vmom%PhaseShiftAnalysis(params)
    call mom%fin()
    call vmom%fin()
  end subroutine phase_shift_analysis

  subroutine add_coulomb_force(vnn)
    use MyLibrary, only: hc, m_proton, alpha, ho_radial_wf_norm, gauss_legendre
    type(NNForceHO), intent(inout) :: vnn
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelChanSpinHOBasis), pointer :: two_ch
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    real(8) :: hw
    integer :: Nmax, lmax
    integer, parameter :: nmesh = 500
    real(8), parameter :: rmax = 40.d0
    real(8), allocatable :: rcoo(:), wcoo(:)
    real(8), allocatable :: howf(:,:,:)
    real(8) :: a, r, factor, coul
    integer :: n, l, i, bra, ket, ich
    integer :: n1, l1, n2, l2

    two => vnn%ms
    hw = two%GetFrequency()
    Nmax = two%GetNmax()
    lmax = min(two%GetJmax()+1, two%GetNmax())
    factor = hc / alpha

    allocate(rcoo(nmesh), wcoo(nmesh))
    allocate(howf(nmesh, 0:nmax/2, 0:lmax))
    call gauss_legendre(0.d0, rmax, rcoo, wcoo, nmesh)
    a = 0.5d0 * m_proton * hw / hc ** 2
    do n = 0, nmax/2
      do l = 0, lmax
        do i = 1, nmesh
          r = rcoo(i)
          howf(i, n, l) = ho_radial_wf_norm(n,l,a,r)
        end do
      end do
    end do

    !$omp parallel
    !$omp do private(ich, two_ch, bra, ho_bra, n1, l1, ket, ho_ket, n2, l2, i, r, coul)
    do ich = 1, two%GetNumberChannels()
      two_ch => two%GetChannel(ich)
      if(two_ch%GetZ() /= -1) cycle
      do bra = 1, two_ch%GetNumberStates()
        ho_bra => two_ch%GetP(bra)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        do ket = 1, two_ch%GetNumberStates()
          ho_ket => two_ch%GetP(ket)
          n2 = ho_ket%GetN()
          l2 = ho_ket%GetL()

          if(l1 /= l2) cycle
          coul = 0.d0
          do i = 1, nmesh
            r = rcoo(i)
            coul = coul + howf(i, n1, l1) * howf(i, n2, l2) * wcoo(i) * &
                &           factor / r
          end do
          vnn%MatCh(ich)%m(bra, ket) = vnn%MatCh(ich)%m(bra, ket) + coul
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(rcoo, wcoo)
    deallocate(howf)
  end subroutine add_coulomb_force

  subroutine trans_mom_to_ho(vho, vmom)
    type(NNForceHO), intent(inout) :: vho
    type(NNFOrceMom), intent(in) :: vmom
    type(TwoBodyRelChanSpinMBasis), pointer :: ch_mesh
    type(TwoBodyRelChanSpinHOBasis), pointer :: ch_ho
    type(HarmonicOscillator), pointer :: ho_bra
    type(Mesh), pointer :: mesh_ket
    integer :: ich, nho, nmom, i, j, ich_mesh
    integer :: n, l, i1
    real(8) :: p, w
    type(DMat) :: ovlap, ovlapt

    do ich = 1, vho%ms%GetNumberChannels()
      ch_ho => vho%ms%GetChannel(ich)
      ich_mesh = vmom%ms%GetIndex( ch_ho%GetJ(), ch_ho%GetParity(), ch_ho%GetS(), ch_ho%GetZ() )
      ch_mesh => vmom%ms%GetChannel(ich_mesh)
      nmom = ch_mesh%GetNumberStates()
      nho = ch_ho%GetNumberStates()
      call ovlap%zeros(nmom, nho); call ovlapt%zeros(nho, nmom)

      do i = 1, nho
        ho_bra => ch_ho%GetP(i)
        n = ho_bra%GetN()
        l = ho_bra%GetL()
        do j = 1, nmom
          mesh_ket => ch_mesh%GetP(j)
          i1= mesh_ket%GetI()
          p = mesh_ket%GetP()
          w = mesh_ket%GetW()
          if(l /= mesh_ket%GetL() ) cycle
          ovlap%m(j, i) = p * w * &
              &           radial(i1, n, l, ch_ho%GetZ()) * (-1.d0) ** n
          ovlapt%m(i, j) = p * w * &
              &           radial(i1, n, l, ch_ho%GetZ()) * (-1.d0) ** n
        end do
      end do
      vho%MatCh(ich) = ovlapt * (vmom%MatCh(ich_mesh) * ovlap)
      call ovlap%fin(); call ovlapt%fin()
    end do
  end subroutine trans_mom_to_ho

  subroutine set_ho_radial(mom, Nmax, Lmax, hw, pn_same_mass)
    use MyLibrary, only: hc, m_proton, m_neutron, ho_radial_wf_norm
    type(TwoBodyRelSpaceSpinMBasis), intent(in) :: mom
    integer, intent(in) :: Nmax, Lmax
    real(8), intent(in) :: hw
    logical, intent(in) :: pn_same_mass
    integer :: NMesh, n, l, i, z
    real(8) :: a, m_red
    real(8), allocatable :: pmom(:)
    NMesh = mom%NMesh
    allocate(pmom(NMesh))
    do i = 1, NMesh
      pmom(i) = mom%jpsz(1)%points(i)%GetP()
    end do

    do z = -1, 1
      if(pn_same_mass) then
        m_red = (m_proton*m_neutron)/(m_proton+m_neutron)
      else
        if(z==-1) m_red = m_proton * 0.5d0
        if(z== 0) m_red = (m_proton*m_neutron) / (m_proton+m_neutron)
        if(z== 1) m_red = m_neutron * 0.5d0
      end if
      a = hc ** 2 / (m_red * hw)
      !$omp parallel
      !$omp do private(n, l, i) schedule(dynamic)
      do n = 0, nmax/2
        do l = 0, lmax
          do i = 1, nmesh
            radial(i, n, l, z) = ho_radial_wf_norm(n,l,a,pmom(i))
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
    deallocate(pmom)
  end subroutine set_ho_radial

  subroutine WriteNNForceHORelative(this, fn)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
    class(NNForceHO), intent(in) :: this
    type(str), intent(in) :: fn
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelChanSpinHOBasis), pointer :: two_ch
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: ich, j, p, s, z, ntot
    integer :: i1, i2
    integer :: n1, n2, l1, l2
    integer :: N2max, jmax
    type(c_ptr) :: f, err
    integer :: length
    character(kind=c_char,len=256) :: buffer
    real(8) :: time

    time = omp_get_wtime()
    two => this%ms
    N2max = two%GetNmax()
    jmax = two%GetJmax()
    ntot = 0
    do ich = 1, two%GetNumberChannels()
      two_ch => two%GetChannel(ich)
      do i1 = 1, two_ch%GetNumberStates()
        do i2 = 1, i1
          ntot = ntot + 1
        end do
      end do
    end do
    f = gzip_open(fn%val, "wt")
    write(buffer,*) ntot, N2max, jmax
    length = len_trim(buffer)
    err = gzip_writeline( f, trim(buffer), length)
    do ich = 1, two%GetNumberChannels()
      two_ch => two%GetChannel(ich)
      j = two_ch%GetJ()
      p = two_ch%GetParity()
      s = two_ch%GetS()
      z = two_ch%GetZ()
      do i1 = 1, two_ch%GetNumberStates()
        ho_bra => two_ch%GetP(i1)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        do i2 = 1, i1
          ho_ket => two_ch%GetP(i2)
          n2 = ho_ket%GetN()
          l2 = ho_ket%GetL()
          write(buffer, '(8i4,es18.6)') j, p, s, z, n1, l1, n2, l2, this%MatCh(ich)%m(i1,i2)
          length = len_trim(buffer)
          err = gzip_writeline( f, trim(buffer), length)
        end do
      end do
    end do
    err = gzip_close(f)
    call timer%add(sy%str('Write to file'), omp_get_wtime() - time)
  end subroutine WriteNNForceHORelative

  subroutine ReadNNForceHORelative(this, fn)
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceSpinHOBasis
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    class(NNForceHO), intent(inout) :: this
    type(str), intent(in) :: fn
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelChanSpinHOBasis), pointer :: two_ch
    integer :: ich, j, p, s, z, n1, n2, ntot, i, l1, l2
    integer :: i1, i2
    integer :: N2max, jmax
    integer :: N2max_r, jmax_r
    real(8) :: v, time
    type(c_ptr) :: f, err
    character(128) :: buffer

    time = omp_get_wtime()
    two => this%ms
    N2max = two%GetNmax()
    jmax = two%GetJmax()
    f = gzip_open(fn%val, "rt")
    buffer = ""
    err = gzip_readline( f, buffer, len(buffer) )
    read(buffer,*) ntot, N2max_r, jmax_r
    if(jmax_r < jmax) then
      write(*,*) "Error in ", __FILE__, " at ", __LINE__
      write(*,*) jmax_r, N2max_r, two%GetJmax(), two%GetNmax()
      stop
    end if

    do i = 1, ntot
      err = gzip_readline( f, buffer, len(buffer) )
      read(buffer,*) j, p, s, z, n1, l1, n2, l2, v
      if(2*n1 + l1 > N2max) cycle
      if(2*n2 + l2 > N2max) cycle
      if(j > Jmax) cycle
      ich = two%GetIndex(j,p,s,z)
      if(ich == 0) cycle
      two_ch => two%GetChannel(ich)
      i1 = two_ch%GetIndex(n1,l1)
      i2 = two_ch%GetIndex(n2,l2)
      this%MatCh(ich)%m(i1,i2) = v
      this%MatCh(ich)%m(i2,i1) = v
    end do
    err = gzip_close(f)
    call timer%add(sy%str('Read from file'), omp_get_wtime() - time)
  end subroutine ReadNNForceHORelative

  function GetFileNameBareNNForceRelative(NNInt,Nmax,hw,path_to_tmp_dir) result(r)
    type(str) :: r
    type(str), intent(in) :: NNInt
    integer, intent(in) :: Nmax
    real(8), intent(in) :: hw
    type(str), intent(in) :: path_to_tmp_dir
    type(sys) :: s
    call s%mkdir(path_to_tmp_dir%val)
    r = path_to_tmp_dir + s%str('/NNint_A2_rel_') + NNInt + s%str('_bare_Nmax') + s%str(Nmax) + &
        & s%str('_hw') + s%str(hw) + s%str(".gz")
  end function GetFileNameBareNNForceRelative

  subroutine check_nnforce_ho_space(vnn, U, hw)
    use MPIFunction, only: myrank
    use MyLibrary, only: hc, m_red_pn
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceSpinHOBasis
    type(NNForceHO), intent(in) :: vnn, U
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelChanSpinHOBasis), pointer :: two_ch
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    real(8), intent(in) :: hw
    integer :: n2max, j, p, z, s, ich, n, m
    integer :: bra, ket, n1, l1, n2, l2
    real(8) :: kine, r, f
    type(DMat) :: h, tkin, rm
    type(DVec) :: v
    type(EigenSolSymD) :: sol
    two => vnn%ms
    N2max = two%GetNmax()
    j = 1
    p = 1
    z = 0
    f = hc ** 2 / (m_red_pn * hw)
    do s = 0, 1
      two_ch => two%GetChannel(j,p,s,z)
      if(.not. associated(two_ch)) cycle
      ich = two%GetIndex(j,p,s,z)
      n = two_ch%GetNumberStates()
      call H%zeros(n,n); call tkin%zeros(n,n)
      call rm%zeros(n,n)
      do bra = 1, n
        ho_bra => two_ch%GetP(bra)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        if(2 * n1 + l1 > n2max) cycle
        do ket = 1, n
          ho_ket => two_ch%GetP(ket)
          n2 = ho_ket%GetN()
          l2 = ho_ket%GetL()
          if(2 * n2 + l2 > n2max) cycle
          kine = 0.d0
          if(l1 /= l2) cycle
          if(abs(n1 - n2) > 1) cycle
          if(n1 == n2) kine = dble(2 * n1 + l1) + 1.5d0
          if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
          if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          tkin%m(bra, ket) = kine * hw * 0.5d0

          r = 0.d0
          if(n1 == n2) r = dble(2 * n1 + l1) + 1.5d0
          if(n1 == n2-1) r = -dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
          if(n1 == n2+1) r = -dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          rm%m(bra, ket) = r * f
        end do
      end do
      H = tkin + vnn%MatCh(ich)
      rm = U%MatCh(ich)%T() * rm * U%MatCh(ich)
      call tkin%fin()
      call sol%init(h)
      call sol%DiagSym(h)
      if(myrank == 0) then
        m = min(5, n)
        call v%zeros(m)
        v%v(:) = sol%eig%v(:m)
        write(*,'(a)') 'd Eigen values in CheckNNForceHOSpace'
        write(*,'(a,i3,a,i3,a,i3,a,i4)') 'J = ', j, ',  P = ', p, ',  Tz = ', Z, ",  Nmax = ", two_ch%GetNmax()
        call v%prt()
        call v%fin()
        write(*, '(a, f12.6)') 'rm = ', &
            & sqrt(dot_product(sol%vec%m(:,1), matmul(rm%m, sol%vec%m(:,1))) / 4.d0)
      end if
      call sol%fin()
      call H%fin()
    end do
  end subroutine check_nnforce_ho_space

  subroutine renorm_ho_space(vnn, Trs, renorm, lambda, A, N_ls, pn_same_mass, generator, Nmax_srg_edge, verbose)
    use MyLibrary, only: hc, m_proton, m_neutron, m_red_pn
    class(NNForceHO), intent(inout) :: vnn
    type(NNForceHO), intent(inout) :: Trs
    type(str), intent(in) :: renorm
    real(8), intent(in) :: lambda
    integer, intent(in) :: A, N_ls, Nmax_srg_edge
    type(str), intent(in) :: generator
    logical, intent(in) :: pn_same_mass
    logical, intent(in), optional :: verbose
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    integer :: ich, Nmax
    real(8) :: lc, alpha, hw, m_red

    if(myrank == 0) then
      if(present(verbose)) write(*,'(a)') '## NN Renormalization in HO space ...'
    end if
    two => vnn%ms
    Nmax = two%GetNmax()
    hw = two%GetFrequency()
    if(renorm%val == 'srg') then
      do ich = 1, two%GetNumberChannels()
        if(pn_same_mass) then
          m_red = m_red_pn
        else
          if(two%jpsz(ich)%GetZ() ==-1) m_red = m_proton * 0.5d0
          if(two%jpsz(ich)%GetZ() == 0) m_red = m_red_pn
          if(two%jpsz(ich)%GetZ() == 1) m_red = m_neutron * 0.5d0
        end if
        lc = lambda * hc / dsqrt(2.d0 * m_red)
        alpha = (1.d0 / lc) ** (4.d0)
        call HOSRGChannel(vnn%MatCh(ich),Trs%MatCh(ich),vnn%ms%jpsz(ich),alpha,hw,generator,Nmax_srg_edge)
      end do
    elseif(renorm%val == 'ls') then
      do ich = 1, two%GetNumberChannels()
        call LeeSuzukiChannel(vnn%MatCh(ich),Trs%MatCh(ich),vnn%ms%jpsz(ich),hw, N_ls, A)
      end do
    end if
  end subroutine renorm_ho_space

  subroutine HOSRGChannel(v, u, ms, alpha, hw, generator_type, Nmax_srg_edge)
    use OperatorDefinitions, only: CalcMERel, OperatorDef
    use Renormalization, only: SRGSolver
    type(DMat), intent(inout) :: v, u
    type(TwoBodyRelChanSpinHOBasis), intent(in) :: ms
    real(8), intent(inout) :: alpha, hw
    type(str), intent(in) :: generator_type
    integer, intent(in) :: Nmax_srg_edge
    integer :: n, n_zero
    integer :: bra, ket, n1, l1, n2, l2
    real(8) :: kine
    integer, allocatable :: rhs_zero_term(:)
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(DMat) :: h, t, generator
    type(SRGSolver) :: sol
    type(OperatorDef) :: Gen

    n = ms%GetNumberStates()
    if(n < 1) return
    call t%zeros(n,n)
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        ho_ket => ms%GetP(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        kine = 0.d0
        if(l1 /= l2) cycle
        if(abs(n1 - n2) > 1) cycle
        if(n1 == n2) kine = dble(2 * n1 + l1) + 1.5d0
        if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
        if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
        t%m(bra, ket) = kine * hw * 0.5d0
      end do
    end do

    ! Generator
    call generator%zeros(n,n)
    call Gen%InitOpDef(generator_type, .true.)
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        ho_ket => ms%GetP(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        generator%m(bra,ket) = CalcMERel(Gen, [n1,l1,ms%GetS(),ms%GetJ(),ms%GetZ()], [n2,l2,ms%GetS(),ms%GetJ(),ms%GetZ()])
      end do
    end do

    ! Neumann boundary for SRG flow
    n_zero = 0
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if( 2*n1+l1 < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
    end do
    allocate(rhs_zero_term(n_zero))
    n_zero = 0
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if( 2*n1+l1 < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
      rhs_zero_term(n_zero) = bra
    end do

    h = t + v
    !call sol%init(h, 'HUflow')
    call sol%init(h, 'Hflow')
    if( n_zero > 0 ) call sol%SRGFlow(h, generator, alpha, rhs_zero_term)
    if( n_zero == 0) call sol%SRGFlow(h, generator, alpha)
    v = sol%h - t
    u = sol%U

    call sol%fin()
    call generator%fin()
    call t%fin()
    call h%fin()
  end subroutine HOSRGChannel

  subroutine LeeSuzukiChannel(vin, u, ms, hw, N_ls, A)
    use Renormalization, only: LeeSuzukiSolver
    type(DMat), intent(inout) :: vin, u
    type(TwoBodyRelChanSpinHOBasis), intent(in) :: ms
    real(8), intent(in) :: hw
    integer, intent(in) :: N_ls, A
    integer :: n
    integer :: n1, l1, bra
    integer :: n2, l2, ket
    integer :: np, nq
    integer :: NN, cnt, cntb, cntk
    integer, allocatable :: reord(:), back(:)
    real(8) :: rpot
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(DMat) :: h, hd, v, ho
    type(LeeSuzukiSolver) :: sol

    n = ms%GetNumberStates()
    if(n < 1) return
    allocate(reord(n), back(n))
    cnt = 0
    do NN = 0, ms%GetNmax()
      do bra = 1, n
        ho_bra => ms%GetP(bra)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        if(NN /= 2 * n1 + l1) cycle
        cnt = cnt + 1
        reord(bra) = cnt
        back(cnt) = bra
      end do
    end do
    np = 0; nq = 0
    do bra = 1, n
      cnt = reord(bra)
      ho_bra => ms%GetP(cnt)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if(2 * n1 + l1 <= N_ls) np = np + 1
      if(2 * n1 + l1 >  N_ls) nq = nq + 1
    end do

    call hd%zeros(n,n); call h%zeros(n,n); call v%zeros(n,n)
    call ho%zeros(n,n)

    do bra = 1, n
      cntb = reord(bra)
      ho_bra => ms%GetP(cntb)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        cntk = reord(ket)
        ho_ket => ms%GetP(cntk)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        rpot = 0.d0
        if(l1 /= l2) then
          if(abs(n1 - n2) > 1) then
            if(n1 == n2) rpot = dble(2 * n1 + l1) + 1.5d0
            if(n1 == n2-1) rpot = -dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
            if(n1 == n2+1) rpot = -dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          end if
        end if
        ho%m(bra, ket) = - rpot * hw / dble(A)
        v%m(bra, ket) = vin%m(cntb, cntk)
      end do
      hd%m(bra, bra) = dble(2 * n1 + l1) * hw
    end do

    h = hd + v + ho
    call sol%init(h)
    call sol%LeeSuzuki(0, np, nq, h)
    v = sol%h - hd - ho
    do bra = 1, n
      cntb = back(bra)
      do ket = 1, n
        cntk = back(ket)
        vin%m(bra,ket) = v%m(cntb, cntk)
        U%m(bra,ket) = sol%expS%m(cntb, cntk)
      end do
    end do
    call sol%fin()
    call h%fin(); call hd%fin(); call v%fin()
    deallocate(reord, back)
  end subroutine LeeSuzukiChannel


  ! Methods for NN force in momentum space
  subroutine FinNNForceMom(this)
    class(NNForceMom), intent(inout) :: this
    integer :: ich
    do ich = 1, this%ms%NChan
      call this%MatCh(ich)%fin()
      call this%UT(ich)%fin()
    end do
    deallocate(this%MatCh)
    deallocate(this%UT)
    this%ms => null()
    this%weighted = .false.
  end subroutine FinNNForceMom

  subroutine InitNNForceMom(this, ms)
    class(NNForceMom), intent(inout) :: this
    type(TwoBodyRelSpaceSpinMBasis), target :: ms
    integer :: ch, n

    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    allocate(this%MatCh(ms%NChan))
    allocate(this%UT(ms%NChan))
    do ch = 1, ms%NChan
      n = ms%jpsz(ch)%GetNumberStates()
      call this%MatCh(ch)%zeros(n,n)
      call this%UT(ch)%eye(n)
    end do
    this%is_init = .true.
    this%weighted = .false.
  end subroutine InitNNForceMom

  subroutine WeightNormalize(this)
    class(NNForceMom), intent(inout) :: this
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: ch
    integer :: ich, bra, ket
    type(Mesh), pointer :: bra_m, ket_m
    real(8) :: p1, w1, p2, w2

    if(this%weighted) then
      write(*,*) "The potential is already normalized..."
      return
    end if
    ms => this%ms
    do ich = 1, ms%GetNumberChannels()
      ch => ms%GetChannel(ich)
      !$omp parallel
      !$omp do private(bra, bra_m, p1, w1, ket, ket_m, p2, w2)
      do bra = 1, ch%GetNumberStates()
        bra_m => ch%GetP(bra)
        p1 = bra_m%GetP()
        w1 = bra_m%GetW()
        do ket = 1, ch%GetNumberStates()
          ket_m => ch%GetP(ket)
          p2 = ket_m%GetP()
          w2 = ket_m%GetW()
          this%MatCh(ich)%m(bra,ket) = this%MatCh(ich)%m(bra,ket) * p1 * p2 * sqrt(w1*w2)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
    this%weighted=.true.
  end subroutine WeightNormalize

  subroutine WeightUnNormalize(this)
    class(NNForceMom), intent(inout) :: this
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: ch
    integer :: ich, bra, ket
    type(Mesh), pointer :: bra_m, ket_m
    real(8) :: p1, w1, p2, w2

    if(.not. this%weighted) then
      write(*,*) "The potential is already unnormalized..."
      return
    end if
    ms => this%ms
    do ich = 1, ms%GetNumberChannels()
      ch => ms%GetChannel(ich)
      !$omp parallel
      !$omp do private(bra, bra_m, p1, w1, ket, ket_m, p2, w2)
      do bra = 1, ch%GetNumberStates()
        bra_m => ch%GetP(bra)
        p1 = bra_m%GetP()
        w1 = bra_m%GetW()
        do ket = 1, ch%GetNumberStates()
          ket_m => ch%GetP(ket)
          p2 = ket_m%GetP()
          w2 = ket_m%GetW()
          this%MatCh(ich)%m(bra,ket) = this%MatCh(ich)%m(bra,ket) / (p1 * p2 * sqrt(w1*w2))
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
    this%weighted=.false.
  end subroutine WeightUnNormalize

  subroutine SetBareNNForceMom(this, NNInt)
    class(NNForceMom), intent(inout) :: this
    type(str), intent(in) :: NNInt
    call read_nn_mom(this, NNInt)
  end subroutine SetBareNNForceMom

  subroutine SetNNForceMom(v, params, fn)
    use NuHamilInput, only: InputParameters
    class(NNForceMom), intent(inout) :: v
    type(InputParameters), intent(in) :: params
    type(str), intent(in), optional :: fn

    if(.not. v%is_init) then
      write(*,*) "Initialize NNForceHO before calling SetNNForceHO"
      return
    end if

    call v%SetBareNNForceMom(params%input_nn_file)
    if(present(fn) .and. params%renorm%val == "bare") then
      call write_nn_mom(v, fn)
    end if
#ifdef NNForceDebug
    call check_nnforce_mom_space(v)
#endif
    if(params%renorm%val /= "bare" .and. params%renorm_space2%val == "mom") then
      call renorm_mom_space(v, params%renorm, params%lambda, params%pn_same_mass)
      if(present(fn)) call write_nn_mom(v, fn)

#ifdef NNForceDebug
      call check_nnforce_mom_space(v)
#endif
    end if
  end subroutine SetNNForceMom

  subroutine PhaseShiftAnalysis(vmom, params)
    use NuHamilInput, only: InputParameters
    class(NNForceMom), intent(in) :: vmom
    type(NNForceMom) :: op, op1, op2
    type(InputParameters), intent(in) :: params
    integer :: ich
    integer :: wunit = 20
    real(8) :: time
    type(str) :: fn_coupled, fn_uncoupled
    type(sys) :: s

    fn_coupled = params%file_name_phase_shift + '_coupled.txt'
    fn_uncoupled = params%file_name_phase_shift + '_uncoupled.txt'
    call op%init(vmom%ms)
    do ich = 1, vmom%ms%GetNumberChannels()
      op%MatCh(ich) = vmom%MatCh(ich)
    end do
    if(params%spin_tensor_decomposition/=-1) then
      select case(params%spin_tensor_decomposition)
      case(0, 1, 2)
        op = spin_tensor_decomposition(vmom, 0)
      case(3)
        op  = spin_tensor_decomposition(vmom, 0)
        op1 = spin_tensor_decomposition(vmom, 1)
        do ich = 1, vmom%ms%GetNumberChannels()
          op%MatCh(ich) = op%MatCh(ich) + op1%MatCh(ich)
        end do
      case(4)
        op  = spin_tensor_decomposition(vmom, 0)
        op2 = spin_tensor_decomposition(vmom, 2)
        do ich = 1, vmom%ms%GetNumberChannels()
          op%MatCh(ich) = op%MatCh(ich) + op2%MatCh(ich)
        end do
      case(-100)
        op  = spin_tensor_decomposition(vmom, 0)
        op1 = spin_tensor_decomposition(vmom, 1)
        op2 = spin_tensor_decomposition(vmom, 2)
        do ich = 1, vmom%ms%GetNumberChannels()
          op%MatCh(ich) = op%MatCh(ich) + 1.1d0*op1%MatCh(ich) + op2%MatCh(ich)
        end do
      end select
    end if
    call check_nnforce_mom_space(op)
    time = omp_get_wtime()
    open(wunit, file = fn_coupled%val, status='replace')
    write(wunit,*) "j, s, p, tz, q, Tlab, delta--, delta++, eps"
    close(wunit)
    open(wunit, file = fn_uncoupled%val, status='replace')
    write(wunit,*) "j, s, p, tz, q, Tlab, delta"
    close(wunit)
    do ich = 1, op%ms%NChan
      call PhaseshiftAnalysisChannel(op%MatCh(ich), op%ms%jpsz(ich), fn_coupled, fn_uncoupled, params%Tlab_phase_shift)
    end do
    call timer%add(sy%str('Phase shift analysis'), omp_get_wtime() - time)
  end subroutine PhaseShiftAnalysis

  subroutine PhaseShiftAnalysisChannel(vint,ms,fn_coupled,fn_uncoupled,T_labs)
    !
    ! phase-shift analysis using K-matrix
    ! K(q',q) = V(q',q) + M P int dk k^2 V(q',k) K(k,q) / (q^2 - k^2)
    ! In the Matrix form, K = V + VGK <=> K = (1-VG)^{-1} V
    ! V: NN interaction
    ! G: propagator
    ! Note that T- and K-matricies are related with
    ! T(E) = K(E) - i \pi K(E) \delta(E-E_0) T(E)
    !
    use MyLibrary, only: hc, m_proton, m_neutron, pi, m_nucleon
    type(DMat), intent(in) :: vint
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(str), intent(in) :: fn_coupled, fn_uncoupled
    real(8), intent(in) :: T_labs(:)
    integer :: j, p, s, tz, i0, n, m
    real(8) :: q0, mu, Tlab, delta, deltas(3)
    type(DMat) :: K, G, V, E
    logical :: couple
    integer :: wunit = 20

    j = ms%GetJ()
    p = ms%GetParity()
    s = ms%GetS()
    tz= ms%GetZ()

    mu = 0.d0
    if(tz == -1) mu = 0.5d0 * m_proton
    if(tz ==  0) mu = m_proton * m_neutron / (m_proton + m_neutron)
    if(tz ==  1) mu = 0.5d0 * m_neutron

    if(ms%GetNumberStates() /= ms%GetNMesh() ) then
      couple = .true.
    else
      couple = .false.
    end if

    m = 0
    if(couple) m = ms%GetNumberStates() + 2
    if(.not. couple) m = ms%GetNumberStates() + 1
    n = ms%GetNumberStates()
    if(couple) open(wunit, file = fn_coupled%val, position='append')
    if(.not. couple) open(wunit, file = fn_uncoupled%val, position='append')
    do i0 = 1, size(T_labs)
      Tlab = T_labs(i0)

      q0 = sqrt(m_proton**2 * Tlab * (Tlab + 2.d0 * m_neutron) / ( (m_proton + m_neutron)**2 + 2.d0 * Tlab * m_proton)) / hc
      if(tz == -1) q0 = sqrt(0.5d0 * m_proton * Tlab) / hc
      if(tz ==  0) q0 = sqrt(m_proton**2 * Tlab * (Tlab + 2.d0 * m_neutron) / ( (2.d0*m_nucleon)**2 + 2.d0 * Tlab * m_proton)) / hc
      if(tz ==  1) q0 = sqrt(0.5d0 * m_neutron * Tlab) / hc

      call G%zeros(m,m)
      call V%zeros(m,m)

      call calc_NNinteraction(vint,ms,q0,V,couple)
      V = V / hc**3 ! unit of MeV^{-2}
      call calc_propagator(ms,q0,mu,G,couple)

      call E%eye(m)
      E = E - (V*G)
      K = E%inv() * V

      if(.not. couple) then
        delta = calc_phase_shift_uncouple(ms,K,mu,q0)
        write(wunit,'(4i4,3f14.8)') j, s, p, tz, q0, Tlab, delta
      end if

      if(couple) then
        deltas = calc_phase_shift_couple(ms,K,mu,q0)
        write(wunit,'(4i4,5f14.8)') j, s, p, tz, q0, Tlab, deltas(:)
      end if

      call G%fin()
      call V%fin()
      call K%fin()
      call E%fin()

    end do
    close(wunit)
  end subroutine PhaseShiftAnalysisChannel

  function calc_phase_shift_uncouple(ms,K,mu,q0) result(delta)
    ! tan (delta) = - pi q M T(q, q) / 2
    use MyLibrary, only: hc, pi
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(DMat), intent(in) :: K
    real(8), intent(in) :: mu, q0
    real(8) :: m
    real(8) :: delta, KL
    m = 2.d0 * mu
    KL = pi * q0 * hc * m * K%m(ms%GetNMesh()+1, ms%GetNMesh()+1)
    delta = atan( - KL * 0.5d0)
  end function calc_phase_shift_uncouple

  function calc_phase_shift_couple(ms,K,mu,q0) result(deltas)
    ! for spin triplet coupled state:
    !        tan (delta_{-/+}) = - pi q M [K^{j-1,j-1} + K^{j+1,j+1} +/- (K^{j-1,j-1} - K^{j+1,j+1} / cos(2 eps))] / 4
    !        tan ( 2 eps     ) =  2 K^{j+1,j-1} / (K^{j-1,j-1} - K^{j+1,j+1})
    !        Blatt and Biedenharn convention
    !        bar{delta}_{+} + bar{delta}_{-} = delta_{+} + delta_{-}
    !        sin (bar{delta}_{+} - bar{delta}_{-} ) = tan (2 bar{eps}) / tan (2 eps)
    !        sin (    delta_{+}  -     delta_{-}  ) = sin (2 bar{eps}) / sin (2 eps)
    !        Stapp convention (usually used)
    use MyLibrary, only: hc, pi
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(DMat), intent(in) :: K
    real(8), intent(in) :: mu, q0
    real(8) :: m
    real(8) :: deltas(3)
    real(8) :: Kpp, Kmm, Kpm
    real(8) :: delta_p, delta_m, eps
    real(8) :: delta_bar_p, delta_bar_m, eps_bar
    logical :: filp
    m = 2.d0 * mu
    Kmm = pi * q0 * hc * m * real(K%m(  ms%GetNMesh()+1,   ms%GetNMesh()+1))
    Kpp = pi * q0 * hc * m * real(K%m(2*ms%GetNMesh()+2, 2*ms%GetNMesh()+2))
    Kpm = pi * q0 * hc * m * real(K%m(2*ms%GetNMesh()+2,   ms%GetNMesh()+1))

    eps = atan(2.d0*Kpm / (Kmm - Kpp)) * 0.5d0
    delta_m = atan( - 0.25d0 * ( (Kmm + Kpp) + (Kmm - Kpp) / (cos( 2.d0 * eps ) ) ) )
    delta_p = atan( - 0.25d0 * ( (Kmm + Kpp) - (Kmm - Kpp) / (cos( 2.d0 * eps ) ) ) )

    eps_bar = 0.5d0 * asin( sin(2.d0 * eps) * sin(delta_m - delta_p) )
    if(abs(eps) > 1.d-16) then
      delta_bar_m = 0.5d0 * (delta_p + delta_m) + 0.5d0 * asin( tan(2.d0*eps_bar) / tan(2.d0*eps) )
      delta_bar_p = 0.5d0 * (delta_p + delta_m) - 0.5d0 * asin( tan(2.d0*eps_bar) / tan(2.d0*eps) )
    else
      delta_bar_m = delta_m
      delta_bar_p = delta_p
    end if

    deltas(1) = delta_bar_m
    deltas(2) = delta_bar_p
    deltas(3) = eps_bar
  end function calc_phase_shift_couple

  subroutine calc_NNinteraction(vint,ms,q0,V,couple)
    use MyLibrary, only: hc
    type(DMat), intent(in) :: vint
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(Mesh), pointer :: bra_m, ket_m
    real(8), intent(in) :: q0
    type(DMat), intent(inout) :: V
    logical, intent(in) :: couple
    integer :: bra, ket, ibra, iket
    integer :: lbra, lket
    integer :: lm, lp
    type(DMat) :: v_sub

    call v_sub%ini(ms%GetNMesh()+1, ms%GetNMesh()+1)
    lm = abs(ms%GetJ() - ms%GetS())
    lp =     ms%GetJ() + ms%GetS()

    ibra = 0
    do bra = 1, ms%GetNumberStates()
      bra_m => ms%GetP(bra)
      lbra = bra_m%GetL()
      if(couple .and. lbra /= lm) cycle
      ibra = ibra + 1

      iket = 0
      do ket = 1, bra
        ket_m => ms%GetP(ket)
        lket = ket_m%GetL()
        if(couple .and. lket /= lm) cycle
        iket = iket + 1

        v_sub%m(ibra,iket) = vint%m(bra,ket)
        v_sub%m(iket,ibra) = v_sub%m(ibra,iket)

      end do
    end do

    call interpolation_vsub(v_sub,ms,q0)
    V%m(:ms%GetNMesh()+1, :ms%GetNMesh()+1) = v_sub%m(:,:)
    call v_sub%fin()

    if(.not. couple) return
    call v_sub%ini(ms%GetNMesh()+1, ms%GetNMesh()+1)
    ibra = 0
    do bra = 1, ms%GetNumberStates()
      bra_m => ms%GetP(bra)
      lbra = bra_m%GetL()
      if(lbra /= lp) cycle
      ibra = ibra + 1

      iket = 0
      do ket = 1, bra
        ket_m => ms%GetP(ket)
        lket = ket_m%GetL()
        if(lket /= lp) cycle
        iket = iket + 1

        v_sub%m(ibra,iket) = vint%m(bra,ket)
        v_sub%m(iket,ibra) = v_sub%m(ibra,iket)

      end do
    end do

    call interpolation_vsub(v_sub,ms,q0)
    V%m(ms%GetNMesh()+2:, ms%GetNMesh()+2:) = v_sub%m(:,:)

    ibra = 0
    do bra = 1, ms%GetNumberStates()
      bra_m => ms%GetP(bra)
      lbra = bra_m%GetL()
      if(lbra /= lm) cycle
      ibra = ibra + 1

      iket = 0
      do ket = 1, ms%GetNumberStates()
        ket_m => ms%GetP(ket)
        lket = ket_m%GetL()
        if(lket /= lp) cycle
        iket = iket + 1

        v_sub%m(ibra,iket) = vint%m(bra,ket)

      end do
    end do

    call interpolation_vsub(v_sub,ms,q0)
    V%m(:ms%GetNMesh()+1, ms%GetNMesh()+2:) = v_sub%m(:,:)
    V%m(ms%GetNMesh()+2:, :ms%GetNMesh()+1) = transpose(v_sub%m(:,:))
    call v_sub%fin()

  end subroutine calc_NNinteraction

  subroutine calc_propagator(ms,q0,mu,G,couple)
    use MyLibrary, only: hc
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(Mesh), pointer :: point
    real(8), intent(in) :: q0, mu
    type(DMat), intent(inout) :: G
    logical, intent(in) :: couple
    integer :: ket, n
    real(8) :: q, w
    ! complex(8) :: i_cmplx = (0.d0,1.d0)

    n = ms%GetNMesh()
    do ket = 1, n
      point => ms%GetP(ket)
      q = point%GetP()
      w = point%GetW()
      G%m(ket,ket) = 2.d0 * mu * w * hc * q**2 / ( q0**2 - q**2)
    end do

    do ket = 1, n
      point => ms%GetP(ket)
      q = point%GetP()
      w = point%GetW()
      G%m(n+1,n+1) = G%m(n+1,n+1) - 2.d0 * mu * w * hc * q0**2 / ( q0**2 - q**2)
    end do
    ! G%m(n+1,n+1) = G%m(n+1,n+1) + pi * mu * q0 * hc * i_cmplx ! if you want to calculate T-matrix

    if(.not. couple) return
    do ket = 1, ms%GetNMesh()+1
      G%m(ket + ms%GetNMesh()+1, ket+ms%GetNMesh()+1) = G%m(ket,ket)
    end do
  end subroutine calc_propagator

  subroutine interpolation_vsub(v, ms, q0)
    use NdSpline
    type(DMat), intent(inout) :: v
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    real(8), intent(in) :: q0
    type(spline) :: sp
    integer, parameter :: k = 4
    real(8), allocatable :: f2din(:,:), f1d(:), pmom(:)
    integer :: i, n

    n = size(v%m,1) - 1
    allocate(f2din(n,n), f1d(n*n))
    allocate(pmom(n))
    do i = 1, n
      pmom(i) = ms%points(i)%GetP()
    end do


    f2din(:,:) = v%m(:n, :n)
    f1d = reshape(f2din, shape(f1d))
    call sp%init([k,k], [n,n], [pmom,pmom], f1d)
    deallocate(f1d)

    call sp%interpolate([1,1], [q0, q0])
    v%m(n+1,n+1) = sp%fs(1)
    call sp%interpolate([n,1], [pmom, q0])
    v%m(1:n,n+1) = sp%fs(:)
    call sp%interpolate([1,n], [q0, pmom])
    v%m(n+1,1:n) = sp%fs(:)

    call sp%fin()

    deallocate(f2din,pmom)
  end subroutine interpolation_vsub


  subroutine check_nnforce_mom_space(vmom)
    use MPIFunction, only: myrank
    use MyLibrary, only: m_proton, m_neutron, hc
    type(NNForceMom), intent(inout) :: vmom
    type(TwoBodyRelSpaceSpinMBasis), pointer :: mom
    type(TwoBodyRelChanSpinMBasis), pointer :: mom_ch
    type(Mesh), pointer :: bra_m
    type(DMat) :: h
    type(DVec) :: v
    real(8) :: rmass
    integer :: ich, n, bra, ket, m
    real(8) :: pbra
    integer :: j, p, s, t
    type(EigenSolSymD) :: sol

    mom => vmom%ms
    rmass = m_proton * m_neutron / (m_proton + m_neutron)
    call vmom%WeightNormalize()
    j = 1
    p = 1
    t = 0
    do s = 0, 1
      ich = mom%jpsz2idx(j, p, s, t)
      if(ich < 1) cycle
      mom_ch => mom%GetChannel(ich)
      n = mom_ch%GetNumberStates()
      call H%zeros(n,n)
      do bra = 1, n
        bra_m => mom_ch%GetP(bra)
        pbra = bra_m%GetP()
        do ket = 1, bra
          H%m(bra, ket) = vmom%MatCh(ich)%m(bra, ket)
          H%m(ket, bra) = H%m(bra, ket)
        end do
        H%m(bra, bra) = (pbra * hc) ** 2 / (2 * rmass) + H%m(bra, bra)
      end do
      call sol%init(h)
      call sol%DiagSym(h)
      if(myrank == 0) then
        m = min(5, n)
        call v%zeros(m)
        v%v(:) = sol%eig%v(:m)
        write(*,'(a)') 'd Eigen values in CheckNNForceMOMSpace'
        write(*,'(a,i3,a,i3,a,i3)') 'J = ', j, ',  P = ', p, ',  T = ', t
        call v%prt()
        call v%fin()
      end if
      call sol%fin()
      call H%fin()
    end do
    call vmom%WeightUnNormalize()
  end subroutine check_nnforce_mom_space

  subroutine renorm_mom_space(vmom, renorm, lambda, pn_same_mass)
    type(NNForceMom), intent(inout) :: vmom
    type(str), intent(in) :: renorm
    real(8), intent(in) :: lambda
    logical, intent(in) :: pn_same_mass
    integer :: ich

    if(myrank == 0) then
      write(*,'(a)') '## NN Renormalization in mom. space ...'
    end if
    call vmom%WeightNormalize()
    if(renorm%val == 'srg') then
      do ich = 1, vmom%ms%NChan
        call MomSRGChannel(vmom%MatCh(ich), vmom%UT(ich), vmom%ms%jpsz(ich),lambda, pn_same_mass)
      end do
    elseif(renorm%val == 'vlowk') then
      do ich = 1, vmom%ms%NChan
        call MomVlowkChannel(vmom%MatCh(ich), vmom%UT(ich), vmom%ms%jpsz(ich),lambda,pn_same_mass)
      end do
    end if
    call vmom%WeightUnNormalize()
  end subroutine renorm_mom_space

  subroutine MomSRGChannel(vmom,UT,ms,lambda,pn_same_mass)
    use MyLibrary, only: hc, m_proton, m_neutron
    use Renormalization, only: SRGSolver
    type(DMat), intent(inout) :: vmom, UT
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    type(Mesh), pointer :: bra_m
    real(8), intent(in) :: lambda
    logical, intent(in) :: pn_same_mass
    integer :: n, itz, bra
    real(8) :: rmass, lc, alpha, p1
    type(DMat) :: t, h
    type(SRGSolver) :: sol
    n = ms%GetNumberStates()
    itz = ms%GetZ()
    if(pn_same_mass) then
      rmass = m_proton * m_neutron / (m_proton + m_neutron)
    else
      if(itz == -1) rmass = m_proton * 0.5d0
      if(itz ==  0) rmass = m_proton * m_neutron / (m_proton + m_neutron)
      if(itz ==  1) rmass = m_neutron * 0.5d0
    end if
    lc = lambda * hc / dsqrt(2.d0 * rmass)
    alpha = (1.d0 / lc) ** 4.d0
    call t%zeros(n,n); call h%zeros(n,n)
    do bra = 1, n
      bra_m => ms%GetP(bra)
      p1 = bra_m%GetP()
      t%m(bra,bra) = (p1 * hc) ** 2 / (2 * rmass)
    end do

    h = t + vmom
    call sol%init(h, 'Hflow')
    call sol%SRGFlow(h,t,alpha)
    vmom = sol%h - t
    UT = sol%u
    call h%fin(); call t%fin()
    call sol%fin()
  end subroutine MomSRGChannel

  subroutine MomVlowkChannel(vmom,UT,ms,lambda,pn_same_mass)
    use MyLibrary, only: hc, m_proton, m_neutron
    use Renormalization, only: LeeSuzukiSolver
    type(DMat), intent(inout) :: vmom, UT
    type(TwoBodyRelChanSpinMBasis), intent(in) :: ms
    logical, intent(in) :: pn_same_mass
    type(Mesh), pointer :: bra_m, ket_m
    real(8), intent(in) :: lambda
    integer :: n, itz, bra, ket
    integer :: i1, i2, cnt, cntb, cntk
    integer :: np, nq
    integer, allocatable :: reord(:), back(:)
    real(8) :: rmass, p1, p2, w1, w2
    type(DMat) :: t, h
    type(LeeSuzukiSolver) :: sol
    n = ms%GetNumberStates()
    itz = ms%GetZ()
    if(pn_same_mass) then
      rmass = m_proton * m_neutron / (m_proton + m_neutron)
    else
      if(itz == -1) rmass = m_proton * 0.5d0
      if(itz ==  0) rmass = m_proton * m_neutron / (m_proton + m_neutron)
      if(itz ==  1) rmass = m_neutron * 0.5d0
    end if
    allocate(reord(n), back(n))
    reord = 0; back = 0
    cnt = 0
    do i1 = 1, ms%GetNMesh()
      do bra = 1, n
        bra_m => ms%GetP(bra)
        i2 = bra_m%GetI()
        if(i1 /= i2) cycle
        cnt = cnt + 1
        reord(bra) = cnt
        back(cnt) = bra
      end do
    end do

    np = 0; nq = 0
    do bra = 1, n
      cnt = reord(bra)
      bra_m => ms%GetP(bra)
      i1 = bra_m%GetI()
      if(bra_m%GetP() <= lambda) np = np + 1
      if(bra_m%GetP() >  lambda) nq = nq + 1
    end do

    call t%zeros(n,n); call h%zeros(n,n)
    do bra = 1, n
      cnt = reord(bra)
      bra_m => ms%GetP(cnt)
      p1 = bra_m%GetP()
      t%m(bra,bra) = (p1 * hc) ** 2 / (2 * rmass)
    end do

    do bra = 1, n
      cntb = reord(bra)
      bra_m => ms%GetP(cntb)
      p1 = bra_m%GetP()
      w1 = bra_m%GetW()
      do ket = 1, n
        cntk = reord(ket)
        ket_m => ms%GetP(cntk)
        p2 = ket_m%GetP()
        w2 = ket_m%GetW()
        h%m(bra,ket) = vmom%m(cntb, cntk)
      end do
    end do
    h = t + h
    call sol%init(h)
    call sol%LeeSuzuki(0, np, nq, h)

    do bra = 1, n
      cntb = back(bra)
      bra_m => ms%GetP(cntb)
      p1 = bra_m%GetP()
      w1 = bra_m%GetW()
      do ket = 1, n
        cntk = back(ket)
        ket_m => ms%GetP(cntk)
        p2 = ket_m%GetP()
        w2 = ket_m%GetW()
        vmom%m(bra,ket) = (sol%h%m(cntb, cntk) - t%m(cntb, cntk))
      end do
    end do
    UT = sol%expS
    call sol%fin()
    call t%fin(); call h%fin()
    deallocate(reord, back)

  end subroutine MomVlowkChannel

  subroutine isospin_symmetric(vnn)
    ! Make the interaction isospin symmetric:
    ! <T=0| V | T=0> = <pn|V|pn>
    ! <T=1| V | T=1> = (<pp|V|pp> + <nn|V|nn> + <pn|V|pn>)/3
    !
    ! In the pn basis
    ! <pp | V | pp > and < nn | V | nn > => ( <pp|V|pp> + <nn|V|nn> + <pn|V|pn> )/3
    ! <pn | V | pn > =>  <pp|V|pp>/6 + <nn|V|nn>/6 + 2<pn|V|pn>/3
    !
    type(NNForceMom), intent(inout) :: vnn
    type(NNForceMom) :: vnn_copy
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: chan
    integer :: ich, j, p, s, z, ich_pp, ich_nn, ich_pn
    ms => vnn%ms
    call vnn_copy%init(ms)
    do ich = 1, ms%GetNumberChannels()
      vnn_copy%MatCh(ich) = vnn%MatCh(ich)
    end do

    do ich = 1, ms%GetNumberChannels()
      chan => ms%GetChannel(ich)
      j = chan%GetJ()
      p = chan%GetParity()
      s = chan%GetS()
      z = chan%GetZ()
      ich_pp = ms%GetIndex(j,p,s,-1)
      ich_nn = ms%GetIndex(j,p,s, 1)
      ich_pn = ms%GetIndex(j,p,s, 0)

      select case(z)
      case(-1, 1)
        vnn%MatCh(ich) = ( vnn_copy%MatCh(ich_pp) + vnn_copy%MatCh(ich_nn) + vnn_copy%MatCh(ich_pn) ) / 3.d0
      case(0)
        if(ich_pp /= 0 .and. ich_nn /= 0) then
          vnn%MatCh(ich) = vnn_copy%MatCh(ich_pp)/6.d0 + vnn_copy%MatCh(ich_nn)/6.d0 + 2.d0*vnn_copy%MatCh(ich_pn)/3.d0
        end if
      end select
    end do
    call vnn_copy%fin()
  end subroutine isospin_symmetric

  function spin_tensor_decomposition(vint, rank) result(op)
    use MyLibrary, only: triag, sjs
    type(NNForceMom), intent(in) :: vint
    type(NNForceMom) :: op
    integer, intent(in) :: rank
    type(TwoBodyRelSpaceSpinMBasis), pointer :: ms
    type(TwoBodyRelChanSpinMBasis), pointer :: relch, relch2
    type(Mesh), pointer :: mbra, mket
    integer :: jj, j, s, par, z, ich, ich2
    integer :: ibra, iket, lbra, lket, ibra2, iket2
    real(8) :: me

    ms => vint%ms
    call op%init(ms)
    do ich = 1, ms%GetNumberChannels()
      relch => ms%GetChannel(ich)
      jj = relch%GetJ()
      par = relch%GetParity()
      s = relch%GetS()
      z = relch%GetZ()
      if(triag(s,s,rank)) cycle
      do ibra = 1, relch%GetNumberStates()
        mbra => relch%GetP(ibra)
        lbra = mbra%GetL()
        do iket = 1, relch%GetNumberStates()
          mket => relch%GetP(iket)
          lket = mket%GetL()
          if(triag(lbra,lket,rank)) cycle
          me = 0.d0
          do j = max(abs(lbra-s), abs(lket-s)), min(lbra+s,lket+s,ms%GetJmax())
            ich2 = ms%GetIndex(j, par, s, z)
            if(ich2 == 0) cycle
            relch2 => ms%GetChannel(j, par, s, z)
            ibra2 = relch2%GetIndex(mbra%GetMeshI(), lbra)
            iket2 = relch2%GetIndex(mket%GetMeshI(), lket)
            me = me + dble( 2*j+1 ) * &
                & sjs(2*lbra, 2*s, 2*j, 2*s, 2*lket, 2*rank) * &
                & (-1.d0)**j * vint%MatCh(ich2)%m(ibra2,iket2)
          end do
          op%MatCh(ich)%m(ibra,iket) = me * (-1.d0)**jj * dble(2*rank+1) * &
                & sjs(2*lbra, 2*s, 2*jj, 2*s, 2*lket, 2*rank)
        end do
      end do
    end do
  end function spin_tensor_decomposition

  subroutine write_nn_mom(this, fn)
    type(NNForceMom), intent(in) :: this
    type(str), intent(in) :: fn
    type(TwoBodyRelSpaceSpinMBasis), pointer :: two
    type(TwoBodyRelChanSpinMBasis), pointer :: two_ch
    type(Mesh), pointer :: m_bra, m_ket
    integer :: ich, j, p, s, z, ntot
    integer :: i, i1, i2
    integer :: n1, n2, l1, l2
    integer :: wunit = 20

    if(sy%find(fn, sy%str('.txt'))) then
      open(wunit, file=fn%val)
      two => this%ms
      write(wunit,*) " j, prty, spin,  Tz, idx', l',  p' (fm-1),  w' (fm-1),   idx,  l,  p (fm-1),  w (fm-1),  ME (MeV fm-3)"
      do ich = 1, two%GetNumberChannels()
        two_ch => two%GetChannel(ich)
        j = two_ch%GetJ()
        p = two_ch%GetParity()
        s = two_ch%GetS()
        z = two_ch%GetZ()
        do i1 = 1, two_ch%GetNumberStates()
          m_bra => two_ch%GetP(i1)
          do i2 = 1, two_ch%GetNumberStates()
            m_ket => two_ch%GetP(i2)
            write(wunit, '(6i4,2f12.6,2i4,2f12.6,es18.6)') j, p, s, z, m_bra%GetI(), m_bra%GetL(), m_bra%GetP(), m_bra%GetW(), &
                & m_ket%GetI(), m_ket%GetL(), m_ket%GetP(), m_bra%GetW(), this%MatCh(ich)%m(i1,i2)
          end do
        end do
      end do
      close(wunit)
    end if

    if(sy%find(fn, sy%str('.bin'))) then
      two => this%ms
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) two%GetNMesh(), two%GetJmax(), two%GetNumberChannels()
      do i = 1, two%GetNMesh()
        write(wunit) two%jpsz(1)%points(i)%GetP()
      end do
      do i = 1, two%GetNMesh()
        write(wunit) two%jpsz(1)%points(i)%GetW()
      end do
      do ich = 1, two%GetNumberChannels()
        two_ch => two%GetChannel(ich)
        j = two_ch%GetJ()
        p = two_ch%GetParity()
        s = two_ch%GetS()
        z = two_ch%GetZ()
        write(wunit) j, p, s, z, this%MatCh(ich)%n_col
        write(wunit) this%MatCh(ich)%m(:,:)
      end do
      close(wunit)
    end if
  end subroutine write_nn_mom

  subroutine read_nn_mom(this, fn)
    type(NNForceMom), intent(inout) :: this
    type(str), intent(in) :: fn
    type(TwoBodyRelSpaceSpinMBasis), pointer :: two
    real(8) :: xmin_read, xmax_read
    real(8), allocatable :: tmp(:,:), pmesh(:), wmesh(:)
    integer :: NMesh_read, Jmax_read, nchan_read
    integer :: ich, j, p, s, z, ndim, idx
    integer :: runit = 21

    two => this%ms
    open(runit, file=fn%val, form='unformatted', access='stream')
    read(runit) NMesh_read, Jmax_read, nchan_read
    if(two%GetNMesh() /= NMesh_read) then
      write(*,*) "Error: NN file quadrature mesh number mismatch"
      write(*,*) __LINE__, __FILE__
      stop
    end if
    allocate(pmesh(NMesh_read), wmesh(NMesh_read))
    read(runit) pmesh
    read(runit) wmesh
    call two%SetMeshWeight(pmesh, wmesh)
    do ich = 1, nchan_read
      read(runit) j, p, s, z, ndim
      idx = two%GetIndex(j,p,s,z)
      allocate(tmp(ndim,ndim))
      read(runit) tmp(:,:)
      if(j > two%GetJmax() .or. idx==0) then
        deallocate(tmp)
        cycle
      end if
      this%MatCh(idx)%m(:,:) = tmp(:,:)
      deallocate(tmp)
    end do
    deallocate(pmesh, wmesh)
    close(runit)
  end subroutine read_nn_mom

  subroutine write_nn_ho(this, fn)
    type(NNForceHO), intent(in) :: this
    type(str), intent(in) :: fn
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: two
    type(TwoBodyRelChanSpinHOBasis), pointer :: two_ch
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: ich, j, p, s, z, ntot
    integer :: i1, i2
    integer :: n1, n2, l1, l2
    integer :: wunit = 20

    two => this%ms
    open(wunit, file=fn%val)
    write(wunit,*) " j, prty, spin,  Tz, n', l',  n,  l,  ME (MeV)"
    do ich = 1, two%GetNumberChannels()
      two_ch => two%GetChannel(ich)
      j = two_ch%GetJ()
      p = two_ch%GetParity()
      s = two_ch%GetS()
      z = two_ch%GetZ()
      do i1 = 1, two_ch%GetNumberStates()
        ho_bra => two_ch%GetP(i1)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        do i2 = 1, two_ch%GetNumberStates()
          ho_ket => two_ch%GetP(i2)
          n2 = ho_ket%GetN()
          l2 = ho_ket%GetL()
          write(wunit, '(8i4,es18.6)') j, p, s, z, n1, l1, n2, l2, this%MatCh(ich)%m(i1,i2)
        end do
      end do
    end do
    close(wunit)
  end subroutine write_nn_ho
end module NNForce

!program test
!  use Profiler, only: timer
!  use TwoBodyRelativeSpace
!  use NNForce
!  use MyLibrary
!  !type(TwoBodyRelSpaceSpinMBasis), target :: two
!  type(TwoBodyRelSpaceSpinHOBasis) :: two
!  type(NNForceHO) :: vho, u
!
!  call timer%init()
!
!  call two%init(20.d0,20,1)
!  call vho%init(two)
!  call vho%SetParams("N3LO_EMN500","srg","ho",2.d0,.True.,200,8)
!  call vho%SetNNForceHO(u, two, 500, 0.d0, 25.d0, .True., 2, 20)
!
!  call timer%fin()
!end program test
