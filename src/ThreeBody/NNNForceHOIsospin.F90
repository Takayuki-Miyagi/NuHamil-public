module NNNForceHOIsospin
  use MPIFunction, only: myrank
  use omp_lib
  use ClassSys
  use LinAlgLib
  use ThreeBodyJacOpsChanIso
  implicit none
  type, extends(ThreeBodyJacOpChanIso) :: NNNForceIsospin
  contains
    procedure :: InitNNNForce
    procedure :: SetNNNForce
    generic :: init => InitNNNForce
  end type NNNForceIsospin
  type(sys), private :: sy
contains

  subroutine InitNNNForce(this,jacobi_ch)
    class(NNNForceIsospin), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: jacobi_ch
    type(sys) :: s
    call this%init(jacobi_ch, jacobi_ch, s%str('NNNint'))
    this%is_init = .true.
  end subroutine InitNNNForce

  subroutine SetNNNForce(this,U, params)
    use NuHamilInput, only: InputParameters
    class(NNNForceIsospin), intent(inout) :: this
    type(ThreeBodyJacOpChanIso), intent(inout) :: U
    type(InputParameters), intent(in) :: params
    type(ThreeBodyJacIsoChan), pointer :: jac
    type(str) :: fv, fut
    integer :: wunit = 300, runit = 301
    type(sys) :: s

    jac => this%jacobi_ch_ket

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling SetNNNForce: NNNForceHOIsospin.F90"
      return
    end if

    if(.not. U%is_init) then
      write(*,*) "Initialize 'U' before calling SetNNNForce: NNNForceHOIsospin.F90"
      return
    end if
    fv =  U%GetFileName(jac,params%hw,s%str('NNNint'),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)
    fut = U%GetFileName(jac,params%hw,s%str('UT'),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)

    if(s%isfile(fv) .and. s%isfile(fut)) then
      open(runit,file=fv%val,status='old',form='unformatted',access='stream')
      call this%readf(jac,runit,jac%GetNmax(),params%hw,&
          & params%hw_target,params%hw_conversion)
      close(runit)

      open(runit,file=fut%val,status='old',form='unformatted',access='stream')
      call U%readf(jac,runit,jac%GetNmax(),params%hw,&
          & params%hw_target,params%hw_conversion)
      close(runit)
      return
    end if

    call set_nnn_force_ho_isospin(this,U, params)
    open(wunit,file=fv%val,status='replace',form='unformatted',access='stream')
    call this%writef(wunit)
    close(wunit)

    open(wunit,file=fut%val,status='replace',form='unformatted',access='stream')
    call U%writef(wunit)
    close(wunit)
  end subroutine SetNNNForce

  subroutine set_nnn_force_ho_isospin(this, U, params)
    use NuHamilInput, only: InputParameters
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    type(NNNForceIsospin), intent(inout) :: this
    type(ThreeBodyJacOpChanIso), intent(inout) :: U
    type(InputParameters), intent(in) :: params ! used only for kinetic term
    type(ThreeBodyJacIsoChan), pointer :: jac
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelOpIso) :: vnn_bare_rel, vnn_sub_rel, kin_rel
    type(NNForceHOIsospin) :: vnn_bare_relspin, vnn_sub_relspin, U_bare_relspin, U_sub_relspin
    type(ThreeBodyJacOpChanIso) :: T_jac, vnn_jac, vnn_sub_jac, Generator
    type(NNNForceIsospin) :: v3n_jac
    type(DMat) :: h
    type(InputParameters) :: input
    type(EigenSolSymD) :: diag
    integer :: ndim
    type(sys) :: s


    jac => this%jacobi_ch_ket
    !call relspin%init(jac%GetFrequency(), jac%GetNmax(), params%J2max_NNint)
    call relspin%init(jac%GetFrequency(), params%N2max, params%J2max_NNint)
    !call rel%init(jac%GetFrequency(), jac%GetNmax(), params%Jmax2)
    call rel%init(jac%GetFrequency(), params%N2max, params%Jmax2)
    ! -----
    call vnn_bare_relspin%init(relspin)
    call U_bare_relspin%init(relspin)
    input = params; input%renorm = s%str("bare")
    call vnn_bare_relspin%SetNNForceHOIsospin(U_bare_relspin, input)
    ! -----
    call vnn_sub_relspin%init(relspin)
    call U_sub_relspin%init(relspin)
    call vnn_sub_relspin%SetNNForceHOIsospin(U_sub_relspin, params)

    call vnn_bare_rel%init(rel,rel,s%str('NNint'))
    call vnn_sub_rel%init(rel,rel,s%str('NNint'))
    call vnn_bare_rel%SetTwoBodyScalarSpinBreaking(vnn_bare_relspin)
    call vnn_sub_rel%SetTwoBodyScalarSpinBreaking(vnn_sub_relspin)

    call vnn_bare_relspin%fin()
    call vnn_sub_relspin%fin()
    call U_bare_relspin%fin()
    call U_sub_relspin%fin()
    call relspin%fin()

    call kin_rel%init(rel,rel,s%str('kinetic'))
    call kin_rel%set()

    call T_jac%init(jac,jac,s%str('kinetic'))
    call vnn_jac%init(jac,jac,s%str('NNint'))
    call vnn_sub_jac%init(jac,jac,s%str('NNint'))
    call v3n_jac%init(jac,jac,s%str('NNNint'))

    call T_jac%set(kin_rel)
    call vnn_jac%set(vnn_bare_rel)
    call vnn_sub_jac%set(vnn_sub_rel)
    if(params%genuine_3bf) then
      call set_genuine_three_body_interaction(v3n_jac,params)
    end if
    h = T_jac%DMat + vnn_jac%DMat + v3n_jac%DMat

    call diag%init(h)
    call diag%DiagSym(h)
    ndim = min(5,size(diag%eig%v))
    write(*,"(a, i4, a, 5f12.6)") "myrank=", myrank, ", Eigen values of original 3-body H:    ", &
        & diag%eig%v(1:ndim)
    call diag%fin()
    if(params%renorm%val == 'srg') then
      select case( params%srg_generator%val )
      case( "kinetic" )
        Generator = T_jac
      case default
        Generator = set_srg_generator( jac, params%srg_generator )
      end select
      call three_body_srg_evolution(h, Generator, U%DMat, params%lambda, params%Nmax_srg_edge)
    end if
    this%DMat = h - T_jac%DMat - vnn_sub_jac%DMat
    call v3n_jac%fin()
    call vnn_jac%fin()
    call vnn_sub_jac%fin()
    call T_jac%fin()

    call kin_rel%fin()
    call vnn_bare_rel%fin()
    call vnn_sub_rel%fin()
    call rel%fin()
  end subroutine set_nnn_force_ho_isospin

  subroutine three_body_srg_evolution(h, gen, U, lambda, Nmax_srg_edge)
    use MyLibrary, only: hc, m_proton, m_neutron
    use renormalization, only: SRGSolver
    type(DMat), intent(inout) :: h, U
    type(ThreeBodyJacOpChanIso), intent(in) :: gen
    integer, intent(in) :: Nmax_srg_edge
    real(8), intent(in) :: lambda
    type(ThreeBodyJacIsoChan), pointer :: ms
    real(8) :: lc, alpha
    type(SRGSolver) :: sol
    type(EigenSolSymD) :: diag
    integer :: ndim, n_zero, ket
    integer, allocatable :: rhs_zero_term(:)

    if(size(h%m) == 0) return
    lc = lambda * hc / dsqrt(2.d0 * m_proton * m_neutron / (m_proton + m_neutron))
    alpha = (1.d0 / lc) ** 4.d0

    ! Neumann boundary for SRG flow
    ms => gen%jacobi_ch_ket
    n_zero = 0
    do ket = 1, ms%GetNumberAStates()
      if( ms%GetEAS(ket) < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
    end do
    allocate(rhs_zero_term(n_zero))
    n_zero = 0
    do ket = 1, ms%GetNumberAStates()
      if( ms%GetEAS(ket) < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
      rhs_zero_term(n_zero) = ket
    end do

    call sol%init(h, 'Hflow')
    !call sol%init(h, 'Hflow', atol=1.d-6, rtol=1.d-7)
    if(n_zero==0) call sol%SRGFlow(h, gen%DMat, alpha)
    if(n_zero> 0) call sol%SRGFlow(h, gen%DMat, alpha, rhs_zero_term)
    h = sol%h
    U = sol%U
    call sol%fin()
    call diag%init(h)
    call diag%DiagSym(h)
    ndim = min(5,size(diag%eig%v))
    write(*,"(a, i4, a, 5f12.6)") "myrank=", myrank, ", Eigen values of 3-body H (after SRG): ", &
        & diag%eig%v(1:ndim)
    call diag%fin()
  end subroutine three_body_srg_evolution

  subroutine set_genuine_three_body_interaction(this, params,verbose)
    use Profiler, only: timer
    use NuHamilInput, only: InputParameters
    class(NNNForceIsospin), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: verbose
    real(8) :: ti
    ti = omp_get_wtime()
    select case(params%NNNInt%val)
    case('ChEFT_N2LO')
      call set_nnn_interaction_chEFT_n2lo(this, params, verbose)
    case('ChEFT_N3LO')
      call set_nnn_interaction_chEFT_n3lo(this, params, verbose)
    case default
      write(*,'(a)') 'Input parameter error:'
      write(*,'(2a)') 'NNNint = ', trim(params%NNNInt%val)
    end select
    call timer%Add(sy%str('SetGenuineThreeBodyForce'),omp_get_wtime()-ti)
  end subroutine set_genuine_three_body_interaction

  subroutine set_nnn_interaction_chEFT_n2lo(this, params, verbose)
    use NuHamilInput, only: InputParameters
    class(NNNForceIsospin), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: verbose
    type(ThreeBodyJacIsoChan), pointer :: jac
    type(DMat) :: cfp, work
    real(8) :: LECs(5)
    integer :: north, nphys

    LECs = [params%c1, params%c3, params%c4, params%cd, params%ce]
    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    call work%zeros(nphys,nphys)
    cfp = jac%GetCFPMat()
    call set_nnn_int_chEFT_n2lo_isospin(work, jac, LECs, params%Regulator, &
        & params%RegulatorPower, params%Lambda_3nf_local, params%Lambda_3nf_nonlocal, &
        & params%Ngrid_3nf, params%J12max, params%save_3nf_before_lec, &
        & params%j3max_initial_3nf, params%path_to_tmp_dir, params%path_to_Hebeler_files, &
        & params%Hebeler_fn_head, params%Hebeler_fn_tails)
    this%DMat = cfp%T() * work * cfp
    call work%fin()
    call cfp%fin()
    select case(params%Regulator%val)
    case("LNL","lnl",'local-non-local',"Local-Non-Local")
      call multiply_non_local_regulator_hospace(this%DMat, this%jacobi_ch_ket, params, verbose)
    case default
    end select
  end subroutine set_nnn_interaction_chEFT_n2lo

  subroutine set_nnn_interaction_chEFT_n3lo(this, params, verbose)
    use NuHamilInput, only: InputParameters
    class(NNNForceIsospin), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: verbose
    type(ThreeBodyJacIsoChan), pointer :: jac
    type(DMat) :: cfp, work
    real(8) :: LECs(8)
    integer :: north, nphys

    LECs = [params%c1, params%c3, params%c4, params%cd, params%ce, 0.d0, 0.d0, 0.d0]
    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    call work%zeros(nphys,nphys)
    cfp = jac%GetCFPMat()
    call set_nnn_int_chEFT_n3lo_isospin(work, jac, LECs, params%Regulator, &
        & params%RegulatorPower, params%Lambda_3nf_local, params%Lambda_3nf_nonlocal, &
        & params%Ngrid_3nf, params%J12max, params%save_3nf_before_lec, &
        & params%j3max_initial_3nf, params%path_to_tmp_dir, params%path_to_Hebeler_files)
    this%DMat = cfp%T() * work * cfp
    call work%fin()
    call cfp%fin()
  end subroutine set_nnn_interaction_chEFT_n3lo

  subroutine multiply_non_local_regulator_hospace(this, jac, params, verbose)
    use MyLibrary, only: gauss_legendre
    use NuHamilInput, only: InputParameters
    type(DMat), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: verbose
    type(DMat) :: cfp, work
    real(8), allocatable :: Radial(:,:,:), p(:), pw(:)
    real(8) :: time
    integer :: north, nphys
    real(8) :: pmin = 0.d0, pmax
    integer :: NMesh = 500
    integer :: wunit=20
    type(str) :: fn
    type(sys) :: s
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    call work%zeros(nphys,nphys)
    cfp = jac%GetCFPMat()
    time = omp_get_wtime()
    pmax = 25.d0 * 6.d0 / sqrt(jac%GetFrequency()) ! estimated such that x^2 = p^2 * m * w / h ~ 25^2
    allocate(p(NMesh), pw(NMesh))
    call gauss_legendre(pmin,pmax,p,pw,Nmesh)
    call store_radial_wf()
    work = non_local_regulator_ho_mat()
    call release_radial_wf()
    deallocate(p, pw)

    work = cfp%T() * work * cfp
    if(present(verbose)) then
      fn = s%str("NNN_Local_") + params%averaged_file_for_test + s%str(".bin")
      open(wunit, file=fn%val, form="unformatted", access="stream")
      write(wunit) this%m(:,:)
      close(wunit)
      fn = s%str("NonLocal_Reg_") + params%averaged_file_for_test + s%str(".bin")
      open(wunit, file=fn%val, form="unformatted", access="stream")
      write(wunit) work%m(:,:)
      close(wunit)
    end if
    this = work%T() * this * work
#ifndef MPI
    write(*,"(a,f12.6,a)") "# Multiply NonLocal regulator: ", omp_get_wtime()-time, " sec"
#endif
    call work%fin()
    call cfp%fin()
  contains
    function non_local_regulator_ho_mat() result(mat)
      type(DMat) :: mat
      integer :: ibra, iket
      type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

      call mat%zeros(jac%GetNumberNAStates(), jac%GetNumberNAStates())
      !$omp parallel
      !$omp do private(ibra, bra, iket, ket) schedule(dynamic)
      do ibra = 1, jac%GetNumberNAStates()
        bra => jac%GetNAS(ibra)
        do iket = 1, ibra
          ket => jac%GetNAS(iket)
          if(bra%GetL12() /= ket%GetL12()) cycle
          if(bra%GetS12() /= ket%GetS12()) cycle
          if(bra%GetJ12() /= ket%GetJ12()) cycle
          if(bra%GetT12() /= ket%GetT12()) cycle
          if(bra%GetL3( ) /= ket%GetL3( )) cycle
          if(bra%GetJ3( ) /= ket%GetJ3( )) cycle
          mat%m(ibra,iket) = non_local_regulator(&
              &bra%GetN12(),ket%GetN12(),ket%GetL12(),bra%GetN3(),ket%GetN3(),ket%GetL3())
          mat%m(iket,ibra) = mat%m(ibra,iket)
          !write(*,"(6i4,es18.8)") n12, n45, n3, n6, l12, l3, FF(bra,ket)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end function non_local_regulator_ho_mat

    function non_local_regulator(n12, n45, LL, n3, n6, l) result(f)
      ! < n12 LL n3 l | F | n45 LL n6 l >, F is non-local regulator.
      use MyLibrary, only: hc
      integer, intent(in) :: n12, n45, LL, n3, n6, l
      real(8) :: f, ex
      integer :: i, k
      real(8) :: pi, pk, wi, wk

      f = 0.d0
      do i = 1, size(p)
        pi = p(i)
        wi = pw(i)
        do k = 1, size(p)
          pk = p(k)
          wk = pw(k)
          ex = - ( 0.5d0 * (pi**2 + pk**2) * hc**2 / params%lambda_3nf_nonlocal**2 )**params%RegulatorPower
          f = f + wi * wk * Radial(i,n12,LL) * Radial(k,n3,l) * &
              &   exp(ex) * Radial(i,n45,LL) * Radial(k,n6,l)
        end do
      end do
    end function non_local_regulator

    subroutine store_radial_wf()
      use MyLibrary, only: m_proton, m_neutron, hc, ho_radial_wf_norm
      integer :: N, L, i
      real(8) :: a

      allocate(Radial(size(p), 0:jac%GetNmax()/2, 0:jac%GetNmax()))
      a = 0.5d0 * (m_proton + m_neutron) * hc**2 / (m_proton * m_neutron * jac%GetFrequency())
      do N = 0, jac%GetNmax()/2
        do L = 0, jac%GetNmax()
          do i = 1, size(p)
            Radial(i,N,L) = ho_radial_wf_norm(N,L,a,p(i)) * (-1.d0)**N ! momentum space phase
          end do
        end do
      end do
    end subroutine store_radial_wf

    subroutine release_radial_wf()
      deallocate(Radial)
    end subroutine release_radial_wf
  end subroutine multiply_non_local_regulator_hospace

  subroutine set_nnn_int_chEFT_n2lo_isospin(v3n, jac, &
        & LECs, Regulator, RegulatorPower, Lambda_local, Lambda_nonlocal, Ngrid, J12max, &
        & save_3nf_before_lec, J3max_initial_3nf, path_to_dir, path_to_hdf5_files, &
        & hdf5_fn_head, hdf5_fn_tails)
    type(DMat), intent(inout) :: v3n
    type(ThreeBodyJacIsoChan), pointer, intent(in) :: jac
    real(8), intent(in) :: LECs(5), Lambda_local, Lambda_nonlocal
    integer, intent(in) :: Ngrid, J12max, RegulatorPower, J3max_initial_3nf
    type(str), intent(in) :: Regulator, path_to_dir, path_to_hdf5_files, hdf5_fn_head, hdf5_fn_tails(:)
    logical, intent(in) :: save_3nf_before_lec
    type(sys) :: s
    type(str) :: fn
    integer :: iunit = 100

    select case(Regulator%val)
    case('Local', 'local','localnonlocal',"LNL","lnl",'local-non-local',"Local-Non-Local")
      fn = get_file_name_ho_partial_wave_n2lo_isospin(jac, jac%GetNmax(), s%str("Local"), RegulatorPower, &
          & Lambda_local, Lambda_nonlocal, path_to_dir)
    case('NonLocal', 'nonlocal', 'Nonlocal', 'SemiLocal', 'semilocal', 'Semilocal')
      fn = get_file_name_ho_partial_wave_n2lo_isospin(jac, jac%GetNmax(), s%str("NonLocal"), RegulatorPower, &
          & Lambda_local, Lambda_nonlocal, path_to_dir)
    case default
      write(*,'(2a)') 'In set_nnn_int_chEFT_n2lo, regulator: ', trim(Regulator%val)
      stop
    end select
    if(s%isfile(fn)) then
      open(iunit, file = fn%val, form='unformatted', access='stream',status='old')
      call read_ho_partial_wave_n2lo(v3n, jac%GetJ(), jac%GetParity(), jac%GetT(), &
          & jac%GetNumberNAStates(), jac%GetNmax(), iunit, LECs)
      close(iunit)
      return
    end if

    select case(Regulator%val)
    case('Local', 'local','localnonlocal',"LNL","lnl",'local-non-local',"Local-Non-Local")
      call set_nnn_int_chEFT_n2lo_isospin_local(v3n, jac, &
          & LECs, Lambda_local, RegulatorPower, save_3nf_before_lec, J3max_initial_3nf,path_to_dir)
    case('NonLocal', 'nonlocal', 'Nonlocal', 'SemiLocal', 'semilocal', 'Semilocal')
      !call set_nnn_int_chEFT_n2lo_isospin_nonlocal(v3n, jac, &
      !    & LECs, Ngrid, J12max, Regulator, RegulatorPower, Lambda_nonlocal)
      call set_nnn_int_ChEFT_n2lo_nonlocal_from_hdf5(v3n, jac, &
        & LECs, Regulator, RegulatorPower, J12max, Lambda_nonlocal, Ngrid, save_3nf_before_lec, J3max_initial_3nf, &
        & path_to_dir, path_to_hdf5_files, hdf5_fn_head, hdf5_fn_tails)
    case default
      write(*,'(2a)') 'In set_nnn_int_chEFT_n2lo, regulator: ', trim(Regulator%val)
    end select
  end subroutine set_nnn_int_chEFT_n2lo_isospin

  subroutine set_nnn_int_chEFT_n3lo_isospin(v3n, jac, &
        & LECs, Regulator, RegulatorPower, Lambda_local, Lambda_nonlocal, Ngrid, J12max, &
        & save_3nf_before_lec, J3max_initial_3nf, path_to_dir, path_to_hdf5_files)
    type(DMat), intent(inout) :: v3n
    type(ThreeBodyJacIsoChan), pointer, intent(in) :: jac
    real(8), intent(in) :: LECs(:), Lambda_local, Lambda_nonlocal
    integer, intent(in) :: Ngrid, J12max, RegulatorPower, J3max_initial_3nf
    type(str), intent(in) :: Regulator, path_to_dir, path_to_hdf5_files
    logical, intent(in) :: save_3nf_before_lec
    type(sys) :: s
    type(str) :: fn
    integer :: iunit = 100

    select case(Regulator%val)
    case('Local', 'local','localnonlocal',"LNL","lnl",'local-non-local',"Local-Non-Local")
      fn = get_file_name_ho_partial_wave_n3lo_isospin(jac, jac%GetNmax(), s%str("Local"), RegulatorPower, &
          & Lambda_local, Lambda_nonlocal, path_to_dir)
    case('NonLocal', 'nonlocal', 'Nonlocal', 'SemiLocal', 'semilocal', 'Semilocal')
      fn = get_file_name_ho_partial_wave_n3lo_isospin(jac, jac%GetNmax(), s%str("NonLocal"), RegulatorPower, &
          & Lambda_local, Lambda_nonlocal, path_to_dir)
    case default
      write(*,'(2a)') 'In set_nnn_int_chEFT_n2lo, regulator: ', trim(Regulator%val)
      stop
    end select
    if(s%isfile(fn)) then
      open(iunit, file = fn%val, form='unformatted', access='stream',status='old')
      call read_ho_partial_wave_n3lo(v3n, jac%GetJ(), jac%GetParity(), jac%GetT(), &
          & jac%GetNumberNAStates(), jac%GetNmax(), iunit, LECs)
      close(iunit)
      return
    end if

    select case(Regulator%val)
    case('Local', 'local','localnonlocal',"LNL","lnl",'local-non-local',"Local-Non-Local")
      write(*,'(3a)') "# Chiral EFT N3LO 3N with ", Regulator%val, " regulator is not implemented."
      write(*,*) __LINE__, __FILE__
      stop
    case('NonLocal', 'nonlocal', 'Nonlocal', 'SemiLocal', 'semilocal', 'Semilocal')
      call set_nnn_int_ChEFT_n3lo_nonlocal_from_hdf5(v3n, jac, &
        & LECs, Regulator, RegulatorPower, J12max, Lambda_nonlocal, Ngrid, save_3nf_before_lec, J3max_initial_3nf, &
        & path_to_dir, path_to_hdf5_files)
    case default
      write(*,'(2a)') 'In set_nnn_int_chEFT_n2lo, regulator: ', trim(Regulator%val)
    end select
  end subroutine set_nnn_int_chEFT_n3lo_isospin

  function get_file_name_ho_partial_wave_n2lo_isospin(jac, Nmax, Rtype, &
        & Rpower, lambda_local, lambda_nonlocal, path_to_dir) result(r)
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    integer, intent(in) :: Nmax, Rpower
    type(str), intent(in) :: Rtype, path_to_dir
    real(8), intent(in) :: lambda_local, lambda_nonlocal
    type(str) :: r
    type(str) :: Ps, dir
    type(sys) :: s

    if(jac%GetParity() == 1) Ps = "+"
    if(jac%GetParity() ==-1) Ps = "-"
    dir = path_to_dir + s%str("/ops")
    call s%mkdir(dir%val)
    r = dir + s%str("/N2LO_3NF_HOPW_J") + s%str(jac%GetJ())
    r = r + s%str("P") + Ps + s%str("T") + s%str(jac%GetT())
    r = r + s%str("_Nmax") + s%str(Nmax)
    r = r + s%str("_") + Rtype + s%str(Rpower)
    select case(Rtype%val)
    case("LNL","lnl",'local-non-local',"Local-Non-Local")
      r = r + s%str("_") + s%str(lambda_local) + s%str("_") + s%str(lambda_nonlocal)
    case('Local', 'local')
      r = r + s%str("_") + s%str(lambda_local)
    case('NonLocal', 'nonlocal', 'Nonlocal')
      r = r + s%str("_") + s%str(lambda_nonlocal)
    end select
    r = r + s%str("_hw") + s%str(jac%GetFrequency())
    r = r + s%str(".bin")
  end function get_file_name_ho_partial_wave_n2lo_isospin

  function get_file_name_ho_partial_wave_n3lo_isospin(jac, Nmax, Rtype, &
        & Rpower, lambda_local, lambda_nonlocal, path_to_dir) result(r)
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    integer, intent(in) :: Nmax, Rpower
    type(str), intent(in) :: Rtype, path_to_dir
    real(8), intent(in) :: lambda_local, lambda_nonlocal
    type(str) :: r
    type(str) :: Ps, dir
    type(sys) :: s

    if(jac%GetParity() == 1) Ps = "+"
    if(jac%GetParity() ==-1) Ps = "-"
    dir = path_to_dir + s%str("/ops")
    call s%mkdir(dir%val)
    r = dir + s%str("/N3LO_3NF_HOPW_J") + s%str(jac%GetJ())
    r = r + s%str("P") + Ps + s%str("T") + s%str(jac%GetT())
    r = r + s%str("_Nmax") + s%str(Nmax)
    r = r + s%str("_") + Rtype + s%str(Rpower)
    select case(Rtype%val)
    case("LNL","lnl",'local-non-local',"Local-Non-Local")
      r = r + s%str("_") + s%str(lambda_local) + s%str("_") + s%str(lambda_nonlocal)
    case('Local', 'local')
      r = r + s%str("_") + s%str(lambda_local)
    case('NonLocal', 'nonlocal', 'Nonlocal')
      r = r + s%str("_") + s%str(lambda_nonlocal)
    end select
    r = r + s%str("_hw") + s%str(jac%GetFrequency())
    r = r + s%str(".bin")
  end function get_file_name_ho_partial_wave_n3lo_isospin

  subroutine read_ho_partial_wave_n2lo(v, j, p, t, n, Nmax, iunit, LECs)
    type(DMat), intent(inout) :: v
    integer, intent(in) :: j, p, t, n, Nmax, iunit
    real(8), intent(in) :: LECs(:)
    integer :: j_in, p_in, t_in, Nmax_in
    real(8), allocatable :: vv(:)
    integer :: bra, ket
    integer(8) :: ndim

    v%m(:,:) = 0.d0
    read(iunit) j_in, p_in, t_in, Nmax_in, ndim
    if(j_in /= j .or. p_in /= p .or. t_in /= t .or. &
        & Nmax_in /= Nmax .or. &
        & ndim /= int(n+1,kind(ndim)) * int(n,kind(ndim)) / int(2,kind(ndim) )) then
      write(*,'(a)') 'File reading error in read_ho_partial_wave'
      stop
    end if
    allocate(vv(ndim))
    read(iunit) vv
    call add_v( LECs(1) )
    read(iunit) vv
    call add_v( LECs(2) )
    read(iunit) vv
    call add_v( LECs(3) )
    read(iunit) vv
    call add_v( LECs(4) )
    read(iunit) vv
    call add_v( LECs(5) )
    deallocate(vv)
    do bra = 1, n
      do ket = 1, bra
        v%m(ket,bra) = v%m(bra,ket)
      end do
    end do
  contains
    subroutine add_v(lec)
      real(8), intent(in) :: lec
      integer :: bra, ket
      integer(8) :: nelm
      do bra = 1, n
        do ket = 1, bra
          nelm = int(bra,kind(nelm)) * int(bra-1,kind(nelm)) / int(2,kind(nelm)) +&
              & int(ket,kind(nelm))
          v%m(bra,ket) = v%m(bra,ket) + 3.d0 * lec * vv(nelm)
        end do
      end do
    end subroutine add_v
  end subroutine read_ho_partial_wave_n2lo

  subroutine read_ho_partial_wave_n3lo(v, j, p, t, n, Nmax, iunit, LECs)
    type(DMat), intent(inout) :: v
    integer, intent(in) :: j, p, t, n, Nmax, iunit
    real(8), intent(in) :: LECs(5)
    integer :: j_in, p_in, t_in, Nmax_in
    real(8), allocatable :: vv(:)
    integer :: bra, ket
    integer(8) :: ndim

    v%m(:,:) = 0.d0
    read(iunit) j_in, p_in, t_in, Nmax_in, ndim
    if(j_in /= j .or. p_in /= p .or. t_in /= t .or. &
        & Nmax_in /= Nmax .or. &
        & ndim /= int(n+1,kind(ndim)) * int(n,kind(ndim)) / int(2,kind(ndim) )) then
      write(*,'(a)') 'File reading error in read_ho_partial_wave'
      stop
    end if
    allocate(vv(ndim))
    read(iunit) vv
    call add_v( LECs(1) )
    read(iunit) vv
    call add_v( LECs(2) )
    read(iunit) vv
    call add_v( LECs(3) )
    read(iunit) vv
    call add_v( LECs(4) )
    read(iunit) vv
    call add_v( LECs(5) )
    read(iunit) vv
    call add_v( 1.d0 )
    read(iunit) vv
    call add_v( 1.d0 )
    read(iunit) vv
    call add_v( 1.d0 )
    read(iunit) vv
    call add_v( 1.d0 )
    read(iunit) vv
    call add_v( 0.d0 )
    read(iunit) vv
    call add_v( 0.d0 )
    read(iunit) vv
    call add_v( 1.d0 )
    deallocate(vv)
    do bra = 1, n
      do ket = 1, bra
        v%m(ket,bra) = v%m(bra,ket)
      end do
    end do
  contains
    subroutine add_v(lec)
      real(8), intent(in) :: lec
      integer :: bra, ket
      integer(8) :: nelm
      do bra = 1, n
        do ket = 1, bra
          nelm = int(bra,kind(nelm)) * int(bra-1,kind(nelm)) / int(2,kind(nelm)) +&
              & int(ket,kind(nelm))
          v%m(bra,ket) = v%m(bra,ket) + 3.d0 * lec * vv(nelm)
        end do
      end do
    end subroutine add_v
  end subroutine read_ho_partial_wave_n3lo

  ! limitation: power of local regulator and non-local regulator have to be same
  subroutine set_nnn_int_chEFT_n2lo_isospin_local(this, jac, LECs, lam_local, RegulatorPower, &
      & save_3nf_before_lec, J3max_initial_3nf, path_to_dir)
    use NNNForceLocal, only: NNNIntLocal
    use MyLibrary, only: gauss_legendre
    type(DMat), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    real(8), intent(in) :: LECs(5), lam_local
    integer, intent(in) :: RegulatorPower, J3max_initial_3nf
    logical, intent(in) :: save_3nf_before_lec
    type(str), intent(in) :: path_to_dir
    integer(8) :: n
    type(NNNIntLocal) :: nnn_force
    real(8), allocatable :: v(:)
    type(str) :: fn
    type(sys) :: s
    !
    if (jac%GetNumberNAStates() < 1) return
    if( jac%GetJ() > J3max_initial_3nf ) return
    fn = get_file_name_ho_partial_wave_n2lo_isospin(jac, jac%GetNmax(), &
        & s%str("Local"), RegulatorPower, lam_local, 0.d0, path_to_dir)

    n = int(jac%GetNumberNAStates(),kind(n)) * int(jac%GetNumberNAStates()+1,kind(n)) / int(2,kind(n))
    allocate(v(n))
    call this%zeros(jac%GetNumberNAStates(), jac%GetNumberNAStates())

    call nnn_force%init(jac%GetJ(), jac%GetParity(), jac%GetT(), jac%GetNmax(), &
        & jac%GetFrequency(), xmin_in=0.d0, xmax_in=15.d0, &
        & lambda_in = lam_local, regulator_power_in = RegulatorPower)

    if(save_3nf_before_lec) open(50, file = fn%val, form = 'unformatted', access = 'stream')
    if(save_3nf_before_lec) write(50) jac%GetJ(), jac%GetParity(), jac%GetT(), jac%GetNmax(), n

    call nnn_force%set("c1")
    call set_inside(this, v, nnn_force, LECs(1))
    if(save_3nf_before_lec) write(50) v
    call nnn_force%release()

    call nnn_force%set("c3")
    call set_inside(this, v, nnn_force, LECs(2))
    if(save_3nf_before_lec) write(50) v
    call nnn_force%release()

    call nnn_force%set("c4")
    call set_inside(this, v, nnn_force, LECs(3))
    if(save_3nf_before_lec) write(50) v
    call nnn_force%release()

    call nnn_force%set("cD")
    call set_inside(this, v, nnn_force, LECs(4))
    if(save_3nf_before_lec) write(50) v
    call nnn_force%release()

    call nnn_force%set("cE")
    call set_inside(this, v, nnn_force, LECs(5))
    if(save_3nf_before_lec) write(50) v
    call nnn_force%release()
    call nnn_force%fin()

    if(save_3nf_before_lec) close(50)
    deallocate(v)
  contains

    subroutine set_inside(that, vtmp, tmp, lec)
      type(DMat), intent(inout) :: that
      real(8), intent(inout) :: vtmp(:)
      type(NNNIntLocal), intent(in) :: tmp
      real(8), intent(in) :: lec
      integer :: ibra, iket
      type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
      integer(8) :: nelm
      real(8) :: v3

      !$omp parallel
      !$omp do private(ibra, bra, iket, ket, nelm, v3) schedule (dynamic)
      do ibra = 1, jac%GetNumberNAStates()
        bra => jac%GetNAS(ibra)
        do iket = 1, ibra
          ket => jac%GetNAS(iket)
          nelm = int(ibra,kind(nelm)) * int(ibra-1,kind(nelm)) / int(2,kind(nelm)) + &
              & int(iket,kind(nelm))
          v3 = tmp%get(bra%GetN12(),bra%GetL12(),bra%GetS12(),bra%GetJ12(),bra%GetT12(),bra%GetN3(),bra%GetL3(),bra%GetJ3(),&
              &        ket%GetN12(),ket%GetL12(),ket%GetS12(),ket%GetJ12(),ket%GetT12(),ket%GetN3(),ket%GetL3(),ket%GetJ3())
          vtmp(nelm) = v3
          that%m(ibra,iket) = this%m(ibra,iket) + 3.d0 * v3 * lec
          that%m(iket,ibra) = this%m(ibra,iket)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine set_inside
  end subroutine set_nnn_int_chEFT_n2lo_isospin_local

  subroutine set_nnn_int_ChEFT_n2lo_nonlocal_from_hdf5(v, jac, &
      & LECs, Regulator, RegulatorPower, J12max, Lambda, NMesh, save_3nf_before_lec, J3max_file, &
      & path_to_dir, path_to_hdf5_files, fn_head, fn_tails)
    type(DMat), intent(inout) :: v
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    real(8), intent(in) :: LECs(5), Lambda
    integer, intent(in) :: RegulatorPower, J12max
    type(str), intent(in) :: Regulator, path_to_dir, path_to_hdf5_files, fn_head, fn_tails(:)
    logical, intent(in) :: save_3nf_before_lec
    integer, intent(in) :: NMesh, J3max_file
    type(str) :: fnho
    integer :: n, iop
    real(8), allocatable :: vtmp(:)
    integer :: wunit=50
    integer(8) :: nelm
    n = jac%GetNumberNAStates()
    call v%zeros(n, n)
    if(jac%GetJ() > J3max_file) then
      v%m(:,:) = 0.d0
      return
    end if

    if (n < 1) return
    nelm = int(n,kind(nelm)) * int(n+1,kind(nelm)) / int(2,kind(nelm))
    allocate(vtmp(nelm))

    fnho = get_file_name_ho_partial_wave_n2lo_isospin(jac, jac%GetNmax(), Regulator, RegulatorPower, 0.d0, Lambda, path_to_dir)
    if(save_3nf_before_lec) open(wunit, file=fnho%val, form="unformatted", access="stream")
    if(save_3nf_before_lec) write(wunit) jac%GetJ(), jac%GetParity(), jac%GetT(), jac%GetNmax(), nelm
    do iop = 1, size(fn_tails)
      call inside(vtmp, fn_tails(iop))
      !if(Ops(iop)=="TwoPionC1") call inside(vtmp, "c1", LECs(1))
      !if(Ops(iop)=="TwoPionC3") call inside(vtmp, "c3", LECs(2))
      !if(Ops(iop)=="TwoPionC4") call inside(vtmp, "c4", LECs(3))
      !if(Ops(iop)=="OnePion"  ) call inside(vtmp, "cD", LECs(4))
      !if(Ops(iop)=="Contact"  ) call inside(vtmp, "cE", LECs(5))
      if(save_3nf_before_lec) write(wunit) vtmp
    end do
    if(save_3nf_before_lec) close(wunit)
    deallocate(vtmp)

  contains
    subroutine inside(vv, fn_tail)
      use NNNFFromFile, only: ChEFT3NIntHO
      real(8), intent(inout) :: vv(:)
      type(str) :: fn_tail
      real(8) :: lec_val
      type(str) :: fn
      type(ChEFT3NIntHO) :: vho
      integer :: ibra, iket
      real(8) :: times(3)
      integer(8) :: nelm
      type(sys) :: s
      type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

      lec_val = 1.d0
      if(s%find(fn_tail, s%str('c1'))) lec_val = LECs(1)
      if(s%find(fn_tail, s%str('c3'))) lec_val = LECs(2)
      if(s%find(fn_tail, s%str('c4'))) lec_val = LECs(3)
      if(s%find(fn_tail, s%str('cD'))) lec_val = LECs(4)
      if(s%find(fn_tail, s%str('cE'))) lec_val = LECs(5)
      if(s%find(fn_head, s%str('1PV1P'))) lec_val = lec_val / 9.d0

      fn = get_hdf5_file_name(jac, fn_tail, path_to_hdf5_files)
      vv(:) = 0.d0
      if(.not. s%isfile(fn)) then
        write(*,*) "Not found file: ", trim(fn%val)
        return
      end if
      times = vho%init(jac%GetJ(), jac%GetParity(), jac%GetT(), &
          & jac%GetNmax(), jac%GetFrequency(), fn%val, &
          &  RegulatorPower=RegulatorPower, Lambda=Lambda, &
          &  J12max=J12max, MeshNum=NMesh)
      ! test mode (w/o interpolation)
      !tiems =  vho%init(jac%j, jac%p, jac%t, jac%Nmax, jac%hw, fn, &
      !  &  RegulatorPower=RegulatorPower, Lambda=Lambda, &
      !  &  J12max=J12max, MeshNum=NMesh, intrp=.false.)
      if(myrank == 0) then
        write(*,"(3a,f12.4,a,f12.4,a,f12.4,a)") "# calculating ", trim(fn_tail%val), " term, reading file:", &
            & times(1), " sec, interpolation:", times(2), " sec, mom->HO:", times(3), " sec"
      end if
      call vho%SetLEC(1.d0)

      !$omp parallel
      !$omp do private(ibra, bra, iket, ket, nelm) schedule(dynamic)
      do ibra = 1, jac%GetNumberNAStates()
        bra => jac%GetNAS(ibra)
        do iket = 1, ibra
          ket => jac%GetNAS(iket)
          nelm = int(ibra,kind(nelm)) * int(ibra-1,kind(nelm)) / int(2,kind(nelm)) +&
              &  int(iket,kind(nelm))
          vv(nelm) = &
              & vho%get(bra%GetN12(),bra%GetL12(),bra%GetS12(),bra%GetJ12(),bra%GetT12(),bra%GetN3(),bra%GetL3(),bra%GetJ3(),&
              &         ket%GetN12(),ket%GetL12(),ket%GetS12(),ket%GetJ12(),ket%GetT12(),ket%GetN3(),ket%GetL3(),ket%GetJ3())
          v%m(ibra,iket) = v%m(ibra,iket) + 3.d0 * vv(nelm) * lec_val
          v%m(iket,ibra) = v%m(ibra,iket)
        end do
      end do
      !$omp end do
      !$omp end parallel
      call vho%fin()
    end subroutine inside

    function get_hdf5_file_name(jac, kw, path) result(fn)
      type(ThreeBodyJacIsoChan), intent(in) :: jac
      type(str), intent(in) :: kw, path
      type(str) :: fn
      type(sys) :: s
      fn = path + s%str("/T3_") + s%str(jac%GetT()) + s%str("/J3_") + s%str(jac%GetJ()) + &
        & s%str("/PAR_") + s%str(jac%GetParity()) + &
        & s%str("/") + fn_head + s%str("_J3_") + s%str(jac%GetJ()) + s%str("_PAR_") + s%str(jac%GetParity()) + &
        & s%str("_T3_") + s%str(jac%GetT()) + s%str("_") + kw + s%str(".h5")
    end function get_hdf5_file_name

  end subroutine set_nnn_int_ChEFT_n2lo_nonlocal_from_hdf5

  subroutine set_nnn_int_ChEFT_n3lo_nonlocal_from_hdf5(v, jac, &
      & LECs, Regulator, RegulatorPower, J12max, Lambda, NMesh, save_3nf_before_lec, J3max_file, &
      & path_to_dir, path_to_hdf5_files)
    type(DMat), intent(inout) :: v
    type(ThreeBodyJacIsoChan), intent(in) :: jac
    real(8), intent(in) :: LECs(5), Lambda
    integer, intent(in) :: RegulatorPower, J12max
    type(str), intent(in) :: Regulator, path_to_dir, path_to_hdf5_files
    logical, intent(in) :: save_3nf_before_lec
    integer, intent(in) :: NMesh, J3max_file
    type(str) :: fnho
    character(20) :: Ops(12)
    integer :: n, iop
    real(8), allocatable :: vtmp(:)
    integer :: wunit=50
    integer(8) :: nelm
    Ops(1) = 'c1'
    Ops(2) = 'c3'
    Ops(3) = 'c4'
    Ops(4) = 'cD'
    Ops(5) = 'cE'
    Ops(6) = '2pi-1pi'
    Ops(7) = '2pi_const'
    Ops(8) = '2pi'
    Ops(9) = 'relcorr_2pi'
    Ops(10) = 'relcorr_CS'
    Ops(11) = 'relcorr_CT'
    Ops(12) = 'rings'
    n = jac%GetNumberNAStates()
    call v%zeros(n, n)
    if(jac%GetJ() > J3max_file) then
      v%m(:,:) = 0.d0
      return
    end if

    if (n < 1) return
    nelm = int(n,kind(nelm)) * int(n+1,kind(nelm)) / int(2,kind(nelm))
    allocate(vtmp(nelm))

    fnho = get_file_name_ho_partial_wave_n3lo_isospin(jac, jac%GetNmax(), Regulator, RegulatorPower, 0.d0, Lambda, path_to_dir)
    if(save_3nf_before_lec) open(wunit, file=fnho%val, form="unformatted", access="stream")
    if(save_3nf_before_lec) write(wunit) jac%GetJ(), jac%GetParity(), jac%GetT(), jac%GetNmax(), nelm
    do iop = 1, 12
      select case(Ops(iop))
      case("c1")
        call inside(vtmp, "N2LO_c1", LECs(1))
      case("c3")
        call inside(vtmp, "N2LO_c3", LECs(2))
      case("c4")
        call inside(vtmp, "N2LO_c4", LECs(3))
      case("cD")
        call inside(vtmp, "N2LO_cD", LECs(4))
      case("cE")
        call inside(vtmp, "N2LO_cE", LECs(5))
      case("2pi-1pi")
        call inside(vtmp, "N3LO_2pi-1pi", 1.d0)
      case("2pi_const")
        call inside(vtmp, "N3LO_2pi_const", 1.d0)
      case("2pi")
        call inside(vtmp, "N3LO_2pi", 1.d0)
      case("relcorr_2pi")
        call inside(vtmp, "N3LO_relcorr_2pi", 1.d0)
      case("relcorr_CS")
        call inside(vtmp, "N3LO_relcorr_CS", 0.d0)
      case("relcorr_CT")
        call inside(vtmp, "N3LO_relcorr_CT", 0.d0)
      case("rings")
        call inside(vtmp, "N3LO_rings", 1.d0)
      end select
      if(save_3nf_before_lec) write(wunit) vtmp
    end do
    if(save_3nf_before_lec) close(wunit)
    deallocate(vtmp)

  contains
    subroutine inside(vv, kwd, lec_val)
      use NNNFFromFile, only: ChEFT3NIntHO
      real(8), intent(inout) :: vv(:)
      character(*), intent(in) :: kwd
      real(8), intent(in) :: lec_val
      type(str) :: fn
      type(ChEFT3NIntHO) :: vho
      integer :: ibra, iket
      real(8) :: times(3)
      integer(8) :: nelm
      type(sys) :: s
      type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

      fn = get_hdf5_file_name(jac, s%str(kwd), path_to_hdf5_files)
      vv(:) = 0.d0
      if(.not. s%isfile(fn)) then
        write(*,*) "Not found file: ", trim(fn%val)
        return
      end if
      times = vho%init(jac%GetJ(), jac%GetParity(), jac%GetT(), &
          & jac%GetNmax(), jac%GetFrequency(), fn%val, &
          &  RegulatorPower=RegulatorPower, Lambda=Lambda, &
          &  J12max=J12max, MeshNum=NMesh)
      if(myrank == 0) then
        write(*,"(3a,f12.4,a,f12.4,a,f12.4,a)") "# calculating ", trim(kwd), " term, reading file:", &
            & times(1), " sec, interpolation:", times(2), " sec, mom->HO:", times(3), " sec"
      end if
      call vho%SetLEC(1.d0)

      !$omp parallel
      !$omp do private(ibra, bra, iket, ket, nelm) schedule(dynamic)
      do ibra = 1, jac%GetNumberNAStates()
        bra => jac%GetNAS(ibra)
        do iket = 1, ibra
          ket => jac%GetNAS(iket)
          nelm = int(ibra,kind(nelm)) * int(ibra-1,kind(nelm)) / int(2,kind(nelm)) +&
              &  int(iket,kind(nelm))
          vv(nelm) = &
              & vho%get(bra%GetN12(),bra%GetL12(),bra%GetS12(),bra%GetJ12(),bra%GetT12(),bra%GetN3(),bra%GetL3(),bra%GetJ3(),&
              &         ket%GetN12(),ket%GetL12(),ket%GetS12(),ket%GetJ12(),ket%GetT12(),ket%GetN3(),ket%GetL3(),ket%GetJ3())
          v%m(ibra,iket) = v%m(ibra,iket) + 3.d0 * vv(nelm) * lec_val
          v%m(iket,ibra) = v%m(ibra,iket)
        end do
      end do
      !$omp end do
      !$omp end parallel
      call vho%fin()
    end subroutine inside

    function get_hdf5_file_name(jac, kwd, path) result(fn)
      type(ThreeBodyJacIsoChan), intent(in) :: jac
      type(str), intent(in) :: kwd, path
      type(str) :: fn
      type(sys) :: s
      fn = path + s%str("/T3_") + s%str(jac%GetT()) + s%str("/J3_") + s%str(jac%GetJ()) + &
        & s%str("/PAR_") + s%str(jac%GetParity()) + &
        & s%str("/3NF_V_J3_") + s%str(jac%GetJ()) + s%str("_PAR_") + s%str(jac%GetParity()) + &
        & s%str("_T3_") + s%str(jac%GetT()) + "_" + kwd + s%str(".h5")
    end function get_hdf5_file_name

  end subroutine set_nnn_int_ChEFT_n3lo_nonlocal_from_hdf5

  function set_srg_generator( modelspace, generator_type ) result(generator)
    use TwoBodyRelativeSpace
    use TwoBodyRelOpsIso
    type(ThreeBodyJacIsoChan), intent(in) :: modelspace
    type(str), intent(in) :: generator_type
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelOpIso) :: gen_rel
    type(ThreeBodyJacOpChanIso) :: generator
    call rel%init(modelspace%GetFrequency(), modelspace%GetNmax(), 0)
    call gen_rel%init(rel,rel,generator_type)
    call gen_rel%set()
    call generator%init(modelspace,modelspace,generator_type)
    call generator%set( gen_rel )
    call gen_rel%fin()
    call rel%fin()
    select case( generator_type%val )
    case("kinetic","Kinetic")
    case  default
      generator%DMat = generator%DMat * 2.d0 / 3.d0
    end select
  end function set_srg_generator

end module NNNForceHOIsospin
