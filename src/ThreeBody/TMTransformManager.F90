module TMTransformManager
  use omp_lib
  use MPIFunction
  use ClassSys
  use MyLibrary
  use half_precision_floating_points, only: my_real16
  use Profiler, only: timer
  use NuHamilInput, only: InputParameters
  use TwoBodyRelativeSpace
  use TwoBodyRelOpsIso
  use NNForceIsospin
  use ThreeBodyJacOpsChanIso
  use ThreeBodyJacOpsIso
  use SingleParticleState
  use ThreeBodyJacobiSpace
  use ThreeBodyLabSpace
  use TransformationCoefficient
  use ThreeBodyLabSpaceIsoMon
  use ThreeBodyLabSpaceIsoNO2B
  use ThreeBodyLabOpsIso
  use ThreeBodyVMonIso
  use ThreeBodyNO2BIsospin
  use ThreeBodyTensorNO2BIsospin
  implicit none
  public :: manage_trans_to_lab
  public :: set_induced_nnn
  public :: set_nnn_component
  private
  type(sys) :: sy
contains
  subroutine manage_trans_to_lab(params)
    type(InputParameters), intent(in) :: params
    if(params%count_memory .and. .not. params%pn_form) then
      if( params%only_hf_monopole ) call count_three_body_monopole_elements_isospin(params)
      if( params%only_no2b_elements ) call count_three_body_no2b_elements_isospin(params)
      if(.not. params%only_hf_monopole .and. .not. params%only_no2b_elements) then
        call count_three_body_matrix_elements_isospin(params)
      end if
      return
    end if

    if(params%lab_3bme_precision%val=="half") call trans_to_lab_basis_isospin_half(params)
    if(params%lab_3bme_precision%val=="single") call trans_to_lab_basis_isospin_single(params)
    if(params%lab_3bme_precision%val=="double") call trans_to_lab_basis_isospin_double(params)

  end subroutine manage_trans_to_lab

  subroutine count_three_body_matrix_elements_isospin(params)
    type(OrbitsIsospin) :: sps
    type(InputParameters), intent(in) :: params
    type(ThreeBodyLabIsoSpace) :: lab
    type(ThreeBodyLabIsoChan), pointer :: ChLab
    integer :: ch
    real(8) :: mem, fact

    write(*,"(a)") "# Memory estimation mode: "
    write(*,"(a)") "# Three-body matrix element w/ isospin formalism (Scalar)"
    write(*,"(a,i3,a,i3,a,i3)") "# emax=", params%emax, ", e2max=", params%e2max, &
        & ", e3max=", params%e3max
    if(params%lab_3bme_precision%val=="half") then
      write(*,"(a)") "# half precision"
      fact = 1.d0
    end if
    if(params%lab_3bme_precision%val=="single") then
      write(*,"(a)") "# single precision"
      fact = 2.d0
    end if
    if(params%lab_3bme_precision%val=="double") then
      write(*,"(a)") "# double precision"
      fact = 4.d0
    end if
    call sps%init(params%emax, params%lmax)
    call lab%init(sps, params%e2max, params%e3max, params%hw_target)

    mem = 0.d0
    do ch = 1, lab%GetNumberChannels()
      ChLab => lab%GetChannel(ch)
      mem = mem + dble(ChLab%GetNumberStates()) * dble(ChLab%GetNumberStates()+1) /1024.d0**3
    end do
    write(*,"(f18.6,a)") mem*fact, " GB"
  end subroutine count_three_body_matrix_elements_isospin

  subroutine count_three_body_monopole_elements_isospin(params)
    type(OrbitsIsospin) :: sps
    type(InputParameters), intent(in) :: params
    type(MonThreeBodyLabIsoSpace) :: lab
    integer :: ch
    real(8) :: mem, fact

    write(*,"(a)") "# Memory estimation mode: "
    write(*,"(a)") "# Three-body matrix element w/ isospin formalism (Scalar, monopole)"
    write(*,"(a,i3,a,i3,a,i3)") "# emax=", params%emax, ", e2max=", params%e2max, &
        & ", e3max=", params%e3max
    if(params%lab_3bme_precision%val=="half") then
      write(*,"(a)") "# half precision"
      fact = 2.d0
    end if
    if(params%lab_3bme_precision%val=="single") then
      write(*,"(a)") "# single precision"
      fact = 4.d0
    end if
    if(params%lab_3bme_precision%val=="double") then
      write(*,"(a)") "# double precision"
      fact = 8.d0
    end if

    call sps%init(params%emax, params%lmax)
    call lab%init(sps, params%e2max, params%e3max)

    mem = 0.d0
    do ch = 1, lab%GetNumberChannels()
      mem = mem + dble(lab%chan(ch)%GetNumberStates()) * dble(lab%chan(ch)%GetNumberStates()) /1024.d0**3
    end do

    write(*,"(a,i14)") "Number of channels: ", lab%GetNumberChannels()
    write(*,"(f18.6,a)") mem*fact, " GB"
    call lab%fin()
    call sps%fin()
  end subroutine count_three_body_monopole_elements_isospin

  subroutine count_three_body_no2b_elements_isospin(params)
    type(OrbitsIsospin) :: sps
    type(InputParameters), intent(in) :: params
    type(ThreeBodyLabIsoSpaceNO2B) :: lab_no2b
    type(NO2BThreeBodyIsoSpace) :: lab
    integer :: ch
    real(8) :: mem, fact

    write(*,"(a)") "# Memory estimation mode: "
    write(*,"(a)") "# Three-body matrix element w/ isospin formalism (Scalar, NO2B)"
    write(*,"(a,i3,a,i3,a,i3)") "# emax=", params%emax, ", e2max=", params%e2max, &
        & ", e3max=", params%e3max
    if(params%lab_3bme_precision%val=="half") then
      write(*,"(a)") "# half precision"
      fact = 1.d0
    end if
    if(params%lab_3bme_precision%val=="single") then
      write(*,"(a)") "# single precision"
      fact = 2.d0
    end if
    if(params%lab_3bme_precision%val=="double") then
      write(*,"(a)") "# double precision"
      fact = 4.d0
    end if

    call sps%init(params%emax, params%lmax)
    call lab%init(sps, params%e2max, params%e3max)
    call lab_no2b%init(sps, params%e2max, params%e3max)

    mem = 0.d0
    do ch = 1, lab%GetNumberChannels()
      mem = mem + dble(lab%chan(ch)%GetNumberStates()) * dble(lab%chan(ch)%GetNumberStates()+1)/1024.d0**3
    end do

    write(*,"(a,i8)") "Number of channels, middle: ", lab_no2b%GetNumberChannels()
    write(*,"(a,i8)") "Number of channels, final:  ", lab%GetNumberChannels()
    write(*,"(f18.6,a)") mem*fact, " GB"
    call lab%fin()
    call sps%fin()
  end subroutine count_three_body_no2b_elements_isospin

  subroutine set_induced_nnn(this, params)
    type(ThreeBodyJacOpIso), intent(inout), target :: this
    type(InputParameters), intent(in) :: params
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(str) :: opname
    integer, allocatable :: slranks(:)
    type(str), allocatable :: strings(:)
    type(InputParameters) :: input
    type(sys) :: s
    real(8) :: time
#ifdef MPI
    integer :: ch, nelms
#endif
    if(params%renorm%val == "bare") return
    opname = this%GetOpName()
    select case(opname%val)
    case("NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
      if(params%NNint%val /= "N3LO_EM500") then
        write(*,*) "Error: inconsitent NN interaction and operator name"
        stop
      end if
      call s%split(opname, s%str("NNNinduced_"), strings)
      input = params
      input%NNint = strings(2); input%renorm="bare"
      call relspin%init(params%hw, params%N2max, params%J2max_NNint)
      call vnn_relspin%init(relspin)
      call U_relspin%init(relspin)
      call vnn_relspin%setNNForceHOIsospin(U_relspin,input)
      call vnn_relspin%fin()
      call U_relspin%fin()
      call relspin%fin()
    end select
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    call parent_child_procedure(set_induced_nnn_ch, this%ms%GetNumberChannels(), slranks, time)
    if(nprocs>1) call timer%add(sy%str('MPI parent-child, nnn_force_induced'), time)
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(slranks(1), this%ms%GetNumberChannels(), mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    do ch = 1, this%ms%GetNumberChannels()
      nelms = this%ms%jpt(ch)%GetNumberAStates()**2
      if(nelms < 1) cycle
      call mpi_bcast(this%MatCh(ch,ch)%m(1,1), nelms, mpi_real8, slranks(ch), mpi_comm_world, ierr)
    end do
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(slranks)
  contains
    subroutine set_induced_nnn_ch(loop)
      integer, intent(in) :: loop
      type(ThreeBodyJacOpChanIso), pointer :: v
      type(TwoBodyRelOpIso) :: vnn, op2, op2_sub, t_rel
      type(TwoBodyRelSpaceIsoHOBasis) :: rel
      type(ThreeBodyJacIsoChan) :: chjac
      type(ThreeBodyJacOpChanIso) :: tmp, ut, t_orig, t_eff, vnn_jac
      integer :: iut = 101, Nmax_file
      logical :: ex
      type(str) :: fut, fn_jacobi

      v => this%MatCh(loop,loop)
      if(v%jacobi_ch_bra%GetNumberAStates() < 1) return
      if(v%jacobi_ch_ket%GetNumberAStates() < 1) return
      Nmax_file = GetRampNmax(v%jacobi_ch_ket%GetJ(), params%ramp)
      fn_jacobi = chjac%GetFileName(v%jacobi_ch_ket%GetJ(), &
          &  v%jacobi_ch_ket%GetParity(), v%jacobi_ch_ket%GetT(), Nmax_file, params%path_to_tmp_dir)
      ex = s%isfile(fn_jacobi, s%str('In set_three_body_jac_operator_isospin'))
      open(iut, file=fn_jacobi%val, form='unformatted', access='stream')
      call chjac%readf(params%hw,iut)
      close(iut)

      call relspin%init(params%hw, chjac%GetNmax(), params%J2max_NNint)
      call vnn_relspin%init(relspin)
      call U_relspin%init(relspin)
      input = params
      input%N_ls = params%N3max; input%renorm="bare"
      call vnn_relspin%setNNForceHOIsospin(U_relspin,input)

      call rel%init(params%hw, chjac%GetNmax(), params%Jmax2)
      call vnn%init(rel,rel,s%str("NNint"))
      call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)

      call vnn_relspin%fin()
      call U_relspin%fin()
      call relspin%fin()

      select case(opname%val)
      case("NNNfrom2N","NNNinduced")
        call op2%init(rel,rel,opname)
        op2 = vnn
      case("NNNfrom2N_central")
        call op2%init(rel,rel,opname)
        op2 = vnn%SpinTensorDecomposition(0)
      case("NNNfrom2N_spinorbit")
        call op2%init(rel,rel,opname)
        op2 = vnn%SpinTensorDecomposition(1)
      case("NNNfrom2N_tensor")
        call op2%init(rel,rel,opname)
        op2 = vnn%SpinTensorDecomposition(2)
      case("NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
        if(params%NNint%val /= "N3LO_EM500") then
          write(*,*) "Error: inconsitent NN interaction and operator name"
          stop
        end if
        call s%split(opname, s%str("NNNinduced_"), strings)
        input = params
        input%NNint = strings(2); input%renorm="bare"
        call relspin%init(params%hw, chjac%GetNmax(), params%J2max_NNint)
        call vnn_relspin%init(relspin)
        call U_relspin%init(relspin)
        call vnn_relspin%setNNForceHOIsospin(U_relspin,input)

        call op2%init(rel,rel,s%str("NNint"))
        call op2%SetTwoBodyScalarSpinBreaking(vnn_relspin)

        call vnn_relspin%fin()
        call U_relspin%fin()
        call relspin%fin()

      case("NNNfromTkin")
        call op2%init(rel,rel,s%str('kinetic'))
        call op2%set()
      case default
        write(*,*) "# Unknown operator at ", trim(opname%val), &
            & __FILE__, __LINE__
        return
      end select
      op2_sub = op2
      input = params
      call op2_sub%evolve(input, chjac%GetNmax(), chjac%GetNmax())
      fut = ut%GetFileName(chjac,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      if(params%hw_conversion) then
        fut = ut%GetFileName(chjac,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
      end if
      ex = s%isfile(fut, s%str('In set_three_body_jac_operator_isospin'))
      open(iut, file=fut%val, form='unformatted', access='stream')
      select case(opname%val)
      case("NNNfromTkin")
        call UT%init(chjac, chjac, s%str('UT'))
        call UT%readf(chjac, iut, chjac%GetNmax(), params%hw, params%hw_target, params%hw_conversion)

        call t_orig%init(chjac, chjac, s%str("kinetic"))
        call t_orig%set(op2)
        t_orig%DMat = UT%DMat%T() * t_orig%DMat * UT%DMat
        call UT%release()

        call t_eff%init(chjac, chjac, s%str("kinetic"))
        call t_eff%set(op2_sub)
        t_orig%DMat = t_orig%DMat - (t_eff%DMat * 1.5d0)
        call t_eff%release()

        call tmp%init(chjac, chjac, s%str("kinetic"))
        call tmp%set(op2)
        t_orig%DMat = t_orig%DMat + 0.5d0 * tmp%DMat
        call tmp%release()
        if(params%hw_conversion) call t_orig%FreqConv(chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target)
        v%m(:,:) = t_orig%m(:v%jacobi_ch_bra%GetNumberAStates(), :v%jacobi_ch_ket%GetNumberAStates())
        call t_orig%release()

      case("NNNinduced")
        call UT%init(chjac, chjac, s%str('UT'))
        call UT%readf(chjac, iut, chjac%GetNmax(), params%hw, params%hw_target, params%hw_conversion)
        call t_rel%init(rel,rel,s%str('kinetic'))
        call t_rel%set()
        ! -- 3N evolution
        call t_orig%init(chjac, chjac, s%str("kinetic"))
        call t_orig%set(t_rel)
        call vnn_jac%init(chjac, chjac, s%str("NNint"))
        call vnn_jac%set(op2)
        call tmp%init(chjac, chjac, s%str("hamil"))
        tmp%DMat = UT%DMat%T() * (t_orig%DMat + vnn_jac%DMat) * UT%DMat
        call vnn_jac%release()
        call vnn_jac%init(chjac, chjac, s%str("NNint"))
        !

        ! -- subtract
        op2_sub = op2 + t_rel
        call op2_sub%evolve(input, chjac%GetNmax(), chjac%GetNmax())
        op2_sub = op2_sub - t_rel
        call vnn_jac%set(op2_sub)
        tmp%DMat = tmp%DMat - vnn_jac%DMat - t_orig%DMat
        !
        call t_orig%release()
        call vnn_jac%release()
        if(params%hw_conversion) call tmp%FreqConv(chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target)
        v%m(:,:) = tmp%m(:v%jacobi_ch_bra%GetNumberAStates(), :v%jacobi_ch_ket%GetNumberAStates())
        call tmp%release()

      case default
        call v%set(chjac, chjac, iut, iut, op2, op2_sub, &
            & chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target, params%hw_conversion)
      end select
      close(iut)
      call vnn%fin()
      call op2%fin()
      call op2_sub%fin()
      call chjac%fin()
    end subroutine set_induced_nnn_ch
  end subroutine set_induced_nnn


  subroutine set_nnn_component(this, params)
    type(ThreeBodyJacOpIso), intent(inout), target :: this
    type(InputParameters), intent(in) :: params
    integer, allocatable :: slranks(:)
    real(8) :: time
#ifdef MPI
    integer :: ch, nelms
#endif
    call parent_child_procedure(set_nnn_component_ch, this%ms%GetNumberChannels(), slranks, time)
    if(nprocs>1) call timer%add(sy%str('MPI parent-child, nnn_force_component'), time)
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(slranks(1), this%ms%GetNumberChannels(), mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    do ch = 1, this%ms%GetNumberChannels()
      nelms = this%ms%jpt(ch)%GetNumberAStates()**2
      if(nelms < 1) cycle
      call mpi_bcast(this%MatCh(ch,ch)%m(1,1), nelms, mpi_real8, slranks(ch), mpi_comm_world, ierr)
    end do
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(slranks)
  contains
    subroutine set_nnn_component_ch(loop)
      use NNNForceHOIsospin
      integer, intent(in) :: loop
      type(ThreeBodyJacOpChanIso), pointer :: v
      type(str) :: opname
      type(ThreeBodyJacIsoChan) :: chjac
      type(ThreeBodyJacOpChanIso) :: tmp, ut
      type(DMat) :: v3n_NAS
      real(8) :: LECs(5)
      integer :: iut = 101, Nmax_file, iut_3nf = 102
      logical :: ex
      type(sys) :: s
      type(str) :: fut, fn_jacobi, fn

      v => this%MatCh(loop,loop)
      opname = this%GetOpName()

      if(.not. params%genuine_3bf) then
        write(*,*) "Warning: Operator", trim(opname%val), " cannot be calclated.", &
            & " Turn genuine_3nf true."
        return
      end if
      if(.not. params%save_3nf_before_lec) then
        write(*,*) "Warning: Operator", trim(opname%val), " cannot be calclated.", &
            & " Turn save_3nf_before_lec true."
        return
      end if
      if(v%jacobi_ch_ket%GetJ() > params%j3max_initial_3nf) return
      if(v%jacobi_ch_bra%GetNumberAStates() < 1) return
      if(v%jacobi_ch_ket%GetNumberAStates() < 1) return
      Nmax_file = GetRampNmax(v%jacobi_ch_ket%GetJ(), params%ramp)
      fn_jacobi = chjac%GetFileName(v%jacobi_ch_ket%GetJ(), &
          &  v%jacobi_ch_ket%GetParity(), v%jacobi_ch_ket%GetT(), Nmax_file, params%path_to_tmp_dir)
      ex = s%isfile(fn_jacobi, s%str('In set_three_body_jac_operator_isospin'))
      open(iut, file=fn_jacobi%val, form='unformatted', access='stream')
      call chjac%readf(params%hw,iut)
      close(iut)

      LECs = 0.d0
      select case(opname%val)
      case("NNN_c1")
        LECs(1) = params%c1
      case("NNN_c3")
        LECs(2) = params%c3
      case("NNN_c4")
        LECs(3) = params%c4
      case("NNN_TPE")
        LECs(1) = params%c1; LECs(2) = params%c3; LECs(3) = params%c4
      case("NNN_cD","NNN_OPE")
        LECs(4) = params%cd
      case("NNN_cE","NNN_Contact")
        LECs(5) = params%ce
      case("NNN_Genuine")
        LECs(1) = params%c1; LECs(2) = params%c3; LECs(3) = params%c4
        LECs(4) = params%cd; LECs(5) = params%ce
      case default
        write(*,*) "Unknown operator name:", trim(opname%val)
        return
      end select

      ! bare operator
      select case(params%Regulator%val)
      case("LNL","lnl","local-non-local","Local-Non-Local")
        if(params%NNNInt%val=="ChEFT_N2LO") then
          fn = get_file_name_ho_partial_wave_n2lo_isospin(chjac, Nmax_file, s%str("Local"), &
              & params%RegulatorPower, params%lambda_3nf_local, params%lambda_3nf_nonlocal, params%path_to_tmp_dir)
        else if (params%NNNInt%val=="ChEFT_N3LO") then
          fn = get_file_name_ho_partial_wave_n3lo_isospin(chjac, Nmax_file, s%str("Local"), &
              & params%RegulatorPower, params%lambda_3nf_local, params%lambda_3nf_nonlocal, params%path_to_tmp_dir)
        end if
      case default
        if(params%NNNInt%val=="ChEFT_N2LO") then
          fn = get_file_name_ho_partial_wave_n2lo_isospin(chjac, Nmax_file, params%Regulator, &
              & params%RegulatorPower, params%lambda_3nf_local, params%lambda_3nf_nonlocal, params%path_to_tmp_dir)
        else if (params%NNNInt%val=="ChEFT_N3LO") then
          fn = get_file_name_ho_partial_wave_n3lo_isospin(chjac, Nmax_file, params%Regulator, &
              & params%RegulatorPower, params%lambda_3nf_local, params%lambda_3nf_nonlocal, params%path_to_tmp_dir)
        end if
      end select
      call v3n_NAS%ini(chjac%GetNumberNAStates(), chjac%GetNumberNAStates())
      open(iut_3nf, file = fn%val, form='unformatted', access='stream',status='old')
      if(params%NNNInt%val=="ChEFT_N2LO") then
        call read_ho_partial_wave_n2lo(v3n_NAS, chjac%GetJ(), chjac%GetParity(), chjac%GetT(), &
            & chjac%GetNumberNAStates(), chjac%GetNmax(), iut_3nf, LECs)
      else if (params%NNNInt%val=="ChEFT_N3LO") then
        call read_ho_partial_wave_n3lo(v3n_NAS, chjac%GetJ(), chjac%GetParity(), chjac%GetT(), &
            & chjac%GetNumberNAStates(), chjac%GetNmax(), iut_3nf, LECs)
      end if
      close(iut_3nf)
      call tmp%Antisymmetrize(v3n_NAS)
      select case(params%Regulator%val)
      case("LNL","lnl","local-non-local","Local-Non-Local")
        call multiply_non_local_regulator_hospace(tmp%DMat, chjac, params)
      case default
      end select
      ! bare operator -> effective operator
      fut = ut%GetFileName(chjac,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      if(params%hw_conversion) then
        fut = ut%GetFileName(chjac,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
      end if
      ex = s%isfile(fut, s%str('In set_three_body_jac_operator_isospin'))
      open(iut, file=fut%val, form='unformatted', access='stream')
      call UT%init(chjac, chjac, s%str('UT'))
      call UT%readf(chjac, iut, chjac%GetNmax(), params%hw, &
          &         params%hw_target, params%hw_conversion)
      close(iut)
      tmp%DMat = UT%DMat%T() * tmp%DMat * UT%DMat
      call UT%release()
      if(params%hw_conversion) call tmp%FreqConv(chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target)
      v%m(:,:) = tmp%m(:v%jacobi_ch_bra%GetNumberAStates(), :v%jacobi_ch_ket%GetNumberAStates())
      call tmp%release()
      call chjac%fin()
    end subroutine set_nnn_component_ch
  end subroutine set_nnn_component


#define PRECISION Half
#define half_precision
#include "TMTransFunctions.inc"
#undef half_precision
#undef PRECISION

#define PRECISION Single
#define single_precision
#include "TMTransFunctions.inc"
#undef single_precision
#undef PRECISION

#define PRECISION Double
#define double_precision
#include "TMTransFunctions.inc"
#undef double_precision
#undef PRECISION
end module TMTransformManager
