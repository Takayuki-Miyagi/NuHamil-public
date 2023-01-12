#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION

subroutine MACRO_ADD_SUFX(trans_to_lab_basis_isospin_)(params)
  type(InputParameters), intent(in) :: params
  integer, allocatable :: jpt(:,:), slranks(:)
  integer :: n, j, p, t, nch, loop
  type(OrbitsIsospin) :: sps
  type(ThreeBodyLabIsoSpace) :: lab
  type(ThreeBodyJacIsoSpace) :: jac
  type(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)) :: jac2lab
  real(8) :: time_ms
  nch = ((params%jmax3 + 1) / 2) * 4
  allocate(jpt(3,nch), slranks(nch))
  n = 0
  do t = 1, 3, 2
    do j = 1, params%jmax3, 2
      do p = 1, -1, -2
        n = n + 1
        jpt(:,n) = [j,p,t]
      end do
    end do
  end do

  call parent_child_procedure(calc_each_channel, nch, slranks, time_ms)
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, three-body force'), time_ms)

  if(params%only_hf_monopole) then

    if(params%no2b_channel_begin == -1) call MACRO_ADD_SUFX(trans_to_lab_monopole_)(params)
    if(params%no2b_channel_begin > 0 ) call MACRO_ADD_SUFX(trans_to_lab_monopole_writing_)(params)
    if(params%no2b_channel_begin == -100 ) call MACRO_ADD_SUFX(merge_trans_to_lab_monopole_)(params)
    return

  end if

  if(params%only_no2b_elements) then
    do loop = 1, size(params%Operators)

      if(params%no2b_channel_begin == -1) then
        select case( params%Operators(loop)%val )
        case("hamil", "NNN_c1", "NNN_c3", "NNN_c4", "NNN_TPE", &
              & "NNN_cD", "NNN_OPE", "NNN_cE", "NNN_Contact", "NNN_Genuine", &
              & "NNNfrom2N", "NNNfrom2N_central", "NNNfrom2N_spinorbit", "NNNfrom2N_tensor", "NNNfromTkin",&
              & "NNNinduced", "NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
          call MACRO_ADD_SUFX(trans_to_lab_no2b_)(params, params%Operators(loop))
        case default
          call MACRO_ADD_SUFX(trans_to_lab_no2b_tensor_)(params, params%Operators(loop))
        end select
      end if
      if(params%no2b_channel_begin > 0 ) call MACRO_ADD_SUFX(trans_to_lab_no2b_writing_)(params, params%Operators(loop))
      if(params%no2b_channel_begin == -100 ) call MACRO_ADD_SUFX(merge_trans_to_lab_no2b_)(params, params%Operators(loop))

    end do
    return
  end if

  do loop = 1, size(params%Operators)

    call trans_to_lab_op(params%Operators(loop))

  end do

  call jac2lab%fin()
  call jac%fin()
  call lab%fin()
  call sps%fin()

contains
  subroutine calc_each_channel(ich)
    use NNNForceHOIsospin
    integer, intent(in) :: ich
    type(ThreeBodyJacIsoChan) :: jac
    type(sys) :: s
    type(ThreeBodyJacOpChanIso) :: U
    type(NNNForceIsospin) :: v3
    type(str) :: f, fv, fut
    integer :: wunit = 21, runit = 22, Nmax

    Nmax = GetRampNmax(jpt(1,ich), params%ramp)
    f = jac%GetFileName(jpt(1,ich),jpt(2,ich),jpt(3,ich),Nmax,params%path_to_tmp_dir)
    fv =  U%GetFileName(jac,params%hw,s%str('NNNint'),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)
    fut = U%GetFileName(jac,params%hw,s%str('UT'),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)

    write(*,'(a, i4, a, i2, a, i2, a, i2, a, i2)') &
        &  'myrank=', myrank, ',  :Calculating for J = ', &
        &   jpt(1,ich), '/2,   P = ', jpt(2,ich), ',   T = ', jpt(3,ich), &
        &   '/2, Nmax = ', Nmax

    if(.not. s%isfile(f)) then
      call jac%init(params%hw,jpt(1,ich),jpt(2,ich),jpt(3,ich),Nmax,params%path_to_tmp_dir)
      ! -- write to file
      open(wunit,file=f%val,status='replace',form='unformatted',access='stream')
      call jac%writef(wunit)
      close(wunit)
      call jac%fin()
    end if

    ! -- read from file
    open(runit,file=f%val,status='old',form='unformatted',access='stream')
    call jac%readf(params%hw,runit,Nmax)
    close(runit)

    if(jac%GetNumberNAStates() < 1 .or. jac%GetNumberAStates() < 1) then
      call jac%fin()
      return
    end if

    if(params%NN_only) then
      call jac%fin()
      write(*,'(a)') 'To calculate the 3BME, turn NN_only to .false.'
      return
    end if

    if(s%isfile(fv) .and. s%isfile(fut)) then
      call jac%fin()
      return
    end if
    call v3%InitNNNForce(jac)
    call U%init(jac,jac,s%str("UT"))
    call v3%SetNNNForce(U, params)

    call v3%fin()
    call U%fin()
    call jac%fin()
    !write(*,*)
  end subroutine calc_each_channel

  subroutine trans_to_lab_op(oprtr)
    type(str), intent(in) :: oprtr
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)) :: oplab
    type(ThreeBodyJacOpIso) :: opjac, optmp
    type(sys) :: s
    type(str) :: filename

    if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
    if(.not. lab%is_Constructed) call lab%init(sps, params%e2max, params%e3max, params%hw_target)

    call oplab%init(lab, oprtr)

    filename = oplab%GetFile(params%file_name_3n, params%NNInt, params%NNNInt, params%renorm, &
        & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
        & params%cd, params%ce, &
        & params%genuine_3bf, params%Regulator, params%RegulatorPower)

    if(s%isfile(filename)) then
      call oplab%fin()
      call lab%fin()
      call sps%fin()
      if(myrank==0) write(*,'(2a)') trim(filename%val), ' already exists.'
      return
    end if

    if(.not. jac%is_Constructed) call jac%init(params%hw,params%jmax3,params%e3max,params%ramp,params%N3max,params%path_to_tmp_dir)
    if(.not. jac2lab%is_Constructed) call jac2lab%init(jac, lab, params%e3max, params%jmax3, params%nblock)

    call opjac%init(jac, oprtr)
    select case(oprtr%val)
    case("NNNfrom2N", "NNNfrom2N_central", "NNNfrom2N_spinorbit", "NNNfrom2N_tensor", "NNNfromTkin","NNNinduced",&
          & "NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
      if(s%isfile(opjac%GetFileName(params))) then
        call opjac%set(params)
      else
        call set_induced_nnn(opjac, params)
        call opjac%writef(params)
        if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
      end if
    case("NNN_c1", "NNN_c3", "NNN_c4", "NNN_TPE", "NNN_cD", "NNN_OPE", "NNN_cE", "NNN_Contact", "NNN_Genuine")
      if(s%isfile(opjac%GetFileName(params))) then
        call opjac%set(params)
      else
        call set_nnn_component(opjac, params)
        call opjac%writef(params)
        if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
      end if
    case default
      if(s%isfile(opjac%GetFileName(params))) then
        call opjac%set(params)
      else
        call opjac%set(params)
        if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
      end if
    end select

    if(params%spin_tensor_decomposition /= -1) then
      select case(params%spin_tensor_decomposition)
      case(0, 1, 2, 3)
        optmp = opjac%SpinTensorDecomposition(params%spin_tensor_decomposition)
        opjac = optmp
      case(4) ! C + LS
        optmp = opjac%SpinTensorDecomposition(0) + opjac%SpinTensorDecomposition(1)
        opjac = optmp
      case(5) ! C + T
        optmp = opjac%SpinTensorDecomposition(0) + opjac%SpinTensorDecomposition(2)
        opjac = optmp
      case(6) ! C + k3
        optmp = opjac%SpinTensorDecomposition(0) + opjac%SpinTensorDecomposition(3)
        opjac = optmp
      case(7) ! C + LS + T
        optmp = opjac%SpinTensorDecomposition(0) + opjac%SpinTensorDecomposition(1) + opjac%SpinTensorDecomposition(2)
        opjac = optmp
      case(-100) ! C + LS + T + (k=3)
        optmp = opjac%SpinTensorDecomposition(0) + opjac%SpinTensorDecomposition(1) + &
            & opjac%SpinTensorDecomposition(2) + opjac%SpinTensorDecomposition(3)
        opjac = optmp
      case default
        write(*,*) 'Unknown value of spin_tensor_decomposition', __LINE__, __FILE__
      end select
    end if

    call oplab%TMtrans(jac2lab,opjac)
    if(myrank == 0) then
      if(params%averaged_file_for_test%val == "none") call oplab%writef(filename)
      if(params%averaged_file_for_test%val /= "none") call oplab%writef(filename, params%averaged_file_for_test)
    end if

    call opjac%fin()
    call oplab%fin()
  end subroutine trans_to_lab_op
end subroutine MACRO_ADD_SUFX(trans_to_lab_basis_isospin_)

subroutine MACRO_ADD_SUFX(trans_to_lab_monopole_)(params)
  type(InputParameters), intent(in) :: params
  type(OrbitsIsospin) :: sps
  type(MonThreeBodyLabIsoSpace) :: spmon
  type(ThreeBodyJacOpIso) :: opjac
  type(ThreeBodyJacIsoSpace) :: jac
  type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)) :: jac2lab
  type(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)) :: opmon
  real(8) :: time_ms
  type(sys) :: s
  type(str) :: filename
#ifdef MPI
#ifdef SPARC
  integer :: n1d
#else
  integer(8) :: n1d
#endif
  integer :: tag_channel=15
  integer :: tag_me_temp=16
#endif

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call spmon%init(sps, params%e2max, params%e3max)
  call opmon%init(spmon)

  filename = opmon%GetFileName(params%file_name_3n, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    call opmon%fin()
    call spmon%fin()
    call sps%fin()
    if(myrank==0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  if(.not. jac%is_Constructed) &
      & call jac%init(params%hw, params%jmax3, params%e3max, params%ramp, params%N3max, params%path_to_tmp_dir)
  call opjac%init(jac, s%str('NNNint'))
  if(s%isfile(opjac%GetFileName(params))) then
    call opjac%set(params)
  else
    call opjac%set(params)
    if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
  end if
  call jac2lab%store_couplings(jac,params%e3max)
#ifndef MPI
  write(*,*)
  write(*,'("############################################################")')
  write(*,'(4x, "doing transformaion to sp basis ...")')
  write(*,'(a,i4)') "    Number of threads for tcoef reordering: ", params%n_threads_tcoef
  write(*,'("############################################################")')
  write(*,*)
#endif

  call parent_child_procedure_data_transfer(trans_to_lab_monopole_channel, receive_data, spmon%GetNumberChannels(), time_ms)
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, TMtrans'), time_ms)

  call jac2lab%release_couplings()
  if(myrank == 0) call opmon%writef(filename)

  call opjac%fin()
contains
  subroutine trans_to_lab_monopole_channel(loop)
    integer, intent(in) :: loop
    type(MonThreeBodyLabIsoChan), pointer :: channel
    integer :: T, j1, p1, j2, p2, j3, p3, ch
#ifdef MPI
    integer :: n_states
#endif
#if defined(half_precision) && defined(MPI)
    integer(2), allocatable :: v(:)
#endif

    ch = spmon%GetSortedChannelIndex(loop)
    channel => spmon%GetChannel(ch)
    T = channel%GetT()
    J1 = channel%GetJ1()
    P1 = channel%GetP1()
    J2 = channel%GetJ2()
    P2 = channel%GetP2()
    J3 = channel%GetJ3()
    P3 = channel%GetP3()
    call jac2lab%init(channel, jac, params%e3max, params%n_threads_tcoef)

    write(*,"(a,i4)",advance='no') "myrank=", myrank
    write(*,"(a,i3)",advance='no') ", T=", T
    write(*,"(a,i3)",advance='no') "/2, J1=", J1
    write(*,"(a,i3)",advance='no') "/2, P1=", P1
    write(*,"(a,i3)",advance='no') ", J2=", J2
    write(*,"(a,i3)",advance='no') "/2, P2=", P2
    write(*,"(a,i3)",advance='no') ", J3=", J3
    write(*,"(a,i3)",advance='no') "/2, P3=", P3
    write(*,"(a,f10.4)",advance='no') ", Mem(ME): ", opmon%MatCh(ch)%GetMemory()*1.d3
    write(*,"(a,f10.4)",advance='no') " MB, Mem(Tc): ", jac2lab%GetMemory()
    write(*,"(a,i8,a,i8)") " GB, ", loop, "/", spmon%GetNumberChannels()

    call opmon%MatCh(ch)%trans(jac2lab, opjac)
#ifdef MPI
    if(myrank /= 0) then
      n_states = channel%GetNumberStates()
      n1d = int( n_states, kind(n1d)) * int( n_states+1, kind(n1d) ) / int(2,kind(n1d))
      call mpi_send(ch,1,mpi_integer,0,tag_channel,mpi_comm_world,ierr)
#if defined(half_precision)
      allocate(v(n1d))
      v(:) = opmon%MatCh(ch)%v(:)%i16
      call mpi_send(v(1), n1d, mpi_integer2, 0, tag_me_temp, mpi_comm_world, ierr)
      deallocate(v)
#elif defined(single_precision)
      call mpi_send(opmon%MatCh(ch)%v(1), n1d, mpi_real4, 0, tag_me_temp, mpi_comm_world, ierr)
#elif defined(double_precision)
      call mpi_send(opmon%MatCh(ch)%v(1), n1d, mpi_real8, 0, tag_me_temp, mpi_comm_world, ierr)
#endif
    end if
#endif
    call jac2lab%fin()
  end subroutine trans_to_lab_monopole_channel

  subroutine receive_data()
#ifdef MPI
    integer :: ch, n_states
    real(8) :: t1
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif
#if defined(half_precision) && defined(MPI)
    integer(2), allocatable :: v(:)
#endif
    call mpi_recv(ch, 1, mpi_integer, mpi_any_source, tag_channel,mpi_comm_world, istatus, ierr)
    t1 = omp_get_wtime()
    n_states = spmon%Chan(ch)%GetNumberStates()
    n1d = int( n_states, kind(n1d)) * int( n_states+1, kind(n1d) ) / int(2,kind(n1d))
#if defined(half_precision)
    allocate(v(n1d))
    call mpi_recv(v(1), n1d, mpi_integer2, istatus(mpi_source), tag_me_temp, mpi_comm_world, istatus, ierr)
    opmon%MatCh(ch)%v(:)%i16 = v(:)
    deallocate(v)
#elif defined(single_precision)
    call mpi_recv(opmon%MatCh(ch)%v(1),n1d,mpi_real4,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#elif defined(double_precision)
    call mpi_recv(opmon%MatCh(ch)%v(1),n1d,mpi_real8,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#endif
    call timer%add(sy%str("Gathering data to parent node"), omp_get_wtime()-t1)
#endif
  end subroutine receive_data
end subroutine MACRO_ADD_SUFX(trans_to_lab_monopole_)

subroutine MACRO_ADD_SUFX(merge_trans_to_lab_monopole_)(params)
  type(InputParameters), intent(in) :: params
  type(OrbitsIsospin) :: sps
  type(MonThreeBodyLabIsoSpace) :: spmon
  type(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)) :: opmon
  type(sys) :: s
  type(str) :: filename
  integer :: i_ch

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call spmon%init(sps, params%e2max, params%e3max)
  if(myrank == 0) then
    call opmon%init(spmon)
    write(*,"(a,f12.6,a)") "# HF relevant matrix elements: ", opmon%GetMemory(), " GB"
  end if

  filename = opmon%GetFileName(params%file_name_3n, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    if(myrank == 0) call opmon%fin()
    call spmon%fin()
    call sps%fin()
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  do i_ch = 1, spmon%GetNumberChannels()
    call merge_channel(i_ch)
  end do
  if(myrank == 0) call opmon%writef(filename)
contains
  subroutine merge_channel(ich_in)
    integer, intent(in) :: ich_in
    integer :: loop
    type(MonThreeBodyLabIsoChan), pointer :: channel
    integer :: runit=22
    type(sys) :: s
    type(str) :: fn

    loop = spmon%GetSortedChannelIndex(ich_in)
    channel => spmon%GetChannel(loop)

    write(*,"(a,i4)",advance='no') "myrank=", myrank
    write(*,"(a,i3)",advance='no') ", T=", channel%GetT()
    write(*,"(a,i3)",advance='no') "/2, J1=", channel%GetJ1()
    write(*,"(a,i3)",advance='no') "/2, P1=", channel%GetP1()
    write(*,"(a,i3)",advance='no') ", J2=", channel%GetJ2()
    write(*,"(a,i3)",advance='no') "/2, P2=", channel%GetP2()
    write(*,"(a,i3)",advance='no') ", J3=", channel%GetJ3()
    write(*,"(a,i3)",advance='no') "/2, P3=", channel%GetP3()
    write(*,"(a,f10.4)",advance='no') ", Mem(ME): ", opmon%MatCh(loop)%GetMemory()*1.d3
    write(*,"(a,i8,a,i8)") " MB, ", ich_in, "/", spmon%GetNumberChannels()

    fn = get_file_name(ich_in)
    if(.not. s%isfile( fn )) then
      write(*,*) "Error, file doesn't exist: ", fn%val
      stop
    end if
    open(runit, file=fn%val, form="unformatted", access="stream")
    read(runit) opmon%MatCh(loop)%v(:)
    close(runit)
  end subroutine merge_channel

  function get_file_name(ch) result(f)
    integer, intent(in) :: ch
    type(str) :: f
    type(sys) :: s
    f = params%path_to_tmp_dir + s%str("/") + params%no2b_temp_dir + s%str("/monopole_temp_ch") + s%str(ch) + s%str(".bin")
  end function get_file_name
end subroutine MACRO_ADD_SUFX(merge_trans_to_lab_monopole_)

subroutine MACRO_ADD_SUFX(trans_to_lab_monopole_writing_)(params)
  type(InputParameters), intent(in) :: params
  type(OrbitsIsospin) :: sps
  type(ThreeBodyJacOpIso) :: opjac
  type(ThreeBodyJacIsoSpace) :: jac
  type(MonThreeBodyLabIsoSpace) :: spmon
  type(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)) :: opmon
  type(sys) :: s
  type(str) :: filename
  real(8) :: time_ms
  type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)) :: jac2lab
  integer :: i_ch

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call spmon%init(sps, params%e2max, params%e3max)
  call opmon%init(spmon)
  if( params%no2b_channel_begin > spmon%GetNumberChannels() ) return

  filename = opmon%GetFileName(params%file_name_3n, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    call spmon%fin()
    call sps%fin()
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  if(.not. jac%is_Constructed) &
      & call jac%init(params%hw, params%jmax3, params%e3max, params%ramp, params%N3max, params%path_to_tmp_dir)
  call opjac%init(jac, s%str('NNNint'))
  if(s%isfile(opjac%GetFileName(params))) then
    call opjac%set(params)
  else
    call opjac%set(params)
    if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
  end if
  call jac2lab%store_couplings(jac,params%e3max)
  if(myrank == 0) then
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'(a,i4)') "    Number of threads for tcoef reordering: ", params%n_threads_tcoef
    write(*,'("############################################################")')
    write(*,*)
  end if

  do i_ch = params%no2b_channel_begin, min( params%no2b_channel_end, spmon%GetNumberChannels() )
    call trans_to_lab_monopole_channel_writing(i_ch)
  end do
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, TMtrans'), time_ms)
  call jac2lab%release_couplings()
  call opjac%fin()
  call spmon%fin()
contains
  subroutine trans_to_lab_monopole_channel_writing(ich_in)
    integer, intent(in) :: ich_in
    integer :: loop
    type(MonThreeBodyLabIsoChan), pointer :: channel
    real(8) :: t1, t2
    integer :: wunit = 25
    type(str) :: fn
    type(sys) :: s

    loop = spmon%GetSortedChannelIndex(ich_in)
    if( s%isfile(get_file_name(ich_in)) ) return
    channel => spmon%GetChannel(loop)
    t1 = omp_get_wtime()
    call jac2lab%init(channel, jac, params%e3max, params%n_threads_tcoef)
    t2 = omp_get_wtime()

    write(*,"(a,i4)",advance='no') "myrank=", myrank
    write(*,"(a,i3)",advance='no') ", T=", channel%GetT()
    write(*,"(a,i3)",advance='no') "/2, J1=", channel%GetJ1()
    write(*,"(a,i3)",advance='no') "/2, P1=", channel%GetP1()
    write(*,"(a,i3)",advance='no') ", J2=", channel%GetJ2()
    write(*,"(a,i3)",advance='no') "/2, P2=", channel%GetP2()
    write(*,"(a,i3)",advance='no') ", J3=", channel%GetJ3()
    write(*,"(a,i3)",advance='no') "/2, P3=", channel%GetP3()
    write(*,"(a,f10.4)",advance='no') ", Mem(ME): ", opmon%MatCh(loop)%GetMemory()*1.d3
    write(*,"(a,f10.4)",advance='no') " MB, Mem(Tc): ", jac2lab%GetMemory()
    write(*,"(a,i8,a,i8)") " GB, ", ich_in, "/", spmon%GetNumberChannels()

    call opmon%MatCh(loop)%trans(jac2lab, opjac)
    fn = get_file_name(ich_in)
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) opmon%MatCh(loop)%v(:)
    close(wunit)
    call jac2lab%fin()
  end subroutine trans_to_lab_monopole_channel_writing

  function get_file_name(ch) result(f)
    integer, intent(in) :: ch
    type(str) :: f, dir
    type(sys) :: s
    dir = params%path_to_tmp_dir + s%str("/") + params%no2b_temp_dir
    call s%mkdir(dir%val)
    f = dir + s%str("/monopole_temp_ch") + s%str(ch) + s%str(".bin")
  end function get_file_name
end subroutine MACRO_ADD_SUFX(trans_to_lab_monopole_writing_)

subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_)(params, opname)
  type(InputParameters), intent(in) :: params
  type(str), intent(in) :: opname
  type(OrbitsIsospin) :: sps
  type(NO2BThreeBodyIsoSpace) :: spno2b
  type(ThreeBodyJacOpIso) :: opjac, optmp
  type(ThreeBodyJacIsoSpace) :: jac
  type(ThreeBodyLabIsoSpaceNO2B) :: lab_no2b
  type(MACRO_ADD_SUFX(ThreeBodyNO2BIso)) :: opno2b
  type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)) :: jac2lab
  type(sys) :: s
  type(str) :: filename
  real(8) :: time_ms
#ifdef MPI
  integer :: tag_channel=15
  integer :: tag_me_temp=16
#endif

  filename = opno2b%GetFileName(params%file_name_3n, opname, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  if(.not. jac%is_Constructed) &
      & call jac%init(params%hw, params%jmax3, params%e3max, params%ramp, params%N3max, params%path_to_tmp_dir)
  if(myrank == 0) then
    write(*,"(a,f12.6,a)") "# 3N Jacobi space: ", jac%GetMemory(), " GB"
  end if

  call opjac%init(jac, OpName)

  select case(opname%val)
  case("NNNfrom2N", "NNNfrom2N_central", "NNNfrom2N_spinorbit", "NNNfrom2N_tensor", "NNNfromTkin","NNNinduced",&
          & "NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call set_induced_nnn(opjac, params)
      call opjac%writef(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case("NNN_c1", "NNN_c3", "NNN_c4", "NNN_TPE", "NNN_cD", "NNN_OPE", "NNN_cE", "NNN_Contact", "NNN_Genuine")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call set_nnn_component(opjac, params)
      call opjac%writef(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case("hamil")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call opjac%set(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case default
    if(myrank==0) write(*,*) "Not supported", __LINE__, " in ", __FILE__
    return
  end select

  if(params%spin_tensor_decomposition /= -1) then
    optmp = opjac%SpinTensorDecomposition(params%spin_tensor_decomposition)
    opjac = optmp
  end if

  call jac2lab%store_couplings(jac,params%e3max)
  if(myrank == 0) then
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'(a,i4)') "    Number of threads for tcoef reordering: ", params%n_threads_tcoef
    write(*,'("############################################################")')
    write(*,*)
  end if

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call lab_no2b%init(sps, params%e2max, params%e3max)
  call spno2b%init(sps, params%e2max, params%e3max)
  if(myrank == 0) then
    call opno2b%init(spno2b)
    write(*,"(a,f12.6,a)") "# NO2B relevant matrix elements: ", opno2b%GetMemory(), " GB"
  end if

  filename = opno2b%GetFileName(params%file_name_3n, opname, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  call parent_child_procedure_data_transfer(trans_to_lab_no2b_channel, receive_data, lab_no2b%GetNumberChannels(), time_ms)
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, TMtrans: '), time_ms)
  call jac2lab%release_couplings()
  if(myrank == 0) call opno2b%writef(filename)
  call opjac%fin()
  if(myrank == 0) call lab_no2b%fin()
contains
  subroutine trans_to_lab_no2b_channel(ich_in)
    integer, intent(in) :: ich_in
    integer :: loop
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)) :: oplab_ch
    type(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)) :: oplab_ch_no2b
    type(ThreeBodyLabIsoChanNO2B), pointer :: ch_original
    type(NO2BThreeBodyIsoChan), pointer :: ch_averaged
    integer :: J, P, T, J12, P12, J3, P3, ch
    real(8) :: t1, t2, t3, t4
#ifdef MPI
    integer :: n_states
#ifdef half_precision
    integer(2), allocatable :: tmp(:)
#endif /* half_precision */
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif /* SPARC */
#endif /* MPI */

    loop = lab_no2b%GetSortedChannelIndex(ich_in)
    ch_original => lab_no2b%GetChannel(loop)
    J = ch_original%GetJ()
    P = ch_original%GetParity()
    T = ch_original%GetT()
    J12 = ch_original%GetJ12()
    P12 = ch_original%GetP12()
    J3  = ch_original%GetJ3()
    P3  = ch_original%GetP3()
    ch = spno2b%GetIndex(T,J12,P12,J3,P3)
    ch_averaged => spno2b%GetChannel(ch)

    call oplab_ch%init( ch_original, ch_original )
    call oplab_ch_no2b%init( ch_averaged, ch_averaged )
    t1 = omp_get_wtime()
    call jac2lab%init(ch_original, jac, params%e3max, params%n_threads_tcoef)
    t2 = omp_get_wtime()

    write(*,"(a,i4,a,i3,a,i3,a,i2,a,i3,a,i3,a,i3,a,i3,a,f10.4,a,f10.4,a,i8,a,i8)") &
        & "myrank=", myrank, ", J=", J, "/2, P=", P, ", T=", T, "/2, J12=",&
        & J12, ", P12=", P12, ", J3=", J3, "/2, P3=", P3, ", Mem(ME): ", &
        & oplab_ch_no2b%GetMemory()*1.d3, " MB, Mem(Tc): ", jac2lab%GetMemory(), &
        & " GB, ", ich_in, "/", lab_no2b%GetNumberChannels()

    call oplab_ch%trans(jac2lab, opjac)
    t3 = omp_get_wtime()
    call oplab_ch_no2b%AngleAverage(oplab_ch)
    t4 = omp_get_wtime()

#ifdef MPI
    if(myrank /= 0) then
      n_states = ch_averaged%GetNumberStates()
      n1d = int( n_states, kind(n1d)) * int( n_states+1, kind(n1d) ) / int(2,kind(n1d))
      call mpi_send(ch,1,mpi_integer,0,tag_channel,mpi_comm_world,ierr)
#if defined(half_precision)
      allocate( tmp( n1d ))
      tmp(:) = oplab_ch_no2b%v(:)%i16
      call mpi_send(tmp(1), n1d, mpi_integer2, 0, tag_me_temp, mpi_comm_world, ierr)
      deallocate(tmp)
#elif defined(single_precision)
      call mpi_send(oplab_ch_no2b%v(1), n1d, mpi_real4, 0, tag_me_temp, mpi_comm_world, ierr)
#elif defined(double_precision)
      call mpi_send(oplab_ch_no2b%v(1), n1d, mpi_real8, 0, tag_me_temp, mpi_comm_world, ierr)
#endif /* precision */
    end if
#endif /* MPI */
    if(myrank == 0) then
      opno2b%MatCh(ch)%v(:) = opno2b%MatCh(ch)%v(:) + oplab_ch_no2b%v(:)
    end if

    call jac2lab%fin()
    call oplab_ch%fin()
    call oplab_ch_no2b%fin()
  end subroutine trans_to_lab_no2b_channel

  subroutine receive_data()
#ifdef MPI
    integer :: ch, n_states
    real(8) :: t1
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif /* SPARC */
#if defined(half_precision)
    integer(2), allocatable :: tmp(:)
    type(my_real16), allocatable :: tmp16(:)
#elif defined(single_precision)
    real(4), allocatable :: tmp(:)
#elif defined(double_precision)
    real(8), allocatable :: tmp(:)
#endif /* precision */
    call mpi_recv(ch,1,mpi_integer,mpi_any_source,tag_channel,mpi_comm_world,istatus,ierr)
    t1 = omp_get_wtime()
    n_states = spno2b%Chan(ch)%GetNumberStates()
    n1d = int( n_states, kind(n1d)) * int( n_states+1, kind(n1d) ) / int(2,kind(n1d))
    allocate(tmp(n1d))
#if defined(half_precision)
    allocate(tmp16(n1d))
    call mpi_recv(tmp,n1d,mpi_real4,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
    tmp16(:)%i16 = tmp(:)
#elif defined(single_precision)
    call mpi_recv(tmp,n1d,mpi_real4,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#elif defined(double_precision)
    call mpi_recv(tmp,n1d,mpi_real8,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#endif /* precision */
#if defined(half_precision)
    opno2b%MatCh(ch)%v(:) = opno2b%MatCh(ch)%v(:) + tmp16(:)
    deallocate(tmp16)
#elif defined(single_precision) || defined(double_precision)
    opno2b%MatCh(ch)%v(:) = opno2b%MatCh(ch)%v(:) + tmp(:)
#endif /* precision */
    deallocate(tmp)
    call timer%add(sy%str("Gathering data to parent node"), omp_get_wtime()-t1)
#endif /* MPI */
  end subroutine receive_data
end subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_)

subroutine MACRO_ADD_SUFX(merge_trans_to_lab_no2b_)(params, opname)
  type(InputParameters), intent(in) :: params
  type(str), intent(in) :: opname
  type(OrbitsIsospin) :: sps
  type(ThreeBodyLabIsoSpaceNO2B) :: lab_no2b
  type(NO2BThreeBodyIsoSpace) :: spno2b
  type(MACRO_ADD_SUFX(ThreeBodyNO2BIso)) :: opno2b
  type(sys) :: s
  type(str) :: filename
  integer :: i_ch

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call lab_no2b%init(sps, params%e2max, params%e3max)
  call spno2b%init(sps, params%e2max, params%e3max)
  if(myrank == 0) then
    call opno2b%init(spno2b)
    write(*,"(a,f12.6,a)") "# NO2B relevant matrix elements: ", opno2b%GetMemory(), " GB"
  end if

  filename = opno2b%GetFileName(params%file_name_3n, opname, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    if(myrank == 0) call opno2b%fin()
    call lab_no2b%fin()
    call spno2b%fin()
    call sps%fin()
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  do i_ch = 1, lab_no2b%GetNumberChannels()
    call merge_channel(i_ch)
  end do
  if(myrank == 0) call opno2b%writef(filename)
  if(myrank == 0) call lab_no2b%fin()
contains
  subroutine merge_channel(ich_in)
    integer, intent(in) :: ich_in
    integer :: loop
    type(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)) :: oplab_ch_no2b
    type(ThreeBodyLabIsoChanNO2B), pointer :: ch_original
    type(NO2BThreeBodyIsoChan), pointer :: ch_averaged
    integer :: J, P, T, J12, P12, J3, P3, ch
    !real(8) :: t1, t2, t3, t4
    integer :: runit=22
    type(sys) :: s
    type(str) :: fn

    loop = lab_no2b%GetSortedChannelIndex(ich_in)
    ch_original => lab_no2b%GetChannel(loop)
    J = ch_original%GetJ()
    P = ch_original%GetParity()
    T = ch_original%GetT()
    J12 = ch_original%GetJ12()
    P12 = ch_original%GetP12()
    J3  = ch_original%GetJ3()
    P3  = ch_original%GetP3()
    ch = spno2b%GetIndex(T,J12,P12,J3,P3)
    ch_averaged => spno2b%GetChannel(ch)

    call oplab_ch_no2b%init( ch_averaged, ch_averaged )

    write(*,"(a,i4,a,i3,a,i3,a,i2,a,i3,a,i3,a,i3,a,i3,a,f10.4,a,i8,a,i8)") &
        & "myrank=", myrank, ", J=", J, "/2, P=", P, ", T=", T, "/2, J12=",&
        & J12, ", P12=", P12, ", J3=", J3, "/2, P3=", P3, ", Mem(ME): ", &
        & oplab_ch_no2b%GetMemory()*1.d3, " MB, ", ich_in, "/", lab_no2b%GetNumberChannels()

    fn = get_file_name(ich_in)
    if(.not. s%isfile( fn )) then
      write(*,*) "Error, file doesn't exist: ", fn%val
      stop
    end if
    open(runit, file=fn%val, form="unformatted", access="stream")
    read(runit) oplab_ch_no2b%v(:)
    close(runit)
    opno2b%MatCh(ch)%v(:) = opno2b%MatCh(ch)%v(:) + oplab_ch_no2b%v(:)
    call oplab_ch_no2b%fin()
  end subroutine merge_channel

  function get_file_name(ch) result(f)
    integer, intent(in) :: ch
    type(str) :: f
    type(sys) :: s
    f = params%path_to_tmp_dir + s%str("/") + params%no2b_temp_dir + s%str("/no2b_temp_ch") + s%str(ch) + s%str(".bin")
  end function get_file_name
end subroutine MACRO_ADD_SUFX(merge_trans_to_lab_no2b_)

subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_writing_)(params, opname)
  type(InputParameters), intent(in) :: params
  type(str), intent(in) :: opname
  type(OrbitsIsospin) :: sps
  type(ThreeBodyJacOpIso) :: opjac
  type(ThreeBodyJacIsoSpace) :: jac
  type(ThreeBodyLabIsoSpaceNO2B) :: lab_no2b
  type(NO2BThreeBodyIsoSpace) :: spno2b
  type(MACRO_ADD_SUFX(ThreeBodyNO2BIso)) :: opno2b
  type(sys) :: s
  type(str) :: filename
  real(8) :: time_ms
  type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)) :: jac2lab
  integer :: i_ch

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call lab_no2b%init(sps, params%e2max, params%e3max)
  call spno2b%init(sps, params%e2max, params%e3max)
  if( params%no2b_channel_begin > lab_no2b%GetNumberChannels() ) return

  filename = opno2b%GetFileName(params%file_name_3n, opname, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, &
      & params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    call lab_no2b%fin()
    call spno2b%fin()
    call sps%fin()
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  if(.not. jac%is_Constructed) &
      & call jac%init(params%hw, params%jmax3, params%e3max, params%ramp, params%N3max, params%path_to_tmp_dir)
  call opjac%init(jac, opname)
  select case(opname%val)
  case("NNNfrom2N", "NNNfrom2N_central", "NNNfrom2N_spinorbit", "NNNfrom2N_tensor", "NNNfromTkin","NNNinduced",&
        & "NNNinduced_N3LO_EM500_OPE", "NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call set_induced_nnn(opjac, params)
      call opjac%writef(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case("NNN_c1", "NNN_c3", "NNN_c4", "NNN_TPE", "NNN_cD", "NNN_OPE", "NNN_cE", "NNN_Contact", "NNN_Genuine")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call set_nnn_component(opjac, params)
      call opjac%writef(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case("hamil", "Hamil")
    if(s%isfile(opjac%GetFileName(params))) then
      call opjac%set(params)
    else
      call opjac%set(params)
      if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
    end if
  case default
    if(myrank==0) write(*,*) "unsupported operator:", trim(opname%val)
    return
  end select
  call jac2lab%store_couplings(jac,params%e3max)
  if(myrank == 0) then
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'(a,i4)') "    Number of threads for tcoef reordering: ", params%n_threads_tcoef
    write(*,'("############################################################")')
    write(*,*)
  end if

  do i_ch = params%no2b_channel_begin, min( params%no2b_channel_end, lab_no2b%GetNumberChannels() )
    call trans_to_lab_no2b_channel_writing(i_ch)
  end do
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, TMtrans'), time_ms)
  call jac2lab%release_couplings()
  call opjac%fin()
  call spno2b%fin()
  if(myrank == 0) call lab_no2b%fin()
contains
  subroutine trans_to_lab_no2b_channel_writing(ich_in)
    integer, intent(in) :: ich_in
    integer :: loop
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)) :: oplab_ch
    type(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)) :: oplab_ch_no2b
    type(ThreeBodyLabIsoChanNO2B), pointer :: ch_original
    type(NO2BThreeBodyIsoChan), pointer :: ch_averaged
    integer :: J, P, T, J12, P12, J3, P3
    integer :: ch
    real(8) :: t1, t2
    integer :: wunit = 25
    type(str) :: fn
!#ifdef half_precision
!    integer(2), allocatable :: tmp(:)
!#endif
    type(sys) :: s
    loop = lab_no2b%GetSortedChannelIndex(ich_in)
    if( s%isfile(get_file_name(ich_in)) ) return
    ch_original => lab_no2b%GetChannel(loop)
    J = ch_original%GetJ()
    P = ch_original%GetParity()
    T = ch_original%GetT()
    J12 = ch_original%GetJ12()
    P12 = ch_original%GetP12()
    J3  = ch_original%GetJ3()
    P3  = ch_original%GetP3()
    ch = spno2b%GetIndex(T,J12,P12,J3,P3)
    ch_averaged => spno2b%GetChannel(ch)

    call oplab_ch%init( ch_original, ch_original )
    call oplab_ch_no2b%init( ch_averaged, ch_averaged )
    t1 = omp_get_wtime()
    call jac2lab%init(ch_original, jac, params%e3max, params%n_threads_tcoef)
    t2 = omp_get_wtime()

    write(*,"(a,i4,a,i3,a,i3,a,i2,a,i3,a,i3,a,i3,a,i3,a,f10.4,a,i8,a,i8)") &
        & "myrank=", myrank, ", J=", J, "/2, P=", P, ", T=", T, "/2, J12=",&
        & J12, ", P12=", P12, ", J3=", J3, "/2, P3=", P3, &
        & ", Mem(Tc): ", jac2lab%GetMemory(), &
        & " GB, ", ich_in, "/", lab_no2b%GetNumberChannels()

    call oplab_ch%trans(jac2lab, opjac)
    call oplab_ch_no2b%AngleAverage(oplab_ch)
    fn = get_file_name(ich_in)
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) oplab_ch_no2b%v(:)
    close(wunit)
    call jac2lab%fin()
    call oplab_ch_no2b%fin()
    call oplab_ch%fin()
  end subroutine trans_to_lab_no2b_channel_writing

  function get_file_name(ch) result(f)
    integer, intent(in) :: ch
    type(str) :: f, dir
    type(sys) :: s
    dir = params%path_to_tmp_dir + s%str("/") + params%no2b_temp_dir
    call s%mkdir(dir%val)
    f = dir + s%str("/no2b_temp_ch") + s%str(ch) + s%str(".bin")
  end function get_file_name
end subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_writing_)

subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_tensor_)(params, opname)
  type(InputParameters), intent(in) :: params
  type(str), intent(in) :: opname
  type(OrbitsIsospin) :: sps
  type(NO2BThreeBodyIsoSpace) :: spno2b
  type(ThreeBodyJacOpIso) :: opjac
  type(ThreeBodyJacIsoSpace) :: jac
  type(ThreeBodyLabIsoSpaceNO2B) :: lab_no2b
  type(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)) :: opno2b
  type(sys) :: s
  type(str) :: filename
  real(8) :: time_ms
#ifdef MPI
  integer :: tag_channel=15
  integer :: tag_me_temp=16
#endif
  type(ThreeBodyLabIsoChanNO2B), pointer :: chbra, chket
  integer :: cnt, nloops, ichbra, ichket
  integer, allocatable :: ibras(:), ikets(:)
  type(OperatorDef) :: op_def
  type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)) :: tr_bra, tr_ket

  filename = opno2b%GetFileName(params%file_name_3n, opname, params%NNInt, params%NNNInt, params%renorm, &
      & params%lambda, params%hw_target, params%emax, params%e2max, params%e3max, &
      & params%cd, params%ce, params%genuine_3bf, params%Regulator, params%RegulatorPower)

  if(s%isfile(filename)) then
    if(myrank == 0) write(*,'(2a)') trim(filename%val), ' already exists.'
    return
  end if

  if(.not. jac%is_Constructed) &
      & call jac%init(params%hw, params%jmax3, params%e3max, params%ramp, params%N3max, params%path_to_tmp_dir)
  if(myrank == 0) then
    write(*,"(a,f12.6,a)") "# 3N Jacobi space: ", jac%GetMemory(), " GB"
  end if

  call opjac%init(jac, OpName)
  if(s%isfile(opjac%GetFileName(params))) then
    call opjac%set(params)
  else
    call opjac%set(params)
    if(params%fname_jacobi_op_for_trans2lab%val /= "none") return
  end if
  call tr_ket%store_couplings(jac,params%e3max)

  if(myrank == 0) then
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'(a,i4)') "    Number of threads for tcoef reordering: ", min(params%n_threads_tcoef,omp_get_max_threads())
    write(*,'("############################################################")')
    write(*,*)
  end if

  if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
  call lab_no2b%init(sps, params%e2max, params%e3max)
  call spno2b%init(sps, params%e2max, params%e3max)
  if(myrank == 0) then
    call opno2b%init(spno2b, OpName)
    write(*,"(a,f12.6,a)") "# NO2B relevant matrix elements: ", opno2b%GetMemory(), " GB"
  end if
  call op_def%InitOpDef( OpName, .false. )

  cnt = 0
  do ichbra = 1, lab_no2b%GetNumberChannels()
    chbra => lab_no2b%GetChannel( lab_no2b%GetSortedChannelIndex(ichbra))
    do ichket = 1, ichbra
      chket => lab_no2b%GetChannel( lab_no2b%GetSortedChannelIndex(ichket))
      if(triag( chbra%GetJ(), chket%GetJ(), 2*op_def%GetOpJ() ) ) cycle
      if(triag( chbra%GetJ12(), chket%GetJ12(), op_def%GetOpJ() ) ) cycle
      if(triag( chbra%GetT(), chket%GetT(), 2*op_def%GetOpT() ) ) cycle
      if(chbra%GetParity() * chket%GetParity() * op_def%GetOpP()==-1 ) cycle
      if(chbra%GetJ3() /= chket%GetJ3() ) cycle
      if(chbra%GetP3() /= chket%GetP3() ) cycle
      cnt = cnt + 1
    end do
  end do
  nloops = cnt
  allocate(ibras(nloops))
  allocate(ikets(nloops))
  cnt = 0
  do ichbra = 1, lab_no2b%GetNumberChannels()
    chbra => lab_no2b%GetChannel( lab_no2b%GetSortedChannelIndex(ichbra))
    do ichket = 1, ichbra
      chket => lab_no2b%GetChannel( lab_no2b%GetSortedChannelIndex(ichket))
      if(triag( chbra%GetJ(), chket%GetJ(), 2*op_def%GetOpJ() ) ) cycle
      if(triag( chbra%GetJ12(), chket%GetJ12(), op_def%GetOpJ() ) ) cycle
      if(triag( chbra%GetT(), chket%GetT(), 2*op_def%GetOpT() ) ) cycle
      if(chbra%GetParity() * chket%GetParity() * op_def%GetOpP()==-1 ) cycle
      if(chbra%GetJ3() /= chket%GetJ3() ) cycle
      if(chbra%GetP3() /= chket%GetP3() ) cycle
      cnt = cnt + 1
      ibras(cnt) = lab_no2b%GetSortedChannelIndex(ichbra)
      ikets(cnt) = lab_no2b%GetSortedChannelIndex(ichket)
    end do
  end do
  call parent_child_procedure_data_transfer(trans_to_lab_no2b_channel, receive_data, nloops, time_ms)
  if(nprocs>1) call timer%add(sy%str('MPI parent-child, TMtrans: '), time_ms)
  call tr_ket%release_couplings()
  if(myrank == 0) call opno2b%writef(filename)
  call opjac%fin()
  if(myrank == 0) call lab_no2b%fin()
contains
  subroutine trans_to_lab_no2b_channel(loop)
    integer, intent(in) :: loop
    type(MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)) :: oplab_ch
    type(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)) :: oplab_ch_no2b
    type(ThreeBodyLabIsoChanNO2B), pointer :: chbra_original, chket_original
    type(NO2BThreeBodyIsoChan), pointer :: chbra_averaged, chket_averaged
    integer :: Jbra, Pbra, Tbra, J12bra, P12bra, J3bra, P3bra, chbra
    integer :: Jket, Pket, Tket, J12ket, P12ket, J3ket, P3ket, chket
    real(8) :: t1, t2, t3, t4
#if defined(half_precision)
    real(4) :: flip
#elif defined(single_precision)
    real(4) :: flip
#elif defined(double_precision)
    real(8) :: flip
#endif
#ifdef MPI
    integer :: chbraket(2)
#if defined(half_precision)
    integer(2), allocatable :: tmp(:,:)
#endif /* half_precision */
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif /* SPARC */
#endif /* MPI */

    chbra_original => lab_no2b%GetChannel( ibras(loop) )
    chket_original => lab_no2b%GetChannel( ikets(loop) )

    Jbra = chbra_original%GetJ()
    Pbra = chbra_original%GetParity()
    Tbra = chbra_original%GetT()
    J12bra = chbra_original%GetJ12()
    P12bra = chbra_original%GetP12()
    J3bra  = chbra_original%GetJ3()
    P3bra  = chbra_original%GetP3()

    Jket = chket_original%GetJ()
    Pket = chket_original%GetParity()
    Tket = chket_original%GetT()
    J12ket = chket_original%GetJ12()
    P12ket = chket_original%GetP12()
    J3ket  = chket_original%GetJ3()
    P3ket  = chket_original%GetP3()

    chbra = spno2b%GetIndex(Tbra,J12bra,P12bra,J3bra,P3bra)
    chbra_averaged => spno2b%GetChannel(chbra)
    chket = spno2b%GetIndex(Tket,J12ket,P12ket,J3ket,P3ket)
    chket_averaged => spno2b%GetChannel(chket)

    call oplab_ch%init( chbra_original, chket_original )
    call oplab_ch_no2b%init( chbra_averaged, chket_averaged )
    t1 = omp_get_wtime()
    call tr_bra%init(chbra_original, jac, params%e3max, params%n_threads_tcoef)
    call tr_ket%init(chket_original, jac, params%e3max, params%n_threads_tcoef)
    t2 = omp_get_wtime()

    write(*,"(a,i4,a,i6,a,i6,a,i6,a,i6,a,f8.4,a,f8.4,a,f8.4,a,i6,a,i8)") &
        & "myrank=", myrank, &
        & ", bra: ", ibras(loop), ", ket: ", ikets(loop), &
        & ", bra(avraged): ", chbra, ", ket(averaged): ", chket, &
        & ", Mem(ME): ", oplab_ch_no2b%GetMemory()*1.d3, " MB, Mem(TcBra): ", tr_bra%GetMemory(), &
        & " GB, Mem(TcKet): ", tr_ket%GetMemory(), " GB, ", loop, "/", nloops

    call oplab_ch%trans(tr_bra, tr_ket, opjac)
    t3 = omp_get_wtime()
    call oplab_ch_no2b%AngleAverage( oplab_ch, op_def )
    t4 = omp_get_wtime()
#ifdef MPI
    if(myrank /= 0) then
      chbraket = [chbra,chket]
      n1d = int(chbra_averaged%GetNumberStates(), kind(n1d)) * int(chket_averaged%GetNumberStates(), kind(n1d))
      call mpi_send(chbraket(1),2,mpi_integer,0,tag_channel,mpi_comm_world,ierr)
#if defined(half_precision)
      allocate( tmp( chbra_averaged%GetNumberStates(), chket_averaged%GetNumberStates() ))
      tmp(:,:) = oplab_ch_no2b%m(:,:)%i16
      call mpi_send(tmp(1,1), n1d, mpi_integer2, 0, tag_me_temp, mpi_comm_world, ierr)
      deallocate(tmp)
#elif defined(single_precision)
      call mpi_send(oplab_ch_no2b%m(1,1), n1d, mpi_real4, 0, tag_me_temp, mpi_comm_world, ierr)
#elif defined(double_precision)
      call mpi_send(oplab_ch_no2b%m(1,1), n1d, mpi_real8, 0, tag_me_temp, mpi_comm_world, ierr)
#endif /* precision */
    end if
#endif /* MPI */
    if(myrank == 0) then
      flip = real( (-1.d0)**( (J12ket-J12bra)+(Tket-Tket)/2 ), kind=kind(flip) )
      if( chbra < chket ) then
        opno2b%MatCh(chket,chbra)%m = opno2b%MatCh(chket,chbra)%m + transpose(oplab_ch_no2b%m) * flip
      else
        opno2b%MatCh(chbra,chket)%m = opno2b%MatCh(chbra,chket)%m + oplab_ch_no2b%m
      end if
    end if

    call tr_bra%fin()
    call tr_ket%fin()
    call oplab_ch%fin()
    call oplab_ch_no2b%fin()
  end subroutine trans_to_lab_no2b_channel

  subroutine receive_data()
#ifdef MPI
    integer :: nbra, nket
    real(8) :: t1
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif /* SPARC */
#if defined(single_precision)
    real(4), allocatable :: tmp(:,:)
    real(4) :: flip
#elif defined(half_precision)
    integer(2), allocatable :: tmp(:,:)
    type(my_real16), allocatable :: tmp16(:,:)
    real(4) :: flip
#elif defined(double_precision)
    real(8), allocatable :: tmp(:,:)
    real(8) :: flip
#endif /* precisoin */
    integer :: chbraket(2), ichbra, ichket
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
    call mpi_recv(chbraket,2,mpi_integer,mpi_any_source,tag_channel,mpi_comm_world,istatus,ierr)
    t1 = omp_get_wtime()
    ichbra = chbraket(1)
    ichket = chbraket(2)
    chbra => spno2b%GetChannel(ichbra)
    chket => spno2b%GetChannel(ichket)
    nbra = chbra%GetNumberStates()
    nket = chket%GetNumberStates()
    n1d = int( nbra, kind(n1d)) * int( nket, kind(n1d) )
    allocate(tmp(nbra,nket))
    flip = real( (-1.d0)**( (2*chket%GetJ12()+chket%GetT()-2*chbra%GetJ12()-chbra%GetT())/2 ), kind=kind(flip) )
#if defined(half_precision)
    allocate(tmp16(nbra,nket))
    call mpi_recv(tmp(1,1),n1d,mpi_real4,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
    tmp16(:,:)%i16 = tmp(:,:)
#elif defined(single_precision)
    call mpi_recv(tmp(1,1),n1d,mpi_real4,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#elif defined(double_precisoin)
    call mpi_recv(tmp(1,1),n1d,mpi_real8,istatus(mpi_source),tag_me_temp,mpi_comm_world,istatus,ierr)
#endif /* precision */

#if defined(half_precision)
    if( ichbra < ichket ) then
      if( flip > 0.0 ) then
        opno2b%MatCh(ichket,ichbra)%m = opno2b%MatCh(ichket,ichbra)%m + transpose(tmp16)
      else
        opno2b%MatCh(ichket,ichbra)%m = opno2b%MatCh(ichket,ichbra)%m - transpose(tmp16)
      end if
    else
      opno2b%MatCh(ichbra,ichket)%m = opno2b%MatCh(ichbra,ichket)%m + tmp16
    end if
    deallocate(tmp16)
#elif defined(single_precision) || defined(double_precision)
    if( ichbra < ichket ) then
      opno2b%MatCh(ichket,ichbra)%m = opno2b%MatCh(ichket,ichbra)%m + transpose(tmp)*flip
    else
      opno2b%MatCh(ichbra,ichket)%m = opno2b%MatCh(ichbra,ichket)%m + tmp
    end if
#endif /* precision */
    deallocate(tmp)
    call timer%add(sy%str("Gathering data to parent node"), omp_get_wtime()-t1)
#endif /* MPI */
  end subroutine receive_data
end subroutine MACRO_ADD_SUFX(trans_to_lab_no2b_tensor_)
