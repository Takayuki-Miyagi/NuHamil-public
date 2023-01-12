module NuHamilInput
  use ClassSys
  implicit none

  public :: InputParameters
  public :: SetInputParameters
  public :: PrintInputParameters
  public :: params
  private
  type :: InputParameters
    ! Input parameters
    ! General
    real(8) :: hw ! \hbar\omega value unit of MeV
    type(str) :: renorm ! renormalization method bare/srg/vlowk/ls are available
    real(8) :: lambda ! momentum cutoff scale in unit of fm-1
    integer :: particle_rank ! number of particles

    ! two-body system
    integer :: N2max
    integer :: jmax2
    integer :: J2max_NNint
    integer :: nmesh2
    integer :: N_ls ! N boundary for LeeSuzuki
    integer :: snt_mass
    integer :: J2maxLab
    integer :: Lcm2Max
    real(8) :: pmax2
    type(str) :: renorm_space2, NNInt, srg_generator
    integer :: Nmax_srg_edge
    type(str) :: NNrel_op_input
    type(str) :: NNrel_op_output
    type(str), allocatable :: Operators(:)
    type(str) :: path_to_tmp_dir
    type(str) :: path_to_NNrel
    type(str) :: path_to_Hebeler_files
    type(str) :: Hebeler_fn_head
    type(str), allocatable :: Hebeler_fn_tails(:)
    type(str) :: fname_jacobi_op_for_trans2lab
    type(str) :: input_nn_file
    logical :: coul
    logical :: trans2lab
    logical :: pn_same_mass
    type(str) :: file_nn_cartesian
    type(str) :: file_name_n
    type(str) :: file_name_nn
    type(str) :: file_name_3n
    type(str) :: averaged_file_for_test
    logical :: evolve_coulomb
    integer :: spin_tensor_decomposition
    real(8) :: LECs_DeltaGO394(17) ! Used if NNint='DeltaGO394_***'
    character(7) :: LECs_names_DeltaGO394(17) = &
        &["Ct1S0pp", "Ct1S0np", "Ct1S0nn", "Ct3S1  ", "C1S0   ", &
        & "C3P0   ", "C1P1   ", "C3P1   ", "C3S1   ", "CE1    ", &
        & "C3P2   ", "c1     ", "c2     ", "c3     ", "c4     ", "cD     ", "cE     " ]


    ! three-body system
    integer :: N3max
    integer :: J12max ! truncation for genuine three-body force
    logical :: genuine_3bf
    real(8) :: c1, c3, c4 ! chiral N2LO 3BF two-pion exchange
    real(8) :: cd         ! chiral N2LO 3BF one-pion exchange
    real(8) :: ce         ! chiral N2LO contact term
    real(8) :: lambda_3nf_local    ! chiral N2LO local regulator cutoff
    real(8) :: lambda_3nf_nonlocal ! chiral N2LO non-local regulator cutoff
    logical :: save_3nf_before_lec
    type(str) :: NNNInt, ramp ! three-body-model space
    type(str) :: Regulator ! Local/NonLocal
    integer :: RegulatorPower
    integer :: Ngrid_3nf
    integer :: jmax3
    integer :: j3max_initial_3nf
    logical :: hw_conversion
    logical :: pn_form
    real(8) :: hw_target
    integer :: nblock
    integer :: n_threads_tcoef
    integer :: svd_rank_tcoef
    integer :: svd_rank_operator
    integer :: svd_rank_op_lab
    type(str) :: lab_3bme_precision

    integer :: NAmax
    integer :: jmaxA, tmaxA
    logical :: NN_only
    integer :: emax, e2max, e3max, lmax
    integer :: jbra, pbra, tbra, tzbra
    integer :: jket, pket, tket, tzket
    logical :: only_hf_monopole
    logical :: only_no2b_elements
    type(str) :: fn_fewbody_wave_function

    ! file format convert mode
    logical :: file_convert
    integer :: emax_convert, e2max_convert, e3max_convert, lmax_convert
    type(str) :: file_name_nn_original
    type(str) :: file_name_3n_original
    type(str) :: file_name_nn_converted
    type(str) :: file_name_3n_converted

    ! file combine mode
    type(str), allocatable :: files_combined(:)
    real(8), allocatable :: weights_combined(:)

    logical :: test_mode = .false.
    logical :: count_memory=.false.
    integer :: no2b_channel_begin
    integer :: no2b_channel_end
    type(str) :: no2b_temp_dir

    type(str) :: file_name_phase_shift
    real(8), allocatable :: Tlab_phase_shift(:)
  contains
    procedure :: PrintParametersALL
    procedure :: CopyInputParameters
    generic :: assignment(=) => CopyInputParameters
  end type InputParameters

  type(InputParameters), protected :: params
contains
  subroutine SetInputParameters()
    use omp_lib
    type(sys) :: s
    ! common
    real(8) :: hw = 20.d0
    character(256) :: renorm = 'bare'
    real(8) :: lambda = 2.d0
    integer :: rank = 2
    integer :: Nmax_srg_edge=-1
    ! 2-body
    integer :: N2max = 100
    integer :: N_ls = 0
    integer :: jmax2 = 12
    integer :: J2max_NNint=8
    integer :: J2maxLab=-1
    integer :: Lcm2Max=-1
    integer :: nmesh2 =100
    real(8) :: pmax2 = 8.d0
    integer :: snt_mass = 4
    integer :: svd_rank_tcoef = -1
    integer :: svd_rank_operator = -1
    integer :: svd_rank_op_lab = -1
    character(256) :: NNint = 'N3LO_EM500'
    character(256) :: srg_generator = "kinetic"
    character(512) :: file_nn_cartesian = "none"
    character(512) :: file_name_n = "default"
    character(512) :: file_name_nn = "default"
    character(512) :: file_name_3n = "default"
    character(256) :: renorm_space2 = 'ho'
    character(2048) :: Operators = 'hamil'
    character(256) :: path_to_tmp_dir = "./"
    character(256) :: path_to_NNrel = "default"
    character(256) :: path_to_hebeler_files = "./"
    character(256) :: Hebeler_fn_head = "3NF_V"
    character(1024) :: Hebeler_fn_tails = "N2LO_c1,N2LO_c3,N2LO_c4,N2LO_cD,N2LO_cE"
    logical :: trans2lab =.true.
    logical :: coul = .true.
    logical :: is_PhaseShift = .false.
    logical :: evolve_coulomb = .true.
    logical :: pn_same_mass = .true. ! mom -> HO, assuming the same pn mass; same for SRG scale
    ! 3-body
    integer :: N3max = 40
    integer :: j3max_initial_3nf=-1
    integer :: jmax3 = 51
    integer :: J12max = -1 ! truncation for genuine three-body force
    real(8) :: hw_target = -1.d0
    real(8) :: LECs(5) = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
    real(8) :: LECs_DeltaGO394(17) = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
    real(8) :: lambda_3nf_local = 400.d0
    real(8) :: lambda_3nf_nonlocal = 500.d0
    character(256) :: Regulator = 'Local'
    character(256) :: NNNint = 'ChEFT_N2LO'
    logical :: save_3nf_before_lec=.false.
    integer :: RegulatorPower = 2
    integer :: Ngrid_3nf = 120
    character(256) :: ramp = 'flat40'
    logical :: genuine_3bf = .true.
    logical :: pn_form = .false.
    logical :: only_hf_monopole=.false.
    logical :: only_no2b_elements=.false.
    integer :: nblock = 400
    integer :: n_threads_tcoef = -1
    character(6) :: lab_3bme_precision="single"
    character(512) :: fname_jacobi_op_for_trans2lab="none"
    !
    integer :: emax = 6
    integer :: e2max=12
    integer :: e3max=6
    integer :: lmax = -1
    ! A-body
    integer :: bra(4) = [-1,-100,-1,-100]
    integer :: ket(4) = [-1,-100,-1,-100]
    integer :: NAmax=10
    logical :: NN_only=.false.
    character(512) :: fn_fewbody_wave_function = "none"

    ! file format convert mode
    logical :: file_convert = .false.
    character(512) :: file_name_nn_original = "none"
    character(512) :: file_name_3n_original = "none"
    character(512) :: file_name_nn_converted = "default"
    character(512) :: file_name_3n_converted = "default"
    character(512) :: NNrel_op_input = "none"
    character(512) :: NNrel_op_output = "none"
    character(512) :: input_nn_file = 'none'
    integer :: emax_convert=-1, e2max_convert=-1, e3max_convert=-1, lmax_convert=-1
    character(2048) :: files_combined = ""
    character(1024) :: weights_combined = "none"
    type(str) :: Operators_str, files_combined_str, del

    logical :: ex
    logical :: test_mode = .false.
    integer :: no2b_channel_begin = -1
    integer :: no2b_channel_end = -1
    character(256) :: no2b_temp_dir = "NO2B_temp"
    character(256) :: averaged_file_for_test="none"
    integer :: spin_tensor_decomposition = -1
    type(str) :: weights_combined_str
    type(str), allocatable :: weights_combined_strs(:)

    character(256) :: file_name_phase_shift = "none"
    character(256) :: Tlab_phase_shift = "1,5,10,25,50,100,150,200,250,300,350"
    type(str) :: Tlab_phase_shift_str
    type(str), allocatable :: Tlab_phase_shift_strs(:)

    ! input file name
    character(256) :: inputfile
    integer :: i

    logical :: count_memory=.false.
    namelist /input/ rank, hw, renorm, lambda, N2max, &
        & jmax2, nmesh2, pmax2, NNInt, trans2lab, &
        & N3max, jmax3, LECs, path_to_tmp_dir, &
        & lambda_3nf_local, lambda_3nf_nonlocal, snt_mass, &
        & genuine_3bf, NAmax, &
        & renorm_space2, hw_target, &
        & NN_only, Operators, coul, &
        & ramp, bra, ket, &
        & pn_form, emax, e2max, e3max, &
        & J12max, Regulator, NNNInt, RegulatorPower, Ngrid_3nf, &
        & nblock, is_PhaseShift, &
        & file_name_nn_original, file_name_3n_original, &
        & file_name_nn_converted, file_name_3n_converted, &
        & file_convert, lmax, &
        & test_mode, save_3nf_before_lec, &
        & only_hf_monopole, only_no2b_elements, count_memory, &
        & file_name_nn, file_name_3n, j3max_initial_3nf, averaged_file_for_test, &
        & srg_generator, J2max_NNint, &
        & no2b_channel_begin, no2b_channel_end, no2b_temp_dir, &
        & emax_convert, e2max_convert, e3max_convert, lmax_convert, &
        & NNrel_op_input, NNrel_op_output, n_threads_tcoef, file_nn_cartesian, evolve_coulomb, Nmax_srg_edge, &
        & spin_tensor_decomposition, LECs_DeltaGO394, files_combined, weights_combined, &
        & lab_3bme_precision, path_to_NNrel, path_to_hebeler_files, fname_jacobi_op_for_trans2lab, file_name_n, &
        & svd_rank_tcoef, svd_rank_operator, svd_rank_op_lab, &
        & file_name_phase_shift, Tlab_phase_shift, J2maxLab, Lcm2Max, pn_same_mass, fn_fewbody_wave_function, &
        & input_nn_file, Hebeler_fn_head, Hebeler_fn_tails


    call getarg(1, inputfile)
    inquire(file = inputfile, exist = ex)
    if(.not.ex) then
      write(*,'(3a)') 'Error: Input file: ', trim(inputfile), ' was not found.'
      stop
    end if
    open(118, file = inputfile, status = 'old')
    read(118, nml=input)
    close(118)
    params%particle_rank = rank
    params%hw = hw
    params%renorm = renorm
    params%lambda = lambda
    params%N_ls = N_ls
    params%pn_same_mass = pn_same_mass

    params%N2max = n2max
    params%jmax2 = jmax2
    params%J2maxLab = J2maxLab
    params%Lcm2Max = Lcm2Max
    params%j2max_NNint = j2max_NNint
    params%nmesh2 = nmesh2
    params%pmax2 = pmax2
    params%NNInt = NNInt
    params%Nmax_srg_edge = Nmax_srg_edge
    params%srg_generator = srg_generator
    params%file_name_n = file_name_n
    params%file_name_nn = file_name_nn
    params%file_name_3n = file_name_3n
    params%renorm_space2 = renorm_space2
    params%path_to_tmp_dir = path_to_tmp_dir
    params%path_to_hebeler_files = path_to_hebeler_files
    params%path_to_NNrel = path_to_NNrel
    params%coul = coul
    params%trans2lab = trans2lab
    params%snt_mass = snt_mass
    params%NNrel_op_input = NNrel_op_input
    params%NNrel_op_output = NNrel_op_output
    params%evolve_coulomb = evolve_coulomb
    params%spin_tensor_decomposition = spin_tensor_decomposition
    params%LECs_DeltaGO394 = LECs_DeltaGO394

    params%N3max = N3max
    params%J12max = J12max
    params%jmax3 = jmax3
    params%j3max_initial_3nf = j3max_initial_3nf
    params%genuine_3bf = genuine_3bf
    params%hw_target = hw_target
    params%c1 = LECs(1)
    params%c3 = LECs(2)
    params%c4 = LECs(3)
    params%cd = LECs(4)
    params%ce = LECs(5)
    params%Ngrid_3nf = Ngrid_3nf
    params%lambda_3nf_local = lambda_3nf_local
    params%lambda_3nf_nonlocal = lambda_3nf_nonlocal
    params%ramp = ramp
    params%Regulator = Regulator
    params%NNNInt = NNNInt
    params%RegulatorPower = RegulatorPower
    params%nblock = nblock
    params%n_threads_tcoef = n_threads_tcoef
    params%save_3nf_before_lec= save_3nf_before_lec
    params%lab_3bme_precision = lab_3bme_precision
    params%fname_jacobi_op_for_trans2lab = fname_jacobi_op_for_trans2lab

    params%NAmax = NAmax
    params%NN_only = NN_only
    params%jbra = bra(1)
    params%pbra = bra(2)
    params%tbra = bra(3)
    params%tzbra = bra(4)
    params%jket = ket(1)
    params%pket = ket(2)
    params%tket = ket(3)
    params%tzket = ket(4)
    params%pn_form = pn_form
    params%emax = emax
    params%lmax = lmax
    params%e2max = e2max
    params%e3max = e3max

    params%file_convert = file_convert
    params%file_name_nn_original = file_name_nn_original
    params%file_name_3n_original = file_name_3n_original
    params%file_name_nn_converted = file_name_nn_converted
    params%file_name_3n_converted = file_name_3n_converted
    params%file_name_phase_shift = file_name_phase_shift
    params%emax_convert = emax_convert
    params%e2max_convert = e2max_convert
    params%e3max_convert = e3max_convert
    params%lmax_convert = lmax_convert
    params%only_hf_monopole = only_hf_monopole
    params%only_no2b_elements = only_no2b_elements
    params%averaged_file_for_test = averaged_file_for_test
    params%test_mode = test_mode
    params%no2b_channel_begin = no2b_channel_begin
    params%no2b_channel_end = no2b_channel_end
    params%no2b_temp_dir = no2b_temp_dir
    params%svd_rank_tcoef = svd_rank_tcoef
    params%svd_rank_operator = svd_rank_operator
    params%svd_rank_op_lab = svd_rank_op_lab
    params%fn_fewbody_wave_function = fn_fewbody_wave_function
    params%input_nn_file = input_nn_file

    params%file_nn_cartesian = file_nn_cartesian
    if(params%lmax==-1) params%lmax = params%emax
    if(params%file_name_nn_converted%val=="default") params%file_name_nn_converted = params%file_name_nn
    if(params%file_name_3n_converted%val=="default") params%file_name_3n_converted = params%file_name_3n
    if(params%path_to_NNrel%val=="default") params%path_to_NNrel = params%path_to_tmp_dir
    if(params%j3max_initial_3nf==-1) params%j3max_initial_3nf = params%jmax3
    if(params%emax_convert == -1) params%emax_convert = params%emax
    if(params%e2max_convert == -1) params%e2max_convert = params%e2max
    if(params%e3max_convert == -1) params%e3max_convert = params%e3max
    if(params%lmax_convert == -1) params%lmax_convert = params%emax_convert
    if(params%n_threads_tcoef == -1)  params%n_threads_tcoef = omp_get_max_threads()

    params%count_memory = count_memory
    if(params%hw_target < 0.d0) then
      params%hw_target=params%hw
      params%hw_conversion = .false.
    else
      params%hw_conversion = .true.
    end if

    params%Hebeler_fn_head = Hebeler_fn_head
    Operators_str = Hebeler_fn_tails
    del = ","
    call s%split(Operators_str,del,params%Hebeler_fn_tails)
    Operators_str = Operators
    files_combined_str = files_combined
    del = ","
    call s%split(Operators_str,del,params%Operators)
    call s%split(files_combined_str,del,params%files_combined)
    Tlab_phase_shift_str = Tlab_phase_shift
    call s%split(Tlab_phase_shift_str, del, Tlab_phase_shift_strs)
    allocate(params%Tlab_phase_shift(size(Tlab_phase_shift_strs)))
    do i = 1, size(Tlab_phase_shift_strs)
      read(Tlab_phase_shift_strs(i)%val,*) params%Tlab_phase_shift(i)
    end do

    if( weights_combined/="none" ) then
      weights_combined_str = weights_combined
      call s%split(weights_combined_str,del,weights_combined_strs)
      allocate(params%weights_combined(size(weights_combined_strs)))
      do i = 1, size(weights_combined_strs)
        read(weights_combined_strs(i)%val,*) params%weights_combined(i)
      end do
    end if
    ex = .false.
#if defined(__GFORTRAN__)
    inquire(file=trim(params%path_to_Hebeler_files%val)//"/.", exist=ex)
#endif
#if defined(__INTEL_COMPILER)
    inquire(directory=params%path_to_Hebeler_files%val, exist=ex)
#endif
    if(.not. ex) then
      write(*,*)
      write(*,"(3a)") "Error, directory ", trim(params%path_to_Hebeler_files%val), " does not exists!"
      write(*,*)
      stop
    end if
  end subroutine SetInputParameters

  subroutine CopyInputParameters(params2, params1)
    class(InputParameters), intent(inout) :: params2
    type(InputParameters), intent(in) :: params1
    integer :: i

    params2%hw =            params1%hw
    params2%renorm =        params1%renorm
    params2%lambda =        params1%lambda
    params2%particle_rank = params1%particle_rank
    params2%N2max =         params1%N2max
    params2%jmax2 =         params1%jmax2
    params2%j2max_NNint = params1%j2max_NNint
    params2%J2maxLab = params1%J2maxLab
    params2%Lcm2Max = params1%Lcm2Max
    params2%nmesh2 =        params1%nmesh2
    params2%pmax2  =        params1%pmax2
    params2%NNInt  =        params1%NNInt
    params2%snt_mass =      params1%snt_mass
    params2%file_name_n =       params1%file_name_n
    params2%file_name_nn=       params1%file_name_nn
    params2%file_name_3n=       params1%file_name_3n
    params2%N_ls   =       params1%N_ls
    params2%Nmax_srg_edge = params1%Nmax_srg_edge
    params2%renorm_space2 = params1%renorm_space2
    params2%coul         =  params1%coul
    params2%count_memory         =  params1%count_memory
    params2%jbra            = params1%jbra
    params2%pbra            = params1%pbra
    params2%tbra            = params1%tbra
    params2%tzbra           = params1%tzbra
    params2%jket            = params1%jket
    params2%pket            = params1%pket
    params2%tket            = params1%tket
    params2%tzket           = params1%tzket
    params2%trans2lab    = params1%trans2lab
    params2%N3max =         params1%N3max
    params2%J12max =         params1%J12max
    params2%jmax3 =         params1%jmax3
    params2%j3max_initial_3nf = params1%j3max_initial_3nf
    params2%genuine_3bf =   params1%genuine_3bf
    params2%pn_same_mass = params1%pn_same_mass
    params2%hw_conversion =params1%hw_conversion
    params2%hw_target =    params1%hw_target
    params2%c1 =            params1%c1
    params2%c3 =            params1%c3
    params2%c4 =            params1%c4
    params2%cd =            params1%cd
    params2%ce =            params1%ce
    params2%lambda_3nf_local =  params1%lambda_3nf_local
    params2%lambda_3nf_nonlocal =  params1%lambda_3nf_nonlocal
    params2%Ngrid_3nf =  params1%Ngrid_3nf
    params2%ramp =          params1%ramp
    params2%Regulator =      params1%Regulator
    params2%NNNInt    =      params1%NNNInt
    params2%srg_generator= params1%srg_generator
    params2%RegulatorPower = params1%RegulatorPower
    params2%nblock = params1%nblock
    params2%NAmax =         params1%NAmax
    params2%jmaxA =         params1%jmaxA
    params2%tmaxA =         params1%tmaxA
    params2%NN_only =     params1%NN_only
    params2%pn_form =       params1%pn_form
    params2%emax =          params1%emax
    params2%lmax =          params1%lmax
    params2%e2max =         params1%e2max
    params2%e3max =         params1%e3max
    params2%save_3nf_before_lec = params1%save_3nf_before_lec
    params2%test_mode = params1%test_mode
    params2%no2b_channel_begin = params1%no2b_channel_begin
    params2%no2b_channel_end = params1%no2b_channel_end
    params2%no2b_temp_dir = params1%no2b_temp_dir
    params2%evolve_coulomb = params1%evolve_coulomb
    params2%spin_tensor_decomposition = params1%spin_tensor_decomposition

    params2%file_convert = params1%file_convert
    params2%file_name_nn_original = params1%file_name_nn_original
    params2%file_name_3n_original = params1%file_name_3n_original
    params2%file_name_nn_converted = params1%file_name_nn_converted
    params2%file_name_3n_converted = params1%file_name_3n_converted
    params2%emax_convert = params1%emax_convert
    params2%e2max_convert = params2%e2max_convert
    params2%e3max_convert = params2%e3max_convert
    params2%lmax_convert = params2%lmax_convert
    params2%averaged_file_for_test = params1%averaged_file_for_test
    params2%NNrel_op_input = params1%NNrel_op_input
    params2%NNrel_op_output = params1%NNrel_op_output

    params2%only_hf_monopole  =         params1%only_hf_monopole
    params2%only_no2b_elements=         params1%only_no2b_elements
    params2%path_to_tmp_dir = params1%path_to_tmp_dir
    params2%path_to_NNrel = params1%path_to_NNrel
    params2%path_to_hebeler_files = params1%path_to_hebeler_files
    params2%file_nn_cartesian = params1%file_nn_cartesian
    params2%LECs_DeltaGO394 = params1%LECs_DeltaGO394
    params2%lab_3bme_precision = params1%lab_3bme_precision
    params2%fname_jacobi_op_for_trans2lab = params1%fname_jacobi_op_for_trans2lab
    params2%n_threads_tcoef = params1%n_threads_tcoef
    params2%count_memory = params1%count_memory
    params2%svd_rank_tcoef = params1%svd_rank_tcoef
    params2%svd_rank_operator = params1%svd_rank_operator
    params2%svd_rank_op_lab = params1%svd_rank_op_lab
    params2%file_name_phase_shift = params1%file_name_phase_shift
    params2%fn_fewbody_wave_function = params1%fn_fewbody_wave_function
    params2%input_nn_file = params1%input_nn_file
    params2%Hebeler_fn_head = params1%Hebeler_fn_head
    if(allocated(params2%Hebeler_fn_tails)) deallocate(params2%Hebeler_fn_tails)
    allocate(params2%Hebeler_fn_tails(size(params1%Hebeler_fn_tails)))
    do i = 1, size(params1%Hebeler_fn_tails)
      params2%Hebeler_fn_tails(i) = params1%Hebeler_fn_tails(i)
    end do

    if(allocated(params2%Operators)) deallocate(params2%Operators)
    allocate(params2%Operators(size(params1%Operators)))
    do i = 1, size(params1%Operators)
      params2%Operators(i) = params1%Operators(i)
    end do

    if(allocated(params1%weights_combined)) then
      if(allocated(params2%weights_combined)) deallocate(params2%weights_combined)
      allocate(params2%weights_combined(size(params1%weights_combined)))
      do i = 1, size(params1%weights_combined)
        params2%weights_combined(i) = params1%weights_combined(i)
      end do
    end if

    if(allocated(params1%files_combined)) then
      if(allocated(params2%files_combined)) deallocate(params2%files_combined)
      allocate(params2%files_combined(size(params1%files_combined)))
      do i = 1, size(params1%files_combined)
        params2%files_combined(i) = params1%files_combined(i)
      end do
    end if

    if(allocated(params2%Tlab_phase_shift))  deallocate(params2%Tlab_phase_shift)
    allocate(params2%Tlab_phase_shift(size(params1%Tlab_phase_shift)))
    do i = 1, size(params1%Tlab_phase_shift)
      params2%Tlab_phase_shift(i) = params1%Tlab_phase_shift(i)
    end do
  end subroutine CopyInputParameters

  subroutine PrintInputParameters(unt)
    use MPIFunction, only: myrank
    integer, intent(in), optional :: unt
    integer :: iunit, i
    type(sys) :: s
    if(present(unt)) then
      iunit = unt
    else
      iunit = 6
    end if
    if(myrank /= 0) return
    write(iunit,*)
    write(iunit,'(a)') '################### Input Parameters #####################'

    if(params%file_convert) then
      write(iunit,'(a)') "# File format convert mode:"
      if(params%particle_rank == 2) then
        write(iunit,"(2a)") "# Input file:  ", trim(params%file_name_nn_original%val)
        write(iunit,"(2a)") "# Output file: ", trim(params%file_name_nn_converted%val)
        return
      end if

      if(params%particle_rank == 3) then
        write(iunit,"(2a)") "# Input file:  ", trim(params%file_name_3n_original%val)
        write(iunit,"(2a)") "# Output file: ", trim(params%file_name_3n_converted%val)
        return
      end if
    end if

    if(allocated(params%weights_combined) .and. size(params%files_combined)==size(params%weights_combined)) then
        write(iunit,"(a)") "# File combine mode:  weight factor,   filename"
        do i = 1, size(params%files_combined)
          write(*,'(a,f12.8,2a)') "# ", params%weights_combined(i), ", ", trim(params%files_combined(i)%val)
        end do
        write(iunit,"(2a)") "# Output file name: ", trim(params%file_name_3n%val)
        return
    end if

    write(iunit,'(a, f6.2, a)') '# hw = ', params%hw, ' MeV'
    if(params%renorm%val == 'bare') then
      write(iunit,'(2a)') '# renorm: ', trim(params%renorm%val)
    else
      write(iunit,'(3a, f6.2)') '# renorm: ', trim(params%renorm%val), ', lambda = ', params%lambda
      if(params%renorm%val == "srg") write(iunit,'(2a)') '# Generator of SRG: ', trim(params%srg_generator%val)
    end if
    write(iunit,"(a)") "# "
    write(iunit, '(a)') '# Parameters for two-body system'
    write(iunit,'(a, i3, a, i3, a, i3)') '# N2max = ', params%N2max, ',  Jmax = ', params%jmax2, &
        &", Jmax(NNint) = ", params%J2max_NNint
    write(iunit, '(2a)') '# NN interaction: ', trim(params%NNInt%val)
    write(iunit,'(a, f6.2, a, i4)') '# Maximum momentum = ', params%pmax2, &
        &' fm-1,  Number of gauss-legendre mesh points', params%nmesh2
    write(iunit,'(a,a)') '# input_nn_file = ', params%input_nn_file%val
    if(params%particle_rank > 2) then
      write(iunit,'(a, i3, a, i3)') '# N3max = ', params%N3max, ',  Jmax = ', params%jmax3
      if(params%genuine_3bf) then
        write(iunit,'(a)') '# Genuine three-body force will be included'
        write(iunit,'(3a,i2,a,f8.4,a,f8.4,a)') '# Regulator: ', trim(params%Regulator%val), &
            & ', Regulator Power: ', params%RegulatorPower, &
            & ", Cutoff Lambda (Local): ", params%lambda_3nf_local, &
            & " MeV, Cutoff Lambda (NonLocal): ", params%lambda_3nf_nonlocal, " MeV"
        write(iunit,'(a, f9.4, a, f9.4, a, f9.4, a)') '# Two-pion exchange: c1 = ', &
            & params%c1, ' GeV-1, c3 = ', params%c3, ' GeV-1, c4 = ', &
            & params%c4, ' GeV-1'
        write(iunit,'(a, f9.4)') '# One-pion exchange: cD = ', &
            & params%cd
        write(iunit,'(a, f9.4)') '# Contact term: cE = ', &
            & params%cE
      end if
    end if

    if(params%particle_rank > 3) then
    end if

    write(iunit,"(a)") "# "
    write(iunit,'(a)') '################### For Output Files #####################'
    if(params%particle_rank == 2) then
      if(params%trans2lab) then
        write(iunit,'(a, i3, a, i3)') '# emax = ', params%emax, ', e2max = ', params%e2max
        write(iunit,'(2a)') "# file name: ", trim(params%file_name_nn%val)
      else
        write(iunit,'(a)') '# Parameters for 2-body calculations'
        write(iunit,'(a, i3, a, i3, a, i3, a, i3)') '# bra : N2max = ', params%N2max, &
            & ',  J = ', params%jbra, ',  Pari = ', params%pbra, ',  Tz = ', params%tzbra
        write(iunit,'(a, i3, a, i3, a, i3, a, i3)') '# ket : N2max = ', params%N2max, &
            & ',  J = ', params%jket, ',  Pari = ', params%pket, ',  Tz = ', params%tzket
      end if
    elseif(params%particle_rank == 3) then
      if(params%trans2lab) then
        write(iunit,'(a, i3, a, i3, a, i3)') '# emax = ', params%emax, ', e2max = ', params%e2max, ", e3max = ", params%e3max
        write(iunit,'(2a)') "# file name: ", trim(params%file_name_3n%val)
        write(iunit,'(3a)') "# ", trim(params%lab_3bme_precision%val), " precision"
      else
        write(iunit,'(a)') '# Parameters for 3-body calculations'
        write(iunit,'(a, i3, a, i3, a, i3, a, i3, a, i3)') '# bra : N3max = ', params%N3max, &
            & ',  J = ', params%jbra, ',  Pari = ', params%pbra, ',  T = ', params%tbra, &
            & ', Tz = ', params%tzbra
        write(iunit,'(a, i3, a, i3, a, i3, a, i3, a, i3)') '# ket : N3max = ', params%N3max, &
            & ',  J = ', params%jket, ',  Pari = ', params%pket, ',  T = ', params%tket, &
            & ', Tz = ', params%tzket
      end if
    elseif(params%particle_rank >= 4) then
      write(iunit,'(a,i1,a)') '# Parameters for ', params%particle_rank, '-body calculations'
      write(iunit,'(a, i3, a, i3, a, i3, a, i3, a, i3)') '# bra : NAmax = ', params%NAmax, &
          & ',  J = ', params%jbra, ',  Pari = ', params%pbra, ',  T = ', params%tbra, &
          & ', Tz = ', params%tzbra
      write(iunit,'(a, i3, a, i3, a, i3, a, i3, a, i3)') '# ket : NAmax = ', params%NAmax, &
          & ',  J = ', params%jket, ',  Pari = ', params%pket, ',  T = ', params%tket, &
          & ', Tz = ', params%tzket
    end if
    write(iunit,'(a)') "# Calculations will be done for Operators:"
    write(iunit,'(a)', advance='no') "# "
    do i = 1, size(params%Operators)
      write(iunit,'(2a)', advance='no') trim(params%Operators(i)%val), ', '
    end do
    write(iunit,*)
    write(iunit,'(2a)') "# directory for temporary files: ", trim(params%path_to_tmp_dir%val)
  end subroutine PrintInputParameters

  subroutine PrintParametersALL(this, iunit)
    ! This method is ugly and is used only for debuging
    class(InputParameters), intent(in) :: this
    integer, intent(in), optional :: iunit
    integer :: i, unt = 6
    if(present(iunit)) unt = iunit

    write(unt,'(a,i2)') '# particle_rank = ', this%particle_rank
    write(unt,'(a,f6.2)') '# hw = ', this%hw
    write(unt,'(a,a)') '# renorm = ', trim(this%renorm%val)
    write(unt,'(a,f8.4)') '# lambda = ', this%lambda
    write(unt,'(a,i4)') '# N2max = ', this%N2max
    write(unt,'(a,i3)') '# jmax2 = ', this%jmax2
    write(unt,'(a,i4)') '# nmesh2 = ', this%nmesh2
    write(unt,'(a,f6.2)') '# pmax2 = ', this%pmax2
    write(unt,'(a,a)') '# NNInt = ', trim(this%NNInt%val)
    write(unt,'(a,a)') '# renorm_space2 = ', trim(this%renorm_space2%val)
    write(unt,'(a,a)') '# path_to_tmp_dir = ', trim(this%path_to_tmp_dir%val)
    write(unt,'(a,a)') '# path_to_NNrel = ', trim(this%path_to_NNrel%val)
    write(unt,'(a,L2)') '# trans2lab = ', this%trans2lab
    write(unt,'(a,L2)') '# coul = ', this%coul
    if(this%particle_rank /= 2) then
      write(unt,'(a,i3)') '# N3max = ', this%N3max
      write(unt,'(a,i3)') '# jmax3 = ', this%jmax3
      write(unt,'(a,L2)') '# hw_conversion = ', this%hw_conversion
      if(this%hw_conversion) write(unt,'(a,f6.2)') '# hw_target = ', this%hw_target
      write(unt,'(a,L2)') '# genuine_3bf = ', this%genuine_3bf
      if(this%genuine_3bf) then
        write(unt,'(a,f8.3)') '# c1 = ', this%c1
        write(unt,'(a,f8.3)') '# c3 = ', this%c3
        write(unt,'(a,f8.3)') '# c4 = ', this%c4
        write(unt,'(a,f8.3)') '# cd = ', this%cd
        write(unt,'(a,f8.3)') '# ce = ', this%ce
        write(unt,'(a,f8.2)') '# lambda_3nf_local = ', this%lambda_3nf_local
        write(unt,'(a,f8.2)') '# lambda_3nf_nonlocal = ', this%lambda_3nf_nonlocal
        write(unt,'(a,a)') '# Regulator = ', trim(this%Regulator%val)
        write(unt,'(a,a)') '# NNNInt = ', trim(this%NNNInt%val)
        write(unt,'(a,i2)') '# RegulatorPower = ', this%RegulatorPower
        write(unt,'(a,i3)') '# Ngrid_3nf = ', this%Ngrid_3nf
        write(unt,'(a,i3)') '# J12max = ', this%J12max
      end if
      write(unt,'(a,i4)') '# nblock = ', this%nblock
      write(unt,'(a,a)') '# ramp = ', trim(this%ramp%val)
      write(unt,'(a,L2)') '# pn_form = ', this%pn_form
      write(unt,'(a,L2)') '# NN_only = ', this%NN_only
    end if

    if(this%trans2lab) then
      write(unt,'(a,i3)') '# emax = ', this%emax
      write(unt,'(a,i3)') '# e2max = ', this%e2max
      if(this%particle_rank ==3) write(unt,'(a,i3)') '# e3max = ', this%e3max
    end if

    if(.not. this%trans2lab) then
      write(unt,'(a,i3)') '# jbra = ', this%jbra
      write(unt,'(a,i3)') '# pbra = ', this%pbra
      write(unt,'(a,i3)') '# tbra = ', this%tbra
      write(unt,'(a,i3)') '# tzbra = ', this%tzbra
      write(unt,'(a,i3)') '# jket = ', this%jket
      write(unt,'(a,i3)') '# pket = ', this%pket
      write(unt,'(a,i3)') '# tket = ', this%tket
      write(unt,'(a,i3)') '# tzket = ', this%tzket
    end if

    if(this%particle_rank > 3) write(unt,'(a,i3)') '# NAmax = ', this%NAmax
    do i = 1, size(this%Operators)
      write(unt,'(2a)') '# Operator: ', trim(this%Operators(i)%val)
    end do
  end subroutine PrintParametersALL

end module NuHamilInput

