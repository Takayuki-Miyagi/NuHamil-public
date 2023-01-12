program main
#ifdef MPI
  use MPIFunction, only: myrank, nprocs, ierr
  use mpi
#else
  use MPIFunction, only: myrank, nprocs
#endif
  use Profiler, only: timer
  use MyLibrary, only: init_dtrinomial, fin_dtrinomial
  use NuHamilInput, only: SetInputParameters, PrintInputParameters, params
  use TwoBodyRelativeSpace
  use OperatorDefinitions, only: SetParameters
  use OneBodyManager
  use TwoBodyManager
  use ThreeBodyManager
  use ABodyManager
  use half_precision_floating_points
  use M1_2b_current, only: init_M1_2b_current, fin_M1_2b_current
  use TwoBodyMultiPoleOperator
  use tests, only: test_main
  implicit none
#ifdef MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
#endif

  if(myrank == 0) then
    write(*,*)
    write(*,'(a)') "############################################################"
    write(*,'(a,10x,a)') "#","NuHamil (Tokyo) code"
    write(*,'(a)') "############################################################"
    write(*,*)
#ifdef VERSION
    write(*,'(2a)') "# Version is ", trim(VERSION)
#endif

#ifdef MPI
    write(*,'(a,i8,a)') "# MPI parallelization: ", nprocs, " process"
#endif
    write(*,'(a,i4)') "# OpenMP number of threads=", omp_get_max_threads()
  end if

  call timer%init()
  call set_table_for_half_precision_floating_point()

  call SetInputParameters()
  call PrintInputParameters()
  call SetFrequency(params%hw) ! pass the frequency to the operator defintion class

  if( nprocs > 1 .and. myrank==0 ) then
    if( params%particle_rank /= 3 .or. .not. params%trans2lab ) then
      write(*,*) "MPI parallelization is not implemented for the selected option"
      write(*,*) "please use with single node"
      stop
    end if
  end if
  call init_dtrinomial()

  call init_M1_2b_current(params%hw, params%N2max, params%Jmax2) ! frequency conversion won't work

  if(params%test_mode) then
    call test_main()
  else

    select case(params%particle_rank)
    case(:0)
      write(*,'(a)') 'Input Error: particle_rank has to be an integer larger than 1.'
    case(1)
      call OneBodyManagers()
    case(2)
      call TwoBodyManagers()
    case(3)
      call ThreeBodyManagers()
    case(4:)
      call ABodyManagers()
    end select

  end if
  call fin_M1_2b_current()

  call fin_dtrinomial()
#ifdef MPI
  call mpi_finalize(ierr)
#endif
  call timer%fin()
end program main
