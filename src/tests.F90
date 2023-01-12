module tests
  use ClassSys
  use MyLibrary
  implicit none
  public :: test_main
  private
  type(sys) :: s
contains

  subroutine test_main()
    use NuHamilInput, only: params
    write(*,*)
    write(*,*) "############################################################"
    write(*,*)
    write(*,*) " Test mode: "
    write(*,*)
    write(*,*) "############################################################"
    select case(params%particle_rank)
    case(:0)
      write(*,'(a)') 'Input Error: particle_rank has to be an integer larger than 1.'
    case(1)
      ! add your test for one-body system
    case(2)
      ! add your test for NN system
    case(3)
      ! add your test for 3N system
    case(4:)
    end select
  end subroutine test_main

end module tests
