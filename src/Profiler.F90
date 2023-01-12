module Profiler
  use omp_lib
  use ClassSys
#ifdef MPI
  use MPIFunction, only: myrank
#endif
  implicit none

  public :: timer
  private
  type, private :: prof
    type(imap) :: counter
    type(dmap) :: timer
    real(8) :: start_time
  contains
    procedure :: init => InitProf
    procedure :: start => StartProf
    procedure :: fin => FinProf
    procedure :: add => AddProf
    procedure :: prt => PrintSummary
    procedure :: PrintMemoryUsage
  end type prof

  type(prof) :: timer
  type(sys) :: s

contains
  subroutine InitProf(this)
    class(prof), intent(inout) :: this

#ifdef MPI
    if(myrank /= 0) return
#endif

    this%start_time = omp_get_wtime()
    call this%counter%append(s%str('Total'),1)
    call this%timer%append(s%str('Total'),0.d0)
  end subroutine InitProf

  subroutine StartProf(this,key)
    class(prof), intent(inout) :: this
    type(str), intent(in) :: key
#ifdef MPI
    if(myrank /= 0) return
#endif
    call this%counter%append(key,0)
    call this%timer%append(key,0.d0)
  end subroutine StartProf

  subroutine AddProf(this,key,time_in)
    class(prof), intent(inout) :: this
    type(str), intent(in) :: key
    real(8), intent(in) :: time_in
    integer :: n
    real(8) :: time
#ifdef MPI
    if(myrank /= 0) return
#endif
    if(.not. this%counter%Search(key)) call this%start(key)
    n = this%counter%Get(key)
    time = this%timer%Get(key)
    call this%counter%append(key,n+1)
    call this%timer%append(key,time+time_in)
  end subroutine AddProf

  subroutine PrintMemoryUsage(this)
    ! only for linux
    class(prof), intent(inout) :: this
    integer :: getpid
    character(512) :: c_num, fn_proc, cmd

#ifdef MPI
    if(myrank /= 0) return
#endif
    write(c_num,'(i0)') getpid()
    fn_proc = '/proc/' // trim(c_num) // '/status'
    cmd = 'cat ' // trim(fn_proc) // ' | grep VmRSS'
    call system( trim(cmd) )
  end subroutine PrintMemoryUsage

  subroutine PrintSummary(this)
    use MPIFunction, only: myrank
    class(prof), intent(inout) :: this
    integer :: n, i
    real(8) :: r, ttotal
#ifdef MPI
    if(myrank /= 0) return
#endif
    call this%timer%append(s%str('Total'),omp_get_wtime() - this%start_time)
    if (myrank==0) write(*,*)
    ttotal = this%timer%get(s%str('Total'))
    n = int(ttotal)
    if (myrank==0) write(*,'(a,i5,a,i2.2,a,i2.2)') &
        '    summary of time, total = ', n/3600, ':', &
        mod(n, 3600)/60, ':', mod(n, 60)
    if (myrank==0) write(*,*)
    if (myrank==0) write(*,'(37x,a)') "time,    ncall, time/ncall,   ratio "
    r = ttotal
    do i = 1, this%counter%GetLen()
      if(this%counter%key(i)%val == 'Total') cycle
      call PrintEach(this%counter%key(i)%val, this%counter%val(i), this%timer%val(i))
      r = r - this%timer%val(i)
    end do
    if (myrank==0) write(*,'(1a30, 1f12.3, 22x, 1f9.4)') &
        "misc", r, r/ttotal
    if (myrank==0)  write(*,*)
  contains
    subroutine PrintEach(title, ncall, time)
      character(*), intent(in) :: title
      integer, intent(in) :: ncall
      real(8), intent(in) :: time
      if (myrank==0) &
          write(*,'(1a30, 1f12.3, 1i10, 1f12.5, 1f9.4)') title,  &
          time, ncall, time/ncall, &
          time/ttotal
    end subroutine PrintEach

  end subroutine PrintSummary

  subroutine FinProf(this)
    class(prof), intent(inout) :: this
#ifdef MPI
    if(myrank /= 0) return
#endif
    call this%prt()
  end subroutine FinProf

end module Profiler

!program test_Profiler
!  use omp_lib
!  use Profiler, only: timer
!  implicit none
!  integer :: i, j, k
!  real(8) :: ti, tj, tk
!  real(8) :: wa
!  real(8), allocatable :: a(:), b(:)
!
!  call timer%init()
!  allocate(a(10000000))
!  a = 0.d0
!  call timer%cmemory(.true.,'a is allocated : ')
!  ti = omp_get_wtime()
!  allocate(b(100000000))
!  call timer%add("Allocation ", omp_get_wtime() - ti)
!  ti = omp_get_wtime()
!  b = 0.d0
!  call timer%add("All zero ", omp_get_wtime() - ti)
!  call timer%tmemory('b')
!
!
!  deallocate(a,b)
!
!  call timer%fin()
!
!end program test_Profiler
