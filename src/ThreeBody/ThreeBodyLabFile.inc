#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION

  public :: MACRO_ADD_SUFX(file_combine)
  private
#if defined(single_precision)
  real(4), allocatable :: MEs(:)
#elif defined(half_precision)
  type(my_real16), allocatable :: MEs(:)
#else
  real(8), allocatable :: MEs(:)
#endif
  integer(8) :: nelms
contains

  subroutine MACRO_ADD_SUFX(file_combine)(params, files_combined, file_output, weights_combined)
    type(InputParameters), intent(in) :: params
    type(str), intent(in) :: files_combined(:)
    type(str), intent(in) :: file_output
    real(8), intent(in) :: weights_combined(:)
    type(sys) :: s
    integer :: i

    if(s%isfile(file_output)) then
      write(*,'(2a)') trim(file_output%val), ' already exists.'
      return
    end if

    if(s%find(file_output,s%str("me3j.gz"))) then
      do i = 1, size(files_combined)
        if(s%find(files_combined(i),s%str("me3j.gz"))) cycle
        write(*,*) "Use the consistent file format!"
        stop
      end do
      call file_combine_gzip(params, files_combined, file_output, weights_combined)
      return
    end if

    call check_file_size(files_combined)

    write(*,'(3a)') "# Now reading ", trim(files_combined(1)%val), " ..."
    call read_lab_file( files_combined(1), weights_combined(1) )
    call print_head(10)
    do i = 2, size(files_combined)
      write(*,'(3a)') "# Now reading ", trim(files_combined(i)%val), " ..."
      call add_MEs(files_combined(i), weights_combined(i))
      call print_head(10)
    end do
    call write_lab_file( file_output )
    deallocate(MEs)
  end subroutine MACRO_ADD_SUFX(file_combine)

  subroutine write_lab_file(filename)
    type(str), intent(in) :: filename
    integer :: funit=20
    open(funit, file=filename%val, form='unformatted', access='stream', action='write')
    write(funit) MEs
    close(funit)
  end subroutine write_lab_file

  subroutine check_file_size(files)
    type(str), intent(in) :: files(:)
    integer :: i, funit=20
    integer(8) :: fsize1, fsize

    open(funit, file=files(1)%val, form='unformatted', access='stream', action='read',position='append')
    inquire( funit, pos=fsize )
    close(funit)

    do i = 1, size(files)
      open(funit, file=files(i)%val, form='unformatted', access='stream', action='read',position='append')
      inquire( funit, pos=fsize1 )
      close(funit)
      if( fsize1 /= fsize ) then
        write(*,'(a)') "  Error: file boundary inconsistency is detected!"
        write(*,'(5a)') "  Make sure the file boundaries of '", trim(files(1)%val), "' and '", trim(files(i)%val), "'"
        stop
      end if
    end do
  end subroutine check_file_size

  subroutine read_lab_file(filename, weight)
    type(str), intent(in) :: filename
    real(8), intent(in) :: weight
    type(sys) :: s
    integer :: funit=20
    integer(8) :: fsize, el_size, i
    real(8) :: time
#if defined(single_precision)
    real(4) :: fact
#elif defined(half_precision)
    real(4) :: fact
#else
    real(8) :: fact
#endif

    fact = real(weight, kind=kind(fact))
    if(.not. s%find(filename,s%str('.stream.bin'))) then
      write(*,*) "Only the binary stream format is supported, convert to the .stream.bin first."
      return
    end if
#if defined(single_precision)
    el_size = 4
#elif defined(half_precision)
    el_size = 2
#else
    el_size = 8
#endif

    open(funit, file=filename%val, form='unformatted', access='stream', action='read',position='append')
    inquire( funit, pos=fsize )
    close(funit)

    fsize = fsize - 1
    write(*,'(a,f12.6,a)') "# Make sure the file size: ", dble(fsize)/1024.d0**3, " GB"
    nelms = fsize / el_size
    write(*,'(a,i13)') "# Number of MEs: ", nelms
    allocate( MEs(nelms) )

    time = omp_get_wtime()
    open(funit, file=filename%val, form='unformatted', access='stream', action='read')
    read(funit) MEs
    close(funit)
    call timer%add(s%str("Read from file"), omp_get_wtime()-time)
    time = omp_get_wtime()
    !$omp parallel do
    do i = 1, nelms
      MEs(i) = MEs(i) * fact
    end do
    !$omp end parallel do
    call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
  end subroutine read_lab_file

  subroutine print_head(n)
    integer, intent(in) :: n
    integer :: i
#if defined(half_precision)
    real(4) :: ME(10)
#endif
    write(*,'(a,i4,a)') "Print first ", n, " lines"
    do i = 1, n
#if defined(half_precision)
      ME(:) = MEs((i-1)*10+1:i*10)
      write(*,'(10f12.6)') ME
#else
      write(*,'(10f12.6)') MEs((i-1)*10+1:i*10)
#endif
    end do
  end subroutine print_head

  subroutine add_MEs(filename, weight)
    type(str), intent(in) :: filename
    real(8), intent(in) :: weight
    integer(8) :: i
    integer :: funit=20
    real(8) :: time
#if defined(single_precision)
    real(4), allocatable :: tmp(:)
    real(4) :: fact
#elif defined(half_precision)
    type(my_real16), allocatable :: tmp(:)
    real(4) :: fact
#else
    real(8), allocatable :: tmp(:)
    real(8) :: fact
#endif
    type(sys) :: s

    allocate( tmp(nelms) )
    open(funit, file=filename%val, form='unformatted', access='stream', action='read')
    time = omp_get_wtime()
    read(funit) tmp
    call timer%add(s%str("Read from file"), omp_get_wtime()-time)
    close(funit)

    time = omp_get_wtime()
    fact = real(weight, kind=kind(fact))
    !$omp parallel do
    do i = 1, nelms
      MEs(i) = MEs(i) + (tmp(i) * fact)
    end do
    !$omp end parallel do
    call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
    deallocate(tmp)
  end subroutine add_MEs

  subroutine add_MEs_memory_save(filename, weight)
    type(str), intent(in) :: filename
    real(8), intent(in) :: weight
    integer(8) :: num_buffer = 1000000000, cnt, i_begin, i_end
    !integer(8) :: num_buffer = 100, cnt, i_begin, i_end
    integer :: funit=20
    real(8) :: time
#if defined(single_precision)
    real(4), allocatable :: MEs_buffer(:)
    real(4) :: fact
#elif defined(half_precision)
    type(my_real16), allocatable :: tmp(:)
    real(4), allocatable :: MEs_buffer(:)
    real(4) :: fact
#else
    real(8), allocatable :: MEs_buffer(:)
    real(8) :: fact
#endif
    type(sys) :: s

    open(funit, file=filename%val, form='unformatted', access='stream', action='read')
    allocate(MEs_buffer(num_buffer))
#if defined(half_precision)
    allocate(tmp(num_buffer))
#endif
    MEs_buffer(:) = real(0.d0, kind=kind(MEs_buffer))
    fact = real(weight, kind=kind(fact))
    if(nelms /num_buffer > 0) then
      do cnt = 1, nelms / num_buffer
        i_begin = (cnt-1)*num_buffer+1
        i_end   = cnt*num_buffer
#if defined(half_precision)
        time = omp_get_wtime()
        read(funit) tmp
        call timer%add(s%str("Read from file"), omp_get_wtime()-time)
        time = omp_get_wtime()
        MEs_buffer(:) = tmp(:)
        call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
#else
        time = omp_get_wtime()
        read(funit) MEs_buffer
        call timer%add(s%str("Read from file"), omp_get_wtime()-time)
#endif
        time = omp_get_wtime()
        MEs(i_begin:i_end) = MEs(i_begin:i_end) + fact * MEs_buffer(:)
        call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
      end do
    end if
    if( nelms - (nelms/num_buffer)*num_buffer > 0 ) then
      i_begin = (nelms/num_buffer)*num_buffer+1
      i_end   = nelms
#if defined(half_precision)
      time = omp_get_wtime()
      read(funit) tmp(: nelms-(nelms/num_buffer)*num_buffer )
      call timer%add(s%str("Read from file"), omp_get_wtime()-time)
      time = omp_get_wtime()
      MEs_buffer(: nelms-(nelms/num_buffer)*num_buffer ) = tmp(: nelms-(nelms/num_buffer)*num_buffer )
      call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
#else
      time = omp_get_wtime()
      read(funit) MEs_buffer(: nelms-(nelms/num_buffer)*num_buffer )
      call timer%add(s%str("Read from file"), omp_get_wtime()-time)
#endif
      time = omp_get_wtime()
      MEs(i_begin:i_end) = MEs(i_begin:i_end) + fact * MEs_buffer(:nelms-(nelms/num_buffer)*num_buffer )
      call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
    end if
    close(funit)
    deallocate(MEs_buffer)
#if defined(half_precision)
    deallocate(tmp)
#endif
  end subroutine add_MEs_memory_save

  subroutine file_combine_gzip(params, files_combined, file_output, weights_combined)
    type(InputParameters), intent(in) :: params
    type(str), intent(in) :: files_combined(:)
    type(str), intent(in) :: file_output
    real(8), intent(in) :: weights_combined(:)
    type(sys) :: s
    integer :: i
    type(str) :: filename
    type(OrbitsIsospin) :: sps
    real(8) :: time
#if defined(single_precision)
    real(4), allocatable :: op1(:), op2(:)
    real(4) :: fact
#elif defined(half_precision)
    real(4), allocatable :: op1(:), op2(:)
    real(4) :: fact
#else
    real(8), allocatable :: op1(:), op2(:)
    real(8) :: fact
#endif


    write(*,*) ""
    write(*,*) "Use this method only for full 3N file!"
    write(*,*) ""

    write(*,'(4(a,i3))') "# emax: ", params%emax, ", e2max: ", params%e2max, &
      & ", e3max: ", params%e3max, ", lmax: ", params%lmax

    call sps%init(params%emax, params%lmax)
    nelms = count_scalar_3bme(sps, params%e2max, params%e3max)
    write(*,'(a,i13)') "# Number of MEs: ", nelms
    allocate(op1(nelms))

    call store_3bme_gzip(op1, files_combined(1))
    time = omp_get_wtime()
    op1(:) = op1(:) * weights_combined(1)
    call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
    do i = 2, size(files_combined)
      allocate(op2(nelms))
      call store_3bme_gzip(op2, files_combined(i))
      time = omp_get_wtime()
      op1(:) = op1(:) + (op2(:) * weights_combined(i))
      call timer%add(s%str("MEs manipulation"), omp_get_wtime()-time)
      deallocate(op2)
    end do
    call write_3bme_gzip(op1, file_output)
    deallocate(op1)
    call sps%fin()
  end subroutine file_combine_gzip

  function count_scalar_3bme(spsf, e2max, e3max) result(r)
    type(OrbitsIsospin), intent(in) :: spsf
    integer, intent(in) :: e2max, e3max
    integer(8) :: r
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456

    r = 0
    do i1 = 1, spsf%norbs
      l1 = spsf%orb(i1)%l
      j1 = spsf%orb(i1)%j
      e1 = spsf%orb(i1)%e
      do i2 = 1, i1
        l2 = spsf%orb(i2)%l
        j2 = spsf%orb(i2)%j
        e2 = spsf%orb(i2)%e
        if(e1 + e2 > e2max) cycle
        do i3 = 1, i2
          l3 = spsf%orb(i3)%l
          j3 = spsf%orb(i3)%j
          e3 = spsf%orb(i3)%e
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = spsf%orb(i4)%l
            j4 = spsf%orb(i4)%j
            e4 = spsf%orb(i4)%e

            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = spsf%orb(i5)%l
              j5 = spsf%orb(i5)%j
              e5 = spsf%orb(i5)%e
              if(e4 + e5 > e2max) cycle

              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3

              do i6 = 1, i6max
                l6 = spsf%orb(i6)%l
                j6 = spsf%orb(i6)%j
                e6 = spsf%orb(i6)%e
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do J = max(abs(2*J12-j3),abs(2*J45-j6)),&
                          &min(   (2*J12+j3),   (2*J45+j6)), 2

                      do T12 = 0, 1
                        do T45 = 0, 1
                          do T = max(abs(2*T12-1),abs(2*T45-1)),&
                                &min(   (2*T12+1),   (2*T45+1)), 2

                            r = r + 1

                          end do

                        end do
                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do


        end do
      end do
    end do
  end function count_scalar_3bme

  subroutine store_3bme_gzip(MEs, f) 
    use MyLibrary, only: gzip_open, gzip_close, gzip_readline
    use, intrinsic :: iso_c_binding
#if defined(half_precision)
    real(4), intent(inout) :: MEs(:)
#elif defined(single_precision)
    real(4), intent(inout) :: MEs(:)
#elif defined(double_precision)
    real(8), intent(inout) :: MEs(:)
#endif
    type(str), intent(in) :: f
    character(256) :: header, buffer
    integer(8) :: line
    type(c_ptr) :: fp, err
    real(8) :: time
    type(sys) :: s

    time = omp_get_wtime()
    fp = gzip_open(f%val, "rt")
    err = gzip_readline(fp, header, len(header))
    do line = 1, nelms/10
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) MEs((line-1)*10+1 : line*10)
    end do
    if(nelms - (nelms/10) * 10 > 0) then
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) MEs((nelms/10)*10+1 : nelms)
    end if
    err = gzip_close(fp)
    call timer%add(s%str("Read gzipped file"), omp_get_wtime()-time)
  end subroutine store_3bme_gzip

  subroutine write_3bme_gzip(MEs, f) 
    use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
    use, intrinsic :: iso_c_binding
#if defined(half_precision)
    real(4), intent(inout) :: MEs(:)
#elif defined(single_precision)
    real(4), intent(inout) :: MEs(:)
#elif defined(double_precision)
    real(8), intent(inout) :: MEs(:)
#endif
    type(str), intent(in) :: f
    character(256) :: header, buffer
    character(12) :: cfmt
    integer(8) :: line, nrest
    type(c_ptr) :: fp, err
    real(8) :: time
    type(sys) :: s

    time = omp_get_wtime()
    fp = gzip_open(f%val, "wt")
    header = 'Created by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    err = gzip_writeline(fp, trim(header), len_trim(header))
    do line = 1, nelms/10
      write(buffer,'(10f16.8)') MEs(10*(line-1)+1:10*line)
      err = gzip_writeline(fp,trim(buffer), len_trim(buffer))
    end do
    nrest = nelms - (nelms/10)*10
    if(nrest > 0) then
      cfmt = '(xf16.8)'
      write(cfmt(2:2),'(i1)') nrest
      write(buffer,cfmt) MEs((nelms/10)*10+1:nelms)
      err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
    end if
    err = gzip_close(fp)
    call timer%add(s%str("Write gzipped file"), omp_get_wtime()-time)
  end subroutine write_3bme_gzip
