!
! IEEE half-precision floating point format (https://en.wikipedia.org/wiki/Half-precision_floating-point_format)
! implementation is done based on the several web pages, which I cannot remenber.
!
module half_precision_floating_points
  use, intrinsic :: iso_fortran_env
  implicit none
  public :: my_real16
  public :: assignment(=)
  public int16_to_real32_  ! test method, don't call
  public real32_to_int16_  ! test method, don't call

  private :: summation
  private :: subtract
  private :: multiply
  private :: division
  private :: summation_real32
  private :: subtract_real32
  private :: multiply_real32
  private :: division_real32
  private :: print_my_real16
  private :: real32_to_int16
  private :: int16_to_real32
  private :: copy_my_real16

  type :: my_real16
    integer(int16) :: i16 = 0
  contains
    procedure :: summation
    procedure :: subtract
    procedure :: multiply
    procedure :: division
    procedure :: summation_real32
    procedure :: subtract_real32
    procedure :: multiply_real32
    procedure :: division_real32
    procedure :: print_my_real16
    generic :: operator(+) => summation, summation_real32
    generic :: operator(-) => subtract, subtract_real32
    generic :: operator(*) => multiply, multiply_real32
    generic :: operator(/) => division, division_real32
    generic :: prnt => print_my_real16
  end type my_real16

  interface assignment(=)
    module procedure :: real32_to_int16
    module procedure :: real64_to_int16
    module procedure :: int16_to_real32
    module procedure :: int16_to_real64
    module procedure :: int16_to_int16
    module procedure :: int16_to_int16_inv
    module procedure :: copy_my_real16
  end interface assignment(=)

  ! loopup tables should be calculated first!
  integer(int32), private :: fraction_table(0:2047)
  integer(int32), private :: exponent_table(0:63)
  integer(int16), private :: offset_table(0:63)
  integer(int16), private :: base_table(0:511)
  integer(int8), private :: shift_table(0:511)
  integer(int32), private, parameter :: z3800_32 = 939524096   ! 0x38000000
  integer(int32), private, parameter :: z8000_32 = ishft(1,31) ! 0x80000000
  integer(int32), private, parameter :: z4780_32 = 1199570944  ! 0x47800000
  integer(int32), private, parameter :: zc780_32 = -947912704  ! 0xc7800000
  integer(int32), private, parameter :: z007f_32 = 8388607     ! 0x007FFFFF
  integer(int32), private, parameter :: z0080_32 = 8388608     ! 0x00800000
  integer(int32), private, parameter :: z3880_32 = 947912704   ! 0x38800000
  integer(int16), private, parameter :: z0000_16 = 0           ! 0x0000
  integer(int16), private, parameter :: z0100_16 = 256         ! 0x0100
  integer(int16), private, parameter :: z8000_16 = -32768      ! 0x8000
  integer(int16), private, parameter :: z0400_16 = 1024        ! 0x0400
  integer(int16), private, parameter :: z03ff_16 = 1023        ! 0x03FF
  integer(int16), private, parameter :: z7c00_16 = 31744       ! 0x7C00
  integer(int16), private, parameter :: zfc00_16 = -1024       ! 0xFC00
  integer(int16), private, parameter :: z01ff_16 = 511         ! 0x01FF
contains

  elemental subroutine real32_to_int16( this, fp32 )
    type(my_real16), intent(out) :: this
    real(real32), intent(in) :: fp32
    integer(int32) :: n32
    n32 = transfer(fp32, mold=n32)
    this%i16 = base_table( iand( int( ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16)) + &
        & int( ishft( iand( n32, z007f_32), &
        & -shift_table( iand( int(ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16 ) ) ), kind=kind(this%i16))
  end subroutine real32_to_int16

  elemental subroutine real64_to_int16( this, fp64 )
    type(my_real16), intent(out) :: this
    real(real64), intent(in) :: fp64
    real(real32) :: fp32
    integer(int32) :: n32
    fp32 = real(fp64, kind=kind(fp32))
    n32 = transfer(fp32, mold=n32)
    this%i16 = base_table( iand( int( ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16)) + &
        & int( ishft( iand( n32, z007f_32), &
        & -shift_table( iand( int(ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16 ) ) ), kind=kind(this%i16))
  end subroutine real64_to_int16

  elemental subroutine int16_to_int16( this, i16 )
    type(my_real16), intent(out) :: this
    integer(int16), intent(in) :: i16
    this%i16 = i16
  end subroutine int16_to_int16

  elemental subroutine int16_to_int16_inv( i16, this )
    integer(int16), intent(out) :: i16
    type(my_real16), intent(in) :: this
    i16 = this%i16
  end subroutine int16_to_int16_inv

  elemental subroutine int16_to_real32( fp32, this  )
    real(real32), intent(out) :: fp32
    type(my_real16), intent(in) :: this
    integer(int32) :: n32
    n32 = fraction_table( offset_table( ishft(this%i16,-10) ) + iand( this%i16,z03ff_16)) + &
        & exponent_table( ishft(this%i16,-10))
    fp32 = transfer(n32,mold=fp32)
  end subroutine int16_to_real32

  elemental subroutine int16_to_real64( fp64, this  )
    real(real64), intent(out) :: fp64
    type(my_real16), intent(in) :: this
    real(real32) :: fp32
    integer(int32) :: n32
    n32 = fraction_table( offset_table( ishft(this%i16,-10) ) + iand( this%i16,z03ff_16)) + &
        & exponent_table( ishft(this%i16,-10))
    fp32 = transfer(n32,mold=fp32)
    fp64 = real( fp32, kind=kind(fp64) )
  end subroutine int16_to_real64

  elemental subroutine copy_my_real16(a, b)
    type(my_real16), intent(out) :: a
    type(my_real16), intent(in) :: b
    a%i16 = b%i16
  end subroutine copy_my_real16

  subroutine print_my_real16( this )
    class(my_real16), intent(in) :: this
    real(real32) :: x
    x = this
    write(*,*) x
  end subroutine print_my_real16

  elemental function summation(a, b) result(c)
    class(my_real16), intent(in) :: a, b
    type(my_real16) :: c
    real(real32) :: aa, bb
    aa = a
    bb = b
    c = aa+bb
  end function summation

  elemental function summation_real32(a, b) result(c)
    class(my_real16), intent(in) :: a
    real(real32), intent(in) :: b
    type(my_real16) :: c
    real(real32) :: aa
    aa = a
    c = aa+b
  end function summation_real32

  elemental function subtract(a, b) result(c)
    class(my_real16), intent(in) :: a, b
    type(my_real16) :: c
    real(real32) :: aa, bb
    aa = a
    bb = b
    c = aa-bb
  end function subtract

  elemental function subtract_real32(a, b) result(c)
    class(my_real16), intent(in) :: a
    real(real32), intent(in) :: b
    type(my_real16) :: c
    real(real32) :: aa
    aa = a
    c = aa-b
  end function subtract_real32

  elemental function multiply(a, b) result(c)
    class(my_real16), intent(in) :: a, b
    type(my_real16) :: c
    real(real32) :: aa, bb
    aa = a
    bb = b
    c = aa*bb
  end function multiply

  elemental function multiply_real32(a, b) result(c)
    class(my_real16), intent(in) :: a
    real(real32), intent(in) :: b
    type(my_real16) :: c
    real(real32) :: aa
    aa = a
    c = aa*b
  end function multiply_real32

  elemental function division(a, b) result(c)
    class(my_real16), intent(in) :: a, b
    type(my_real16) :: c
    real(real32) :: aa, bb
    aa = a
    bb = b
    c = aa/bb
  end function division

  elemental function division_real32(a, b) result(c)
    class(my_real16), intent(in) :: a
    real(real32), intent(in) :: b
    type(my_real16) :: c
    real(real32) :: aa
    aa = a
    c = aa/b
  end function division_real32


  subroutine real32_to_int16_( fp32, i16 )
    integer(int16), intent(out) :: i16
    real(real32), intent(in) :: fp32
    integer(int32) :: n32
    n32 = transfer(fp32, mold=n32)
    i16 = base_table( iand( int( ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16)) + &
        & int( ishft( iand( n32, z007f_32), &
        & -shift_table( iand( int(ishft(n32,-23),kind=kind(z01ff_16)), z01ff_16 ) ) ), kind=kind(i16))
  end subroutine real32_to_int16_

  subroutine int16_to_real32_( i16, fp32 )
    integer(int16), intent(in) :: i16
    real(real32), intent(out) :: fp32
    integer(int32) :: n32
    n32 = fraction_table( offset_table( ishft(i16,-10) ) + iand( i16,z03ff_16)) + &
        & exponent_table( ishft(i16,-10))
    fp32 = transfer(n32,mold=fp32)
  end subroutine int16_to_real32_

  subroutine set_table_for_half_precision_floating_point()
    integer(int32) :: i
    integer(int16) :: ii, ee, i_1, i_2

    fraction_table(:) = 0
    do i = 1, 1023
      fraction_table(i) = fraction_convert(i)
    end do
    do i = 1024, 2047
      fraction_table(i) = z3800_32 + ishft((i-1024), 13)
    end do

    exponent_table(:) = 0
    do i = 1, 30
      exponent_table(i) = ishft(i, 23)
    end do
    exponent_table(31) = z4780_32
    do i = 32, 62
      exponent_table(i) = z8000_32 + ishft( i-32, 23 )
    end do
    exponent_table(63) = zc780_32

    offset_table(:) = 1024
    offset_table(0)  = 0
    offset_table(32) = 0

    do ii = 0, 255
      ee = ii-127
      i_1 = int(ior(ii, z0000_16), kind=kind(i_1))
      i_2 = int(ior(ii, z0100_16), kind=kind(i_2))
      if(ee < -24) then
        base_table( i_1 ) = z0000_16
        base_table( i_2 ) = z8000_16
        shift_table( i_1 ) = 24
        shift_table( i_2 ) = 24
      elseif(ee < -14) then
        base_table( i_1 ) = ishft( z0400_16, ee+14 )
        base_table( i_2 ) = ior( z8000_16, ishft( z0400_16, ee+14 ))
        shift_table( i_1 ) = int(-ee-1, kind(shift_table(0)))
        shift_table( i_2 ) = int(-ee-1, kind(shift_table(0)))
      elseif(ee <= 15) then
        base_table( i_1 ) = int(ishft( ee+15, 10 ), kind=kind(i_1))
        base_table( i_2 ) = ior( z8000_16, int(ishft( ee+15, 10 ), kind=kind(i_2)))
        shift_table( i_1 ) = 13
        shift_table( i_2 ) = 13
      elseif(ee < 128) then
        base_table( i_1 ) = z7c00_16
        base_table( i_2 ) = zfc00_16
        shift_table( i_1 ) = 24
        shift_table( i_2 ) = 24
      else
        base_table( i_1 ) = z7c00_16
        base_table( i_2 ) = zfc00_16
        shift_table( i_1 ) = 13
        shift_table( i_2 ) = 13
      end if
    end do
    !do i = 0, 2047
    !  write(*,"(i6,4x,i10,4x,b32,4x,z8)") i, fraction_table(i), fraction_table(i), fraction_table(i)
    !end do
    !do i = 0, 63
    !  write(*,"(i6,4x,i16,4x,b32,4x,z8)") i, exponent_table(i), exponent_table(i), exponent_table(i)
    !end do
    !do i = 0, 511
    !  write(*,"(i6,4x,i16,4x,b32,4x,z8)") i, base_table(i), base_table(i), base_table(i)
    !end do
    !do i = 0, 511
    !  write(*,"(i6,4x,i16,4x,b32,4x,z8)") i, shift_table(i), shift_table(i), shift_table(i)
    !end do
  end subroutine set_table_for_half_precision_floating_point

  function fraction_convert( i ) result(r)
    integer(int32), intent(in) :: i
    integer(int32) :: r
    integer(int32) :: m, e
    m = ishft(i, 13)
    e = 0
    do
      if( iand(m,z0080_32) /= 0 ) exit
      e = e - z0080_32
      m = ishft(m, 1)
    end do
    m = iand( m, not(z0080_32))
    e = e + z3880_32
    r = ior(m,e)
  end function fraction_convert
end module half_precision_floating_points

!program test
!  use omp_lib
!  use half_precision_floating_points
!  implicit none
!
!  call set_table_for_half_precision_floating_point()
!  !call test_type()
!  !call test_array_operation()
!  call test_covered_range()
!contains
!
!  subroutine test_type()
!    type(my_real16) :: a, b, c
!    real(4) :: aa, bb, cc
!    real(8) :: d
!
!    d = 1.d0
!    aa = 0.001
!    bb = 0.0
!
!    a = real(d)
!    b = bb
!    c = a+b
!    cc = aa+bb
!    call c%prnt()
!    write(*,*) cc
!
!    c = a-b
!    cc = aa-bb
!    call c%prnt()
!    write(*,*) cc
!  end subroutine test_type
!
!  subroutine test_array_operation()
!    type(my_real16), allocatable :: x16(:,:), y16(:,:), z16(:,:)
!    real(4), allocatable :: x32(:,:), y32(:,:), z32(:,:)
!    real(4), allocatable :: test_array(:)
!    integer :: n = 30000
!    real(8) :: time
!    allocate(x32(n,n), y32(n,n), z32(n,n))
!    time = omp_get_wtime()
!    x32(:,:) = 3.0
!    y32(:,:) = 4.0
!    z32(:,:) = 5.0
!    write(*,"(a, f12.6, a)") "32-bit real 10^4 x 10^4 array settting x=3, y=4, z=5 took ", omp_get_wtime()-time, " sec"
!    test_array = transfer( x32, test_array )
!    write(*,"(a, f12.6, a)") "32-bit real 10^4 x 10^4 array needs ", size(test_array) * 4.d0/1024.d0**3, " GB memory"
!    time = omp_get_wtime()
!    z32 = x32 + y32
!    write(*,"(a, f12.6, a)") "32-bit real 10^4 x 10^4 array x+y took ", omp_get_wtime()-time, " sec"
!
!    write(*,*)
!    allocate(x16(n,n), y16(n,n), z16(n,n))
!    time = omp_get_wtime()
!    x16 = x32
!    y16 = y32
!    z16 = z32
!    write(*,"(a, f12.6, a)") "16-bit real 10^4 x 10^4 array settting x=x32, y=y32, z=z32 took ", omp_get_wtime()-time, " sec"
!    test_array = transfer( x16, test_array )
!    write(*,"(a, f12.6, a)") "16-bit real 10^4 x 10^4 array needs ", size(test_array) * 4.d0/1024.d0**3, " GB memory"
!    time = omp_get_wtime()
!    z16 = x16 + y16
!    write(*,"(a, f12.6, a)") "16-bit real 10^4 x 10^4 array x+y took ", omp_get_wtime()-time, " sec"
!
!
!    deallocate(x16,y16,z16)
!    deallocate(x32,y32,z32)
!  end subroutine test_array_operation
!
!  subroutine test_covered_range()
!    integer :: iloop
!    integer(2) :: i, j
!    real(4) :: r
!    do iloop = -32768, 32767
!    i = int(iloop,kind=kind(i))
!    call int16_to_real32_(i,r)
!    call real32_to_int16_(r,j)
!    if( i/=j ) write(*,*) "Error: ", i, j, r
!    write(*,*) i, j, r
!    !if(mod(i, 10000)==0) write(*,*) i,r
!    end do
!  end subroutine test_covered_range
!end program test
