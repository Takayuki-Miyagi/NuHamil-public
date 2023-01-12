module StoreCouplings
  use omp_lib
  use ClassSys
  use Profiler, only: timer
  implicit none

  public :: CGzeroStore
  public :: CGsStore
  public :: SixJsStore
  public :: NineJsStore
  public :: TMbracketStore

  ! methods for CG coefficient
  private :: InitCGzero_final
  private :: FinCGzero_final
  private :: InitCGzeroStore
  private :: FinCGzeroStore
  private :: GetStoredCGzero
  private :: InitCGs_final
  private :: FinCGs_final
  private :: InitCGsStore
  private :: FinCGsStore
  private :: GetStoredCG
  private :: GetMemoryCG
  private :: GetMemoryCGZero

  ! methods for 6-j symbol
  private :: InitSixJsStore
  private :: FinSixJsStore
  private :: GetStoredSixJ
  private :: InitSixJs_intermediate
  private :: FinSixJs_intermediate
  private :: InitSixJs_final
  private :: FinSixJs_final
  private :: GetMemorySixJ

  ! methods for 9-j symbol
  private :: InitNineJsStore
  private :: FinNineJsStore
  private :: GetStoredNineJ
  private :: InitNineJs_intermediate
  private :: FinNineJs_intermediate
  private :: InitNineJs_final
  private :: FinNineJs_final
  private :: GetMemoryNineJ

  ! methods for Talmi-Moshinsky bracket
  private :: InitTMbracketStore
  private :: FinTMbracketStore
  private :: GetStoredTMbracket
  private :: InitTMbk_intermediate
  private :: FinTMbk_intermediate
  private :: GetMemoryTMbk

  ! CG coefficient storing
  !
  ! CG%v(j1i,j2i)%v(m1i,m2i,j12i)
  !   (  j1  j2 | j12 }
  ! = (   0   0 |   0 }
  !
  !   j1i = j1       (j1: integer)
  type :: CGzero_final
    real(8), allocatable :: v(:)
    real(8) :: mem = 0.d0
  contains
    procedure :: InitCGzero_final
    procedure :: FinCGzero_final

    generic :: init => InitCGzero_final
    generic :: fin => FinCGzero_final
  end type CGzero_final

  type :: CGzeroStore
    type(CGzero_final), allocatable :: v(:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    real(8) :: mem = 0.d0
  contains
    procedure :: InitCGzeroStore
    procedure :: FinCGzeroStore
    procedure :: GetStoredCGzero
    procedure :: GetMemoryCGZero

    generic :: init => InitCGzeroStore
    generic :: fin => FinCGzeroStore
    generic :: get => GetStoredCGzero
    generic :: GetMemory => GetMemoryCGZero
  end type CGzeroStore

  ! CG coefficient storing
  !
  ! CG%v(j1i,j2i)%v(m1i,m2i,j12i)
  !   (  j1  j2 | j12 }
  ! = (  m1  m2 | m12 }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: CGs_final
    real(8), allocatable :: v(:,:,:)
    real(8) :: mem = 0.d0
  contains
    procedure :: InitCGs_final
    procedure :: FinCGs_final

    generic :: init => InitCGs_final
    generic :: fin => FinCGs_final
  end type CGs_final

  type :: CGsStore
    type(CGs_final), allocatable :: v(:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j12
    real(8) :: mem = 0.d0
  contains
    procedure :: InitCGsStore
    procedure :: FinCGsStore
    procedure :: GetStoredCG
    procedure :: GetMemoryCG

    generic :: init => InitCGsStore
    generic :: fin => FinCGsStore
    generic :: get => GetStoredCG
    generic :: GetMemory => GetMemoryCG
  end type CGsStore

  ! 6-j symbol storing
  !
  ! sixj%v(j1i,j2i,j3i)%v(j12i,j23i)%v(ji)
  !   {  j1  j2 j12 }
  ! = {  j3   J j23 }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: SixJs_final
    real(8), allocatable :: v(:)
    real(8) :: mem = 0.d0
  contains
    procedure :: InitSixJs_final
    procedure :: FinSixJs_final

    generic :: init => InitSixJs_final
    generic :: fin => FinSixJs_final
  end type SixJs_final

  type :: SixJs_intermediate
    type(SixJs_final), allocatable :: v(:,:)
    integer :: j12min, j12max
    integer :: j23min, j23max
    real(8) :: mem = 0.d0
  contains
    procedure :: InitSixJs_intermediate
    procedure :: FinSixJs_intermediate

    generic :: init => InitSixJs_intermediate
    generic :: fin => FinSixJs_intermediate
  end type SixJs_intermediate

  type :: SixJsStore
    type(SixJs_intermediate), allocatable :: v(:,:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    integer :: j3min, j3max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j3
    logical :: half_j12
    logical :: half_j23
    logical :: half_j
    real(8) :: mem = 0.d0
  contains
    procedure :: InitSixJsStore
    procedure :: FinSixJsStore
    procedure :: GetStoredSixJ
    procedure :: GetMemorySixJ

    generic :: init => InitSixJsStore
    generic :: fin => FinSixJsStore
    generic :: get => GetStoredSixJ
    generic :: GetMemory => GetMemorySixJ
  end type SixJsStore

  ! 9-j symbol storing
  !
  ! ninej%v(j1i,j2i,j3i,j4i)%v(j12i,j34i,j13i,j24i)%v(ji)
  !   {  j1  j2 j12 }
  ! = {  j3  j4 j34 }
  !   { j13 j24   J }
  !
  !   j1i = j1       (j1: integer)
  !   j1i = (j1+1)/2 (j1: half integer)
  type :: NineJs_final
    real(8), allocatable :: v(:)
    real(8) :: mem = 0.d0
  contains
    procedure :: InitNineJs_final
    procedure :: FinNineJs_final

    generic :: init => InitNineJs_final
    generic :: fin => FinNineJs_final
  end type NineJs_final

  type :: NineJs_intermediate
    type(NineJs_final), allocatable :: v(:,:)
    integer, allocatable :: jds2i_1234(:,:), jds2i_1324(:,:)
    integer, allocatable :: i2j12d(:), i2j34d(:), i2j13d(:), i2j24d(:)
    integer :: imax_1234, imax_1324
    real(8) :: mem = 0.d0
  contains
    procedure :: InitNineJs_intermediate
    procedure :: FinNineJs_intermediate

    generic :: init => InitNineJs_intermediate
    generic :: fin => FinNineJs_intermediate
  end type NineJs_intermediate

  type :: NineJsStore
    type(NineJs_intermediate), allocatable :: v(:,:,:,:)
    integer :: j1min, j1max
    integer :: j2min, j2max
    integer :: j3min, j3max
    integer :: j4min, j4max
    logical :: half_j1
    logical :: half_j2
    logical :: half_j3
    logical :: half_j4
    logical :: half_j
    real(8) :: mem = 0.d0
  contains
    procedure :: InitNineJsStore
    procedure :: FinNineJsStore
    procedure :: GetStoredNineJ
    procedure :: GetMemoryNineJ

    generic :: init => InitNineJsStore
    generic :: fin => FinNineJsStore
    generic :: get => GetStoredNineJ
    generic :: GetMemory => GetMemoryNineJ
  end type NineJsStore

  ! Talmi-Moshinsky bracket storing
  !
  type :: TMbk_final
    real(8), allocatable :: v(:)
    real(8) :: mem = 0.d0
  end type TMbk_final

  type :: TMbk_intermediate
    type(TMbk_final), allocatable :: bk(:,:)
    integer, allocatable :: idx(:,:,:,:)
    integer, allocatable :: ll1(:)
    integer, allocatable :: ll2(:)
    integer, allocatable :: nn1(:)
    integer, allocatable :: nn2(:)
    integer :: n_idx
    real(8) :: mem = 0.d0
  contains
    procedure :: InitTMbk_intermediate
    procedure :: FinTMbk_intermediate

    generic :: init => InitTMbk_intermediate
    generic :: fin => FinTMbk_intermediate
  end type TMbk_intermediate

  type :: TMbracketStore
    type(TMbk_intermediate), allocatable :: N(:)
    real(8) :: mass_ratio
    integer :: Nmax
    real(8) :: mem = 0.d0
  contains
    procedure :: InitTMbracketStore
    procedure :: FinTMbracketStore
    procedure :: GetStoredTMbracket
    procedure :: GetMemoryTMbk

    generic :: init => InitTMbracketStore
    generic :: fin => FinTMbracketStore
    generic :: get => GetStoredTMbracket
    generic :: GetMemory => GetMemoryTMbk
  end type TMbracketStore

  type(sys), private :: sy
contains
  !
  !
  !  CG zero storing
  !
  !
  function GetMemoryCGZero(this) result(r)
    class(CGzeroStore), intent(in) :: this
    real(8) :: r
    r = this%mem
  end function GetMemoryCGZero

  subroutine FinCGzeroStore(this)
    class(CGzeroStore), intent(inout) :: this
    integer :: j1, j2
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        call this%v(j1,j2)%fin()
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinCGzeroStore

  subroutine InitCGzeroStore(this, j1min, j1max, j2min, j2max, j3min, j3max)
    class(CGzeroStore), intent(inout) :: this
    integer, intent(in) :: j1min, j1max, j2min, j2max
    integer, intent(in), optional :: j3min, j3max
    integer :: j1, j2
    real(8) :: ti

    ti = omp_get_wtime()
    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max

    allocate(this%v(j1min:j1max,j2min:j2max))

    do j1 = j1min, j1max
      do j2 = j2min, j2max
        call this%v(j1,j2)%init(j1,j2,j3min_in=j3min,j3max_in=j3max)
        this%mem = this%mem + this%v(j1,j2)%mem
      end do
    end do
  end subroutine InitCGzeroStore

  pure function GetStoredCGzero(this, j1, j2, j12) result(r)
    ! inputs are double j
    use MyLibrary, only: triag
    class(CGzeroStore), intent(in) :: this
    integer, intent(in) :: j1, j2, j12
    real(8) :: r

    r = 0.d0
    if(triag(j1,j2,j12)) return
    r = this%v(j1,j2)%v(j12)
  end function GetStoredCGzero

  subroutine FinCGzero_final(this)
    class(CGzero_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinCGzero_final

  subroutine InitCGzero_final(this,j1,j2,j3min_in,j3max_in)
    use MyLibrary, only: dcg
    class(CGzero_final), intent(inout) :: this
    integer, intent(in) :: j1, j2
    integer, intent(in), optional :: j3min_in, j3max_in
    integer :: j12
    integer :: j3min, j3max

    j3min = abs(j1-j2)
    j3max = j1+j2
    if(present(j3min_in)) j3min = j3min_in
    if(present(j3max_in)) j3max = j3max_in

    allocate(this%v( max(abs(j1-j2), j3min): min((j1+j2), j3max) ))
    do j12 = max(abs(j1-j2),j3min), min(j1+j2,j3max)
      this%v(j12) = dcg(2*j1,0,2*j2,0,2*j12,0)
    end do
    this%mem = this%mem + 8.d0 * dble(size(this%v)) / 1024.d0**3
  end subroutine InitCGzero_final
  !
  !
  !  end CG zero storing
  !
  !
  !
  !
  !  CG storing
  !
  !
  function GetMemoryCG(this) result(r)
    class(CGsStore), intent(in) :: this
    real(8) :: r
    r = this%mem
  end function GetMemoryCG

  subroutine FinCGsStore(this)
    class(CGsStore), intent(inout) :: this
    integer :: j1, j2
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        call this%v(j1,j2)%fin()
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinCGsStore

  subroutine InitCGsStore(this, j1dmin, j1dmax, half_j1, &
        & j2dmin, j2dmax, half_j2, j3dmin, j3dmax)
    class(CGsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax
    logical, intent(in) :: half_j1, half_j2
    integer, intent(in), optional :: j3dmin, j3dmax
    logical :: half_j12
    integer :: j1min, j1max, j2min, j2max
    integer :: j1, j2, j1d, j2d
    real(8) :: ti

    ti = omp_get_wtime()
    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%half_j1 = half_j1
    this%half_j2 = half_j2

    half_j12 = half_j1 .neqv. half_j2
    this%half_j12 = half_j12

    allocate(this%v(j1min:j1max,j2min:j2max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        call this%v(j1,j2)%init(j1d,j2d,half_j1,half_j2,half_j12,j3dmin,j3dmax)
        this%mem = this%mem + this%v(j1,j2)%mem
      end do
    end do
    call timer%add(sy%str('Storing CG coefs'), omp_get_wtime() - ti)
  end subroutine InitCGsStore

#ifdef debug
  function GetStoredCG(this, j1d, m1d, j2d, m2d, jd, md) result(r)
#else
  pure function GetStoredCG(this, j1d, m1d, j2d, m2d, jd, md) result(r)
#endif
    ! inputs are double j
    class(CGsStore), intent(in) :: this
    integer, intent(in) :: j1d, m1d, j2d, m2d, jd, md
    integer :: j1, j2, m1, m2, j
    real(8) :: r

    r = 0.d0
#ifdef debug
    if(m1d + m2d /= md) then
      write(*,"(a)") "Warning, GetStoredCG: m1+m2 != m3"
      return
    end if
#endif
    j1 = j1d / 2
    j2 = j2d / 2
    m1 = m1d / 2
    m2 = m2d / 2
    j = jd / 2

    if(this%half_j1 ) then
      j1  = (j1d +1) / 2
      m1  = (m1d +1) / 2
    end if
    if(this%half_j2 ) then
      j2  = (j2d +1) / 2
      m2  = (m2d +1) / 2
    end if

    if(this%half_j12 ) then
      j  = (jd +1) / 2
    end if

    r = this%v(j1,j2)%v(m1,m2,j)
  end function GetStoredCG

  subroutine FinCGs_final(this)
    class(CGs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinCGs_final

  subroutine InitCGs_final(this,j1d,j2d,half_j1,half_j2,half_j12,j3dmin,j3dmax)
    use MyLibrary, only: dcg
    class(CGs_final), intent(inout) :: this
    integer, intent(in) :: j1d, j2d
    logical, intent(in) :: half_j1, half_j2, half_j12
    integer, intent(in), optional :: j3dmin, j3dmax
    integer :: jdmax, jdmin, jmax, jmin
    integer :: m1dmin, m1dmax, m1min, m1max
    integer :: m2dmin, m2dmax, m2min, m2max
    integer :: j, jd, m1, m1d, m2, m2d

    jdmin = abs(j1d-j2d)
    jdmax =     j1d+j2d
    if(present(j3dmin)) jdmin = max(jdmin, j3dmin)
    if(present(j3dmax)) jdmax = min(jdmax, j3dmax)
    m1dmin = -j1d; m1dmax = j1d
    m2dmin = -j2d; m2dmax = j2d

    jmin = jdmin / 2
    jmax = jdmax / 2
    if(half_j12) then
      jmin = (jdmin+1) / 2
      jmax = (jdmax+1) / 2
    end if

    m1min = m1dmin / 2
    m1max = m1dmax / 2
    if(half_j1) then
      m1min = (m1dmin+1)/2
      m1max = (m1dmax+1)/2
    end if

    m2min = m2dmin / 2
    m2max = m2dmax / 2
    if(half_j2) then
      m2min = (m2dmin+1)/2
      m2max = (m2dmax+1)/2
    end if

    allocate(this%v(m1min:m1max,m2min:m2max,jmin:jmax))
    do m1 = m1min, m1max
      m1d = 2*m1
      if(half_j1) m1d = 2*m1-1
      do m2 = m2min, m2max
        m2d = 2*m2
        if(half_j2) m2d = 2*m2-1
        do j = jmin, jmax
          jd = 2*j
          if(half_j12) jd = 2*j-1
          this%v(m1,m2,j) = dcg(j1d,m1d,j2d,m2d,jd,m1d+m2d)
        end do
      end do
    end do
    this%mem = this%mem + 8.d0 * dble(size(this%v)) / 1024.d0**3
  end subroutine InitCGs_final

  !
  !
  !  end CG storing
  !
  !
  !
  !
  !  6-j storing
  !
  !
  function GetMemorySixJ(this) result(r)
    class(SixJsStore), intent(in) :: this
    real(8) :: r
    r = this%mem
  end function GetMemorySixJ

  subroutine FinSixJsStore(this)
    class(SixJsStore), intent(inout) :: this
    integer :: j1, j2, j3
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        do j3 = this%j3min, this%j3max
          call this%v(j1,j2,j3)%fin()
        end do
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinSixJsStore

  subroutine InitSixJsStore(this, j1dmin, j1dmax, half_j1, &
        & j2dmin, j2dmax, half_j2, j3dmin, j3dmax, half_j3, &
        & j12dmin_in, j12dmax_in, j23dmin_in, j23dmax_in, j123dmin_in, j123dmax_in)
    class(SixJsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax, j3dmin, j3dmax
    logical, intent(in) :: half_j1, half_j2, half_j3
    integer, intent(in), optional :: j12dmin_in, j12dmax_in, j23dmin_in, j23dmax_in, j123dmin_in, j123dmax_in
    logical :: half_j12, half_j23, half_j123, half_j231
    integer :: j1min, j1max, j2min, j2max, j3min, j3max
    integer :: j1, j2, j3, j1d, j2d, j3d
    real(8) :: ti

    ti = omp_get_wtime()

    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    j3min = j3dmin/2
    j3max = j3dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2
    if(half_j3) j3min = (j3dmin+1)/2
    if(half_j3) j3max = (j3dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%j3min = j3min
    this%j3max = j3max
    this%half_j1 = half_j1
    this%half_j2 = half_j2
    this%half_j3 = half_j3

    half_j12 = half_j1 .neqv. half_j2
    half_j23 = half_j2 .neqv. half_j3
    half_j123 = half_j12 .neqv. half_j3
    half_j231 = half_j23 .neqv. half_j1
    if(half_j123 .neqv. half_j231) then
      write(*,*) "Error occures in 6-j symbol storing"
      return
    end if

    this%half_j12 = half_j12
    this%half_j23 = half_j23
    this%half_j = half_j123

    allocate(this%v(j1min:j1max,j2min:j2max,j3min:j3max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        do j3 = j3min, j3max
          j3d = 2*j3
          if(half_j3) j3d = 2*j3-1
          call this%v(j1,j2,j3)%init(j1d,j2d,j3d, &
              & half_j1, half_j3, half_j12,half_j23,&
              & j12dmin_in=j12dmin_in, j12dmax_in=j12dmax_in, &
              & j23dmin_in=j23dmin_in, j23dmax_in=j23dmax_in, &
              & j123dmin_in=j123dmin_in, j123dmax_in=j123dmax_in)
          this%mem = this%mem + this%v(j1,j2,j3)%mem
        end do
      end do
    end do
    call timer%add(sy%str('Storing 6-j couplings'), omp_get_wtime() - ti)
  end subroutine InitSixJsStore

#ifdef debug
  function GetStoredSixJ(this, j1d, j2d, j12d, j3d, jd, j23d) result(r)
#else
  pure function GetStoredSixJ(this, j1d, j2d, j12d, j3d, jd, j23d) result(r)
#endif
    ! inputs are double j
    class(SixJsStore), intent(in) :: this
    integer, intent(in) :: j1d, j2d, j12d, j3d, jd, j23d
    integer :: j1, j2, j12, j3, j, j23
    real(8) :: r

    j1 = j1d / 2
    j2 = j2d / 2
    j3 = j3d / 2
    j12 = j12d / 2
    j23 = j23d / 2
    j = jd / 2
    if(this%half_j1 ) j1  = (j1d +1) / 2
    if(this%half_j2 ) j2  = (j2d +1) / 2
    if(this%half_j3 ) j3  = (j3d +1) / 2
    if(this%half_j12) j12 = (j12d+1) / 2
    if(this%half_j23) j23 = (j23d+1) / 2
    if(this%half_j  ) j   = (jd  +1) / 2

    r = this%v(j1,j2,j3)%v(j12,j23)%v(j)
  end function GetStoredSixJ

  subroutine FinSixJs_intermediate(this)
    class(SixJs_intermediate), intent(inout) :: this
    integer :: j12, j23
    do j12 = this%j12min, this%j12max
      do j23 = this%j23min, this%j23max
        call this%v(j12,j23)%fin()
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinSixJs_intermediate

  subroutine InitSixJs_intermediate(this,j1d,j2d,j3d,half_j1,half_j3,&
        & half_j12,half_j23, j12dmin_in, j12dmax_in, j23dmin_in, j23dmax_in, j123dmin_in, j123dmax_in)
    class(Sixjs_intermediate), intent(inout) :: this
    integer, intent(in) :: j1d, j2d, j3d
    logical, intent(in) :: half_j12, half_j23, half_j1, half_j3
    integer, intent(in), optional :: j12dmin_in, j12dmax_in, j23dmin_in, j23dmax_in, j123dmin_in, j123dmax_in
    integer :: j12, j23, j12d, j23d
    integer :: j12dmin, j12dmax, j12min, j12max
    integer :: j23dmin, j23dmax, j23min, j23max
    logical :: half_j, half_j123, half_j231

    j12dmin = abs(j1d - j2d)
    j12dmax =     j1d + j2d
    j23dmin = abs(j2d - j3d)
    j23dmax =     j2d + j3d
    if(present(j12dmin_in)) j12dmin = max(j12dmin, j12dmin_in)
    if(present(j12dmax_in)) j12dmax = min(j12dmax, j12dmax_in)
    if(present(j23dmin_in)) j23dmin = max(j23dmin, j23dmin_in)
    if(present(j23dmax_in)) j23dmax = min(j23dmax, j23dmax_in)

    j12min = j12dmin/2
    j12max = j12dmax/2
    j23min = j23dmin/2
    j23max = j23dmax/2
    if(half_j12) then
      j12min = (j12dmin+1)/2
      j12max = (j12dmax+1)/2
    end if

    if(half_j23) then
      j23min = (j23dmin+1)/2
      j23max = (j23dmax+1)/2
    end if

    half_j123 = half_j12 .neqv. half_j3
    half_j231 = half_j23 .neqv. half_j1
    if(half_j123 .neqv. half_j231) then
      write(*,*) "Error occures in 6-j symbol storing"
      return
    end if
    half_j = half_j123

    this%j12min = j12min
    this%j12max = j12max
    this%j23min = j23min
    this%j23max = j23max

    allocate(this%v(j12min:j12max,j23min:j23max))
    do j12 = this%j12min, this%j12max
      j12d = 2*j12
      if(half_j12) j12d = 2*j12-1
      do j23 = this%j23min, this%j23max
        j23d = 2*j23
        if(half_j23) j23d = 2*j23-1
        call this%v(j12,j23)%init(j1d,j2d,j3d,j12d,j23d,half_j,&
            & j123dmin_in=j123dmin_in, j123dmax_in=j123dmax_in)
        this%mem = this%v(j12,j23)%mem
      end do
    end do
  end subroutine InitSixJs_intermediate

  subroutine FinSixJs_final(this)
    class(Sixjs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinSixJs_final

  subroutine InitSixJs_final(this,j1d,j2d,j3d,j12d,j23d,half_j,&
        & j123dmin_in, j123dmax_in)
    use MyLibrary, only: sjs
    class(Sixjs_final), intent(inout) :: this
    integer, intent(in) :: j1d, j2d, j3d, j12d, j23d
    logical, intent(in) :: half_j
    integer, intent(in), optional :: j123dmin_in, j123dmax_in
    integer :: jdmax, jdmin, jd, jmax, jmin, j
    jdmin = max(abs(j12d-j3d), abs(j1d-j23d))
    jdmax = min(   (j12d+j3d),    (j1d+j23d))
    if(present(j123dmin_in)) jdmin = max(jdmin, j123dmin_in)
    if(present(j123dmax_in)) jdmax = min(jdmax, j123dmax_in)
    if(jdmin > jdmax) return
    jmin = jdmin / 2
    jmax = jdmax / 2
    if(half_j) then
      jmin = (jdmin+1) / 2
      jmax = (jdmax+1) / 2
    end if

    allocate(this%v(jmin:jmax))
    do j = jmin, jmax
      jd = 2*j
      if(half_j) jd = 2*j - 1
      this%v(j) = sjs(j1d,j2d,j12d,j3d,jd,j23d)
    end do
    this%mem = this%mem + 8.d0 * dble(size(this%v)) / 1024.d0**3
  end subroutine InitSixJs_final
  !
  !
  !  end 6-j storing
  !
  !

  !
  !
  !  9-j storing
  !
  !
  function GetMemoryNineJ(this) result(r)
    class(NineJsStore), intent(in) :: this
    real(8) :: r
    r = this%mem
  end function GetMemoryNineJ

  subroutine FinNineJsStore(this)
    class(NineJsStore), intent(inout) :: this
    integer :: j1, j2, j3, j4
    do j1 = this%j1min, this%j1max
      do j2 = this%j2min, this%j2max
        do j3 = this%j3min, this%j3max
          do j4 = this%j4min, this%j4max
            call this%v(j1,j2,j3,j4)%fin()
          end do
        end do
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinNineJsStore

  subroutine InitNineJsStore(this, j1dmin, j1dmax, half_j1, j2dmin, j2dmax, half_j2, &
        & j3dmin, j3dmax, half_j3, j4dmin, j4dmax, half_j4, &
        & j12dmin_in, j12dmax_in, j34dmin_in, j34dmax_in, &
        & j13dmin_in, j13dmax_in, j24dmin_in, j24dmax_in, jdmin_in, jdmax_in)
    class(NineJsStore), intent(inout) :: this
    integer, intent(in) :: j1dmin, j1dmax, j2dmin, j2dmax, j3dmin, j3dmax, j4dmin, j4dmax
    logical, intent(in) :: half_j1, half_j2, half_j3, half_j4
    integer, intent(in), optional :: j12dmin_in, j12dmax_in, j34dmin_in, j34dmax_in
    integer, intent(in), optional :: j13dmin_in, j13dmax_in, j24dmin_in, j24dmax_in, jdmin_in, jdmax_in
    logical :: half_j12, half_j34, half_j13, half_j24, half_j1234, half_j1324
    integer :: j1min, j1max, j2min, j2max, j3min, j3max, j4min, j4max
    integer :: j1d, j2d, j3d, j4d
    integer :: j1, j2, j3, j4
    real(8) :: ti

    ti = omp_get_wtime()

    j1min = j1dmin/2
    j1max = j1dmax/2
    j2min = j2dmin/2
    j2max = j2dmax/2
    j3min = j3dmin/2
    j3max = j3dmax/2
    j4min = j4dmin/2
    j4max = j4dmax/2
    if(half_j1) j1min = (j1dmin+1)/2
    if(half_j1) j1max = (j1dmax+1)/2
    if(half_j2) j2min = (j2dmin+1)/2
    if(half_j2) j2max = (j2dmax+1)/2
    if(half_j3) j3min = (j3dmin+1)/2
    if(half_j3) j3max = (j3dmax+1)/2
    if(half_j4) j4min = (j4dmin+1)/2
    if(half_j4) j4max = (j4dmax+1)/2

    this%j1min = j1min
    this%j1max = j1max
    this%j2min = j2min
    this%j2max = j2max
    this%j3min = j3min
    this%j3max = j3max
    this%j4min = j4min
    this%j4max = j4max
    this%half_j1 = half_j1
    this%half_j2 = half_j2
    this%half_j3 = half_j3
    this%half_j4 = half_j4

    half_j12 = half_j1 .neqv. half_j2
    half_j34 = half_j3 .neqv. half_j4
    half_j13 = half_j1 .neqv. half_j3
    half_j24 = half_j2 .neqv. half_j4

    half_j1234 = half_j12 .neqv. half_j34
    half_j1324 = half_j13 .neqv. half_j24
    if(half_j1234 .neqv. half_j1324) then
      write(*,*) "Error occures in 9-j symbol storing"
      return
    end if
    this%half_j = half_j1234

    allocate(this%v(j1min:j1max,j2min:j2max,j3min:j3max,j4min:j4max))

    do j1 = j1min, j1max
      j1d = 2*j1
      if(half_j1) j1d = 2*j1-1
      do j2 = j2min, j2max
        j2d = 2*j2
        if(half_j2) j2d = 2*j2-1
        do j3 = j3min, j3max
          j3d = 2*j3
          if(half_j3) j3d = 2*j3-1
          do j4 = j4min, j4max
            j4d = 2*j4
            if(half_j4) j4d = 2*j4-1
            call this%v(j1,j2,j3,j4)%init(j1d,j2d,j3d,j4d, &
                & half_j12,half_j34,&
                & j12dmin_in=j12dmin_in, j12dmax_in=j12dmax_in,&
                & j34dmin_in=j34dmin_in, j34dmax_in=j34dmax_in,&
                & j13dmin_in=j13dmin_in, j13dmax_in=j13dmax_in,&
                & j24dmin_in=j24dmin_in, j24dmax_in=j24dmax_in,&
                & jdmin_in=jdmin_in, jdmax_in=jdmax_in)
            this%mem = this%mem + this%v(j1,j2,j3,j4)%mem
          end do
        end do
      end do
    end do

    call timer%add(sy%str('Storing 9-j couplings'), omp_get_wtime() - ti)
  end subroutine InitNineJsStore

#ifdef debug
  function GetStoredNineJ(this, j1d, j2d, j12d, j3d, j4d, j34d, j13d, j24d, jd) result(r)
#else
  pure function GetStoredNineJ(this, j1d, j2d, j12d, j3d, j4d, j34d, j13d, j24d, jd) result(r)
#endif
    ! inputs are double j
    class(NineJsStore), intent(in) :: this
    integer, intent(in) :: j1d, j2d, j3d, j4d, j12d, j34d, j13d, j24d, jd
    integer :: j1, j2, j3, j4, j
    integer :: i1234, i1324
    real(8) :: r
    j1 = j1d / 2
    j2 = j2d / 2
    j3 = j3d / 2
    j4 = j4d / 2
    j = jd / 2

    if(this%half_j1 ) j1  = (j1d +1) / 2
    if(this%half_j2 ) j2  = (j2d +1) / 2
    if(this%half_j3 ) j3  = (j3d +1) / 2
    if(this%half_j4 ) j4  = (j4d +1) / 2
    if(this%half_j  ) j   = (jd  +1) / 2

    i1234 = this%v(j1,j2,j3,j4)%jds2i_1234(j12d,j34d)
    i1324 = this%v(j1,j2,j3,j4)%jds2i_1324(j13d,j24d)
#ifdef debug
    if(i1234 * i1324 == 0) then
      write(*,"(a)") "Warning: GetStoredNineJ"
      r = 0.d0
      return
    end if
#endif

    r = this%v(j1,j2,j3,j4)%v(i1234,i1324)%v(j)
  end function GetStoredNineJ

  subroutine FinNineJs_intermediate(this)
    class(NineJs_intermediate), intent(inout) :: this
    integer :: i1234, i1324
    do i1234 = 1, this%imax_1234
      do i1324 = 1, this%imax_1324
        call this%v(i1234,i1324)%fin()
      end do
    end do
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinNineJs_intermediate

  subroutine InitNineJs_intermediate(this,j1d,j2d,j3d,j4d,&
        & half_j12,half_j34,&
        & j12dmin_in, j12dmax_in, j34dmin_in, j34dmax_in, &
        & j13dmin_in, j13dmax_in, j24dmin_in, j24dmax_in, jdmin_in, jdmax_in)
    use MyLibrary, only: triag
    class(NineJs_intermediate), intent(inout) :: this
    integer, intent(in) :: j1d,j2d,j3d,j4d
    logical, intent(in) :: half_j12, half_j34
    integer, intent(in), optional :: j12dmin_in, j12dmax_in, j34dmin_in, j34dmax_in
    integer, intent(in), optional :: j13dmin_in, j13dmax_in, j24dmin_in, j24dmax_in, jdmin_in, jdmax_in
    logical :: half_j
    integer :: j12d,j34d,j13d,j24d
    integer :: j12dmin,j12dmax
    integer :: j34dmin,j34dmax
    integer :: j13dmin,j13dmax
    integer :: j24dmin,j24dmax
    integer :: cnt, i1234, i1324

    j12dmin = abs(j1d-j2d)
    j34dmin = abs(j3d-j4d)
    j13dmin = abs(j1d-j3d)
    j24dmin = abs(j2d-j4d)
    j12dmax =     j1d+j2d
    j34dmax =     j3d+j4d
    j13dmax =     j1d+j3d
    j24dmax =     j2d+j4d
    if(present(j12dmin_in)) j12dmin = max(j12dmin, j12dmin_in)
    if(present(j12dmax_in)) j12dmax = min(j12dmax, j12dmax_in)
    if(present(j34dmin_in)) j34dmin = max(j34dmin, j34dmin_in)
    if(present(j34dmax_in)) j34dmax = min(j34dmax, j34dmax_in)
    if(present(j13dmin_in)) j13dmin = max(j13dmin, j13dmin_in)
    if(present(j13dmax_in)) j13dmax = min(j13dmax, j13dmax_in)
    if(present(j24dmin_in)) j24dmin = max(j24dmin, j24dmin_in)
    if(present(j24dmax_in)) j24dmax = min(j24dmax, j24dmax_in)

    half_j = half_j12 .neqv. half_j34
    cnt = 0
    do j12d = j12dmin, j12dmax, 2
      do j34d = j34dmin, j34dmax, 2
        cnt = cnt + 1
      end do
    end do
    this%imax_1234 = cnt
    allocate(this%jds2i_1234(j12dmin:j12dmax, j34dmin:j34dmax))
    allocate(this%i2j12d(this%imax_1234))
    allocate(this%i2j34d(this%imax_1234))
    this%jds2i_1234(:,:) = 0
    cnt = 0
    do j12d = j12dmin, j12dmax, 2
      do j34d = j34dmin, j34dmax, 2
        cnt = cnt + 1
        this%jds2i_1234(j12d,j34d) = cnt
        this%i2j12d(cnt) = j12d
        this%i2j34d(cnt) = j34d
      end do
    end do

    cnt = 0
    do j13d = j13dmin, j13dmax, 2
      do j24d = j24dmin, j24dmax, 2
        cnt = cnt + 1
      end do
    end do
    this%imax_1324 = cnt
    allocate(this%jds2i_1324(j13dmin:j13dmax, j24dmin:j24dmax))
    allocate(this%i2j13d(this%imax_1324))
    allocate(this%i2j24d(this%imax_1324))
    this%jds2i_1324(:,:) = 0
    cnt = 0
    do j13d = j13dmin, j13dmax, 2
      do j24d = j24dmin, j24dmax, 2
        cnt = cnt + 1
        this%jds2i_1324(j13d,j24d) = cnt
        this%i2j13d(cnt) = j13d
        this%i2j24d(cnt) = j24d
      end do
    end do

    allocate(this%v(this%imax_1234,this%imax_1324))
    do i1234 = 1, this%imax_1234
      j12d = this%i2j12d(i1234)
      j34d = this%i2j34d(i1234)
      do i1324 = 1, this%imax_1324
        j13d = this%i2j13d(i1324)
        j24d = this%i2j24d(i1324)
        call this%v(i1234,i1324)%init(j1d,j2d,j3d,j4d,&
            & j12d,j34d,j13d,j24d,half_j,jdmin_in=jdmin_in, jdmax_in=jdmax_in)
        this%mem = this%mem + this%v(i1234,i1324)%mem
      end do
    end do
    this%mem = this%mem + 4.d0 * &
      & dble(size(this%jds2i_1234) + 2*size(this%i2j12d) + &
      & size(this%jds2i_1324) + 2*size(this%i2j13d)) / 1024.d0**3
  end subroutine InitNineJs_intermediate

  subroutine FinNineJs_final(this)
    class(NineJs_final), intent(inout) :: this
    if(.not. allocated(this%v)) return
    deallocate(this%v)
    this%mem = 0.d0
  end subroutine FinNineJs_final

  subroutine InitNineJs_final(this,j1d,j2d,j3d,j4d,j12d,j34d,j13d,j24d,half_j,jdmin_in,jdmax_in)
    use MyLibrary, only: snj
    class(NineJs_final), intent(inout) :: this
    integer, intent(in) :: j1d,j2d,j3d,j4d,j12d,j34d,j13d,j24d
    logical, intent(in) :: half_j
    integer, intent(in), optional :: jdmin_in, jdmax_in
    integer :: j, jj
    integer :: jmin=-1, jmax=-1
    integer :: jdmin, jdmax


    jdmin = max(abs(j12d-j34d),abs(j13d-j24d))
    jdmax = min(   (j12d+j34d),   (j13d+j24d))
    if(present(jdmin_in)) jdmin = max(jdmin, jdmin_in)
    if(present(jdmax_in)) jdmax = min(jdmax, jdmax_in)

    if(jdmin > jdmax) then
      return
    end if

    jmin = jdmin/2
    jmax = jdmax/2
    if(half_j) then
      jmin = (jdmin + 1)/2
      jmax = (jdmax + 1)/2
    end if

    allocate(this%v(jmin:jmax))
    do j = jmin, jmax
      jj = 2*j
      if(half_j) jj = 2*j-1
      this%v(j) = snj(j1d,j2d,j12d,j3d,j4d,j34d,j13d,j24d,jj)
    end do
    this%mem = this%mem + 8.d0 * dble(size(this%v)) / 1024.d0**3
  end subroutine InitNineJs_final

  !
  !
  !  end 9-j storing
  !
  !

  !
  !
  !  Talmi-Moshinsky bracket storing
  !
  !
  function GetMemoryTMbk(this) result(r)
    class(TMbracketStore), intent(in) :: this
    real(8) :: r
    r = this%mem
  end function GetMemoryTMbk

  subroutine FinTMbracketStore(this)
    class(TMbracketStore), intent(inout) :: this
    integer :: N
    do N = 0, this%Nmax
      call this%N(N)%fin()
    end do
    deallocate(this%N)
    this%mem = 0.d0
  end subroutine FinTMbracketStore

  subroutine InitTMbracketStore(this, Nmax, d)
    class(TMbracketStore), intent(inout) :: this
    integer, intent(in) :: Nmax
    real(8), intent(in) :: d
    integer :: Ntot
    real(8) :: ti

    ti = omp_get_wtime()

    this%mass_ratio = d
    this%Nmax = Nmax
    allocate(this%N(0:Nmax))
    do ntot = 0, Nmax
      call this%N(ntot)%init(ntot, d)
      this%mem = this%mem + this%N(ntot)%mem
    end do

    call timer%add(sy%str('Storing Talmi-Moshinsky brackets'), omp_get_wtime() - ti)
  end subroutine InitTMbracketStore

#ifdef debug
  function GetStoredTMbracket(this, n1, l1, n2, l2, n3, l3, n4, l4, lam) result(r)
#else
  pure function GetStoredTMbracket(this, n1, l1, n2, l2, n3, l3, n4, l4, lam) result(r)
#endif
    class(TMbracketStore), intent(in) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, lam
    integer :: Nmax, idx_bra, idx_ket
    real(8) :: r

    r = 0.d0
    Nmax = 2*n1 + l1 + 2*n2 + l2
    if(Nmax /= 2*n3 + l3 + 2*n4 + l4) return

    idx_bra = this%N(Nmax)%idx(n1,l1,n2,l2)
    idx_ket = this%N(Nmax)%idx(n3,l3,n4,l4)
#ifdef debug
    if(idx_bra * idx_ket == 0) then
      write(*,*) "error : GetStoreadTMbracket !"
      return
    end if
#endif
    r = this%N(Nmax)%bk(idx_bra,idx_ket)%v(lam)
  end function GetStoredTMbracket

  subroutine FinTMbk_intermediate(this)
    class(TMbk_intermediate), intent(inout) :: this
    integer :: bra, ket
    do bra = 1, this%n_idx
      do ket = 1, this%n_idx
        deallocate(this%bk(bra,ket)%v)
      end do
    end do
    deallocate(this%bk)
    deallocate(this%ll1)
    deallocate(this%ll2)
    deallocate(this%nn1)
    deallocate(this%nn2)
    deallocate(this%idx)
    this%mem = 0.d0
  end subroutine FinTMbk_intermediate

  subroutine InitTMbk_intermediate(this, Nmax, d)
    use MyLibrary, only: gmosh
    class(TMbk_intermediate), intent(inout) :: this
    integer, intent(in) :: Nmax
    real(8), intent(in) :: d
    integer :: idx, n1, l1, n2, l2, n
    integer :: bra, ket, n3, l3, n4, l4
    integer :: lammin, lammax, lam
    real(8) :: mem

    allocate(this%idx(0:Nmax/2, 0:Nmax, 0:Nmax/2, 0:Nmax))
    idx = 0
    do n1 = 0, Nmax / 2
      do l1 = 0, Nmax - 2 * n1
        do n2 = 0, (Nmax - 2 * n1 - l1) / 2
          l2 = Nmax - 2 * n1 - 2 * n2 - l1
          idx = idx + 1
        end do
      end do
    end do
    n = idx
    this%n_idx = n
    allocate(this%bk(n,n))
    allocate(this%ll1(n))
    allocate(this%ll2(n))
    allocate(this%nn1(n))
    allocate(this%nn2(n))

    idx = 0
    do n1 = 0, Nmax / 2
      do l1 = 0, Nmax - 2 * n1
        do n2 = 0, (Nmax - 2 * n1 - l1) / 2
          l2 = Nmax - 2 * n1 - 2 * n2 - l1
          idx = idx + 1
          this%ll1(idx) = l1
          this%ll2(idx) = l2
          this%nn1(idx) = n1
          this%nn2(idx) = n2
          this%idx(n1, l1, n2, l2) = idx
        end do
      end do
    end do

    mem = 0.d0
    n = this%n_idx
    !$omp parallel
    !$omp do private(bra, l1, l2, n1, n2, ket, l3, l4, n3, n4, &
    !$omp &          lammin, lammax, lam) reduction(+:mem) schedule(dynamic)
    do bra = 1, n
      l1 = this%ll1(bra)
      l2 = this%ll2(bra)
      n1 = this%nn1(bra)
      n2 = this%nn2(bra)
      do ket = 1, n
        l3 = this%ll1(ket)
        l4 = this%ll2(ket)
        n3 = this%nn1(ket)
        n4 = this%nn2(ket)

        lammin = max(abs(l1 - l2), abs(l3 - l4))
        lammax = min(l1 + l2, l3 + l4)
        allocate(this%bk(bra, ket)%v(lammin:lammax))
        this%bk(bra, ket)%v(:) = 0.d0
        do lam = lammin, lammax
          this%bk(bra, ket)%v(lam)    = gmosh(n1, l1, n2, l2, n3, l3, n4, l4, lam, d)
        end do
        mem = mem + 8.d0 * dble(size(this%bk(bra,ket)%v)) / 1024.d0**3
      end do
    end do
    !$omp end do
    !$omp end parallel
    this%mem = mem + 4.d0 * dble(size(this%idx) + 4*size(this%ll1)) / 1024.d0**3
  end subroutine InitTMbk_intermediate

  function two_body_tcoef_pn(na, la, ja, za, nb, lb, jb, zb, Ncm, Lcm, nr, lr, s, jr, J, Tz) result(r)
    use MyLibrary, only: sjs, snj, gmosh
    integer, intent(in) :: na, la, ja, za, nb, lb, jb, zb
    integer, intent(in) :: Ncm, Lcm, nr, lr, s, jr
    integer, intent(in) :: J, Tz
    real(8) :: r, ft, delab
    ! loop variables
    integer :: lam, lammin, lammax
    r = 0.d0
    ft = 1.d0
    if(za+zb /= 2*Tz) return
    if(abs(Tz)==1) ft = (1.d0 + (-1.d0)**(la + lb + s - Lcm)) / sqrt(2.d0)
    if(abs(ft) < 1.d-8) return
    delab = 1.d0
    if(na==nb .and. la==lb .and. ja==jb .and. za==zb) delab = 1.d0/sqrt(2.d0)
    lammin = max(abs(lcm - lr), abs(J - s), abs(la - lb))
    lammax = min(lcm + lr, J + s, la + lb)
    r = 0.d0
    do lam = lammin, lammax
      r = r + (-1.d0) **(lcm + lr + s + J) * &
          & sqrt(dble(2*lam+1) * dble(2*s+1) * dble(ja+1) * dble(jb+1)) * &
          & snj(2*la, 2*lb, 2*lam, 1, 1, 2*s, ja, jb, 2*J) * &
          & sqrt(dble(2*lam+1) * dble(2*jr+1)) * &
          & sjs(2*lcm, 2*lr, 2*lam, 2*s, 2*J, 2*jr) * &
          & gmosh(ncm, lcm, nr, lr, na, la, nb, lb, lam, 1.d0) * &
          & ft * delab
    end do
  end function two_body_tcoef_pn

  function three_body_tcoef_ls_iso(na, la, ja, nb, lb, jb, nc, lc, jc, tab, jab, &
        & Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jr, J) result(r)
    use MyLibrary, only: sjs, snj, gmosh
    integer, intent(in) :: na, la, ja, nb, lb, jb, nc, lc, jc, tab, jab
    integer, intent(in) :: Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jr
    integer, intent(in) :: J
    real(8) :: r
    real(8) :: snj_ab, rel_gmosh, cm_gmosh, sjs_lc, snj_c
    ! loop variables
    integer :: lab, lab_min, lab_max
    integer :: lit, lit_min, lit_max, nit
    integer :: lam, lam_min, lam_max
    integer :: ltot, ltot_min, ltot_max
    integer :: stot, stot_min, stot_max
    integer :: lrtot, lrtot_min, lrtot_max
    r = 0.d0
    if(tab /= t12) return
    if(2*(na+nb+nc)+la+lb+lc /= 2*(Ncm+n12+n3)+Lcm+l12+l3) return
    if(mod(la+lb+lc+Lcm+l12+l3,2)==1) return

    r = 0.d0
    lab_min = max(iabs(jab-s12), iabs(la-lb))
    lab_max = min(     jab+s12,       la+lb)
    do lab = lab_min, lab_max
      snj_ab = snj(2*la, 1, ja, 2*lb, 1, jb, 2*lab, 2*s12, 2*jab)
      if(abs(snj_ab) < 1.d-16) cycle

      lit_min = max(iabs(lab-l12), iabs(jab-j12))
      lit_max = min(     lab+l12,       jab+j12, 2*na + la + 2*nb + lb - 2*n12 - l12)
      do lit = lit_min, lit_max
        nit = (2*na + la + 2*nb + lb - 2*n12 - l12 - lit) / 2
        rel_gmosh = gmosh(nit, lit, n12, l12, na, la, nb, lb, lab, 1.d0)
        if(abs(rel_gmosh) < 1.d-16) cycle

        lam_min = max(iabs(lit-lc), iabs(lcm-l3))
        lam_max = min(     lit+lc,       lcm+l3)
        do lam = lam_min, lam_max
          cm_gmosh = gmosh(ncm, lcm, n3, l3, nit, lit, nc, lc, lam, 2.d0)
          if(abs(cm_gmosh) < 1.d-16) cycle

          ltot_min = max(iabs(lam-l12), iabs(lab-lc))
          ltot_max = min(     lam+l12,       lab+lc)
          do ltot = ltot_min, ltot_max
            sjs_lc = sjs(2*l12, 2*lit, 2*lab, 2*lc, 2*ltot, 2*lam)
            if(abs(sjs_lc) < 1.d-16) cycle

            stot_min = max(iabs(J-2*ltot), iabs(2*s12-1))
            stot_max = min(     J+2*ltot,       2*s12+1)
            do stot = stot_min, stot_max, 2
              snj_c = snj(2*lab,2*s12,2*jab,2*lc,1,jc,2*ltot,stot,J )
              if(abs(snj_c) < 1.d-16) cycle

              lrtot_min = max(iabs(l12-l3), iabs(stot-jr)/2, iabs(ltot-lcm))
              lrtot_max = min(     l12+l3,      (stot+jr)/2,      ltot+lcm)
              do lrtot = lrtot_min, lrtot_max
                r = r + (-1.d0) ** (lc + lam + lab + ltot + l12 + (stot + J)/2) * &
                    & sqrt(dble((ja+1)*(jb+1)*(jc+1)*(2*jab+1)*(jr+1)*(2*j12+1)*(j3+1)*(2*s12+1))) * &
                    & dble((stot+1)*(2*lab+1)*(2*ltot+1)*(2*lrtot+1)*(2*lam+1)) * &
                    & rel_gmosh * cm_gmosh * snj_ab * snj_c * sjs_lc * &
                    & snj(2*l12,2*s12,2*j12,2*l3,1,j3,2*lrtot,stot,jr) * &
                    & sjs(2*lcm, 2*l3, 2*lam, 2*l12, 2*ltot, 2*lrtot) * &
                    & sjs(2*lcm, 2*lrtot, 2*ltot, stot, J, jr)
              end do
            end do
          end do
        end do
      end do
    end do
  end function three_body_tcoef_ls_iso

  function three_body_tcoef_jj_iso(na, la, ja, nb, lb, jb, nc, lc, jc, tab, jab, &
        & Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jr, J) result(r)
    use MyLibrary, only: sjs, snj, gmosh
    integer, intent(in) :: na, la, ja, nb, lb, jb, nc, lc, jc, tab, jab
    integer, intent(in) :: Ncm, Lcm, n12, l12, s12, j12, t12, n3, l3, j3, jr
    integer, intent(in) :: J
    real(8) :: r
    real(8) :: flab, flit, fllam, flam
    ! loop variables
    integer :: lab, lab_min, lab_max
    integer :: lit, lit_min, lit_max, nit
    integer :: LLam, LLam_min, LLam_max
    integer :: lam, lam_min, lam_max
    r = 0.d0
    if(tab /= t12) return
    if(2*(na+nb+nc)+la+lb+lc /= 2*(Ncm+n12+n3)+Lcm+l12+l3) return
    if(mod(la+lb+lc+Lcm+l12+l3,2)==1) return

    r = 0.d0
    lab_min = max(iabs(jab-s12), iabs(la-lb))
    lab_max = min(     jab+s12,       la+lb)
    do lab = lab_min, lab_max
      flab = dble(2*lab+1) * snj(2*la, 1, ja, 2*lb, 1, jb, 2*lab, 2*s12, 2*jab)
      if(abs(flab) < 1.d-16) cycle

      lit_min = max(iabs(lab-l12), iabs(jab-j12))
      lit_max = min(     lab+l12,       jab+j12, 2*na+la+2*nb+lb-2*n12-l12, 2*Ncm+lcm+2*n3+l3-2*nc-lc)
      do lit = lit_min, lit_max
        nit = (2*na + la + 2*nb + lb - 2*n12 - l12 - lit) / 2
        flit = gmosh(nit, lit, n12, l12, na, la, nb, lb, lab, 1.d0) * &
            & sjs(2*lit, 2*l12, 2*lab, 2*s12, 2*jab, 2*j12)
        if(abs(flit) < 1.d-16) cycle

        LLam_min = max(iabs(2*j12-J), iabs(2*lit-jc), iabs(2*Lcm-j3))
        LLam_max = min(     2*j12+J,       2*lit+jc,       2*Lcm+j3)
        do LLam = LLam_min, LLam_max, 2
          fLLam = (-1.d0)**( (LLam-Jr)/2 ) * &
              & dble(LLam+1) * &
              & sjs(2*j12, 2*lit, 2*jab, jc, J, LLam) * &
              & sjs(2*Lcm, j3, LLam, 2*j12, J, jr)
          if(abs(FLLam) < 1.d-16) cycle

          lam_min = max(iabs(lit-lc), iabs(1-LLam)/2, iabs(Lcm-l3))
          lam_max = min(     lit+lc,      (1+LLam)/2,     (Lcm+l3))
          do lam = lam_min, lam_max
            flam = dble(2*lam+1) * &
                & sjs(2*lit, 2*lc, 2*lam, 1, LLam, jc) * &
                & sjs(2*Lcm, 2*l3, 2*lam, 1, LLam, j3) * &
                & gmosh(ncm, lcm, n3, l3, nit, lit, nc, lc, lam, 2.d0)
            if(abs(flam) < 1.d-16) cycle
            r = r + flam * fLLam * flab * flit

          end do
        end do
      end do
    end do
    r = r * sqrt(dble((ja+1)*(jb+1)*(jc+1)*(2*jab+1)*(2*s12+1)*(2*j12+1)*(j3+1)*(jr+1))) * &
        & (-1.d0)**(s12+l12+1+j12+lc+l3+(jc+J)/2)
  end function three_body_tcoef_jj_iso

end module StoreCouplings
