#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION
!
! Three-body matrix elements in laboratory frame
! it is possible to calculate e3max = 40 or so, hopefully.
!

  public :: MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)
  public :: MACRO_ADD_SUFX(ThreeBodyMonopoleIso)
  private

  ! Methods for writing to file
  private :: WriteThreeBodyMonopoleIso
  !private :: ReadThreeBodyMonopole
  ! actual methods
  private :: write_three_body_monopole_iso_ascii
  private :: write_three_body_monopole_iso_ascii_me3j
  private :: write_three_body_monopole_iso_gzip_me3j
  !private :: write_scalar_operator_binary_stream_isospin
  !private :: write_scalar_operator_binary_isospin

  type :: MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)
#if defined(half_precision)
    type(my_real16), allocatable :: v(:)
#elif defined(single_precision)
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    type(MonThreeBodyLabIsoChan), pointer :: chbra, chket
  contains
    procedure :: InitThreeBodyMonopoleChIso
    procedure :: FinThreeBodyMonopoleChIso
    procedure :: GetMemoryCh
    procedure :: TMTransMonChannel
    generic :: init => InitThreeBodyMonopoleChIso
    generic :: fin => FinThreeBodyMonopoleChIso
    generic :: GetMemory => GetMemoryCh
    generic :: trans => TMTransMonChannel
  end type MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)

  type :: MACRO_ADD_SUFX(ThreeBodyMonopoleIso)
    type(MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)), allocatable :: MatCh(:)
    type(MonThreeBodyLabIsoSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitThreeBodyMonopoleIso
    procedure :: FinThreeBodyMonopoleIso
    procedure :: WriteThreeBodyMonopoleIso
    procedure :: GetFileName
    procedure :: GetThBMEMonopole
    procedure :: GetMemory
    generic :: init => InitThreeBodyMonopoleIso
    generic :: fin => FinThreeBodyMonopoleIso
    generic :: writef => WriteThreeBodyMonopoleIso
  end type MACRO_ADD_SUFX(ThreeBodyMonopoleIso)

contains

  function idx1d(bra,ket) result(n1d)
    integer, intent(in) :: bra, ket
    integer :: imax, imin
    integer(8) :: n1d
    imax = max(bra,ket)
    imin = min(bra,ket)
    n1d = (int(imax,kind(n1d))*int(imax-1,kind(n1d)))/int(2,kind(n1d)) + int(imin,kind(n1d))
  end function idx1d

  function GetMemory(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    real(8) :: memory
    integer :: ch
    memory = 0.d0
    do ch = 1, this%ms%GetNumberChannels()
      memory = memory + this%MatCh(ch)%GetMemory()
    end do
  end function GetMemory

  function GetMemoryCh(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)), intent(in) :: this
    integer(8) :: tmp
    real(8) :: memory
#if defined(half_precision)
    memory = 2.d0 * size( this%v, kind=kind(tmp) ) / 1024.d0**3
#elif defined(single_precision)
    memory = 4.d0 * size( this%v, kind=kind(tmp) ) / 1024.d0**3
#elif defined(double_precision)
    memory = 8.d0 * size( this%v, kind=kind(tmp) ) / 1024.d0**3
#endif
  end function GetMemoryCh

  function GetThBMEMonopole(this, i1, i2, i3, t12, i4, i5, i6, t45, T) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    integer, intent(in) :: i1, i2, i3, i4, i5, i6, t12, t45, T
    integer :: ch, ch1, ch2, ch3, bra, ket
    type(OrbitsIsospin), pointer :: isps
    type(MonThreeBodyLabIsoChan), pointer :: chan
    type(SingleParticleOrbitIsospin), pointer :: o1, o2, o3, o4, o5, o6
    integer :: ph
#if defined(single_precision) || defined(half_precision)
    real(4) :: r
    r = 0.0
#else
    real(8) :: r
    r = 0.d0
#endif

    isps => this%ms%sps
    o1 => isps%GetOrbit(i1)
    o2 => isps%GetOrbit(i2)
    o3 => isps%GetOrbit(i3)
    o4 => isps%GetOrbit(i4)
    o5 => isps%GetOrbit(i5)
    o6 => isps%GetOrbit(i6)
    if(o1%l /= o4%l) return
    if(o1%j /= o4%j) return
    if(o2%l /= o5%l) return
    if(o2%j /= o5%j) return
    if(o3%l /= o6%l) return
    if(o3%j /= o6%j) return

    ch1 = this%ms%one%GetIndex(o1%j, (-1)**o1%l)
    ch2 = this%ms%one%GetIndex(o2%j, (-1)**o2%l)
    ch3 = this%ms%one%GetIndex(o3%j, (-1)**o3%l)
    if(ch1 * ch2 * ch3 == 0) return
    ch = this%ms%TMon2Ch(T,ch1,ch2,ch3)
    if(ch == 0) return
    chan => this%ms%GetChannel(ch)
    if(this%ms%Flip(T,ch1,ch2,ch3)) then
      bra = chan%TNs2Index(t12,o2%n,o1%n,o3%n)
      ket = chan%TNs2Index(t45,o5%n,o4%n,o6%n)
      ph = (-1)**(t12+t45) * chan%TNs2Phase(t12,o2%n,o1%n,o3%n) * chan%TNs2Phase(t45,o5%n,o4%n,o6%n)
    else
      bra = chan%TNs2Index(t12,o1%n,o2%n,o3%n)
      ket = chan%TNs2Index(t45,o4%n,o5%n,o6%n)
      ph =  chan%TNs2Phase(t12,o1%n,o2%n,o3%n) * chan%TNs2Phase(t45,o4%n,o5%n,o6%n)
    end if
    if(bra * ket == 0) return
    r = this%MatCh(ch)%v( idx1d(bra,ket) ) * real(ph,kind=kind(r))
  end function GetThBMEMonopole

  subroutine TMTransMonChannel(this,jac2lab,vjac)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)), intent(inout), target :: this
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinMon)), intent(in), target :: jac2lab
    type(MonThreeBodyLabIsoChan), pointer :: lab
    type(ThreeBodyJacOpIso), intent(in), target :: vjac
    type(MACRO_ADD_SUFX(TCoefs)), pointer :: tcoef
    integer :: nl, ch, nj, trans_index, bra, ket, J
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v3(:,:)
    real(4), allocatable :: work(:,:)
    real(4) :: fac
#else
    real(8), allocatable :: v3(:,:)
    real(8), allocatable :: work(:,:)
    real(8) :: fac
#endif
#if defined(half_precision)
    real(4), allocatable :: t(:,:)
#endif
    real(8), pointer :: vv(:,:)
    real(8) :: ti
    integer(8) :: n1d
    type(sys) :: s
    lab => this%chket
    nl = lab%GetNumberStates()
    if(nl < 1) return
    allocate(v3(nl,nl))
#if defined(single_precision) || defined(half_precision)
    v3 = 0.0
#else
    v3 = 0.d0
#endif
    ti = omp_get_wtime()
    do trans_index = 1, jac2lab%GetNumberChannels()
      ch = jac2lab%GetJacobiIndex(trans_index)
      J = jac2lab%GetJ(trans_index)
      tcoef => jac2lab%Chan(trans_index)
      nj = tcoef%GetNumberASJacCMStates()
      allocate(work(nl,nj))
      vv => vjac%MatCh(ch,ch)%m(:nj,:nj)
#if defined(half_precision)
      allocate(t(nl,nj))
      t = tcoef%mat
      fac = 6.0 * real(J+1)
      call sgemm( 'n','n',nl,nj,nj,1.0,t,nl,real(vv),nj,0.0,work,nl )
      call sgemm( 'n','t',nl,nl,nj, fac,work,nl,t,nl,1.0,v3,nl )
      deallocate(t)
#elif defined(single_precision)
      fac = 6.0 * real(J+1)
      call sgemm( 'n','n',nl,nj,nj,1.0,tcoef%mat,nl,real(vv),nj,0.0,work,nl )
      call sgemm( 'n','t',nl,nl,nj, fac,work,nl,tcoef%mat,nl,1.0,v3,nl )
#else
      fac = 6.d0 * dble(J+1)
      call dgemm( 'n','n',nl,nj,nj,1.d0,tcoef%mat,nl,vv,nj,0.d0,work,nl )
      call dgemm( 'n','t',nl,nl,nj, fac,work,nl,tcoef%mat,nl,1.d0,v3,nl )
#endif
      deallocate(work)
    end do

    !$omp parallel
    !$omp do private(bra,ket,n1d)
    do bra = 1, lab%GetNumberStates()
      do ket = 1, bra
        n1d = idx1d(bra,ket)
        this%v(n1d) = v3(bra,ket)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(v3)
    call timer%Add(s%str('TMTransChannel'),omp_get_wtime()-ti)
  end subroutine TMTransMonChannel

  subroutine FinThreeBodyMonopoleChIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)), intent(inout) :: this
    deallocate(this%v)
    this%chbra => null()
    this%chket => null()
  end subroutine FinThreeBodyMonopoleChIso

  subroutine InitThreeBodyMonopoleChIso(this, chbra, chket)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleChIso)), intent(inout) :: this
    type(MonThreeBodyLabIsoChan), intent(in), target :: chbra, chket
    this%chbra => chbra
    this%chket => chket
    if( loc(chbra) /= loc(chket) ) then
      write(*,*) "Tensor is not implemented yet."
      return
    end if
    allocate( this%v( idx1d(chket%GetNumberStates(), chket%GetNumberStates()) ) )
#if defined(single_precision) || defined(half_precision)
    this%v(:) = 0.0
#else
    this%v(:) = 0.d0
#endif
  end subroutine InitThreeBodyMonopoleChIso

  subroutine FinThreeBodyMonopoleIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(inout) :: this
    integer :: ch
    do ch = 1, this%ms%GetNumberChannels()
      call this%MatCh(ch)%fin()
    end do
    deallocate(this%MatCh)
    this%ms => null()
  end subroutine FinThreeBodyMonopoleIso

  subroutine InitThreeBodyMonopoleIso(this, ms)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(inout) :: this
    type(MonThreeBodyLabIsoSpace), intent(in), target :: ms
    integer :: ch
    if(allocated(this%MatCh)) call this%fin()
    allocate(this%MatCh(ms%GetNumberChannels()))
    this%ms => ms
    do ch = 1, ms%GetNumberChannels()
      call this%MatCh(ch)%init( ms%GetChannel(ch), ms%GetChannel(ch) )
    end do
    this%is_init = .true.
  end subroutine InitThreeBodyMonopoleIso

  function GetFileName(this, filename, &
        &  NNInt, NNNInt, renorm, lambda, hw, emax, e2max, e3max, &
        &  cd, ce, genuine3bf, Rtype, Rpower) result(f)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: filename, NNInt, NNNInt, renorm, Rtype
    real(8), intent(in) :: lambda, hw, cd, ce
    integer, intent(in) :: emax, e2max, e3max, Rpower
    logical, intent(in) :: genuine3bf
    type(str) :: f
    type(sys) :: s
    logical :: dummy

    dummy = this%is_init
    if(filename%val /= "default") then
      f = filename
      return
    end if

    if(NNNInt%val /= 'ChEFT_N2LO') then
      write(*,'(a)') 'NNNInt other than ChEFT_N2LO has not been implemented'
      stop
    end if
    f = 'Monopole_ThBME'
    if(renorm%val == 'bare') then
      if(genuine3bf) then
        f = f + s%str('_') + NNNInt + s%str('_cD') + s%str(cd)
        f = f + s%str('cE') + s%str(ce) + s%str('_')
        f = f + Rtype + s%str(Rpower)
      else
        write(*,'(a)') 'Error in GetFileNameThreeBodyLabOp'
        write(*,'(a)') 'renorm = "bare" and genuine_3bf = False'
        stop
      end if
    else
        f = f + s%str('_') + renorm + s%str(lambda)
        f = f + s%str('_') + NNInt
      if(genuine3bf) then
        f = f + s%str('_') + NNNInt + s%str('_cD') + s%str(cd)
        f = f + s%str('cE') + s%str(ce) + s%str('_')
        f = f + Rtype + s%str(Rpower)
      end if
    end if
    f = f + s%str('_IS')
    f = f + s%str('_hw')  + s%str(hw)
    f = f + s%str('_ms')  + s%str(emax) + s%str('_') + s%str(e2max)
    f = f + s%str('_') + s%str(e3max) + s%str('.me3j.gz')
  end function GetFileName

  subroutine WriteThreeBodyMonopoleIso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti
    write(*,*) trim(f%val)
    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling WriteThreeBodyMonopole"
      return
    end if
    ti = omp_get_wtime()
    call write_three_body_monopole_iso(this, f)
    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyMonopoleIso

  subroutine write_three_body_monopole_iso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str('.txt'))) then
      call write_three_body_monopole_iso_ascii(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j.gz'))) then
      call write_three_body_monopole_iso_gzip_me3j(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j'))) then
      call write_three_body_monopole_iso_ascii_me3j(this,f)
      return
    end if

    write(*,*) "Unknown format, ascii is assumed: ", trim(f%val)
    call write_three_body_monopole_iso_ascii(this, f)
    return

  end subroutine write_three_body_monopole_iso

  subroutine write_three_body_monopole_iso_ascii(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(MonThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1, ch1
    integer :: i2, n2, l2, j2, ch2
    integer :: i3, n3, l3, j3, ch3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: t12, t45
    integer :: cntbf, total_num=0
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    header = 'NNN Monopole int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    write(wunit,'(a)') trim(header)
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      ch1 = lab%one%GetIndex(j1,(-1)**l1)
      do i2 = 1, sps%norbs
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        ch2 = lab%one%GetIndex(j2,(-1)**l2)
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          ch3 = lab%one%GetIndex(j3,(-1)**l3)
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            if(j1 /= j4) cycle
            if(l1 /= l4) cycle

            do i5 = 1, sps%norbs
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              if(j2 /= j5) cycle
              if(l2 /= l5) cycle
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, sps%norbs
                n6 = sps%orb(i6)%n
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do t12 = 0, 1
                  do t45 = 0, 1
                    do t123 = 1, min(2*t12+1,2*t45+1), 2

                      total_num = total_num + 1
                      write(wunit,'(9i4,f12.6)') i1, i2, i3, t12, i4, i5, i6, t45, T123, &
                          & this%GetThBMEMonopole(i1,i2,i3,t12,i4,i5,i6,t45,T123)
                    end do

                  end do
                end do

              end do

            end do

          end do
        end do
      end do
    end do

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    close(wunit)
  end subroutine write_three_body_monopole_iso_ascii

  subroutine write_three_body_monopole_iso_gzip_me3j(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(MonThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1, ch1
    integer :: i2, n2, l2, j2, ch2
    integer :: i3, n3, l3, j3, ch3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: t12, t45, nrest, cnt
    integer :: numbf = 10000000, cntbf, total_num=0
    integer :: wunit=25
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#else
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    type(c_ptr) :: fp, err
    character(256) :: header, buffer
    character(12) :: cfmt

    fp = gzip_open(f%val, "wt")
    header = 'NNN Monopole int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    err = gzip_writeline(fp, trim(header), len_trim(header))

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      ch1 = lab%one%GetIndex(j1,(-1)**l1)
      do i2 = 1, sps%norbs
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        ch2 = lab%one%GetIndex(j2,(-1)**l2)
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          ch3 = lab%one%GetIndex(j3,(-1)**l3)
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            if(j1 /= j4) cycle
            if(l1 /= l4) cycle

            do i5 = 1, sps%norbs
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              if(j2 /= j5) cycle
              if(l2 /= l5) cycle
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, sps%norbs
                n6 = sps%orb(i6)%n
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle


                do t12 = 0, 1
                  do t45 = 0, 1
                    do t123 = 1, min(2*t12+1,2*t45+1), 2

                      total_num = total_num + 1
                      cntbf = cntbf + 1
                      v123 = this%GetThBMEMonopole(i1,i2,i3,t12,i4,i5,i6,t45,t123)
                      v_buffer(cntbf) = v123
                      if(cntbf == numbf) then
                        do cnt = 1, cntbf/10
                          write(buffer,'(10f16.8)') v_buffer(10*(cnt-1)+1:10*cnt)
                          err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
                        end do
                        write(wunit,'(10f16.8)') v_buffer(:)
                        cntbf = 0
                      end if
                    end do

                  end do
                end do

              end do

            end do

          end do
        end do
      end do
    end do

    if(cntbf > 1) then
      do cnt = 1, cntbf/10
        write(buffer,'(10f16.8)') v_buffer(10*(cnt-1)+1:10*cnt)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end do

      nrest = cntbf - 10 * (cntbf/10)
      if(nrest > 1) then
        cfmt = '(xf16.8)'
        write(cfmt(2:2),'(i1)') nrest
        write(buffer,cfmt) v_buffer((cntbf/10)*10+1:cntbf)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end if
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    err = gzip_close(fp)
  end subroutine write_three_body_monopole_iso_gzip_me3j

  subroutine write_three_body_monopole_iso_ascii_me3j(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyMonopoleIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(MonThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: wunit=25
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1, ch1
    integer :: i2, n2, l2, j2, ch2
    integer :: i3, n3, l3, j3, ch3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: t12, t45
    integer :: numbf = 10000000, cntbf, total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#else
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    header = 'NNN Monopole int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    write(wunit,'(a)') trim(header)
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      ch1 = lab%one%GetIndex(j1,(-1)**l1)
      do i2 = 1, sps%norbs
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        ch2 = lab%one%GetIndex(j2,(-1)**l2)
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          ch3 = lab%one%GetIndex(j3,(-1)**l3)
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            if(j1 /= j4) cycle
            if(l1 /= l4) cycle

            do i5 = 1, sps%norbs
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              if(j2 /= j5) cycle
              if(l2 /= l5) cycle
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, sps%norbs
                n6 = sps%orb(i6)%n
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle


                do t12 = 0, 1
                  do t45 = 0, 1
                    do t123 = 1, min(2*t12+1,2*t45+1), 2

                      total_num = total_num + 1
                      cntbf = cntbf + 1
                      v123 = this%GetThBMEMonopole(i1,i2,i3,t12,i4,i5,i6,t45,t123)
                      v_buffer(cntbf) = v123
                      if(cntbf == numbf) then
                        write(wunit,'(10f16.8)') v_buffer(:)
                        cntbf = 0
                      end if
                    end do

                  end do
                end do

              end do

            end do

          end do
        end do
      end do
    end do

    if(cntbf > 1) then
      write(wunit,'(10f16.8)') v_buffer(:cntbf)
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_three_body_monopole_iso_ascii_me3j
