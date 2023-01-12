#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION
!
! Three-body matrix elements in laboratory frame
!
  public :: MACRO_ADD_SUFX(ThreeBodyNO2BChIso)
  public :: MACRO_ADD_SUFX(ThreeBodyNO2BIso)
  public :: MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)

  private
  private :: GetMemory
  private :: GetMemoryCh
  private :: FinThreeBodyNO2BIso
  private :: InitThreeBodyNO2BIso
  private :: AngleAverage
  private :: GetThBMENO2B
  private :: SetThBMENO2B
  private :: TMTransNO2BChannel
  private :: GetFileName
  private :: WriteThreeBodyNO2BIso
  private :: ReadThreeBodyNO2BIso
  private :: write_three_body_NO2B_iso
  private :: write_three_body_NO2B_iso_ascii
  private :: write_three_body_NO2B_iso_gzip_me3j
  private :: write_three_body_NO2B_iso_ascii_me3j
  private :: write_three_body_NO2B_binary_stream
  private :: read_three_body_NO2B_iso
  private :: read_three_body_NO2B_iso_ascii
  private :: read_three_body_NO2B_iso_gzip_me3j
  private :: read_three_body_NO2B_iso_ascii_me3j
  private :: read_three_body_NO2B_binary_stream

  type :: MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)
#if defined(half_precision)
    type(my_real16), allocatable :: v(:)
#elif defined(single_precision)
    real(4), allocatable :: v(:)
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif
    type(ThreeBodyLabIsoChanNO2B), pointer :: chbra, chket
  contains
    procedure :: InitThreeBodyLabOpIsoChan
    procedure :: FinThreeBodyLabOpIsoChan
    procedure :: TMTransNO2BChannel
    generic :: init => InitThreeBodyLabOpIsoChan
    generic :: fin => FinThreeBodyLabOpIsoChan
    generic :: trans => TMTransNO2BChannel
  end type MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)

  type :: MACRO_ADD_SUFX(ThreeBodyNO2BChIso)
#if defined(half_precision)
    type(my_real16), allocatable :: v(:)
#elif defined(single_precision)
    real(4), allocatable :: v(:)
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
  contains
    procedure :: InitThreeBodyNO2BChIso
    procedure :: FinThreeBodyNO2BChIso
    procedure :: AngleAverage
    procedure :: GetMemoryCh
    generic :: init => InitThreeBodyNO2BChIso
    generic :: fin => FinThreeBodyNO2BChIso
    generic :: GetMemory => GetMemoryCh
  end type MACRO_ADD_SUFX(ThreeBodyNO2BChIso)

  type :: MACRO_ADD_SUFX(ThreeBodyNO2BIso)
    type(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)), allocatable :: MatCh(:)
    type(NO2BThreeBodyIsoSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitThreeBodyNO2BIso
    procedure :: FinThreeBodyNO2BIso
    procedure :: WriteThreeBodyNO2BIso
    procedure :: ReadThreeBodyNO2BIso
    procedure :: GetThBMENO2B
    procedure :: SetThBMENO2B
    procedure :: GetFileName
    procedure :: GetMemory
    generic :: init => InitThreeBodyNO2BIso
    generic :: fin => FinThreeBodyNO2BIso
    generic :: writef => WriteThreeBodyNO2BIso
    generic :: readf => ReadThreeBodyNO2BIso
  end type MACRO_ADD_SUFX(ThreeBodyNO2BIso)
  integer, private :: buffer_length = 100000000

contains

  function GetMemory(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    real(8) :: memory
    integer :: ch
    memory = 0.d0
    do ch = 1, this%ms%GetNumberChannels()
      memory = memory + this%MatCh(ch)%GetMemory()
    end do
  end function GetMemory

  function GetMemoryCh(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)), intent(in) :: this
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

  function idx1d(bra,ket) result(n1d)
    integer, intent(in) :: bra, ket
    integer :: imax, imin
    integer(8) :: n1d
    imax = max(bra,ket)
    imin = min(bra,ket)
    n1d = (int(imax,kind(n1d))*int(imax-1,kind(n1d)))/int(2,kind(n1d)) + int(imin,kind(n1d))
  end function idx1d

  subroutine FinThreeBodyLabOpIsoChan(this)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)), intent(inout) :: this
    deallocate(this%v)
    this%chbra => null()
    this%chket => null()
  end subroutine FinThreeBodyLabOpIsoChan

  subroutine InitThreeBodyLabOpIsoChan(this, chbra, chket)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)), intent(inout) :: this
    type(ThreeBodyLabIsoChanNO2B), intent(in), target :: chbra, chket
    this%chbra => chbra
    this%chket => chket
    if( loc(chbra) /= loc(chket) ) then
      write(*,*) "Tensor is not implemented yet."
      return
    end if
    allocate( this%v( idx1d(chket%GetNumberStates(), chket%GetNumberStates()) ) )
#if defined(single_precision) || defined(half_precision)
    this%v(:) = 0.0
#elif defined(double_precision)
    this%v(:) = 0.d0
#endif
  end subroutine InitThreeBodyLabOpIsoChan

  subroutine FinThreeBodyNO2BChIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)), intent(inout) :: this
    deallocate(this%v)
    this%chbra => null()
    this%chket => null()
  end subroutine FinThreeBodyNO2BChIso

  subroutine InitThreeBodyNO2BChIso(this, chbra, chket)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)), intent(inout) :: this
    type(NO2BThreeBodyIsoChan), intent(in), target :: chbra, chket
    this%chbra => chbra
    this%chket => chket
    if( loc(chbra) /= loc(chket) ) then
      write(*,*) "Tensor is not implemented yet."
      return
    end if

    allocate( this%v( idx1d(chket%GetNumberStates(), chket%GetNumberStates()) ) )
#if defined(single_precision) || defined(half_precision)
    this%v(:) = 0.0
#elif defined(double_precision)
    this%v(:) = 0.d0
#endif
  end subroutine InitThreeBodyNO2BChIso

  subroutine FinThreeBodyNO2BIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    integer :: ch
    do ch = 1, this%ms%GetNumberChannels()
      deallocate(this%MatCh(ch)%v)
    end do
    deallocate(this%MatCh)
    this%ms => null()
  end subroutine FinThreeBodyNO2BIso

  subroutine InitThreeBodyNO2BIso(this, ms)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(NO2BThreeBodyIsoSpace), intent(in), target :: ms
    integer :: ch
    if(allocated(this%MatCh)) call this%fin()
    allocate(this%MatCh(ms%GetNumberChannels() ))
    this%ms => ms
    do ch = 1, ms%GetNumberChannels()
      call this%MatCh(ch)%init( ms%GetChannel(ch), ms%GetChannel(ch) )
    end do
    this%is_init = .true.
  end subroutine InitThreeBodyNO2BIso

  subroutine AngleAverage(this, v_ch)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BChIso)), intent(inout) :: this
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)), intent(in) :: v_ch
    type(NO2BThreeBodyIsoChan), pointer :: ch_averaged
    type(ThreeBodyLabIsoChanNO2B), pointer :: ch_original
    integer :: bra, ket, ibra, iket, jtot
    integer :: i1, i2, i3, i4, i5, i6, t12, t45
    real(8) :: phase
    integer :: QNsBra(4), QNsKet(4)
    integer(8) :: n1d
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: tmp(:,:)
#elif defined(double_precision)
    real(8), allocatable :: tmp(:,:)
#endif
    real(8) :: ti
    type(sys) :: s

    ti = omp_get_wtime()

    ch_original => v_ch%chket
    ch_averaged => this%chket
    Jtot = ch_original%GetJ()
    allocate(tmp(ch_averaged%GetNumberStates(), ch_averaged%GetNumberStates()))
#if defined(single_precision) || defined(half_precision)
    tmp(:,:) = 0.0
#elif defined(double_precision)
    tmp(:,:) = 0.d0
#endif

    !$omp parallel
    !$omp do private(bra,QNsBra,i1,i2,i3,t12,ket,QNsKet,i4,i5,i6,t45,&
    !$omp &  ibra,iket,phase)
    do bra = 1, ch_original%GetNumberStates()
      QNsBra = ch_original%GetABCT(bra)
      i1 = QNsBra(1)
      i2 = QNsBra(2)
      i3 = QNsBra(3)
      t12= QNsBra(4)
      do ket = 1, ch_original%GetNumberStates()
        QNsKet = ch_original%GetABCT(ket)
        i4 = QNsKet(1)
        i5 = QNsKet(2)
        i6 = QNsKet(3)
        t45= QNsKet(4)

        ibra = ch_averaged%GetIndex(i1, i2, i3, t12)
        iket = ch_averaged%GetIndex(i4, i5, i6, t45)
        phase = ch_averaged%GetPhase(i1,i2,i3,t12) * ch_averaged%GetPhase(i4,i5,i6,t45)
        if(ibra * iket == 0) cycle
#if defined(half_precision)
        tmp(ibra, iket) = v_ch%v(idx1d(bra,ket)) * real(phase) * real(jtot+1) + tmp(ibra,iket)
#elif defined(single_precision)
        tmp(ibra, iket) = tmp(ibra,iket) + real(jtot+1)  * v_ch%v(idx1d(bra,ket)) * real(phase)
#elif defined(double_precision)
        tmp(ibra, iket) = tmp(ibra, iket) + dble(jtot+1) * v_ch%v(idx1d(bra,ket)) * phase
#endif
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(ibra,iket,n1d)
    do ibra = 1, ch_averaged%GetNumberStates()
      do iket = 1, ibra
        n1d = idx1d(ibra,iket)
        this%v(n1d) = tmp(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(tmp)
    call timer%Add(s%str("Angle average"), omp_get_wtime() - ti)
  end subroutine AngleAverage

  subroutine TMTransNO2BChannel(this,jac2lab,vjac)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIsoChan)), intent(inout), target :: this
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in), target :: jac2lab
    type(ThreeBodyLabIsoChanNO2B), pointer :: lab
    type(ThreeBodyJacOpIso), intent(in), target :: vjac
    type(MACRO_ADD_SUFX(TCoefs)), pointer :: tcoef
    integer :: nl, ch, nj, trans_index, bra, ket
#if defined(half_precision)
    real(4), allocatable :: v3(:,:), work(:,:), t(:,:)
    real(4) :: fac
#elif defined(single_precision)
    real(4), allocatable :: v3(:,:), work(:,:)
    real(4) :: fac
#elif defined(double_precision)
    real(8), allocatable :: v3(:,:), work(:,:)
    real(8) :: fac
#endif
    real(8), pointer :: vv(:,:)
    real(8) :: ti
    integer(8) :: n1d
    type(sys) :: s
    lab => this%chket
    nl = lab%GetNumberStates()
    if(nl < 1) return
    allocate(v3(nl,nl))
    v3 = 0.d0
    ti = omp_get_wtime()
    do trans_index = 1, jac2lab%GetNumberChannels()
      ch = jac2lab%GetJacobiIndex(trans_index)
      tcoef => jac2lab%Chan(trans_index)
      nj = tcoef%GetNumberASJacCMStates()
      vv => vjac%MatCh(ch,ch)%m(:nj,:nj)
      allocate(work(nl,nj))
#if defined(half_precision)
      allocate(t(nl,nj))
      t = tcoef%mat
      fac = 6.0
      call sgemm( 'n','n',nl,nj,nj,1.0,t,nl,real(vv),nj,0.0,work,nl )
      call sgemm( 'n','t',nl,nl,nj, fac,work,nl,t,nl,1.0,v3,nl )
      deallocate(t)
#elif defined(single_precision)
      fac = 6.0
      call sgemm( 'n','n',nl,nj,nj,1.0,tcoef%mat,nl,real(vv),nj,0.0,work,nl )
      call sgemm( 'n','t',nl,nl,nj, fac,work,nl,tcoef%mat,nl,1.0,v3,nl )
#elif defined(double_precision)
      fac = 6.d0
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
  end subroutine TMTransNO2BChannel

  function GetFileName(this, filename, opname, &
        &  NNInt, NNNInt, renorm, lambda, hw, emax, e2max, e3max, &
        &  cd, ce, genuine3bf, Rtype, Rpower) result(f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: filename, opname, NNInt, NNNInt, renorm, Rtype
    real(8), intent(in) :: lambda, hw, cd, ce
    integer, intent(in) :: emax, e2max, e3max, Rpower
    logical, intent(in) :: genuine3bf
    type(str) :: f
    logical :: dummy
    type(sys) :: s

    dummy = this%is_init
    if(filename%val /= "default") then
      f = opname + s%str('_') + filename
      if(opname%val == 'hamil' .or. opname%val == 'NNNint') f = filename
      return
    end if

    if(NNNInt%val /= 'ChEFT_N2LO') then
      write(*,'(a)') 'NNNInt other than ChEFT_N2LO has not been implemented'
      stop
    end if
    f = opname + s%str('_NO2B_ThBME')
    if(opname%val == 'hamil' .or. opname%val == 'NNNint') f = s%str('NO2B_ThBME')
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

  function GetThBMENO2B(this,i1,i2,i3,T12,i4,i5,i6,T45,J,T) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: T12,T45,T,J
    type(OrbitsIsospin), pointer :: isps
    type(NO2BThreeBodyIsoChan), pointer :: chan
    integer :: ch, bra, ket
    real(8) :: phase
    integer :: P12, j3, p3
#if defined(single_precision) || defined(half_precision)
    real(4) :: r
    r = 0.0
#elif defined(double_precision)
    real(8) :: r
    r = 0.d0
#endif
    isps => this%ms%sps
    if(isps%orb(i3)%j /= isps%orb(i6)%j) return
    if(isps%orb(i3)%l /= isps%orb(i6)%l) return

    if((-1)**(isps%orb(i1)%l+ isps%orb(i2)%l+ &
        & isps%orb(i4)%l+ isps%orb(i5)%l) /= 1) then
      write(*,*) 'Warning: in GetThBMEIso_scalar: P'
      return
    end if

    p12 = (-1)**(isps%orb(i1)%l + isps%orb(i2)%l)
    j3 = isps%orb(i3)%j
    p3 = (-1)**(isps%orb(i3)%l)
    ch = this%ms%GetIndex(T,J,P12,j3,p3)
    if(ch == 0) return
    chan => this%ms%GetChannel(ch)
    bra = chan%TNcAB2Index(t12,isps%orb(i3)%n,i1,i2)
    ket = chan%TNcAB2Index(t45,isps%orb(i6)%n,i4,i5)
    if(bra * ket == 0) return
    phase = dble( chan%iphase(t12,isps%orb(i3)%n,i1,i2) * chan%iphase(t45,isps%orb(i6)%n,i4,i5) )
#if defined(single_precision) || defined(half_precision)
    r = this%MatCh(ch)%v(idx1d(bra,ket)) * real(phase)
#elif defined(double_precision)
    r = this%MatCh(ch)%v(idx1d(bra,ket)) * dble(phase)
#endif
  end function GetThBMENO2B

  subroutine SetThBMENO2B(this,i1,i2,i3,T12,i4,i5,i6,T45,J,T,me)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: T12,T45,T,J
#if defined(single_precision) || defined(half_precision)
    real(4) :: me
#elif defined(double_precision)
    real(8) :: me
#endif
    type(OrbitsIsospin), pointer :: isps
    type(NO2BThreeBodyIsoChan), pointer :: chan
    integer :: ch, bra, ket
    real(8) :: phase
    integer :: P12, j3, p3
    isps => this%ms%sps
    if(isps%orb(i3)%j /= isps%orb(i6)%j) return
    if(isps%orb(i3)%l /= isps%orb(i6)%l) return

    if((-1)**(isps%orb(i1)%l+ isps%orb(i2)%l + &
        & isps%orb(i4)%l+ isps%orb(i5)%l ) /= 1) then
      write(*,*) 'Warning: in SetThBMEIso_scalar: P'
      return
    end if

    p12 = (-1)**(isps%orb(i1)%l + isps%orb(i2)%l)
    j3 = isps%orb(i3)%j
    p3 = (-1)**(isps%orb(i3)%l)
    ch = this%ms%GetIndex(T,J,P12,j3,p3)
    if(ch == 0) return
    chan => this%ms%GetChannel(ch)
    bra = chan%TNcAB2Index(t12,isps%orb(i3)%n,i1,i2)
    ket = chan%TNcAB2Index(t45,isps%orb(i6)%n,i4,i5)
    if(bra * ket == 0) return
    phase = dble( chan%iphase(t12,isps%orb(i3)%n,i1,i2) * chan%iphase(t45,isps%orb(i6)%n,i4,i5) )
#if defined(single_precision) || defined(half_precision)
    this%MatCh(ch)%v(idx1d(bra,ket)) = real(phase) * me
#elif defined(double_precision)
    this%MatCh(ch)%v(idx1d(bra,ket)) = dble(phase) * me
#endif
  end subroutine SetThBMENO2B

  subroutine WriteThreeBodyNO2BIso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti
    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling WriteThreeBodyNO2B"
      return
    end if
    ti = omp_get_wtime()
    call write_three_body_NO2B_iso(this, f)
    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyNO2BIso

  subroutine ReadThreeBodyNO2BIso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti
    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling ReadThreeBodyNO2B"
      return
    end if
    ti = omp_get_wtime()
    call read_three_body_NO2B_iso(this, f)
    call timer%Add(s%str('Read from file'), omp_get_wtime() - ti)
  end subroutine ReadThreeBodyNO2BIso

  subroutine write_three_body_NO2B_iso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str('.txt'))) then
      call write_three_body_NO2B_iso_ascii(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j.gz'))) then
      call write_three_body_NO2B_iso_gzip_me3j(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j'))) then
      call write_three_body_NO2B_iso_ascii_me3j(this,f)
      return
    end if

    if(s%find(f,s%str('.stream.bin'))) then
      call write_three_body_NO2B_binary_stream(this, f)
      return
    end if

    write(*,'(a)') ' Unknown format, writting file assuming readable me3j'
    call write_three_body_NO2B_iso_ascii_me3j(this,f)
    return
  end subroutine write_three_body_NO2B_iso

  subroutine write_three_body_NO2B_iso_ascii(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf
    integer(8) :: total_num=0
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    header = 'NNN NO2B int.'
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
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        total_num = total_num + 1
                        write(wunit,'(10i4,f12.6)') i1, i2, i3, t12, i4, i5, i6, t45, J, T123, &
                            & this%GetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123)
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    close(wunit)
  end subroutine write_three_body_NO2B_iso_ascii

  subroutine write_three_body_NO2B_iso_gzip_me3j(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: cnt, nrest
    integer :: J, t12, t45
    integer :: cntbf
    integer(8) :: total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    type(c_ptr) :: fp, err
    character(256) :: header, buffer
    character(12) :: cfmt

    fp = gzip_open(f%val, "wt")
    header = 'NNN NO2B int.'
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
    allocate(v_buffer(buffer_length))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        total_num = total_num + 1
                        cntbf = cntbf + 1
                        v123 = this%GetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123)
                        v_buffer(cntbf) = v123
                        if(cntbf == buffer_length) then
                          do cnt = 1, cntbf/10
                            write(buffer,'(10f16.8)') v_buffer(10*(cnt-1)+1:10*cnt)
                            err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
                          end do
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
    end do

    if(cntbf > 1) then
      do cnt = 1, cntbf/10
        write(buffer,'(10f16.8)') v_buffer(10*(cnt-1)+1:10*cnt)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end do

      nrest = cntbf - 10 * (cntbf/10)
      if(nrest > 0) then
        cfmt = '(xf16.8)'
        write(cfmt(2:2),'(i1)') nrest
        write(buffer,cfmt) v_buffer((cntbf/10)*10+1:cntbf)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end if
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    err = gzip_close(fp)
  end subroutine write_three_body_NO2B_iso_gzip_me3j

  subroutine write_three_body_NO2B_iso_ascii_me3j(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: wunit=25
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf
    integer(8) :: total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    header = 'NNN NO2B int.'
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
    allocate(v_buffer(buffer_length))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        total_num = total_num + 1
                        cntbf = cntbf + 1
                        v123 = this%GetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123)
                        v_buffer(cntbf) = v123
                        if(cntbf == buffer_length) then
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
    end do

    if(cntbf > 0) then
      write(wunit,'(10f16.8)') v_buffer(:cntbf)
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_three_body_NO2B_iso_ascii_me3j

  subroutine write_three_body_NO2B_binary_stream(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf
    integer(8) :: total_num=0
#if defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(half_precision)
    integer(2), allocatable :: v_buffer(:)
    real(4) :: tmp
    type(my_real16) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif

    open(wunit, file=f%val, form='unformatted',access='stream',status='replace', action='write')
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(buffer_length))
    total_num = 0; cntbf = 0

    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2
                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        total_num = total_num + 1
                        cntbf = cntbf + 1
#ifdef half_precision
                        tmp = this%GetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,t123)
                        v123 = tmp
#else
                        v123 = this%GetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,t123)
#endif

                        v_buffer(cntbf) = v123
                        if(cntbf == buffer_length) then
                          write(wunit) v_buffer
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
    end do

    if(cntbf > 0) then
      write(wunit) v_buffer(1:cntbf)
    end if
    close(wunit)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
  end subroutine write_three_body_NO2B_binary_stream

  subroutine read_three_body_NO2B_iso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    integer(8) :: nelms

    nelms = count_three_body_no2b_me(this)

    if(s%find(f,s%str('.txt'))) then
      call read_three_body_NO2B_iso_ascii(this, f, nelms)
      return
    end if

    if(s%find(f,s%str('.me3j.gz'))) then
      call read_three_body_NO2B_iso_gzip_me3j(this, f, nelms)
      return
    end if

    if(s%find(f,s%str('.me3j'))) then
      call read_three_body_NO2B_iso_ascii_me3j(this,f, nelms)
      return
    end if

    if(s%find(f,s%str('.stream.bin'))) then
      call read_three_body_NO2B_binary_stream(this, f, nelms)
      return
    end if

    write(*,'(a)') ' Unknown format, writting file assuming readable me3j'
    call read_three_body_NO2B_iso_ascii_me3j(this,f,nelms)
    return
  end subroutine read_three_body_NO2B_iso

  subroutine read_three_body_NO2B_iso_ascii(this, f, nelms)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    integer(8), intent(in) :: nelms
    integer :: runit = 24
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf, n
    integer(8) :: total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v(:)
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif
    open(runit, file=f%val, status='old', action='read')
    read(runit,*)
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v(buffer_length))

    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        if(cntbf==buffer_length) cntbf = 0
                        if(cntbf == 0) then
#if defined(single_precision) || defined(half_precision)
                          v(:) = 0.0
#elif defined(double_precision)
                          v(:) = 0.d0
#endif
                          if(nelms - total_num >= buffer_length) then
                            do n = 1, buffer_length
                              read(runit,*) v( n )
                            end do
                          else
                            do n = 1, int((nelms - total_num)/10,kind(n))
                              read(runit,*) v( n )
                            end do
                          end if
                        end if
                        total_num = total_num + 1
                        cntbf = cntbf + 1

                        if(i1==i2 .and. mod(J+T12,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        if(i4==i5 .and. mod(J+T45,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        call this%SetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)
  end subroutine read_three_body_NO2B_iso_ascii

  subroutine read_three_body_NO2B_iso_gzip_me3j(this, f, nelms)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    integer(8), intent(in) :: nelms
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf, n
    integer(8) :: total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v(:)
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif
    type(c_ptr) :: fp, err
    character(256) :: header="", buffer=""

    fp = gzip_open(f%val, "rt")
    err = gzip_readline(fp, header, len(header))
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v(buffer_length))

    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        if(cntbf==buffer_length) cntbf = 0
                        if(cntbf == 0) then
#if defined(single_precision) || defined(half_precision)
                          v(:) = 0.0
#elif defined(double_precision)
                          v(:) = 0.d0
#endif
                          if(nelms - total_num >= buffer_length) then
                            do n = 1, buffer_length/10
                              err = gzip_readline(fp, buffer, len(buffer))
                              read(buffer,*) v( (n-1)*10+1 : n*10 )
                            end do
                          else
                            do n = 1, int((nelms - total_num)/10,kind(n))
                              err = gzip_readline(fp, buffer, len(buffer))
                              read(buffer,*) v( (n-1)*10+1 : n*10 )
                            end do
                            n = mod(int(nelms-total_num,kind(n)),10)
                            if(n > 0) then
                              buffer = ""
                              err = gzip_readline(fp, buffer, len(buffer))
                              read(buffer,*) v( nelms-n+1 : nelms )
                            end if
                          end if
                        end if
                        total_num = total_num + 1
                        cntbf = cntbf + 1

                        if(i1==i2 .and. mod(J+T12,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        if(i4==i5 .and. mod(J+T45,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        call this%SetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    err = gzip_close(fp)
  end subroutine read_three_body_NO2B_iso_gzip_me3j

  subroutine read_three_body_NO2B_iso_ascii_me3j(this, f, nelms)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    integer(8), intent(in) :: nelms
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: runit=24
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf, n
    integer(8) :: total_num=0
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: v(:)
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif

    open(runit, file=f%val, status='old', action='read')
    read(runit,*)
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v(buffer_length))

    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        if(cntbf==buffer_length) cntbf = 0
                        if(cntbf == 0) then
#if defined(single_precision) || defined(half_precision)
                          v(:) = 0.0
#elif defined(double_precision)
                          v(:) = 0.d0
#endif
                          if(nelms - total_num >= buffer_length) then
                            do n = 1, buffer_length/10
                              read(runit,*) v( (n-1)*10+1 : n*10 )
                            end do
                          else
                            do n = 1, int((nelms - total_num)/10,kind(n))
                              read(runit,*) v( (n-1)*10+1 : n*10 )
                              if(nelms-total_num-10*n - ((nelms-total_num-10*n)/10)*10 > 0) then
                                read(runit,*) v(((nelms-total_num-10*n)/10)*10+1 : nelms)
                              end if
                            end do
                          end if
                        end if
                        total_num = total_num + 1
                        cntbf = cntbf + 1

                        if(i1==i2 .and. mod(J+T12,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        if(i4==i5 .and. mod(J+T45,2)==0) then
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
                          cycle
                        end if

                        call this%SetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)
  end subroutine read_three_body_NO2B_iso_ascii_me3j

  subroutine read_three_body_NO2B_binary_stream(this, f, nelms)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    integer :: runit = 24
    integer(8), intent(in) :: nelms
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    integer :: cntbf
    integer(8) :: total_num=0
#if defined(single_precision)
    real(4), allocatable :: v(:)
#elif defined(half_precision)
    integer(2), allocatable :: v(:)
    type(my_real16) :: tmp
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v(:)
#endif

    open(runit, file=f%val, action='read', form='unformatted',access='stream')
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v(buffer_length))

    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
                        if(cntbf==buffer_length) cntbf = 0
                        if(cntbf == 0) then
#if defined(single_precision)
                          v(:) = 0.0
#elif defined(half_precision)
                          v(:) = 0
#elif defined(double_precision)
                          v(:) = 0.d0
#endif
                          if(nelms - total_num >= buffer_length) then
                            read(runit) v(:)
                          else
                            read(runit) v(:nelms-total_num)
                          end if
                        end if
                        total_num = total_num + 1
                        cntbf = cntbf + 1

                        if(i1==i2 .and. mod(J+T12,2)==0) then
#if defined(half_precision)
                          if(v(cntbf) /= 0 .and. v(cntbf) /= -32768) then
                            call show_warning_reading_me_int( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
#else
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
#endif
                          cycle
                        end if

                        if(i4==i5 .and. mod(J+T45,2)==0) then
#if defined(half_precision)
                          if(v(cntbf) /= 0 .and. v(cntbf) /= -32768) then
                            call show_warning_reading_me_int( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
#else
                          if(abs(v(cntbf)) > 1.d-6) then
                            call show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf) )
                          end if
#endif
                          cycle
                        end if

#if defined(half_precision)
                        tmp = v(cntbf)
                        v123 = tmp
                        call this%SetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v123)
#else
                        call this%SetThBMENO2B(i1,i2,i3,t12,i4,i5,i6,t45,J,T123,v(cntbf))
#endif
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)

  end subroutine read_three_body_NO2B_binary_stream

  function count_three_body_no2b_me(this) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyNO2BIso)), intent(in), target :: this
    integer(8) :: r
    type(NO2BThreeBodyIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, n1, l1, j1
    integer :: i2, n2, l2, j2
    integer :: i3, n3, l3, j3
    integer :: i4, n4, l4, j4
    integer :: i5, n5, l5, j5
    integer :: i6, n6, l6, j6
    integer :: p123, p456, t123
    integer :: J, t12, t45
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    r = 0
    do i1 = 1, sps%norbs
      n1 = sps%orb(i1)%n
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        n2 = sps%orb(i2)%n
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, sps%norbs
          n3 = sps%orb(i3)%n
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle
          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            n4 = sps%orb(i4)%n
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            do i5 = 1, i4
              n5 = sps%orb(i5)%n
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
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

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2),(j4+j5))/2

                  do t12 = 0, 1
                    do t45 = 0, 1
                      do t123 = 1, min(2*t12+1, 2*t45+1), 2
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
  end function count_three_body_no2b_me

  subroutine show_warning_reading_me( i1,i2,i3,t12,i4,i5,i6,t45,J,T,v )
    integer, intent(in) :: i1,i2,i3,i4,i5,i6,t12,t45,J,T
#if defined(single_precision) || defined(half_precision)
    real(4) :: v
#elif defined(double_precision)
    real(8) :: v
#endif
    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
    write(*,"(10i4,f12.6)") i1,i2,i3,T12,i4,i5,i6,T45,J,T,v
  end subroutine show_warning_reading_me

#ifdef half_precision
  subroutine show_warning_reading_me_int( i1,i2,i3,t12,i4,i5,i6,t45,J,T,v )
    integer, intent(in) :: i1,i2,i3,i4,i5,i6,t12,t45,J,T
    integer(2) :: v
    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
    write(*,"(10i4,i8)") i1,i2,i3,T12,i4,i5,i6,T45,J,T,v
  end subroutine show_warning_reading_me_int
#endif
