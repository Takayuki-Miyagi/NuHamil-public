#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION
!
! Three-body matrix elements in laboratory frame
!

  public :: MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)
  public :: MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)
  public :: MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)

  private
  private :: GetMemory
  private :: GetMemoryCh
  private :: FinThreeBodyTensorNO2BIso
  private :: InitThreeBodyTensorNO2BIso
  private :: AngleAverage
  private :: GetThBMENO2B
  private :: SetThBMENO2B
  private :: TMTransNO2BChannel
  private :: GetFileName
  private :: WriteThreeBodyTensorNO2BIso
  private :: ReadThreeBodyTensorNO2BIso
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

  type :: MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)
#if defined(half_precision)
    type(my_real16), allocatable :: m(:,:)
#elif defined(single_precision)
    real(4), allocatable :: m(:,:)
#elif defined(double_precision)
    real(8), allocatable :: m(:,:)
#endif
    type(ThreeBodyLabIsoChanNO2B), pointer :: chbra, chket
    logical :: is_zero=.true.
  contains
    procedure :: InitThreeBodyLabTensorIsoChan
    procedure :: FinThreeBodyLabTensorIsoChan
    procedure :: TMTransNO2BChannel
    generic :: init => InitThreeBodyLabTensorIsoChan
    generic :: fin => FinThreeBodyLabTensorIsoChan
    generic :: trans => TMTransNO2BChannel
  end type MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)

  type :: MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)
#if defined(half_precision)
    type(my_real16), allocatable :: m(:,:)
#elif defined(single_precision)
    real(4), allocatable :: m(:,:)
#elif defined(double_precision)
    real(8), allocatable :: m(:,:)
#endif
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
    logical :: is_zero=.true.
  contains
    procedure :: InitThreeBodyTensorNO2BChIso
    procedure :: FinThreeBodyTensorNO2BChIso
    procedure :: AngleAverage
    procedure :: GetMemoryCh
    generic :: init => InitThreeBodyTensorNO2BChIso
    generic :: fin => FinThreeBodyTensorNO2BChIso
    generic :: GetMemory => GetMemoryCh
  end type MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)

  type, extends(OperatorDef) :: MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)
    type(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)), allocatable :: MatCh(:,:)
    type(NO2BThreeBodyIsoSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitThreeBodyTensorNO2BIso
    procedure :: FinThreeBodyTensorNO2BIso
    procedure :: WriteThreeBodyTensorNO2BIso
    procedure :: ReadThreeBodyTensorNO2BIso
    procedure :: GetThBMENO2B
    procedure :: SetThBMENO2B
    procedure :: GetFileName
    procedure :: GetMemory
    generic :: init => InitThreeBodyTensorNO2BIso
    generic :: fin => FinThreeBodyTensorNO2BIso
    generic :: writef => WriteThreeBodyTensorNO2BIso
    generic :: readf => ReadThreeBodyTensorNO2BIso
  end type MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)
  integer, private :: buffer_length = 100000000

contains

  function GetMemory(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
    real(8) :: memory
    integer :: ichbra, ichket
    memory = 0.d0
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        memory = memory + this%MatCh(ichbra,ichket)%GetMemory()
      end do
    end do
  end function GetMemory

  function GetMemoryCh(this) result(memory)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)), intent(in) :: this
    integer(8) :: tmp
    real(8) :: memory
#if defined(half_precision)
    memory = 2.d0 * size( this%m, kind=kind(tmp) ) / 1024.d0**3
#elif defined(single_precision)
    memory = 4.d0 * size( this%m, kind=kind(tmp) ) / 1024.d0**3
#elif defined(double_precision)
    memory = 8.d0 * size( this%m, kind=kind(tmp) ) / 1024.d0**3
#endif
  end function GetMemoryCh

  subroutine FinThreeBodyLabTensorIsoChan(this)
    class(MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)), intent(inout) :: this
    deallocate(this%m)
    this%chbra => null()
    this%chket => null()
    this%is_zero=.true.
  end subroutine FinThreeBodyLabTensorIsoChan

  subroutine InitThreeBodyLabTensorIsoChan(this, chbra, chket)
    class(MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)), intent(inout) :: this
    type(ThreeBodyLabIsoChanNO2B), intent(in), target :: chbra, chket
    if( chbra%GetNumberStates() * chket%GetNumberStates() == 0 ) return
    this%chbra => chbra
    this%chket => chket
    allocate( this%m( chbra%GetNumberStates(), chket%GetNumberStates()) )
#if defined(single_precision) || defined(half_precision)
    this%m(:,:) = 0.0
#elif defined(double_precision)
    this%m(:,:) = 0.d0
#endif
    this%is_zero = .false.
  end subroutine InitThreeBodyLabTensorIsoChan

  subroutine FinThreeBodyTensorNO2BChIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)), intent(inout) :: this
    deallocate(this%m)
    this%chbra => null()
    this%chket => null()
    this%is_zero=.true.
  end subroutine FinThreeBodyTensorNO2BChIso

  subroutine InitThreeBodyTensorNO2BChIso(this, chbra, chket)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)), intent(inout) :: this
    type(NO2BThreeBodyIsoChan), intent(in), target :: chbra, chket
    if( chbra%GetNumberStates() * chket%GetNumberStates() == 0 ) return
    this%chbra => chbra
    this%chket => chket

    allocate( this%m( chbra%GetNumberStates(), chket%GetNumberStates()) )
#if defined(single_precision) || defined(half_precision)
    this%m(:,:) = 0.0
#elif defined(double_precision)
    this%m(:,:) = 0.d0
#endif
    this%is_zero = .false.
  end subroutine InitThreeBodyTensorNO2BChIso

  subroutine FinThreeBodyTensorNO2BIso(this)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
    integer :: ichbra, ichket
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        call this%MatCh(ichbra,ichket)%fin()
      end do
    end do
    deallocate(this%MatCh)
    this%ms => null()
  end subroutine FinThreeBodyTensorNO2BIso

  subroutine InitThreeBodyTensorNO2BIso(this, ms, OpName)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
    type(NO2BThreeBodyIsoSpace), intent(in), target :: ms
    type(str), intent(in) :: OpName
    integer :: ichbra, ichket
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
    if(allocated(this%MatCh)) call this%fin()
    call this%InitOpDef(OpName, .false.)
    allocate(this%MatCh(ms%GetNumberChannels(), ms%GetNumberChannels() ))
    this%ms => ms
    do ichbra = 1, ms%GetNumberChannels()
      chbra => ms%GetChannel(ichbra)
      do ichket = 1, ichbra
        chket => ms%GetChannel(ichket)
        if( triag( chbra%GetJ12(), chket%GetJ12(), this%GetOpJ() )) cycle
        if( triag( chbra%GetT(), chket%GetT(), 2*this%GetOpT() )) cycle
        if( chbra%GetP12() * chket%GetP12() * this%GetOpP() == -1 ) cycle
        call this%MatCh(ichbra,ichket)%init( ms%GetChannel(ichbra), ms%GetChannel(ichket) )
      end do
    end do
    this%is_init = .true.
  end subroutine InitThreeBodyTensorNO2BIso

  subroutine AngleAverage(this, v_ch, op_def)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BChIso)), intent(inout) :: this
    type(MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)), intent(in) :: v_ch
    type(OperatorDef), intent(in) :: op_def
    type(NO2BThreeBodyIsoChan), pointer :: chbra_averaged, chket_averaged
    type(ThreeBodyLabIsoChanNO2B), pointer :: chbra_original, chket_original
    integer :: bra, ket, ibra, iket, Jbra, Jket, J12bra, J12ket, Jh
    integer :: i1, i2, i3, i4, i5, i6, t12, t45
    real(8) :: phase
    integer :: QNsBra(4), QNsKet(4)
    integer(8) :: n1d
#if defined(single_precision) || defined(half_precision)
    real(4), allocatable :: tmp(:,:)
    real(4) :: factor
#elif defined(double_precision)
    real(8), allocatable :: tmp(:,:)
    real(4) :: factor
#endif
    real(8) :: ti
    type(sys) :: s

    ti = omp_get_wtime()

    chbra_original => v_ch%chbra
    chket_original => v_ch%chket

    chbra_averaged => this%chbra
    chket_averaged => this%chket

    Jbra = chbra_original%GetJ()
    Jket = chket_original%GetJ()
    J12bra = chbra_averaged%GetJ12()
    J12ket = chket_averaged%GetJ12()
    Jh = chket_averaged%GetJ3()
    factor = real( sqrt( dble(Jbra+1)*dble(Jket+1))*&
        & (-1.d0)**(J12bra+op_def%GetOpJ()+(Jh+Jket)/2 ) * &
        & sjs(2*J12bra, 2*J12ket, 2*op_def%GetOpJ(), Jket, Jbra, Jh), kind=kind(factor))
    if( op_def%GetOpJ()==0 .and. op_def%GetOpP()==1 .and. op_def%GetOpT()==0 ) then
      factor = real( sqrt(dble(Jbra+1)*dble(Jket+1)), kind=kind(factor))
    end if
    allocate(tmp(chbra_averaged%GetNumberStates(), chket_averaged%GetNumberStates()))
#if defined(single_precision) || defined(half_precision)
    tmp(:,:) = 0.0
#elif defined(double_precision)
    tmp(:,:) = 0.d0
#endif

    !$omp parallel
    !$omp do private(bra,QNsBra,i1,i2,i3,t12,ket,QNsKet,i4,i5,i6,t45,&
    !$omp &  ibra,iket,phase)
    do bra = 1, chbra_original%GetNumberStates()
      QNsBra = chbra_original%GetABCT(bra)
      i1 = QNsBra(1)
      i2 = QNsBra(2)
      i3 = QNsBra(3)
      t12= QNsBra(4)
      do ket = 1, chket_original%GetNumberStates()
        QNsKet = chket_original%GetABCT(ket)
        i4 = QNsKet(1)
        i5 = QNsKet(2)
        i6 = QNsKet(3)
        t45= QNsKet(4)

        ibra = chbra_averaged%GetIndex(i1, i2, i3, t12)
        iket = chket_averaged%GetIndex(i4, i5, i6, t45)
        phase = chbra_averaged%GetPhase(i1,i2,i3,t12) * chket_averaged%GetPhase(i4,i5,i6,t45)
        if(ibra * iket == 0) cycle
#if defined(single_precision)
        tmp(ibra, iket) = tmp(ibra,iket) + factor * v_ch%m(bra,ket) * real(phase)
#elif defined(half_precision)
        tmp(ibra, iket) = v_ch%m(bra,ket) * real(phase) * factor + tmp(ibra,iket)
#elif defined(double_precision)
        tmp(ibra, iket) = tmp(ibra, iket) + factor * v_ch%m(bra,ket) * phase
#endif
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(ibra,iket,n1d)
    do ibra = 1, chbra_averaged%GetNumberStates()
      do iket = 1, chket_averaged%GetNumberStates()
        this%m(ibra,iket) = tmp(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(tmp)
    call timer%Add(s%str("Angle average"), omp_get_wtime() - ti)
  end subroutine AngleAverage

  subroutine TMTransNO2BChannel(this,tr_bra,tr_ket,vjac)
    class(MACRO_ADD_SUFX(ThreeBodyLabTensorIsoChan)), intent(inout), target :: this
    type(MACRO_ADD_SUFX(TransJac2LabChanIsospinNO2B)), intent(in), target :: tr_bra, tr_ket
    type(ThreeBodyLabIsoChanNO2B), pointer :: lab_bra, lab_ket
    type(ThreeBodyJacOpIso), intent(in), target :: vjac
    type(MACRO_ADD_SUFX(TCoefs)), pointer :: tcoef_bra, tcoef_ket
    integer :: nlbra, nlket, trans_index_bra, trans_index_ket, bra, ket
    integer :: jchbra, jchket, njbra, njket, Lcm, jjbra, jjket, Jbra, Jket
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
    lab_bra => this%chbra
    lab_ket => this%chket
    nlbra = lab_bra%GetNumberStates()
    nlket = lab_ket%GetNumberStates()
    Jbra = lab_bra%GetJ()
    Jket = lab_ket%GetJ()
    if(nlbra*nlket < 1) return
    allocate(v3(nlbra,nlket))
    v3 = 0.d0
    ti = omp_get_wtime()
    do trans_index_bra = 1, tr_bra%GetNumberChannels()
      jchbra = tr_bra%GetJacobiIndex(trans_index_bra)
      do trans_index_ket = 1, tr_ket%GetNumberChannels()
        jchket = tr_ket%GetJacobiIndex(trans_index_ket)
        if( tr_bra%GetNcm( trans_index_bra ) /= tr_ket%GetNcm( trans_index_ket ) ) cycle
        if( tr_bra%GetLcm( trans_index_bra ) /= tr_ket%GetLcm( trans_index_ket ) ) cycle
        if( vjac%MatCh(jchbra,jchket)%is_zero ) cycle
        Lcm = tr_ket%GetLcm( trans_index_ket )
        jjbra = vjac%MatCh(jchbra,jchket)%jacobi_ch_bra%GetJ()
        jjket = vjac%MatCh(jchbra,jchket)%jacobi_ch_ket%GetJ()
        tcoef_bra => tr_bra%Chan(trans_index_bra)
        tcoef_ket => tr_ket%Chan(trans_index_ket)
        njbra = tcoef_bra%GetNumberASJacCMStates()
        njket = tcoef_ket%GetNumberASJacCMStates()
        vv => vjac%MatCh(jchbra,jchket)%m(:njbra,:njket)
        allocate(work( nlbra,njket ))
#if defined(half_precision)
        fac = 6.0 * real( (-1.d0)**( (2*Lcm+2*vjac%GetOpJ()+jjket+Jbra)/2 ) * &
            & sqrt( dble(Jbra+1) * dble(Jket+1) ) * &
            & sjs(jjbra, jjket, 2*vjac%GetOpJ(), Jket, Jket, 2*Lcm) )
        allocate(t(nlbra,njbra))
        t = tcoef_bra%mat
        call sgemm( 'n','n',nlbra,njket,njbra,1.0,t,nlbra,real(vv),njbra,0.0,work,nlbra )
        deallocate(t)

        allocate(t(nlket,njket))
        t = tcoef_ket%mat
        call sgemm( 'n','t',nlbra,nlket,njket, fac,work,nlbra,t,nlket,1.0,v3,nlbra )
        deallocate(t)
#elif defined(single_precision)
        fac = real(6.d0 * (-1.d0)**( (2*Lcm+2*vjac%GetOpJ()+jjket+Jbra)/2 ) * &
            & sqrt( dble(Jbra+1) * dble(Jket+1) ) * &
            & sjs(jjbra, jjket, 2*vjac%GetOpJ(), Jket, Jket, 2*Lcm))
        call sgemm( 'n','n',nlbra,njket,njbra,1.0,tcoef_bra%mat,nlbra,real(vv),njbra,0.0,work,nlbra )
        call sgemm( 'n','t',nlbra,nlket,njket, fac,work,nlbra,tcoef_ket%mat,nlket,1.0,v3,nlbra )
#elif defined(double_precision)
        fac = 6.d0 * (-1.d0)**( (2*Lcm+2*vjac%GetOpJ()+jjket+Jbra)/2 ) * &
            & sqrt( dble(Jbra+1) * dble(Jket+1) ) * &
            & sjs(jjbra, jjket, 2*vjac%GetOpJ(), Jket, Jket, 2*Lcm)
        call dgemm( 'n','n',nlbra,njket,njbra,1.d0,tcoef_bra%mat,nlbra,vv,njbra,0.d0,work,nlbra )
        call dgemm( 'n','t',nlbra,nlket,njket, fac,work,nlbra,tcoef_ket%mat,nlket,1.d0,v3,nlbra )
#endif
        deallocate(work)
      end do
    end do

    !$omp parallel
    !$omp do private(bra,ket,n1d)
    do bra = 1, lab_bra%GetNumberStates()
      do ket = 1, lab_ket%GetNumberStates()
#if defined(single_precision) || defined(half_precision)
        this%m(bra,ket) = real(v3(bra,ket))
#elif defined(double_precision)
        this%m(bra,ket) = v3(bra,ket)
#endif
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
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
    type(str), intent(in) :: filename, opname, NNInt, NNNInt, renorm, Rtype
    real(8), intent(in) :: lambda, hw, cd, ce
    integer, intent(in) :: emax, e2max, e3max, Rpower
    logical, intent(in) :: genuine3bf
    type(str) :: f
    type(sys) :: s
    logical :: dummy

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
        f = f + s%str('_') + renorm  + s%str(lambda)
        f = f + s%str('_') + NNInt
      if(genuine3bf) then
        f = f + s%str('_') + NNNInt + s%str('_cD') + s%str(cd)
        f = f + s%str('cE') + s%str(ce) + s%str('_')
        f = f + Rtype + s%str(Rpower)
      end if
    end if
    f = f + s%str('_IS')
    f = f + s%str('_hw') + s%str(hw)
    f = f + s%str('_ms') + s%str(emax) + s%str('_')  + s%str(e2max)
    f = f + s%str('_') + s%str(e3max) + s%str('.me3j.gz')
  end function GetFileName

  function GetThBMENO2B(this,i1,i2,i3,J12,T12,T123,i4,i5,i6,J45,T45,T456) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,T123,T456
    type(OrbitsIsospin), pointer :: isps
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
    integer :: ichbra, ichket, bra, ket
    real(8) :: phase
    integer :: P12, P45, j3, p3
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od, oe, of
#if defined(single_precision) || defined(half_precision)
    real(4) :: r
    r = 0.0
#elif defined(double_precision)
    real(8) :: r
    r = 0.d0
#endif
    if( triag( J12, J45, this%GetOpJ() )) return
    if( triag( T123, T456, 2*this%GetOpT() )) return
    isps => this%ms%sps
    oa => isps%GetOrbit(i1)
    ob => isps%GetOrbit(i2)
    oc => isps%GetOrbit(i3)
    od => isps%GetOrbit(i4)
    oe => isps%GetOrbit(i5)
    of => isps%GetOrbit(i6)
    if(oc%j /= of%j) return
    if(oc%l /= of%l) return
    P12 = (-1)**(oa%l+ob%l)
    P45 = (-1)**(od%l+oe%l)

    if( P12 * P45 * this%GetOpP() /= 1) then
      write(*,*) 'Warning: ', __LINE__, " in ", __FILE__
      return
    end if

    j3 = oc%j
    p3 = (-1)**(oc%l)
    ichbra = this%ms%GetIndex(T123,J12,P12,j3,p3)
    ichket = this%ms%GetIndex(T456,J45,P45,j3,p3)
    if(ichbra*ichket == 0) return
    chbra => this%ms%GetChannel(ichbra)
    chket => this%ms%GetChannel(ichket)
    bra = chbra%TNcAB2Index(t12,oc%n,i1,i2)
    ket = chket%TNcAB2Index(t45,of%n,i4,i5)
    if(bra * ket == 0) return
    phase = dble( chbra%iphase(t12,oc%n,i1,i2) * chket%iphase(t45,of%n,i4,i5) )
    if( ichbra < ichket ) then
#if defined(single_precision) || defined(half_precision)
      r = this%MatCh(ichket,ichbra)%m(ket,bra) * real(phase)
#elif defined(double_precision)
      r = this%MatCh(ichket,ichbra)%m(ket,bra) * dble(phase)
#endif
      return
    end if
#if defined(single_precision) || defined(half_precision)
    r = this%MatCh(ichbra,ichket)%m(bra,ket) * real(phase)
#elif defined(double_precision)
    r = this%MatCh(ichbra,ichket)%m(bra,ket) * dble(phase)
#endif
  end function GetThBMENO2B

  subroutine SetThBMENO2B(this,i1,i2,i3,J12,T12,T123,i4,i5,i6,J45,T45,T456,me)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    integer, intent(in) :: J12,T12,J45,T45,T123,T456
#if defined(single_precision) || defined(half_precision)
    real(4) :: me
#elif defined(double_precision)
    real(8) :: me
#endif
    type(OrbitsIsospin), pointer :: isps
    type(NO2BThreeBodyIsoChan), pointer :: chbra, chket
    integer :: ichbra, ichket, bra, ket
    real(8) :: phase
    integer :: P12, P45, j3, p3
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od, oe, of

    if( triag( J12, J45, this%GetOpJ() )) return
    if( triag( T123, T456, 2*this%GetOpT() )) return
    isps => this%ms%sps
    oa => isps%GetOrbit(i1)
    ob => isps%GetOrbit(i2)
    oc => isps%GetOrbit(i3)
    od => isps%GetOrbit(i4)
    oe => isps%GetOrbit(i5)
    of => isps%GetOrbit(i6)
    if(oc%j /= of%j) return
    if(oc%l /= of%l) return
    P12 = (-1)**(oa%l+ob%l)
    P45 = (-1)**(od%l+oe%l)

    if( P12 * P45 * this%GetOpP() /= 1) then
      write(*,*) 'Warning: ', __LINE__, " in ", __FILE__
      return
    end if

    j3 = oc%j
    p3 = (-1)**(oc%l)
    ichbra = this%ms%GetIndex(T123,J12,P12,j3,p3)
    ichket = this%ms%GetIndex(T456,J45,P45,j3,p3)
    if(ichbra*ichket == 0) return
    chbra => this%ms%GetChannel(ichbra)
    chket => this%ms%GetChannel(ichket)
    bra = chbra%TNcAB2Index(t12,oc%n,i1,i2)
    ket = chket%TNcAB2Index(t45,of%n,i4,i5)
    if(bra * ket == 0) return
    phase = dble( chbra%iphase(t12,oc%n,i1,i2) * chket%iphase(t45,of%n,i4,i5) )
    if( ichbra < ichket ) then
#if defined(single_precision) || defined(half_precision)
      this%MatCh(ichket,ichbra)%m(ket,bra) = me * real(phase)
#elif defined(double_precision)
      this%MatCh(ichket,ichbra)%m(ket,bra) = me * dble(phase)
#endif
      return
    end if
#if defined(single_precision) || defined(half_precision)
    this%MatCh(ichbra,ichket)%m(bra,ket) = me * real(phase)
#elif defined(double_precision)
    this%MatCh(ichbra,ichket)%m(bra,ket) = me * dble(phase)
#endif
    if(ichbra==ichket) then
#if defined(single_precision) || defined(half_precision)
      this%MatCh(ichbra,ichket)%m(ket,bra) = me * real(phase)
#elif defined(double_precision)
      this%MatCh(ichbra,ichket)%m(ket,bra) = me * dble(phase)
#endif
    end if
  end subroutine SetThBMENO2B

  subroutine WriteThreeBodyTensorNO2BIso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti
    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling WriteThreeBodyTensorNO2B"
      return
    end if
    ti = omp_get_wtime()
    call write_three_body_NO2B_iso(this, f)
    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyTensorNO2BIso

  subroutine ReadThreeBodyTensorNO2BIso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti
    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling ReadThreeBodyTensorNO2B"
      return
    end if
    ti = omp_get_wtime()
    call read_three_body_NO2B_iso(this, f)
    call timer%Add(s%str('Read from file'), omp_get_wtime() - ti)
  end subroutine ReadThreeBodyTensorNO2BIso

  subroutine write_three_body_NO2B_iso(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
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

    !if(s%find(f,'.bin')) then
    !  call write_three_body_NO2B_binary(this, f)
    !  return
    !end if

    write(*,'(a)') ' Unknown format, writting file assuming readable me3j'
    call write_three_body_NO2B_iso_ascii_me3j(this,f)
    return
  end subroutine write_three_body_NO2B_iso

  subroutine write_three_body_NO2B_iso_ascii(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
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
    integer :: p123, p456, t123, t456
    integer :: cnt
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() == -1 ) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBMENO2B(i1,i2,i3,j12,t12,T123,i4,i5,i6,j45,t45,T456)
                            v_buffer(cntbf) = v123
                            if(cntbf == buffer_length) then
                              do cnt = 1, cntbf
                                write(wunit,*) v_buffer(cnt)
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
      end do
    end do

    if(cntbf > 0) then
      do cnt = 1, cntbf
        write(wunit,*) v_buffer(cnt)
      end do
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_three_body_NO2B_iso_ascii

  subroutine write_three_body_NO2B_iso_gzip_me3j(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
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
    integer :: p123, p456, t123, t456
    integer :: cnt, nrest
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() == -1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBMENO2B(i1,i2,i3,j12,t12,T123,i4,i5,i6,j45,t45,T456)
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
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
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
    integer :: p123, p456, t123, t456
    integer :: j12, j45, t12, t45
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
                if( p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBMENO2B(i1,i2,i3,j12,t12,T123,i4,i5,i6,j45,t45,T456)
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
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in) :: this
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
    integer :: p123, p456, t123, t456
    integer :: j12, j45, t12, t45
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
                if(p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBMENO2B(i1,i2,i3,j12,t12,T123,i4,i5,i6,j45,t45,T456)
#ifdef half_precision
                            v123 = tmp
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
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
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
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
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
    integer :: p123, p456, t123, t456
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
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

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            call this%SetThBMENO2B(i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)
  end subroutine read_three_body_NO2B_iso_ascii

  subroutine read_three_body_NO2B_iso_gzip_me3j(this, f, nelms)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
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
    integer :: p123, p456, t123, t456
    integer :: j12, j45, t12, t45
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
                if(p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
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

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            call this%SetThBMENO2B(i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    err = gzip_close(fp)
  end subroutine read_three_body_NO2B_iso_gzip_me3j

  subroutine read_three_body_NO2B_iso_ascii_me3j(this, f, nelms)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
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
    integer :: p123, p456, t123, t456
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() == -1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
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

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf) )
                              end if
                              cycle
                            end if

                            call this%SetThBMENO2B(i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v(cntbf))
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

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)
  end subroutine read_three_body_NO2B_iso_ascii_me3j

  subroutine read_three_body_NO2B_binary_stream(this, f, nelms)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(inout) :: this
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
    integer :: p123, p456, t123, t456
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
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

                            if(i1==i2 .and. mod(J12+T12,2)==0) then
#if defined(half_precision)
                              if(v(cntbf) /= 0 .and. v(cntbf) /= -32768) then
                                call show_warning_reading_me_int( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v(cntbf) )
                              end if
#elif defined(single_precision) || defined(double_precision)
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v(cntbf) )
                              end if
#endif
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
#ifdef half_precision
                              if(v(cntbf) /= 0 .and. v(cntbf) /= -32768) then
                                call show_warning_reading_me_int( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v(cntbf) )
                              end if
#elif defined(single_precision) || defined(double_precision)
                              if(abs(v(cntbf)) > 1.d-6) then
                                call show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v(cntbf) )
                              end if
#endif
                              cycle
                            end if

#ifdef half_precision
                            tmp = v(cntbf)
                            v123 = tmp
                            call this%SetThBMENO2B(i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v123)
#elif defined(single_precision) || defined(double_precision)
                            call this%SetThBMENO2B(i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,t456,v(cntbf))
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
      end do
    end do

    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v)
    close(runit)

  end subroutine read_three_body_NO2B_binary_stream

  function count_three_body_no2b_me(this) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyTensorNO2BIso)), intent(in), target :: this
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
    integer :: p123, p456, t123, t456
    integer :: J12, J45, t12, t45
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
                if(p123 * p456 * this%GetOpP() ==-1) cycle

                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2

                    do t12 = 0, 1
                      do t45 = 0, 1
                        do t123 = 1, 2*t12+1, 2
                          do t456 = 1, 2*t45+1, 2
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
  end function count_three_body_no2b_me

  subroutine show_warning_reading_me( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v )
    integer, intent(in) :: i1,i2,i3,i4,i5,i6,j12,t12,t123,j45,t45,t456
#if defined(single_precision) || defined(half_precision)
    real(4) :: v
#elif defined(double_precision)
    real(8) :: v
#endif
    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
    write(*,"(12i4,f12.6)") i1,i2,i3,j12,T12,t123,i4,i5,i6,j45,T45,t456,v
  end subroutine show_warning_reading_me

#ifdef half_precision
  subroutine show_warning_reading_me_int( i1,i2,i3,j12,t12,t123,i4,i5,i6,j45,t45,T456,v )
    integer, intent(in) :: i1,i2,i3,i4,i5,i6,j12,t12,t123,j45,t45,t456
    integer(2) :: v
    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
    write(*,"(12i4,i8)") i1,i2,i3,j12,T12,t123,i4,i5,i6,j45,T45,t456,v
  end subroutine show_warning_reading_me_int
#endif

