#define MACRO_PASTE(A) A
#define MACRO_ADD_SUFX(STR) MACRO_PASTE(STR)PRECISION

  public :: MACRO_ADD_SUFX(ThreeBodyLabOpIso)
  private
  private :: FinThreeBodyLabOp ! destructor
  private :: InitThreeBodyLabOp ! constructor
  private :: InitThreeBodyLabOpFromString ! constructor (more convenient)
  private :: ReducedToNonReduced
  private :: NonReducedToReduced
  private :: GetThBME_scalar
  private :: GetThBME_tensor
  private :: SetThBME_scalar
  private :: SetThBME_tensor

  ! Methods for TalmiMoshinskyTransformation, TMTrans is interface
  private :: TalmiMoshinskyTransformation
  ! actual methods
  private :: TMTransScalar
  private :: TMTransTensor

  ! Methods for FileName
  private :: GetFileNameThreeBodyLabOp

  ! Methods for writing to file
  private :: WriteThreeBodyLabOp
  private :: ReadThreeBodyLabOp
  ! interfaces
  private :: write_scalar_operator_isospin
  private :: write_tensor_operator_isospin
  ! actual methods
  private :: write_scalar_operator_ascii_isospin
  private :: write_scalar_operator_ascii_me3j_isospin
  private :: write_scalar_operator_gzip_me3j_isospin
  private :: write_scalar_operator_binary_stream_isospin
  private :: write_scalar_operator_binary_isospin
  private :: write_tensor_operator_binary_stream_isospin
  private :: write_tensor_operator_gzip_isospin

  type, private :: ThreeBodyLabOpChanIso
#if defined(half_precision)
    type(my_real16), allocatable :: m(:,:)
#elif defined(single_precision)
    real(4), allocatable :: m(:,:)
#elif defined(double_precision)
    real(8), allocatable :: m(:,:)
#endif
    logical :: is_Zero = .true.
  end type ThreeBodyLabOpChanIso

  type, extends(OperatorDef) :: MACRO_ADD_SUFX(ThreeBodyLabOpIso)
    type(ThreeBodyLabOpChanIso), allocatable :: MatCh(:,:)
    type(ThreeBodyLabIsoSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitThreeBodyLabOpFromString
    procedure :: InitThreeBodyLabOp
    procedure :: FinThreeBodyLabOp
    procedure :: GetFileNameThreeBodyLabOp
    procedure :: TalmiMoshinskyTransformation
    procedure :: WriteThreeBodyLabOp
    procedure :: ReadThreeBodyLabOp
    procedure :: Truncate
    procedure :: ReducedToNonReduced
    procedure :: NonReducedToReduced
    procedure :: GetThBME_scalar
    procedure :: GetThBME_tensor
    procedure :: SetThBME_scalar
    procedure :: SetThBME_tensor
    generic :: init => InitThreeBodyLabOpFromString, InitThreeBodyLabOp
    generic :: fin => FinThreeBodyLabOp
    generic :: GetFile => GetFileNameThreeBodyLabOp
    generic :: TMTrans => TalmiMoshinskyTransformation
    generic :: writef => WriteThreeBodyLabOp
    generic :: readf => ReadThreeBodyLabOp
    generic :: GetThBME => GetThBME_scalar, GetThBME_tensor
    generic :: SetThBME => SetThBME_scalar, SetThBME_tensor

    procedure :: CopyThreeBodyLabOp
    procedure :: SumThreeBodyLabOp
    procedure :: SubtractThreeBodyLabOp
    procedure :: ScaleThreeBodyLabOp
    generic :: assignment(=) => CopyThreeBodyLabOp
    generic :: operator(+) => SumThreeBodyLabOp
    generic :: operator(-) => SubtractThreeBodyLabOp
    generic :: operator(*) => ScaleThreeBodyLabOp
  end type MACRO_ADD_SUFX(ThreeBodyLabOpIso)
contains

  subroutine CopyThreeBodyLabOp(a, b)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: a
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: b
    integer :: ichbra, ichket, n
    type(ThreeBodyLabIsoChan), pointer :: bra, ket
    if(allocated(a%MatCh)) call a%fin()
    n = size(b%MatCh, 1)
    call a%CopyOperatorDef(b%OperatorDef)
    a%ms => b%ms
    allocate(a%MatCh(n,n))
    do ichbra = 1, n
      do ichket = 1, ichbra
        bra => a%ms%GetChannel(ichbra)
        ket => a%ms%GetChannel(ichket)
        if(b%MatCh(ichbra,ichket)%is_Zero) then
          a%MatCh(ichbra,ichket)%is_zero = .true.
          cycle
        end if
        a%MatCh(ichbra,ichket)%is_Zero = .false.
        allocate(a%MatCh(ichbra,ichket)%m(bra%GetNumberStates(),ket%GetNumberStates()))
        a%MatCh(ichbra,ichket)%m(:,:) = b%MatCh(ichbra, ichket)%m(:,:)
      end do
    end do
    a%is_init = b%is_init
  end subroutine CopyThreeBodyLabOp

  function SumThreeBodyLabOp(a, b) result(c)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: a, b
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)) :: c
    integer :: ichbra, ichket, n
    c = a
    n = size(a%MatCh, 1)
    do ichbra = 1, n
      do ichket = 1, ichbra
        if(b%MatCh(ichbra,ichket)%is_Zero) then
          c%MatCh(ichbra,ichket)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichbra,ichket)%is_Zero = .false.
        c%MatCh(ichbra, ichket)%m = a%MatCh(ichbra, ichket)%m + b%MatCh(ichbra, ichket)%m
      end do
    end do
  end function SumThreeBodyLabOp

  function SubtractThreeBodyLabOp(a, b) result(c)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: a, b
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)) :: c
    integer :: ichbra, ichket, n
    c = a
    n = size(a%MatCh, 1)
    do ichbra = 1, n
      do ichket = 1, ichbra
        if(b%MatCh(ichbra,ichket)%is_Zero) then
          c%MatCh(ichbra,ichket)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichbra,ichket)%is_Zero = .false.
        c%MatCh(ichbra, ichket)%m = a%MatCh(ichbra, ichket)%m - b%MatCh(ichbra, ichket)%m
      end do
    end do
  end function SubtractThreeBodyLabOp

  function ScaleThreeBodyLabOp(a, b) result(c)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: a
    real(8), intent(in) :: b
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)) :: c
    integer :: ichbra, ichket, n
    c = a
    n = size(a%MatCh, 1)
    do ichbra = 1, n
      do ichket = 1, ichbra
        if(a%MatCh(ichbra,ichket)%is_Zero) then
          c%MatCh(ichbra,ichket)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichbra,ichket)%is_Zero = .false.
#if defined(half_precision) || defined(single_precision)
        c%MatCh(ichbra, ichket)%m = a%MatCh(ichbra, ichket)%m * real(b)
#else
        c%MatCh(ichbra, ichket)%m = a%MatCh(ichbra, ichket)%m * b
#endif
      end do
    end do
  end function ScaleThreeBodyLabOp

  subroutine FinThreeBodyLabOp(this)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    integer :: chbra, chket, n
    n = size(this%MatCh, 1)
    do chbra = 1, n
      do chket = 1, n
        if(.not. allocated(this%MatCh(chbra,chket)%m)) cycle
        deallocate(this%MatCh(chbra,chket)%m)
      end do
    end do
    call this%FinOperatorDef()
    deallocate(this%MatCh)
    this%is_init = .false.
  end subroutine FinThreeBodyLabOp

  subroutine InitThreeBodyLabOpFromString(this, ms, oprtr)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), intent(in) :: ms
    type(str), intent(in) :: oprtr

    call this%InitOpDef(oprtr, .false.)
    call this%init(ms,this%GetOpJ(),this%GetOpP(),this%GetOpT())
  end subroutine InitThreeBodyLabOpFromString

  subroutine InitThreeBodyLabOp(this, ms, jr, pr, tr)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), intent(in), target :: ms
    type(ThreeBodyLabIsoChan), pointer :: bra, ket
    integer, intent(in) :: jr, pr, tr
    integer :: chbra, chket
    integer :: nb, jb, pb, tb
    integer :: nk, jk, pk, tk

    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%SetOpJ(jr)
    call this%SetOpP(pr)
    call this%SetOpT(tr)

    allocate(this%MatCh(ms%GetNumberChannels() , ms%GetNumberChannels() ))
    do chbra = 1, ms%GetNumberChannels()
      bra => ms%GetChannel(chbra)
      nb = bra%GetNumberStates()
      jb = bra%GetJ()
      pb = bra%GetParity()
      tb = bra%GetT()
      !do chket = 1, ms%GetNumberChannels()
      do chket = 1, chbra
        ket => ms%GetChannel(chket)
        nk = ket%GetNumberStates()
        jk = ket%GetJ()
        pk = ket%GetParity()
        tk = ket%GetT()
        if(triag(jb, jk, 2*jr)) cycle
        if(pb * pr /= pk) cycle
        if(triag(tb, tk, 2*tr)) cycle
        this%MatCh(chbra, chket)%is_Zero = .false.
        allocate(this%MatCh(chbra,chket)%m(nb,nk))
#if defined(half_precision) || defined(single_precision)
        this%MatCh(chbra,chket)%m(:,:) = 0.0
#elif defined(double_precision)
        this%MatCh(chbra,chket)%m(:,:) = 0.0d0
#endif
      end do
    end do
    this%is_init = .true.
  end subroutine InitThreeBodyLabOp

  function GetThBME_scalar(this, a, b, c, Jab, Tab, d, e, f, Jde, Tde, J, T) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    integer, intent(in) :: a, b, c, Jab, Tab
    integer, intent(in) :: d, e, f, Jde, Tde
    integer, intent(in) :: J, T
#if defined(half_precision) || defined(single_precision)
    real(4) :: r
#elif defined(double_precision)
    real(8) :: r
#endif
    integer :: pari, ich, bra, ket
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(ThreeBodyLabIsoChan), pointer :: ChLab
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc

    r = real(0.d0, kind=kind(r))
    if( a<b .or. b<c ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( d<e .or. e<f ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    lab => this%ms
    sps => this%ms%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    pari = (-1)**( oa%l + ob%l + oc%l )
    ich = lab%GetIndex(J, pari, T)
    if( ich==0 ) return
    ChLab => lab%GetChannel(ich)
    bra = ChLab%GetIndex(a,b,c,Jab,Tab)
    ket = ChLab%GetIndex(d,e,f,Jde,Tde)
    if( bra*ket == 0 ) return
    r = this%MatCh(ich,ich)%m(bra,ket)
  end function GetThBME_scalar

  function GetThBME_tensor(this, a, b, c, Jab, Tab, Jbra, Tbra, d, e, f, Jde, Tde, Jket, Tket) result(r)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    integer, intent(in) :: a, b, c, Jab, Tab, Jbra, Tbra
    integer, intent(in) :: d, e, f, Jde, Tde, Jket, Tket
#if defined(half_precision) || defined(single_precision)
    real(4) :: r
#elif defined(double_precision)
    real(8) :: r
#endif
    integer :: ichbra, ichket, Pbra, Pket, bra, ket
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(ThreeBodyLabIsoChan), pointer :: LabBra, LabKet
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od, oe, of

    r = real(0.d0, kind=kind(r))
    if( a<b .or. b<c ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( d<e .or. e<f ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( triag(Jbra, Jket, 2*this%GetOpJ()) ) return
    if( triag(Tbra, Tket, 2*this%GetOpT()) ) return
    lab => this%ms
    sps => this%ms%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    oe => sps%GetOrbit(e)
    of => sps%GetOrbit(f)
    pbra = (-1)**( oa%l + ob%l + oc%l )
    pket = (-1)**( od%l + oe%l + of%l )
    if( Pbra * Pket * this%GetOpP() /= 1 ) return
    ichbra = lab%GetIndex(Jbra, Pbra, Tbra)
    ichket = lab%GetIndex(Jket, Pket, Tket)
    if( ichbra*ichket==0 ) return
    LabBra => lab%GetChannel(ichbra)
    LabKet => lab%GetChannel(ichket)
    bra = LabBra%GetIndex(a,b,c,Jab,Tab)
    ket = LabKet%GetIndex(d,e,f,Jde,Tde)
    if( bra*ket == 0 ) return
    if( ichbra < ichket ) then
      r = this%MatCh(ichket,ichbra)%m(ket,bra) * real((-1.d0)**( (Jket+Tket-Jbra-Tbra)/2 ), kind=kind(r))
      return
    end if
    r = this%MatCh(ichbra,ichket)%m(bra,ket)
  end function GetThBME_tensor

  subroutine SetThBME_scalar(this, a, b, c, Jab, Tab, d, e, f, Jde, Tde, J, T, ME)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    integer, intent(in) :: a, b, c, Jab, Tab
    integer, intent(in) :: d, e, f, Jde, Tde
    integer, intent(in) :: J, T
#if defined(half_precision) || defined(single_precision)
    real(4), intent(in) :: ME
#elif defined(double_precision)
    real(8), intent(in) :: ME
#endif
    integer :: pari, ich, bra, ket
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(ThreeBodyLabIsoChan), pointer :: ChLab
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc

    if( a<b .or. b<c ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( d<e .or. e<f ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    lab => this%ms
    sps => this%ms%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    pari = (-1)**( oa%l + ob%l + oc%l )
    ich = lab%GetIndex(J, pari, T)
    if( ich==0 ) return
    ChLab => lab%GetChannel(ich)
    bra = ChLab%GetIndex(a,b,c,Jab,Tab)
    ket = ChLab%GetIndex(d,e,f,Jde,Tde)
    if( bra*ket == 0 ) return
    this%MatCh(ich,ich)%m(bra,ket) = ME
    this%MatCh(ich,ich)%m(ket,bra) = ME
  end subroutine SetThBME_scalar

  subroutine SetThBME_tensor(this, a, b, c, Jab, Tab, Jbra, Tbra, d, e, f, Jde, Tde, Jket, Tket, ME)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    integer, intent(in) :: a, b, c, Jab, Tab, Jbra, Tbra
    integer, intent(in) :: d, e, f, Jde, Tde, Jket, Tket
#if defined(half_precision) || defined(single_precision)
    real(4), intent(in) :: ME
#elif defined(double_precision)
    real(8), intent(in) :: ME
#endif
    integer :: ichbra, ichket, Pbra, Pket, bra, ket
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(ThreeBodyLabIsoChan), pointer :: LabBra, LabKet
    type(OrbitsIsospin), pointer :: sps
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od, oe, of

    if( a<b .or. b<c ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( d<e .or. e<f ) then
      write(*,*) "Error, ", __LINE__, " in ", __FILE__
      return
    end if
    if( triag(Jbra, Jket, 2*this%GetOpJ()) ) return
    if( triag(Tbra, Tket, 2*this%GetOpT()) ) return
    lab => this%ms
    sps => this%ms%sps
    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    oe => sps%GetOrbit(e)
    of => sps%GetOrbit(f)
    pbra = (-1)**( oa%l + ob%l + oc%l )
    pket = (-1)**( od%l + oe%l + of%l )
    if( Pbra * Pket * this%GetOpP() /= 1 ) return
    ichbra = lab%GetIndex(Jbra, Pbra, Tbra)
    ichket = lab%GetIndex(Jket, Pket, Tket)
    if( ichbra*ichket==0 ) return
    LabBra => lab%GetChannel(ichbra)
    LabKet => lab%GetChannel(ichket)
    bra = LabBra%GetIndex(a,b,c,Jab,Tab)
    ket = LabKet%GetIndex(d,e,f,Jde,Tde)
    if( bra*ket == 0 ) return
    if( ichbra < ichket ) then
      this%MatCh(ichket,ichbra)%m(ket,bra) = ME * real((-1.d0)**( (Jket+Tket-Jbra-Tbra)/2 ),kind=kind(ME))
      return
    end if
    this%MatCh(ichbra,ichket)%m(bra,ket) = ME
    if( ichbra == ichket ) then
      this%MatCh(ichbra,ichket)%m(ket,bra) = ME
    end if
  end subroutine SetThBME_tensor

  function Truncate( this, ms_new ) result(op)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(ThreeBodyLabIsoSpace), intent(in), target :: ms_new
    type(MACRO_ADD_SUFX(ThreeBodyLabOpIso)) :: op
    type(ThreeBodyLabIsoSpace), pointer :: ms_old
    type(ThreeBodyLabIsoChan), pointer :: chbra_new, chket_new, chbra_old, chket_old
    integer :: ichbra, jbra, pbra, tbra, bra, abcJT(5)
    integer :: ichket, jket, pket, tket, ket, defJT(5)
    integer :: ichbra_old, ichket_old, bra_old, ket_old
    real(8) :: time
    type(sys) :: s

    time = omp_get_wtime()
    call op%init( ms_new, this%GetOpName() )
    ms_old => this%ms
    !$omp parallel
    !$omp do private( ichbra, ichket, chbra_new, chket_new, &
    !$omp &  jbra, pbra, tbra, jket, pket, tket, ichbra_old, ichket_old, &
    !$omp &  chbra_old, chket_old, bra, abcJT, bra_old, ket, defJT, ket_old ) schedule(dynamic)
    do ichbra = 1, ms_new%GetNumberChannels()
      do ichket = 1, ms_new%GetNumberChannels()
        if( op%MatCh(ichbra,ichket)%is_Zero ) cycle
        chbra_new => ms_new%GetChannel( ichbra )
        chket_new => ms_new%GetChannel( ichket )
        jbra = chbra_new%GetJ()
        pbra = chbra_new%GetParity()
        tbra = chbra_new%GetT()
        jket = chket_new%GetJ()
        pket = chket_new%GetParity()
        tket = chket_new%GetT()
        ichbra_old = ms_old%GetIndex( jbra, pbra, tbra )
        ichket_old = ms_old%GetIndex( jket, pket, tket )
        if( ichbra_old * ichket_old == 0) stop "Error in ThreeBodyLabOpsIso: Truncate, chan"
        chbra_old => ms_old%GetChannel( ichbra_old )
        chket_old => ms_old%GetChannel( ichket_old )

        do bra = 1, chbra_new%GetNumberStates()
          abcJT = chbra_new%GetABCJT( bra )
          bra_old = chbra_old%GetIndex(abcJT)
          do ket = 1, chket_new%GetNumberStates()
            defJT = chket_new%GetABCJT( ket )
            ket_old = chket_old%GetIndex(defJT)
            if(bra_old * ket_old == 0) stop "Error in ThreeBodyLabOpsIso: Truncate, sub"
            op%MatCh(ichbra,ichket)%m(bra,ket) = &
                & this%MatCh(ichbra_old,ichket_old)%m(bra_old,ket_old)
          end do
        end do

      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add(s%str("ThreeBodyLabOp: Truncate"), omp_get_wtime()-time)
  end function Truncate

  subroutine ReducedToNonReduced( this, convJ, convT )
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    logical, intent(in), optional :: convJ, convT
    logical :: cJ=.true.
    logical :: cT=.true.
    integer :: ichbra, ichket
    type(ThreeBodyLabIsoChan), pointer :: chbra, chket
#if defined(half_precision) || defined(single_precision)
    real(4) :: fact
#elif defined(double_precision)
    real(8) :: fact
#endif

    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. this%GetOpJ()/=0 ) then
      write(*,*) "Error, at ", __LINE__, " in ", __FILE__
    end if
    if( cT .and. this%GetOpT()/=0 ) then
      write(*,*) "Error, at ", __LINE__, " in ", __FILE__
    end if
    do ichbra = 1, this%ms%GetNumberChannels()
      chbra => this%ms%GetChannel( ichbra )
      do ichket = 1, this%ms%GetNumberChannels()
        chket => this%ms%GetChannel( ichket )
        if( this%MatCh(ichbra,ichket)%is_Zero ) cycle
        fact = 1.d0
        if( cJ ) fact = fact / real(sqrt( dble( chket%GetJ()+1 )), kind=kind(fact))
        if( cT ) fact = fact / real(sqrt( dble( chket%GetT()+1 )), kind=kind(fact))
        this%MatCh(ichbra,ichket)%m = this%MatCh(ichbra,ichket)%m * real(fact,kind=kind( fact ))
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReduced( this, convJ, convT )
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    logical, intent(in), optional :: convJ, convT
    logical :: cJ=.true.
    logical :: cT=.true.
    integer :: ichbra, ichket
    type(ThreeBodyLabIsoChan), pointer :: chbra, chket
#if defined(half_precision) || defined(single_precision)
    real(4) :: fact
#elif defined(double_precision)
    real(8) :: fact
#endif

    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. this%GetOpJ()/=0 ) then
      write(*,*) "Error, at ", __LINE__, " in ", __FILE__
    end if
    if( cT .and. this%GetOpT()/=0 ) then
      write(*,*) "Error, at ", __LINE__, " in ", __FILE__
    end if

    do ichbra = 1, this%ms%GetNumberChannels()
      chbra => this%ms%GetChannel( ichbra )
      do ichket = 1, this%ms%GetNumberChannels()
        chket => this%ms%GetChannel( ichket )
        if( this%MatCh(ichbra,ichket)%is_Zero ) cycle
        fact = 1.d0
        if( cJ ) fact = fact * real(sqrt( dble( chket%GetJ()+1 )), kind=kind(fact))
        if( cT ) fact = fact * real(sqrt( dble( chket%GetT()+1 )), kind=kind(fact))
        this%MatCh(ichbra,ichket)%m = this%MatCh(ichbra,ichket)%m * real(fact,kind=kind( fact ))
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NonReducedToReduced

  subroutine TalmiMoshinskyTransformation(this,jac2lab,opjac)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: jac2lab
    type(ThreeBodyJacOpIso), intent(in) :: opjac

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling TalmiMoshinskyTransformation"
      return
    end if

    if(.not. jac2lab%is_constructed) then
      write(*,*) "Initialize 'jac2lab' before calling TalmiMoshinskyTransformation"
      return
    end if

    if(.not. opjac%is_init) then
      write(*,*) "Initialize 'opjac' before calling TalmiMoshinskyTransformation"
      return
    end if


    if( .not. this%reduced_me() .and. .not. opjac%reduced_me() ) call TMTransScalar(this,jac2lab,opjac)
    if(       this%reduced_me() .and.       opjac%reduced_me() ) call TMTransTensor(this,jac2lab,opjac)
    if(       this%reduced_me() .and. .not. opjac%reduced_me() ) then
      write(*,*) "Error at", __LINE__, " in ", __FILE__
      stop
    end if
    if( .not. this%reduced_me() .and.       opjac%reduced_me() ) then
      write(*,*) "Error at", __LINE__, " in ", __FILE__
      stop
    end if
  end subroutine TalmiMoshinskyTransformation

  subroutine WriteThreeBodyLabOp(this, f, fn_average)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(str), intent(in), optional :: fn_average
    type(sys) :: s
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling WriteThreeBodyLabOp"
      return
    end if

    ti = omp_get_wtime()
    if(.not. this%reduced_me() ) call write_scalar_operator_isospin(this, f)
    if(      this%reduced_me() ) call write_tensor_operator_isospin(this, f)
    if( present(fn_average ))  call write_average_norm(this, fn_average )

    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyLabOp

  subroutine write_average_norm(this, fn)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: fn
    type(ThreeBodyLabIsoSpace), pointer :: ms
    type(ThreeBodyLabIsoChan), pointer :: chbra, chket
    type(OrbitsIsospin), pointer :: sps
    integer :: abcjt_bra(5), abcjt_ket(5)
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od, oe, of
    integer :: ichbra, ichket, bra, ket, ebra, eket
    integer :: wunit=22
    real(8), allocatable :: v_array(:,:), n_array(:,:), s_array(:,:)
    real(8), allocatable, save :: v_array_threads(:,:), n_array_threads(:,:)
    !$omp threadprivate( v_array_threads, n_array_threads )
    type(sys) :: s
    logical :: is_norm, is_averaged_norm, is_mean
#if defined(half_precision)
    real(4) :: vtmp
#endif

    ms => this%ms
    sps => this%ms%sps
    is_norm=.false.
    is_averaged_norm=.false.
    is_mean=.false.
    if( s%find( fn, s%str("AveragedNorm") ) ) is_averaged_norm = .true.
    if( s%find( fn, s%str("Norm") ) .and. (.not. is_averaged_norm) ) is_norm = .true.
    if( s%find( fn, s%str("Mean") ) ) is_mean = .true.
    if( (.not. is_norm) .and. (.not. is_averaged_norm) .and. (.not. is_mean) ) is_averaged_norm = .true. ! default

    if( is_norm ) write(*,"(2a)") "# Writing the norm file: ", trim(fn%val)
    if( is_averaged_norm ) write(*,"(2a)") "# Writing the averaged norm file: ", trim(fn%val)
    if( is_mean ) write(*,"(2a)") "# Writing the mean-value file: ", trim(fn%val)
    open(wunit, file=fn%val, action="write", status="replace")
    write(wunit,'(a)') "Jbra Pbra Tbra Ebra Jket Pket Tket Eket Mean Std dev."
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)


        allocate(v_array(0:ms%GetE3max(), 0:ms%GetE3max()) )
        allocate(n_array(0:ms%GetE3max(), 0:ms%GetE3max()) )
        allocate(s_array(0:ms%GetE3max(), 0:ms%GetE3max()) )
        v_array(:,:) = 0.d0
        n_array(:,:) = 0.d0
        s_array(:,:) = 0.d0
        !$omp parallel
        allocate(v_array_threads(0:ms%GetE3max(), 0:ms%GetE3max()) )
        allocate(n_array_threads(0:ms%GetE3max(), 0:ms%GetE3max()) )
        v_array_threads(:,:) = 0.d0
        n_array_threads(:,:) = 0.d0
        !$omp do private(bra, abcjt_bra, oa, ob, oc, ebra, ket, abcjt_ket, od, oe, of, eket)
        do bra = 1, chbra%GetNumberStates()
          abcjt_bra = chbra%GetABCJT( bra )
          oa => sps%GetOrbit( abcjt_bra(1) )
          ob => sps%GetOrbit( abcjt_bra(2) )
          oc => sps%GetOrbit( abcjt_bra(3) )
          ebra =  oa%e + ob%e + oc%e
          do ket = 1, chbra%GetNumberStates()
            abcjt_ket = chket%GetABCJT( ket )
            od => sps%GetOrbit( abcjt_ket(1) )
            oe => sps%GetOrbit( abcjt_ket(2) )
            of => sps%GetOrbit( abcjt_ket(3) )
            eket = od%e + oe%e + of%e
            n_array_threads(ebra,eket) = n_array_threads(ebra,eket) + 1.d0
#if defined(half_precision)
            vtmp = this%MatCh(ichbra,ichket)%m(bra,ket)
            if(is_norm .or. is_averaged_norm) then
              v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + vtmp**2
            end if
            if(is_mean) then
              v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + vtmp
            end if
#elif defined(single_precision) || defined(double_precision)
            if(is_norm .or. is_averaged_norm) then
              v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + this%MatCh(ichbra,ichket)%m(bra,ket)**2
            end if
            if(is_mean) then
              v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + this%MatCh(ichbra,ichket)%m(bra,ket)
            end if
#endif
          end do
        end do
        !$omp end do
        !$omp barrier

        !$omp critical
        n_array(:,:) = n_array(:,:) + n_array_threads(:,:)
        v_array(:,:) = v_array(:,:) + v_array_threads(:,:)
        !$omp end critical
        deallocate(n_array_threads)
        deallocate(v_array_threads)
        !$omp end parallel

        if( is_mean .or. is_averaged_norm ) then
          do ebra = 0, ms%GetE3max()
            do eket = 0, ms%GetE3max()
              v_array(ebra,eket) = v_array(ebra,eket) / n_array(ebra,eket)
            end do
          end do
          !$omp parallel
          allocate(v_array_threads(0:ms%GetE3max(), 0:ms%GetE3max()) )
          v_array_threads(:,:) = 0.d0
          !$omp do private(bra, abcjt_bra, oa, ob, oc, ebra, ket, abcjt_ket, od, oe, of, eket)
          do bra = 1, chbra%GetNumberStates()
            abcjt_bra = chbra%GetABCJT( bra )
            oa => sps%GetOrbit( abcjt_bra(1) )
            ob => sps%GetOrbit( abcjt_bra(2) )
            oc => sps%GetOrbit( abcjt_bra(3) )
            ebra =  oa%e + ob%e + oc%e
            do ket = 1, chbra%GetNumberStates()
              abcjt_ket = chket%GetABCJT( ket )
              od => sps%GetOrbit( abcjt_ket(1) )
              oe => sps%GetOrbit( abcjt_ket(2) )
              of => sps%GetOrbit( abcjt_ket(3) )
              eket = od%e + oe%e + of%e
#if defined(half_precision)
              vtmp = this%MatCh(ichbra,ichket)%m(bra,ket)
              if(is_averaged_norm) then
                v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + &
                    &( v_array(ebra,eket) - vtmp**2 )**2
              end if
              if(is_mean) then
                v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + &
                    &( v_array(ebra,eket) - vtmp )**2
              end if
#elif defined(single_precision) || defined(double_precision)
              if(is_averaged_norm) then
                v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + &
                    &( v_array(ebra,eket) - this%MatCh(ichbra,ichket)%m(bra,ket)**2 )**2
              end if
              if(is_mean) then
                v_array_threads(ebra,eket) = v_array_threads(ebra,eket) + &
                    &( v_array(ebra,eket) - this%MatCh(ichbra,ichket)%m(bra,ket) )**2
              end if
#endif
            end do
          end do
          !$omp end do
          !$omp barrier

          !$omp critical
          s_array(:,:) = s_array(:,:) + v_array_threads(:,:)
          !$omp end critical
          deallocate(v_array_threads)
          !$omp end parallel

          do ebra = 0, ms%GetE3max()
            if( chbra%GetParity() /= (-1)**ebra ) cycle
            do eket = 0, ms%GetE3max()
              if( chket%GetParity() /= (-1)**eket ) cycle
              if(n_array(ebra,eket) < 1.d-1) cycle

              write(wunit,'(8i3,2es18.8)') chbra%GetJ(), chbra%GetParity(), chbra%GetT(), ebra, &
                  & chket%GetJ(), chket%GetParity(), chket%GetT(), eket, v_array(ebra,eket), &
                  & sqrt(s_array(ebra,eket)/n_array(ebra,eket))

            end do
          end do
        else
          do ebra = 0, ms%GetE3max()
            if( chbra%GetParity() /= (-1)**ebra ) cycle
            do eket = 0, ms%GetE3max()
              if( chket%GetParity() /= (-1)**eket ) cycle
              if(n_array(ebra,eket) < 1.d-1) cycle
              write(wunit,'(8i3,es18.8)') chbra%GetJ(), chbra%GetParity(), chbra%GetT(), ebra, &
                  & chket%GetJ(), chket%GetParity(), chket%GetT(), eket, v_array(ebra,eket)

            end do
          end do
        end if

        deallocate(v_array, n_array, s_array)

      end do
    end do
    close(wunit)
  end subroutine write_average_norm

  subroutine write_scalar_operator_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s


    if(s%find(f,s%str('.me3j.gz'))) then
      call write_scalar_operator_gzip_me3j_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j'))) then
      call write_scalar_operator_ascii_me3j_isospin(this,f)
      return
    end if

    if(s%find(f,s%str('.stream.bin'))) then
      call write_scalar_operator_binary_stream_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.bin'))) then
      call write_scalar_operator_binary_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.readable.txt'))) then
      call write_operator_humanreadable_ascii_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.txt'))) then
      call write_scalar_operator_ascii_isospin(this, f)
      return
    end if

    write(*,'(a)') ' Unknown format, writting file assuming readable me3j'
    call write_scalar_operator_ascii_me3j_isospin(this,f)
    return
  end subroutine write_scalar_operator_isospin

  subroutine write_tensor_operator_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str('.bin.stream'))) then
      call write_tensor_operator_binary_stream_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.me3j.gz'))) then
      call write_tensor_operator_gzip_isospin(this, f)
      return
    end if

    if(s%find(f,s%str('.readable.txt'))) then
      call write_operator_humanreadable_ascii_isospin(this, f)
      return
    end if

    write(*,'(a)') ' Unknown format, writting file assuming binary stream i/o format'
    call write_tensor_operator_binary_stream_isospin(this, f)
    return

  end subroutine write_tensor_operator_isospin

  subroutine TMTransScalar(this,jac2lab,opjac)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: jac2lab
    type(ThreeBodyJacOpIso), intent(in) :: opjac
    type(ThreeBodyJacIsoSpace), pointer :: jac
    type(OrbitsIsospin), pointer :: sps
    type(ThreeBodyLabIsoSpace), pointer :: lab
    integer :: e3max
    integer :: jtot, ptot, ttot, nch, ich
    integer :: jrel, prel, ncm, lcm
    integer, allocatable :: jpnl(:,:)
    integer, allocatable :: jpttemp(:,:)
    integer, allocatable :: slranks(:)
    type(sys) :: s
#ifdef MPI
    type(ThreeBodyLabIsoChan), pointer :: ChLab
    integer :: ch
#ifdef SPARC
      integer :: n1d
#else
      integer(8) :: n1d
#endif
#endif
    real(8) :: time_ms

    jac => opjac%ms
    lab => this%ms
    sps => this%ms%sps
    e3max = lab%GetE3max()
#ifndef MPI
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis tuned using DGEMM ...")')
    write(*,'(4x, a, i4)') "Block number is", jac2lab%GetNumberBlock()
    write(*,'("############################################################")')
    write(*,*)
#endif

    nch = min(e3max + 2, (jac%GetJmax()+1)/2) * 2 * (e3max/2 + 1) * (e3max + 1)
    allocate( jpnl(4, nch) )
    ich = 0
    do lcm = 0, e3max
      do ncm = 0, e3max / 2
        do prel = -1, 1, 2
          do jrel = 1, min(2*e3max + 3, jac%GetJmax()), 2
            ich = ich + 1
            jpnl(:, ich) = [jrel, prel, ncm, lcm]
          end do
        end do
      end do
    end do

    allocate(jpttemp(3, (e3max + 2) * 2 * 2))
    ich = 0
    do ttot = 1, 3, 2
      do jtot = 1,2 * e3max + 3, 2
        do ptot = 1, -1, -2
          ich = ich + 1
          jpttemp(:,ich) = [jtot, ptot, ttot]
        end do
      end do
    end do

    call parent_child_procedure(TMTransChannel, (e3max+2)*4, slranks, time_ms)
#ifdef MPI
    call timer%Add(s%str('TMTransScalarIsospin'),time_ms)
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(slranks(1), (e3max+2)*4, mpi_integer, 0, mpi_comm_world, ierr)
    do ich = 1, (e3max+2)*4
      jtot = jpttemp(1, ich)
      ptot = jpttemp(2, ich)
      ttot = jpttemp(3, ich)
      ch = lab%GetIndex(Jtot,Ptot,Ttot)
      if(ch == 0) cycle
      ChLab => lab%GetChannel(ch)
      n1d = int(ChLab%GetNumberStates(),kind(n1d)) * int(ChLab%GetNumberStates(),kind(n1d))
      if(n1d < 1) cycle
#if defined(half_precision)
      call mpi_bcast(this%MatCh(ch,ch)%m(1,1)%i16,n1d,mpi_integer2,slranks(ich),mpi_comm_world,ierr)
#elif defined(single_precision)
      call mpi_bcast(this%MatCh(ch,ch)%m(1,1),n1d,mpi_real4,slranks(ich),mpi_comm_world,ierr)
#elif defined(double_precision)
      call mpi_bcast(this%MatCh(ch,ch)%m(1,1),n1d,mpi_real8,slranks(ich),mpi_comm_world,ierr)
#endif
    end do
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(slranks)
  contains
    subroutine TMTransChannel(loop)
      integer, intent(in) :: loop
      type r8
#if defined(half_precision) || defined(single_precision)
        real(4), allocatable :: v(:,:)
#elif defined(double_precision)
        real(8), allocatable :: v(:,:)
#endif
        logical :: is
        integer :: jpnl(4)
      end type r8
      type(r8), allocatable :: tc(:,:,:,:)
      type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), pointer :: ChJac2Lab
      type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: asqns
      type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
      integer :: nblock
      integer :: nlab, njac, iblock, nbc, ichrel, ichlab
      integer :: ndim, i, m, k
      integer :: vbra, vket, tcbra, tcket, nabc
#if defined(half_precision) || defined(single_precision)
      real(4) :: fac
      real(4), allocatable :: work(:,:), v3(:,:), opjacm(:,:)
#elif defined(double_precision)
      real(8) :: fac
      real(8), allocatable :: work(:,:), v3(:,:)
#endif
#ifndef MPI
      real(8) :: ti
      ti = omp_get_wtime()
#endif
      nblock = jac2lab%GetNumberBlock()
      jtot = jpttemp(1, loop)
      ptot = jpttemp(2, loop)
      ttot = jpttemp(3, loop)
      ichlab = lab%GetIndex(jtot,ptot,ttot)
      if(ichlab == 0) return
      write(*,'(a, i4, a, i2, a, i2, a, i2, a)') 'myrank=', myrank, ',  :Calculating for J = ', &
          & jtot, '/2,   P = ', ptot, ',   T = ', ttot, '/2'
      nch = min(e3max + 2, (jac%GetJmax()+1)/2) * 2 * (e3max/2 + 1) * (e3max + 1)
      ChJac2Lab => jac2lab%GetChannel(ichlab)
      nlab = ChJac2Lab%GetNumberLabStates()
      njac = ChJac2Lab%GetNumberASJacCMStates()
      if(nlab < 1) return
      if(njac < 1) return
      allocate(v3(nlab,nlab))
      allocate(tc(1:2*e3max+3,-1:1,0:e3max/2, 0:e3max))
      v3(:,:) = 0.d0
      do iblock = 1, nch, nblock
        do nbc = 1, min(nblock, nch - iblock + 1)
          jrel = jpnl(1, iblock+nbc-1)
          prel = jpnl(2, iblock+nbc-1)
          ncm  = jpnl(3, iblock+nbc-1)
          lcm  = jpnl(4, iblock+nbc-1)
          tc(jrel, prel, ncm, lcm)%is = .false.
          ichrel = jac%GetIndex(jrel, prel, ttot)
          if(ichrel == 0) cycle
          if(.not. allocated(opjac%MatCh(ichrel,ichrel)%m)) cycle
          ndim = jac%jpt(ichrel)%GetNumberAStates()
          allocate(tc(jrel,prel,ncm,lcm)%v(ndim,nlab))
          tc(jrel, prel, ncm, lcm)%v(:,:) = 0.d0
        end do

        !$omp parallel
        !$omp do private(i, asqns, nabc, ichrel, vbra, tcbra, vket, labqns, tcket)
        do i = 1, njac
          asqns => ChJac2Lab%GetASJacCM( i )
          nabc = asqns%Nmax
          if(nabc > e3max) cycle
          ichrel = jac%GetIndex(asqns%jrel, (-1)**(asqns%E), ttot)
          if(ichrel == 0) cycle
          if(jac%jpt(ichrel)%GetNumberAStates() < 1) cycle
          vbra = jac%jpt(ichrel)%GetASIndex(asqns%E,asqns%i)
          !tcbra = jac2lab%jpt(ichlab)%n2jac(i)
          tcbra = asqns%idx
          if(.not. allocated(ChJac2Lab%N3max(nabc)%mat)) cycle
          if(.not. allocated(tc(asqns%jrel, (-1)**asqns%E, asqns%ncm, asqns%lcm)%v) ) cycle

          do vket = 1, nlab
            labqns => ChJac2Lab%GetLab( vket )
            if(labqns%Nmax /= nabc) cycle
            !tcket = jac2lab%jpt(ichlab)%n2lab(vket)
            tcket = labqns%idx
            tc(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%is = .true.
            tc(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%v( vbra, vket ) = ChJac2Lab%N3max(nabc)%mat( tcbra, tcket )
          end do
        end do
        !$omp end do
        !$omp end parallel

        do nbc = 1, min(nblock, nch - iblock + 1)
          jrel = jpnl(1, iblock+nbc-1)
          prel = jpnl(2, iblock+nbc-1)
          ncm  = jpnl(3, iblock+nbc-1)
          lcm  = jpnl(4, iblock+nbc-1)
          ichrel = jac%GetIndex(jrel, prel, ttot)
          if(ichrel == 0) cycle
          if (.not. allocated(opjac%MatCh(ichrel,ichrel)%m)) cycle
          if ( .not. tc( jrel, prel, ncm, lcm )%is ) cycle
          m = size( opjac%MatCh(ichrel,ichrel)%m, 1)
          k = size( tc(jrel,prel,ncm,lcm)%v, 2)
          if (m < 1 .or. k < 1) cycle
#if defined(half_precision) || defined(single_precision)
          fac = 6.0
#elif defined(double_precision)
          fac = 6.d0
#endif
          allocate(work(m,k))
#if defined(half_precision) || defined(single_precision)
          allocate(opjacm(m,m))
          opjacm = real(opjac%MatCh(ichrel,ichrel)%m)
          call sgemm( 'n', 'n', m, k, m, 1.0, opjacm, m, tc(jrel,prel,ncm,lcm)%v, m, 0.0, work, m )
          call sgemm( 't', 'n', k, k, m, fac, tc(jrel,prel,ncm,lcm)%v, m, work, m, 1.0, v3, k )
          deallocate(opjacm)
#elif defined(double_precision)
          call dgemm( 'n', 'n', m, k, m, 1.d0, opjac%MatCh(ichrel,ichrel)%m, m, tc(jrel,prel,ncm,lcm)%v, m, 0.d0, work, m )
          call dgemm( 't', 'n', k, k, m, fac, tc(jrel,prel,ncm,lcm)%v, m, work, m, 1.d0, v3, k )
#endif
          deallocate(work)
        end do

        do nbc = 1, min(nblock, nch - iblock + 1)
          jrel = jpnl(1, iblock+nbc-1)
          prel = jpnl(2, iblock+nbc-1)
          ncm  = jpnl(3, iblock+nbc-1)
          lcm  = jpnl(4, iblock+nbc-1)
          if(allocated(tc(jrel,prel,ncm,lcm)%v)) then
            deallocate(tc(jrel,prel,ncm,lcm)%v)
          end if
        end do

      end do
      deallocate(tc)

      this%MatCh(ichlab,ichlab)%m(:,:) =  v3(:,:)
      deallocate(v3)
#ifndef MPI
      call timer%Add(s%str('TMTransChannel'),omp_get_wtime()-ti)
#endif
    end subroutine TMTransChannel
  end subroutine TMTransScalar

  subroutine TMTransTensor(this,jac2lab,opjac)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(MACRO_ADD_SUFX(TransJac2LabSpaceIsospin)), intent(in) :: jac2lab
    type(ThreeBodyJacOpIso), intent(in) :: opjac
    type(ThreeBodyJacIsoSpace), pointer :: jac
    type(OrbitsIsospin), pointer :: sps
    type(ThreeBodyLabIsoSpace), pointer :: lab
    integer :: e3max
    integer :: ncm, lcm, cnt
    integer :: nchloop
    integer :: jbrel, pbrel, jkrel, pkrel
    integer :: ichb, ichk, ch1dim
    integer :: ttotb, ttotk
    integer, allocatable :: slranks(:)
    integer, allocatable :: jpnl_loop(:,:), loop2chb(:), loop2chk(:)
    type(sys) :: s
#ifdef MPI
    type(ThreeBodyLabIsoChan), pointer :: Ch_bra, Ch_ket
    integer :: ich, chbra, chket
#ifdef SPARC
    integer :: n1d
#else
    integer(8) :: n1d
#endif
#endif
    real(8) :: ti, time_ms

    lab => this%ms
    jac => opjac%ms
    sps => this%ms%sps
    e3max = lab%GetE3max()
#ifndef MPI
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis tuned using DGEMM ...")')
    write(*,'(4x, a, i4)') "Block number is", jac2lab%GetNumberBlock()
    write(*,'("############################################################")')
    write(*,*)
#endif

    ti =  omp_get_wtime()

    nchloop = min(e3max+2, (jac%GetJmax()+1)/2) ** 2 * 4 * (e3max/2+1) * (e3max+1)
    allocate( jpnl_loop(6,nchloop) )
    cnt = 0
    do lcm = 0, e3max
      do ncm = 0, e3max/2
        do pbrel = 1, -1, -2
          do jbrel = 1, min(2*e3max + 3, jac%GetJmax()), 2
            do pkrel = 1, -1, -2
              do jkrel = 1, min(2*e3max + 3, jac%GetJmax()), 2
                cnt = cnt + 1
                jpnl_loop(:,cnt) = (/jbrel, pbrel, jkrel, pkrel, ncm, lcm/)
              end do
            end do
          end do
        end do
      end do
    end do

    ch1dim = 0
    do ichb = 1, lab%GetNumberChannels()
      do ichk = 1, lab%GetNumberChannels()
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        ch1dim = ch1dim + 1
      end do
    end do

    allocate(loop2chb(ch1dim), loop2chk(ch1dim))
    ch1dim = 0
    do ichb = 1, lab%GetNumberChannels()
      do ichk = 1, lab%GetNumberChannels()
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        ch1dim = ch1dim + 1
        loop2chb(ch1dim) = ichb
        loop2chk(ch1dim) = ichk
      end do
    end do

    call parent_child_procedure(TMTransChannel, ch1dim, slranks, time_ms)
#ifdef MPI
    call timer%Add(s%str('TMTransTensorIsospin'),time_ms)
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(slranks(1), ch1dim, mpi_integer, 0, mpi_comm_world, ierr)
    do ich = 1, ch1dim
      chbra = loop2chb(ich)
      chket = loop2chk(ich)
      Ch_bra => lab%GetChannel(chbra)
      Ch_ket => lab%GetChannel(chket)
      n1d = int(Ch_bra%GetNumberStates(),kind(n1d)) * int(Ch_ket%GetNumberStates(),kind(n1d))
      if(n1d < 1) cycle
#if defined(half_precision)
      call mpi_bcast(this%MatCh(chbra,chket)%m(1,1)%i16,n1d,mpi_integer2,slranks(ich),mpi_comm_world,ierr)
#elif defined(single_precision)
      call mpi_bcast(this%MatCh(chbra,chket)%m(1,1),n1d,mpi_real4,slranks(ich),mpi_comm_world,ierr)
#elif defined(double_precision)
      call mpi_bcast(this%MatCh(chbra,chket)%m(1,1),n1d,mpi_real8,slranks(ich),mpi_comm_world,ierr)
#endif
    end do
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(slranks)
    deallocate(jpnl_loop, loop2chb, loop2chk)
  contains
    subroutine TMTransChannel(loop)
      integer, intent(in) :: loop
      type r8
#if defined(half_precision) || defined(single_precision)
        real(4), allocatable :: v(:,:)
#elif defined(double_precision)
        real(8), allocatable :: v(:,:)
#endif
        logical :: is
        integer :: jpnl(4)
      end type r8
      type(r8), allocatable :: tcb(:,:,:,:), tck(:,:,:,:)
      type(ThreeBodyLabIsoChan), pointer :: ChLabBra, ChLabKet
      type(MACRO_ADD_SUFX(TransJac2LabChannelIsospin)), pointer :: Jac2LabBra, Jac2LabKet
      type(MACRO_ADD_SUFX(ASJacobiCMQNs)), pointer :: asqns
      type(MACRO_ADD_SUFX(ThreeBodyLabQNs)), pointer :: labqns
      integer :: nblock
      integer :: nblab, nbjac, nklab, nkjac
      integer :: iblock, nbc, ichbrel, ichkrel
      integer :: i, nbdim, nkdim
      integer :: vbra, vket, tcbra, tcket, nabc
      integer :: nloop, nb, nk, mb, mk, jk, jb
#if defined(half_precision) || defined(single_precision)
      real(4), allocatable :: work(:,:), v3(:,:), opjacm(:,:)
      real(4) :: fac
#elif defined(double_precision)
      real(8), allocatable :: work(:,:), v3(:,:)
      real(8) :: fac
#endif

#ifndef MPI
      ti = omp_get_wtime()
#endif
      nblock = jac2lab%GetNumberBlock()
      ichb = loop2chb(loop)
      ichk = loop2chk(loop)

      ChLabBra => lab%GetChannel(ichb)
      ChLabKet => lab%GetChannel(ichk)
      Jac2LabBra => jac2lab%GetChannel(ichb)
      Jac2LabKet => jac2lab%GetChannel(ichk)
      jb = ChLabBra%GetJ()
      jk = ChLabKet%GetJ()
      ttotb = ChLabBra%GetT()
      ttotk = ChLabKet%GetT()
      write(*,'(a, i4, a, i2, a, i2)') 'myrank=', myrank, ',  :Calculating for bra = ', &
          & loop2chb(loop), ',   ket = ', loop2chk(loop)
      nloop = min(e3max + 2, (jac%GetJmax()+1)/2)**2 * 4 * (e3max/2 + 1) * (e3max + 1)
      nblab = Jac2LabBra%GetNumberLabStates()
      nklab = Jac2Labket%GetNumberLabStates()
      nbjac = Jac2LabBra%GetNumberASJacCMStates()
      nkjac = Jac2Labket%GetNumberASJacCMStates()
      if(nblab < 1) return
      if(nbjac < 1) return
      if(nklab < 1) return
      if(nkjac < 1) return
      allocate(v3(nblab, nklab))
      allocate(tcb(1:2*e3max+3, -1:1, 0:e3max/2, 0:e3max))
      allocate(tck(1:2*e3max+3, -1:1, 0:e3max/2, 0:e3max))
      v3 = 0.d0

      do iblock = 1, nloop, nblock
        do nbc = 1, min(nblock, nloop - iblock + 1)
          jbrel = jpnl_loop(1, iblock+nbc-1)
          pbrel = jpnl_loop(2, iblock+nbc-1)
          jkrel = jpnl_loop(3, iblock+nbc-1)
          pkrel = jpnl_loop(4, iblock+nbc-1)
          ncm   = jpnl_loop(5, iblock+nbc-1)
          lcm   = jpnl_loop(6, iblock+nbc-1)
          tcb(jbrel,pbrel,ncm,lcm)%is = .false.
          tck(jkrel,pkrel,ncm,lcm)%is = .false.
          ichbrel = jac%GetIndex(jbrel,pbrel,ttotb)
          ichkrel = jac%GetIndex(jkrel,pkrel,ttotk)
          if(opjac%MatCh(ichbrel,ichkrel)%is_Zero) cycle
          nbdim = jac%jpt(ichbrel)%GetNumberAStates()
          nkdim = jac%jpt(ichkrel)%GetNumberAStates()

          if(.not. allocated(tcb(jbrel,pbrel,ncm,lcm)%v)) allocate(tcb(jbrel,pbrel,ncm,lcm)%v(nbdim,nblab))
          if(.not. allocated(tck(jkrel,pkrel,ncm,lcm)%v)) allocate(tck(jkrel,pkrel,ncm,lcm)%v(nkdim,nklab))

          tcb(jbrel,pbrel,ncm,lcm)%v(:,:) = 0.d0
          tck(jkrel,pkrel,ncm,lcm)%v(:,:) = 0.d0

        end do

        !$omp parallel
        !$omp do private(i, asqns, nabc, ichbrel, vbra, tcbra, vket, labqns, tcket)
        do i = 1, nbjac
          asqns => Jac2LabBra%GetASJacCM( i )
          nabc = asqns%Nmax
          if(nabc > e3max) cycle
          ichbrel = jac%GetIndex(asqns%jrel, (-1)**asqns%E, ttotb)
          if(ichbrel == 0) cycle
          if(jac%jpt(ichbrel)%GetNumberAStates() < 1) cycle
          vbra = jac%jpt(ichbrel)%GetASIndex(asqns%E,asqns%i)
          tcbra = asqns%idx
          if(.not. allocated(Jac2LabBra%N3max(nabc)%mat)) cycle
          if(.not. allocated(tcb(asqns%jrel, (-1)**asqns%E, asqns%ncm, asqns%lcm)%v) ) cycle

          do vket = 1, nblab
            labqns => Jac2LabBra%GetLab( vket )
            if(labqns%Nmax /= nabc) cycle
            tcket = labqns%idx
            tcb(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%is = .true.
            tcb(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%v( vbra, vket) = Jac2LabBra%N3max(nabc)%mat( tcbra, tcket )
          end do
        end do
        !$omp end do
        !$omp end parallel

        !$omp parallel
        !$omp do private(i, asqns, nabc, ichkrel, vbra, tcbra, vket, labqns, tcket)
        do i = 1, nkjac
          asqns => Jac2LabKet%GetASJacCM( i )
          nabc = asqns%Nmax
          if(nabc > e3max) cycle
          ichkrel = jac%GetIndex(asqns%jrel, (-1)**asqns%E, ttotk)
          if(ichkrel == 0) cycle
          if(jac%jpt(ichkrel)%GetNumberAStates() < 1) cycle
          vbra = jac%jpt(ichkrel)%GetASIndex(asqns%E, asqns%i)
          tcbra = asqns%idx
          if(.not. allocated(Jac2LabKet%N3max(nabc)%mat)) cycle
          if(.not. allocated(tck(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%v) ) cycle

          do vket = 1, nklab
            labqns => Jac2LabKet%GetLab( vket )
            if(labqns%Nmax /= nabc) cycle
            tcket = labqns%idx
            tck(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%is = .true.
            tck(asqns%jrel,(-1)**asqns%E,asqns%ncm,asqns%lcm)%v( vbra, vket ) = Jac2LabKet%N3max(nabc)%mat( tcbra, tcket )
          end do
        end do
        !$omp end do
        !$omp end parallel

        do nbc = 1, min(nblock, nloop - iblock + 1)
          jbrel = jpnl_loop(1, iblock+nbc-1)
          pbrel = jpnl_loop(2, iblock+nbc-1)
          jkrel = jpnl_loop(3, iblock+nbc-1)
          pkrel = jpnl_loop(4, iblock+nbc-1)
          ncm   = jpnl_loop(5, iblock+nbc-1)
          lcm   = jpnl_loop(6, iblock+nbc-1)
          ichbrel = jac%GetIndex(jbrel,pbrel,ttotb)
          ichkrel = jac%GetIndex(jkrel,pkrel,ttotk)
          if(opjac%MatCh(ichbrel,ichkrel)%is_Zero) cycle
          if(.not. tcb(jbrel,pbrel,ncm,lcm)%is) cycle
          if(.not. tck(jbrel,pbrel,ncm,lcm)%is) cycle

          nb = size(tcb(jbrel, pbrel, ncm, lcm)%v, 2)
          mb = size(tcb(jbrel, pbrel, ncm, lcm)%v, 1)

          nk = size(tck(jkrel, pkrel, ncm, lcm)%v, 2)
          mk = size(tck(jkrel, pkrel, ncm, lcm)%v, 1)
          if(mk < 1 .or. nk < 1) cycle
          allocate(work(mb, nk))
#if defined(half_precision) || defined(single_precision)
          fac = real(6.d0 * (-1.d0) ** (lcm + (jkrel + jb)/2 + opjac%GetOpJ()) * hat(jb) * hat(jk) * &
              & sjs(jbrel, jkrel, 2*opjac%GetOpJ(), jk, jb, 2*lcm))
#elif defined(double_precision)
          fac = 6.d0 * (-1.d0) ** (lcm + (jkrel + jb)/2 + opjac%GetOpJ()) * hat(jb) * hat(jk) * &
              & sjs(jbrel, jkrel, 2*opjac%GetOpJ(), jk, jb, 2*lcm)
#endif

#if defined(half_precision) || defined(single_precision)
          allocate(opjacm(mb,mk))
          opjacm = real(opjac%MatCh(ichbrel,ichkrel)%m)
          call sgemm('n','n',mb,nk,mk,1.0,opjacm,mb,tck(jkrel,pkrel,ncm,lcm)%v,mk,0.0,work,mb)
          call sgemm('t', 'n', nb, nk, mb, fac, tcb(jbrel, pbrel, ncm, lcm)%v, mb, work, mb, 1.0, v3, nb)
          deallocate(opjacm)
#elif defined(double_precision)
          call dgemm('n','n',mb,nk,mk,1.d0,opjac%MatCh(ichbrel,ichkrel)%m,mb,tck(jkrel,pkrel,ncm,lcm)%v,mk,0.d0,work,mb)
          call dgemm('t', 'n', nb, nk, mb, fac, tcb(jbrel, pbrel, ncm, lcm)%v, mb, work, mb, 1.d0, v3, nb)
#endif
          deallocate(work)
        end do

        do nbc = 1, min(nblock, nloop - iblock + 1)
          jbrel = jpnl_loop(1, iblock+nbc-1)
          pbrel = jpnl_loop(2, iblock+nbc-1)
          jkrel = jpnl_loop(3, iblock+nbc-1)
          pkrel = jpnl_loop(4, iblock+nbc-1)
          ncm   = jpnl_loop(5, iblock+nbc-1)
          lcm   = jpnl_loop(6, iblock+nbc-1)
          if(allocated(tcb(jbrel,pbrel,ncm,lcm)%v)) then
            deallocate(tcb(jbrel,pbrel,ncm,lcm)%v)
          end if
          if(allocated(tck(jkrel,pkrel,ncm,lcm)%v)) then
            deallocate(tck(jkrel,pkrel,ncm,lcm)%v)
          end if
        end do

      end do
      this%MatCh(ichb, ichk)%m(:,:) =  v3(:,:)
      deallocate(v3)
      deallocate(tcb)
      deallocate(tck)
#ifndef MPI
      call timer%Add(s%str('TMTransChannel'),omp_get_wtime()-ti)
#endif
    end subroutine TMTransChannel
  end subroutine TMTransTensor

  function GetFileNameThreeBodyLabOp(this, filename, &
        &  NNInt, NNNInt, renorm, lambda, hw, emax, e2max, e3max, &
        &  cd, ce, genuine3bf, Rtype, Rpower) result(f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: NNInt, NNNInt, renorm, Rtype, filename
    real(8), intent(in) :: lambda, hw, cd, ce
    integer, intent(in) :: emax, e2max, e3max, Rpower
    logical, intent(in) :: genuine3bf
    type(str) :: f, OpName
    type(sys) :: s

    OpName = this%GetOpName()
    if(filename%val /= "default") then
      select case(OpName%val)
      case("hamil","Hamil","NNNint")
        f = filename
      case default
        f = this%GetOpName() + s%str("_") + filename
      end select
      return
    end if

    if(NNNInt%val /= 'ChEFT_N2LO') then
      write(*,'(a)') 'NNNInt other than ChEFT_N2LO has not been implemented'
      stop
    end if
    f = this%GetOpName() + s%str('_ThBME')
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') f = 'ThBME'
    if(renorm%val == 'bare') then
      if(genuine3bf) then
        f = f + s%str('_')  + NNNInt + s%str('_cD') + s%str(cd)
        f = f + s%str('cE')  + s%str(ce) + s%str('_')
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
  end function GetFileNameThreeBodyLabOp

  subroutine write_operator_humanreadable_ascii_isospin( this, f )
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: ms
    type(ThreeBodyLabIsoChan), pointer :: chbra, chket
    integer :: ichbra, ichket, bra, ket, ketmax
    character(256) :: header
    integer :: abcJT(5), defJT(5)
    type(str) :: OpName
#if defined(half_precision)
    real(4) :: vtmp
#endif


    ms => this%ms
    OpName = this%GetOpName()
    header = trim(OpName%val)
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NNN int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    open(wunit, file=f%val, status='replace', action='write')
    write(wunit,'(a)') trim(header)
    write(wunit,'(a)') "   a,  b,  c,Jab,Tab,  d,  e,  f,Jde,Tde,Jabc,Tabc,Jdef,Tdef,               ME"
    do ichbra = 1, ms%GetNumberChannels()
      chbra => ms%GetChannel(ichbra)
      do ichket = 1, ichbra
        chket => ms%GetChannel(ichket)
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        do bra = 1, chbra%GetNumberStates()
          ketmax = chket%GetNumberStates()
          if( ichbra == ichket ) ketmax = bra
          do ket = 1, ketmax
            abcJT = chbra%GetABCJT(bra)
            defJT = chket%GetABCJT(ket)
#if defined(half_precision)
            vtmp = this%MatCh(ichbra,ichket)%m(bra,ket)
            if( abs(vtmp) < 1.d-8) cycle
            write(wunit,"(10i4,4i5,f18.10)")  abcJT, defJT, &
                &chbra%GetJ(), chbra%GetT(), chket%GetJ(), chket%GetT(),vtmp
#elif defined(single_precision) || defined(double_precision)
            if( abs(this%MatCh(ichbra,ichket)%m(bra,ket)) < 1.d-8) cycle
            write(wunit,"(10i4,4i5,f18.10)")  abcJT, defJT, &
                &chbra%GetJ(), chbra%GetT(), chket%GetJ(), chket%GetT(),this%MatCh(ichbra,ichket)%m(bra,ket)
#endif
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine write_operator_humanreadable_ascii_isospin

  subroutine write_scalar_operator_ascii_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(str) :: OpName
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, j123
    integer :: cnt
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf
    integer(8) :: total_num=0
#if defined(half_precision) || defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    OpName = this%GetOpName()
    header = trim(OpName%val)
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NNN int.'
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
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = max(abs(2*j12-j3),abs(2*j45-j6)), min(2*j12+j3,2*j45+j6), 2

                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t123 = 1, min(2*t12+1,2*t45+1), 2


                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J123,T123)
                            v_buffer(cntbf) = v123
                            if(cntbf == numbf) then
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
  end subroutine write_scalar_operator_ascii_isospin

  subroutine write_scalar_operator_ascii_me3j_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(str) :: OpName
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, j123
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf
    integer(8) :: total_num=0
#if defined(half_precision) || defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    character(256) :: header

    open(wunit, file=f%val, status='replace', action='write')
    OpName = this%GetOpName()
    header = trim(OpName%val)
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NNN int.'
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
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = max(abs(2*j12-j3),abs(2*j45-j6)), min(2*j12+j3,2*j45+j6), 2

                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t123 = 1, min(2*t12+1,2*t45+1), 2


                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J123,T123)
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
        end do
      end do
    end do

    if(cntbf > 0) then
      write(wunit,'(10f16.8)') v_buffer(:cntbf)
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_scalar_operator_ascii_me3j_isospin

  subroutine write_scalar_operator_gzip_me3j_isospin(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, j123
    integer :: cnt
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf, nrest
    integer(8) :: total_num=0
#if defined(half_precision) || defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    type(c_ptr) :: fp, err
    character(256) :: header, buffer
    character(12) :: cfmt
    type(str) :: OpName

    fp = gzip_open(f%val, "wt")
    OpName = this%GetOpName()
    header = trim(OpName%val)
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NNN int.'
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
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = max(abs(2*j12-j3),abs(2*j45-j6)), min(2*j12+j3,2*j45+j6), 2

                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t123 = 1, min(2*t12+1,2*t45+1), 2


                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J123,T123)
                            v_buffer(cntbf) = v123
                            if(cntbf == numbf) then
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

    if(cntbf > 0) then
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
  end subroutine write_scalar_operator_gzip_me3j_isospin

  subroutine write_scalar_operator_binary_stream_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, j123
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf
    integer(8) :: total_num=0
#if defined(half_precision)
    integer(2), allocatable :: v_buffer(:)
    type(my_real16) :: v123
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif

    open(wunit,file=f%val,form='unformatted',access='stream',status='replace')
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = max(abs(2*j12-j3),abs(2*j45-j6)), min(2*j12+j3,2*j45+j6), 2

                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t123 = 1, min(2*t12+1,2*t45+1), 2


                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J123,T123)
                            v_buffer(cntbf) = v123
                            if(cntbf == numbf) then
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
      write(wunit) v_buffer(:cntbf)
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_scalar_operator_binary_stream_isospin

  subroutine write_scalar_operator_binary_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, j123
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf
    integer(8) :: total_num=0, cnt
#if defined(half_precision)
    integer(2), allocatable :: v_buffer(:)
    type(my_real16) :: v123
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif

    open(wunit,file=f%val,form='unformatted',status='replace')
    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 /= p456) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = max(abs(2*j12-j3),abs(2*j45-j6)), min(2*j12+j3,2*j45+j6), 2

                      do t12 = 0, 1
                        do t45 = 0, 1
                          do t123 = 1, min(2*t12+1,2*t45+1), 2

                            total_num = total_num + 1
                            cntbf = cntbf + 1
                            v123 = this%GetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J123,T123)
                            v_buffer(cntbf) = v123
                            if(cntbf == numbf) then
                              do cnt = 1, cntbf
                                write(wunit) v_buffer(cnt)
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
        write(wunit) v_buffer(cnt)
      end do
    end if
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
    close(wunit)
  end subroutine write_scalar_operator_binary_isospin

  subroutine write_tensor_operator_binary_stream_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    integer :: wunit = 25
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, t456, j123, j456
    integer :: j12, j45, t12, t45
    integer :: numbf = 10000000, cntbf
    integer(8) :: total_num=0
#if defined(half_precision)
    integer(2), allocatable :: v_buffer(:)
    type(my_real16) :: v123
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif

    lab => this%ms
    sps => this%ms%sps
    open(wunit,file=f%val,form='unformatted',access='stream',status='replace')
    write(wunit) this%GetOpJ(), this%GetOpP(), this%GetOpT(), emax, e2max, e3max, sps%lmax
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 * p456 /= this%GetOpP()) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = abs(2*j12-j3), 2*j12+j3, 2
                      do j456 = abs(2*j45-j6), 2*j45+j6, 2
                        if(triag(j123,j456,2*this%GetOpJ())) cycle

                        do t12 = 0, 1
                          do t45 = 0, 1
                            do t123 = abs(2*t12-1), 2*t12+1, 2
                              do t456 = abs(2*t45-1), 2*t45+1, 2
                                if( triag(t123, t456, 2*this%GetOpT()) ) cycle

                                total_num = total_num + 1
                                cntbf = cntbf + 1
                                v123 = this%GetThBME(i1,i2,i3,J12,T12,J123,T123,i4,i5,i6,J45,T45,J456,T456)
                                v_buffer(cntbf) = v123
                                if(cntbf == numbf) then
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
      end do
    end do

    if(cntbf > 0) then
      write(wunit) v_buffer(1:cntbf)
    end if
    close(wunit)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', total_num
    deallocate(v_buffer)
  end subroutine write_tensor_operator_binary_stream_isospin

  subroutine write_tensor_operator_gzip_isospin(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(in) :: this
    type(str), intent(in) :: f
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
    integer :: i1, l1, j1
    integer :: i2, l2, j2
    integer :: i3, l3, j3
    integer :: i4, l4, j4
    integer :: i5, l5, j5
    integer :: i6, l6, j6
    integer :: i5max, i6max
    integer :: p123, p456, t123, t456, j123, j456
    integer :: j12, j45, t12, t45
    integer :: numbf = 100000000, cntbf, cnt, nrest
    integer(8) :: total_num=0
#if defined(half_precision) || defined(single_precision)
    real(4), allocatable :: v_buffer(:)
    real(4) :: v123
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
    real(8) :: v123
#endif
    type(c_ptr) :: fp, err
    character(256) :: header, buffer
    character(12) :: cfmt
    type(str) :: OpName

    fp = gzip_open(f%val, "wt")
    OpName = this%GetOpName()
    header = trim(OpName%val)
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NNN int.'
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
    write(buffer,"(7i4)") this%GetOpJ(), this%GetOpP(), this%GetOpT(), emax, e2max, e3max, sps%lmax
    err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
    allocate(v_buffer(numbf))
    total_num = 0; cntbf = 0
    do i1 = 1, sps%norbs
      l1 = sps%orb(i1)%l
      j1 = sps%orb(i1)%j
      do i2 = 1, i1
        l2 = sps%orb(i2)%l
        j2 = sps%orb(i2)%j
        if(sps%orb(i1)%e + sps%orb(i2)%e > e2max) cycle
        do i3 = 1, i2
          l3 = sps%orb(i3)%l
          j3 = sps%orb(i3)%j
          if(sps%orb(i1)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i2)%e + sps%orb(i3)%e > e2max) cycle
          if(sps%orb(i1)%e + sps%orb(i2)%e + sps%orb(i3)%e > e3max) cycle

          p123 = (-1) ** (l1+l2+l3)

          do i4 = 1, i1
            l4 = sps%orb(i4)%l
            j4 = sps%orb(i4)%j
            i5max = i4
            if(i1 == i4) i5max = i2

            do i5 = 1, i5max
              l5 = sps%orb(i5)%l
              j5 = sps%orb(i5)%j
              i6max = i5
              if(i1 == i4 .and. i2 == i5) i6max = i3
              if(sps%orb(i4)%e + sps%orb(i5)%e > e2max) cycle

              do i6 = 1, i6max
                l6 = sps%orb(i6)%l
                j6 = sps%orb(i6)%j
                if(sps%orb(i4)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i5)%e + sps%orb(i6)%e > e2max) cycle
                if(sps%orb(i4)%e + sps%orb(i5)%e + sps%orb(i6)%e > e3max) cycle

                p456 = (-1)**(l4+l5+l6)
                if(p123 * p456 /= this%GetOpP()) cycle

                do j12 = abs(j1-j2)/2, (j1+j2)/2
                  do j45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = abs(2*j12-j3), 2*j12+j3, 2
                      do j456 = abs(2*j45-j6), 2*j45+j6, 2
                        if(triag(j123,j456,2*this%GetOpJ())) cycle

                        do t12 = 0, 1
                          do t45 = 0, 1
                            do t123 = abs(2*t12-1), 2*t12+1, 2
                              do t456 = abs(2*t45-1), 2*t45+1, 2
                                if( triag(t123, t456, 2*this%GetOpT()) ) cycle

                                total_num = total_num + 1
                                cntbf = cntbf + 1
                                v123 = this%GetThBME(i1,i2,i3,J12,T12,J123,T123,i4,i5,i6,J45,T45,J456,T456)
                                v_buffer(cntbf) = v123
                                if(cntbf == numbf) then
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
      end do
    end do

    if(cntbf > 0) then
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
  end subroutine write_tensor_operator_gzip_isospin

  subroutine ReadThreeBodyLabOp(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling WriteThreeBodyLabOp"
      return
    end if

    ti = omp_get_wtime()
    if(.not. this%reduced_me() ) call read_scalar_operator_isospin(this, f)
    if(      this%reduced_me() ) call read_tensor_operator_isospin(this, f)
    call timer%Add(s%str('Reading from file'), omp_get_wtime() - ti)
  end subroutine ReadThreeBodyLabOp

  subroutine read_scalar_operator_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f, s%str(".txt"))) then
      call read_scalar_operator_ascii_isospin(this, f)
      return
    end if

    if(s%find(f, s%str(".me3j.gz"))) then
      call read_scalar_operator_gzip_me3j_isospin(this, f)
      return
    end if

    if(s%find(f, s%str(".me3j"))) then
      call read_scalar_operator_ascii_me3j_isospin(this,f)
      return
    end if

    if(s%find(f, s%str(".stream.bin"))) then
      call read_scalar_operator_binary_stream_isospin(this, f)
      return
    end if

    if(s%find(f, s%str(".bin"))) then
      call read_scalar_operator_binary_isospin(this, f)
      return
    end if

    write(*,"(a)") "# Unknown input file format"
    stop
  end subroutine read_scalar_operator_isospin

  subroutine read_tensor_operator_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    if(s%find(f, s%str(".bin.stream"))) then
      call read_tensor_operator_binary_stream_isospin(this, f)
      return
    end if

    if(s%find(f, s%str(".me3j.gz"))) then
      call read_tensor_operator_gzip_me3j_isospin(this, f)
      return
    end if

    write(*,"(a)") "# Unknown input file format"
    stop
  end subroutine read_tensor_operator_isospin

  subroutine read_scalar_operator_ascii_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    type(str), intent(in) :: f
    integer :: runit = 25
    integer :: emax, e2max, e3max
#if defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    type(my_real16), allocatable :: v_tmp(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    character(256) :: header
    integer(8) :: nelms, line

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_scalar_3bme(sps, e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))
    open(runit, file=f%val, status='replace', action='write')
    read(runit,*) header
    do line = 1, nelms
      read(runit,*) v_buffer(line)
    end do
    close(runit)
#if defined(half_precision)
    allocate(v_tmp(nelms))
    v_tmp = v_buffer
    call store_scalar_3bme(this,v_tmp,e2max,e3max)
    deallocate(v_tmp)
#elif defined(single_precision) || defined(double_precision)
    call store_scalar_3bme(this,v_buffer,e2max,e3max)
#endif
    deallocate(v_buffer)
  end subroutine read_scalar_operator_ascii_isospin

  subroutine read_scalar_operator_gzip_me3j_isospin(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    type(str), intent(in) :: f
    integer :: emax, e2max, e3max
#if defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    type(my_real16), allocatable :: v_tmp(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    character(256) :: header, buffer
    integer(8) :: nelms, line
    type(c_ptr) :: fp, err

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_scalar_3bme(sps, e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    fp = gzip_open(f%val, "rt")
    err = gzip_readline(fp, header, len(header))
    do line = 1, nelms/10
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v_buffer((line-1)*10+1 : line*10)
    end do

    ! basically, this is not needed
    if(nelms - (nelms/10) * 10 > 0) then
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v_buffer((nelms/10)*10+1 : nelms)
    end if
    err = gzip_close(fp)

#if defined(half_precision)
    allocate(v_tmp(nelms))
    v_tmp = v_buffer
    call store_scalar_3bme(this,v_tmp,e2max,e3max)
    deallocate(v_tmp)
#elif defined(single_precision) || defined(double_precision)
    call store_scalar_3bme(this,v_buffer,e2max,e3max)
#endif
    deallocate(v_buffer)
  end subroutine read_scalar_operator_gzip_me3j_isospin

  subroutine read_scalar_operator_ascii_me3j_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    type(str), intent(in) :: f
    integer :: runit = 25
    integer :: emax, e2max, e3max
#if defined(half_precision)
    real(4), allocatable :: v_buffer(:)
    type(my_real16), allocatable :: v_tmp(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    integer(8) :: nelms, line
    character(256) :: header

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_scalar_3bme(sps, e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    open(runit, file=f%val, action='read')
    read(runit,*) header
    do line = 1, nelms/10
      read(runit,*) v_buffer((line-1)*10+1 : line*10)
    end do

    ! basically, this is not needed
    if(nelms - (nelms/10) * 10 > 0) then
      read(runit,*) v_buffer((nelms/10)*10+1 : nelms)
    end if
    close(runit)
#if defined(half_precision)
    allocate(v_tmp(nelms))
    v_tmp = v_buffer
    call store_scalar_3bme(this,v_tmp,e2max,e3max)
    deallocate(v_tmp)
#elif defined(single_precision) || defined(double_precision)
    call store_scalar_3bme(this,v_buffer,e2max,e3max)
#endif
    deallocate(v_buffer)
  end subroutine read_scalar_operator_ascii_me3j_isospin

  subroutine read_scalar_operator_binary_stream_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    type(str), intent(in) :: f
    integer :: runit = 25
    integer :: emax, e2max, e3max
#if defined(half_precision)
    type(my_real16), allocatable :: v_buffer(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    integer(8) :: nelms

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_scalar_3bme(sps, e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    open(runit, file=f%val, form="unformatted", access="stream", action='read')
    read(runit) v_buffer
    close(runit)

    call store_scalar_3bme(this,v_buffer,e2max,e3max)
    deallocate(v_buffer)
  end subroutine read_scalar_operator_binary_stream_isospin

  subroutine read_scalar_operator_binary_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    type(str), intent(in) :: f
    integer :: runit = 25
    integer :: emax, e2max, e3max
#if defined(half_precision)
    type(my_real16), allocatable :: v_buffer(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    integer(8) :: nelms, line

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_scalar_3bme(sps, e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    open(runit, file=f%val, form="unformatted", action='read')
    do line = 1, nelms
      read(runit) v_buffer(line)
    end do
    close(runit)

    call store_scalar_3bme(this,v_buffer,e2max,e3max)
    deallocate(v_buffer)
  end subroutine read_scalar_operator_binary_isospin

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

  subroutine store_scalar_3bme(this,v,e2max,e3max)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(ThreeBodyLabIsoSpace), pointer :: ms
    type(OrbitsIsospin), pointer :: spsf
#if defined(half_precision)
    type(my_real16), intent(in) :: v(:)
    real(4) :: v123
#elif defined(single_precision)
    real(4), intent(in) :: v(:)
#elif defined(double_precision)
    real(8), intent(in) :: v(:)
#endif
    integer, intent(in) :: e2max, e3max
    integer(8) :: cnt
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J, T, P123, P456

    ms => this%ms
    spsf => ms%sps

    cnt = 0
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
                            cnt = cnt + 1

                            if(e1 > ms%GetEmax()) cycle
                            if(e2 > ms%GetEmax()) cycle
                            if(e3 > ms%GetEmax()) cycle

                            if(e4 > ms%GetEmax()) cycle
                            if(e5 > ms%GetEmax()) cycle
                            if(e6 > ms%GetEmax()) cycle

                            if(e1 + e2 > ms%GetE2max()) cycle
                            if(e2 + e3 > ms%GetE2max()) cycle
                            if(e3 + e1 > ms%GetE2max()) cycle

                            if(e4 + e5 > ms%GetE2max()) cycle
                            if(e5 + e6 > ms%GetE2max()) cycle
                            if(e6 + e4 > ms%GetE2max()) cycle

                            if(e1 + e2 + e3 > ms%GetE3max()) cycle
                            if(e4 + e5 + e6 > ms%GetE3max()) cycle
#if defined(half_precision)
                            v123 = v(cnt)
                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v123) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v123) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            call this%SetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J,T,v123)
#elif defined(single_precision) || defined(double_precision)
                            if(i1==i2 .and. mod(J12+T12,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            if(i4==i5 .and. mod(J45+T45,2)==0) then
                              if(abs(v(cnt)) > 1.d-6) then
                                write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                              end if
                              cycle
                            end if

                            call this%SetThBME(i1,i2,i3,J12,T12,i4,i5,i6,J45,T45,J,T,v(cnt))
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
  end subroutine store_scalar_3bme

  subroutine read_tensor_operator_binary_stream_isospin(this, f)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: runit = 25
    integer :: emax, e2max, e3max
#if defined(half_precision)
    !integer(2), allocatable :: v_buffer(:)
    type(my_real16), allocatable :: v_buffer(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    integer(8) :: nelms

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_tensor_3bme(sps, this%GetOpJ(), this%GetOpP(), this%GetOpT(), e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    open(runit, file=f%val, form="unformatted", access="stream", action='read')
    read(runit) v_buffer
    close(runit)

    call store_tensor_3bme(this,v_buffer,e2max,e3max)
    deallocate(v_buffer)
  end subroutine read_tensor_operator_binary_stream_isospin

  subroutine read_tensor_operator_gzip_me3j_isospin(this, f)
    use, intrinsic :: iso_c_binding
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
    type(str), intent(in) :: f
    type(ThreeBodyLabIsoSpace), pointer :: lab
    type(OrbitsIsospin), pointer :: sps
    integer :: emax, e2max, e3max
#if defined(half_precision)
    type(my_real16), allocatable :: v_tmp(:)
    real(4), allocatable :: v_buffer(:)
#elif defined(single_precision)
    real(4), allocatable :: v_buffer(:)
#elif defined(double_precision)
    real(8), allocatable :: v_buffer(:)
#endif
    character(256) :: header, buffer
    integer(8) :: nelms, line
    type(c_ptr) :: fp, err

    lab => this%ms
    sps => this%ms%sps
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    e3max = lab%GetE3max()
    nelms = count_tensor_3bme(sps, this%GetOpJ(), this%GetOpP(), this%GetOpT(), e2max, e3max)
    write(*,'(a,i13)') 'Number of three-body matrix elements is ', nelms
    allocate(v_buffer(nelms))

    fp = gzip_open(f%val, "rt")
    err = gzip_readline(fp, header, len(header))
    err = gzip_readline(fp, header, len(header))
    do line = 1, nelms/10
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v_buffer((line-1)*10+1 : line*10)
    end do
    ! basically, this is not needed
    if(nelms - (nelms/10) * 10 > 0) then
      err = gzip_readline(fp, buffer, len(buffer))
      read(buffer,*) v_buffer((nelms/10)*10+1 : nelms)
    end if

    err = gzip_close(fp)
#if defined(half_precision)
    allocate(v_tmp(nelms))
    v_tmp = v_buffer
    call store_tensor_3bme(this,v_tmp,e2max,e3max)
    deallocate(v_tmp)
#elif defined(single_precision) || defined(double_precision)
    call store_tensor_3bme(this,v_buffer,e2max,e3max)
#endif
    deallocate(v_buffer)
  end subroutine read_tensor_operator_gzip_me3j_isospin

  function count_tensor_3bme(spsf, rankJ, rankP, rankT, e2max, e3max) result(r)
    type(OrbitsIsospin), intent(in) :: spsf
    integer, intent(in) :: e2max, e3max
    integer, intent(in) :: rankJ, rankP, rankT
    integer(8) :: r
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer :: J12, T12, J45, T45, J123, J456, T123, T456, P123, P456

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

                if(P123 * P456 /= rankP ) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = abs(2*j12-j3), 2*j12+j3, 2
                      do j456 = abs(2*j45-j6), 2*j45+j6, 2
                        if(triag(j123,j456,2*rankJ)) cycle

                        do T12 = 0, 1
                          do T45 = 0, 1
                            do t123 = abs(2*t12-1), 2*t12+1, 2
                              do t456 = abs(2*t45-1), 2*t45+1, 2
                                if( triag(t123, t456, 2*rankT) ) cycle
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
      end do
    end do
  end function count_tensor_3bme

  subroutine store_tensor_3bme(this, v, e2max, e3max)
    class(MACRO_ADD_SUFX(ThreeBodyLabOpIso)), intent(inout) :: this
#if defined(half_precision)
    type(my_real16), intent(in) :: v(:)
    real(4) :: v123
#elif defined(single_precision)
    real(4), intent(in) :: v(:)
#elif defined(double_precision)
    real(8), intent(in) :: v(:)
#endif
    integer, intent(in) :: e2max, e3max
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5, i5max
    integer :: i6, l6, j6, e6, i6max
    integer(8) :: cnt
    integer :: J12, T12, J45, T45, J123, J456, T123, T456, P123, P456
    type(ThreeBodyLabIsoSpace), pointer :: ms
    type(OrbitsIsospin), pointer :: spsf

    ms => this%ms
    spsf => ms%sps
    cnt = 0
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

                if(P123 * P456 /= this%GetOpP() ) cycle
                do J12 = abs(j1-j2)/2, (j1+j2)/2
                  do J45 = abs(j4-j5)/2, (j4+j5)/2
                    do j123 = abs(2*j12-j3), 2*j12+j3, 2
                      do j456 = abs(2*j45-j6), 2*j45+j6, 2
                        if(triag(j123,j456,2*this%GetOpJ())) cycle

                        do T12 = 0, 1
                          do T45 = 0, 1
                            do t123 = abs(2*t12-1), 2*t12+1, 2
                              do t456 = abs(2*t45-1), 2*t45+1, 2
                                if( triag(t123, t456, 2*this%GetOpT()) ) cycle
                                cnt = cnt + 1

                                if(e1 > ms%GetEmax()) cycle
                                if(e2 > ms%GetEmax()) cycle
                                if(e3 > ms%GetEmax()) cycle

                                if(e4 > ms%GetEmax()) cycle
                                if(e5 > ms%GetEmax()) cycle
                                if(e6 > ms%GetEmax()) cycle

                                if(e1 + e2 > ms%GetE2max()) cycle
                                if(e2 + e3 > ms%GetE2max()) cycle
                                if(e3 + e1 > ms%GetE2max()) cycle

                                if(e4 + e5 > ms%GetE2max()) cycle
                                if(e5 + e6 > ms%GetE2max()) cycle
                                if(e6 + e4 > ms%GetE2max()) cycle

                                if(e1 + e2 + e3 > ms%GetE3max()) cycle
                                if(e4 + e5 + e6 > ms%GetE3max()) cycle

#if defined(half_precision)
                                v123 = v(cnt)
                                if(i1==i2 .and. mod(J12+T12,2)==0) then
                                  if(abs(v123) > 1.d-6) then
                                    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                                  end if
                                  cycle
                                end if

                                if(i4==i5 .and. mod(J45+T45,2)==0) then
                                  if(abs(v123) > 1.d-6) then
                                    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                                  end if
                                  cycle
                                end if
                                call this%SetThBME( i1,i2,i3,J12,T12,J123,T123, i4,i5,i6,J45,T45,J456,T456, v123 )
#elif defined(single_precision) || defined(double_precision)
                                if(i1==i2 .and. mod(J12+T12,2)==0) then
                                  if(abs(v(cnt)) > 1.d-6) then
                                    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                                  end if
                                  cycle
                                end if

                                if(i4==i5 .and. mod(J45+T45,2)==0) then
                                  if(abs(v(cnt)) > 1.d-6) then
                                    write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                                  end if
                                  cycle
                                end if
                                call this%SetThBME( i1,i2,i3,J12,T12,J123,T123, i4,i5,i6,J45,T45,J456,T456, v(cnt) )
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
      end do
    end do
  end subroutine store_tensor_3bme
