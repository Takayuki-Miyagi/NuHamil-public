module TwoBodyLabOps
  use omp_lib
  use ClassSys
  use Profiler, only: timer
  use LinAlgLib
  use TwoBodyRelativeSpace
  use TwoBodyLabSpacePN
  use SingleParticleState
  use OperatorDefinitions
  implicit none

  public :: TwoBodyLabOp
  private :: TwoBodyLabOpChan

  ! Methods
  private :: InitTwoBodyLabOp ! constructor
  private :: InitTwoBodyLabOpFromString ! constructor (more convenient)
  private :: FinTwoBodyLabOp
  ! Talmi-Moshinsky transformation
  private :: TalmiMoshinskyTransformation
  private :: TMTransScalar
  private :: TMTransTensor
  private :: GetFileNameTwoBodyLabOp ! Getter for file name
  private :: SetTwBME_tensor
  private :: SetTwBME_scalar
  private :: AddToTwBME_tensor
  private :: AddToTwBME_scalar
  ! Writting to file name
  ! interfaces
  private :: WriteOperator
  private :: ReadOperator
  private :: write_scalar_operator
  private :: write_general_operator
  private :: read_scalar_operator
  ! actual routines
  private :: write_scalar_operator_ascii_myg
  private :: write_scalar_operator_binary_myg
  private :: write_general_operator_ascii
  private :: write_scalar_operator_ascii_snt
  private :: write_scalar_operator_kshell_snt
  private :: write_scalar_operator_kshell_A_snt
  private :: write_general_operator_gzip
  private :: write_scalar_operator_ascii_me2j
  private :: write_scalar_operator_gzip_me2j
  private :: write_scalar_operator_binary_me2j
  private :: write_scalar_operator_ascii_mfd

  private :: read_scalar_operator_ascii_myg
  private :: read_scalar_operator_binary_myg
  private :: read_scalar_operator_ascii_snt
  private :: read_scalar_operator_kshell_snt
  private :: read_scalar_operator_kshell_A_snt
  !private :: read_scalar_operator_ascii_me2j
  !private :: read_scalar_operator_gzip_me2j
  !private :: read_scalar_operator_binary_me2j

  private :: p_dot_p
  private :: r_dot_r
  private :: red_nab_j
  private :: red_r_j
  private :: GetMemory
  private :: ReducedToNonReduced
  private :: NonReducedToReduced
  private :: SingularValueDecomposition
  private :: SpinTensorDecomposition

  type, extends(DMat) :: TwoBodyLabOpChan
    logical :: Zero = .true.
  contains
    procedure :: GetScalar
    procedure :: GetTensor
  end type TwoBodyLabOpChan

  type, extends(OperatorDef) :: TwoBodyLabOp
    type(TwoBodyLabOpChan), allocatable :: MatCh(:,:)
    type(DMat) :: OneBody
    type(TwoBodyLabPNSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyLabOpFromString
    procedure :: InitTwoBodyLabOp
    procedure :: CopyTwoBodyLabOp
    procedure :: SumTwoBodyLabOp
    procedure :: SubtractTwoBodyLabOp
    procedure :: ScaleTwoBodyLabOp
    procedure :: GetMemory
    procedure :: Truncate
    generic :: init => InitTwoBodyLabOpFromString, InitTwoBodyLabOp
    procedure :: fin => FinTwoBodyLabOp
    procedure :: writef => WriteOperator
    procedure :: TalmiMoshinskyTransformation
    procedure :: TMTransCMTensor
    procedure :: TMTransCMTensor_explicit
    procedure :: GetFile => GetFileNameTwoBodyLabOp
    procedure :: readf => ReadOperator
    procedure :: GetTBME
    procedure :: SetTwBME_tensor
    procedure :: SetTwBME_scalar
    procedure :: AddToTwBME_tensor
    procedure :: AddToTwBME_scalar
    procedure :: ReducedToNonReduced
    procedure :: NonReducedToReduced
    procedure :: SingularValueDecomposition
    procedure :: SpinTensorDecomposition
    generic :: TMTrans => TalmiMoshinskyTransformation, TMTransCMTensor
    generic :: SetTBME => SetTwBME_tensor, SetTwBME_scalar
    generic :: AddToTBME => AddToTwBME_tensor, AddToTwBME_scalar
    generic :: SVD => SingularValueDecomposition
    generic :: assignment(=) => CopyTwoBodyLabOp
    generic :: operator(+) => SumTwoBodyLabOp
    generic :: operator(-) => SubtractTwoBodyLabOp
    generic :: operator(*) => ScaleTwoBodyLabOp
  end type TwoBodyLabOp
contains

  subroutine CopyTwoBodyLabOp(a, b)
    class(TwoBodyLabOp), intent(inout) :: a
    type(TwoBodyLabOp), intent(in) :: b
    integer :: ichb, ichk, n
    if(allocated(a%MatCh)) call a%fin()
    call a%OneBody%fin()
    n = size(b%MatCh, 1)
    call a%CopyOperatorDef(b%OperatorDef)
    a%ms => b%ms
    a%OneBody = b%OneBody
    allocate(a%MatCh(n,n))
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%Zero) then
          a%MatCh(ichb,ichk)%zero = .true.
          cycle
        end if
        a%MatCh(ichb,ichk)%Zero = .false.
        a%MatCh(ichb,ichk)%DMat = b%MatCh(ichb, ichk)%DMat
      end do
    end do
    a%is_init = b%is_init
  end subroutine CopyTwoBodyLabOp

  function SumTwoBodyLabOp(a, b) result(c)
    class(TwoBodyLabOp), intent(in) :: a, b
    type(TwoBodyLabOp) :: c
    integer :: ichb, ichk, n
    c = a
    n = size(b%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%Zero) then
          c%MatCh(ichb,ichk)%zero = .true.
          cycle
        end if
        c%MatCh(ichb,ichk)%Zero = .false.
        c%MatCh(ichb,ichk)%DMat = a%MatCh(ichb, ichk)%DMat + b%MatCh(ichb, ichk)%DMat
      end do
    end do
  end function SumTwoBodyLabOp

  function SubtractTwoBodyLabOp(a, b) result(c)
    class(TwoBodyLabOp), intent(in) :: a, b
    type(TwoBodyLabOp) :: c
    integer :: ichb, ichk, n
    c = a
    n = size(b%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%Zero) then
          c%MatCh(ichb,ichk)%zero = .true.
          cycle
        end if
        c%MatCh(ichb,ichk)%Zero = .false.
        c%MatCh(ichb,ichk)%DMat = a%MatCh(ichb, ichk)%DMat - b%MatCh(ichb, ichk)%DMat
      end do
    end do
  end function SubtractTwoBodyLabOp

  function ScaleTwoBodyLabOp(a, b) result(c)
    class(TwoBodyLabOp), intent(in) :: a
    type(TwoBodyLabOp) :: c
    real(8), intent(in) :: b
    integer :: ichb, ichk, n
    c = a
    n = size(a%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(a%MatCh(ichb,ichk)%Zero) then
          c%MatCh(ichb,ichk)%zero = .true.
          cycle
        end if
        c%MatCh(ichb,ichk)%Zero = .false.
        c%MatCh(ichb,ichk)%DMat = a%MatCh(ichb, ichk)%DMat * b
      end do
    end do
  end function ScaleTwoBodyLabOp

  subroutine FinTwoBodyLabOp(this)
    class(TwoBodyLabOp), intent(inout) :: this
    integer :: chbra, chket, n
    n = size(this%MatCh, 1)
    do chbra = 1, n
      do chket = 1, n
        call this%MatCh(chbra,chket)%fin()
      end do
    end do
    call this%OneBody%fin()
    call this%FinOperatorDef()
    deallocate(this%MatCh)
    this%ms => null()
    this%is_init = .false.
  end subroutine FinTwoBodyLabOp

  subroutine InitTwoBodyLabOpFromString(this, two, oprtr)
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), intent(in) :: two
    type(str), intent(in) :: oprtr

    call this%InitOpDef(oprtr, .true.)
    call this%init(two, this%GetOpJ(), this%GetOpP(), this%GetOpZ())
  end subroutine InitTwoBodyLabOpFromString

  function GetMemory(this) result(r)
    class(TwoBodyLabOp), intent(in) :: this
    real(8) :: r
    integer :: chbra, chket
    r = 0.d0
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%zero) cycle
        r = r + 8.d0 * dble(size(this%MatCh(chbra,chket)%m)) / 1024.d0**3
      end do
    end do
  end function GetMemory

  function Truncate( this, ms_new ) result(op)
    class(TwoBodyLabOp), intent(in) :: this
    type(TwoBodyLabPNSpace), intent(in), target :: ms_new
    type(TwoBodyLabOp) :: op
    type(TwoBodyLabPNChan), pointer :: chbra, chket
    integer :: ichbra, jbra, bra, a, b
    integer :: ichket, jket, ket, c, d
    real(8) :: time
    type(sys) :: s

    time = omp_get_wtime()
    call op%init( ms_new, this%GetOpName() )

    op%OneBody%m(:,:) = this%OneBody%m(:ms_new%sps%norbs, :ms_new%sps%norbs)

    !$omp parallel
    !$omp do private( ichbra, ichket, chbra, chket, &
    !$omp &  jbra, jket, bra, a, b, ket, c, d ) schedule(dynamic)
    do ichbra = 1, ms_new%GetNumberChannels()
      do ichket = 1, ms_new%GetNumberChannels()
        if( op%MatCh(ichbra,ichket)%Zero ) cycle
        chbra => ms_new%GetChannel( ichbra )
        chket => ms_new%GetChannel( ichket )
        jbra = chbra%GetJ()
        jket = chket%GetJ()

        do bra = 1, chbra%GetNumberStates()
          a = chbra%n2label1( bra )
          b = chbra%n2label2( bra )
          do ket = 1, chket%GetNumberStates()
            c = chket%n2label1( ket )
            d = chket%n2label2( ket )
            op%MatCh(ichbra,ichket)%m(bra,ket) = this%GetTBME(a,b,c,d,Jbra,Jket)
          end do
        end do

      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%Add(s%str("TwoBodyLabOp: Truncate"), omp_get_wtime()-time)
  end function Truncate

  subroutine SingularValueDecomposition(this, svd_rank)
    class(TwoBodyLabOp), intent(inout) :: this
    integer, intent(in) :: svd_rank
    integer :: ichb, ichk, ndim
    type(DMat) :: A, tmp
    type(DVec) :: Sigma
    type(DSingularValueDecomposition) :: sol
    type(TwoBodyLabPNChan), pointer :: chbra, chket
    type(str) :: fn
    if(svd_rank==-1) return
    do ichb = 1, this%ms%GetNumberChannels()
      do ichk = 1, ichb
        if(this%MatCh(ichb, ichk)%Zero) cycle
        chbra => this%ms%GetChannel(ichb)
        chket => this%ms%GetChannel(ichk)
        A = this%MatCh(ichb,ichk)%DMat
        call sol%init(A)
        call sol%svd(A)
        call Sigma%zeros(sol%Sigma%n_size)
        if(chbra%GetJ() == 1 .and. chket%GetJ()==1 .and. &
            & chbra%GetParity()==1 .and. chket%GetParity()==1 .and. &
            & chbra%GetZ()==0 .and. chket%GetZ()==0 ) then
          fn = "svalues_J1P+Z0.bin"
          open(15, file=fn%val, form="unformatted", access="stream")
          write(15) sol%sigma%v(:)
          close(15)
        end if
        ndim = min(svd_rank, Sigma%n_size)
        Sigma%v(:ndim) = sol%Sigma%v(:ndim)
        call tmp%DiagMat(Sigma)
        this%MatCh(ichb,ichk)%DMat = sol%U * tmp * sol%V
        call sol%fin()
      end do
    end do
  end subroutine SingularValueDecomposition

  function SpinTensorDecomposition(this, rank) result(op)
    use MyLibrary, only: triag, sjs, snj
    class(TwoBodyLabOp), intent(inout) :: this
    integer, intent(in) :: rank
    type(TwoBodyLabOp) :: op
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: tbc
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ich, Lab, Sab, Lcd, Scd, JJ, J, a, b, c, d
    integer :: jaa, jbb, jcc, jdd, ibra, iket, aa, bb, cc, dd
    real(8) :: sum1, sum2, sum3, norm_1, norm_2

    ms => this%ms
    sps => ms%sps
    call op%init(ms, this%GetOpName())
    do ich = 1, ms%GetNumberChannels()
      op%MatCh(ich,ich) = this%MatCh(ich,ich)
    end do
    if( this%GetOpJ() /= 0 ) then
      write(*,*) "SpinTensorDecompsition is impolemented only for scalar operator"
      return
    end if
    if( rank > 2) then
      write(*,*) "Two-body spin-tensor decomposition is valid only up to rank=2"
      return
    end if

    do ich = 1, ms%GetNumberChannels()
      tbc => ms%GetChannel(ich)
      J = tbc%GetJ()

      do ibra = 1, tbc%GetNumberStates()
        a = tbc%n2label1(ibra)
        b = tbc%n2label2(ibra)
        oa => sps%GetOrbit(a)
        ob => sps%GetOrbit(b)
        do iket = 1, ibra
          c = tbc%n2label1(iket)
          d = tbc%n2label2(iket)
          oc => sps%GetOrbit(c)
          od => sps%GetOrbit(d)
          norm_1 = 1.d0
          if(a==b) norm_1 = norm_1 * sqrt(2.d0)
          if(c==d) norm_1 = norm_1 * sqrt(2.d0)

          sum3 = 0.d0
          do Sab = 0, 1
            do Lab = max(abs(J-Sab), abs(oa%l-ob%l)), min(J+Sab, oa%l+ob%l)

              do Scd = 0, 1
                do Lcd = max(abs(J-Scd), abs(oc%l-od%l)), min(J+Scd, oc%l+od%l)

                  if(triag(Sab,Scd,rank)) cycle
                  if(triag(Lab,Lcd,rank)) cycle

                  sum2 = 0.d0
                  do JJ = max(abs(Lab-Sab), abs(Lcd-Scd)), min(Lab+Sab, Lcd+Scd)

                    sum1 = 0.d0
                    do jaa = abs(2*oa%l-1), 2*oa%l+1, 2
                      do jbb = abs(2*ob%l-1), 2*ob%l+1, 2
                        do jcc = abs(2*oc%l-1), 2*oc%l+1, 2
                          do jdd = abs(2*od%l-1), 2*od%l+1, 2
                            if(triag(jaa, jbb, 2*JJ)) cycle
                            if(triag(jcc, jdd, 2*JJ)) cycle
                            aa = sps%nljz2idx(oa%n, oa%l, jaa, oa%z)
                            bb = sps%nljz2idx(ob%n, ob%l, jbb, ob%z)
                            cc = sps%nljz2idx(oc%n, oc%l, jcc, oc%z)
                            dd = sps%nljz2idx(od%n, od%l, jdd, od%z)
                            norm_2 = 1.d0
                            if(aa==bb) norm_2 = norm_2 * sqrt(2.d0)
                            if(cc==dd) norm_2 = norm_2 * sqrt(2.d0)
                            sum1 = sum1 + this%GetTBME(aa,bb,cc,dd,JJ,JJ) * norm_2 * &
                                & ls_to_jj(oa%l, jaa, ob%l, jbb, Lab, Sab, JJ) * &
                                & ls_to_jj(oc%l, jcc, od%l, jdd, Lcd, Scd, JJ)
                          end do
                        end do
                      end do
                    end do
                    sum2 = sum2 + sum1 * dble(2*JJ+1) * (-1.d0)**(JJ+J) * &
                        & sjs(2*Lab, 2*Sab, 2*JJ, 2*Scd, 2*Lcd, 2*rank)

                  end do

                  sum3 = sum3 + sum2 * dble(2*rank+1) * &
                      & sjs(2*Lab, 2*Sab, 2*J, 2*Scd, 2*Lcd, 2*rank) * &
                      & ls_to_jj(oa%l, oa%j, ob%l, ob%j, Lab, Sab, J) * &
                      & ls_to_jj(oc%l, oc%j, od%l, od%j, Lcd, Scd, J)

                end do
              end do

            end do
          end do
          call op%SetTBME(a, b, c, d, J, sum3/norm_1)
        end do
      end do

    end do
  contains
    function ls_to_jj(la, ja, lb, jb, Lab, Sab, J) result(res)
      integer, intent(in) :: la, ja, lb, jb, Lab, Sab, J
      real(8) :: res
      res = sqrt(dble((ja+1)*(jb+1)*(2*Lab+1)*(2*Sab+1))) * &
          & snj(2*la, 1, ja, 2*lb, 1, jb, 2*Lab, 2*Sab, 2*J)
    end function ls_to_jj
  end function SpinTensorDecomposition

  subroutine ReducedToNonReduced( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: ichbra, ichket, J

    if( this%GetOpJ() /= 0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
      return
    end if
    if( .not. this%reduced_me() ) return
    ms => this%ms
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ms%GetNumberChannels()
        if( this%MatCh(ichbra,ichket)%zero) cycle
        J = ms%jpz(ichket)%GetJ()
        this%MatCh(ichbra,ichket)%DMat = this%MatCh(ichbra,ichket)%DMat / sqrt(dble(2*J+1))
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReduced( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: ichbra, ichket, J

    if( this%GetOpJ() /= 0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
      return
    end if
    if( this%reduced_me() ) return
    ms => this%ms
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ms%GetNumberChannels()
        if( this%MatCh(ichbra,ichket)%zero) cycle
        J = ms%jpz(ichket)%GetJ()
        this%MatCh(ichbra,ichket)%DMat = this%MatCh(ichbra,ichket)%DMat * sqrt(dble(2*J+1))
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NonReducedToReduced

  subroutine InitTwoBodyLabOp(this, two, jr, pr, tzr)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), target, intent(in) :: two
    integer, intent(in) :: jr, pr, tzr
    integer :: chbra, chket
    integer :: nb, jb, pb, tzb
    integer :: nk, jk, pk, tzk

    if(allocated(this%MatCh)) call this%fin()
    this%ms => two
    call this%InitOpDef(.true., jr, pr, tzr)
    allocate(this%MatCh(two%GetNumberChannels(), two%GetNumberChannels()))
    do chbra = 1, two%GetNumberChannels()
      nb = this%ms%jpz(chbra)%GetNumberStates()
      jb = this%ms%jpz(chbra)%GetJ()
      pb = this%ms%jpz(chbra)%GetParity()
      tzb= this%ms%jpz(chbra)%GetZ()
      !do chket = 1, two%GetNumberChannels()
      do chket = 1, chbra
        nk = this%ms%jpz(chket)%GetNumberStates()
        jk = this%ms%jpz(chket)%GetJ()
        pk = this%ms%jpz(chket)%GetParity()
        tzk= this%ms%jpz(chket)%GetZ()
        if(triag(jb, jk, jr)) cycle
        if(pb * pr /= pk) cycle
        if(abs(tzk-tzb) /= tzr) cycle
        this%MatCh(chbra,chket)%Zero = .false.
        call this%MatCh(chbra, chket)%zeros(nb,nk)
      end do
    end do
    call this%OneBody%zeros( this%ms%sps%norbs, this%ms%sps%norbs )
    this%is_init = .true.
  end subroutine InitTwoBodyLabOp

  subroutine TalmiMoshinskyTransformation(this,rel2lab,oprel,svd_rank_tcoef)
    use TwoBodyTransCoef, only: TransRel2LabSpace
    use TwoBodyRelOps
    class(TwoBodyLabOp), intent(inout) :: this
    type(TransRel2LabSpace), intent(in) :: rel2lab
    type(TwoBodyRelOp), intent(in) :: oprel
    integer, intent(in), optional :: svd_rank_tcoef

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyLabOp before calling TalmiMoshinskyTransofrmation"
      return
    end if

    if(.not. rel2lab%is_constructed) then
      write(*,*) "Initialize TransRel2LabSapce before calling TalmiMoshinskyTransofrmation"
      return
    end if

    if(.not. oprel%is_init) then
      write(*,*) "Initialize TwoBodyRelOperator before calling TalmiMoshinskyTransofrmation"
      return
    end if

    if( .not. this%reduced_me() .and. .not. oprel%reduced_me() ) call TMTransScalar(this,rel2lab,oprel,svd_rank_tcoef)
    if(       this%reduced_me() .and.       oprel%reduced_me() ) call TMTransTensor(this,rel2lab,oprel,svd_rank_tcoef)
    if(       this%reduced_me() .and. .not. oprel%reduced_me() ) then
      write(*,*) "Error at", __LINE__, " in ", __FILE__
      stop
    end if
    if( .not. this%reduced_me() .and.       oprel%reduced_me() ) then
      write(*,*) "Error at", __LINE__, " in ", __FILE__
      stop
    end if
  end subroutine TalmiMoshinskyTransformation

  function GetTBME(this,a,b,c,d,Jab,Jcd) result(r)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(in) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer, intent(in) :: a, b, c, d, Jab, Jcd
    integer :: Pab, Pcd, Tzab, Tzcd
    integer :: chbra, chket
    real(8) :: r

    r = 0.d0
    ms => this%ms
    sps => this%ms%sps
    Pab = (-1)**(sps%orb(a)%l + sps%orb(b)%l)
    Pcd = (-1)**(sps%orb(c)%l + sps%orb(d)%l)

    Tzab = (sps%orb(a)%z + sps%orb(b)%z) / 2
    Tzcd = (sps%orb(c)%z + sps%orb(d)%z) / 2

    if(triag(Jab,Jcd,this%GetOpJ())) then
      write(*,*) 'Error, at ', __LINE__, " in ", __FILE__
      stop
    end if
    if(Pab * Pcd * this%GetOpP() /= 1) then
      write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      stop
    end if
    if(abs(Tzab-Tzcd) /= this%GetOpZ() ) return

    chbra = ms%jpz2idx(Jab,Pab,Tzab)
    chket = ms%jpz2idx(Jcd,Pcd,Tzcd)
    if(chbra * chket == 0) return

    if( .not. this%reduced_me() ) r = this%MatCh(chbra,chket)%GetScalar(ms%jpz(chket),a,b,c,d)
    if(       this%reduced_me() ) then
      if( chbra < chket ) then
        r = this%MatCh(chket,chbra)%GetTensor(ms%jpz(chket),ms%jpz(chbra),c,d,a,b) * (-1.d0)**(Jcd-Jab)
        return
      end if
      r = this%MatCh(chbra,chket)%GetTensor(ms%jpz(chbra),ms%jpz(chket),a,b,c,d)
    end if
  end function GetTBME

  function GetScalar(this,ms,a,b,c,d) result(r)
    class(TwoBodyLabOpChan), intent(in) :: this
    type(TwoBodyLabPNChan), intent(in) :: ms
    integer, intent(in) :: a, b, c, d
    integer :: bra, ket
    real(8) :: r, ph

    r = 0.d0
    ph = ms%GetPhase(a,b) * ms%GetPhase(c,d)
    bra = ms%GetIndex(a,b)
    ket = ms%GetIndex(c,d)
    if(bra * ket == 0) return
    r = this%m(bra,ket) * ph

  end function GetScalar

  function GetTensor(this,chbra,chket,a,b,c,d) result(r)
    class(TwoBodyLabOpChan), intent(in) :: this
    type(TwoBodyLabPNChan), intent(in) :: chbra,chket
    integer, intent(in) :: a, b, c, d
    integer :: bra, ket
    real(8) :: r, ph

    r = 0.d0
    ph = chbra%GetPhase(a,b) * chket%GetPhase(c,d)
    bra = chbra%GetIndex(a,b)
    ket = chket%GetIndex(c,d)
    if(bra * ket == 0) return
    r = this%m(bra,ket) * ph
  end function GetTensor

  subroutine TMTransScalar(this,rel2lab,oprel,svd_rank_tcoef)
    use TwoBodyTransCoef
    use TwoBodyRelOps
    class(TwoBodyLabOp), intent(inout) :: this
    type(TransRel2LabSpace), intent(in) :: rel2lab
    type(TwoBodyRelOp), intent(in) :: oprel
    integer, intent(in), optional :: svd_rank_tcoef
    type(TwoBodyLabPNSpace), pointer :: lab
    type(TwoBodyRelSpaceHOBasis), pointer :: rel
    integer :: emax, e2max, Jmax, nch
    integer :: lcm, ncm, prel, jrel, ich
    integer, allocatable :: jpnl_rel(:,:)
    integer :: svd_rank_tcoef_

#ifndef MPI
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'("############################################################")')
    write(*,*)
#endif
    svd_rank_tcoef_ = -1
    if(present(svd_rank_tcoef)) svd_rank_tcoef_ = svd_rank_tcoef
    lab => this%ms
    rel => oprel%rel_sp_ket
    emax = lab%GetEmax()
    e2max = lab%GetE2max()
    Jmax = rel%GetJmax()
    nch = (Jmax + 1) * 2 * (e2max/2 + 1) * (e2max + 1)
    allocate(jpnl_rel(4,nch))
    ich = 0
    do lcm = 0, e2max
      do ncm = 0, e2max/2
        do prel = 1, -1, -2
          do jrel = 0, Jmax
            ich = ich + 1
            jpnl_rel(:,ich) = [jrel, prel, ncm, lcm]
          end do
        end do
      end do
    end do

    do ich = 1, lab%GetNumberChannels()
      call TMTransChannel(ich)
    end do
  contains
    subroutine TMTransChannel(loop)
      integer, intent(in) :: loop
      type :: trans
        real(8), allocatable :: v(:,:)
        logical :: is
      end type trans
      type(trans), allocatable :: tc(:,:,:,:)
      type(TwoBodyLabPNChan), pointer :: ChLab
      type(TransRel2LabChannel), pointer :: ChRel2Lab
      type(RelativeCMQNs), pointer :: QNRelCM
      type(TwoBodyLabQNs), pointer :: QN2BLab
      type(DSingularValueDecomposition) :: DSVD
      type(DMat) :: Tcoef_svd
      integer :: n, m, nrel
      integer :: j, p, z
      integer :: looprel, ichrel
      integer :: jrel, prel, s, ncm, lcm, lrel
      integer :: nn, kk, k, num, i, ii
      real(8), allocatable :: v2(:,:), work(:,:)
      real(8) :: fac
      real(8) :: ti
      type(sys) :: sy

      ti = omp_get_wtime()

      ChLab => lab%GetChannel(loop)
      j = ChLab%GetJ()
      p = ChLab%GetParity()
      z = ChLab%GetZ()
      ChRel2Lab => rel2lab%GetChannel( j, p, z )

      n = ChLab%GetNumberStates()
      m = ChRel2Lab%GetNumberRelCMStates()
      if(n < 1) return
      if(m < 1) return

      allocate(tc(0:Jmax, -1:1, 0:e2max/2, 0:e2max), v2(n,n))
      v2 = 0.d0
      do looprel = 1, nch
        jrel = jpnl_rel(1, looprel)
        prel = jpnl_rel(2, looprel)
        ncm  = jpnl_rel(3, looprel)
        lcm  = jpnl_rel(4, looprel)
        tc(jrel, prel, ncm, lcm)%is = .false.
        ichrel = rel%GetIndex(jrel, prel, z)

        if(ichrel == 0) cycle
        if( oprel%MatCh(ichrel,ichrel)%is_Zero ) cycle

        nrel = rel%jpz(ichrel)%GetNumberStates()
        allocate(tc(jrel, prel, ncm, lcm)%v(nrel, n))
        tc(jrel, prel, ncm, lcm)%v(:,:) = 0.d0
      end do

      !$omp parallel
      !$omp do private(k, QNRelCM, nrel, lrel, s, jrel, ncm, lcm, prel, &
      !$omp & ichrel, nn, kk, num, i, QN2BLab, ii) schedule(dynamic)
      do k = 1, ChRel2Lab%GetNumberRelCMStates()
        QNRelCM => ChRel2Lab%GetRelCM(k)
        nrel = QNRelCM%nrel
        lrel = QNRelCM%lrel
        s    = QNRelCM%s
        jrel = QNRelCM%jrel
        ncm  = QNRelCM%Ncm
        lcm  = QNRelCM%Lcm
        prel = (-1) ** (lrel)
        ichrel = rel%GetIndex(jrel, prel, z)
        if(ichrel == 0) cycle

        nn = QNRelCM%Nmax
        kk = QNRelCM%idx
        num = rel%jpz(ichrel)%GetIndex(nrel, lrel, s)
        tc(jrel, prel, ncm, lcm)%is = .true.

        do i = 1, ChRel2Lab%GetNumberLabStates()
          QN2BLab => ChRel2Lab%GetLab(i)
          if(QN2BLab%Nmax /= nn) cycle
          ii = QN2BLab%idx
          tc(jrel, prel, ncm, lcm)%v(num, i) = &
            & ChRel2Lab%N2max(NN)%mat(kk, ii)
        end do
      end do
      !$omp end do
      !$omp end parallel

      do looprel = 1, nch

        jrel = jpnl_rel(1, looprel)
        prel = jpnl_rel(2, looprel)
        ncm  = jpnl_rel(3, looprel)
        lcm  = jpnl_rel(4, looprel)
        ichrel = rel%GetIndex(jrel, prel, z)
        if(ichrel == 0) cycle
        if( .not. allocated(oprel%MatCh(ichrel,ichrel)%m)) cycle
        if( .not. tc(jrel, prel, ncm, lcm)%is) cycle
        fac = 1.d0
        m = size(oprel%MatCh(ichrel,ichrel)%m, 1)
        k = size(tc(jrel, prel, ncm, lcm)%v, 2)
        if(m < 1 .or. k < 1) cycle
        allocate(work(m, k))
        call dgemm('n', 'n', m, k, m, 1.d0, &
            &  oprel%MatCh(ichrel,ichrel)%m, m, &
            &  tc(jrel, prel, ncm, lcm)%v, m, &
            &  0.d0, work, m)
        call dgemm('t', 'n', k, k, m, fac, &
            &  tc(jrel, prel, ncm, lcm)%v, m, &
            &  work, m, 1.d0, v2, k)
        deallocate(work)
      end do

      do looprel = 1, nch
        jrel = jpnl_rel(1, looprel)
        prel = jpnl_rel(2, looprel)
        ncm  = jpnl_rel(3, looprel)
        lcm  = jpnl_rel(4, looprel)
        if(allocated(tc(jrel, prel, ncm, lcm)%v)) then
          deallocate(tc(jrel, prel, ncm, lcm)%v)
        end if
      end do
      deallocate(tc)
      this%MatCh(loop,loop)%m = v2
      call timer%Add(sy%str('TMTransChannel'), omp_get_wtime() - ti)
    end subroutine TMTransChannel
  end subroutine TMTransScalar

  subroutine TMTransTensor(this,rel2lab,oprel,svd_rank_tcoef)
    use TwoBodyTransCoef
    use TwoBodyRelOps
    class(TwoBodyLabOp), intent(inout) :: this
    type(TransRel2LabSpace), intent(in) :: rel2lab
    type(TwoBodyRelOp), intent(in) :: oprel
    integer, intent(in), optional :: svd_rank_tcoef
    type(TwoBodyLabPNSpace), pointer :: lab
    type(TwoBodyRelSpaceHOBasis), pointer :: rel
    integer :: chbra, chket
    integer :: emax, e2max, Jmax, nch, nchloop, ch1dim, ch
    integer :: lcm, ncm, pbrel, jbrel
    integer :: pkrel, jkrel
    integer :: svd_rank_tcoef_
    integer, allocatable :: jpnl_rel(:,:)
    integer, allocatable :: jpnl_loop(:,:)
    integer, allocatable :: loop2chb(:), loop2chk(:)

#ifndef MPI
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'("############################################################")')
    write(*,*)
#endif
    svd_rank_tcoef_ = -1
    if(present(svd_rank_tcoef)) svd_rank_tcoef_ = svd_rank_tcoef
    rel => oprel%rel_sp_ket
    lab => this%ms
    emax = lab%GetEMax()
    e2max = lab%GetE2max()
    Jmax = rel%GetJmax()
    nch = (Jmax + 1) * 2 * (e2max/2 + 1) * (e2max + 1)
    nchloop = (Jmax + 1) * (Jmax + 1) * 4 * (e2max/2 + 1) * (e2max + 1)
    allocate(jpnl_rel(4,nch))
    ch = 0
    do lcm = 0, e2max
      do ncm = 0, e2max/2
        do pbrel = 1, -1, -2
          do jbrel = 0, Jmax
            ch = ch + 1
            jpnl_rel(:,ch) = (/jbrel, pbrel, ncm, lcm/)
          end do
        end do
      end do
    end do

    allocate(jpnl_loop(6,nchloop))
    ch = 0
    do lcm = 0, e2max
      do ncm = 0, e2max/2
        do pbrel = 1, -1, -2
          do jbrel = 0, Jmax
            do pkrel = 1, -1, -2
              do jkrel = 0, Jmax
                ch = ch + 1
                jpnl_loop(:,ch) = (/jbrel, pbrel, jkrel, pkrel, ncm, lcm/)
              end do
            end do
          end do
        end do
      end do
    end do

    ch1dim = 0
    do chbra = 1, lab%GetNumberChannels()
      do chket = 1, lab%GetNumberChannels()
        if(this%MatCh(chbra,chket)%Zero) cycle
        ch1dim = ch1dim + 1
      end do
    end do
    allocate(loop2chb(ch1dim), loop2chk(ch1dim))

    ch1dim = 0
    do chbra = 1, lab%GetNumberChannels()
      do chket = 1, lab%GetNumberChannels()
        if(this%MatCh(chbra,chket)%Zero) cycle
        ch1dim = ch1dim + 1
        loop2chb(ch1dim) = chbra
        loop2chk(ch1dim) = chket
      end do
    end do

    do ch = 1, ch1dim
      call TMTransChannel(ch)
    end do

  contains
    subroutine TMTransChannel(loop)
      use MyLibrary, only: sjs, hat
      integer, intent(in) :: loop
      integer :: chb, chk
      type :: trans
        real(8), allocatable :: v(:,:)
        logical :: is
      end type trans
      type(trans), allocatable :: tcb(:,:,:,:), tck(:,:,:,:)
      type(TwoBodyLabPNChan), pointer :: ChLabBra, ChLabKet
      type(TransRel2LabChannel), pointer :: ChRel2LabBra, ChRel2LabKet
      type(RelativeCMQNs), pointer :: QNRelCM
      type(TwoBodyLabQNs), pointer :: QN2BLab
      integer :: nb, mb, nk, mk, nbrel, nkrel
      integer :: jb, pb, zb
      integer :: jk, pk, zk
      integer :: looprel, ichbrel, ichkrel
      integer :: nrel, lrel, jrel, prel, s, ncm, lcm
      integer :: jbrel, pbrel, jkrel, pkrel
      integer :: nn, kk, k, num, i, ii
      real(8), allocatable :: v2(:,:), work(:,:)
      real(8) :: fac
      real(8) :: ti
      type(sys) :: sy

      ti = omp_get_wtime()

      chb = loop2chb(loop)
      chk = loop2chk(loop)
      ChLabBra => lab%GetChannel(chb)
      ChLabKet => lab%GetChannel(chk)

      jb = ChLabBra%GetJ()
      pb = ChLabBra%GetParity()
      zb = ChLabBra%GetZ()

      jk = ChLabKet%GetJ()
      pk = ChLabKet%GetParity()
      zk = ChLabKet%GetZ()

      ChRel2LabBra => rel2lab%GetChannel(jb, pb, zb)
      ChRel2LabKet => rel2lab%GetChannel(jk, pk, zk)

      nb = ChLabBra%GetNumberStates()
      nk = ChLabKet%GetNumberStates()
      mb = ChRel2LabBra%GetNumberRelCMStates()
      mk = ChRel2LabBra%GetNumberRelCMStates()

      if(nb < 1) return
      if(mb < 1) return
      if(nk < 1) return
      if(mk < 1) return

      allocate(tcb(0:Jmax, -1:1, 0:e2max/2, 0:e2max))
      allocate(tck(0:Jmax, -1:1, 0:e2max/2, 0:e2max))
      allocate(v2(nb,nk))
      v2 = 0.d0

      do looprel = 1, nch
        jrel = jpnl_rel(1, looprel)
        prel = jpnl_rel(2, looprel)
        ncm  = jpnl_rel(3, looprel)
        lcm  = jpnl_rel(4, looprel)
        tcb(jrel, prel, ncm, lcm)%is = .false.
        tck(jrel, prel, ncm, lcm)%is = .false.
        ichbrel = rel%GetIndex(jrel, prel, zb)
        ichkrel = rel%GetIndex(jrel, prel, zk)

        if(ichbrel /= 0) then
          nbrel = rel%jpz(ichbrel)%GetNumberStates()
          allocate(tcb(jrel, prel, ncm, lcm)%v(nbrel, nb))
          tcb(jrel, prel, ncm, lcm)%v(:,:) = 0.d0
        end if

        if(ichkrel /= 0) then
          nkrel = rel%jpz(ichkrel)%GetNumberStates()
          allocate(tck(jrel, prel, ncm, lcm)%v(nkrel, nk))
          tck(jrel, prel, ncm, lcm)%v(:,:) = 0.d0
        end if
      end do

      !$omp parallel
      !$omp do private(k, QNRelCM, nrel, lrel, s, jrel, ncm, lcm, prel, &
      !$omp &  ichbrel, nn, kk, num, i, QN2BLab, ii) schedule(dynamic)
      do k = 1, ChRel2LabBra%GetNumberRelCMStates()
        QNRelCM => ChRel2LabBra%GetRelCM(k)
        nrel = QNRelCM%nrel
        lrel = QNRelCM%lrel
        s    = QNRelCM%s
        jrel = QNRelCM%jrel
        ncm  = QNRelCM%Ncm
        lcm  = QNRelCM%Lcm
        prel = (-1) ** (lrel)
        ichbrel = rel%GetIndex(jrel, prel, zb)
        if(ichbrel == 0) cycle

        nn = QNRelCM%Nmax
        kk = QNRelCM%idx
        num = rel%jpz(ichbrel)%GetIndex(nrel, lrel, s)
        tcb(jrel, prel, ncm, lcm)%is = .true.

        do i = 1, ChRel2LabBra%GetNumberLabStates()
          QN2BLab => ChRel2LabBra%GetLab(i)
          if(QN2BLab%Nmax /= nn) cycle
          ii = QN2BLab%idx
          tcb(jrel, prel, ncm, lcm)%v(num, i) =&
              & ChRel2LabBra%N2max(nn)%mat(kk,ii)
        end do
      end do
      !$omp end do
      !$omp end parallel

      !$omp parallel
      !$omp do private(k, QNRelCM, nrel, lrel, s, jrel, ncm, lcm, prel, &
      !$omp &  ichkrel, nn, kk, num, i, QN2BLab, ii) schedule(dynamic)
      do k = 1, ChRel2LabKet%GetNumberRelCMStates()
        QNRelCM => ChRel2LabKet%GetRelCM(k)
        nrel = QNRelCM%nrel
        lrel = QNRelCM%lrel
        s    = QNRelCM%s
        jrel = QNRelCM%jrel
        ncm  = QNRelCM%Ncm
        lcm  = QNRelCM%Lcm
        prel = (-1) ** (lrel)
        ichkrel = rel%GetIndex(jrel, prel, zk)
        if(ichkrel == 0) cycle

        nn = QNRelCM%Nmax
        kk = QNRelCM%idx
        num = rel%jpz(ichkrel)%GetIndex(nrel, lrel, s)
        tck(jrel, prel, ncm, lcm)%is = .true.

        do i = 1, ChRel2LabKet%GetNumberLabStates()
          QN2BLab => ChRel2LabKet%GetLab(i)
          if(QN2BLab%Nmax /= nn) cycle
          ii = QN2BLab%idx
          tck(jrel, prel, ncm, lcm)%v(num, i) =&
              & ChRel2LabKet%N2max(nn)%mat(kk,ii)
        end do
      end do
      !$omp end do
      !$omp end parallel

      do looprel = 1, nchloop
        jbrel= jpnl_loop(1, looprel)
        pbrel= jpnl_loop(2, looprel)
        jkrel= jpnl_loop(3, looprel)
        pkrel= jpnl_loop(4, looprel)
        ncm  = jpnl_loop(5, looprel)
        lcm  = jpnl_loop(6, looprel)
        ichbrel = rel%GetIndex(jbrel, pbrel, zb)
        ichkrel = rel%GetIndex(jkrel, pkrel, zk)
        if(ichbrel == 0) cycle
        if(ichkrel == 0) cycle
        if(oprel%MatCh(ichbrel,ichkrel)%is_Zero) cycle
        if( .not. tcb(jbrel, pbrel, ncm, lcm)%is) cycle
        if( .not. tck(jkrel, pkrel, ncm, lcm)%is) cycle

        fac = (-1.d0) ** (lcm + jkrel + jb + oprel%GetOpJ()) * hat(2*jb) * hat(2*jk) * &
            & sjs(2*jbrel, 2*jkrel, 2*oprel%GetOpJ(), 2*jk, 2*jb, 2*lcm)

        nb = size(tcb(jbrel, pbrel, ncm, lcm)%v, 2)
        mb = size(tcb(jbrel, pbrel, ncm, lcm)%v, 1)
        if(mb < 1 .or. nb < 1) cycle

        nk = size(tck(jkrel, pkrel, ncm, lcm)%v, 2)
        mk = size(tck(jkrel, pkrel, ncm, lcm)%v, 1)
        if(mk < 1 .or. nk < 1) cycle
        allocate(work(mb, nk))

        call dgemm('n', 'n', mb, nk, mk, 1.d0, &
          &  oprel%MatCh(ichbrel,ichkrel)%m, mb, &
          &  tck(jkrel, pkrel, ncm, lcm)%v, mk, &
          &  0.d0, work, mb)

        call dgemm('t', 'n', nb, nk, mb, fac, &
          &  tcb(jbrel, pbrel, ncm, lcm)%v, mb, &
          &  work, mb, 1.d0, v2, nb)
        deallocate(work)
      end do
      this%MatCh(chb, chk)%m = v2
      do looprel = 1, nch
        jbrel= jpnl_rel(1, looprel)
        pbrel= jpnl_rel(2, looprel)
        ncm  = jpnl_rel(3, looprel)
        lcm  = jpnl_rel(4, looprel)
        if(allocated(tcb(jbrel, pbrel, ncm, lcm)%v)) then
          deallocate(tcb(jbrel, pbrel, ncm, lcm)%v)
        end if

        if(allocated(tck(jbrel, pbrel, ncm, lcm)%v)) then
          deallocate(tck(jbrel, pbrel, ncm, lcm)%v)
        end if
      end do
      deallocate(tcb, tck)
      call timer%Add(sy%str('TMTransChannel'), omp_get_wtime() - ti)
    end subroutine TMTransChannel
  end subroutine TMTransTensor

  subroutine TMTransCMTensor(this,rel2lab,oprelcm)
    use TwoBodyTransCoef
    use TwoBodyRelCMOps
    class(TwoBodyLabOp), intent(inout) :: this
    type(TransRel2LabSpace), intent(in) :: rel2lab
    type(TwoBodyRelCMOp), intent(in) :: oprelcm
    type(TwoBodyLabPNSpace), pointer :: lab
    type(TwoBodyRelCMSpaceHOBasis), pointer :: rel
    integer :: chbra, chket
    integer :: emax, e2max, JRelmax, LcmMax, nch, nchloop, ch1dim, ch
    integer :: lcm_bra, lcm_ket, pbrel, jbrel
    integer :: pkrel, jkrel
    integer, allocatable :: jpnl_rel(:,:)
    integer, allocatable :: jpnl_loop(:,:)
    integer, allocatable :: loop2chb(:), loop2chk(:)

#ifndef MPI
    write(*,*)
    write(*,'("############################################################")')
    write(*,'(4x, "Doing transformaion to sp basis ...")')
    write(*,'("############################################################")')
    write(*,*)
#endif
    rel => oprelcm%ms
    lab => this%ms
    emax = lab%GetEMax()
    e2max = lab%GetE2max()
    JRelmax = rel%GetJRelmax()
    LcmMax = rel%GetLcmMax()

    nch = (JRelmax + 1) * 2 * (LcmMax + 1)
    nchloop = nch**2
    allocate(jpnl_rel(3,nch))
    ch = 0
    do lcm_bra = 0, LcmMax
      do pbrel = 1, -1, -2
        do jbrel = 0, JRelmax
          ch = ch + 1
          jpnl_rel(:,ch) = (/jbrel, pbrel, lcm_bra/)
        end do
      end do
    end do

    allocate(jpnl_loop(6,nchloop))
    ch = 0
    do lcm_bra = 0, LcmMax
      do pbrel = 1, -1, -2
        do jbrel = 0, JRelmax
          do lcm_ket = 0, LcmMax
            do pkrel = 1, -1, -2
              do jkrel = 0, JRelmax
                ch = ch + 1
                jpnl_loop(:,ch) = (/jbrel, pbrel, lcm_bra, jkrel, pkrel, lcm_ket/)
              end do
            end do
          end do
        end do
      end do
    end do

    ch1dim = 0
    do chbra = 1, lab%GetNumberChannels()
      do chket = 1, lab%GetNumberChannels()
        if(this%MatCh(chbra,chket)%Zero) cycle
        ch1dim = ch1dim + 1
      end do
    end do
    allocate(loop2chb(ch1dim), loop2chk(ch1dim))

    ch1dim = 0
    do chbra = 1, lab%GetNumberChannels()
      do chket = 1, lab%GetNumberChannels()
        if(this%MatCh(chbra,chket)%Zero) cycle
        ch1dim = ch1dim + 1
        loop2chb(ch1dim) = chbra
        loop2chk(ch1dim) = chket
      end do
    end do

    do ch = 1, ch1dim
      call TMTransChannel(ch)
    end do

    if(this%GetOpJ()==0 .and. this%GetOpP()==0 .and. this%GetOpZ()==0 .and. this%reduced_me() .and. oprelcm%reduced_me()) then
      call this%ReducedToNonReduced()
    end if

  contains
    subroutine TMTransChannel(loop)
      use MyLibrary, only: sjs, hat
      integer, intent(in) :: loop
      integer :: chb, chk
      type :: trans
        real(8), allocatable :: v(:,:)
        logical :: is
      end type trans
      type(trans), allocatable :: tcb(:,:,:), tck(:,:,:)
      type(TwoBodyLabPNChan), pointer :: ChLabBra, ChLabKet
      type(TransRel2LabChannel), pointer :: ChRel2LabBra, ChRel2LabKet
      type(TwoBodyRelCMChanHOBasis), pointer :: ChRelCMBra, ChRelCMKet
      type(RelativeCMQNs), pointer :: QNRelCM
      type(TwoBodyLabQNs), pointer :: QN2BLab
      integer :: nb, mb, nk, mk, nbrel, nkrel
      integer :: jb, pb, zb
      integer :: jk, pk, zk
      integer :: looprel, ichbrel, ichkrel
      integer :: nrel, lrel, jrel, prel, s, ncm, lcm, lcm_bra, lcm_ket
      integer :: jbrel, pbrel, jkrel, pkrel
      integer :: nn, kk, k, num, i, ii
      real(8), allocatable :: v2(:,:), work(:,:)
      real(8) :: fac
      real(8) :: ti
      type(sys) :: sy

      ti = omp_get_wtime()

      chb = loop2chb(loop)
      chk = loop2chk(loop)
      ChLabBra => lab%GetChannel(chb)
      ChLabKet => lab%GetChannel(chk)

      jb = ChLabBra%GetJ()
      pb = ChLabBra%GetParity()
      zb = ChLabBra%GetZ()


      jk = ChLabKet%GetJ()
      pk = ChLabKet%GetParity()
      zk = ChLabKet%GetZ()

      ChRel2LabBra => rel2lab%GetChannel(jb, pb, zb)
      ChRel2LabKet => rel2lab%GetChannel(jk, pk, zk)

      nb = ChLabBra%GetNumberStates()
      nk = ChLabKet%GetNumberStates()
      mb = ChRel2LabBra%GetNumberRelCMStates()
      mk = ChRel2LabBra%GetNumberRelCMStates()

      if(jb > rel%GetJmax()) return
      if(jk > rel%GetJmax()) return
      if(nb < 1) return
      if(mb < 1) return
      if(nk < 1) return
      if(mk < 1) return

      allocate(tcb(0:JRelMax, -1:1, 0:LcmMax))
      allocate(tck(0:JRelMax, -1:1, 0:LcmMax))
      allocate(v2(nb,nk))
      do looprel = 1, nch
        jrel = jpnl_rel(1, looprel)
        prel = jpnl_rel(2, looprel)
        lcm  = jpnl_rel(3, looprel)
        tcb(jrel, prel, lcm)%is = .false.
        tck(jrel, prel, lcm)%is = .false.
        ichbrel = rel%GetIndex(jb, pb, zb, jrel, lcm)
        ichkrel = rel%GetIndex(jk, pk, zk, jrel, lcm)
        ChRelCMBra => rel%GetChannel(ichbrel)
        ChRelCMKet => rel%GetChannel(ichkrel)

        if(ichbrel /= 0) then
          nbrel = ChRelCMBra%GetNumberStates()
          allocate(tcb(jrel, prel, lcm)%v(nbrel, nb))
          tcb(jrel, prel, lcm)%v(:,:) = 0.d0
        end if

        if(ichkrel /= 0) then
          nkrel = ChRelCMKet%GetNumberStates()
          allocate(tck(jrel, prel, lcm)%v(nkrel, nk))
          tck(jrel, prel, lcm)%v(:,:) = 0.d0
        end if
      end do

      !$omp parallel
      !$omp do private(k, QNRelCM, nrel, lrel, s, jrel, ncm, lcm, prel, &
      !$omp &  ichbrel, ChRelCMBra, nn, kk, num, i, QN2BLab, ii) schedule(dynamic)
      do k = 1, ChRel2LabBra%GetNumberRelCMStates()
        QNRelCM => ChRel2LabBra%GetRelCM(k)
        nrel = QNRelCM%nrel
        lrel = QNRelCM%lrel
        s    = QNRelCM%s
        jrel = QNRelCM%jrel
        ncm  = QNRelCM%Ncm
        lcm  = QNRelCM%Lcm
        prel = (-1)**lrel
        if(lcm > rel%GetLcmMax()) cycle
        if(2*(nrel+ncm)+lrel+lcm > rel%GetNmax()) cycle
        ichbrel = rel%GetIndex(jb, pb, zb, jrel, lcm)
        if(ichbrel == 0) cycle
        ChRelCMBra => rel%GetChannel(ichbrel)

        nn = QNRelCM%Nmax
        kk = QNRelCM%idx
        num = ChRelCMBra%GetIndex(nrel, lrel, s, ncm)
        tcb(jrel, prel, lcm)%is = .true.

        do i = 1, ChRel2LabBra%GetNumberLabStates()
          QN2BLab => ChRel2LabBra%GetLab(i)
          if(QN2BLab%Nmax /= nn) cycle
          ii = QN2BLab%idx
          tcb(jrel, prel, lcm)%v(num, i) = ChRel2LabBra%N2max(nn)%mat(kk,ii)
        end do
      end do
      !$omp end do
      !$omp end parallel

      !$omp parallel
      !$omp do private(k, QNRelCM, nrel, lrel, s, jrel, ncm, lcm, prel, &
      !$omp &  ichkrel, ChRelCMKet, nn, kk, num, i, QN2BLab, ii) schedule(dynamic)
      do k = 1, ChRel2LabKet%GetNumberRelCMStates()
        QNRelCM => ChRel2LabKet%GetRelCM(k)
        nrel = QNRelCM%nrel
        lrel = QNRelCM%lrel
        s    = QNRelCM%s
        jrel = QNRelCM%jrel
        ncm  = QNRelCM%Ncm
        lcm  = QNRelCM%Lcm
        prel = (-1)**lrel
        if(lcm > rel%GetLcmMax()) cycle
        if(2*(nrel+ncm)+lrel+lcm > rel%GetNmax()) cycle
        ichkrel = rel%GetIndex(jk, pk, zk, jrel, lcm)
        if(ichkrel == 0) cycle
        ChRelCMKet => rel%GetChannel(ichkrel)

        nn = QNRelCM%Nmax
        kk = QNRelCM%idx
        num = ChRelCMKet%GetIndex(nrel, lrel, s, ncm)
        tck(jrel, prel, lcm)%is = .true.

        do i = 1, ChRel2LabKet%GetNumberLabStates()
          QN2BLab => ChRel2LabKet%GetLab(i)
          if(QN2BLab%Nmax /= nn) cycle
          ii = QN2BLab%idx
          tck(jrel, prel, lcm)%v(num, i) = ChRel2LabKet%N2max(nn)%mat(kk,ii)
        end do
      end do
      !$omp end do
      !$omp end parallel

      v2 = 0.d0
      do looprel = 1, nchloop
        jbrel=   jpnl_loop(1, looprel)
        pbrel=   jpnl_loop(2, looprel)
        lcm_bra= jpnl_loop(3, looprel)
        jkrel=   jpnl_loop(4, looprel)
        pkrel=   jpnl_loop(5, looprel)
        lcm_ket= jpnl_loop(6, looprel)
        if(pbrel * (-1)**lcm_bra /= pb) cycle
        if(pkrel * (-1)**lcm_ket /= pk) cycle
        ichbrel = rel%GetIndex(jb, pb, zb, jbrel, lcm_bra)
        ichkrel = rel%GetIndex(jk, pk, zk, jkrel, lcm_ket)
        if(ichbrel == 0) cycle
        if(ichkrel == 0) cycle
        if(oprelcm%MatCh(ichbrel,ichkrel)%is_Zero) cycle
        if( .not. tcb(jbrel, pbrel, lcm_bra)%is) cycle
        if( .not. tck(jkrel, pkrel, lcm_ket)%is) cycle

        fac = 1.d0
        nb = size(tcb(jbrel, pbrel, lcm_bra)%v, 2)
        mb = size(tcb(jbrel, pbrel, lcm_bra)%v, 1)
        if(mb < 1 .or. nb < 1) cycle

        nk = size(tck(jkrel, pkrel, lcm_ket)%v, 2)
        mk = size(tck(jkrel, pkrel, lcm_ket)%v, 1)
        if(mk < 1 .or. nk < 1) cycle
        allocate(work(mb, nk))

        call dgemm('n', 'n', mb, nk, mk, 1.d0, &
            & oprelcm%MatCh(ichbrel,ichkrel)%m, mb, &
            & tck(jkrel, pkrel, lcm_ket)%v, mk, &
            & 0.d0, work, mb)
        call dgemm('t', 'n', nb, nk, mb, fac, &
            & tcb(jbrel, pbrel, lcm_bra)%v, mb, &
            & work, mb, 1.d0, v2, nb)
        deallocate(work)
      end do
      this%MatCh(chb, chk)%m = v2
      do looprel = 1, nch
        jrel= jpnl_rel(1, looprel)
        prel= jpnl_rel(2, looprel)
        lcm = jpnl_rel(3, looprel)
        if(allocated(tcb(jrel, prel, lcm)%v)) then
          deallocate(tcb(jrel, prel, lcm)%v)
        end if

        if(allocated(tck(jrel, prel, lcm)%v)) then
          deallocate(tck(jrel, prel, lcm)%v)
        end if
      end do
      deallocate(tcb, tck)
      call timer%Add(sy%str('TMTransChannel'), omp_get_wtime() - ti)
    end subroutine TMTransChannel
  end subroutine TMTransCMTensor

  subroutine TMTransCMTensor_explicit(this,rel2lab,oprelcm)
    use TwoBodyTransCoef
    use TwoBodyRelCMOps
    class(TwoBodyLabOp), intent(inout) :: this
    type(TransRel2LabSpace), intent(in) :: rel2lab
    type(TwoBodyRelCMOp), intent(in) :: oprelcm
    type(TwoBodyLabPNSpace), pointer :: lab
    type(TwoBodyLabPNChan), pointer :: chbra, chket
    type(Orbits), pointer :: sps
    integer :: ichbra, ichket, ibra, iket
    integer :: jbra, pbra, zbra, jket, pket, zket, ebra, eket
    integer :: ncmbra, lcmbra, nrelbra, lrelbra, sbra, jrelbra
    integer :: ncmket, lcmket, nrelket, lrelket, sket, jrelket
    integer :: bras(8), kets(8)
    integer :: emax, e2max, JRelmax, LcmMax
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: me

    lab => this%ms
    sps => lab%sps
    emax = lab%GetEMax()
    e2max = lab%GetE2max()
    JRelmax = oprelcm%ms%GetJRelmax()
    LcmMax = oprelcm%ms%GetLcmMax()

    do ichbra = 1, lab%GetNumberChannels()
      do ichket = 1, lab%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%Zero) cycle
        chbra => lab%GetChannel(ichbra)
        chket => lab%GetChannel(ichket)
        jbra = chbra%GetJ()
        pbra = chbra%GetParity()
        zbra = chbra%GetZ()
        jket = chket%GetJ()
        pket = chket%GetParity()
        zket = chket%GetZ()

        do ibra = 1, chbra%GetNumberStates()
          do iket = 1, chket%GetNumberStates()
            a = chbra%n2label1( ibra )
            b = chbra%n2label2( ibra )
            c = chket%n2label1( iket )
            d = chket%n2label2( iket )
            oa => sps%GetOrbit(a)
            ob => sps%GetOrbit(b)
            oc => sps%GetOrbit(c)
            od => sps%GetOrbit(d)
            ebra = 2*oa%n+oa%l+2*ob%n+ob%l
            eket = 2*oc%n+oc%l+2*od%n+od%l

            me = 0.d0
            do lrelbra = 0, e2max
              do lcmbra = 0, e2max
                if(lcmbra > LcmMax) cycle
                do nrelbra = 0, e2max/2
                  do ncmbra = 0, e2max/2
                    if(2*(nrelbra+ncmbra)+lrelbra+lcmbra /= ebra) cycle

                    do lrelket = 0, e2max
                      do lcmket = 0, e2max
                        if(lcmbra > LcmMax) cycle
                        do nrelket = 0, e2max/2
                          do ncmket = 0, e2max/2
                            if(2*(nrelket+ncmket)+lrelket+lcmket /= eket) cycle

                            do sbra = 0, 1
                              if(abs(Zbra)==1 .and. (-1)**(lrelbra+sbra)==-1) cycle
                              do sket = 0, 1
                                if(abs(Zket)==1 .and. (-1)**(lrelket+sket)==-1) cycle
                                do jrelbra = max(abs(Jbra-lcmbra), abs(sbra-lrelbra)), min(Jbra+lcmbra, sbra+lrelbra)
                                  do jrelket = max(abs(Jket-lcmket), abs(sket-lrelket)), min(Jket+lcmket, sket+lrelket)
                                    if(jrelbra > JrelMax) cycle
                                    if(jrelket > JrelMax) cycle

                                    bras = [ncmbra, lcmbra, nrelbra, lrelbra, sbra, jrelbra, Jbra, Zbra]
                                    kets = [ncmket, lcmket, nrelket, lrelket, sket, jrelket, Jket, Zket]
                                    me = me + oprelcm%Get2ME(bras, kets) * &
                                        & rel2lab%GetME(ncmbra, lcmbra, nrelbra, lrelbra, sbra, jrelbra, Jbra, oa, ob) * &
                                        & rel2lab%GetME(ncmket, lcmket, nrelket, lrelket, sket, jrelket, Jket, oc, od)
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
            this%MatCh(ichbra,ichket)%m(ibra,iket) = me
          end do
        end do

      end do
    end do
  end subroutine TMTransCMTensor_explicit

  function GetFileNameTwoBodyLabOp(this, filename, &
        &  pot, renorm, lambda, hw, emax, e2max, coul, &
        &  snt_mass) result(f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: filename, pot, renorm
    real(8), intent(in) :: lambda, hw
    integer, intent(in) :: emax, e2max, snt_mass
    logical, intent(in) :: coul
    type(sys) :: s
    type(str) :: f, OpName

    OpName = this%GetOpName()
    if(filename%val /= "default") then
      f = this%GetOpName() + "_" + filename
      if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') f = filename
      return
    end if

    f = OpName + '_TwBME'
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') f = "TwBME"
    f = f + '-HO_NN-only_' + pot + '_' + renorm
    if(renorm%val == 'srg' .or. renorm%val == 'vlowk') f = f + s%str(lambda)
    f = f + '_hw' + s%str(hw) + '_emax' + s%str(emax)
    f = f + '_e2max' + s%str(e2max)
    if(.not. coul) f = f + '_wocoul'
    f = f + '_A' + s%str(snt_mass) + '.snt'
  end function GetFileNameTwoBodyLabOp

  subroutine set_two_body_hcm( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: tbc
    type(Orbits), pointer :: sps
    type(OperatorDef) :: opdef
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ch, bra, ket, a, b, c, d, J
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: me
    type(sys) :: s

    ms => this%ms
    sps => ms%sps
    call opdef%InitOpDef(s%str("HOHamil"), .true.)

    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      do b =1, sps%norbs
        ob => sps%GetOrbit(b)

        this%OneBody%m(a,b) = CalcMEOneBody( opdef, &
            & [oa%n,oa%l,oa%j,oa%z], [ob%n,ob%l,ob%j,ob%z] )
      end do
    end do

    do ch = 1, ms%GetNumberChannels()
      tbc => ms%GetChannel(ch)
      J = tbc%GetJ()
      do bra = 1, tbc%GetNumberStates()
        a = tbc%n2label1(bra)
        b = tbc%n2label2(bra)
        oa => sps%GetOrbit(a)
        ob => sps%GetOrbit(b)
        aa = [ oa%n, oa%l, oa%j, oa%z ]
        bb = [ ob%n, ob%l, ob%j, ob%z ]
        do ket = 1, bra
          c = tbc%n2label1(ket)
          d = tbc%n2label2(ket)
          oc => sps%GetOrbit(c)
          od => sps%GetOrbit(d)
          cc = [ oc%n, oc%l, oc%j, oc%z ]
          dd = [ od%n, od%l, od%j, od%z ]
          me = ( p_dot_p(aa,bb,cc,dd,J) + r_dot_r(aa,bb,cc,dd,J) ) * ms%GetFrequency()
          this%MatCh(ch,ch)%m(bra,ket) = me
          this%MatCh(ch,ch)%m(ket,bra) = me
        end do
      end do
    end do
  end subroutine set_two_body_hcm

  subroutine set_two_body_hcm2( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: tbc
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    type(TwoBodyLabOp) :: op
    integer :: ch, bra, ket, a, b, c, d, J
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: me
    type(sys) :: s

    ms => this%ms
    sps => ms%sps
    call op%init(ms,s%str("HOHamil"))
    do ch = 1, ms%GetNumberChannels()
      tbc => ms%GetChannel(ch)
      J = tbc%GetJ()
      do bra = 1, tbc%GetNumberStates()
        a = tbc%n2label1(bra)
        b = tbc%n2label2(bra)
        oa => sps%GetOrbit(a)
        ob => sps%GetOrbit(b)
        aa = [ oa%n, oa%l, oa%j, oa%z ]
        bb = [ ob%n, ob%l, ob%j, ob%z ]
        do ket = 1, bra
          c = tbc%n2label1(ket)
          d = tbc%n2label2(ket)
          oc => sps%GetOrbit(c)
          od => sps%GetOrbit(d)
          cc = [ oc%n, oc%l, oc%j, oc%z ]
          dd = [ od%n, od%l, od%j, od%z ]
          me = ( p_dot_p(aa,bb,cc,dd,J) + r_dot_r(aa,bb,cc,dd,J) ) * ms%GetFrequency() / 2.d0
          this%MatCh(ch,ch)%m(bra,ket) = me
          this%MatCh(ch,ch)%m(ket,bra) = me
        end do
      end do
    end do
    call SubtractOneBody(op)
    do ch = 1, ms%GetNumberChannels()
      this%MatCh(ch,ch)%DMat = this%MatCh(ch,ch)%DMat - (1.d0/2.d0)*op%MatCh(ch,ch)%DMat
    end do
  end subroutine set_two_body_hcm2

  subroutine set_two_body_hrel2( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: tbc
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    type(TwoBodyLabOp) :: op
    integer :: ch, bra, ket, a, b, c, d, J
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: me
    type(sys) :: s

    ms => this%ms
    sps => ms%sps
    call op%init(ms,s%str("HOHamil"))
    do ch = 1, ms%GetNumberChannels()
      tbc => ms%GetChannel(ch)
      J = tbc%GetJ()
      do bra = 1, tbc%GetNumberStates()
        a = tbc%n2label1(bra)
        b = tbc%n2label2(bra)
        oa => sps%GetOrbit(a)
        ob => sps%GetOrbit(b)
        aa = [ oa%n, oa%l, oa%j, oa%z ]
        bb = [ ob%n, ob%l, ob%j, ob%z ]
        do ket = 1, bra
          c = tbc%n2label1(ket)
          d = tbc%n2label2(ket)
          oc => sps%GetOrbit(c)
          od => sps%GetOrbit(d)
          cc = [ oc%n, oc%l, oc%j, oc%z ]
          dd = [ od%n, od%l, od%j, od%z ]
          me = -( p_dot_p(aa,bb,cc,dd,J) + r_dot_r(aa,bb,cc,dd,J) ) * ms%GetFrequency() / 2.d0
          this%MatCh(ch,ch)%m(bra,ket) = me
          this%MatCh(ch,ch)%m(ket,bra) = me
        end do
      end do
    end do
    call SubtractOneBody(op)
    do ch = 1, ms%GetNumberChannels()
      this%MatCh(ch,ch)%DMat = this%MatCh(ch,ch)%DMat - (1.d0/2.d0)*op%MatCh(ch,ch)%DMat
    end do
  end subroutine set_two_body_hrel2

  subroutine set_two_body_kin2( this )
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: tbc
    type(Orbits), pointer :: sps
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    type(TwoBodyLabOp) :: op
    integer :: ch, bra, ket, a, b, c, d, J
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: me
    type(sys) :: s

    ms => this%ms
    sps => ms%sps
    call op%init(ms,s%str("Kinetic"))
    do ch = 1, ms%GetNumberChannels()
      tbc => ms%GetChannel(ch)
      J = tbc%GetJ()
      do bra = 1, tbc%GetNumberStates()
        a = tbc%n2label1(bra)
        b = tbc%n2label2(bra)
        oa => sps%GetOrbit(a)
        ob => sps%GetOrbit(b)
        aa = [ oa%n, oa%l, oa%j, oa%z ]
        bb = [ ob%n, ob%l, ob%j, ob%z ]
        do ket = 1, bra
          c = tbc%n2label1(ket)
          d = tbc%n2label2(ket)
          oc => sps%GetOrbit(c)
          od => sps%GetOrbit(d)
          cc = [ oc%n, oc%l, oc%j, oc%z ]
          dd = [ od%n, od%l, od%j, od%z ]
          me = -p_dot_p(aa,bb,cc,dd,J) * ms%GetFrequency() / 2.d0
          this%MatCh(ch,ch)%m(bra,ket) = me
          this%MatCh(ch,ch)%m(ket,bra) = me
        end do
      end do
    end do
    call SubtractOneBody(op)
    do ch = 1, ms%GetNumberChannels()
      this%MatCh(ch,ch)%DMat = this%MatCh(ch,ch)%DMat - (1.d0-1.d0/2.d0)*op%MatCh(ch,ch)%DMat
    end do
  end subroutine set_two_body_kin2

  subroutine WriteOperator(this, f, fn_averaged)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(str), intent(in), optional :: fn_averaged
    type(sys) :: s
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Error:", __LINE__, __FILE__
      return
    end if

    if( this%GetSubtract() ) then
      write(*,*) "One-Body term is subtracted"
      call this%OneBody%zeros( this%ms%sps%norbs, this%ms%sps%norbs )
      call SubtractOneBody(this)
    end if
    ti = omp_get_wtime()
    if( this%GetOpJ()==0 .and. this%GetOpP()==1 .and. this%GetOpZ()==0) call write_scalar_operator(this, f)
    if( this%GetOpJ()/=0 .or.  this%GetOpP()/=1 .or.  this%GetOpZ()/=0) call write_general_operator(this, f)
    if( present(fn_averaged) ) call write_averaged_operator(this, fn_averaged)

    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteOperator

  subroutine write_averaged_operator( this, fn )
    type(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: fn
    type(TwoBodyLabPNSpace), pointer :: ms
    type(TwoBodyLabPNChan), pointer :: chbra, chket
    type(Orbits), pointer :: sps
    integer :: wunit = 22
    integer :: ichbra, ichket, iket, ibra, ebra, eket, cnt
    integer :: a, b, c, d
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: vsum, ss
    logical :: is_norm, is_averaged_norm, is_mean
    type(sys) :: s

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
    ms => this%ms
    sps => ms%sps
    open(wunit, file=fn%val, action="write")
    write(wunit,'(a)') "Jbra Pbra Tbra Ebra Jket Pket Tket Eket Mean Std dev."
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%zero ) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        do ebra = 0, ms%GetE2max()
          do eket = 0, ms%GetE2max()

            cnt = 0
            vsum = 0.d0
            !$omp parallel
            !$omp do private(ibra, a, b, oa, ob, iket, c, d, oc, od) reduction(+:cnt,vsum)
            do ibra = 1, chbra%GetNumberStates()
              a = chbra%n2label1( ibra )
              b = chbra%n2label2( ibra )
              oa => sps%GetOrbit(a)
              ob => sps%GetOrbit(b)
              if( oa%e + ob%e /= ebra ) cycle
              do iket = 1, chket%GetNumberStates()
                c = chket%n2label1( iket )
                d = chket%n2label2( iket )
                oc => sps%GetOrbit(c)
                od => sps%GetOrbit(d)
                if( oc%e + od%e /= eket ) cycle
                cnt = cnt + 1
                if(is_norm .or. is_averaged_norm) vsum = vsum + this%MatCh(ichbra,ichket)%m(ibra,iket)**2
                if(is_mean) vsum = vsum + this%MatCh(ichbra,ichket)%m(ibra,iket)
              end do
            end do
            !$omp end do
            !$omp end parallel
            if(cnt == 0) cycle
            if( is_mean .or. is_averaged_norm ) then
              vsum = vsum / dble(cnt)
              ss = 0.d0
              !$omp parallel
              !$omp do private(ibra, a, b, oa, ob, iket, c, d, oc, od) reduction(+:ss)
              do ibra = 1, chbra%GetNumberStates()
                a = chbra%n2label1( ibra )
                b = chbra%n2label2( ibra )
                oa => sps%GetOrbit(a)
                ob => sps%GetOrbit(b)
                if( oa%e + ob%e /= ebra ) cycle
                do iket = 1, chket%GetNumberStates()
                  c = chket%n2label1( iket )
                  d = chket%n2label2( iket )
                  oc => sps%GetOrbit(c)
                  od => sps%GetOrbit(d)
                  if( oc%e + od%e /= eket ) cycle
                  if(is_averaged_norm) ss = ss + ( vsum - this%MatCh(ichbra,ichket)%m(ibra,iket)**2 )**2
                  if(is_mean) ss = ss + ( vsum - this%MatCh(ichbra,ichket)%m(ibra,iket) )**2
                end do
              end do
              !$omp end do
              !$omp end parallel
              ss = sqrt( ss / dble(cnt) )
              write(wunit,'(8i3,2f16.8)') chbra%GetJ(), chbra%GetParity(), chbra%GetZ(), ebra, &
                  & chket%GetJ(), chket%GetParity(), chket%GetZ(), eket, vsum, ss
            else
              write(wunit,'(8i3,f16.8)') chbra%GetJ(), chbra%GetParity(), chbra%GetZ(), ebra, &
                  & chket%GetJ(), chket%GetParity(), chket%GetZ(), eket, vsum
            end if

          end do
        end do

      end do
    end do
    close(wunit)
  end subroutine write_averaged_operator

  subroutine write_scalar_operator(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(str) :: OpName
    type(sys) :: s

    if(s%find(f,s%str(".myg.bin"))) then
      call write_scalar_operator_binary_myg(this, f)
      return
    end if

    if(s%find(f,s%str(".myg"))) then
      call write_scalar_operator_ascii_myg(this, f)
      return
    end if

    if(s%find(f,s%str(".Akshell.snt"))) then
      call write_scalar_operator_kshell_A_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".kshell.snt"))) then
      call write_scalar_operator_kshell_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".op.me2j.gz"))) then
      call write_general_operator_gzip(this, f)
      return
    end if

    if(s%find(f,s%str(".op.snt"))) then
      call write_general_operator_kshell_ascii(this, f)
      return
    end if

    if(s%find(f,s%str(".snt"))) then
      call write_scalar_operator_ascii_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j.bin"))) then
      call write_scalar_operator_binary_me2j(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j.gz"))) then
      OpName = this%GetOpName()
      select case(OpName%val)
      case('hamil', 'NNint')
        call write_scalar_operator_gzip_me2j(this, f)
      case default
        call write_general_operator_gzip(this, f)
      end select
      return
    end if

    if(s%find(f,s%str(".me2j"))) then
      call write_scalar_operator_ascii_me2j(this, f)
      return
    end if

    if(s%find(f,s%str(".MFDn"))) then
      call write_scalar_operator_ascii_mfd(this, f)
      return
    end if

  end subroutine write_scalar_operator

  subroutine write_general_operator(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str(".op.me2j.gz"))) then
      call write_general_operator_gzip(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j.gz"))) then
      call write_general_operator_gzip(this, f)
      return
    end if

    if(s%find(f,s%str(".op.snt"))) then
      call write_general_operator_kshell_ascii(this, f)
      return
    end if

    if(s%find(f,s%str(".snt"))) then
      call write_general_operator_kshell_ascii(this, f)
      return
    end if


    if(s%find(f,s%str(".op.me2j"))) then
      call write_general_operator_ascii(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j"))) then
      call write_general_operator_ascii(this, f)
      return
    end if

  end subroutine write_general_operator

  subroutine write_scalar_operator_ascii_myg(this, f)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, bra, ket
    integer :: j, ich ,nn
    integer :: cnt
    character(255) :: header
    type(str) :: OpName
    ms => this%ms
    cnt = 0
    do ich = 1, ms%GetNumberChannels()
      cnt = cnt + ms%jpz(ich)%GetNumberStates() * (ms%jpz(ich)%GetNumberStates() + 1) / 2
    end do
    open(15, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15,*) cnt
    write(15, '(4a)') "#### a, b, c, d, J,  <ab| scalar |cd> "
    do ich = 1, ms%GetNumberChannels()
      j    = ms%jpz(ich)%GetJ()
      nn   = ms%jpz(ich)%GetNumberStates()
      if(nn < 1) cycle
      do bra = 1, nn
        a = ms%jpz(ich)%n2label1(bra)
        b = ms%jpz(ich)%n2label2(bra)
        do ket = 1, bra
          c = ms%jpz(ich)%n2label1(ket)
          d = ms%jpz(ich)%n2label2(ket)
          write(15, '(5i4, f16.8)') a, b, c, d, &
              & J, this%MatCh(ich,ich)%m(bra, ket)
        end do
      end do
    end do
    close(15)
  end subroutine write_scalar_operator_ascii_myg

  subroutine write_scalar_operator_binary_myg(this, f)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, bra, ket
    integer :: j, ipar, itz, ich, nn
    integer, allocatable :: ia(:), ib(:), ic(:), id(:), jj(:)
    real(8), allocatable :: v2(:)
    integer :: cnt
    character(255) :: header
    type(str) :: OpName
    ms => this%ms
    cnt = 0
    do ich = 1, ms%GetNumberChannels()
      cnt = cnt + ms%jpz(ich)%GetNumberStates() * (ms%jpz(ich)%GetNumberStates() + 1) / 2
    end do
    allocate(ia(cnt))
    allocate(ib(cnt))
    allocate(ic(cnt))
    allocate(id(cnt))
    allocate(jj(cnt), v2(cnt))
    cnt = 0
    do ich = 1, ms%GetNumberChannels()
      j    = ms%jpz(ich)%GetJ()
      ipar = ms%jpz(ich)%GetParity()
      itz  = ms%jpz(ich)%GetZ()
      nn   = ms%jpz(ich)%GetNumberStates()
      if(nn < 1) cycle
      do bra = 1, nn
        a = ms%jpz(ich)%n2label1(bra)
        b = ms%jpz(ich)%n2label2(bra)
        do ket = 1, bra
          c = ms%jpz(ich)%n2label1(ket)
          d = ms%jpz(ich)%n2label2(ket)
          cnt = cnt + 1
          ia(cnt)   = a
          ib(cnt)   = b
          ic(cnt)   = c
          id(cnt)   = d
          jj(cnt)   = j
          v2(cnt)   = this%MatCh(ich,ich)%m(bra,ket)
        end do
      end do
    end do
    open(15, file = f%val, status = 'replace', form = 'unformatted', access = 'stream')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15) cnt
    write(15) ia, ib, ic, id, jj, v2
    close(15)
    deallocate(ia)
    deallocate(ib)
    deallocate(ic)
    deallocate(id)
    deallocate(jj, v2)
  end subroutine write_scalar_operator_binary_myg

  subroutine write_scalar_operator_ascii_snt(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    integer :: emax
    integer :: a, b, c, d, bra, ket
    integer :: j, n, cnt, ich
    character(255) :: header
    type(str) :: OpName
    two => this%ms
    sps => this%ms%sps
    emax = sps%emax

    open(15, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15, '(4i5)') sps%norbs/2, sps%norbs/2, 0, 0
    do bra = 1, sps%norbs
      write(15, '(5i5)') bra, sps%orb(bra)%n, sps%orb(bra)%l, sps%orb(bra)%j, sps%orb(bra)%z
    end do
    write(15, '(a)') '####### one-body term'
    write(15, '(i5, i3)') 0, 0
    write(15, '(a)') '### a, b, <a|t|b>'

    cnt = 0
    do ich = 1, two%GetNumberChannels()
      cnt = cnt + two%jpz(ich)%GetNumberStates() * (two%jpz(ich)%GetNumberStates() + 1) / 2
    end do
    write(15, '(i15,i3)') cnt, 0
    write(15, '(a)') '##### two-body term #####'
    write(15, '(a)') '## a, b, c, d, J, <ab|v|cd>'
    do ich = 1, two%GetNumberChannels()
      j = two%jpz(ich)%GetJ()
      n = two%jpz(ich)%GetNumberStates()
      if(n < 1) cycle
      do bra = 1, n
        a = two%jpz(ich)%n2label1(bra)
        b = two%jpz(ich)%n2label2(bra)
        do ket = 1, bra
          c = two%jpz(ich)%n2label1(ket)
          d = two%jpz(ich)%n2label2(ket)
          write(15, '(5i5, f16.8)') a, b, c, d, j, this%MatCh(ich,ich)%m(bra, ket)
        end do
      end do
    end do
    close(15)
  end subroutine write_scalar_operator_ascii_snt

  subroutine write_scalar_operator_kshell_snt(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    type(Orbits) :: sps_k
    integer :: emax
    integer :: a, b, c, d, bra, ket
    integer :: a_k, b_k, c_k, d_k
    integer :: j, n, cnt, ich, l, n1, n2
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: kine, Tcm2
    character(255) :: header
    type(str) :: OpName
    two => this%ms
    sps => this%ms%sps
    emax = sps%emax

    call sps_k%init(emax, mode="kshell")

    open(15, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15, '(4i5)') sps_k%norbs/2, sps_k%norbs/2, 0, 0
    do bra = 1, sps%norbs
      write(15, '(5i5)') bra, sps_k%orb(bra)%n, sps_k%orb(bra)%l, sps_k%orb(bra)%j, sps_k%orb(bra)%z
    end do

    cnt = 0
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(abs(sps_k%orb(bra)%n - sps_k%orb(ket)%n) > 1) cycle
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        if(sps_k%orb(bra)%z /= sps_k%orb(ket)%z) cycle
        cnt = cnt + 1
      end do
    end do

    write(15, '(a)') '####### one-body term'
    write(15, '(i5, i3, f6.2)') cnt, 10, two%GetFrequency()
    write(15, '(a)') '### a, b, <a|t|b>'
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(abs(sps_k%orb(bra)%n - sps_k%orb(ket)%n) > 1) cycle
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        if(sps_k%orb(bra)%z /= sps_k%orb(ket)%z) cycle
        l = sps_k%orb(bra)%l
        n1= sps_k%orb(bra)%n
        n2= sps_k%orb(ket)%n
        kine = 0.d0
        if(n1 == n2) kine = dble(2 * n1 + l) + 1.5d0
        if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l) + 1.5d0))
        if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l) + 1.5d0))
        write(15,'(2i5,f16.8)') bra,ket, 0.5d0 * kine
      end do
    end do

    cnt = 0
    do ich = 1, two%GetNumberChannels()
      cnt = cnt + two%jpz(ich)%GetNumberStates() * (two%jpz(ich)%GetNumberStates() + 1) / 2
    end do
    write(15, '(a)') '##### two-body term #####'
    write(15, '(i15, i3, f6.2)') cnt, 10, two%GetFrequency()
    write(15, '(a)') '## a, b, c, d, J, <ab|v|cd>'
    do ich = 1, two%GetNumberChannels()
      j = two%jpz(ich)%GetJ()
      n = two%jpz(ich)%GetNumberStates()
      if(n < 1) cycle
      do bra = 1, n
        a = two%jpz(ich)%n2label1(bra)
        b = two%jpz(ich)%n2label2(bra)
        a_k = sps_k%nljz2idx(sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j, sps%orb(a)%z)
        b_k = sps_k%nljz2idx(sps%orb(b)%n, sps%orb(b)%l, sps%orb(b)%j, sps%orb(b)%z)
        do ket = 1, bra
          c = two%jpz(ich)%n2label1(ket)
          d = two%jpz(ich)%n2label2(ket)
          c_k = sps_k%nljz2idx(sps%orb(c)%n, sps%orb(c)%l, sps%orb(c)%j, sps%orb(c)%z)
          d_k = sps_k%nljz2idx(sps%orb(d)%n, sps%orb(d)%l, sps%orb(d)%j, sps%orb(d)%z)
          aa = [sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j, sps%orb(a)%z]
          bb = [sps%orb(b)%n, sps%orb(b)%l, sps%orb(b)%j, sps%orb(b)%z]
          cc = [sps%orb(c)%n, sps%orb(c)%l, sps%orb(c)%j, sps%orb(c)%z]
          dd = [sps%orb(d)%n, sps%orb(d)%l, sps%orb(d)%j, sps%orb(d)%z]
          Tcm2 = p_dot_p(aa,bb,cc,dd,J)
          write(15, '(5i5, 2f16.8)') a_k, b_k, c_k, d_k, j, this%MatCh(ich,ich)%m(bra, ket), -Tcm2
        end do
      end do
    end do
    close(15)
    call sps_k%fin()
  end subroutine write_scalar_operator_kshell_snt

  subroutine write_scalar_operator_kshell_A_snt(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    type(Orbits) :: sps_k
    integer :: emax
    integer :: a, b, c, d, bra, ket
    integer :: a_k, b_k, c_k, d_k
    integer :: j, n, cnt, ich, l, n1, n2
    integer :: aa(4), bb(4), cc(4), dd(4)
    real(8) :: kine, Tcm2
    character(255) :: header
    type(str) :: OpName
    two => this%ms
    sps => this%ms%sps
    emax = sps%emax

    call sps_k%init(emax, mode="kshell")

    open(15, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15, '(4i5)') sps_k%norbs/2, sps_k%norbs/2, 0, 0
    do bra = 1, sps%norbs
      write(15, '(5i5)') bra, sps_k%orb(bra)%n, sps_k%orb(bra)%l, sps_k%orb(bra)%j, sps_k%orb(bra)%z
    end do

    cnt = 0
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(abs(sps_k%orb(bra)%n - sps_k%orb(ket)%n) > 1) cycle
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        if(sps_k%orb(bra)%z /= sps_k%orb(ket)%z) cycle
        cnt = cnt + 1
      end do
    end do

    write(15, '(a)') '####### one-body term'
    write(15, '(i5, i2)') cnt, 0
    write(15, '(a)') '### a, b, <a|t|b>'
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(abs(sps_k%orb(bra)%n - sps_k%orb(ket)%n) > 1) cycle
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        if(sps_k%orb(bra)%z /= sps_k%orb(ket)%z) cycle
        l = sps_k%orb(bra)%l
        n1= sps_k%orb(bra)%n
        n2= sps_k%orb(ket)%n
        kine = 0.d0
        if(n1 == n2) kine = dble(2 * n1 + l) + 1.5d0
        if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l) + 1.5d0))
        if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l) + 1.5d0))
        write(15,'(2i5,f16.8)') bra,ket, 0.5d0 * kine * two%GetFrequency() * (1.d0 - 1.d0 / dble(two%mass_snt))
      end do
    end do

    cnt = 0
    do ich = 1, two%GetNumberChannels()
      cnt = cnt + two%jpz(ich)%GetNumberStates() * (two%jpz(ich)%GetNumberStates() + 1) / 2
    end do
    write(15, '(i15, i2)') cnt, 0
    write(15, '(a)') '##### two-body term #####'
    write(15, '(a)') '## a, b, c, d, J, <ab|v|cd>'
    do ich = 1, two%GetNumberChannels()
      j = two%jpz(ich)%GetJ()
      n = two%jpz(ich)%GetNumberStates()
      if(n < 1) cycle
      do bra = 1, n
        a = two%jpz(ich)%n2label1(bra)
        b = two%jpz(ich)%n2label2(bra)
        a_k = sps_k%nljz2idx(sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j, sps%orb(a)%z)
        b_k = sps_k%nljz2idx(sps%orb(b)%n, sps%orb(b)%l, sps%orb(b)%j, sps%orb(b)%z)
        do ket = 1, bra
          c = two%jpz(ich)%n2label1(ket)
          d = two%jpz(ich)%n2label2(ket)
          c_k = sps_k%nljz2idx(sps%orb(c)%n, sps%orb(c)%l, sps%orb(c)%j, sps%orb(c)%z)
          d_k = sps_k%nljz2idx(sps%orb(d)%n, sps%orb(d)%l, sps%orb(d)%j, sps%orb(d)%z)
          aa = [sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j, sps%orb(a)%z]
          bb = [sps%orb(b)%n, sps%orb(b)%l, sps%orb(b)%j, sps%orb(b)%z]
          cc = [sps%orb(c)%n, sps%orb(c)%l, sps%orb(c)%j, sps%orb(c)%z]
          dd = [sps%orb(d)%n, sps%orb(d)%l, sps%orb(d)%j, sps%orb(d)%z]
          Tcm2 = p_dot_p(aa,bb,cc,dd,J) * two%GetFrequency() / dble(two%mass_snt)
          write(15, '(5i5, f16.8)') a_k, b_k, c_k, d_k, j, this%MatCh(ich,ich)%m(bra, ket) - Tcm2
        end do
      end do
    end do
    close(15)
    call sps_k%fin()
  end subroutine write_scalar_operator_kshell_A_snt

  subroutine write_scalar_operator_ascii_me2j(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: remap
    integer :: a, na, la, ja
    integer :: b, nb, lb, jb
    integer :: c, nc, lc, jc
    integer :: d, nd, ld, jd
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: dmax
    integer :: Jmin, Jmax, J
    real(8) :: norm_fact, me_00, me_pp, me_10, me_nn
    real(8), allocatable :: v_buffer(:)
    integer :: cnt_buf
    integer :: nbuf = 100000000
    character(255) :: header
    type(str) :: OpName

    ms => this%ms
    sps => this%ms%sps
    OpName = this%GetOpName()
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NN int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    allocate(v_buffer(4*nbuf))
    call remap%init(ms%GetEmax(), ms%GetLmax() )
    open(15, file = f%val, status = 'replace', form='formatted',access='stream')
    write(15,'(a)') trim(header)
    cnt_buf = 0
    do a = 1, remap%norbs
      na = remap%orb(a)%n; la = remap%orb(a)%l; ja = remap%orb(a)%j
      ap = sps%nljz2idx(na,la,ja,-1)
      an = sps%nljz2idx(na,la,ja, 1)
      do b = 1, a
        if(remap%orb(a)%e +remap%orb(b)%e > ms%GetE2max() ) cycle
        nb = remap%orb(b)%n; lb = remap%orb(b)%l; jb = remap%orb(b)%j
        bp = sps%nljz2idx(nb,lb,jb,-1)
        bn = sps%nljz2idx(nb,lb,jb, 1)
        do c = 1, a
          nc = remap%orb(c)%n; lc = remap%orb(c)%l; jc = remap%orb(c)%j
          cp = sps%nljz2idx(nc,lc,jc,-1)
          cn = sps%nljz2idx(nc,lc,jc, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            nd = remap%orb(d)%n; ld = remap%orb(d)%l; jd = remap%orb(d)%j
            dp = sps%nljz2idx(nd,ld,jd,-1)
            dn = sps%nljz2idx(nd,ld,jd, 1)
            if(remap%orb(c)%e + remap%orb(d)%e > ms%GetE2max() ) cycle
            if((-1)**(la+lb) /= (-1)**(lc+ld)) cycle
            Jmin = max(abs(remap%orb(a)%j-remap%orb(b)%j),abs(remap%orb(c)%j-remap%orb(d)%j))/2
            Jmax = min(remap%orb(a)%j+remap%orb(b)%j,remap%orb(c)%j+remap%orb(d)%j)/2
            if(Jmin > Jmax) cycle

            norm_fact = 1.d0
            if(a == b) norm_fact = dsqrt(2.d0) * norm_fact
            if(c == d) norm_fact = dsqrt(2.d0) * norm_fact
            do J = Jmin, Jmax
              me_pp = this%GetTBME(ap,bp,cp,dp,J,J) * norm_fact
              me_nn = this%GetTBME(an,bn,cn,dn,J,J) * norm_fact
              me_00 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) - ( &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) ))
              me_10 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) + &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) )
              v_buffer(cnt_buf+1) = me_00
              v_buffer(cnt_buf+2) = me_nn
              v_buffer(cnt_buf+3) = me_10
              v_buffer(cnt_buf+4) = me_pp
              cnt_buf = cnt_buf + 4
              if(cnt_buf == 4*nbuf) then
                write(15,'(10f16.8)') v_buffer(:)
                !write(15,'(4f16.8)') v_buffer(:)
                cnt_buf = 0
              end if
            end do
          end do
        end do
      end do
    end do
    if(cnt_buf /= 0) then
      write(15,'(10f16.8)') v_buffer(:cnt_buf)
      !write(15,'(4f16.8)') v_buffer(:cnt_buf)
    end if
    close(15)
    call remap%fin()
    deallocate(v_buffer)
  end subroutine write_scalar_operator_ascii_me2j

  subroutine write_scalar_operator_ascii_mfd(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: remap
    integer :: a, na, la, ja
    integer :: b, nb, lb, jb
    integer :: c, nc, lc, jc
    integer :: d, nd, ld, jd
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: dmax
    integer :: Jmin, Jmax, J, cnt, loop
    real(8) :: norm_fact, me_00, me_pp, me_10, me_nn, t_rel, h_ho_rel
    type(str) :: OpName
    type(sys) :: s
    type(TwoBodyLabOp) :: trel, hrel


    ms => this%ms
    sps => this%ms%sps
    OpName = this%GetOpName()
    call trel%init(ms,s%str("Kinetic"))
    call set_two_body_kin2(trel)
    call hrel%init(ms,s%str("HOHamil"))
    call set_two_body_hrel2(hrel)

    call remap%init(ms%GetEmax(), ms%GetLmax() )
    open(15, file = f%val, status = 'replace', form='formatted',access='stream')
    cnt = 0
    do loop = 1, 2
      if(loop==2) write(15,'(i12,i4,i4,f6.2)') cnt, ms%GetEmax(), ms%GetE2max(), ms%GetFrequency()
      cnt = 0
      do a = 1, remap%norbs
        na = remap%orb(a)%n; la = remap%orb(a)%l; ja = remap%orb(a)%j
        ap = sps%nljz2idx(na,la,ja,-1)
        an = sps%nljz2idx(na,la,ja, 1)
        do b = 1, a
          if(remap%orb(a)%e +remap%orb(b)%e > ms%GetE2max() ) cycle
          nb = remap%orb(b)%n; lb = remap%orb(b)%l; jb = remap%orb(b)%j
          bp = sps%nljz2idx(nb,lb,jb,-1)
          bn = sps%nljz2idx(nb,lb,jb, 1)
          do c = 1, a
            nc = remap%orb(c)%n; lc = remap%orb(c)%l; jc = remap%orb(c)%j
            cp = sps%nljz2idx(nc,lc,jc,-1)
            cn = sps%nljz2idx(nc,lc,jc, 1)
            dmax = c
            if(a == c) dmax = b
            do d = 1, dmax
              nd = remap%orb(d)%n; ld = remap%orb(d)%l; jd = remap%orb(d)%j
              dp = sps%nljz2idx(nd,ld,jd,-1)
              dn = sps%nljz2idx(nd,ld,jd, 1)
              if(remap%orb(c)%e + remap%orb(d)%e > ms%GetE2max() ) cycle
              if((-1)**(la+lb) /= (-1)**(lc+ld)) cycle
              Jmin = max(abs(remap%orb(a)%j-remap%orb(b)%j),abs(remap%orb(c)%j-remap%orb(d)%j))/2
              Jmax = min(remap%orb(a)%j+remap%orb(b)%j,remap%orb(c)%j+remap%orb(d)%j)/2
              if(Jmin > Jmax) cycle
              norm_fact = 1.d0
              if(a == b) norm_fact = 1.d0/sqrt(2.d0) * norm_fact
              if(c == d) norm_fact = 1.d0/sqrt(2.d0) * norm_fact

              do J = Jmin, Jmax
                me_pp = this%GetTBME(ap,bp,cp,dp,J,J) 
                me_nn = this%GetTBME(an,bn,cn,dn,J,J) 
                me_00 = 0.5d0 * norm_fact * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                    & this%GetTBME(an,bp,cn,dp,J,J) - ( &
                    & this%GetTBME(ap,bn,cn,dp,J,J) + &
                    & this%GetTBME(an,bp,cp,dn,J,J) ))
                me_10 = 0.5d0 * norm_fact * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                    & this%GetTBME(an,bp,cn,dp,J,J) + &
                    & this%GetTBME(ap,bn,cn,dp,J,J) + &
                    & this%GetTBME(an,bp,cp,dn,J,J) )
                if(abs(me_00) > 1.d-16) then
                  ! a, b, c, d, J, T, t_rel, ho_rel, vcoul, me_pn, me_pp, me_nn
                  t_rel = 0.5d0 * norm_fact * (trel%GetTBME(ap,bn,cp,dn,J,J) + &
                    & trel%GetTBME(an,bp,cn,dp,J,J) - ( &
                    & trel%GetTBME(ap,bn,cn,dp,J,J) + &
                    & trel%GetTBME(an,bp,cp,dn,J,J) )) / ms%GetFrequency()
                  h_ho_rel = 0.5d0 * norm_fact * (hrel%GetTBME(ap,bn,cp,dn,J,J) + &
                    & hrel%GetTBME(an,bp,cn,dp,J,J) - ( &
                    & hrel%GetTBME(ap,bn,cn,dp,J,J) + &
                    & hrel%GetTBME(an,bp,cp,dn,J,J) )) / ms%GetFrequency()
                  if(loop==2) write(15,'(6i4,6f16.8)') a, b, c, d, J, 0, &
                    & t_rel, h_ho_rel, 0.d0, me_00, 0.d0, 0.d0
                  cnt = cnt + 1
                end if

                if(abs(me_pp) > 1.d-16 .or. abs(me_nn) > 1.d-16 .or. abs(me_10) > 1.d-16) then
                  t_rel = 0.5d0 * norm_fact * (trel%GetTBME(ap,bn,cp,dn,J,J) + &
                    & trel%GetTBME(an,bp,cn,dp,J,J) + ( &
                    & trel%GetTBME(ap,bn,cn,dp,J,J) + &
                    & trel%GetTBME(an,bp,cp,dn,J,J) )) / ms%GetFrequency()
                  h_ho_rel = 0.5d0 * norm_fact * (hrel%GetTBME(ap,bn,cp,dn,J,J) + &
                    & hrel%GetTBME(an,bp,cn,dp,J,J) + ( &
                    & hrel%GetTBME(ap,bn,cn,dp,J,J) + &
                    & hrel%GetTBME(an,bp,cp,dn,J,J) )) / ms%GetFrequency()
                  if(loop==2) write(15,'(6i4,6f16.8)') a, b, c, d, J, 1, &
                    & t_rel, h_ho_rel, 0.d0, me_10, me_pp, me_nn
                  cnt = cnt + 1
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    close(15)
    call remap%fin()
  end subroutine write_scalar_operator_ascii_mfd

  subroutine write_scalar_operator_gzip_me2j(this, f)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: remap
    integer :: a, na, la, ja
    integer :: b, nb, lb, jb
    integer :: c, nc, lc, jc
    integer :: d, nd, ld, jd
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: dmax
    integer :: Jmin, Jmax, J
    real(8) :: norm_fact, me_00, me_pp, me_10, me_nn
    real(8), allocatable :: v_buffer(:)
    integer :: cnt_buf, line, nrest
    !integer :: nbuf = 100000000
    integer :: nbuf = 100
    character(256) :: header
    type(c_ptr) :: fp, err
    character(kind=c_char, len=160) :: buffer
    character(12) :: cfmt
    type(str) :: OpName

    ms => this%ms
    sps => this%ms%sps
    OpName = this%GetOpName()
    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') header = 'NN int.'
    header = trim(header) // ' calculated by NuHamil (Tokyo code)'
#ifdef VERSION
    header = trim(header) // ', ' // trim(VERSION)
#endif
    allocate(v_buffer(4*nbuf))
    call remap%init(ms%GetEmax(), ms%GetLmax() )
    fp = gzip_open(f%val,"wt")
    err = gzip_writeline(fp, trim(header), len_trim(header))
    cnt_buf = 0
    do a = 1, remap%norbs
      na = remap%orb(a)%n; la = remap%orb(a)%l; ja = remap%orb(a)%j
      ap = sps%nljz2idx(na,la,ja,-1)
      an = sps%nljz2idx(na,la,ja, 1)
      do b = 1, a
        if(remap%orb(a)%e +remap%orb(b)%e > ms%GetE2max() ) cycle
        nb = remap%orb(b)%n; lb = remap%orb(b)%l; jb = remap%orb(b)%j
        bp = sps%nljz2idx(nb,lb,jb,-1)
        bn = sps%nljz2idx(nb,lb,jb, 1)
        do c = 1, a
          nc = remap%orb(c)%n; lc = remap%orb(c)%l; jc = remap%orb(c)%j
          cp = sps%nljz2idx(nc,lc,jc,-1)
          cn = sps%nljz2idx(nc,lc,jc, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            nd = remap%orb(d)%n; ld = remap%orb(d)%l; jd = remap%orb(d)%j
            dp = sps%nljz2idx(nd,ld,jd,-1)
            dn = sps%nljz2idx(nd,ld,jd, 1)
            if(remap%orb(c)%e + remap%orb(d)%e > ms%GetE2max() ) cycle
            if((-1)**(la+lb) /= (-1)**(lc+ld)) cycle
            Jmin = max(abs(remap%orb(a)%j-remap%orb(b)%j),abs(remap%orb(c)%j-remap%orb(d)%j))/2
            Jmax = min(remap%orb(a)%j+remap%orb(b)%j,remap%orb(c)%j+remap%orb(d)%j)/2
            if(Jmin > Jmax) cycle

            norm_fact = 1.d0
            if(a == b) norm_fact = dsqrt(2.d0) * norm_fact
            if(c == d) norm_fact = dsqrt(2.d0) * norm_fact
            do J = Jmin, Jmax
              me_pp = this%GetTBME(ap,bp,cp,dp,J,J) * norm_fact
              me_nn = this%GetTBME(an,bn,cn,dn,J,J) * norm_fact
              me_00 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) - ( &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) ))
              me_10 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) + &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) )
              v_buffer(cnt_buf+1) = me_00
              v_buffer(cnt_buf+2) = me_nn
              v_buffer(cnt_buf+3) = me_10
              v_buffer(cnt_buf+4) = me_pp
              cnt_buf = cnt_buf + 4
              if(cnt_buf == 4*nbuf) then
                do line = 1, 4*nbuf/10
                  write(buffer,'(10f16.8)') v_buffer( 10*(line-1)+1:10*line)
                  err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
                end do
                cnt_buf = 0
              end if
            end do
          end do
        end do
      end do
    end do
    if(cnt_buf /= 0) then

      do line = 1, cnt_buf/10
        write(buffer,'(10f16.8)') v_buffer( 10*(line-1)+1 : 10*line )
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end do

      nrest = cnt_buf - (cnt_buf/10) * 10
      if(nrest > 0) then
        cfmt = '(xf16.8)'
        write(cfmt(2:2),'(i1)') nrest
        write(buffer,cfmt) v_buffer((cnt_buf/10)*10+1:cnt_buf)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end if

    end if
    err = gzip_close(fp)
    call remap%fin()
    deallocate(v_buffer)
  end subroutine write_scalar_operator_gzip_me2j

  subroutine write_scalar_operator_binary_me2j(this, f)
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: remap
    integer :: a, na, la, ja
    integer :: b, nb, lb, jb
    integer :: c, nc, lc, jc
    integer :: d, nd, ld, jd
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: dmax
    integer :: Jmin, Jmax, J
    real(8) :: norm_fact, me_00, me_pp, me_10, me_nn
#ifdef single_precision
    real(4), allocatable :: v_buffer(:)
#else
    real(8), allocatable :: v_buffer(:)
#endif
    integer :: cnt_buf
    integer :: nbuf = 100000000
    !character(255) :: header

    ms => this%ms
    sps => this%ms%sps
    allocate(v_buffer(4*nbuf))
    call remap%init(ms%GetEmax(), ms%GetLmax() )
    open(15, file = f%val, status = 'replace', form = 'unformatted', access = 'stream')
    cnt_buf = 0
    do a = 1, remap%norbs
      na = remap%orb(a)%n; la = remap%orb(a)%l; ja = remap%orb(a)%j
      ap = sps%nljz2idx(na,la,ja,-1)
      an = sps%nljz2idx(na,la,ja, 1)
      do b = 1, a
        if(remap%orb(a)%e +remap%orb(b)%e > ms%GetE2max() ) cycle
        nb = remap%orb(b)%n; lb = remap%orb(b)%l; jb = remap%orb(b)%j
        bp = sps%nljz2idx(nb,lb,jb,-1)
        bn = sps%nljz2idx(nb,lb,jb, 1)
        do c = 1, a
          nc = remap%orb(c)%n; lc = remap%orb(c)%l; jc = remap%orb(c)%j
          cp = sps%nljz2idx(nc,lc,jc,-1)
          cn = sps%nljz2idx(nc,lc,jc, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            nd = remap%orb(d)%n; ld = remap%orb(d)%l; jd = remap%orb(d)%j
            dp = sps%nljz2idx(nd,ld,jd,-1)
            dn = sps%nljz2idx(nd,ld,jd, 1)
            if(remap%orb(c)%e +remap%orb(d)%e > ms%GetE2max() ) cycle
            if((-1)**(la+lb) /= (-1)**(lc+ld)) cycle
            Jmin = max(abs(remap%orb(a)%j-remap%orb(b)%j),abs(remap%orb(c)%j-remap%orb(d)%j))/2
            Jmax = min(remap%orb(a)%j+remap%orb(b)%j,remap%orb(c)%j+remap%orb(d)%j)/2
            if(Jmin > Jmax) cycle

            norm_fact = 1.d0
            if(a == b) norm_fact = dsqrt(2.d0) * norm_fact
            if(c == d) norm_fact = dsqrt(2.d0) * norm_fact
            do J = Jmin, Jmax
              me_pp = this%GetTBME(ap,bp,cp,dp,J,J) * norm_fact
              me_nn = this%GetTBME(an,bn,cn,dn,J,J) * norm_fact
              me_00 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) - ( &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) ))
              me_10 = 0.5d0 * (this%GetTBME(ap,bn,cp,dn,J,J) + &
                  & this%GetTBME(an,bp,cn,dp,J,J) + &
                  & this%GetTBME(ap,bn,cn,dp,J,J) + &
                  & this%GetTBME(an,bp,cp,dn,J,J) )
#ifdef single_precision
              v_buffer(cnt_buf+1) = real(me_00)
              v_buffer(cnt_buf+2) = real(me_nn)
              v_buffer(cnt_buf+3) = real(me_10)
              v_buffer(cnt_buf+4) = real(me_pp)
#else
              v_buffer(cnt_buf+1) = me_00
              v_buffer(cnt_buf+2) = me_nn
              v_buffer(cnt_buf+3) = me_10
              v_buffer(cnt_buf+4) = me_pp
#endif
              cnt_buf = cnt_buf + 4
              if(cnt_buf == 4*nbuf) then
                write(15) v_buffer
                cnt_buf = 0
              end if
            end do
          end do
        end do
      end do
    end do
    if(cnt_buf/=0) then
      write(15) v_buffer(:cnt_buf)
    end if
    close(15)
    call remap%fin()
    deallocate(v_buffer)
  end subroutine write_scalar_operator_binary_me2j

  subroutine write_general_operator_ascii(this, f)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(OrbitsIsospin) :: remap
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, Jab, Jcd
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od
    character(255) :: header
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, me_pnnn, me_npnp, me_npnn, me_nnnn
    type(str) :: OpName
    ms => this%ms
    call remap%init(ms%GetEmax(), ms%GetLmax() )
    open(15, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // trim(OpName%val) // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(15,"(a)") trim(header)
    write(15,"(a)") "# J, Parity, Z"
    write(15,"(3i3)") this%GetOpJ(), this%GetOpP(), this%GetOpZ()
    write(15,"(a)") "# Model space definition:"
    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      write(15,"(4i4)") a, oa%n, oa%l, oa%j
    end do
    write(15,"(a)") "# Zero Body term: 0.00000000000"
    write(15,"(a)") "# Number of one-body ME"
    write(15,"(a)") "# a, b, pp, nn, np, pn"
    write(15,*) count_general_me_onebody(this, remap)
    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      do b = 1, remap%norbs
        ob => remap%GetOrbit(b)
        bp = remap%iso2pn(ms%sps,b,-1)
        bn = remap%iso2pn(ms%sps,b, 1)
        if((-1)**(oa%l+ob%l) * this%GetOpP() /= 1) cycle
        if( triag(oa%j,ob%j,2*this%GetOpJ() ) ) cycle
        write(15,"(2i4,4f16.8)") a, b, this%OneBody%m(ap,bp), &
            & this%OneBody%m(an,bn), this%OneBody%m(an,bp), this%OneBody%m(ap,bn)
      end do
    end do

    write(15,"(a)") "# Number of two-body ME"
    write(15,"(a)") "# a, b, c, d, Jab, Jcd, pppp, pppn, ppnp, ppnn, pnpn, pnnp, pnnn, npnp, npnn, nnnn"
    write(15,*) count_general_me_twobody(this, remap)
    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      do b = 1, a
        ob => remap%GetOrbit(b)
        bp = remap%iso2pn(ms%sps,b,-1)
        bn = remap%iso2pn(ms%sps,b, 1)
        if( oa%e + ob%e > ms%GetE2max() ) cycle

        do c = 1, remap%norbs
          oc => remap%GetOrbit(c)
          cp = remap%iso2pn(ms%sps,c,-1)
          cn = remap%iso2pn(ms%sps,c, 1)
          do d = 1, c
            od => remap%GetOrbit(d)
            dp = remap%iso2pn(ms%sps,d,-1)
            dn = remap%iso2pn(ms%sps,d, 1)
            if( oc%e + od%e > ms%GetE2max() ) cycle
            if((-1)**(oa%l+ob%l+oc%l+od%l) * this%GetOpP() /= 1) cycle

            do Jab = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              do Jcd = abs(oc%j-od%j)/2, (oc%j+od%j)/2
                if( triag(Jab,Jcd,this%GetOpJ() ) ) cycle
                me_pppp = this%GetTBME(ap,bp,cp,dp,Jab,Jcd)
                me_pppn = this%GetTBME(ap,bp,cp,dn,Jab,Jcd)
                me_ppnp = this%GetTBME(ap,bp,cn,dp,Jab,Jcd)
                me_ppnn = this%GetTBME(ap,bp,cn,dn,Jab,Jcd)
                me_pnpn = this%GetTBME(ap,bn,cp,dn,Jab,Jcd)
                me_pnnp = this%GetTBME(ap,bn,cn,dp,Jab,Jcd)
                me_pnnn = this%GetTBME(ap,bn,cn,dn,Jab,Jcd)
                me_npnp = this%GetTBME(an,bp,cn,dp,Jab,Jcd)
                me_npnn = this%GetTBME(an,bp,cn,dn,Jab,Jcd)
                me_nnnn = this%GetTBME(an,bn,cn,dn,Jab,Jcd)
                write(15,"(6i4,10es16.8)") a, b, c, d, Jab, Jcd, &
                    & me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, &
                    & me_pnnn, me_npnp, me_npnn, me_nnnn
              end do
            end do
          end do
        end do
      end do
    end do
    close(15)
    call remap%fin()
  end subroutine write_general_operator_ascii

  subroutine write_general_operator_kshell_ascii(this, f)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits) :: sps_k
    integer :: a, b, c, d, Jab, Jcd, a_k, b_k, c_k, d_k
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    character(255) :: header
    integer :: wunit=25, cnt=0
    integer :: chbra, chket, bra, ket, max_ket
    type(str) :: OpName
    logical :: is_scalar
    real(8) :: me
    ms => this%ms

    is_scalar = .false.
    if(this%GetOpJ()==0 .and. this%GetOpP()==1 .and. this%GetOpZ()==0) is_scalar=.true.
    call sps_k%init(ms%sps%emax, mode="kshell")
    open(wunit, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // trim(OpName%val) // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(wunit,"(a)") trim(header)
    write(wunit,"(a)") "# Zero Body term: 0.00000000000"
    write(wunit,"(a)") "# Model space definition:"
    write(wunit,'(4i5)') sps_k%norbs/2, sps_k%norbs/2, 0, 0
    do a = 1, sps_k%norbs
      oa => sps_k%GetOrbit(a)
      write(wunit,"(5i4)") a, oa%n, oa%l, oa%j, oa%z
    end do
    write(wunit,"(a)") "# reduced matrix elements "
    cnt = 0
    do a = 1, sps_k%norbs
      do b = 1, sps_k%norbs
        if( abs( this%OneBody%m(a,b) ) < 1.d-16 ) cycle
        cnt = cnt + 1
      end do
    end do
    write(wunit,"(2i5,f6.2)") cnt, 0, ms%GetFrequency()
    do a = 1, sps_k%norbs
      oa => ms%sps%GetOrbit(a)
      a_k = sps_k%nljz2idx(oa%n, oa%l, oa%j, oa%z)
      do b = 1, sps_k%norbs
        ob => ms%sps%GetOrbit(b)
        b_k = sps_k%nljz2idx(ob%n, ob%l, ob%j, ob%z)
        if( abs( this%OneBody%m(a,b) ) < 1.d-16 ) cycle
        cnt = cnt + 1
        write(wunit,"(2i3,es16.8)") a_k, b_k, this%OneBody%m(a,b)
      end do
    end do

    cnt = 0
    do chbra = 1, ms%GetNumberChannels()
      do chket = 1, ms%GetNumberChannels()
        if( this%MatCh(chbra,chket)%zero ) cycle

        do bra = 1, ms%jpz(chbra)%GetNumberStates()
          max_ket = ms%jpz(chket)%GetNumberStates()
          if(chbra==chket) max_ket = bra
          do ket = 1, max_ket
            if( abs( this%MatCh(chbra,chket)%m(bra,ket) ) < 1.d-16 ) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do

    write(wunit,"(i8,i3,f6.2)") cnt, 0, ms%GetFrequency()
    do chbra = 1, ms%GetNumberChannels()
      Jab = ms%jpz(chbra)%GetJ()
      do chket = 1, ms%GetNumberChannels()
        Jcd = ms%jpz(chket)%GetJ()
        max_ket = ms%jpz(chket)%GetNumberStates()
        if( this%MatCh(chbra,chket)%zero ) cycle
        do bra = 1, ms%jpz(chbra)%GetNumberStates()

          if(chbra==chket) max_ket = bra
          do ket = 1, max_ket
            if( abs( this%MatCh(chbra,chket)%m(bra,ket) ) < 1.d-16 ) cycle
            a = ms%jpz(chbra)%n2label1(bra)
            b = ms%jpz(chbra)%n2label2(bra)
            c = ms%jpz(chket)%n2label1(ket)
            d = ms%jpz(chket)%n2label2(ket)
            oa => ms%sps%GetOrbit(a)
            ob => ms%sps%GetOrbit(b)
            oc => ms%sps%GetOrbit(c)
            od => ms%sps%GetOrbit(d)
            me = this%GetTBME(a,b,c,d,Jab,Jcd)
            a_k = sps_k%nljz2idx(oa%n, oa%l, oa%j, oa%z)
            b_k = sps_k%nljz2idx(ob%n, ob%l, ob%j, ob%z)
            c_k = sps_k%nljz2idx(oc%n, oc%l, oc%j, oc%z)
            d_k = sps_k%nljz2idx(od%n, od%l, od%j, od%z)
            if(is_scalar .and. .not. this%reduced_me()) then
              write(wunit,"(5i4,es16.8)") a_k, b_k, c_k, d_k, Jab, me
            else
              write(wunit,"(6i4,es16.8)") a_k, b_k, c_k, d_k, Jab, Jcd, me
            end if
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine write_general_operator_kshell_ascii

  subroutine write_general_operator_gzip(this, f)
    use MyLibrary, only: triag
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_writeline, gzip_close
    class(TwoBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(OrbitsIsospin) :: remap
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, Jab, Jcd
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od
    character(255) :: header
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, me_pnnn, me_npnp, me_npnn, me_nnnn
    character(kind=c_char, len=100) :: bufferm
    character(kind=c_char, len=80) :: buffer1
    character(kind=c_char, len=184) :: buffer2
    type(c_ptr) :: fp, err
    type(str) :: OpName
    ms => this%ms
    call remap%init(ms%GetEmax(), ms%GetLmax() )
    OpName = this%GetOpName()
    header = "# " // trim(OpName%val) // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    fp = gzip_open(f%val,"wt")
    err = gzip_writeline(fp, trim(header), len_trim(header))

    write(bufferm,"(6i4)") this%GetOpJ(), this%GetOpP(), this%GetOpZ(), ms%GetEmax(), ms%GetE2max(), ms%GetLmax()
    err = gzip_writeline(fp, trim(bufferm), len_trim(bufferm))

    write(bufferm,'(f16.8)') 0.d0
    err = gzip_writeline(fp, trim(bufferm), len_trim(bufferm))

    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      do b = 1, remap%norbs
        ob => remap%GetOrbit(b)
        bp = remap%iso2pn(ms%sps,b,-1)
        bn = remap%iso2pn(ms%sps,b, 1)
        if((-1)**(oa%l+ob%l) * this%GetOpP() /= 1) cycle
        if( triag(oa%j,ob%j,2*this%GetOpJ() ) ) cycle
        write(buffer1,"(4es16.6)") this%OneBody%m(ap,bp), &
            & this%OneBody%m(an,bn), this%OneBody%m(an,bp), this%OneBody%m(ap,bn)
        err = gzip_writeline(fp, trim(buffer1), len_trim(buffer1))
      end do
    end do

    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      do b = 1, a
        ob => remap%GetOrbit(b)
        bp = remap%iso2pn(ms%sps,b,-1)
        bn = remap%iso2pn(ms%sps,b, 1)
        if( oa%e + ob%e > ms%GetE2max() ) cycle

        do c = 1, remap%norbs
          oc => remap%GetOrbit(c)
          cp = remap%iso2pn(ms%sps,c,-1)
          cn = remap%iso2pn(ms%sps,c, 1)
          do d = 1, c
            od => remap%GetOrbit(d)
            dp = remap%iso2pn(ms%sps,d,-1)
            dn = remap%iso2pn(ms%sps,d, 1)
            if( oc%e + od%e > ms%GetE2max() ) cycle
            if((-1)**(oa%l+ob%l+oc%l+od%l) * this%GetOpP() /= 1) cycle

            do Jab = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              do Jcd = abs(oc%j-od%j)/2, (oc%j+od%j)/2
                if( triag(Jab,Jcd,this%GetOpJ() ) ) cycle
                me_pppp = this%GetTBME(ap,bp,cp,dp,Jab,Jcd)
                me_pppn = this%GetTBME(ap,bp,cp,dn,Jab,Jcd)
                me_ppnp = this%GetTBME(ap,bp,cn,dp,Jab,Jcd)
                me_ppnn = this%GetTBME(ap,bp,cn,dn,Jab,Jcd)
                me_pnpn = this%GetTBME(ap,bn,cp,dn,Jab,Jcd)
                me_pnnp = this%GetTBME(ap,bn,cn,dp,Jab,Jcd)
                me_pnnn = this%GetTBME(ap,bn,cn,dn,Jab,Jcd)
                me_npnp = this%GetTBME(an,bp,cn,dp,Jab,Jcd)
                me_npnn = this%GetTBME(an,bp,cn,dn,Jab,Jcd)
                me_nnnn = this%GetTBME(an,bn,cn,dn,Jab,Jcd)
                write(buffer2,"(10es16.6)") &
                    & me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, &
                    & me_pnnn, me_npnp, me_npnn, me_nnnn
                err = gzip_writeline(fp, trim(buffer2), len_trim(buffer2))
              end do
            end do
          end do
        end do
      end do
    end do
    err = gzip_close(fp)
    call remap%fin()
  end subroutine write_general_operator_gzip

  !
  !  (1/2m) * p_{i} \cdot p_{j} / hw
  !
  function p_dot_p(a, b, c, d, J) result(r)
    use MyLibrary, only: sjs
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ja, jb, jc, jd
    real(8) :: fact

    r = 0.d0
    ja = a(3)
    jb = b(3)
    jc = c(3)
    jd = d(3)

    fact = 1.d0
    if(a(1) == b(1) .and. a(2) == b(2) .and. a(3) == b(3) .and. a(4) == b(4)) fact = fact / dsqrt(2.d0)
    if(c(1) == d(1) .and. c(2) == d(2) .and. c(3) == d(3) .and. c(4) == d(4)) fact = fact / dsqrt(2.d0)

    r = - (-1.d0) ** ((jb+jc)/2+J) * &
        & sjs(ja,jb,2*J,jd,jc,2) * red_nab_j(a,c) * red_nab_j(b,d) - &
        & (-1.d0) ** ((jb+jc)/2) * &
        & sjs(ja,jb,2*J,jc,jd,2) * red_nab_j(a,d) * red_nab_j(b,c)

    r = r * fact
  end function p_dot_p

  !
  !  (1/2) * m * omega**2 * r_{i} \cdot r_{j} / hw
  !
  function r_dot_r(a, b, c, d, J) result(r)
    use MyLibrary, only: sjs
    integer, intent(in) :: a(4), b(4), c(4), d(4), J
    real(8) :: r
    integer :: ja, jb, jc, jd
    integer :: ea, eb, ec, ed
    real(8) :: fact

    r = 0.d0
    ea = 2*a(1) + a(2)
    eb = 2*b(1) + b(2)
    ec = 2*c(1) + c(2)
    ed = 2*d(1) + d(2)
    if(abs(ea+eb-ec-ed) > 2 .or. mod(abs(ea+eb-ec-ed),2)==1) return

    ja = a(3)
    jb = b(3)
    jc = c(3)
    jd = d(3)

    fact = 1.d0
    if(a(1) == b(1) .and. a(2) == b(2) .and. a(3) == b(3) .and. a(4) == b(4)) fact = fact / dsqrt(2.d0)
    if(c(1) == d(1) .and. c(2) == d(2) .and. c(3) == d(3) .and. c(4) == d(4)) fact = fact / dsqrt(2.d0)

    r = (-1.d0) ** ((jb+jc)/2+J) * &
        & sjs(ja,jb,2*J,jd,jc,2) * red_r_j(a,c) * red_r_j(b,d) + &
        & (-1.d0) ** ((jb+jc)/2) * &
        & sjs(ja,jb,2*J,jc,jd,2) * red_r_j(a,d) * red_r_j(b,c)

    r = r * fact
  end function r_dot_r

  function red_nab_j(a, b) result(r)
    use MyLibrary, only: red_nab_l, sjs
    real(8) :: r
    integer, intent(in) :: a(4), b(4)
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb

    r = 0.d0

    na = a(1)
    la = a(2)
    ja = a(3)
    za = a(4)

    nb = b(1)
    lb = b(2)
    jb = b(3)
    zb = b(4)

    if(za /= zb) return

    r = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_nab_l(na, la, nb, lb)
  end function red_nab_j

  ! < a || r || b >
  function red_r_j(a, b) result(r)
    use MyLibrary, only: red_r_l, sjs
    real(8) :: r
    integer, intent(in) :: a(4), b(4)
    integer :: na, la, ja, za
    integer :: nb, lb, jb, zb

    r = 0.d0

    na = a(1)
    la = a(2)
    ja = a(3)
    za = a(4)

    nb = b(1)
    lb = b(2)
    jb = b(3)
    zb = b(4)

    if(za /= zb) return

    r = (-1.d0) ** ((3+2*la+jb)/2) * dsqrt(dble(ja + 1) * dble(jb + 1)) * &
        &  sjs(ja, 2, jb, 2 * lb, 1, 2 * la) * red_r_l(na, la, nb, lb)
  end function red_r_j

  subroutine ReadOperator(this, f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyLabOp before calling TalmiMoshinskyTransofrmation"
      return
    end if

    ti = omp_get_wtime()

    if(.not. this%reduced_me() ) call read_scalar_operator(this, f)
    if( this%reduced_me() ) call read_general_operator(this, f)
    call timer%Add(s%str('Read from file'), omp_get_wtime() - ti)
  end subroutine ReadOperator

  subroutine read_scalar_operator(this, f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str(".myg.bin"))) then
      call read_scalar_operator_binary_myg(this, f)
      return
    end if

    if(s%find(f,s%str(".myg"))) then
      call read_scalar_operator_ascii_myg(this, f)
      return
    end if

    if(s%find(f,s%str(".Akshell.snt"))) then
      call read_scalar_operator_kshell_A_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".kshell.snt"))) then
      call read_scalar_operator_kshell_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".snt"))) then
      call read_scalar_operator_ascii_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j.bin"))) then
      call read_scalar_operator_binary_me2j(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j"))) then
      call read_scalar_operator_ascii_me2j(this, f)
      return
    end if

    write(*,"(2a)") "# Unknown input file format: ", trim(f%val)
    stop
  end subroutine read_scalar_operator

  subroutine read_general_operator(this, f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str(".op.me2j.gz"))) then
      call read_general_operator_gzip(this, f)
      return
    end if

    if(s%find(f,s%str(".op.me2j"))) then
      call read_general_operator_ascii(this, f)
      return
    end if

    if(s%find(f,s%str(".op.snt"))) then
      call read_general_operator_kshell_ascii(this, f)
      return
    end if

    write(*,"(2a)") "# Unknown input file format: ", trim(f%val)
    stop
  end subroutine read_general_operator

  subroutine read_scalar_operator_ascii_myg(this, f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, bra, ket
    integer :: j, ch
    integer :: cnt, i
    real(8) :: v, ph
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    ms => this%ms
    open(15, file = f%val, action = 'read')
    read(15,*) cnt
    read(15,*)
    do i = 1, cnt
      read(15,*) a, b, c, d, J, v
      oa => ms%sps%orb(a)
      ob => ms%sps%orb(b)
      oc => ms%sps%orb(c)
      od => ms%sps%orb(d)
      if(oa%z + ob%z /= oc%z + od%z) cycle
      if((-1) ** (oa%l + ob%l + oc%l + od%l) /= 1) cycle
      ch = ms%jpz2idx(J, (-1)**(oa%l+ob%l), (oa%z+ob%z)/2)
      bra = ms%jpz(ch)%GetIndex(a,b)
      ket = ms%jpz(ch)%GetIndex(c,d)
      ph = ms%jpz(ch)%GetPhase(a,b) * ms%jpz(ch)%GetPhase(c,d)
      this%MatCh(ch,ch)%m(bra,ket) = v * ph
    end do
    close(15)
  end subroutine read_scalar_operator_ascii_myg

  subroutine read_scalar_operator_binary_myg(this, f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: cnt
    integer, allocatable :: ia(:), ib(:), ic(:), id(:), jj(:)
    real(8), allocatable :: v2(:)
    integer :: i, a, b, c, d, J, ch, bra, ket
    real(8) :: v, ph
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    ms => this%ms
    open(15, file = f%val, action = 'read', form = 'unformatted', access = 'stream')
    read(15) cnt
    allocate(ia(cnt), ib(cnt), ic(cnt), id(cnt), jj(cnt), v2(cnt))
    read(15) ia, ib, ic, id, jj, v2
    close(15)

    do i = 1, cnt
      a = ia(i)
      b = ib(i)
      c = ic(i)
      d = id(i)
      J = jj(i)
      v = v2(i)
      oa => ms%sps%orb(a)
      ob => ms%sps%orb(b)
      oc => ms%sps%orb(c)
      od => ms%sps%orb(d)
      if(oa%z + ob%z /= oc%z + od%z) cycle
      if((-1) ** (oa%l + ob%l + oc%l + od%l) /= 1) cycle
      ch = ms%jpz2idx(J, (-1)**(oa%l+ob%l), (oa%z+ob%z)/2)
      bra = ms%jpz(ch)%GetIndex(a,b)
      ket = ms%jpz(ch)%GetIndex(c,d)
      ph = ms%jpz(ch)%GetPhase(a,b) * ms%jpz(ch)%GetPhase(c,d)
      this%MatCh(ch,ch)%m(bra,ket) = v * ph
    end do

    deallocate(ia)
    deallocate(ib)
    deallocate(ic)
    deallocate(id)
    deallocate(jj, v2)
  end subroutine read_scalar_operator_binary_myg

  subroutine read_scalar_operator_ascii_snt(this, f)
    use MyLibrary, only: skip_comment
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    type(Orbits), target :: sps_snt
    integer :: porbs, norbs, pc, nc
    integer :: line, lines
    integer :: idx, n, l, j, z
    real(8) :: zero_body, me
    integer :: aa, bb, cc, dd
    integer :: a, b, c, d, ch, bra, ket
    real(8) :: ph
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    two => this%ms
    sps => this%ms%sps

    call sps_snt%init(sps%emax,sps%lmax)
    open(15, file = f%val, action = 'read')
    call skip_comment(15,'#')
    read(15,*) porbs, norbs, pc, nc

    do line = 1, porbs+norbs
      read(15,*) idx, n, l, j, z
      call sps_snt%orb(idx)%set(n,l,j,z,idx)
    end do

    call skip_comment(15,'#')
    read(15,*) zero_body
    call skip_comment(15,'#')
    read(15,*) lines
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) a, b, me ! one-body part
    end do

    call skip_comment(15,'#')
    read(15,*) lines
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) aa, bb, cc, dd, J, me
      oa => sps_snt%orb(aa)
      ob => sps_snt%orb(bb)
      oc => sps_snt%orb(cc)
      od => sps_snt%orb(dd)
      if(oa%z + ob%z /= oc%z + od%z) cycle
      if((-1) ** (oa%l + ob%l + oc%l + od%l) /= 1) cycle
      a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
      b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
      c = sps%nljz2idx(oc%n, oc%l, oc%j, oc%z)
      d = sps%nljz2idx(od%n, od%l, od%j, od%z)
      ch = two%jpz2idx(J, (-1)**(oa%l+ob%l), (oa%z+ob%z)/2)
      bra = two%jpz(ch)%GetIndex(a,b)
      ket = two%jpz(ch)%GetIndex(c,d)
      ph = two%jpz(ch)%GetPhase(a,b) * two%jpz(ch)%GetPhase(c,d)
      this%MatCh(ch,ch)%m(bra,ket) = me * ph
    end do
    close(15)

    call sps_snt%fin()
  end subroutine read_scalar_operator_ascii_snt

  subroutine read_scalar_operator_kshell_snt(this, f)
    use MyLibrary, only: skip_comment
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    type(Orbits), target :: sps_snt
    integer :: porbs, norbs, pc, nc
    integer :: line, lines
    integer :: idx, n, l, j, z
    real(8) :: zero_body, me, pdp
    integer :: aa, bb, cc, dd
    integer :: a, b, c, d, ch, bra, ket
    real(8) :: ph
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    two => this%ms
    sps => this%ms%sps

    call sps_snt%init(sps%emax,sps%lmax)
    open(15, file = f%val, action = 'read')
    call skip_comment(15,'#')
    read(15,*) porbs, norbs, pc, nc

    do line = 1, porbs+norbs
      read(15,*) idx, n, l, j, z
      call sps_snt%orb(idx)%set(n,l,j,z,idx)
    end do

    call skip_comment(15,'#')
    read(15,*) zero_body
    call skip_comment(15,'#')
    read(15,*) lines
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) a, b, me ! one-body part
    end do

    call skip_comment(15,'#')
    read(15,*) lines
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) aa, bb, cc, dd, J, me, pdp
      oa => sps_snt%orb(aa)
      ob => sps_snt%orb(bb)
      oc => sps_snt%orb(cc)
      od => sps_snt%orb(dd)
      if(oa%z + ob%z /= oc%z + od%z) cycle
      if((-1) ** (oa%l + ob%l + oc%l + od%l) /= 1) cycle
      a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
      b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
      c = sps%nljz2idx(oc%n, oc%l, oc%j, oc%z)
      d = sps%nljz2idx(od%n, od%l, od%j, od%z)
      ch = two%jpz2idx(J, (-1)**(oa%l+ob%l), (oa%z+ob%z)/2)
      bra = two%jpz(ch)%GetIndex(a,b)
      ket = two%jpz(ch)%GetIndex(c,d)
      ph = two%jpz(ch)%GetPhase(a,b) * two%jpz(ch)%GetPhase(c,d)
      this%MatCh(ch,ch)%m(bra,ket) = me * ph
    end do
    close(15)
    call sps_snt%fin()
  end subroutine read_scalar_operator_kshell_snt

  subroutine read_scalar_operator_kshell_A_snt(this, f)
    use MyLibrary, only: skip_comment
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: two
    type(Orbits), pointer :: sps
    type(Orbits), target :: sps_snt
    integer :: porbs, norbs, pc, nc
    integer :: line, lines
    integer :: idx, n, l, j, z, id
    real(8) :: me, v, Tcm2, ph
    integer :: aa, bb, cc, dd
    integer :: a, b, c, d, ch, bra, ket
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    two => this%ms
    sps => this%ms%sps

    call sps_snt%init(sps%emax,sps%lmax)
    open(15, file = f%val, action = 'read')
    call skip_comment(15,'#')
    read(15,*) porbs, norbs, pc, nc

    do line = 1, porbs+norbs
      read(15,*) idx, n, l, j, z
      call sps_snt%orb(idx)%set(n,l,j,z,idx)
    end do

    call skip_comment(15,'#')
    read(15,*) lines, id
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) a, b, me ! one-body part
    end do

    call skip_comment(15,'#')
    read(15,*) lines, id
    call skip_comment(15,'#')
    do line = 1, lines
      read(15,*) aa, bb, cc, dd, J, me
      oa => sps_snt%orb(aa)
      ob => sps_snt%orb(bb)
      oc => sps_snt%orb(cc)
      od => sps_snt%orb(dd)
      if(oa%z + ob%z /= oc%z + od%z) cycle
      if((-1) ** (oa%l + ob%l + oc%l + od%l) /= 1) cycle
      a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
      b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
      c = sps%nljz2idx(oc%n, oc%l, oc%j, oc%z)
      d = sps%nljz2idx(od%n, od%l, od%j, od%z)
      ch = two%jpz2idx(J, (-1)**(oa%l+ob%l), (oa%z+ob%z)/2)
      bra = two%jpz(ch)%GetIndex(a,b)
      ket = two%jpz(ch)%GetIndex(c,d)
      ph = two%jpz(ch)%GetPhase(a,b) * two%jpz(ch)%GetPhase(c,d)
      Tcm2 = p_dot_p([oa%n, oa%l, oa%j, oa%z], [ob%n, ob%l, ob%j, ob%z], &
          & [oc%n, oc%l, oc%j, oc%z], [od%n, od%l, od%j, od%z], j) * dble(two%mass_snt) / two%GetFrequency()
      v = me + Tcm2
      this%MatCh(ch,ch)%m(bra,ket) = v * ph
    end do
    close(15)
    call sps_snt%fin()
  end subroutine read_scalar_operator_kshell_A_snt

  subroutine read_scalar_operator_ascii_me2j(this,f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: sps_me2j
    integer :: nelm
    integer :: io, runit = 22
    real(8), allocatable :: v(:)
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact
    type(sys) :: s

    ms => this%ms
    sps => this%ms%sps

    call sps_me2j%init(sps%emax, sps%lmax)
    nelm = count_scalar_me2j(sps_me2j,ms%GetE2max() )
    allocate(v(nelm))
    if( s%find(f, s%str('.gz')) ) then
      call get_vector_me2j_gz(f,v)
    end if

    if( .not. s%find(f, s%str('.gz')) ) then
      open(runit, file=f%val, action='read',iostat=io)
      if(io /= 0) then
        write(*,'(2a)') 'File open error: ', trim(f%val)
        return
      end if
      call get_vector_me2j_formatted(runit,v)
      close(runit)
    end if

    icnt = 0
    do a = 1, sps_me2j%norbs
      la = sps_me2j%orb(a)%l
      ja = sps_me2j%orb(a)%j
      ea = sps_me2j%orb(a)%e
      ap = sps_me2j%iso2pn(sps,a,-1)
      an = sps_me2j%iso2pn(sps,a, 1)
      do b = 1, a
        lb = sps_me2j%orb(b)%l
        jb = sps_me2j%orb(b)%j
        eb = sps_me2j%orb(b)%e
        bp = sps_me2j%iso2pn(sps,b,-1)
        bn = sps_me2j%iso2pn(sps,b, 1)
        if(ea + eb > ms%GetE2max() ) cycle
        do c = 1, a
          lc = sps_me2j%orb(c)%l
          jc = sps_me2j%orb(c)%j
          ec = sps_me2j%orb(c)%e
          cp = sps_me2j%iso2pn(sps,c,-1)
          cn = sps_me2j%iso2pn(sps,c, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps_me2j%orb(d)%l
            jd = sps_me2j%orb(d)%j
            ed = sps_me2j%orb(d)%e
            dp = sps_me2j%iso2pn(sps,d,-1)
            dn = sps_me2j%iso2pn(sps,d, 1)
            if(ec + ed > ms%GetE2max() ) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              me_00 = v(icnt+1)
              me_nn = v(icnt+2)
              me_10 = v(icnt+3)
              me_pp = v(icnt+4)
              icnt = icnt + 4

              if(a == b .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(a == b .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_pp)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_pp)

              if(ap == 0) cycle
              if(bp == 0) cycle
              if(cp == 0) cycle
              if(dp == 0) cycle

              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%GetE2max() ) cycle
              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%GetE2max() ) cycle
#ifdef NOperatorsDebug
              write(*,'(5i3,4f12.6)') a, b, c, d, J, &
                  &  me_00, me_nn, me_10, me_pp
#endif
              fact = 1.d0
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                call this%SetTBME(ap,bp,cp,dp,J,me_pp*fact)
                call this%SetTBME(an,bn,cn,dn,J,me_nn*fact)
                call this%AddToTBME(ap,bn,cp,dn,J,0.5d0*me_10) ! pnpn
                if(c/=d) call this%AddToTBME(ap,bn,cn,dp,J,0.5d0*me_10) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call this%AddToTBME(an,bp,cn,dp,J,0.5d0*me_10) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call this%AddToTBME(an,bp,cp,dn,J,0.5d0*me_10) ! nppn
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call this%AddToTBME(ap,bn,cp,dn,J,0.5d0*me_00) ! pnpn
                if(c/=d) call this%AddToTBME(ap,bn,cn,dp,J,-0.5d0*me_00) ! pnnp
                if(a/=b .and. c/=d) &
                    &    call this%AddToTBME(an,bp,cn,dp,J, 0.5d0*me_00) ! npnp
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call this%AddToTBME(an,bp,cp,dn,J,-0.5d0*me_00) ! nppn
              end if
            end do
          end do
        end do
      end do
    end do

    deallocate(v)
    call sps_me2j%fin()
    return
  end subroutine read_scalar_operator_ascii_me2j

  subroutine read_scalar_operator_binary_me2j(this,f)
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer ::sps
    type(OrbitsIsospin) :: sps_me2j
    integer :: nelm, io, runit = 22
#ifdef single_precision
    real(4), allocatable :: v(:)
#else
    real(8), allocatable :: v(:)
#endif
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax, icnt
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    real(8) :: me_00, me_pp, me_10, me_nn, fact

    ms => this%ms
    sps => this%ms%sps

    call sps_me2j%init(sps%emax, sps%lmax)
    nelm = count_scalar_me2j(sps_me2j,ms%GetE2max())
    allocate(v(nelm))
    open(runit, file=f%val, action='read',iostat=io, &
        & form='unformatted', access='stream')
    if(io /= 0) then
      write(*,'(2a)') 'File open error: ', trim(f%val)
      return
    end if
    read(runit) v
    close(runit)

    icnt = 0
    do a = 1, sps_me2j%norbs
      la = sps_me2j%orb(a)%l
      ja = sps_me2j%orb(a)%j
      ea = sps_me2j%orb(a)%e
      ap = sps_me2j%iso2pn(sps,a,-1)
      an = sps_me2j%iso2pn(sps,a, 1)
      do b = 1, a
        lb = sps_me2j%orb(b)%l
        jb = sps_me2j%orb(b)%j
        eb = sps_me2j%orb(b)%e
        bp = sps_me2j%iso2pn(sps,b,-1)
        bn = sps_me2j%iso2pn(sps,b, 1)
        if(ea + eb > ms%GetE2max()) cycle
        do c = 1, a
          lc = sps_me2j%orb(c)%l
          jc = sps_me2j%orb(c)%j
          ec = sps_me2j%orb(c)%e
          cp = sps_me2j%iso2pn(sps,c,-1)
          cn = sps_me2j%iso2pn(sps,c, 1)
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps_me2j%orb(d)%l
            jd = sps_me2j%orb(d)%j
            ed = sps_me2j%orb(d)%e
            dp = sps_me2j%iso2pn(sps,d,-1)
            dn = sps_me2j%iso2pn(sps,d, 1)
            if(ec + ed > ms%GetE2max()) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              me_00 = v(icnt+1)
              me_nn = v(icnt+2)
              me_10 = v(icnt+3)
              me_pp = v(icnt+4)
              icnt = icnt + 4

              if(a == b .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_nn) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(a == b .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)
              if(c == d .and. mod(J,2) == 1 .and. abs(me_pp) > 1.d-8) call me2j_read_warning(a,b,c,d,J,me_nn)

              if(ap == 0) cycle
              if(bp == 0) cycle
              if(cp == 0) cycle
              if(dp == 0) cycle

              if(sps%orb(ap)%e + sps%orb(bp)%e > ms%GetE2max() ) cycle
              if(sps%orb(cp)%e + sps%orb(dp)%e > ms%GetE2max() ) cycle

              fact = 1.d0
              if(a == b) fact = fact / dsqrt(2.d0)
              if(c == d) fact = fact / dsqrt(2.d0)

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 0) then
                call this%SetTBME(ap,bp,cp,dp,J,me_pp*fact)
                call this%SetTBME(an,bn,cn,dn,J,me_nn*fact)
                call this%AddToTBME(ap,bn,cp,dn,J,0.5d0*me_10)
                if(c/=d) call this%AddToTBME(ap,bn,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. c/=d) &
                    &    call this%AddToTBME(an,bp,cn,dp,J,0.5d0*me_10)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call this%AddToTBME(an,bp,cp,dn,J,0.5d0*me_10)
              end if

              if(abs(1.d0-fact) < 1.d-4 .or. mod(J,2) == 1) then
                call this%AddToTBME(ap,bn,cp,dn,J,0.5d0*me_00)
                if(c/=d) call this%AddToTBME(ap,bn,cn,dp,J,-0.5d0*me_00)
                if(a/=b .and. c/=d) &
                    &    call this%AddToTBME(an,bp,cn,dp,J, 0.5d0*me_00)
                if(a/=b .and. (a/=c .or. b/=d)) &
                    &    call this%AddToTBME(an,bp,cp,dn,J,-0.5d0*me_00)
              end if
            end do
          end do
        end do
      end do
    end do
    deallocate(v)
    call sps_me2j%fin()
    return
  end subroutine read_scalar_operator_binary_me2j

  subroutine me2j_read_warning(a,b,c,d,J,me)
    integer, intent(in) :: a, b, c, d, J
    real(8), intent(in) :: me
    write(*,'(a,5i3,f12.6)') "Warning: this TBME should be zero: ", a, b, c, d, J, me
  end subroutine me2j_read_warning

  function count_scalar_me2j(sps,e2max) result(r)
    type(OrbitsIsospin), intent(in) :: sps
    integer, intent(in) :: e2max
    integer :: r
    integer :: a, b, c, d, dmax
    integer :: la, ja, ea
    integer :: lb, jb, eb
    integer :: lc, jc, ec
    integer :: ld, jd, ed
    integer :: J, Jmin, Jmax
    r = 0
    do a = 1, sps%norbs
      la = sps%orb(a)%l
      ja = sps%orb(a)%j
      ea = sps%orb(a)%e
      do b = 1, a
        lb = sps%orb(b)%l
        jb = sps%orb(b)%j
        eb = sps%orb(b)%e
        if(ea + eb > e2max) cycle
        do c = 1, a
          lc = sps%orb(c)%l
          jc = sps%orb(c)%j
          ec = sps%orb(c)%e
          dmax = c
          if(a == c) dmax = b
          do d = 1, dmax
            ld = sps%orb(d)%l
            jd = sps%orb(d)%j
            ed = sps%orb(d)%e
            if(ec + ed > e2max) cycle
            if( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
            Jmin = max(abs(ja-jb), abs(jc-jd))/2
            Jmax = min(   (ja+jb),    (jc+jd))/2
            if(Jmin > Jmax) cycle
            do J = Jmin, Jmax
              r = r + 4
            end do
          end do
        end do
      end do
    end do
  end function count_scalar_me2j

  subroutine get_vector_me2j_formatted(ut,v)
    integer, intent(in) :: ut
    real(8), intent(inout) :: v(:)
    integer :: nelm, lines, i
    nelm = size(v)
    lines = nelm/10
    read(ut,*) ! header
    do i = 1, lines
      read(ut,*) v( (i-1)*10+1 : i*10)
    end do
    if(mod(nelm,10) == 0) return
    read(ut,*) v(lines*10+1:nelm)
  end subroutine get_vector_me2j_formatted

  subroutine get_vector_me2j_gz(f,v)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    type(str), intent(in) :: f
    real(8), intent(inout) :: v(:)
    integer :: nelm, line, lines, nelm_tail, cnt
    character(512) :: buffer = ''
    type(c_ptr) :: p, buf
    nelm = size(v)
    lines = nelm/10
    nelm_tail = nelm - 10 * lines
    p = gzip_open(f%val,'r')
    buf = gzip_readline(p, buffer, len(buffer))
    cnt = 0
    do line = 1, lines
      buf = gzip_readline(p, buffer, len(buffer))
      read(buffer,*) v(cnt+1:cnt+10)
      cnt = cnt + 10
    end do
    if(nelm_tail == 0) then
      buf = gzip_close(p)
      return
    end if

    buf = gzip_readline(p, buffer, len(buffer))
    read(buffer,*) v(lines*10+1:nelm)
    buf = gzip_close(p)
  end subroutine get_vector_me2j_gz

  subroutine SetTwBME_tensor(this,i1,i2,i3,i4,J12,J34,me)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(inout) :: this
    type(Orbits), pointer :: sps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket
    real(8) :: phase

    ms => this%ms
    sps => ms%sps

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(triag(J12, J34, this%GetOpJ() )) then
      write(*,*) "Warning: in SetTwBME_general: J"
      return
    end if

    if(P12 * P34 * this%GetOpP() /= 1) then
      write(*,*) "Warning: in SetTwBME_general: P"
      return
    end if

    if(abs(Z12 - Z34) - this%GetOpZ() /= 0) then
      write(*,*) "Warning: in SetTwBME_general: Tz"
      return
    end if

    chbra = ms%GetIndex(J12,P12,Z12)
    chket = ms%GetIndex(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%GetIndex(i1,i2)
    ket = ms%jpz(chket)%GetIndex(i3,i4)
    if(bra * ket == 0) return

    phase = ms%jpz(chbra)%GetPhase(i1,i2) * &
        &   ms%jpz(chket)%GetPhase(i3,i4)
    if( chbra < chket ) then
      this%MatCh(chket,chbra)%m(ket,bra) = phase * me
      return
    end if
    this%MatCh(chbra,chket)%m(bra,ket) = phase * me
  end subroutine SetTwBME_tensor

  subroutine SetTwBME_scalar(this,i1,i2,i3,i4,J,me)
    class(TwoBodyLabOp), intent(inout) :: this
    type(Orbits), pointer :: sps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket
    real(8) :: phase

    ms => this%ms
    sps => ms%sps

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in SetTwBME_general: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in SetTwBME_general: Tz"
      return
    end if

    ch = ms%GetIndex(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%GetIndex(i1,i2)
    ket = ms%jpz(ch)%GetIndex(i3,i4)
    if(bra * ket == 0) return

    phase = ms%jpz(ch)%GetPhase(i1,i2) * &
        &   ms%jpz(ch)%GetPhase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = phase * me
    this%MatCh(ch,ch)%m(ket,bra) = phase * me
  end subroutine SetTwBME_scalar

  subroutine AddToTwBME_tensor(this,i1,i2,i3,i4,J12,J34,me)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(inout) :: this
    type(Orbits), pointer :: sps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer, intent(in) :: i1,i2,i3,i4,J12,J34
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: chbra, chket, bra, ket
    real(8) :: phase

    ms => this%ms
    sps => ms%sps

    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(triag(J12, J34, this%GetOpJ())) then
      write(*,*)"Warning: in AddToTwBME_general: J"
      return
    end if

    if(P12 * P34 * this%GetOpP() /= 1) then
      write(*,*) "Warning: in AddToTwBME_general: P"
      return
    end if

    if(abs(Z12 - Z34) - this%GetOpZ() /= 0) then
      write(*,*) "Warning: in AddToTwBME_general: Tz"
      return
    end if

    chbra = ms%GetIndex(J12,P12,Z12)
    chket = ms%GetIndex(J34,P34,Z34)
    if(chbra * chket == 0) return

    bra = ms%jpz(chbra)%GetIndex(i1,i2)
    ket = ms%jpz(chket)%GetIndex(i3,i4)
    if(bra * ket == 0) return

    phase = ms%jpz(chbra)%GetPhase(i1,i2) * &
        &   ms%jpz(chket)%GetPhase(i3,i4)

    if( chbra < chket ) then
      ! Maybe double count --
      write(*,*) "Maybe better to make sure not double-counting: ", __LINE__, " at ", __FILE__
      this%MatCh(chket,chbra)%m(ket,bra) = &
          & this%MatCh(chket,chbra)%m(ket,bra) + phase * me
      return
    end if
    this%MatCh(chbra,chket)%m(bra,ket) = &
        & this%MatCh(chbra,chket)%m(bra,ket) + phase * me
  end subroutine AddToTwBME_tensor

  subroutine AddToTwBME_scalar(this,i1,i2,i3,i4,J,me)
    class(TwoBodyLabOp), intent(inout) :: this
    type(Orbits), pointer :: sps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer, intent(in) :: i1,i2,i3,i4,J
    real(8), intent(in) :: me
    integer :: P12, P34
    integer :: Z12, Z34
    integer :: ch, bra, ket
    real(8) :: phase

    ms => this%ms
    sps => ms%sps
    P12 = (-1) ** (sps%orb(i1)%l + sps%orb(i2)%l)
    P34 = (-1) ** (sps%orb(i3)%l + sps%orb(i4)%l)
    Z12 = (sps%orb(i1)%z + sps%orb(i2)%z) / 2
    Z34 = (sps%orb(i3)%z + sps%orb(i4)%z) / 2

    if(P12 * P34 /= 1) then
      write(*,*) "Warning: in AddToTwBME_general: P"
      return
    end if

    if(Z12 - Z34 /= 0) then
      write(*,*) "Warning: in AddToTwBME_general: Tz"
      return
    end if

    ch = ms%GetIndex(J,P12,Z12)
    if(ch == 0) return

    bra = ms%jpz(ch)%GetIndex(i1,i2)
    ket = ms%jpz(ch)%GetIndex(i3,i4)
    if(bra * ket == 0) return

    phase = ms%jpz(ch)%GetPhase(i1,i2) * &
        &   ms%jpz(ch)%GetPhase(i3,i4)

    this%MatCh(ch,ch)%m(bra,ket) = this%MatCh(ch,ch)%m(bra,ket) + &
        & me * phase
    this%MatCh(ch,ch)%m(ket,bra) = this%MatCh(ch,ch)%m(bra,ket)
  end subroutine AddToTwBME_scalar

  function count_general_me_onebody(this, isps) result(r)
    use MyLibrary, only: triag
    type(TwoBodyLabOp), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: isps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: r
    integer :: a, b
    type(SingleParticleOrbitIsospin), pointer :: oa, ob
    ms => this%ms
    r = 0
    do a = 1, isps%norbs
      oa => isps%GetOrbit(a)
      do b = 1, isps%norbs
        ob => isps%GetOrbit(b)
        if((-1)**(oa%l+ob%l) * this%GetOpP() /= 1) cycle
        if( triag(oa%j,ob%j,2*this%GetOpJ() ) ) cycle
        r = r + 1
      end do
    end do
  end function count_general_me_onebody

  function count_general_me_twobody(this, isps) result(r)
    use MyLibrary, only: triag
    type(TwoBodyLabOp), intent(in) :: this
    type(OrbitsIsospin), intent(in) :: isps
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: r
    integer :: a, b, c, d, Jab, Jcd
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od
    ms => this%ms
    r = 0
    do a = 1, isps%norbs
      oa => isps%GetOrbit(a)
      do b = 1, a
        ob => isps%GetOrbit(b)

        do c = 1, isps%norbs
          oc => isps%GetOrbit(c)
          do d = 1, c
            od => isps%GetOrbit(d)
            if((-1)**(oa%l+ob%l+oc%l+od%l) * this%GetOpP() /= 1) cycle

            do Jab = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              do Jcd = abs(oc%j-od%j)/2, (oc%j+od%j)/2
                if( triag(Jab,Jcd,this%GetOpJ() ) ) cycle
                r = r + 1
              end do
            end do
          end do
        end do
      end do
    end do
  end function count_general_me_twobody

  subroutine read_general_operator_ascii(this, f)
    use MyLibrary, only: triag
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(OrbitsIsospin) :: remap
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: rankJ, rankP, rankZ
    integer :: a, b, c, d, Jab, Jcd
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od
    character(255) :: line
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: nline1, nline2, iline
    real(8) :: me_pp, me_nn, me_np, me_pn
    real(8) :: me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, me_pnnn, me_npnp, me_npnn, me_nnnn
    ms => this%ms
    call remap%init( ms%GetEmax(), ms%GetLmax() )
    open(15, file = f%val, status = 'replace')
    read(15,*) line
    read(15,*) line
    read(15,*) rankJ, rankP, rankZ
    if( this%GetOpJ() /= rankJ ) stop "Error: rank mismatch J"
    if( this%GetOpP() /= rankP ) stop "Error: rank mismatch P"
    if( this%GetOpZ() /= rankZ ) stop "Error: rank mismatch Z"
    read(15,*) line
    do a = 1, remap%norbs
      read(15,*) line
    end do
    read(15,*) line
    read(15,*) line
    read(15,*) line
    read(15,*) nline1
    do iline = 1, nline1
      read(15,*) a, b, me_pp, me_nn, me_np, me_pn
      oa => remap%GetOrbit(a)
      ob => remap%GetOrbit(b)
      if(oa%e > ms%GetEmax() ) cycle
      if(ob%e > ms%GetEmax() ) cycle
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      bp = remap%iso2pn(ms%sps,b,-1)
      bn = remap%iso2pn(ms%sps,b, 1)
      this%OneBody%m(ap,bp) = me_pp
      this%OneBody%m(an,bn) = me_nn
      this%OneBody%m(an,bp) = me_np
      this%OneBody%m(ap,bn) = me_pn
    end do

    read(15,*) line
    read(15,*) line
    read(15,*) nline2
    do iline = 1, nline2
      read(15,*) a, b, c, d, Jab, Jcd, &
          & me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, &
          & me_pnnn, me_npnp, me_npnn, me_nnnn
      oa => remap%GetOrbit(a)
      ob => remap%GetOrbit(b)
      oc => remap%GetOrbit(c)
      od => remap%GetOrbit(d)
      if( oa%e > ms%GetEmax() ) cycle
      if( ob%e > ms%GetEmax() ) cycle
      if( oc%e > ms%GetEmax() ) cycle
      if( od%e > ms%GetEmax() ) cycle
      if( oa%e+ob%e > ms%GetE2max() ) cycle
      if( oc%e+od%e > ms%GetE2max() ) cycle
      ap = remap%iso2pn(ms%sps,a,-1)
      an = remap%iso2pn(ms%sps,a, 1)
      bp = remap%iso2pn(ms%sps,b,-1)
      bn = remap%iso2pn(ms%sps,b, 1)
      cp = remap%iso2pn(ms%sps,c,-1)
      cn = remap%iso2pn(ms%sps,c, 1)
      dp = remap%iso2pn(ms%sps,d,-1)
      dn = remap%iso2pn(ms%sps,d, 1)
      if(abs(me_pppp) > 1.d-10) call this%SetTBME(ap,bp,cp,dp,Jab,Jcd,me_pppp)
      if(abs(me_pppn) > 1.d-10) call this%SetTBME(ap,bp,cp,dn,Jab,Jcd,me_pppn)
      if(abs(me_ppnp) > 1.d-10) call this%SetTBME(ap,bp,cn,dp,Jab,Jcd,me_ppnp)
      if(abs(me_ppnn) > 1.d-10) call this%SetTBME(ap,bp,cn,dn,Jab,Jcd,me_ppnn)
      if(abs(me_pnpn) > 1.d-10) call this%SetTBME(ap,bn,cp,dn,Jab,Jcd,me_pnpn)
      if(abs(me_pnnp) > 1.d-10) call this%SetTBME(ap,bn,cn,dp,Jab,Jcd,me_pnnp)
      if(abs(me_pnnn) > 1.d-10) call this%SetTBME(ap,bn,cn,dn,Jab,Jcd,me_pnnn)
      if(abs(me_npnp) > 1.d-10) call this%SetTBME(an,bp,cn,dp,Jab,Jcd,me_npnp)
      if(abs(me_npnn) > 1.d-10) call this%SetTBME(an,bp,cn,dn,Jab,Jcd,me_npnn)
      if(abs(me_nnnn) > 1.d-10) call this%SetTBME(an,bn,cn,dn,Jab,Jcd,me_nnnn)
    end do
    close(15)
    call remap%fin()
  end subroutine read_general_operator_ascii

  subroutine read_general_operator_gzip(this, f)
    use MyLibrary, only: triag
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_readline, gzip_close
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(OrbitsIsospin) :: remap
    type(TwoBodyLabPNSpace), pointer :: ms
    integer :: a, b, c, d, Jab, Jcd, iline, cnt
    type(SingleParticleOrbitIsospin), pointer :: oa, ob, oc, od
    integer :: ap, bp, cp, dp
    integer :: an, bn, cn, dn
    integer :: rankJ, rankP, rankZ, emax_read, e2max_read, lmax_read
    real(8) :: me_pp, me_nn, me_np, me_pn
    real(8) :: me_pppp, me_pppn, me_ppnp, me_ppnn, me_pnpn, me_pnnp, me_pnnn, me_npnp, me_npnn, me_nnnn
    character(512) :: line = ""
    type(c_ptr) :: fp, err
    real(8), allocatable :: me1(:), me2(:)
    ms => this%ms

    fp = gzip_open(f%val,"r")
    line = ""
    err = gzip_readline(fp, line, len(line)) ! header

    line = ""
    err = gzip_readline(fp, line, len(line)) ! J, P, Z, emax, e2max
    read(line,*) rankJ, rankP, rankZ, emax_read, e2max_read, lmax_read
    if( this%GetOpJ() /= rankJ ) stop "Error: rank mismatch J"
    if( this%GetOpP() /= rankP ) stop "Error: rank mismatch P"
    if( this%GetOpZ() /= rankZ ) stop "Error: rank mismatch Z"
    call remap%init( emax_read, lmax_read )

    err = gzip_readline(fp, line, len(line)) ! Zero-body
    allocate(me1( 4*count_general_me_onebody( this, remap )))
    allocate(me2( 10*count_general_me_twobody( this, remap )))
    do iline = 1, count_general_me_onebody( this, remap )
      err = gzip_readline(fp, line, len(line))
      read(line,*) me1(4*(iline-1)+1:  4*iline)
    end do
    do iline = 1, count_general_me_twobody( this, remap )
      err = gzip_readline(fp, line, len(line))
      read(line,*) me2(10*(iline-1)+1:  10*iline)
    end do
    err = gzip_close(fp)


    cnt = 0
    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      do b = 1, remap%norbs
        ob => remap%GetOrbit(b)
        if((-1)**(oa%l+ob%l) * this%GetOpP() /= 1) cycle
        if( triag(oa%j,ob%j,2*this%GetOpJ() ) ) cycle
        me_pp = me1(cnt+1)
        me_nn = me1(cnt+2)
        me_np = me1(cnt+3)
        me_pn = me1(cnt+4)
        cnt = cnt + 4
        if(oa%e > ms%GetEmax() ) cycle
        if(ob%e > ms%GetEmax() ) cycle
        ap = remap%iso2pn(ms%sps,a,-1)
        an = remap%iso2pn(ms%sps,a, 1)
        bp = remap%iso2pn(ms%sps,b,-1)
        bn = remap%iso2pn(ms%sps,b, 1)
        this%OneBody%m(ap,bp) = me_pp
        this%OneBody%m(an,bn) = me_nn
        this%OneBody%m(an,bp) = me_np
        this%OneBody%m(ap,bn) = me_pn
      end do
    end do

    cnt = 0
    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      if( oa%e > ms%GetEmax() ) exit
      do b = 1, a
        ob => remap%GetOrbit(b)
        if( oa%e + ob%e > e2max_read ) cycle

        do c = 1, remap%norbs
          oc => remap%GetOrbit(c)
          do d = 1, c
            od => remap%GetOrbit(d)
            if( oc%e + od%e > e2max_read ) cycle
            if((-1)**(oa%l+ob%l+oc%l+od%l) * this%GetOpP() /= 1) cycle

            do Jab = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              do Jcd = abs(oc%j-od%j)/2, (oc%j+od%j)/2
                if( triag(Jab,Jcd,this%GetOpJ() ) ) cycle
                me_pppp = me2(cnt+1)
                me_pppn = me2(cnt+2)
                me_ppnp = me2(cnt+3)
                me_ppnn = me2(cnt+4)
                me_pnpn = me2(cnt+5)
                me_pnnp = me2(cnt+6)
                me_pnnn = me2(cnt+7)
                me_npnp = me2(cnt+8)
                me_npnn = me2(cnt+9)
                me_nnnn = me2(cnt+10)
                cnt = cnt + 10
                if( oa%e > ms%GetEmax() ) cycle
                if( ob%e > ms%GetEmax() ) cycle
                if( oc%e > ms%GetEmax() ) cycle
                if( od%e > ms%GetEmax() ) cycle
                if( oa%e+ob%e > ms%GetE2max() ) cycle
                if( oc%e+od%e > ms%GetE2max() ) cycle
                ap = remap%iso2pn(ms%sps,a,-1)
                an = remap%iso2pn(ms%sps,a, 1)
                bp = remap%iso2pn(ms%sps,b,-1)
                bn = remap%iso2pn(ms%sps,b, 1)
                cp = remap%iso2pn(ms%sps,c,-1)
                cn = remap%iso2pn(ms%sps,c, 1)
                dp = remap%iso2pn(ms%sps,d,-1)
                dn = remap%iso2pn(ms%sps,d, 1)
                if(abs(me_pppp) > 1.d-10) call this%SetTBME(ap,bp,cp,dp,Jab,Jcd,me_pppp)
                if(abs(me_pppn) > 1.d-10) call this%SetTBME(ap,bp,cp,dn,Jab,Jcd,me_pppn)
                if(abs(me_ppnp) > 1.d-10) call this%SetTBME(ap,bp,cn,dp,Jab,Jcd,me_ppnp)
                if(abs(me_ppnn) > 1.d-10) call this%SetTBME(ap,bp,cn,dn,Jab,Jcd,me_ppnn)
                if(abs(me_pnpn) > 1.d-10) call this%SetTBME(ap,bn,cp,dn,Jab,Jcd,me_pnpn)
                if(abs(me_pnnp) > 1.d-10) call this%SetTBME(ap,bn,cn,dp,Jab,Jcd,me_pnnp)
                if(abs(me_pnnn) > 1.d-10) call this%SetTBME(ap,bn,cn,dn,Jab,Jcd,me_pnnn)
                if(abs(me_npnp) > 1.d-10) call this%SetTBME(an,bp,cn,dp,Jab,Jcd,me_npnp)
                if(abs(me_npnn) > 1.d-10) call this%SetTBME(an,bp,cn,dn,Jab,Jcd,me_npnn)
                if(abs(me_nnnn) > 1.d-10) call this%SetTBME(an,bn,cn,dn,Jab,Jcd,me_nnnn)
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate( me1, me2 )
    call remap%fin()
  end subroutine read_general_operator_gzip

  subroutine read_general_operator_kshell_ascii( this, f )
    class(TwoBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    write(*,*) trim(f%val)
    write(*,*) this%is_init
    write(*,*) "not implemented yet: return 0"
    return
  end subroutine read_general_operator_kshell_ascii

  subroutine SubtractOneBody(this)
    class(TwoBodyLabOp), intent(inout) :: this
    type(TwoBodyLabPNSpace), pointer :: ms
    type(Orbits), pointer :: sps
    integer :: chbra, chket, bra, ket
    integer :: a, b, c, d, Jab, Jcd
    type(SingleParticleOrbit), pointer :: oa, ob
    real(8) :: me
    ms => this%ms
    sps => ms%sps
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      do b =1, sps%norbs
        ob => sps%GetOrbit(b)

        this%OneBody%m(a,b) = CalcMEOneBody( this%OperatorDef, &
            & [oa%n,oa%l,oa%j,oa%z], [ob%n,ob%l,ob%j,ob%z] )
      end do
    end do

    do chbra = 1, ms%GetNumberChannels()
      Jab = ms%jpz(chbra)%GetJ()
      do chket = 1, ms%GetNumberChannels()
        Jcd = ms%jpz(chket)%GetJ()
        if( this%MatCh(chbra,chket)%zero ) cycle
        do bra = 1, ms%jpz(chbra)%GetNumberStates()
          a = ms%jpz(chbra)%n2label1(bra)
          b = ms%jpz(chbra)%n2label2(bra)
          do ket = 1, ms%jpz(chket)%GetNumberStates()
            c = ms%jpz(chket)%n2label1(ket)
            d = ms%jpz(chket)%n2label2(ket)
            me = embedding_onebody_in_twobody( this%OneBody, sps, a, b, c, d, Jab, Jcd, this%GetOpJ(), this%reduced_me() )
            !if(abs(this%MatCh(chbra,chket)%m(bra,ket)-me) > 1.d-4 ) write(*,*) this%MatCh(chbra,chket)%m(bra,ket), me
            this%MatCh(chbra,chket)%m(bra,ket) = this%MatCh(chbra,chket)%m(bra,ket) - me
          end do
        end do
      end do
    end do

  end subroutine SubtractOneBody

  function embedding_onebody_in_twobody( OneBody, sps, a, b, c, d, Jab, Jcd, Lambda, reduced ) result(r)
    use MyLibrary, only: sjs
    type(DMat), intent(in) :: OneBody
    type(Orbits), intent(in), target :: sps
    integer, intent(in) :: a, b, c, d, Jab, Jcd, Lambda
    logical, intent(in) :: reduced
    integer :: ja, jb, jc, jd
    type(SingleParticleOrbit), pointer :: oa, ob, oc, od
    real(8) :: r

    oa => sps%GetOrbit(a)
    ob => sps%GetOrbit(b)
    oc => sps%GetOrbit(c)
    od => sps%GetOrbit(d)
    ja = oa%j
    jb = ob%j
    jc = oc%j
    jd = od%j

    r = 0.d0
    if( (.not. reduced) .and. Lambda==0 ) then
      if(b==d) r = r + OneBody%m(a,c)
      if(a==c) r = r + OneBody%m(b,d)
      if(a==d) r = r - OneBody%m(b,c) * (-1.d0)**( (ja+jb)/2-Jab )
      if(b==c) r = r - OneBody%m(a,d) * (-1.d0)**( (ja+jb)/2-Jab )
      if(a==b) r = r / sqrt(2.d0)
      if(c==d) r = r / sqrt(2.d0)
      return
    end if

    if(b==d) r = r + OneBody%m(a,c) * (-1.d0)**( (ja+jb)/2+Jcd ) * sjs(2*Jab,2*Jcd,2*Lambda,jc,ja,jb)
    if(a==c) r = r + OneBody%m(b,d) * (-1.d0)**( (jc+jd)/2-Jab ) * sjs(2*Jab,2*Jcd,2*Lambda,jd,jb,ja)
    if(b==c) r = r - OneBody%m(a,d) * (-1.d0)**( (ja+jb+jc+jd)/2 )*sjs(2*Jab,2*Jcd,2*Lambda,jd,ja,jb)
    if(a==d) r = r - OneBody%m(b,c) * (-1.d0)**( Jcd-Jab )       * sjs(2*Jab,2*Jcd,2*Lambda,jc,jb,ja)
    r = r * sqrt( dble( (2*Jab+1)*(2*Jcd+1) )) * (-1.d0)**Lambda
    if(a==b) r = r / sqrt(2.d0)
    if(c==d) r = r / sqrt(2.d0)
  end function embedding_onebody_in_twobody
end module TwoBodyLabOps
