module TwoBodyRelOps
  use LinAlgLib
  use omp_lib
  use ClassSys
  use TwoBodyRelativeSpace
  use OperatorDefinitions

  implicit none
  public :: TwoBodyRelOpChan
  public :: TwoBodyRelOp
  private :: InitTwoBodyRelOp
  private :: InitTwoBodyRelOpFromString
  private :: FinTwoBodyRelOp
  private :: CopyTwoBodyRelOp
  private :: SumTwoBodyRelOp
  private :: SubtractTwoBodyRelOp
  private :: ScaleTwoBodyRelOp
  private :: SetTwoBodyScalarSpinBreaking
  private :: SetTwoBodyScalarChanSpinBreaking
  private :: SetTwoBodyRelOp
  private :: SetNNSpecialOperators
  private :: SetNNChannelSpecialOperators
  private :: UnitaryTransformationTensor
  private :: SetTwoBodyRelOpChan
  private :: PrintTwoBodyRelOp
  private :: PrintTwoBodyRelOpChan
  private :: InitTwoBodyRelOpChan
  private :: FinTwoBodyRelOpChan
  private :: OperatorEvolution
  private :: SpinTensorDecomposition
  private :: GetTwoBodyRelOpME
  private :: SetTwoBodyRelOpME
  private :: TransformFromMomToHO
  private :: TransformFromMomToHOChan
  private :: get_mom_ho_overlap
  private :: ReadFile
  private :: WriteFile
  private :: ReadUTFile
  private :: WriteUTFile
  private :: NonReducedToReducedChan
  private :: ReducedToNonReducedChan
  private :: NonReducedToReduced
  private :: ReducedToNonReduced
  private :: TruncateModelSpace
  private :: SingularValueDecomposition

  type, extends(DMat) :: TwoBodyRelOpChan
    type(TwoBodyRelChanHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(str) :: OpName
    logical :: is_zero = .true.
    logical :: pn_same_mass = .true.
  contains
    procedure :: InitTwoBodyRelOpChan
    procedure :: FinTwoBodyRelOpChan
    procedure :: SetTwoBodyScalarChanSpinBreaking
    procedure :: SetTwoBodyRelOpChan
    procedure :: PrintTwoBodyRelOpChan
    procedure :: SetNNChannelSpecialOperators
    procedure :: TransformFromMomToHOChan
    procedure :: NonReducedToReducedChan
    procedure :: ReducedToNonReducedChan

    generic :: init => InitTwoBodyRelOpChan
    generic :: release => FinTwoBodyRelOpChan
    generic :: prtm => PrintTwoBodyRelOpChan
    generic :: set => SetTwoBodyRelOpChan
  end type TwoBodyRelOpChan

  type, extends(OperatorDef) :: TwoBodyRelOp
    type(TwoBodyRelOpChan), allocatable :: MatCh(:,:)
    type(TwoBodyRelSpaceHOBasis), pointer :: rel_sp_bra=>null(), rel_sp_ket=>null()
    logical :: is_init = .false.
    logical :: pn_same_mass = .true.
  contains
    procedure :: InitTwoBodyRelOp
    procedure :: InitTwoBodyRelOpFromString
    procedure :: FinTwoBodyRelOp
    procedure :: CopyTwoBodyRelOp
    procedure :: SumTwoBodyRelOp
    procedure :: SubtractTwoBodyRelOp
    procedure :: ScaleTwoBodyRelOp
    procedure :: SetTwoBodyScalarSpinBreaking
    procedure :: UnitaryTransformationTensor
    procedure :: SetTwoBodyRelOp
    procedure :: SetNNSpecialOperators
    procedure :: PrintTwoBodyRelOp
    procedure :: TruncateModelSpace
    procedure :: OperatorEvolution
    procedure :: SpinTensorDecomposition
    procedure :: SetTwoBodyRelOpME
    procedure :: GetTwoBodyRelOpME
    procedure :: TransformFromMomToHO
    procedure :: ReadFile
    procedure :: WriteFile
    procedure :: ReadUTFile
    procedure :: WriteUTFile
    procedure :: NonReducedToReduced
    procedure :: ReducedToNonReduced
    procedure :: SingularValueDecomposition

    generic :: assignment(=) => CopyTwoBodyRelOp
    generic :: operator(+) => SumTwoBodyRelOp
    generic :: operator(-) => SubtractTwoBodyRelOp
    generic :: operator(*) => ScaleTwoBodyRelOp

    generic :: init => InitTwoBodyRelOp, InitTwoBodyRelOpFromString
    generic :: fin => FinTwoBodyRelOp
    generic :: prt => PrintTwoBodyRelOp
    generic :: set => SetTwoBodyRelOp
    generic :: UT => UnitaryTransformationTensor
    generic :: truncate => TruncateModelSpace
    generic :: evolve => OperatorEvolution
    generic :: Get2ME => GetTwoBodyRelOpME
    generic :: SVD => SingularValueDecomposition
  end type TwoBodyRelOp

contains
  subroutine InitTwoBodyRelOpChan(this,OpName,ch_bra,ch_ket,pn_same_mass)
    class(TwoBodyRelOpChan), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelChanHOBasis), target, intent(in) :: ch_bra, ch_ket
    logical, intent(in) :: pn_same_mass
    this%is_zero = .false.
    this%rel_ch_bra => ch_bra
    this%rel_ch_ket => ch_ket
    this%OpName = OpName
    this%pn_same_mass = pn_same_mass
    if( OpName%val == "UT" .or. OpName%val == "UnitaryTransformation" ) then
      call this%eye( ch_ket%GetNumberStates() )
    else
      call this%zeros( ch_bra%GetNumberStates() ,ch_ket%GetNumberStates() )
    end if
  end subroutine InitTwoBodyRelOpChan

  subroutine FinTwoBodyRelOpChan(this)
    class(TwoBodyRelOpChan), intent(inout) :: this
    if(.not. allocated(this%m)) return
    this%is_zero = .true.
    this%pn_same_mass = .true.
    this%rel_ch_bra => null()
    this%rel_ch_ket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelOpChan

  subroutine InitTwoBodyRelOpFromString(this, rel_sp_bra, rel_sp_ket, oprtr, pn_same_mass)
    class(TwoBodyRelOp), intent(inout) :: this
    type(TwoBodyRelSpaceHOBasis), intent(in) :: rel_sp_bra, rel_sp_ket
    type(str), intent(in) :: oprtr
    logical, intent(in) :: pn_same_mass

    if(allocated(this%MatCh)) call this%fin()
    call this%InitOpDef(oprtr, .true. )
    call this%init(rel_sp_bra,rel_sp_ket,this%GetOpJ(),this%GetOpP(),this%GetOpZ(), pn_same_mass)
  end subroutine InitTwoBodyRelOpFromString

  subroutine InitTwoBodyRelOp(this, rel_sp_bra, rel_sp_ket, jr, pr, zr, pn_same_mass)
    use MyLibrary, only: triag
    class(TwoBodyRelOp), intent(inout) :: this
    type(TwoBodyRelSpaceHOBasis), target, intent(in) :: rel_sp_bra, rel_sp_ket
    type(TwoBodyRelChanHOBasis), pointer :: chbra, chket
    integer, intent(in) :: jr, pr, zr
    logical, intent(in) :: pn_same_mass
    integer :: ichb, ichk, nb, nk
    integer :: jb, pb, tzb
    integer :: jk, pk, tzk
    if(allocated(this%MatCh)) call this%fin()
    this%rel_sp_bra => rel_sp_bra
    this%rel_sp_ket => rel_sp_ket
    this%pn_same_mass = pn_same_mass
    call this%InitOpDef(.true., jr, pr, zr )
    allocate(this%MatCh(rel_sp_bra%GetNumberChannels(), rel_sp_ket%GetNumberChannels() ))
    do ichb = 1, this%rel_sp_bra%GetNumberChannels()
      chbra => this%rel_sp_bra%GetChannel(ichb)
      nb = chbra%GetNumberStates()
      jb = chbra%GetJ()
      pb = chbra%GetParity()
      tzb= chbra%GetZ()

      do ichk = 1, this%rel_sp_ket%GetNumberChannels()
        chket => this%rel_sp_ket%GetChannel(ichk)
        nk = chket%GetNumberStates()
        jk = chket%GetJ()
        pk = chket%GetParity()
        tzk= chket%GetZ()
        if(triag(jb, jk, jr)) cycle
        if(pb * pr /= pk) cycle
        if(abs(tzk-tzb)/=zr) cycle
        call this%MatCh(ichb,ichk)%init(this%GetOpName(),chbra,chket,pn_same_mass)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelOp

  subroutine FinTwoBodyRelOp(this)
    class(TwoBodyRelOp), intent(inout) :: this
    integer :: ichb, ichk
    if(.not. this%is_init) return
    do ichb = 1, size(this%MatCh, 1)
      do ichk = 1, size(this%MatCh, 2)
        call this%MatCh(ichb, ichk)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%rel_sp_bra => null()
    this%rel_sp_ket => null()
    call this%FinOperatorDef()
    this%is_init = .false.
  end subroutine FinTwoBodyRelOp

  subroutine CopyTwoBodyRelOp(a, b)
    class(TwoBodyRelOp), intent(inout) :: a
    class(TwoBodyRelOp), intent(in) :: b
    integer :: ichb, ichk, n
    if(allocated(a%MatCh)) call a%fin()
    n = size(b%MatCh, 1)
    call a%CopyOperatorDef(b%OperatorDef)
    a%rel_sp_bra => b%rel_sp_bra
    a%rel_sp_ket => b%rel_sp_ket
    a%pn_same_mass = b%pn_same_mass
    allocate(a%MatCh(n,n))
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%is_Zero) then
          a%MatCh(ichb,ichk)%is_zero = .true.
          cycle
        end if
        a%MatCh(ichb,ichk)%is_Zero = .false.
        a%MatCh(ichb,ichk)%pn_same_mass = b%MatCh(ichb, ichk)%pn_same_mass
        a%MatCh(ichb,ichk)%DMat = b%MatCh(ichb, ichk)%DMat
        a%MatCh(ichb,ichk)%rel_ch_bra => b%MatCh(ichb,ichk)%rel_ch_bra
        a%MatCh(ichb,ichk)%rel_ch_ket => b%MatCh(ichb,ichk)%rel_ch_ket
        a%MatCh(ichb,ichk)%OpName = b%MatCh(ichb,ichk)%OpName
      end do
    end do
    a%is_init = b%is_init
  end subroutine CopyTwoBodyRelOp

  function SumTwoBodyRelOp(a, b) result(c)
    use LinAlgLib
    type(TwoBodyRelOp) :: c
    class(TwoBodyRelOp), intent(in) :: a, b
    integer :: ichb, ichk, n
    if(a%pn_same_mass .neqv. b%pn_same_mass) then
      write(*,*) 'Error:', __LINE__, __FILE__
      stop
    end if
    c = a
    n = size(a%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%is_Zero) then
          c%MatCh(ichb,ichk)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichb,ichk)%is_Zero = .false.
        c%MatCh(ichb, ichk)%DMat = a%MatCh(ichb, ichk)%DMat + b%MatCh(ichb, ichk)%DMat
      end do
    end do
  end function SumTwoBodyRelOp

  function SubtractTwoBodyRelOp(a, b) result(c)
    type(TwoBodyRelOp) :: c
    class(TwoBodyRelOp), intent(in) :: a, b
    integer :: ichb, ichk, n
    if(a%pn_same_mass .neqv. b%pn_same_mass) then
      write(*,*) 'Error:', __LINE__, __FILE__
      stop
    end if
    c = a
    n = size(a%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%is_Zero) then
          c%MatCh(ichb,ichk)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichb, ichk)%is_Zero = .false.
        c%MatCh(ichb, ichk)%DMat = a%MatCh(ichb, ichk)%DMat - b%MatCh(ichb, ichk)%DMat
      end do
    end do
  end function SubtractTwoBodyRelOp

  function ScaleTwoBodyRelOp(a, b) result(c)
    type(TwoBodyRelOp) :: c
    class(TwoBodyRelOp), intent(in) :: a
    real(8), intent(in) :: b
    integer :: ichb, ichk, n
    c = a
    n = size(a%MatCh, 1)
    do ichb = 1, n
      do ichk = 1, n
        if(a%MatCh(ichb,ichk)%is_Zero) then
          c%MatCh(ichb,ichk)%is_Zero = .true.
          cycle
        end if
        c%MatCh(ichb, ichk)%is_Zero = .false.
        c%MatCh(ichb, ichk)%DMat = a%MatCh(ichb, ichk)%DMat * b
      end do
    end do
  end function ScaleTwoBodyRelOp

  subroutine SetTwoBodyScalarSpinBreaking(this, vho)
    use NNForce, only: NNForceHO
    class(TwoBodyRelOp), intent(inout) :: this
    type(NNForceHO), intent(in) :: vho
    type(TwoBodyRelChanHOBasis), pointer :: tbc
    integer :: ich
    type(str) :: OpName

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOp before calling SetTwoBodyScalarSpinBreaking"
      return
    end if

    if(.not. vho%is_init) then
      write(*,*) "Initialize NNForceHO before calling SetTwoBodyScalarSpinBreaking"
      return
    end if

    if(this%GetOpJ() /= 0 .or. this%GetOpP() /= 1 .or. this%GetOpZ() /= 0) then
      OpName = this%GetOpName()
      write(*,*) this%GetOpJ(), this%GetOpP(), this%GetOpZ(), OpName%val
      write(*,*) "Operator rank is wrong: NN interaction"
      return
    end if

    if(this%rel_sp_bra%GetNmax() /= vho%ms%GetNmax() ) then
      write(*,'(a)') 'In SetNNForce: Nmax of two and relspin does not match!'
      return
    end if

    if(this%rel_sp_ket%GetNmax() /= vho%ms%GetNmax() ) then
      write(*,'(a)') 'In SetNNForce: Nmax of two and relspin does not match!'
      return
    end if

    do ich = 1, this%rel_sp_ket%GetNumberChannels()
      tbc => this%rel_sp_ket%GetChannel(ich)
      call this%MatCh(ich,ich)%SetTwoBodyScalarChanSpinBreaking(vho, tbc%GetJ(), tbc%GetParity(), tbc%GetZ())
    end do
  end subroutine SetTwoBodyScalarSpinBreaking

  subroutine SetTwoBodyScalarChanSpinBreaking(this, vho, j, p, tz)
    use NNForce, only: NNForceHO
    class(TwoBodyRelOpChan), intent(inout) :: this
    type(NNForceHO), intent(in) :: vho
    type(TwoBodyRelChanHOBasis), pointer :: two
    type(TwoBodyRelSpaceSpinHOBasis), pointer :: relspin
    type(TwoBodyRelChanSpinHOBasis), pointer :: spin_bra, spin_ket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer, intent(in) :: j, p, tz
    integer :: n
    integer :: bra, ket
    integer :: n1, l1, s1, n2, l2, s2
    integer :: chb, chk, ib, ik

    if( j > vho%ms%GetJmax() ) return
    two => this%rel_ch_ket
    relspin => vho%ms
    n = two%GetNumberStates()
    call this%zeros(n,n)
    do bra = 1, n
      ho_bra => two%getp(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      s1 = ho_bra%GetS()
      chb = relspin%GetIndex(j,p,s1,tz)
      if(chb == 0) cycle
      spin_bra => relspin%GetChannel(chb)
      ib = spin_bra%GetIndex(n1,l1)
      do ket = 1, n
        ho_ket => two%getp(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        s2 = ho_ket%GetS()
        chk = relspin%GetIndex(j,p,s2,tz)
        if(chk == 0) cycle
        spin_ket => relspin%GetChannel(chk)
        ik = spin_ket%GetIndex(n2,l2)
        if(chb /= chk) cycle
        this%m(bra,ket) = vho%MatCh(chb)%m(ib,ik)
      end do
    end do
    this%is_Zero = .false.

  end subroutine SetTwoBodyScalarChanSpinBreaking

  function SpinTensorDecomposition(this, rank) result(op)
    use MyLibrary, only: triag, sjs
    class(TwoBodyRelOp), intent(in) :: this
    integer, intent(in) :: rank
    type(TwoBodyRelOp) :: op
    type(TwoBodyRelChanHOBasis), pointer :: tbc, tbc_j
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: bra, ket, nbra, nket, lbra, lket, sbra, sket
    integer :: j, p, z, jj
    integer :: ich, ich_j, bra_j, ket_j, Jmax
    real(8) :: me
    type(str) :: opname

    op = this
    if( this%GetOpJ() /= 0 ) then
      write(*,*) "SpinTensorDecompsition is impolemented only for scalar operator"
      write(*,*) __LINE__, __FILE__
      opname = this%GetOpName()
      write(*,*) "OpName=", trim(opname%val), "J=", this%GetOpJ(), ", Parity=", this%GetOpP(), ", Z=", this%GetOpZ()
      return
    end if
    if( rank > 2) then
      write(*,*) "Two-body spin-tensor decomposition is valid only up to rank=2"
      return
    end if

    Jmax = op%rel_sp_ket%GetJmax()
    do ich = 1, op%rel_sp_ket%GetNumberChannels()
      tbc => op%rel_sp_ket%GetChannel(ich)
      jj = tbc%GetJ()
      p  = tbc%GetParity()
      z  = tbc%GetZ()
      op%MatCh(ich,ich)%m = 0.d0
      do bra = 1, tbc%GetNumberStates()
        ho_bra => tbc%getp(bra)
        nbra = ho_bra%GetN()
        lbra = ho_bra%GetL()
        sbra = ho_bra%GetS()
        do ket = 1, tbc%GetNumberStates()
          ho_ket => tbc%getp(ket)
          nket = ho_ket%GetN()
          lket = ho_ket%GetL()
          sket = ho_ket%GetS()
          if(triag(lbra,lket,rank)) cycle
          if(triag(sbra,sket,rank)) cycle
          me = 0.d0
          do j = max(abs(lbra-sbra), abs(lket-sket)), min(lbra+sbra,lket+sket,Jmax)
            ich_j = this%rel_sp_ket%GetIndex(j,p,z)
            if(ich_j == 0) cycle
            tbc_j => this%rel_sp_ket%GetChannel(ich_j)
            bra_j = tbc_j%GetIndex(nbra,lbra,sbra)
            ket_j = tbc_j%GetIndex(nket,lket,sket)
            if(bra_j * ket_j == 0) cycle
            me = me + dble( 2*j+1 ) * &
                & sjs(2*lbra, 2*sbra, 2*j, 2*sket, 2*lket, 2*rank) * &
                & (-1.d0)**j * this%MatCh(ich_j,ich_j)%m(bra_j,ket_j)
          end do
          op%MatCh(ich,ich)%m(bra,ket) = me * (-1.d0)**jj * dble(2*rank+1) * &
                & sjs(2*lbra, 2*sbra, 2*jj, 2*sket, 2*lket, 2*rank)
        end do
      end do
    end do
  end function SpinTensorDecomposition

  subroutine UnitaryTransformationTensor(this, Ubra, Uket)
    class(TwoBodyRelOp), intent(inout) :: this
    type(TwoBodyRelOp), intent(in) :: Ubra, Uket
    type(DMat) :: tmp
    integer :: ichb, ichk, nbra, nket
    integer :: ichb_U, ichk_U
    integer :: jb, pb, zb, jk, pk, zk, bra, ket
    type(TwoBodyRelSpaceHOBasis), pointer :: rel_bra, rel_ket
    type(TwoBodyRelChanHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket

    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    do ichb = 1, rel_bra%GetNumberChannels()
      rel_ch_bra => rel_bra%GetChannel(ichb)
      jb = rel_ch_bra%GetJ()
      pb = rel_ch_bra%GetParity()
      zb = rel_ch_bra%GetZ()
      ichb_U = Ubra%rel_sp_bra%GetIndex(jb,pb,zb)
      do ichk = 1, rel_bra%GetNumberChannels()
        rel_ch_ket => rel_ket%GetChannel(ichk)
        jk = rel_ch_ket%GetJ()
        pk = rel_ch_ket%GetParity()
        zk = rel_ch_ket%GetZ()
        ichk_U = Uket%rel_sp_ket%GetIndex(jk,pk,zk)
        if(ichb_U * ichk_U==0) cycle
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        nbra = Ubra%MatCh(ichb_U,ichb_U)%rel_ch_ket%GetNumberStates()
        nket = Uket%MatCh(ichk_U,ichk_U)%rel_ch_ket%GetNumberStates()
        if(nbra*nket==0) cycle
        call tmp%ini(nbra,nket)
        !$omp parallel
        !$omp do private(bra, ho_bra, ket, ho_ket)
        do bra = 1, nbra
          ho_bra => Ubra%MatCh(ichb_U, ichb_U)%rel_ch_ket%GetP(bra)
          do ket = 1, nket
            ho_ket => Uket%MatCh(ichk_U, ichk_U)%rel_ch_ket%GetP(ket)
            tmp%m(bra,ket) = this%MatCh(ichb,ichk)%m(rel_ch_bra%GetIndex(ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS()), &
                &                                    rel_ch_ket%GetIndex(ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS()))
          end do
        end do
        !$omp end do
        !$omp end parallel
        tmp = Ubra%MatCh(ichb_U,ichb_U)%DMat%T() * tmp * Uket%MatCh(ichk_U,ichk_U)%DMat
        !$omp parallel
        !$omp do private(bra, ho_bra, ket, ho_ket)
        do bra = 1, nbra
          ho_bra => Ubra%MatCh(ichb_U, ichb_U)%rel_ch_ket%GetP(bra)
          do ket = 1, nket
            ho_ket => Uket%MatCh(ichk_U, ichk_U)%rel_ch_ket%GetP(ket)
            this%MatCh(ichb,ichk)%m(rel_ch_bra%GetIndex(ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS()), &
                &                   rel_ch_ket%GetIndex(ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS())) = tmp%m(bra,ket)
          end do
        end do
        !$omp end do
        !$omp end parallel
        call tmp%fin()
      end do
    end do
  end subroutine UnitaryTransformationTensor

  subroutine SetTwoBodyRelOp(this, NNint, pmax, NMesh)
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceHOBasis
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in), optional :: NNint
    integer, intent(in), optional :: NMesh
    real(8), intent(in), optional :: pmax
    type(TwoBodyRelChanHOBasis), pointer :: chbra, chket
    integer :: ichb, ichk
    integer :: bra(3), ket(3)
    type(str) :: opname
    type(sys) :: s

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOp before calling SetTwoBodyRelOp"
      return
    end if

    if(present( NNint ) .and. present( pmax ) .and. present( NMesh ) ) then
      opname = this%GetOpName()
      if( s%find(this%GetOpName(), s%str("PVTC")) .or. &
          & s%find(this%GetOpName(), s%str("PVTV")) .or. &
          & s%find(this%GetOpName(), s%str("AxialV")) .or. &
          & s%find(this%GetOpName(), s%str("0vbb")) .or. &
          & s%find(this%GetOpName(), s%str("0veeCap"))) then
        call this%SetNNSpecialOperators( NNint, pmax, NMesh )
        return
      end if
    end if

    do ichb = 1, this%rel_sp_bra%GetNumberChannels()
      chbra => this%rel_sp_bra%GetChannel(ichb)
      bra(1) = chbra%GetJ()
      bra(2) = chbra%GetParity()
      bra(3) = chbra%GetZ()
      do ichk = 1, this%rel_sp_ket%GetNumberChannels()
        chket => this%rel_sp_ket%GetChannel(ichk)
        ket(1) = chket%GetJ()
        ket(2) = chket%GetParity()
        ket(3) = chket%GetZ()
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        call this%MatCh(ichb, ichk)%set(bra, ket)
      end do
    end do
  end subroutine SetTwoBodyRelOp

  subroutine SetTwoBodyRelOpChan(this, bra, ket, NNint, pmax, NMesh)
    use Profiler, only: timer
    use TwoBodyRelativeSpace, only: TwoBodyRelChanHOBasis
    class(TwoBodyRelOpChan), intent(inout) :: this
    integer, intent(in) :: bra(3), ket(3)
    type(str), intent(in), optional :: NNint
    real(8), intent(in), optional :: pmax
    integer, intent(in), optional :: NMesh
    type(TwoBodyRelChanHOBasis), pointer :: ch_bra, ch_ket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: n1, n2
    real(8) :: ti
    type(sys) :: s
    type(str) :: OpName
    type(OperatorDef) :: op

    ch_bra => this%rel_ch_bra
    ch_ket => this%rel_ch_ket

    OpName = this%OpName
    if(present( NNint ) .and. present( pmax ) .and. present( NMesh ) ) then
      if( s%find(OpName, s%str("PVTC")) .or. &
          & s%find(OpName, s%str("PVTV")) .or. &
          & s%find(OpName, s%str("AxialV")) .or. &
          & s%find(OpName, s%str("0vbb")) .or. &
          & s%find(OpName, s%str("0veeCap"))) then
          call this%SetNNChannelSpecialOperators( NNint, pmax, NMesh )
          return
        end if
      end if

    ti = omp_get_wtime()
    call this%zeros( ch_bra%GetNumberStates(), ch_ket%GetNumberStates() )
    call op%InitOpDef(this%OpName, .true.)
    !$omp parallel
    !$omp do private(n1,ho_bra,n2,ho_ket)
    do n1 = 1, ch_bra%GetNumberStates()
      ho_bra => ch_bra%getp(n1)
      do n2 = 1, ch_ket%GetNumberStates()
        ho_ket => ch_ket%getp(n2)
        this%m(n1, n2) = CalcMERel( op, &
            & [ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), bra(1), bra(3)],&
            & [ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), ket(1), ket(3)])
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%add(s%str('SetTwoBodyRelOpChan'),omp_get_wtime()-ti)
  end subroutine SetTwoBodyRelOpChan

  subroutine SetNNSpecialOperators( this, NNint, pmax, NMesh )
    use MyLibrary, only: gauss_legendre
    use TwoBodyRelOpsMesh
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in) :: NNint
    integer, intent(in), optional :: NMesh
    real(8), intent(in), optional :: pmax
    type(TwoBodyRelSpaceMBasis) :: msMom
    type(TwoBodyRelOpMesh) :: OpMom
    real(8), allocatable :: p(:), w(:)
    real(8) :: pp
    integer :: NN
    pp = 8.d0
    NN = 100
    if( present(pmax) ) pp = pmax
    if( present(NMesh)) NN = NMesh
    call msMom%init(0.d0, pp, NN, this%rel_sp_ket%GetJmax())
    call gauss_legendre(0.d0, pp, p, w, NN)
    call msMom%setp(p,w)
    call OpMom%init( msMom, this%GetOpName() )
    call OpMom%set( NNint )
    call this%TransformFromMomToHO( OpMom )
    call OpMom%fin()
    call msMom%fin()
    deallocate(p)
    deallocate(w)
  end subroutine SetNNSpecialOperators

  subroutine SetNNChannelSpecialOperators( this, NNint, pmax, NMesh )
    use MyLibrary, only: gauss_legendre
    use TwoBodyRelOpsMesh
    class(TwoBodyRelOpChan), intent(inout) :: this
    type(str), intent(in) :: NNint
    integer, intent(in) :: NMesh
    real(8), intent(in) :: pmax
    type(TwoBodyRelChanHOBasis), pointer :: bra, ket
    type(TwoBodyRelChanMBasis) :: mbra, mket
    type(TwoBodyRelOpMeshChan) :: OpMom
    real(8), allocatable :: p(:), w(:)

    bra => this%rel_ch_bra
    ket => this%rel_ch_ket
    call mbra%init( NMesh, bra%GetJ(), bra%GetParity(), bra%GetZ() )
    call mket%init( NMesh, ket%GetJ(), ket%GetParity(), ket%GetZ() )
    call gauss_legendre(0.d0, pmax, p, w, NMesh)
    call mbra%setp(p, w)
    call mket%setp(p, w)
    call OpMom%init( this%OpName, mbra, mket )
    call OpMom%SetTwoBodyRelOpMeshChannel(NNint )
    call this%TransformFromMomToHOChan( OpMom )
    call OpMom%fin()
    call mbra%fin()
    call mket%fin()
    deallocate(p,w)
  end subroutine SetNNChannelSpecialOperators

  function TruncateModelSpace(this, new_bra, new_ket) result(op)
    class(TwoBodyRelOp), intent(in) :: this
    type(TwoBodyRelSpaceHOBasis), intent(in), target :: new_bra, new_ket
    type(twoBodyRelSpaceHOBasis), pointer :: old_bra, old_ket
    type(TwoBodyRelChanHOBasis), pointer :: ch_bra, ch_ket, ch_bra_old, ch_ket_old
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(TwoBodyRelOp) :: op
    integer :: chbra, jbra, pbra, zbra, bra, nbra, lbra, sbra
    integer :: chket, jket, pket, zket, ket, nket, lket, sket
    integer :: chbra_old, chket_old, bra_old, ket_old

    call op%init(new_bra, new_ket, this%GetOpName(), op%pn_same_mass )
    old_bra => this%rel_sp_bra
    old_ket => this%rel_sp_ket
    do chbra = 1, new_bra%GetNumberChannels()
      ch_bra => new_bra%GetChannel(chbra)
      jbra = ch_bra%GetJ()
      pbra = ch_bra%GetParity()
      zbra = ch_bra%GetZ()
      do chket = 1, new_ket%GetNumberChannels()
        ch_ket => new_ket%GetChannel(chket)
        jket = ch_ket%GetJ()
        pket = ch_ket%GetParity()
        zket = ch_ket%GetZ()
        if(op%MatCh(chbra,chket)%is_zero) cycle
        chbra_old = old_bra%GetIndex(jbra,pbra,zbra)
        chket_old = old_ket%GetIndex(jket,pket,zket)
        if(chbra_old * chket_old == 0) cycle
        ch_bra_old => old_bra%GetChannel(chbra_old)
        ch_ket_old => old_ket%GetChannel(chket_old)

        do bra = 1, ch_bra%GetNumberStates()
          ho_bra => ch_bra%getp(bra)
          nbra = ho_bra%GetN()
          lbra = ho_bra%GetL()
          sbra = ho_bra%GetS()
          if(2*nbra + lbra > ch_bra_old%GetNmax() ) cycle
          do ket = 1, ch_ket%GetNumberStates()
            ho_ket => ch_ket%getp(ket)
            nket = ho_ket%GetN()
            lket = ho_ket%GetL()
            sket = ho_ket%GetS()
            if(2*nket + lket > ch_ket_old%GetNmax() ) cycle
            bra_old = ch_bra_old%GetIndex(nbra,lbra,sbra)
            ket_old = ch_ket_old%GetIndex(nket,lket,sket)
            op%MatCh(chbra,chket)%m(bra,ket) = &
              & this%MatCh(chbra_old, chket_old)%m(bra_old, ket_old)
          end do
        end do

      end do
    end do
  end function TruncateModelSpace

  subroutine SetSymmetryViolatingForce(this)
    class(TwoBodyRelOp), intent(inout) :: this

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOp before calling SetTwoBodyRelOp"
      return
    end if
  end subroutine SetSymmetryViolatingForce

  subroutine SingularValueDecomposition(this, svd_rank)
    class(TwoBodyRelOp), intent(inout) :: this
    integer, intent(in) :: svd_rank
    integer :: ichb, ichk
    type(DMat) :: A, tmp
    type(DVec) :: Sigma
    type(DSingularValueDecomposition) :: sol
    if(svd_rank==-1) return
    do ichb = 1, this%rel_sp_bra%GetNumberChannels()
      do ichk = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        A = this%MatCh(ichb,ichk)%DMat
        call sol%init(A)
        call sol%svd(A)
        call Sigma%ini(sol%Sigma%n_size)
        Sigma%v(:svd_rank) = sol%Sigma%v(:svd_rank)
        call tmp%DiagMat(Sigma)
        this%MatCh(ichb,ichk)%DMat = sol%U * tmp * sol%V
        call sol%fin()
      end do
    end do
  end subroutine SingularValueDecomposition

  subroutine OperatorEvolution(this, params, Nmax_evolution_bra, Nmax_evolution_ket)
    use NuHamilInput, only: InputParameters
    use NNForce, only: NNForceHO
    class(TwoBodyRelOp), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    integer, intent(in), optional :: Nmax_evolution_bra, Nmax_evolution_ket
    integer :: Nmax_bra, Nmax_ket
    type(TwoBodyRelSpaceSpinHOBasis) :: relspin
    type(NNForceHO) :: vnnspin, Uspin
    type(TwoBodyRelSpaceHOBasis) :: rel_bra, rel_ket
    type(TwoBodyRelOp) :: Ubra, Uket
    type(sys) :: s

    Nmax_bra = params%N2max
    Nmax_ket = params%N2max
    if(present(Nmax_evolution_bra)) Nmax_bra=Nmax_evolution_bra
    if(present(Nmax_evolution_ket)) Nmax_ket=Nmax_evolution_ket

    call relspin%init(params%hw, Nmax_bra, params%Jmax2)
    call vnnspin%init(relspin)
    call Uspin%init(relspin)
    call vnnspin%setNNForceHO(Uspin, params)

    call rel_bra%init(params%hw, Nmax_bra, params%Jmax2)
    call Ubra%init(rel_bra,rel_bra,s%str('UT'),params%pn_same_mass)
    call Ubra%SetTwoBodyScalarSpinBreaking(Uspin)

    call vnnspin%fin()
    call Uspin%fin()
    call relspin%fin()

    if( Nmax_bra == Nmax_ket ) then
      rel_ket = rel_bra
      Uket = Ubra
    end if

    if( Nmax_bra /= Nmax_ket ) then
      call relspin%init(params%hw, Nmax_ket, params%Jmax2)
      call vnnspin%init(relspin)
      call Uspin%init(relspin)
      call vnnspin%setNNForceHO(Uspin, params)

      call rel_ket%init(params%hw, Nmax_ket, params%Jmax2)
      call Uket%init(rel_ket,rel_ket,s%str('UT'),params%pn_same_mass)
      call Uket%SetTwoBodyScalarSpinBreaking(Uspin)

      call vnnspin%fin()
      call Uspin%fin()
      call relspin%fin()
    end if

    call this%UT(Ubra,Uket)
    call Ubra%fin()
    call Uket%fin()
    call rel_bra%fin()
    call rel_ket%fin()
  end subroutine OperatorEvolution

  subroutine PrintTwoBodyRelOp(this,iunit)
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceHOBasis
    class(TwoBodyRelOp), intent(in) :: this
    integer, optional, intent(in) :: iunit
    integer :: unt = 6
    integer :: chbra, chket

    if(present(iunit)) unt = iunit

    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, chbra
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%prtm(unt)
      end do
    end do
  end subroutine PrintTwoBodyRelOp

  subroutine PrintTwoBodyRelOpChan(this,iunit)
    use TwoBodyRelativeSpace, only: TwoBodyRelChanHOBasis
    class(TwoBodyRelOpChan), intent(in), target :: this
    type(TwoBodyRelChanHOBasis), pointer :: chb, chk
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer, intent(in) :: iunit
    integer :: bra, ket, Nbra, Nket
    chb => this%rel_ch_bra
    chk => this%rel_ch_ket
    write(iunit,"(a)") "# n', l', s', j',Tz',  n,  l,  s,  j, Tz,          ME"
    do bra = 1, chb%GetNumberStates()
      ho_bra => chb%getp(bra)
      Nbra = 2*ho_bra%GetN() + ho_bra%GetL()
      do ket = 1, chk%GetNumberStates()
        ho_ket => chk%getp(ket)
        Nket = 2*ho_ket%GetN() + ho_ket%GetL()
        if(abs(this%m(bra,ket))<1.d-8) cycle
        write(iunit,'(10i4, es18.6)') &
            & ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), chb%GetJ(), chb%GetZ(), &
            & ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), chk%GetJ(), chk%GetZ(), this%m(bra,ket)
      end do
    end do
  end subroutine PrintTwoBodyRelOpChan

  subroutine NonReducedToReduced(this)
    class(TwoBodyRelOp), intent(inout) :: this
    integer :: chbra, chket

    if( this%GetOpJ() /= 0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
      return
    end if

    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%NonReducedToReducedChan()
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NonReducedToReduced

  subroutine ReducedToNonReduced(this)
    class(TwoBodyRelOp), intent(inout) :: this
    integer :: chbra, chket

    if( this%GetOpJ() /= 0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
      return
    end if

    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%ReducedToNonReducedChan()
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReducedChan(this)
    class(TwoBodyRelOpChan), intent(inout), target :: this
    type(TwoBodyRelChanHOBasis), pointer :: chbra, chket
    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    if( chbra%GetJ() /= chket%GetJ() ) then
      write(*,*) "Error, J' and J have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    this%m(:,:) = this%m(:,:) * sqrt( dble(2*chket%GetJ() + 1 ))
  end subroutine NonReducedToReducedChan

  subroutine ReducedToNonReducedChan(this)
    class(TwoBodyRelOpChan), intent(inout), target :: this
    type(TwoBodyRelChanHOBasis), pointer :: chbra, chket
    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    if( chbra%GetJ() /= chket%GetJ() ) then
      write(*,*) "Error, J' and J have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    this%m(:,:) = this%m(:,:) / sqrt( dble(2*chket%GetJ() + 1 ))
  end subroutine ReducedToNonReducedChan

  function GetTwoBodyRelOpME(this, bra, ket) result(r)
    !
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,z1]
    ! ket = [n2,l2,s2,j2,z2]
    !
    class(TwoBodyRelOp), intent(in) :: this
    integer, intent(in) :: bra(5), ket(5)
    real(8) :: r
    integer :: chbra, chket, ibra, iket
    type(TwoBodyRelSpaceHOBasis), pointer :: rel_bra, rel_ket

    r = 0.d0
    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    if(bra(2) > rel_bra%GetJmax()+1) return
    if(ket(2) > rel_ket%GetJmax()+1) return
    if(bra(4) > rel_bra%GetJmax()) return
    if(ket(4) > rel_ket%GetJmax()) return
    if(2*bra(1)+bra(2) > rel_bra%GetNmax()) return
    if(2*ket(1)+ket(2) > rel_ket%GetNmax()) return
    chbra = rel_bra%GetIndex(bra(4), (-1) ** bra(2), bra(5))
    chket = rel_ket%GetIndex(ket(4), (-1) ** ket(2), ket(5))
    if(chbra*chket == 0) return
    if(this%MatCh(chbra,chket)%is_zero) return
    ibra = rel_bra%jpz(chbra)%GetIndex(bra(1), bra(2), bra(3))
    iket = rel_ket%jpz(chket)%GetIndex(ket(1), ket(2), ket(3))
    if(ibra*iket == 0) return
    r = this%MatCh(chbra,chket)%m(ibra,iket)
  end function GetTwoBodyRelOpME

  subroutine SetTwoBodyRelOpME(this, bra, ket, me)
    !
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,z1]
    ! ket = [n2,l2,s2,j2,z2]
    !
    class(TwoBodyRelOp), intent(inout) :: this
    integer, intent(in) :: bra(5), ket(5)
    real(8), intent(in) :: me
    integer :: chbra, chket, ibra, iket
    type(TwoBodyRelSpaceHOBasis), pointer :: rel_bra, rel_ket

    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    if(bra(2) > rel_bra%GetJmax()+1) return
    if(ket(2) > rel_ket%GetJmax()+1) return
    if(bra(4) > rel_bra%GetJmax()) return
    if(ket(4) > rel_ket%GetJmax()) return
    if(2*bra(1)+bra(2) > rel_bra%GetNmax()) return
    if(2*ket(1)+ket(2) > rel_ket%GetNmax()) return
    chbra = rel_bra%GetIndex(bra(4), (-1) ** bra(2), bra(5))
    chket = rel_ket%GetIndex(ket(4), (-1) ** ket(2), ket(5))
    if(chbra*chket == 0) return
    if(this%MatCh(chbra,chket)%is_zero) return
    ibra = rel_bra%jpz(chbra)%GetIndex(bra(1), bra(2), bra(3))
    iket = rel_ket%jpz(chket)%GetIndex(ket(1), ket(2), ket(3))
    if(ibra*iket == 0) return
    this%MatCh(chbra,chket)%m(ibra,iket) = me
  end subroutine SetTwoBodyRelOpME

  subroutine TransformFromMomToHO(OpHO, OpMom)
    use TwoBodyRelOpsMesh
    class(TwoBodyRelOp), intent(inout), target :: OpHO
    type(TwoBodyRelOpMesh), intent(in), target :: OpMom
    type(TwoBodyRelSpaceHOBasis), pointer :: msHO_bra, msHO_ket
    type(TwoBodyRelSpaceMBasis), pointer :: msMom
    type(TwoBodyRelChanHOBasis), pointer :: HObra, HOket
    integer :: ichbra_ho, ichbra_mom, jbra, pbra, zbra
    integer :: ichket_ho, ichket_mom, jket, pket, zket

    msHO_bra => OpHO%rel_sp_bra
    msHO_ket => OpHO%rel_sp_ket
    msMom => OpMom%ms
    do ichbra_ho = 1, msHO_bra%GetNumberChannels()
      HObra => msHO_bra%GetChannel(ichbra_ho)
      jbra = HObra%GetJ()
      pbra = HObra%GetParity()
      zbra = HObra%GetZ()
      ichbra_mom = msMom%GetIndex(jbra,pbra,zbra)

      do ichket_ho = 1, msHO_ket%GetNumberChannels()
        if( OpHO%MatCh(ichbra_ho,ichket_ho)%is_zero ) cycle
        HOket => msHO_ket%GetChannel(ichket_ho)
        jket = HOket%GetJ()
        pket = HOket%GetParity()
        zket = HOket%GetZ()
        ichket_mom = msMom%GetIndex(jket,pket,zket)
        call OpHO%MatCh(ichbra_ho, ichket_ho)%TransformFromMomToHOChan( OpMom%MatCh(ichbra_mom,ichket_mom) )
      end do
    end do
  end subroutine TransformFromMomToHO

  subroutine TransformFromMomToHOChan( OpHO, OpMom )
    use TwoBodyRelOpsMesh
    class(TwoBodyRelOpChan), intent(inout), target :: OpHO
    type(TwoBodyRelOpMeshChan), intent(in), target :: OpMom
    type(DMat) :: ovlap_bra, ovlap_ket
    type(TwoBodyRelChanHOBasis), pointer :: HObra, HOket
    type(TwoBodyRelChanMBasis), pointer :: Mombra, Momket

    HObra => OpHO%rel_ch_bra
    HOket => OpHO%rel_ch_ket
    Mombra => OpMom%chbra
    Momket => OpMom%chket
    ovlap_bra = get_mom_ho_overlap( HObra, Mombra, OpHO%pn_same_mass )
    ovlap_ket = get_mom_ho_overlap( HOket, Momket, OpHO%pn_same_mass )
    OpHO%DMat = ovlap_bra%t() * OpMom%DMat * ovlap_ket
    !call OpHO%DMat%prt("")
  end subroutine TransformFromMomToHOChan

  function get_mom_ho_overlap( HO, Mom, pn_same_mass ) result(Mat)
    use MyLibrary, only: hc, m_proton, m_neutron, ho_radial_wf_norm
    type(TwoBodyRelChanHOBasis), intent(in) :: HO
    type(TwoBodyRelChanMBasis), intent(in) :: Mom
    logical, intent(in) :: pn_same_mass
    type(DMat) :: Mat
    integer :: iho, n, l, imom, s
    type(HarmonicOscillator), pointer :: ho_set
    type(Mesh), pointer :: mesh_set
    real(8) :: p, w, a, m_red
    m_red = 0.d0
    if(pn_same_mass) then
      m_red = (m_proton*m_neutron)/(m_proton+m_neutron)
    else
      if(HO%GetZ()==-1) m_red = m_proton * 0.5d0
      if(HO%GetZ()== 0) m_red = (m_proton*m_neutron)/(m_proton+m_neutron)
      if(HO%GetZ()== 1) m_red = m_neutron * 0.5d0
    end if
    a = hc ** 2 / (m_red * HO%GetFrequency())
    call Mat%zeros( Mom%GetNumberStates(), HO%GetNumberStates() )
    !$omp parallel
    !$omp do private(imom, mesh_set, p, w, l, s, iho, ho_set, n)
    do imom = 1, Mom%GetNumberStates()
      mesh_set => Mom%GetP(imom)
      p  = mesh_set%GetP()
      w  = mesh_set%GetW()
      l  = mesh_set%GetL()
      s  = mesh_set%GetS()
      do iho = 1, HO%GetNumberStates()
        ho_set => HO%GetP(iho)
        n = ho_set%GetN()
        if( l /= ho_set%GetL() ) cycle
        if( s /= ho_set%GetS() ) cycle
        Mat%m(imom,iho) = p * w * ho_radial_wf_norm(n,l,a,p) * (-1.d0)**n
      end do
    end do
    !$omp end do
    !$omp end parallel
  end function get_mom_ho_overlap

  subroutine ReadFile( this, filename )
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in) :: filename
    integer :: runit=20, io
    integer :: nbra, lbra, sbra, jbra, zbra
    integer :: nket, lket, sket, jket, zket
    real(8) :: v

    open(runit, file=filename%val, action="read")
    read(runit,*)
    read(runit,*)
    do
      read(runit,*,iostat=io) nbra, lbra, sbra, jbra, zbra, nket, lket, sket, jket, zket, v
      if(io == 0) then
        call this%SetTwoBodyRelOpME( [nbra,lbra,sbra,jbra,zbra], [nket,lket,sket,jket,zket], v )
        call this%SetTwoBodyRelOpME( [nket,lket,sket,jket,zket], [nbra,lbra,sbra,jbra,zbra], v*(-1.d0)**(jbra-jket) )
      end if
      if(io < 0) exit
      if(io > 0) then
        write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      end if
    end do
    close(runit)
  end subroutine ReadFile

  subroutine WriteFile( this, filename )
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in) :: filename
    type(str) :: OpName
    integer :: wunit=20
    integer :: ichbra, ichket, ibra, iket, iketmax
    integer :: jbra, jket, zbra, zket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    character(255) :: header

    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    open(wunit, file=filename%val, action="write")
    write(wunit,'(a)') trim(header)
    write(wunit,'(a)') "#  n'  l'  s'  j' tz'   n   l  s   j  tz           ME "
    do ichbra = 1, this%rel_sp_bra%GetNumberCHannels()
      do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        jbra = this%MatCh(ichbra,ichket)%rel_ch_bra%GetJ()
        jket = this%MatCh(ichbra,ichket)%rel_ch_ket%GetJ()
        zbra = this%MatCh(ichbra,ichket)%rel_ch_bra%GetZ()
        zket = this%MatCh(ichbra,ichket)%rel_ch_ket%GetZ()
        do ibra = 1, this%MatCh(ichbra,ichket)%rel_ch_bra%GetNumberStates()
          iketmax = this%MatCh(ichbra,ichket)%rel_ch_ket%GetNumberStates()
          if(ichbra == ichket) iketmax = ibra
          do iket = 1, iketmax
            ho_bra => this%MatCh(ichbra,ichket)%rel_ch_bra%getp(ibra)
            ho_ket => this%MatCh(ichbra,ichket)%rel_ch_ket%getp(iket)
            if(abs(this%MatCh(ichbra,ichket)%m(ibra,iket)) < 1.d-16) cycle
            write(wunit,'(10i4,es18.8)') ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), jbra, zbra, &
                & ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), jket, zket, this%MatCh(ichbra,ichket)%m(ibra,iket)
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine WriteFile

  subroutine ReadUTFile( this, filename )
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in) :: filename
    integer :: runit=20, io
    integer :: nbra, lbra, sbra, jbra, zbra
    integer :: nket, lket, sket, jket, zket
    real(8) :: v

    open(runit, file=filename%val, action="read")
    read(runit,*)
    read(runit,*)
    do
      read(runit,*,iostat=io) nbra, lbra, sbra, jbra, zbra, nket, lket, sket, jket, zket, v
      if(io == 0) then
        call this%SetTwoBodyRelOpME( [nbra,lbra,sbra,jbra,zbra], [nket,lket,sket,jket,zket], v )
      end if
      if(io < 0) exit
      if(io > 0) then
        write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      end if
    end do
    close(runit)
  end subroutine ReadUTFile

  subroutine WriteUTFile( this, filename )
    class(TwoBodyRelOp), intent(inout) :: this
    type(str), intent(in) :: filename
    type(str) :: OpName
    integer :: wunit=20
    integer :: ich, ibra, iket
    integer :: j, z
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    character(255) :: header

    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    open(wunit, file=filename%val, action="write")
    write(wunit,'(a)') trim(header)
    write(wunit,'(a)') "#  n'  l'  s'  j' tz'   n   l  s   j  tz           ME "
    do ich = 1, this%rel_sp_ket%GetNumberChannels()
      if( this%MatCh(ich,ich)%is_zero ) cycle
      j = this%MatCh(ich,ich)%rel_ch_ket%GetJ()
      z = this%MatCh(ich,ich)%rel_ch_ket%GetZ()
      do ibra = 1, this%MatCh(ich,ich)%rel_ch_bra%GetNumberStates()
        do iket = 1, this%MatCh(ich,ich)%rel_ch_ket%GetNumberStates()
          ho_bra => this%MatCh(ich,ich)%rel_ch_bra%getp(ibra)
          ho_ket => this%MatCh(ich,ich)%rel_ch_ket%getp(iket)
          write(wunit,'(10i4,es18.8)') ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), j, z, &
              & ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), j, z, this%MatCh(ich,ich)%m(ibra,iket)
        end do
      end do
    end do
    close(wunit)
  end subroutine WriteUTFile
end module TwoBodyRelOps
