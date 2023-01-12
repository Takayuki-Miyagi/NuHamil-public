module TwoBodyRelOpsIso
  use LinAlgLib
  use omp_lib
  use ClassSys
  use OperatorDefinitions
  use TwoBodyRelativeSpace
  implicit none
  public :: TwoBodyRelOpChanIso
  public :: TwoBodyRelOpIso

  private :: FinTwoBodyRelOp
  private :: InitTwoBodyRelOpFromString
  private :: InitTwoBodyRelOp
  private :: CopyTwoBodyRelOp
  private :: SumTwoBodyRelOp
  private :: SubtractTwoBodyRelOp
  private :: ScaleTwoBodyRelOp
  private :: SetTwoBodyScalarSpinBreaking
  private :: SetTwoBodyScalarChanSpinBreaking
  private :: MakeIsospinSymmetric
  private :: MakeIsospinSymmetricChannel
  private :: UnitaryTransformationTensor
  private :: SetTwoBodyRelOp
  private :: SetTwoBodyRelOpChan
  private :: PrintTwoBodyRelOp
  private :: PrintTwoBodyRelOpChan
  private :: InitTwoBodyRelOpChan
  private :: FinTwoBodyRelOpChan
  private :: OperatorEvolution
  private :: SpinTensorDecomposition
  private :: GetTwoBodyRelOpIsoME
  private :: SetNNSpecialOperators
  private :: SetNNChannelSpecialOperators
  private :: TransformFromMomToHO
  private :: TransformFromMomToHOChan
  private :: NonReducedToReducedChan
  private :: ReducedToNonReducedChan
  private :: NonReducedToReduced
  private :: ReducedToNonReduced
  private :: TruncateModelSpace

  type, extends(DMat) :: TwoBodyRelOpChanIso
    type(TwoBodyRelChanIsoHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(str) :: OpName
    logical :: is_Zero = .true.
  contains
    procedure :: InitTwoBodyRelOpChan
    procedure :: FinTwoBodyRelOpChan
    procedure :: SetTwoBodyScalarChanSpinBreaking
    procedure :: MakeIsospinSymmetricChannel
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
  end type TwoBodyRelOpChanIso

  type, extends(OperatorDef) :: TwoBodyRelOpIso
    type(TwoBodyRelOpChanIso), allocatable :: MatCh(:,:)
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: rel_sp_bra, rel_sp_ket
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyRelOp
    procedure :: InitTwoBodyRelOpFromString
    procedure :: FinTwoBodyRelOp
    procedure :: CopyTwoBodyRelOp
    procedure :: SumTwoBodyRelOp
    procedure :: SubtractTwoBodyRelOp
    procedure :: ScaleTwoBodyRelOp
    procedure :: SetTwoBodyScalarSpinBreaking
    procedure :: MakeIsospinSymmetric
    procedure :: UnitaryTransformationTensor
    procedure :: SetTwoBodyRelOp
    procedure :: PrintTwoBodyRelOp
    procedure :: OperatorEvolution
    procedure :: SpinTensorDecomposition
    procedure :: GetTwoBodyRelOpIsoME
    procedure :: SetTwoBodyRelOpIsoME
    procedure :: SetNNSpecialOperators
    procedure :: TransformFromMomToHO
    procedure :: TruncateModelSpace
    procedure :: NonReducedToReduced
    procedure :: ReducedToNonReduced
    procedure :: ReadFile
    procedure :: WriteFile

    generic :: assignment(=) => CopyTwoBodyRelOp
    generic :: operator(+) => SumTwoBodyRelOp
    generic :: operator(-) => SubtractTwoBodyRelOp
    generic :: operator(*) => ScaleTwoBodyRelOp

    generic :: fin => FinTwoBodyRelOp
    generic :: init => InitTwoBodyRelOp, InitTwoBodyRelOpFromString
    generic :: prt => PrintTwoBodyRelOp
    generic :: UT => UnitaryTransformationTensor
    generic :: truncate => TruncateModelSpace
    generic :: set => SetTwoBodyRelOp
    generic :: evolve => OperatorEvolution
    generic :: Get2ME => GetTwoBodyRelOpIsoME
    generic :: Set2ME => SetTwoBodyRelOpIsoME
  end type TwoBodyRelOpIso

contains
  subroutine FinTwoBodyRelOp(this)
    class(TwoBodyRelOpIso), intent(inout) :: this
    integer :: ichb, ichk
    do ichb = 1, size(this%MatCh, 1)
      do ichk = 1, size(this%MatCh, 2)
        this%MatCh(ichb,ichk)%rel_ch_bra => null()
        this%MatCh(ichb,ichk)%rel_ch_ket => null()
        call this%MatCh(ichb, ichk)%fin()
      end do
    end do
    deallocate(this%MatCh)
    this%rel_sp_bra => null()
    this%rel_sp_ket => null()
    call this%FinOperatorDef()
    this%is_init = .false.
  end subroutine FinTwoBodyRelOp

  subroutine InitTwoBodyRelOpFromString(this, rel_sp_bra, rel_sp_ket, oprtr)
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(TwoBodyRelSpaceIsoHOBasis), intent(in) :: rel_sp_bra, rel_sp_ket
    type(str), intent(in) :: oprtr

    if(allocated(this%MatCh)) call this%fin()
    call this%InitOpDef(oprtr, .false. )
    call this%init(rel_sp_bra, rel_sp_ket, this%GetOpJ(), this%GetOpP(), this%GetOpT())
  end subroutine InitTwoBodyRelOpFromString

  subroutine InitTwoBodyRelOp(this, rel_sp_bra, rel_sp_ket, jr, pr, tr)
    use MyLibrary, only: triag
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(TwoBodyRelSpaceIsoHOBasis), target, intent(in) :: rel_sp_bra, rel_sp_ket
    type(TwoBodyRelChanIsoHOBasis), pointer :: chbra, chket
    integer, intent(in) :: jr, pr, tr
    integer :: ichb, ichk, nb, nk
    integer :: jb, pb, tb
    integer :: jk, pk, tk
    if(allocated(this%MatCh)) call this%fin()
    this%rel_sp_bra => rel_sp_bra
    this%rel_sp_ket => rel_sp_ket
    call this%InitOpDef(.false., jr, pr, tr )
    allocate(this%MatCh(rel_sp_bra%GetNumberChannels(), rel_sp_ket%GetNumberChannels() ))
    do ichb = 1, rel_sp_bra%GetNumberChannels()
      chbra => rel_sp_bra%GetChannel(ichb)
      nb = chbra%GetNumberStates()
      jb = chbra%GetJ()
      pb = chbra%GetParity()
      tb = chbra%GetT()
      do ichk = 1, rel_sp_ket%GetNumberChannels()
        chket => rel_sp_ket%GetChannel(ichk)
        nk = chket%GetNumberStates()
        jk = chket%GetJ()
        pk = chket%GetParity()
        tk = chket%GetT()
        if(triag(jb, jk, jr)) cycle
        if(pb * pr /= pk) cycle
        if(triag(tb, tk, tr)) cycle
        call this%MatCh(ichb,ichk)%init(this%GetOpName(),chbra,chket)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelOp

  subroutine InitTwoBodyRelOpChan(this,OpName,ch_bra,ch_ket)
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelChanIsoHOBasis), target, intent(in) :: ch_bra, ch_ket
    this%is_zero = .false.
    this%rel_ch_bra => ch_bra
    this%rel_ch_ket => ch_ket
    this%OpName = OpName
    if( OpName%val == "UT" .or. OpName%val == "UnitaryTransformation" ) then
      call this%eye( ch_ket%GetNumberStates() )
    else
      call this%zeros( ch_bra%GetNumberStates() ,ch_ket%GetNumberStates() )
    end if
  end subroutine InitTwoBodyRelOpChan

  subroutine FinTwoBodyRelOpChan(this)
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    if(.not. allocated(this%m)) return
    this%is_zero = .true.
    this%rel_ch_bra => null()
    this%rel_ch_ket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelOpChan

  subroutine CopyTwoBodyRelOp(a, b)
    class(TwoBodyRelOpIso), intent(inout) :: a
    class(TwoBodyRelOpIso), intent(in) :: b
    integer :: ichb, ichk, n
    if(allocated(a%MatCh)) call a%fin()
    call a%CopyOperatorDef(b%OperatorDef)
    n = size(b%MatCh, 1)
    a%rel_sp_bra => b%rel_sp_bra
    a%rel_sp_ket => b%rel_sp_ket
    allocate(a%MatCh(n,n))
    do ichb = 1, n
      do ichk = 1, n
        if(b%MatCh(ichb,ichk)%is_Zero) then
          a%MatCh(ichb,ichk)%is_zero = .true.
          cycle
        end if
        a%MatCh(ichb,ichk)%is_Zero = .false.
        a%MatCh(ichb,ichk)%DMat = b%MatCh(ichb, ichk)%DMat
        a%MatCh(ichb,ichk)%rel_ch_bra => b%MatCh(ichb,ichk)%rel_ch_bra
        a%MatCh(ichb,ichk)%rel_ch_ket => b%MatCh(ichb,ichk)%rel_ch_ket
      end do
    end do
    a%OperatorDef = b%OperatorDef
    a%is_init = b%is_init
  end subroutine CopyTwoBodyRelOp

  function SumTwoBodyRelOp(a, b) result(c)
    use LinAlgLib
    type(TwoBodyRelOpIso) :: c
    class(TwoBodyRelOpIso), intent(in) :: a, b
    integer :: ichb, ichk, n
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
    type(TwoBodyRelOpIso) :: c
    class(TwoBodyRelOpIso), intent(in) :: a, b
    integer :: ichb, ichk, n
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
    type(TwoBodyRelOpIso) :: c
    class(TwoBodyRelOpIso), intent(in) :: a
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
    use NNForceIsospin, only: NNForceHOIsospin
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(NNForceHOIsospin), intent(in) :: vho
    type(TwoBodyRelChanIsoHOBasis), pointer :: tbc
    integer :: ich

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOpIso before calling SetTwoBodyScalarSpinBreaking"
      return
    end if

    if(.not. vho%is_init) then
      write(*,*) "Initialize NNForceIsospin before calling SetTwoBodyScalarSpinBreaking"
      return
    end if

    if(this%GetOpJ() /= 0 .or. this%GetOpP() /= 1 .or. this%GetOpT() /= 0) then
      write(*,*) "Operator rank is wrong: NN interaction"
      return
    end if

    if(this%rel_sp_bra%GetNmax() /= vho%ms%GetNmax()) then
      write(*,'(a)') 'In SetNNForce: Nmax of two and relspin does not match!'
      return
    end if

    if(this%rel_sp_ket%GetNmax() /= vho%ms%GetNmax()) then
      write(*,'(a)') 'In SetNNForce: Nmax of two and relspin does not match!'
      return
    end if

    do ich = 1, this%rel_sp_ket%GetNumberChannels()
      tbc => this%rel_sp_ket%GetChannel(ich)
      call this%MatCh(ich,ich)%SetTwoBodyScalarChanSpinBreaking(vho, tbc%GetJ(), tbc%GetParity(), tbc%GetT())
    end do
  end subroutine SetTwoBodyScalarSpinBreaking

  subroutine SetTwoBodyScalarChanSpinBreaking(this, vho, j, p, t)
    use NNForceIsospin, only: NNForceHOIsospin
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    type(TwoBodyRelChanIsoHOBasis), pointer :: two
    type(TwoBodyRelSpaceSpinIsoHOBasis), pointer :: relspin
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(NNForceHOIsospin), intent(in) :: vho
    integer, intent(in) :: j, p, t
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
      chb = relspin%GetIndex(j,p,s1,t)
      ib = relspin%jpst(chb)%GetIndex(n1,l1)
      do ket = 1, n
        ho_ket => two%getp(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        s2 = ho_ket%GetS()
        chk = relspin%GetIndex(j,p,s2,t)
        ik = relspin%jpst(chk)%GetIndex(n2,l2)

        if(chb /= chk) cycle
        this%m(bra,ket) = vho%MatCh(chb)%m(ib,ik)
      end do
    end do
    this%is_Zero = .false.
  end subroutine SetTwoBodyScalarChanSpinBreaking

  subroutine MakeIsospinSymmetric(this, op, flip_projection)
    !
    ! Operator should have good isospin symmetry
    !
    use TwoBodyRelOps, only: TwoBodyRelOp
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(TwoBodyRelOp), intent(in) :: op
    integer, intent(in), optional :: flip_projection
    integer :: chbra, chket

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOpIso before calling ImposingIsospinTwoBodyRelOp"
      return
    end if

    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%MakeIsospinSymmetricChannel(op, this%GetOpT(), flip_projection )
      end do
    end do
  end subroutine MakeIsospinSymmetric

  subroutine MakeIsospinSymmetricChannel(this, op, Trank, flip_projection)
    !
    ! Operator should have good isospin symmetry
    !
    use MyLibrary, only: tjs, triag
    use TwoBodyRelOps, only: TwoBodyRelOp
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    type(TwoBodyRelOp), intent(in) :: op
    integer, intent(in) :: Trank
    integer, intent(in), optional :: flip_projection
    type(TwoBodyRelChanIsoHOBasis), pointer :: chbra, chket
    type(TwoBodyRelChanHOBasis), pointer :: chbra_op, chket_op
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: ichbra_op, ichket_op
    integer :: ibra, nbra, lbra, sbra, tbra
    integer :: iket, nket, lket, sket, tket
    integer :: idx_op_bra, idx_op_ket
    integer :: zbra, zket, z
    integer :: flip
    real(8) :: v, fact

    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    tbra = chbra%GetT()
    tket = chket%GetT()
    if(triag(tbra,tket,Trank)) return
    if(abs(op%GetOpZ()) > Trank) return
    flip = 1
    if(present(flip_projection)) flip = flip * flip_projection

    zbra = 0
    zket = 0
    do z = -tket, tket
      zket = z
      zbra = z + op%GetOpZ()
      if(abs(zbra) <= tbra) exit
      if(z==tket) return
    end do
    fact = (-1.d0)**(tbra-zbra) * tjs(2*tbra, 2*Trank, 2*tket, -2*zbra, 2*op%GetOpZ(), 2*zket)
    if(abs(fact) < 1.d-16) return
    fact = 1.d0 / fact
    ichbra_op = op%rel_sp_bra%GetIndex(chbra%GetJ(), chbra%GetParity(), zbra)
    ichket_op = op%rel_sp_ket%GetIndex(chket%GetJ(), chket%GetParity(), zket)
    if(ichbra_op * ichket_op ==0 ) return
    if(op%MatCh(ichbra_op,ichket_op)%is_zero) return
    chbra_op => op%MatCh(ichbra_op,ichket_op)%rel_ch_bra
    chket_op => op%MatCh(ichbra_op,ichket_op)%rel_ch_ket
    do ibra = 1, chbra%GetNumberStates()
      ho_bra => chbra%getp(ibra)
      nbra = ho_bra%GetN()
      lbra = ho_bra%GetL()
      sbra = ho_bra%GetS()
      do iket = 1, chket%GetNumberStates()
        ho_ket => chket%getp(iket)
        nket = ho_ket%GetN()
        lket = ho_ket%GetL()
        sket = ho_ket%GetS()
        idx_op_bra = chbra_op%GetIndex(nbra,lbra,sbra)
        idx_op_ket = chket_op%GetIndex(nket,lket,sket)
        v = op%MatCh(ichbra_op, ichket_op)%m(idx_op_bra, idx_op_ket)
        this%m(ibra,iket) = v * fact
      end do
    end do
  end subroutine MakeIsospinSymmetricChannel

  function SpinTensorDecomposition(this, rank) result(op)
    use MyLibrary, only: triag, sjs
    class(TwoBodyRelOpIso), intent(in) :: this
    integer, intent(in) :: rank
    type(TwoBodyRelOpIso) :: op
    type(TwoBodyRelChanIsoHOBasis), pointer :: tbc
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: bra, ket, nbra, nket, lbra, lket, sbra, sket
    integer :: j, p, t, jj
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
      jj= tbc%GetJ()
      p = tbc%GetParity()
      t = tbc%GetT()
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
            ich_j = this%rel_sp_ket%GetIndex(j,p,t)
            if(ich_j == 0) cycle
            bra_j = this%MatCh(ich_j,ich_j)%rel_ch_bra%GetIndex(nbra,lbra,sbra)
            ket_j = this%MatCh(ich_j,ich_j)%rel_ch_ket%GetIndex(nket,lket,sket)
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

  function TruncateModelSpace(this, new_bra, new_ket) result(op)
    class(TwoBodyRelOpIso), intent(in) :: this
    type(TwoBodyRelSpaceIsoHOBasis), intent(in), target :: new_bra, new_ket
    type(twoBodyRelSpaceIsoHOBasis), pointer :: old_bra, old_ket
    type(TwoBodyRelChanIsoHOBasis), pointer :: ch_bra, ch_ket, ch_bra_old, ch_ket_old
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(TwoBodyRelOpIso) :: op
    integer :: chbra, jbra, pbra, tbra, bra, nbra, lbra, sbra
    integer :: chket, jket, pket, tket, ket, nket, lket, sket
    integer :: chbra_old, chket_old, bra_old, ket_old

    call op%init(new_bra, new_ket, this%GetOpName() )
    old_bra => this%rel_sp_bra
    old_ket => this%rel_sp_ket
    do chbra = 1, new_bra%GetNumberChannels()
      ch_bra => new_bra%GetChannel(chbra)
      jbra = ch_bra%GetJ()
      pbra = ch_bra%GetParity()
      tbra = ch_bra%GetT()
      do chket = 1, new_ket%GetNumberChannels()
        ch_ket => new_ket%GetChannel(chket)
        jket = ch_ket%GetJ()
        pket = ch_ket%GetParity()
        tket = ch_ket%GetT()
        if(op%MatCh(chbra,chket)%is_zero) cycle
        chbra_old = old_bra%GetIndex(jbra,pbra,tbra)
        chket_old = old_ket%GetIndex(jket,pket,tket)
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

  subroutine UnitaryTransformationTensor(this, Ubra, Uket)
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: Ubra, Uket
    type(DMat) :: tmp
    integer :: ichb, ichk, nbra, nket
    integer :: ichb_U, ichk_U
    integer :: jb, pb, tb, jk, pk, tk, bra, ket
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: rel_bra, rel_ket
    type(TwoBodyRelChanIsoHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket

    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    do ichb = 1, rel_bra%GetNumberChannels()
      rel_ch_bra => rel_bra%GetChannel(ichb)
      jb = rel_ch_bra%GetJ()
      pb = rel_ch_bra%GetParity()
      tb = rel_ch_bra%GetT()
      ichb_U = Ubra%rel_sp_ket%GetIndex(jb,pb,tb)
      do ichk = 1, rel_ket%GetNumberChannels()
        rel_ch_ket => rel_ket%GetChannel(ichk)
        jk = rel_ch_ket%GetJ()
        pk = rel_ch_ket%GetParity()
        tk = rel_ch_ket%GetT()
        ichk_U = Uket%rel_sp_bra%GetIndex(jk,pk,tk)
        if(ichb_U * ichk_U==0) cycle
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        nbra = Ubra%MatCh(ichb_U,ichb_U)%rel_ch_ket%GetNumberStates()
        nket = Uket%MatCh(ichk_U,ichk_U)%rel_ch_bra%GetNumberStates()
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

  subroutine SetTwoBodyRelOp(this)
    use NuHamilInput, only: params
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceIsoHOBasis
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(str) :: opname
    type(TwoBodyRelChanIsoHOBasis), pointer :: chbra, chket
    integer :: ichb, ichk
    integer :: bra(3), ket(3)
    real(8) :: hw
    type(sys) :: s

    hw = this%rel_sp_bra%GetFrequency()

    if(.not. this%is_init) then
      write(*,*) "Initialize TwoBodyRelOpIso before calling SetTwoBodyRelOp"
      return
    end if

    opname = this%GetOpName()
    if( s%find(this%GetOpName(), s%str("PVTC")) .or. &
        & s%find(this%GetOpName(), s%str("PVTV")) .or. &
        & s%find(this%GetOpName(), s%str("AxialV")) .or. &
        & s%find(this%GetOpName(), s%str("0vbb")) .or. &
        & s%find(this%GetOpName(), s%str("0veeCap"))) then
      call this%SetNNSpecialOperators( params%NNint, params%pmax2, params%NMesh2 )
      return
    end if

    do ichb = 1, this%rel_sp_bra%GetNumberChannels()
      chbra => this%rel_sp_bra%GetChannel(ichb)
      bra(1) = chbra%GetJ()
      bra(2) = chbra%GetParity()
      bra(3) = chbra%GetT()
      do ichk = 1, this%rel_sp_ket%GetNumberChannels()
        chket => this%rel_sp_ket%GetChannel(ichk)
        ket(1) = chket%GetJ()
        ket(2) = chket%GetParity()
        ket(3) = chket%GetT()
        if(this%MatCh(ichb, ichk)%is_Zero) cycle
        call this%MatCh(ichb, ichk)%set(bra, ket)
      end do
    end do
  end subroutine SetTwoBodyRelOp

  subroutine SetTwoBodyRelOpChan(this, bra, ket)
    use Profiler, only: timer
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    type(TwoBodyRelChanIsoHOBasis), pointer :: chb, chk
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer, intent(in) :: bra(3), ket(3)
    integer :: n1, n2
    real(8) :: ti
    type(sys) :: s
    type(OperatorDef) :: op


    chb => this%rel_ch_bra
    chk => this%rel_ch_ket

    ti = omp_get_wtime()
    call this%zeros(chb%GetNumberStates(), chk%GetNumberStates())
    call op%InitOpDef(this%OpName, .false.)
    !$omp parallel
    !$omp do private(n1,ho_bra,n2,ho_ket)
    do n1 = 1, chb%GetNumberStates()
      ho_bra => chb%getp(n1)
      do n2 = 1, chk%GetNumberStates()
        ho_ket => chk%getp(n2)
        this%m(n1, n2) = CalcMERel( op, &
            & [ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), bra(1), bra(3)],&
            & [ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), ket(1), ket(3)])
      end do
    end do
    !$omp end do
    !$omp end parallel
    call timer%add(s%str('SetTwoBodyRelOpChannel'),omp_get_wtime()-ti)
  end subroutine SetTwoBodyRelOpChan

  subroutine SetNNSpecialOperators( this, NNint, pmax, NMesh )
    use MyLibrary, only: gauss_legendre
    use TwoBodyRelOpsMeshIso
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(str), intent(in), optional :: NNint
    integer, intent(in), optional :: NMesh
    real(8), intent(in), optional :: pmax
    type(TwoBodyRelSpaceIsoMBasis) :: msMom
    type(TwoBodyRelOpMeshIso) :: OpMom
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
    use TwoBodyRelOpsMeshIso
    class(TwoBodyRelOpChanIso), intent(inout) :: this
    type(str), intent(in), optional :: NNint
    integer, intent(in) :: NMesh
    real(8), intent(in) :: pmax
    type(TwoBodyRelChanIsoHOBasis), pointer :: bra, ket
    type(TwoBodyRelChanIsoMBasis) :: mbra, mket
    type(TwoBodyRelOpMeshChanIso) :: OpMom
    real(8), allocatable :: p(:), w(:)

    bra => this%rel_ch_bra
    ket => this%rel_ch_ket
    call mbra%init( NMesh, bra%GetJ(), bra%GetParity(), bra%GetT() )
    call mket%init( NMesh, ket%GetJ(), ket%GetParity(), ket%GetT() )
    call gauss_legendre(0.d0, pmax, p, w, NMesh)
    call mbra%setp(p, w)
    call mket%setp(p, w)
    call OpMom%init( this%OpName, mbra, mket )
    call OpMom%SetTwoBodyRelOpMeshChannel(NNint)
    call this%TransformFromMomToHOChan( OpMom )
    call OpMom%fin()
    call mbra%fin()
    call mket%fin()
    deallocate(p,w)
  end subroutine SetNNChannelSpecialOperators

  subroutine OperatorEvolution(this, params, Nmax_evolution_bra, Nmax_evolution_ket)
    use NuHamilInput, only: InputParameters
    use NNForceIsospin, only: NNForceHOIsospin
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    integer, intent(in), optional :: Nmax_evolution_bra, Nmax_evolution_ket
    integer :: Nmax_bra, Nmax_ket
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel_bra, rel_ket
    type(NNForceHOIsospin) :: vnnspin, Uspin
    type(TwoBodyRelOpIso) :: Ubra, Uket
    type(sys) :: s

    Nmax_bra = params%N2max
    Nmax_ket = params%N2max
    if(present(Nmax_evolution_bra)) Nmax_bra=Nmax_evolution_bra
    if(present(Nmax_evolution_ket)) Nmax_ket=Nmax_evolution_ket

    call relspin%init(params%hw, Nmax_bra, params%J2max_NNint)
    call vnnspin%init(relspin)
    call Uspin%init(relspin)
    call vnnspin%setNNForceHOIsospin(Uspin, params)

    call rel_bra%init(params%hw, Nmax_bra, params%Jmax2)
    call Ubra%init(rel_bra,rel_bra,s%str('UT'))
    call Ubra%SetTwoBodyScalarSpinBreaking(Uspin)

    call vnnspin%fin()
    call Uspin%fin()
    call relspin%fin()

    if( Nmax_bra == Nmax_ket ) then
      rel_ket = rel_bra
      Uket = Ubra
    end if

    if( Nmax_bra /= Nmax_ket ) then
      call relspin%init(params%hw, Nmax_ket, params%J2max_NNint)
      call vnnspin%init(relspin)
      call Uspin%init(relspin)
      call vnnspin%setNNForceHOIsospin(Uspin, params)

      call rel_ket%init(params%hw, Nmax_ket, params%Jmax2)
      call Uket%init(rel_ket,rel_ket,s%str('UT'))
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
    class(TwoBodyRelOpIso), intent(in) :: this
    integer, optional, intent(in) :: iunit
    integer :: unt = 6
    integer :: chbra, chket
    if(present(iunit)) unt = 6
    do chbra = 1, this%rel_sp_ket%GetNumberChannels()
      do chket = 1, chbra
        if( this%MatCh(chbra,chket)%is_zero ) cycle
        call this%MatCh(chbra,chket)%prtm(unt)
      end do
    end do
  end subroutine PrintTwoBodyRelOp

  subroutine PrintTwoBodyRelOpChan(this, iunit)
    class(TwoBodyRelOpChanIso), intent(in) :: this
    integer, intent(in) :: iunit
    type(TwoBodyRelChanIsoHOBasis), pointer :: chb, chk
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    integer :: bra, ket, Nbra, Nket
    chb => this%rel_ch_bra
    chk => this%rel_ch_ket
    write(*,"(a, 3i3)") "Jbra, Pbra, Tbra: ", chb%GetJ(), chb%GetParity(), chb%GetT()
    write(*,"(a, 3i3)") "Jket, Pket, Tket: ", chk%GetJ(), chk%GetParity(), chk%GetT()
    do bra = 1, chb%GetNumberStates()
      ho_bra => chb%getp(bra)
      Nbra = 2*ho_bra%GetN() + ho_bra%GetL()
      do ket = 1, chk%GetNumberStates()
        ho_ket => chk%getp(ket)
        Nket = 2*ho_ket%GetN() + ho_ket%GetL()
        write(iunit,'(6i4, f18.6)') ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(),&
            & ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), this%m(bra,ket)
      end do
    end do
  end subroutine PrintTwoBodyRelOpChan

  function GetTwoBodyRelOpIsoME(this, bra, ket) result(r)
    !
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,t1]
    ! ket = [n2,l2,s2,j2,t2]
    !
    class(TwoBodyRelOpIso), intent(in) :: this
    integer, intent(in) :: bra(5), ket(5)
    real(8) :: r
    integer :: chbra, chket, ibra, iket
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: rel_bra, rel_ket

    r = 0.d0
    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    if(bra(2) > rel_bra%GetJmax()+1) return
    if(ket(2) > rel_ket%GetJmax()+1) return
    if(bra(4) > rel_bra%GetJmax()) return
    if(ket(4) > rel_ket%GetJmax()) return
    chbra = rel_bra%GetIndex(bra(4), (-1) ** bra(2), bra(5))
    chket = rel_ket%GetIndex(ket(4), (-1) ** ket(2), ket(5))
    if(chbra*chket == 0) return
    if(this%MatCh(chbra,chket)%is_zero) return
    ibra = rel_bra%jpt(chbra)%GetIndex(bra(1), bra(2), bra(3))
    iket = rel_ket%jpt(chket)%GetIndex(ket(1), ket(2), ket(3))
    if(ibra*iket == 0) return
    r = this%MatCh(chbra,chket)%m(ibra,iket)
  end function GetTwoBodyRelOpIsoME

  subroutine SetTwoBodyRelOpIsoME(this, bra, ket, me)
    !
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,t1]
    ! ket = [n2,l2,s2,j2,t2]
    !
    class(TwoBodyRelOpIso), intent(inout) :: this
    integer, intent(in) :: bra(5), ket(5)
    real(8), intent(in) :: me
    integer :: chbra, chket, ibra, iket
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: rel_bra, rel_ket

    rel_bra => this%rel_sp_bra
    rel_ket => this%rel_sp_ket
    if(bra(2) > rel_bra%GetJmax()+1) return
    if(ket(2) > rel_ket%GetJmax()+1) return
    if(bra(4) > rel_bra%GetJmax()) return
    if(ket(4) > rel_ket%GetJmax()) return
    chbra = rel_bra%GetIndex(bra(4), (-1) ** bra(2), bra(5))
    chket = rel_ket%GetIndex(ket(4), (-1) ** ket(2), ket(5))
    if(chbra*chket == 0) return
    if(this%MatCh(chbra,chket)%is_zero) return
    ibra = rel_bra%jpt(chbra)%GetIndex(bra(1), bra(2), bra(3))
    iket = rel_ket%jpt(chket)%GetIndex(ket(1), ket(2), ket(3))
    if(ibra*iket == 0) return
    this%MatCh(chbra,chket)%m(ibra,iket) = me
  end subroutine SetTwoBodyRelOpIsoME

  subroutine TransformFromMomToHO(OpHO, OpMom)
    use TwoBodyRelOpsMeshIso
    class(TwoBodyRelOpIso), intent(inout), target :: OpHO
    type(TwoBodyRelOpMeshIso), intent(in), target :: OpMom
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: msHO_bra, msHO_ket
    type(TwoBodyRelSpaceIsoMBasis), pointer :: msMom
    type(TwoBodyRelChanIsoHOBasis), pointer :: HObra, HOket
    integer :: ichbra_ho, ichbra_mom, jbra, pbra, tbra
    integer :: ichket_ho, ichket_mom, jket, pket, tket

    msHO_bra => OpHO%rel_sp_bra
    msHO_ket => OpHO%rel_sp_ket
    msMom => OpMom%ms
    do ichbra_ho = 1, msHO_bra%GetNumberChannels()
      HObra => msHO_bra%GetChannel(ichbra_ho)
      jbra = HObra%GetJ()
      pbra = HObra%GetParity()
      tbra = HObra%GetT()
      ichbra_mom = msMom%GetIndex(jbra,pbra,tbra)

      do ichket_ho = 1, msHO_ket%GetNumberChannels()
        if( OpHO%MatCh(ichbra_ho,ichket_ho)%is_zero ) cycle
        HOket => msHO_ket%GetChannel(ichket_ho)
        jket = HOket%GetJ()
        pket = HOket%GetParity()
        tket = HOket%GetT()
        ichket_mom = msMom%GetIndex(jket,pket,tket)
        call OpHO%MatCh(ichbra_ho, ichket_ho)%TransformFromMomToHOChan( OpMom%MatCh(ichbra_mom,ichket_mom) )
      end do
    end do
  end subroutine TransformFromMomToHO

  subroutine TransformFromMomToHOChan( OpHO, OpMom )
    use TwoBodyRelOpsMeshIso
    class(TwoBodyRelOpChanIso), intent(inout), target :: OpHO
    type(TwoBodyRelOpMeshChanIso), intent(in), target :: OpMom
    type(DMat) :: ovlap_bra, ovlap_ket
    type(TwoBodyRelChanIsoHOBasis), pointer :: HObra, HOket
    type(TwoBodyRelChanIsoMBasis), pointer :: Mombra, Momket

    HObra => OpHO%rel_ch_bra
    HOket => OpHO%rel_ch_ket
    Mombra => OpMom%chbra
    Momket => OpMom%chket
    ovlap_bra = get_mom_ho_overlap( HObra, Mombra )
    ovlap_ket = get_mom_ho_overlap( HOket, Momket )
    OpHO%DMat = ovlap_bra%t() * OpMom%DMat * ovlap_ket
  end subroutine TransformFromMomToHOChan

  function get_mom_ho_overlap( HO, Mom ) result(Mat)
    use MyLibrary, only: hc, m_red_pn, ho_radial_wf_norm
    type(TwoBodyRelChanIsoHOBasis), intent(in) :: HO
    type(TwoBodyRelChanIsoMBasis), intent(in) :: Mom
    type(DMat) :: Mat
    integer :: iho, n, l, s, imom
    type(HarmonicOscillator), pointer :: ho_set
    type(Mesh), pointer :: mesh_set
    real(8) :: p, w, a
    a = hc ** 2  / (m_red_pn * HO%GetFrequency())
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

  subroutine ReducedToNonReduced(this, convJ, convT)
    class(TwoBodyRelOpIso), intent(inout) :: this
    logical, intent(in), optional :: convJ, convT
    logical :: cJ = .true.
    logical :: cT = .true.
    integer :: chbra, chket
    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. this%GetOpJ()/=0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    if( cT .and. this%GetOpT()/=0 ) then
      write(*,*) "Error, T has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%ReducedToNonReducedChan(convJ,convT)
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReduced(this, convJ, convT)
    class(TwoBodyRelOpIso), intent(inout) :: this
    logical, intent(in), optional :: convJ, convT
    logical :: cJ = .true.
    logical :: cT = .true.
    integer :: chbra, chket
    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. this%GetOpJ()/=0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    if( cT .and. this%GetOpT()/=0 ) then
      write(*,*) "Error, T has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    do chbra = 1, this%rel_sp_bra%GetNumberChannels()
      do chket = 1, this%rel_sp_ket%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%ReducedToNonReducedChan(convJ,convT)
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NonReducedToReduced

  subroutine NonReducedToReducedChan(this, convJ, convT)
    class(TwoBodyRelOpChanIso), intent(inout), target :: this
    logical, intent(in), optional :: convJ, convT
    type(TwoBodyRelChanIsoHOBasis), pointer :: chbra, chket
    logical :: cJ=.true.
    logical :: cT=.true.
    real(8) :: fact
    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. chbra%GetJ() /= chket%GetJ() ) then
      write(*,*) "Error, J' and J have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    if( cT .and. chbra%GetT() /= chket%GetT() ) then
      write(*,*) "Error, T' and T have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    fact = 1.d0
    if( cJ )  fact = fact * sqrt(dble( chket%GetJ() + 1 ))
    if( cT )  fact = fact * sqrt(dble( chket%GetT() + 1 ))
    this%m(:,:) = this%m(:,:) * fact
  end subroutine NonReducedToReducedChan

  subroutine ReducedToNonReducedChan(this, convJ, convT)
    class(TwoBodyRelOpChanIso), intent(inout), target :: this
    logical, intent(in), optional :: convJ, convT
    type(TwoBodyRelChanIsoHOBasis), pointer :: chbra, chket
    logical :: cJ=.true.
    logical :: cT=.true.
    real(8) :: fact
    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    if( present(convJ) ) cJ = convJ
    if( present(convT) ) cT = convT
    if( cJ .and. chbra%GetJ() /= chket%GetJ() ) then
      write(*,*) "Error, J' and J have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    if( cT .and. chbra%GetT() /= chket%GetT() ) then
      write(*,*) "Error, T' and T have to be equal at ", __LINE__, " in ", __FILE__
      return
    end if
    fact = 1.d0
    if( cJ )  fact = fact / sqrt(dble( chket%GetJ() + 1 ))
    if( cT )  fact = fact / sqrt(dble( chket%GetT() + 1 ))
    this%m(:,:) = this%m(:,:) * fact
  end subroutine ReducedToNonReducedChan

  subroutine ReadFile( this, filename )
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(str), intent(in) :: filename
    integer :: runit=20, io
    integer :: nbra, lbra, sbra, jbra, tbra
    integer :: nket, lket, sket, jket, tket
    real(8) :: v

    open(runit, file=filename%val, action="read")
    read(runit,*)
    read(runit,*)
    do
      read(runit,*,iostat=io) nbra, lbra, sbra, jbra, tbra, nket, lket, sket, jket, tket, v
      if(io == 0) then
        call this%Set2ME( [nbra,lbra,sbra,jbra,tbra], [nket,lket,sket,jket,tket], v )
        call this%Set2ME( [nket,lket,sket,jket,tket], [nbra,lbra,sbra,jbra,tbra], v*(-1.d0)**(jbra-jket+tbra-tket) )
      end if
      if(io < 0) exit
      if(io > 0) then
        write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      end if
    end do
    close(runit)
  end subroutine ReadFile

  subroutine WriteFile( this, filename )
    class(TwoBodyRelOpIso), intent(inout) :: this
    type(str), intent(in) :: filename
    type(str) :: OpName
    integer :: wunit=20
    integer :: ichbra, ichket, ibra, iket, iketmax
    integer :: jbra, jket, tbra, tket
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    character(255) :: header

    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    open(wunit, file=filename%val, action="write")
    write(wunit,'(a)') trim(header)
    write(wunit,'(a)') "#  n'  l'  s'  j'  t'   n   l  s   j   t           ME "
    do ichbra = 1, this%rel_sp_bra%GetNumberCHannels()
      do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        jbra = this%MatCh(ichbra,ichket)%rel_ch_bra%GetJ()
        jket = this%MatCh(ichbra,ichket)%rel_ch_ket%GetJ()
        tbra = this%MatCh(ichbra,ichket)%rel_ch_bra%GetT()
        tket = this%MatCh(ichbra,ichket)%rel_ch_ket%GetT()
        do ibra = 1, this%MatCh(ichbra,ichket)%rel_ch_bra%GetNumberStates()
          iketmax = this%MatCh(ichbra,ichket)%rel_ch_ket%GetNumberStates()
          if(ichbra == ichket) iketmax = ibra
          do iket = 1, iketmax
            ho_bra => this%MatCh(ichbra,ichket)%rel_ch_bra%getp(ibra)
            ho_ket => this%MatCh(ichbra,ichket)%rel_ch_ket%getp(iket)
            if(abs(this%MatCh(ichbra,ichket)%m(ibra,iket)) < 1.d-16) cycle
            write(wunit,'(10i4,es18.8)') ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS(), jbra, tbra, &
                & ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS(), jket, tket, this%MatCh(ichbra,ichket)%m(ibra,iket)
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine WriteFile
end module TwoBodyRelOpsIso
