module TwoBodyRelCMIsoOps
  use LinAlgLib
  use omp_lib
  use ClassSys
  use TwoBodyRelativeSpace
  use OperatorDefinitions
  implicit none
  public :: TwoBodyRelCMIsoOpChan
  public :: TwoBodyRelCMIsoOp
  private

  type, extends(DMat) :: TwoBodyRelCMIsoOpChan
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(str) :: OpName
    logical :: is_zero = .true.
  contains
    procedure :: InitTwoBodyRelCMIsoOpChan
    procedure :: FinTwoBodyRelCMIsoOpChan
    procedure :: PrintTwoBodyRelCMIsoOpChan
    procedure :: MakeIsospinSymmetricChannel
    generic :: init => InitTwoBodyRelCMIsoOpChan
    generic :: release => FinTwoBodyRelCMIsoOpChan
    generic :: prtm => PrintTwoBodyRelCMIsoOpChan
  end type TwoBodyRelCMIsoOpChan

  type, extends(OperatorDef) :: TwoBodyRelCMIsoOp
    type(TwoBodyRelCMIsoOpChan), allocatable :: MatCh(:,:)
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms=>null()
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyRelCMIsoOp
    procedure :: InitTwoBodyRelCMIsoOpFromString
    procedure :: FinTwoBodyRelCMIsoOp
    procedure :: CopyTwoBodyRelCMIsoOp
    procedure :: SumTwoBodyRelCMIsoOp
    procedure :: SubtractTwoBodyRelCMIsoOp
    procedure :: ScaleTwoBodyRelCMIsoOp
    procedure :: PrintTwoBodyRelCMIsoOp
    procedure :: GetTwoBodyRelCMIsoOpME
    procedure :: SetTwoBodyRelCMIsoOpME
    procedure :: SetTwoBodyRelCMIsoOp
    procedure :: MakeIsospinSymmetric
    procedure :: SetFromRelOp
    procedure :: GetMemory
    procedure :: IsZero
    procedure :: ReadFile
    procedure :: WriteFile

    generic :: assignment(=) => CopyTwoBodyRelCMIsoOp
    generic :: operator(+) => SumTwoBodyRelCMIsoOp
    generic :: operator(-) => SubtractTwoBodyRelCMIsoOp
    generic :: operator(*) => ScaleTwoBodyRelCMIsoOp

    generic :: init => InitTwoBodyRelCMIsoOp, InitTwoBodyRelCMIsoOpFromString
    generic :: fin => FinTwoBodyRelCMIsoOp
    generic :: prt => PrintTwoBodyRelCMIsoOp
    generic :: Get2ME => GetTwoBodyRelCMIsoOpME
    generic :: Set2ME => SetTwoBodyRelCMIsoOpME
    generic :: set => SetTwoBodyRelCMIsoOp, SetFromRelOp
  end type TwoBodyRelCMIsoOp

contains

  function GetMemory(this) result(mem)
    class(TwoBodyRelCMIsoOp), intent(in) :: this
    real(8) :: mem
    integer :: ichbra, ichket
    mem = 0.d0
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        mem = mem + this%MatCh(ichbra,ichket)%n_row * this%MatCh(ichbra,ichket)%n_col * 8.d0 / (1024.d0**3)
      end do
    end do
  end function GetMemory

  subroutine IsZero(this)
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chb, chk
    integer :: ichbra, ichket
    integer :: bra, ket
    real(8) :: norm, me_scale
    me_scale = -1.d0
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        me_scale = max(me_scale, abs(maxval(this%MatCh(ichbra,ichket)%m)))
      end do
    end do

    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        chb => this%ms%GetChannel(ichbra)
        chk => this%ms%GetChannel(ichket)
        norm = 0.d0
        do bra = 1, chb%GetNumberStates()
          do ket = 1, chk%GetNumberStates()
            norm = norm + this%MatCh(ichbra,ichket)%m(bra,ket)**2
          end do
        end do
        if(sqrt(norm/(dble(this%MatCh(ichbra,ichket)%n_row * this%MatCh(ichbra,ichket)%n_col))) > 1.d-10*me_scale) cycle
        call this%MatCh(ichbra,ichket)%release()
      end do
    end do
  end subroutine IsZero

  subroutine InitTwoBodyRelCMIsoOpChan(this,OpName,ch_bra,ch_ket)
    class(TwoBodyRelCMIsoOpChan), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelCMChanIsoHOBasis), target, intent(in) :: ch_bra, ch_ket
    this%is_zero = .false.
    this%rel_ch_bra => ch_bra
    this%rel_ch_ket => ch_ket
    this%OpName = OpName
    if( OpName%val == "UT" .or. OpName%val == "UnitaryTransformation" ) then
      call this%eye( ch_ket%GetNumberStates() )
    else
      call this%zeros( ch_bra%GetNumberStates() ,ch_ket%GetNumberStates() )
    end if
  end subroutine InitTwoBodyRelCMIsoOpChan

  subroutine FinTwoBodyRelCMIsoOpChan(this)
    class(TwoBodyRelCMIsoOpChan), intent(inout) :: this
    if(.not. allocated(this%m)) return
    this%is_zero = .true.
    this%rel_ch_bra => null()
    this%rel_ch_ket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelCMIsoOpChan

  subroutine InitTwoBodyRelCMIsoOpFromString(this, ms, oprtr)
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceIsoHOBasis), intent(in) :: ms
    type(str), intent(in) :: oprtr

    if(allocated(this%MatCh)) call this%fin()
    call this%InitOpDef(oprtr, .false. )
    call this%init(ms,this%GetOpJ(),this%GetOpP(),this%GetOpT())
  end subroutine InitTwoBodyRelCMIsoOpFromString

  subroutine InitTwoBodyRelCMIsoOp(this, ms, jr, pr, zr)
    use MyLibrary, only: triag
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceIsoHOBasis), target, intent(in) :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket
    integer, intent(in) :: jr, pr, zr
    integer :: ichb, ichk, nb, nk
    integer :: jb, pb, tzb, lcmb, jrelb
    integer :: jk, pk, tzk, lcmk, jrelk
    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%InitOpDef(.false., jr, pr, zr )
    allocate(this%MatCh(ms%GetNumberChannels(), ms%GetNumberChannels() ))
    do ichb = 1, this%ms%GetNumberChannels()
      chbra => this%ms%GetChannel(ichb)
      nb = chbra%GetNumberStates()
      jb = chbra%GetJ()
      pb = chbra%GetParity()
      tzb= chbra%GetT()
      jrelb = chbra%GetJRel()
      lcmb = chbra%GetLCM()

      do ichk = 1, this%ms%GetNumberChannels()
        chket => this%ms%GetChannel(ichk)
        nk = chket%GetNumberStates()
        jk = chket%GetJ()
        pk = chket%GetParity()
        tzk= chket%GetT()
        jrelk = chket%GetJRel()
        lcmk = chket%GetLCM()
        if(triag(jb, jk, jr)) cycle
        if(pb * pr /= pk) cycle
        if(triag(tzk,tzb,zr)) cycle
        call this%MatCh(ichb,ichk)%init(this%GetOpName(),chbra,chket)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelCMIsoOp

  subroutine FinTwoBodyRelCMIsoOp(this)
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    integer :: ichb, ichk
    if(.not. this%is_init) return
    do ichb = 1, size(this%MatCh, 1)
      do ichk = 1, size(this%MatCh, 2)
        call this%MatCh(ichb, ichk)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%ms => null()
    this%ms => null()
    call this%FinOperatorDef()
    this%is_init = .false.
  end subroutine FinTwoBodyRelCMIsoOp

  subroutine CopyTwoBodyRelCMIsoOp(a, b)
    class(TwoBodyRelCMIsoOp), intent(inout) :: a
    class(TwoBodyRelCMIsoOp), intent(in) :: b
    integer :: ichb, ichk, n
    if(allocated(a%MatCh)) call a%fin()
    n = size(b%MatCh, 1)
    call a%CopyOperatorDef(b%OperatorDef)
    a%ms => b%ms
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
    a%is_init = b%is_init
  end subroutine CopyTwoBodyRelCMIsoOp

  function SumTwoBodyRelCMIsoOp(a, b) result(c)
    use LinAlgLib
    type(TwoBodyRelCMIsoOp) :: c
    class(TwoBodyRelCMIsoOp), intent(in) :: a, b
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
        c%MatCh(ichb,ichk)%DMat = a%MatCh(ichb, ichk)%DMat + b%MatCh(ichb, ichk)%DMat
      end do
    end do
  end function SumTwoBodyRelCMIsoOp

  function SubtractTwoBodyRelCMIsoOp(a, b) result(c)
    type(TwoBodyRelCMIsoOp) :: c
    class(TwoBodyRelCMIsoOp), intent(in) :: a, b
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
  end function SubtractTwoBodyRelCMIsoOp

  function ScaleTwoBodyRelCMIsoOp(a, b) result(c)
    type(TwoBodyRelCMIsoOp) :: c
    class(TwoBodyRelCMIsoOp), intent(in) :: a
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
  end function ScaleTwoBodyRelCMIsoOp

  subroutine PrintTwoBodyRelCMIsoOp(this,iunit)
    class(TwoBodyRelCMIsoOp), intent(in) :: this
    integer, optional, intent(in) :: iunit
    integer :: unt = 6
    integer :: chbra, chket

    if(present(iunit)) unt = iunit

    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%prtm(unt)
      end do
    end do
  end subroutine PrintTwoBodyRelCMIsoOp

  subroutine PrintTwoBodyRelCMIsoOpChan(this,iunit)
    class(TwoBodyRelCMIsoOpChan), intent(in), target :: this
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chb, chk
    type(RelCMHarmonicOscillator), pointer :: ho_bra, ho_ket
    integer, intent(in) :: iunit
    integer :: bra, ket
    chb => this%rel_ch_bra
    chk => this%rel_ch_ket
    write(iunit,"(a)") "#  n', l', s', j', N', L', J', z', n,  l,  s,  j,  N,  L,  J,  z,                 ME"
    do bra = 1, chb%GetNumberStates()
      ho_bra => chb%GetHOs(bra)
      do ket = 1, chk%GetNumberStates()
        ho_ket => chk%GetHOs(ket)
        write(iunit,'(16i4, es18.6)') &
            & ho_bra%GetNrel(), ho_bra%GetLrel(), ho_bra%GetSpin(), chb%GetJrel(), &
            & ho_bra%GetNCM(), chb%GetLCM(), chb%GetJ(), chb%GetT(), &
            & ho_ket%GetNrel(), ho_ket%GetLrel(), ho_ket%GetSpin(), chk%GetJrel(), &
            & ho_ket%GetNCM(), chk%GetLCM(), chk%GetJ(), chk%GetT(), &
            & this%m(bra,ket)
      end do
    end do
  end subroutine PrintTwoBodyRelCMIsoOpChan

  subroutine SetTwoBodyRelCMIsoOpME(this, bra, ket, me)
    !
    ! < bra || Op || ket >
    ! bra = [ncm1,lcm1,nr1,lr1,s1,jr1,j1,z1]
    ! ket = [ncm2,lcm2,nr2,lr2,s2,jr2,j2,z2]
    !
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    integer, intent(in) :: bra(8), ket(8)
    real(8), intent(in) :: me
    integer :: ichbra, ichket, ibra, iket
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket

    ms => this%ms
    if(2*bra(1)+bra(2) > ms%GetNmaxCM()) return
    if(2*ket(1)+ket(2) > ms%GetNmaxCM()) return
    if(2*bra(3)+bra(4) > ms%GetNmaxRel()) return
    if(2*ket(3)+ket(4) > ms%GetNmaxRel()) return
    if(2*(bra(1)+bra(3))+bra(2)+bra(4) > ms%GetNmax()) return
    if(2*(ket(1)+ket(3))+ket(2)+ket(4) > ms%GetNmax()) return
    if(bra(2) > ms%GetLcmMax()) return
    if(ket(2) > ms%GetLcmMax()) return
    if(bra(7) > ms%GetJMax()) return
    if(ket(7) > ms%GetJMax()) return
    if(bra(6) > ms%GetJRelMax()) return
    if(ket(6) > ms%GetJRelMax()) return
    if(bra(4) > ms%GetJRelMax()+1) return
    if(ket(4) > ms%GetJRelMax()+1) return
    ichbra = ms%GetIndex(bra(7), (-1) ** (bra(2)+bra(4)), bra(8), bra(6), bra(2))
    ichket = ms%GetIndex(ket(7), (-1) ** (ket(2)+ket(4)), ket(8), ket(6), ket(2))
    if(ichbra*ichket == 0) return
    if(this%MatCh(ichbra,ichket)%is_zero) return
    chbra => ms%GetChannel(ichbra)
    chket => ms%GetChannel(ichket)
    ibra = chbra%GetIndex(bra(3), bra(4), bra(5), bra(1))
    iket = chket%GetIndex(ket(3), ket(4), ket(5), ket(1))
    if(ibra*iket == 0) return
    this%MatCh(ichbra,ichket)%m(ibra,iket) = me
  end subroutine SetTwoBodyRelCMIsoOpME

  function GetTwoBodyRelCMIsoOpME(this, bra, ket) result(r)
    !
    ! < bra || Op || ket >
    ! bra = [ncm1,lcm1,nr1,lr1,s1,jr1,j1,z1]
    ! ket = [ncm2,lcm2,nr2,lr2,s2,jr2,j2,z2]
    !
    class(TwoBodyRelCMIsoOp), intent(in) :: this
    integer, intent(in) :: bra(8), ket(8)
    real(8) :: r
    integer :: ichbra, ichket, ibra, iket
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket

    r = 0.d0
    ms => this%ms
    if(2*bra(1)+bra(2) > ms%GetNmaxCM()) return
    if(2*ket(1)+ket(2) > ms%GetNmaxCM()) return
    if(2*bra(3)+bra(4) > ms%GetNmaxRel()) return
    if(2*ket(3)+ket(4) > ms%GetNmaxRel()) return
    if(2*(bra(1)+bra(3))+bra(2)+bra(4) > ms%GetNmax()) return
    if(2*(ket(1)+ket(3))+ket(2)+ket(4) > ms%GetNmax()) return
    if(bra(2) > ms%GetLcmMax()) return
    if(ket(2) > ms%GetLcmMax()) return
    if(bra(7) > ms%GetJMax()) return
    if(ket(7) > ms%GetJMax()) return
    if(bra(6) > ms%GetJRelMax()) return
    if(ket(6) > ms%GetJRelMax()) return
    if(bra(4) > ms%GetJRelMax()+1) return
    if(ket(4) > ms%GetJRelMax()+1) return
    ichbra = ms%GetIndex(bra(7), (-1) ** (bra(2)+bra(4)), bra(8), bra(6), bra(2))
    ichket = ms%GetIndex(ket(7), (-1) ** (ket(2)+ket(4)), ket(8), ket(6), ket(2))
    if(ichbra*ichket == 0) return
    if(this%MatCh(ichbra,ichket)%is_zero) return
    chbra => ms%GetChannel(ichbra)
    chket => ms%GetChannel(ichket)
    ibra = chbra%GetIndex(bra(3), bra(4), bra(5), bra(1))
    iket = chket%GetIndex(ket(3), ket(4), ket(5), ket(1))
    if(ibra*iket == 0) return
    r = this%MatCh(ichbra,ichket)%m(ibra,iket)
  end function GetTwoBodyRelCMIsoOpME

  subroutine SetTwoBodyRelCMIsoOp(this)
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    integer :: ichbra, ichket
    integer :: ibra, iket, iket_max
    type(OperatorDef) :: op
    real(8) :: me

    call op%InitOpDef(this%GetOpName(), .false.)
    ms => this%ms
    !$omp parallel
    !$omp do private(ichbra, ichket, chbra, chket, ibra, iket_max, iket, hobra, hoket, me) schedule(dynamic)
    do ichbra = 1, ms%GetNumberChannels()
      !do ichket = 1, ichbra
      do ichket = 1, ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        do ibra = 1, chbra%GetNumberStates()
          iket_max = chket%GetNumberStates()
          !if(ichbra==ichket) iket_max = ibra
          do iket = 1, iket_max
            hobra => chbra%GetHOs(ibra)
            hoket => chket%GetHOs(iket)
            me = CalcMERelCM(op, &
                & [hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), chbra%GetJrel(), &
                & hobra%GetNCM(), chbra%GetLCM(), chbra%GetJ(), chbra%GetT()], &
                & [hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), chket%GetJrel(), &
                & hoket%GetNCM(), chket%GetLCM(), chket%GetJ(), chket%GetT()])
            this%MatCh(ichbra,ichket)%m(ibra,iket) = me
            !this%MatCh(ichket,ichbra)%m(iket,ibra) = me * &
            !    & (-1.d0)**(chbra%GetJ()-chket%GetJ()+chbra%GetT()-chket%GetT())
            ! write(*,*) ichbra, ichket, ibra, iket, me, &
            !     & me*(-1.d0)**(chbra%GetJ()-chket%GetJ()+chbra%GetT()-chket%GetT())
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetTwoBodyRelCMIsoOp

  subroutine SetFromRelOp(this, op)
    use MyLibrary, only: sjs
    use TwoBodyRelOpsIso
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    integer :: ichbra, ichket
    integer :: ibra, iket, iket_max
    integer :: rankJ
    real(8) :: me

    rankJ = op%GetOpJ()
    ms => this%ms
    !$omp parallel
    !$omp do private(ichbra, ichket, chbra, chket, ibra, iket_max, iket, hobra, hoket, me) schedule(dynamic)
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ichbra
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        this%MatCh(ichbra,ichket)%m(:,:) = 0.d0
        this%MatCh(ichket,ichbra)%m(:,:) = 0.d0
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        if(chbra%GetLCM() /= chket%GetLCM()) cycle
        do ibra = 1, chbra%GetNumberStates()
          iket_max = chket%GetNumberStates()
          if(ichbra==ichket) iket_max = ibra
          do iket = 1, iket_max
            hobra => chbra%GetHOs(ibra)
            hoket => chket%GetHOs(iket)
            if(hobra%GetNCM() /= hoket%GetNCM()) cycle
            me = op%Get2ME([hobra%GetNRel(), hobra%GetLRel(), hobra%GetSpin(), chbra%GetJrel(), chbra%GetT()], &
                & [hoket%GetNRel(), hoket%GetLRel(), hoket%GetSpin(), chket%GetJrel(), chket%GetT()]) * &
                & (-1.d0)**(chbra%GetJ()+chket%GetLCM()+chket%GetJrel()+rankJ) * &
                & sqrt(dble((2*chbra%GetJ()+1)*(2*chket%GetJ()+1))) * &
                & sjs(2*chbra%GetJRel(), 2*chbra%GetJ(), 2*chket%GetLCM(), 2*chket%GetJ(), 2*chket%GetJRel(), 2*rankJ)
            this%MatCh(ichbra,ichket)%m(ibra,iket) = me
            this%MatCh(ichket,ichbra)%m(iket,ibra) = me * &
                & (-1.d0)**(chbra%GetJ()-chket%GetJ()+chbra%GetT()-chket%GetT())
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetFromRelOp

  subroutine ReadFile( this, filename )
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(str), intent(in) :: filename
    integer :: runit=20, io
    integer :: ncm_bra, lcm_bra, nrel_bra, lrel_bra, s_bra, jrel_bra, Jbra, zbra
    integer :: ncm_ket, lcm_ket, nrel_ket, lrel_ket, s_ket, jrel_ket, Jket, zket
    real(8) :: v, phase

    open(runit, file=filename%val, action="read")
    read(runit,*)
    read(runit,*)
    do
      !read(runit,*,iostat=io) ncm_bra, lcm_bra, nrel_bra, lrel_bra, s_bra, jrel_bra, Jbra, zbra, &
      !    & ncm_ket, lcm_ket, nrel_ket, lrel_ket, s_ket, jrel_ket, Jket, zket, v
      !read(runit,*,iostat=io) nrel_bra, lrel_bra, s_bra, jrel_bra, ncm_bra, lcm_bra, Jbra, zbra, &
      !    & nrel_ket, lrel_ket, s_ket, jrel_ket, ncm_ket, lcm_ket, Jket, zket, v
      read(runit,*,iostat=io) nrel_bra, ncm_bra, lrel_bra, s_bra, jrel_bra, zbra, lcm_bra, Jbra, &
          & nrel_ket, ncm_ket, lrel_ket, s_ket, jrel_ket, zket, lcm_ket, Jket, v
      phase = 1.d0
      phase = (-1.d0)**(lcm_bra+jrel_bra-Jbra+lcm_ket+jrel_ket-Jket)
      if(io == 0) then
        call this%SetTwoBodyRelCMIsoOpME( [ncm_bra,lcm_bra,nrel_bra,lrel_bra,s_bra,jrel_bra,Jbra,zbra], &
            & [ncm_ket,lcm_ket,nrel_ket,lrel_ket,s_ket,jrel_ket,Jket,zket], v*phase)
        !call this%SetTwoBodyRelCMIsoOpME( [ncm_ket,lcm_ket,nrel_ket,lrel_ket,s_ket,jrel_ket,Jket,zket], &
        !    & [ncm_bra,lcm_bra,nrel_bra,lrel_bra,s_bra,jrel_bra,Jbra,zbra], v*phase*(-1.d0)**(Jbra-Jket+zbra-zket))
      end if
      if(io < 0) exit
      if(io > 0) then
        write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      end if
    end do
    close(runit)
  end subroutine ReadFile

  subroutine WriteFile( this, filename )
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(str), intent(in) :: filename
    type(str) :: OpName
    integer :: wunit=20
    integer :: ichbra, ichket, ibra, iket, iketmax
    type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    character(255) :: header

    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    ms => this%ms
    open(wunit, file=filename%val, action="write")
    write(wunit,'(a)') trim(header)
    write(wunit,"(a)") "#  n', l', s', j', N', L', J', t', n,  l,  s,  j,  N,  L,  J,  t,                 ME"
    do ichbra = 1, this%ms%GetNumberCHannels()
      do ichket = 1, this%ms%GetNumberCHannels()
      !do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        do ibra = 1, chbra%GetNumberStates()
          hobra => chbra%GetHOs(ibra)
          iketmax = chket%GetNumberStates()
          !if(ichbra == ichket) iketmax = ibra
          do iket = 1, iketmax
            hoket => chket%GetHOs(iket)
            if(abs(this%MatCh(ichbra,ichket)%m(ibra,iket)) < 1.d-8) cycle
            write(wunit,'(16i4,es18.8)') &
                & hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), chbra%GetJrel(), &
                & hobra%GetNCM(), chbra%GetLCM(), chbra%GetJ(), chbra%GetT(), &
                & hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), chket%GetJrel(), &
                & hoket%GetNCM(), chket%GetLCM(), chket%GetJ(), chket%GetT(), &
                & this%MatCh(ichbra,ichket)%m(ibra,iket)
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine WriteFile

  subroutine MakeIsospinSymmetric(this, op, flip_projection)
    !
    ! Operator should have good isospin symmetry
    !
    use TwoBodyRelCMOps, only: TwoBodyRelCMOp
    class(TwoBodyRelCMIsoOp), intent(inout) :: this
    type(TwoBodyRelCMOp), intent(in) :: op
    integer, intent(in), optional :: flip_projection
    integer :: chbra, chket

    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
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
    use TwoBodyRelCMOps, only: TwoBodyRelCMOp
    class(TwoBodyRelCMIsoOpChan), intent(inout) :: this
    type(TwoBodyRelCMOp), intent(in) :: op
    integer, intent(in) :: Trank
    integer, intent(in), optional :: flip_projection
    type(TwoBodyRelCMChanIsoHOBasis), pointer :: chbra, chket
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra_op, chket_op
    type(RelCMHarmonicOscillator), pointer :: ho_bra, ho_ket
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
    ichbra_op = op%ms%GetIndex(chbra%GetJ(), chbra%GetParity(), zbra, chbra%GetJRel(), chbra%GetLCM())
    ichket_op = op%ms%GetIndex(chket%GetJ(), chket%GetParity(), zket, chket%GetJRel(), chket%GetLCM())
    if(ichbra_op * ichket_op ==0 ) return
    if(op%MatCh(ichbra_op,ichket_op)%is_zero) return
    chbra_op => op%MatCh(ichbra_op,ichket_op)%rel_ch_bra
    chket_op => op%MatCh(ichbra_op,ichket_op)%rel_ch_ket
    !$omp parallel
    !$omp do private(ibra, ho_bra, iket, ho_ket, idx_op_bra, idx_op_ket, v)
    do ibra = 1, chbra%GetNumberStates()
      ho_bra => chbra%GetHOs(ibra)
      do iket = 1, chket%GetNumberStates()
        ho_ket => chket%GetHOs(iket)
        idx_op_bra = chbra_op%GetIndex(ho_bra)
        idx_op_ket = chket_op%GetIndex(ho_ket)
        v = op%MatCh(ichbra_op, ichket_op)%m(idx_op_bra, idx_op_ket)
        this%m(ibra,iket) = v * fact
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine MakeIsospinSymmetricChannel
end module TwoBodyRelCMIsoOps
