module TwoBodyRelCMOps
  use LinAlgLib
  use Profiler, only: timer
  use omp_lib
  use ClassSys
  use TwoBodyRelCMChanHO
  use TwoBodyRelativeSpace
  use OperatorDefinitions
  implicit none
  public :: TwoBodyRelCMOpChan
  public :: TwoBodyRelCMOp
  private

  type, extends(DMat) :: TwoBodyRelCMOpChan
    type(TwoBodyRelCMChanHOBasis), pointer :: rel_ch_bra, rel_ch_ket
    type(str) :: OpName
    logical :: is_zero = .true.
  contains
    procedure :: InitTwoBodyRelCMOpChan
    procedure :: FinTwoBodyRelCMOpChan
    procedure :: PrintTwoBodyRelCMOpChan
    procedure :: OperatorEvolutionChannel
    procedure :: TruncateChannel

    generic :: init => InitTwoBodyRelCMOpChan
    generic :: release => FinTwoBodyRelCMOpChan
    generic :: prtm => PrintTwoBodyRelCMOpChan
    generic :: evolve => OperatorEvolutionChannel
    generic :: truncate => TruncateChannel
  end type TwoBodyRelCMOpChan

  type, extends(OperatorDef) :: TwoBodyRelCMOp
    type(TwoBodyRelCMOpChan), allocatable :: MatCh(:,:)
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms=>null()
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyRelCMOp
    procedure :: InitTwoBodyRelCMOpFromString
    procedure :: FinTwoBodyRelCMOp
    procedure :: CopyTwoBodyRelCMOp
    procedure :: SumTwoBodyRelCMOp
    procedure :: SubtractTwoBodyRelCMOp
    procedure :: ScaleTwoBodyRelCMOp
    procedure :: PrintTwoBodyRelCMOp
    procedure :: GetTwoBodyRelCMOpME
    procedure :: SetTwoBodyRelCMOpME
    procedure :: SetTwoBodyRelCMOp
    procedure :: GetMemory
    procedure :: IsZero
    procedure :: ReadFile
    procedure :: WriteFile
    procedure :: Truncate
    procedure :: OperatorEvolution

    generic :: assignment(=) => CopyTwoBodyRelCMOp
    generic :: operator(+) => SumTwoBodyRelCMOp
    generic :: operator(-) => SubtractTwoBodyRelCMOp
    generic :: operator(*) => ScaleTwoBodyRelCMOp

    generic :: init => InitTwoBodyRelCMOp, InitTwoBodyRelCMOpFromString
    generic :: fin => FinTwoBodyRelCMOp
    generic :: prt => PrintTwoBodyRelCMOp
    generic :: Get2ME => GetTwoBodyRelCMOpME
    generic :: Set2ME => SetTwoBodyRelCMOpME
    generic :: set => SetTwoBodyRelCMOp
    generic :: evolve => OperatorEvolution
  end type TwoBodyRelCMOp

contains

  function GetMemory(this) result(mem)
    class(TwoBodyRelCMOp), intent(in) :: this
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
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(TwoBodyRelCMChanHOBasis), pointer :: chb, chk
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

  subroutine InitTwoBodyRelCMOpChan(this,OpName,ch_bra,ch_ket)
    class(TwoBodyRelCMOpChan), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelCMChanHOBasis), target, intent(in) :: ch_bra, ch_ket
    this%is_zero = .false.
    this%rel_ch_bra => ch_bra
    this%rel_ch_ket => ch_ket
    this%OpName = OpName
    if( OpName%val == "UT" .or. OpName%val == "UnitaryTransformation" ) then
      call this%eye( ch_ket%GetNumberStates() )
    else
      call this%zeros( ch_bra%GetNumberStates() ,ch_ket%GetNumberStates() )
    end if
  end subroutine InitTwoBodyRelCMOpChan

  subroutine FinTwoBodyRelCMOpChan(this)
    class(TwoBodyRelCMOpChan), intent(inout) :: this
    if(.not. allocated(this%m)) return
    this%is_zero = .true.
    this%rel_ch_bra => null()
    this%rel_ch_ket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelCMOpChan

  subroutine InitTwoBodyRelCMOpFromString(this, ms, oprtr)
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceHOBasis), intent(in) :: ms
    type(str), intent(in) :: oprtr

    if(allocated(this%MatCh)) call this%fin()
    call this%InitOpDef(oprtr, .true. )
    call this%init(ms,this%GetOpJ(),this%GetOpP(),this%GetOpZ())
  end subroutine InitTwoBodyRelCMOpFromString

  subroutine InitTwoBodyRelCMOp(this, ms, jr, pr, zr)
    use MyLibrary, only: triag
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceHOBasis), target, intent(in) :: ms
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
    integer, intent(in) :: jr, pr, zr
    integer :: ichb, ichk, nb, nk
    integer :: jb, pb, tzb, lcmb, jrelb
    integer :: jk, pk, tzk, lcmk, jrelk
    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%InitOpDef(.true., jr, pr, zr )
    allocate(this%MatCh(ms%GetNumberChannels(), ms%GetNumberChannels() ))
    do ichb = 1, this%ms%GetNumberChannels()
      chbra => this%ms%GetChannel(ichb)
      nb = chbra%GetNumberStates()
      jb = chbra%GetJ()
      pb = chbra%GetParity()
      tzb= chbra%GetZ()
      jrelb = chbra%GetJRel()
      lcmb = chbra%GetLCM()

      do ichk = 1, this%ms%GetNumberChannels()
        chket => this%ms%GetChannel(ichk)
        nk = chket%GetNumberStates()
        jk = chket%GetJ()
        pk = chket%GetParity()
        tzk= chket%GetZ()
        jrelk = chket%GetJRel()
        lcmk = chket%GetLCM()
        if(triag(jb, jk, jr)) cycle
        if(pb * pr /= pk) cycle
        if(abs(tzk-tzb)/=zr) cycle
        call this%MatCh(ichb,ichk)%init(this%GetOpName(),chbra,chket)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelCMOp

  subroutine FinTwoBodyRelCMOp(this)
    class(TwoBodyRelCMOp), intent(inout) :: this
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
  end subroutine FinTwoBodyRelCMOp

  subroutine CopyTwoBodyRelCMOp(a, b)
    class(TwoBodyRelCMOp), intent(inout) :: a
    class(TwoBodyRelCMOp), intent(in) :: b
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
  end subroutine CopyTwoBodyRelCMOp

  function SumTwoBodyRelCMOp(a, b) result(c)
    use LinAlgLib
    type(TwoBodyRelCMOp) :: c
    class(TwoBodyRelCMOp), intent(in) :: a, b
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
  end function SumTwoBodyRelCMOp

  function SubtractTwoBodyRelCMOp(a, b) result(c)
    type(TwoBodyRelCMOp) :: c
    class(TwoBodyRelCMOp), intent(in) :: a, b
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
  end function SubtractTwoBodyRelCMOp

  function ScaleTwoBodyRelCMOp(a, b) result(c)
    type(TwoBodyRelCMOp) :: c
    class(TwoBodyRelCMOp), intent(in) :: a
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
  end function ScaleTwoBodyRelCMOp

  subroutine PrintTwoBodyRelCMOp(this,iunit)
    class(TwoBodyRelCMOp), intent(in) :: this
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
  end subroutine PrintTwoBodyRelCMOp

  subroutine PrintTwoBodyRelCMOpChan(this,iunit)
    class(TwoBodyRelCMOpChan), intent(in), target :: this
    type(TwoBodyRelCMChanHOBasis), pointer :: chb, chk
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
            & ho_bra%GetNCM(), chb%GetLCM(), chb%GetJ(), chb%GetZ(), &
            & ho_ket%GetNrel(), ho_ket%GetLrel(), ho_ket%GetSpin(), chk%GetJrel(), &
            & ho_ket%GetNCM(), chk%GetLCM(), chk%GetJ(), chk%GetZ(), &
            & this%m(bra,ket)
      end do
    end do
  end subroutine PrintTwoBodyRelCMOpChan

  subroutine SetTwoBodyRelCMOpME(this, bra, ket, me)
    !
    ! < bra || Op || ket >
    ! bra = [ncm1,lcm1,nr1,lr1,s1,jr1,j1,z1]
    ! ket = [ncm2,lcm2,nr2,lr2,s2,jr2,j2,z2]
    !
    class(TwoBodyRelCMOp), intent(inout) :: this
    integer, intent(in) :: bra(8), ket(8)
    real(8), intent(in) :: me
    integer :: ichbra, ichket, ibra, iket
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket

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
  end subroutine SetTwoBodyRelCMOpME

  function GetTwoBodyRelCMOpME(this, bra, ket) result(r)
    !
    ! < bra || Op || ket >
    ! bra = [ncm1,lcm1,nr1,lr1,s1,jr1,j1,z1]
    ! ket = [ncm2,lcm2,nr2,lr2,s2,jr2,j2,z2]
    !
    class(TwoBodyRelCMOp), intent(in) :: this
    integer, intent(in) :: bra(8), ket(8)
    real(8) :: r
    integer :: ichbra, ichket, ibra, iket
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket

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
  end function GetTwoBodyRelCMOpME

  subroutine SetTwoBodyRelCMOp(this)
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    integer :: ichbra, ichket
    integer :: ibra, iket, iket_max
    type(OperatorDef) :: op
    real(8) :: me

    call op%InitOpDef(this%GetOpName(), .true.)
    ms => this%ms
    !$omp parallel
    !$omp do private(ichbra, ichket, chbra, chket, ibra, iket_max, iket, hobra, hoket, me) schedule(dynamic)
    do ichbra = 1, ms%GetNumberChannels()
      do ichket = 1, ichbra
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        do ibra = 1, chbra%GetNumberStates()
          iket_max = chket%GetNumberStates()
          if(ichbra==ichket) iket_max = ibra
          do iket = 1, iket_max
            hobra => chbra%GetHOs(ibra)
            hoket => chket%GetHOs(iket)
            me = CalcMERelCM(op, &
                & [hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), chbra%GetJrel(), &
                & hobra%GetNCM(), chbra%GetLCM(), chbra%GetJ(), chbra%GetZ()], &
                & [hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), chket%GetJrel(), &
                & hoket%GetNCM(), chket%GetLCM(), chket%GetJ(), chket%GetZ()])
            this%MatCh(ichbra,ichket)%m(ibra,iket) = me
            this%MatCh(ichket,ichbra)%m(iket,ibra) = me * (-1.d0)**(chbra%GetJ()-chket%GetJ())
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine SetTwoBodyRelCMOp

  function TruncateChannel(this, chbra_new, chket_new) result(op)
    class(TwoBodyRelCMOpChan), intent(in) :: this
    type(TwoBodyRelCMChanHOBasis), intent(in) :: chbra_new, chket_new
    type(TwoBodyRelCMOpChan) :: op
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra_old, chket_old
    integer :: ibra, iket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket
    integer :: ibra_old, iket_old

    call op%init(this%OpName, chbra_new, chket_new)
    chbra_old => this%rel_ch_bra
    chket_old => this%rel_ch_ket
    do ibra = 1, chbra_new%GetNumberStates()
      do iket = 1, chket_new%GetNumberStates()
        hobra => chbra_new%GetHOs(ibra)
        hoket => chket_new%GetHOs(iket)
        ibra_old = chbra_old%GetIndex(hobra)
        iket_old = chket_old%GetIndex(hoket)
        if(ibra_old * iket_old == 0) cycle
        op%m(ibra,iket) = this%m(ibra_old,iket_old)
      end do
    end do
  end function TruncateChannel

  function Truncate(this, ms_new) result(op)
    class(TwoBodyRelCMOp), intent(in) :: this
    type(TwoBodyRelCMSpaceHOBasis), target, intent(in) :: ms_new
    type(TwoBodyRelCMOp) :: op
    type(TwoBodyRelCMOpChan) :: opch
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
    integer :: ibra, iket
    type(RelCMHarmonicOscillator) :: hobra, hoket
    integer :: ibra_old, iket_old
    type(str) :: tmp

    call op%init(ms_new, this%GetOpName())
    do ibra = 1, ms_new%GetNumberChannels()
      do iket = 1, ms_new%GetNumberChannels()
        chbra => ms_new%GetChannel(ibra)
        chket => ms_new%GetChannel(iket)
        ibra_old = this%ms%GetIndex(chbra%GetJ(), chbra%GetParity(), chbra%GetZ(), chbra%GetJrel(), chbra%GetLcm())
        iket_old = this%ms%GetIndex(chket%GetJ(), chket%GetParity(), chket%GetZ(), chket%GetJrel(), chket%GetLcm())
        if(ibra_old * iket_old == 0) cycle
        if(this%MatCh(ibra_old, iket_old)%is_zero) cycle
        opch = this%MatCh(ibra_old, iket_old)%Truncate(chbra, chket)
        op%MatCh(ibra,iket)%DMat = opch%DMat
      end do
    end do
  end function Truncate

  subroutine OperatorEvolutionChannel(this, U)
    use TwoBodyRelOps
    class(TwoBodyRelCMOpChan), intent(inout) :: this
    type(TwoBodyRelOp), intent(in) :: U
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
    type(TwoBodyRelCMOpChan) :: Ubra, Uket
    type(sys) :: s
    integer :: jrel_bra, jrel_ket, prel_bra, prel_ket, ibra, iket, zbra, zket
    type(RelCMHarmonicOscillator), pointer :: hobra, hoket

    chbra => this%rel_ch_bra
    chket => this%rel_ch_ket
    call Ubra%init(s%str("UT"), chbra, chbra)
    call Uket%init(s%str("UT"), chket, chket)

    jrel_bra = chbra%GetJRel()
    prel_bra = chbra%GetParity() * (-1)**chbra%GetLCM()
    zbra = chbra%GetZ()

    jrel_ket = chket%GetJRel()
    prel_ket = chket%GetParity() * (-1)**chket%GetLCM()
    zket = chket%GetZ()

    do ibra = 1, chbra%GetNumberStates()
      do iket = 1, chbra%GetNumberStates()
        hobra => chbra%GetHOs(ibra)
        hoket => chbra%GetHOs(iket)
        if(hobra%GetNcm() /= hoket%GetNcm()) cycle
        Ubra%m(ibra,iket) = U%Get2ME([hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), jrel_bra, zbra], &
          & [hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), jrel_bra, zbra])
      end do
    end do

    do ibra = 1, chket%GetNumberStates()
      do iket = 1, chket%GetNumberStates()
        hobra => chket%GetHOs(ibra)
        hoket => chket%GetHOs(iket)
        if(hobra%GetNcm() /= hoket%GetNcm()) cycle
        Uket%m(ibra,iket) = U%Get2ME([hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), jrel_ket, zket], &
          & [hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), jrel_ket, zket])
      end do
    end do
    this%DMat = Ubra%DMat%t() * this%DMat * Uket%DMat
    call Ubra%fin()
    call Uket%fin()
  end subroutine OperatorEvolutionChannel

  subroutine OperatorEvolution(this, U)
    use TwoBodyRelOps
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(TwoBodyRelOp), intent(in) :: U
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms
    integer :: ibra, iket
    real(8) :: ti
    type(sys) :: s
    ti = omp_get_wtime()
    ms => this%ms
    do ibra = 1, ms%GetNumberChannels()
      do iket = 1, ms%GetNumberChannels()
        if(this%MatCh(ibra, iket)%is_zero) cycle
        call this%MatCh(ibra, iket)%evolve(U)
      end do
    end do
    call timer%add(s%str('cm-dep op evolution'), omp_get_wtime() - ti)
  end subroutine OperatorEvolution

  subroutine ReadFile( this, filename )
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(str), intent(in) :: filename
    integer :: runit=20, io
    integer :: ncm_bra, lcm_bra, nrel_bra, lrel_bra, s_bra, jrel_bra, Jbra, zbra
    integer :: ncm_ket, lcm_ket, nrel_ket, lrel_ket, s_ket, jrel_ket, Jket, zket
    real(8) :: v

    open(runit, file=filename%val, action="read")
    read(runit,*)
    read(runit,*)
    do
      !read(runit,*,iostat=io) ncm_bra, lcm_bra, nrel_bra, lrel_bra, s_bra, jrel_bra, Jbra, zbra, &
      !    & ncm_ket, lcm_ket, nrel_ket, lrel_ket, s_ket, jrel_ket, Jket, zket, v
      read(runit,*,iostat=io) nrel_bra, lrel_bra, s_bra, jrel_bra, ncm_bra, lcm_bra, Jbra, zbra, &
          & nrel_ket, lrel_ket, s_ket, jrel_ket, ncm_ket, lcm_ket, Jket, zket, v
      if(io == 0) then
        call this%SetTwoBodyRelCMOpME( [ncm_bra,lcm_bra,nrel_bra,lrel_bra,s_bra,jrel_bra,Jbra,zbra], &
            & [ncm_ket,lcm_ket,nrel_ket,lrel_ket,s_ket,jrel_ket,Jket,zket], v)
        call this%SetTwoBodyRelCMOpME( [ncm_ket,lcm_ket,nrel_ket,lrel_ket,s_ket,jrel_ket,Jket,zket], &
            & [ncm_bra,lcm_bra,nrel_bra,lrel_bra,s_bra,jrel_bra,Jbra,zbra], v*(-1.d0)**(Jbra-Jket))
      end if
      if(io < 0) exit
      if(io > 0) then
        write(*,*) "Error, at ", __LINE__, " in ", __FILE__
      end if
    end do
    close(runit)
  end subroutine ReadFile

  subroutine WriteFile( this, filename )
    class(TwoBodyRelCMOp), intent(inout) :: this
    type(str), intent(in) :: filename
    type(str) :: OpName
    integer :: wunit=20
    integer :: ichbra, ichket, ibra, iket, iketmax
    type(TwoBodyRelCMSpaceHOBasis), pointer :: ms
    type(TwoBodyRelCMChanHOBasis), pointer :: chbra, chket
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
    write(wunit,"(a)") "#  n', l', s', j', N', L', J', z', n,  l,  s,  j,  N,  L,  J,  z,                 ME"
    do ichbra = 1, this%ms%GetNumberCHannels()
      do ichket = 1, ichbra
        if( this%MatCh(ichbra,ichket)%is_zero ) cycle
        chbra => ms%GetChannel(ichbra)
        chket => ms%GetChannel(ichket)
        do ibra = 1, chbra%GetNumberStates()
          hobra => chbra%GetHOs(ibra)
          iketmax = chket%GetNumberStates()
          if(ichbra == ichket) iketmax = ibra
          do iket = 1, iketmax
            hoket => chket%GetHOs(iket)
            if(abs(this%MatCh(ichbra,ichket)%m(ibra,iket)) < 1.d-16) cycle
            write(wunit,'(16i4,es18.8)') &
                & hobra%GetNrel(), hobra%GetLrel(), hobra%GetSpin(), chbra%GetJrel(), &
                & hobra%GetNCM(), chbra%GetLCM(), chbra%GetJ(), chbra%GetZ(), &
                & hoket%GetNrel(), hoket%GetLrel(), hoket%GetSpin(), chket%GetJrel(), &
                & hoket%GetNCM(), chket%GetLCM(), chket%GetJ(), chket%GetZ(), &
                & this%MatCh(ichbra,ichket)%m(ibra,iket)
          end do
        end do
      end do
    end do
    close(wunit)
  end subroutine WriteFile
end module TwoBodyRelCMOps
