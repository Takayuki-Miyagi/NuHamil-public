module TwoBodyRelCMOpsMesh
  use LinAlgLib
  use omp_lib
  use ClassSys
  use TwoBodyRelativeSpace
  use OperatorDefinitions
  implicit none
  public :: TwoBodyRelCMOpMeshChan
  public :: TwoBodyRelCMOpMesh
  private

  type, extends(DMat) :: TwoBodyRelCMOpMeshChan
    type(TwoBodyRelCMChanMBasis), pointer :: chbra, chket
    type(str) :: OpName
    logical :: is_zero = .true.
  contains
    procedure :: InitTwoBodyRelCMOpMeshChan
    procedure :: FinTwoBodyRelCMOpMeshChan
    !procedure :: SetMultiPoleOperatorChannel
    generic :: init => InitTwoBodyRelCMOpMeshChan
    generic :: release => FinTwoBodyRelCMOpMeshChan
    !generic :: SetMultiPoleOperator => SetMultiPoleOperatorChannel
  end type TwoBodyRelCMOpMeshChan

  type, extends(OperatorDef) :: TwoBodyRelCMOpMesh
    type(TwoBodyRelCMOpMeshChan), allocatable :: MatCh(:,:)
    type(TwoBodyRelCMSpaceMBasis), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyRelCMOpMeshFromString
    procedure :: InitTwoBodyRelCMOpMesh
    procedure :: FinTwoBodyRelCMOpMesh
    procedure :: CopyTwoBodyRelCMOpMesh
    procedure :: SumTwoBodyRelCMOpMesh
    procedure :: SubtractTwoBodyRelCMOpMesh
    procedure :: ScaleTwoBodyRelCMOpMesh
    procedure :: PrintTwoBodyRelCMOpMesh
    generic :: assignment(=) => CopyTwoBodyRelCMOpMesh
    generic :: operator(+) => SumTwoBodyRelCMOpMesh
    generic :: operator(-) => SubtractTwoBodyRelCMOpMesh
    generic :: operator(*) => ScaleTwoBodyRelCMOpMesh
    generic :: init => InitTwoBodyRelCMOpMeshFromString, InitTwoBodyRelCMOpMesh
    generic :: fin => FinTwoBodyRelCMOpMesh
    generic :: prt => PrintTwoBodyRelCMOpMesh
  end type TwoBodyRelCMOpMesh
contains
  subroutine InitTwoBodyRelCMOpMeshChan(this,OpName,chbra,chket)
    class(TwoBodyRelCMOpMeshChan), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelCMChanMBasis), target, intent(in) :: chbra, chket
    this%is_zero = .false.
    this%chbra => chbra
    this%chket => chket
    this%OpName = OpName
    call this%zeros( chbra%GetNumberStates() ,chket%GetNumberStates() )
  end subroutine InitTwoBodyRelCMOpMeshChan

  subroutine FinTwoBodyRelCMOpMeshChan(this)
    class(TwoBodyRelCMOpMeshChan), intent(inout) :: this
    if( this%is_zero ) return
    this%is_zero = .true.
    this%chbra => null()
    this%chket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelCMOpMeshChan

  subroutine InitTwoBodyRelCMOpMeshFromString(this, ms, OpName)
    class(TwoBodyRelCMOpMesh), intent(inout) :: this
    type(TwoBodyRelCMSpaceMBasis), target, intent(in) :: ms
    type(str), intent(in) :: OpName

    if( allocated(this%MatCh) ) call this%fin()
    call this%InitOpDef( OpName, .true. )
    call this%init( ms, this%GetOpJ(), this%GetOpP(), this%GetOpZ() )
  end subroutine InitTwoBodyRelCMOpMeshFromString

  subroutine InitTwoBodyRelCMOpMesh(this, ms, jr, pr, zr)
    use MyLibrary, only: triag
    class(TwoBodyRelCMOpMesh), intent(inout) :: this
    type(TwoBodyRelCMSpaceMBasis), target, intent(in) :: ms
    integer, intent(in) :: jr, pr, zr
    integer :: ichbra, jbra, pbra, zbra, jrelbra, lcmbra
    integer :: ichket, jket, pket, zket, jrelket, lcmket
    type(TwoBodyRelCMChanMBasis), pointer :: chbra, chket
    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%InitOpDef( .true., jr, pr, zr )
    allocate( this%MatCh( ms%GetNumberChannels(), ms%GetNumberChannels() ))
    do ichbra = 1, ms%GetNumberChannels()
      chbra => ms%GetChannel(ichbra)
      jbra = chbra%GetJ()
      pbra = chbra%GetParity()
      zbra = chbra%GetZ()
      jrelbra = chbra%GetJRel()
      lcmbra = chbra%GetLCM()
      do ichket = 1, ms%GetNumberChannels()
        chket => ms%GetChannel(ichket)
        jket = chket%GetJ()
        pket = chket%GetParity()
        zket = chket%GetZ()
        jrelket = chket%GetJRel()
        lcmket = chket%GetLCM()
        if(triag( jbra, jket, jr )) cycle
        if(pbra * pket * pr == -1) cycle
        if(abs(zket-zbra) /= zr) cycle
        call this%MatCh(ichbra,ichket)%init(this%GetOpName(),chbra,chket)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelCMOpMesh

  subroutine FinTwoBodyRelCMOpMesh(this)
    use MyLibrary, only: triag
    class(TwoBodyRelCMOpMesh), intent(inout) :: this
    integer :: ichbra, ichket
    if(.not. this%is_init ) return
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        call this%MatCh(ichbra,ichket)%release()
      end do
    end do
    deallocate( this%MatCh )
    this%ms => null()
    call this%FinOperatorDef()
    this%is_init = .false.
  end subroutine FinTwoBodyRelCMOpMesh

  subroutine CopyTwoBodyRelCMOpMesh(a, b)
    class(TwoBodyRelCMOpMesh), intent(inout) :: a
    type(TwoBodyRelCMOpMesh), intent(in) :: b
    integer :: ichbra, ichket
    if( allocated(a%MatCh) ) call a%fin()
    call a%CopyOperatorDef( b%OperatorDef )
    a%ms => b%ms
    allocate( a%MatCh( a%ms%GetNumberChannels(), a%ms%GetNumberChannels() ) )
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( b%MatCh(ichbra,ichket)%is_zero ) then
          a%MatCh(ichbra,ichket)%is_zero = .true.
          cycle
        end if
        a%MatCh(ichbra,ichket)%is_zero = .false.
        a%MatCh(ichbra,ichket)%DMat = b%MatCh(ichbra,ichket)%DMat
        a%MatCh(ichbra,ichket)%chbra => b%MatCh(ichbra,ichket)%chbra
        a%MatCh(ichbra,ichket)%chket => b%MatCh(ichbra,ichket)%chket
      end do
    end do
    a%is_init = b%is_init
  end subroutine CopyTwoBodyRelCMOpMesh

  function SumTwoBodyRelCMOpMesh(a, b) result(c)
    class(TwoBodyRelCMOpMesh), intent(in) :: a
    type(TwoBodyRelCMOpMesh), intent(in) :: b
    type(TwoBodyRelCMOpMesh) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat + b%MatCh(ichbra,ichket)%DMat
      end do
    end do
  end function SumTwoBodyRelCMOpMesh

  function SubtractTwoBodyRelCMOpMesh(a, b) result(c)
    class(TwoBodyRelCMOpMesh), intent(in) :: a
    type(TwoBodyRelCMOpMesh), intent(in) :: b
    type(TwoBodyRelCMOpMesh) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat - b%MatCh(ichbra,ichket)%DMat
      end do
    end do
  end function SubtractTwoBodyRelCMOpMesh

  function ScaleTwoBodyRelCMOpMesh(a, b) result(c)
    class(TwoBodyRelCMOpMesh), intent(in) :: a
    real(8), intent(in) :: b
    type(TwoBodyRelCMOpMesh) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat * b
      end do
    end do
  end function ScaleTwoBodyRelCMOpMesh

  subroutine PrintTwoBodyRelCMOpMesh(this,iunit)
    class(TwoBodyRelCMOpMesh), intent(in) :: this
    integer, optional, intent(in) :: iunit
    integer :: chbra, chket, unt=6

    if(present(iunit)) unt = iunit
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        call this%MatCh(chbra,chket)%prt()
      end do
    end do
  end subroutine PrintTwoBodyRelCMOpMesh

!  subroutine SetMultiPoleOperator(this, ms, Q, rankJ, rankZ, OpName)
!    class(TwoBodyRelCMOpMesh), intent(inout) :: this
!    type(TwoBodyRelCMSpaceMBasis), intent(in) :: ms
!    real(8), intent(in) :: Q
!    integer, intent(in) :: rankJ, rankZ
!    type(str), intent(in) :: OpName
!    integer :: ichbra, ichket
!
!    select case(OpName%val)
!    case("L5", "Tel5", "Tmag")
!      call this%InitOpDef(.true., rankJ, (-1)**(rankJ+1), rankZ)
!      call this%SetOpName(OpName)
!      call this%SetReduced(.true.)
!    case("L", "Tmag5", "Tel")
!      call this%InitOpDef(.true., rankJ, (-1)**rankJ, rankZ)
!      call this%SetOpName(OpName)
!      call this%SetReduced(.true.)
!    end select
!    call this%init(ms, this%GetOpJ(), this%GetOpP(), this%GetOpZ())
!
!    do ichbra = 1, ms%GetNumberChannels()
!      do ichket = 1, ichbra
!        if(this%MatCh(ichbra,ichket)%is_zero) cycle
!        call this%MatCh(ichbra,ichket)%SetMultiPoleOperator(Q, rankJ, rankZ, OpName)
!      end do
!    end do
!  end subroutine SetMultiPoleOperator

!  subroutine SetMultiPoleOperatorChannel(this, Q, rankJ, rankZ, OpName)
!    use MyLibrary, only: hc
!    use MultiPoleResponse2
!    class(TwoBodyRelCMOpMeshChan), intent(inout) :: this
!    real(8), intent(in) :: Q
!    integer, intent(in) :: rankJ, rankZ
!    type(str), intent(in) :: OpName
!    type(TwoBodyRelCMChanMBasis), pointer :: chbra, chket
!    integer :: ibra, iket, op_type
!    type(RelCMMesh), pointer :: mbra, mket
!    real(8) :: me
!
!    if(rankZ/=0) then
!      write(*,*) "Not implemented; the Rank of Tz has to be 0."
!      return
!    end if
!    if(OpName%val=="L5")    op_type = 0
!    if(OpName%val=="Tel5")  op_type = 1
!    if(OpName%val=="Tmag5") op_type = 2
!    chbra => this%chbra
!    chket => this%chket
!    do ibra = 1, chbra%GetNumberStates()
!      mbra => chbra%GetRelCMMesh(ibra)
!      do iket = 1, chket%GetNumberStates()
!        mket => chket%GetRelCMMesh(iket)
!        if(abs(mbra%GetPCM()-mket%GetPCM())*hc>Q) cycle
!        if(   (mbra%GetPCM()+mket%GetPCM())*hc<Q) cycle
!        me = calc_me(op_type, rankJ, &
!            & mbra%GetPCM()*hc, chbra%GetLCM(), mbra%GetPRel()*hc, mbra%GetLRel(), mbra%GetSpin(), chbra%GetJRel(), chbra%GetJ(), &
!            & mket%GetPCM()*hc, chket%GetLCM(), mket%GetPRel()*hc, mket%GetLRel(), mket%GetSpin(), chket%GetJRel(), chket%GetJ())
!        this%m(ibra,iket) = me * hc**3
!      end do
!    end do
!    !call this%prt('')
!  end subroutine SetMultiPoleOperatorChannel
!
!  subroutine TestTwoBodyRelCMOpMesh()
!    use MyLibrary, only: gauss_legendre
!    type(TwoBodyRelCMOpMesh) :: op
!    type(TwoBodyRelCMSpaceMBasis) :: ms
!    type(sys) :: s
!    real(8), allocatable :: x_rel(:), w_rel(:), x_cm(:), w_cm(:)
!    call ms%init(0.d0, 1.d0, 10, 0.d0, 1.d0, 10, 1, 1, 0)
!    call gauss_legendre(0.d0, 1.d0, x_rel, w_rel, 10)
!    call gauss_legendre(0.d0, 1.d0, x_cm, w_cm, 10)
!    call ms%SetMeshWeight(x_rel, w_rel, x_cm, w_cm)
!    call op%SetMultiPoleOperator(ms, 10.d0, 1, 0, s%str('L5'))
!  end subroutine TestTwoBodyRelCMOpMesh

end module TwoBodyRelCMOpsMesh
