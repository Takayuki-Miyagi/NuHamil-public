module TwoBodyRelOpsMeshIso
  use LinAlgLib
  use omp_lib
  use ClassSys
  use TwoBodyRelativeSpace
  use OperatorDefinitions
  implicit none
  public :: TwoBodyRelOpMeshChanIso, TwoBodyRelOpMeshIso
  private

  type, extends(DMat) :: TwoBodyRelOpMeshChanIso
    type(TwoBodyRelChanIsoMBasis), pointer :: chbra, chket
    type(str) :: OpName
    logical :: is_zero = .true.
  contains
    procedure :: InitTwoBodyRelOpMeshChan
    procedure :: FinTwoBodyRelOpMeshChan
    procedure :: SetTwoBodyRelOpMeshChannel

    generic :: init => InitTwoBodyRelOpMeshChan
    generic :: release => FinTwoBodyRelOpMeshChan
  end type TwoBodyRelOpMeshChanIso

  type, extends(OperatorDef) :: TwoBodyRelOpMeshIso
    type(TwoBodyRelOpMeshChanIso), allocatable :: MatCh(:,:)
    type(TwoBodyRelSpaceIsoMBasis), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitTwoBodyRelOpMeshFromString
    procedure :: InitTwoBodyRelOpMesh
    procedure :: FinTwoBodyRelOpMesh
    procedure :: CopyTwoBodyRelOpMesh
    procedure :: SumTwoBodyRelOpMesh
    procedure :: SubtractTwoBodyRelOpMesh
    procedure :: ScaleTwoBodyRelOpMesh
    procedure :: SetTwoBodyRelOpMesh
    procedure :: PrintTwoBodyRelOpMesh

    generic :: assignment(=) => CopyTwoBodyRelOpMesh
    generic :: operator(+) => SumTwoBodyRelOpMesh
    generic :: operator(-) => SubtractTwoBodyRelOpMesh
    generic :: operator(*) => ScaleTwoBodyRelOpMesh
    generic :: init => InitTwoBodyRelOpMeshFromString, InitTwoBodyRelOpMesh
    generic :: fin => FinTwoBodyRelOpMesh
    generic :: set => SetTwoBodyRelOpMesh
    generic :: prt => PrintTwoBodyRelOpMesh
  end type TwoBodyRelOpMeshIso
contains
  subroutine InitTwoBodyRelOpMeshChan(this,OpName,chbra,chket)
    class(TwoBodyRelOpMeshChanIso), intent(inout) :: this
    type(str), intent(in) :: OpName
    type(TwoBodyRelChanIsoMBasis), target, intent(in) :: chbra, chket
    this%is_zero = .false.
    this%chbra => chbra
    this%chket => chket
    this%OpName = OpName
    call this%zeros( chbra%GetNumberStates() ,chket%GetNumberStates() )
  end subroutine InitTwoBodyRelOpMeshChan

  subroutine FinTwoBodyRelOpMeshChan(this)
    class(TwoBodyRelOpMeshChanIso), intent(inout) :: this
    if( this%is_zero ) return
    this%is_zero = .true.
    this%chbra => null()
    this%chket => null()
    call this%DMat%fin()
  end subroutine FinTwoBodyRelOpMeshChan

  subroutine SetTwoBodyRelOpMeshChannel(this, NNint)
    use OperatorDefinitions, only: parse_operator_string, parse_opname, set_pv_couplings, set_pvtv_couplings
    use LeptonNumberViolation, only: LNViolation
    use ParityViolation, only: PViolation
    use ParityTimeReViolation, only: PTViolation
    use AxialVectorQ0, only: AxialCurrentQ0
    use VectorQ0, only: VectorCurrentQ0
    use NuHamilInput, only: params
    class(TwoBodyRelOpMeshChanIso), intent(inout) :: this
    type(str), intent(in) :: NNint
    type(LNViolation) :: lnv
    type(PViolation) :: pv
    type(PTViolation) :: ptv
    type(AxialCurrentQ0) :: acq0
    type(VectorCurrentQ0) :: vcq0
    integer :: ch_order, reg_pow, i
    logical :: meson_exchange
    real(8) :: Ec, lam, pv_couplings(11), pvtv_couplings(10)
    type(str) :: OpName, op, regulator, op_string
    type(str), allocatable :: par_names(:), par_vals(:)
    type(sys) :: s

    op_string = this%OpName
    call parse_operator_string(op_string, OpName, ch_order, regulator, reg_pow, lam)
    call parse_opname(OpName, op, par_names, par_vals)

    if(s%find(OpName,s%str("0vbb")) .or. s%find(OpName,s%str("0veeCap"))) then
      Ec = 0.d0
      do i = 1, size(par_names)
        if(par_names(i)%val=="Ec") Ec = s%dbl(par_vals(i))
      end do
      call lnv%init(op%val, NNint%val, rank=0, Ec=Ec, ch_order=ch_order, &
          & Jmax=max(this%chbra%GetJ(), this%chket%GetJ()), Regulator=regulator%val, RegulatorPower=reg_pow, RegulatorLambda=lam)
      call lnv%show_lnv_options()
      call inside(lnv)
      call lnv%fin()
      return
    end if

    if(s%find(OpName,s%str("PVTC"))) then
      pv_couplings(:) = 0.d0
      meson_exchange = .false.
      do i = 1, size(par_names)
        if(par_names(i)%val=="ChEFT" .or. par_names(i)%val=="DDH") then
          call set_pv_couplings(pv_couplings, par_names(i)%val, par_vals(i)%val)
          if(par_names(i)%val=="DDH") meson_exchange=.true.
        end if
      end do
      call pv%init(Op%val, NNint%val, pv_couplings, rank=0, ch_order=ch_order, &
          & Jmax=max(this%chbra%GetJ(), this%chket%GetJ()), Regulator=regulator%val, RegulatorPower=reg_pow, &
          & RegulatorLambda=lam, meson_exchange=meson_exchange)
      call pv%show_pv_options()
      call inside(pv)
      call pv%fin()
      return
    end if

    if(s%find(OpName,s%str("PVTV"))) then
      pvtv_couplings(:) = 0.d0
      meson_exchange = .false.
      do i = 1, size(par_names)
        if(par_names(i)%val=="ChEFT" .or. par_names(i)%val=="OBE") then
          call set_pvtv_couplings(pvtv_couplings, par_names(i)%val, par_vals(i)%val)
          if(par_names(i)%val=="OBE") meson_exchange=.true.
        end if
      end do
      call ptv%init(Op%val, NNint%val, pvtv_couplings, rank=0, ch_order=ch_order, &
          & Jmax=max(this%chbra%GetJ(), this%chket%GetJ()), Regulator=regulator%val, RegulatorPower=reg_pow, &
          & RegulatorLambda=lam, meson_exchange=meson_exchange)
      call ptv%show_pvtv_options()
      call inside(ptv)
      call ptv%fin()
      return
    end if

    if(s%find(OpName,s%str("AxialV"))) then
      call acq0%init(op%val, [params%c1, params%c3, params%c4, params%cD], ch_order=ch_order, &
          & Jmax=max(this%chbra%GetJ(), this%chket%GetJ()), Regulator=regulator%val, &
          & RegulatorPower=reg_pow, RegulatorLambda=lam)
      call acq0%show_options()
      call inside(acq0)
      call acq0%fin()
      return
    end if

    if(s%find(OpName,s%str("Vector"))) then
      call vcq0%init(op%val, NNint%val, ch_order=ch_order, &
          & Jmax=max(this%chbra%GetJ(), this%chket%GetJ()), Regulator=regulator%val, &
          & RegulatorPower=reg_pow, RegulatorLambda=lam)
      call vcq0%show_options()
      call inside(vcq0)
      call vcq0%fin()
      return
    end if
  contains
    subroutine inside(op_class)
      use MyLibrary, only: hc
      class(*), intent(in) :: op_class
      type(TwoBodyRelChanIsoMBasis), pointer :: chbra, chket
      integer :: jbra, tbra, lbra, sbra, bra
      integer :: jket, tket, lket, sket, ket
      type(Mesh), pointer :: mbra, mket
      real(8) :: pbra, pket, me
      select type(op_class)
      type is (LNViolation)
#include "include/block_for_rel_op_mesh_channel_iso.inc"
      type is (PViolation)
#include "include/block_for_rel_op_mesh_channel_iso.inc"
      type is (PTViolation)
#include "include/block_for_rel_op_mesh_channel_iso.inc"
      type is (AxialCurrentQ0)
#include "include/block_for_rel_op_mesh_channel_iso.inc"
      type is (VectorCurrentQ0)
#include "include/block_for_rel_op_mesh_channel_iso.inc"
      end select
    end subroutine inside
  end subroutine SetTwoBodyRelOpMeshChannel

  subroutine InitTwoBodyRelOpMeshFromString(this, ms, OpName)
    class(TwoBodyRelOpMeshIso), intent(inout) :: this
    type(TwoBodyRelSpaceIsoMBasis), target, intent(in) :: ms
    type(str), intent(in) :: OpName

    if( allocated(this%MatCh) ) call this%fin()
    call this%InitOpDef( OpName, .false. )
    call this%init( ms, this%GetOpJ(), this%GetOpP(), this%GetOpT() )
  end subroutine InitTwoBodyRelOpMeshFromString

  subroutine InitTwoBodyRelOpMesh(this, ms, jr, pr, tr)
    use MyLibrary, only: triag
    class(TwoBodyRelOpMeshIso), intent(inout) :: this
    type(TwoBodyRelSpaceIsoMBasis), target, intent(in) :: ms
    integer, intent(in) :: jr, pr, tr
    integer :: ichbra, jbra, pbra, tbra
    integer :: ichket, jket, pket, tket
    type(TwoBodyRelChanIsoMBasis), pointer :: chbra, chket
    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%InitOpDef( .false., jr, pr, tr )
    allocate( this%MatCh( ms%GetNumberChannels(), ms%GetNumberChannels() ))
    do ichbra = 1, ms%GetNumberChannels()
      chbra => ms%GetChannel(ichbra)
      jbra = chbra%GetJ()
      pbra = chbra%GetParity()
      tbra = chbra%GetT()
      do ichket = 1, ms%GetNumberChannels()
        chket => ms%GetChannel(ichket)
        jket = chket%GetJ()
        pket = chket%GetParity()
        tket = chket%GetT()
        if(triag( jbra, jket, jr )) cycle
        if(pbra * pket * pr == -1) cycle
        if(triag( tbra, tket, tr )) cycle
        call this%MatCh(ichbra,ichket)%init(this%GetOpName(),chbra,chket)
      end do
    end do
    this%is_init = .true.
  end subroutine InitTwoBodyRelOpMesh

  subroutine FinTwoBodyRelOpMesh(this)
    class(TwoBodyRelOpMeshIso), intent(inout) :: this
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
  end subroutine FinTwoBodyRelOpMesh

  subroutine CopyTwoBodyRelOpMesh(a, b)
    class(TwoBodyRelOpMeshIso), intent(inout) :: a
    type(TwoBodyRelOpMeshIso), intent(in) :: b
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
  end subroutine CopyTwoBodyRelOpMesh

  function SumTwoBodyRelOpMesh(a, b) result(c)
    class(TwoBodyRelOpMeshIso), intent(in) :: a
    type(TwoBodyRelOpMeshIso), intent(in) :: b
    type(TwoBodyRelOpMeshIso) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat + b%MatCh(ichbra,ichket)%DMat
      end do
    end do
  end function SumTwoBodyRelOpMesh

  function SubtractTwoBodyRelOpMesh(a, b) result(c)
    class(TwoBodyRelOpMeshIso), intent(in) :: a
    type(TwoBodyRelOpMeshIso), intent(in) :: b
    type(TwoBodyRelOpMeshIso) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat - b%MatCh(ichbra,ichket)%DMat
      end do
    end do
  end function SubtractTwoBodyRelOpMesh

  function ScaleTwoBodyRelOpMesh(a, b) result(c)
    class(TwoBodyRelOpMeshIso), intent(in) :: a
    real(8), intent(in) :: b
    type(TwoBodyRelOpMeshIso) :: c
    integer :: ichbra, ichket
    c = a
    do ichbra = 1, a%ms%GetNumberChannels()
      do ichket = 1, a%ms%GetNumberChannels()
        if( c%MatCh(ichbra,ichket)%is_zero ) cycle
        c%MatCh(ichbra,ichket)%DMat = a%MatCh(ichbra,ichket)%DMat * b
      end do
    end do
  end function ScaleTwoBodyRelOpMesh

  subroutine SetTwoBodyRelOpMesh(this, NNint)
    use OperatorDefinitions, only: parse_operator_string, parse_opname, set_pv_couplings, set_pvtv_couplings
    use LeptonNumberViolation, only: LNViolation
    use ParityViolation, only: PViolation
    use ParityTimeReViolation, only: PTViolation
    use AxialVectorQ0, only: AxialCurrentQ0
    use VectorQ0, only: VectorCurrentQ0
    use NuHamilInput, only: params
    class(TwoBodyRelOpMeshIso), intent(inout), target :: this
    type(str), intent(in) :: NNint
    type(LNViolation) :: lnv
    type(PViolation) :: pv
    type(PTViolation) :: ptv
    type(AxialCurrentQ0) :: acq0
    type(VectorCurrentQ0) :: vcq0
    integer :: ch_order, reg_pow, i
    logical :: meson_exchange
    real(8) :: Ec, lam, pv_couplings(11), pvtv_couplings(10)
    type(str) :: OpName, op, regulator, op_string
    type(str), allocatable :: par_names(:), par_vals(:)
    type(sys) :: s

    op_string = this%GetOpName()
    call parse_operator_string(op_string, OpName, ch_order, regulator, reg_pow, lam)
    call parse_opname(OpName, op, par_names, par_vals)

    if(s%find(OpName,s%str("0vbb")) .or. s%find(OpName,s%str("0veeCap"))) then
      Ec = 0.d0
      do i = 1, size(par_names)
        if(par_names(i)%val=="Ec") Ec = s%dbl(par_vals(i))
      end do
      call lnv%init(op%val, NNint%val, rank=0, Ec=Ec, ch_order=ch_order, &
          & Jmax=this%ms%GetJmax(), Regulator=regulator%val, RegulatorPower=reg_pow, RegulatorLambda=lam)
      call lnv%show_lnv_options()
      call inside(lnv)
      call lnv%fin()
    end if

    if(s%find(OpName,s%str("PVTC"))) then
      pv_couplings(:) = 0.d0
      meson_exchange = .false.
      do i = 1, size(par_names)
        if(par_names(i)%val=="ChEFT" .or. par_names(i)%val=="DDH") then
          call set_pv_couplings(pv_couplings, par_names(i)%val, par_vals(i)%val)
          if(par_names(i)%val=="DDH") meson_exchange=.true.
        end if
      end do
      call pv%init(Op%val, NNint%val, pv_couplings, rank=0, ch_order=ch_order, &
          & Jmax=this%ms%GetJmax(), Regulator=regulator%val, RegulatorPower=reg_pow, &
          & RegulatorLambda=lam, meson_exchange=meson_exchange)
      call pv%show_pv_options()
      call inside(pv)
      call pv%fin()
      return
    end if

    if(s%find(OpName,s%str("PVTV"))) then
      pvtv_couplings(:) = 0.d0
      meson_exchange = .false.
      do i = 1, size(par_names)
        if(par_names(i)%val=="ChEFT" .or. par_names(i)%val=="OBE") then
          call set_pvtv_couplings(pvtv_couplings, par_names(i)%val, par_vals(i)%val)
          if(par_names(i)%val=="OBE") meson_exchange=.true.
        end if
      end do
      call ptv%init(Op%val, NNint%val, pvtv_couplings, rank=0, ch_order=ch_order, &
          & Jmax=this%ms%GetJmax(), Regulator=regulator%val, RegulatorPower=reg_pow, &
          & RegulatorLambda=lam, meson_exchange=meson_exchange)
      call ptv%show_pvtv_options()
      call inside(ptv)
      call ptv%fin()
      return
    end if

    if(s%find(OpName,s%str("AxialV"))) then
      call acq0%init(op%val, [params%c1, params%c3, params%c4, params%cD], ch_order=ch_order, &
          & Jmax=this%ms%GetJmax(), Regulator=regulator%val, RegulatorPower=reg_pow, RegulatorLambda=lam)
      call acq0%show_options()
      call inside(acq0)
      call acq0%fin()
      return
    end if

    if(s%find(OpName,s%str("Vector"))) then
      call vcq0%init(op%val, NNint%val, ch_order=ch_order, &
          & Jmax=this%ms%GetJmax(), Regulator=regulator%val, RegulatorPower=reg_pow, RegulatorLambda=lam)
      call vcq0%show_options()
      call inside(vcq0)
      call vcq0%fin()
      return
    end if

  contains
    subroutine inside(op_class)
      use MyLibrary, only: hc
      class(*), intent(in) :: op_class
      class(TwoBodyRelOpMeshChanIso), pointer :: op_ch
      type(TwoBodyRelChanIsoMBasis), pointer :: chbra, chket
      integer :: ichbra, ichket
      integer :: jbra, tbra, lbra, sbra, bra
      integer :: jket, tket, lket, sket, ket
      type(Mesh), pointer :: mbra, mket
      real(8) :: pbra, pket, me
      select type(op_class)
      type is (LNViolation)
#include "include/block_for_rel_op_mesh_iso.inc"
      type is (PViolation)
#include "include/block_for_rel_op_mesh_iso.inc"
      type is (PTViolation)
#include "include/block_for_rel_op_mesh_iso.inc"
      type is (AxialCurrentQ0)
#include "include/block_for_rel_op_mesh_iso.inc"
      type is (VectorCurrentQ0)
#include "include/block_for_rel_op_mesh_iso.inc"
      end select
    end subroutine inside
  end subroutine SetTwoBodyRelOpMesh

  subroutine PrintTwoBodyRelOpMesh(this,iunit)
    class(TwoBodyRelOpMeshIso), intent(in) :: this
    integer, optional, intent(in) :: iunit
    integer :: chbra, chket, unt=6

    if(present(iunit)) unt = iunit
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        call this%MatCh(chbra,chket)%prt()
      end do
    end do
  end subroutine PrintTwoBodyRelOpMesh
end module TwoBodyRelOpsMeshIso
