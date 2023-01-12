module ABodyJacOpsIso
  use omp_lib
  use LinAlgLib
  use ClassSys
  use OperatorDefinitions
  use ABodyJacSpaceIso
  implicit none

  public :: ABodyJacOpChanIso
  public :: ABodyJacOpIso

  private :: FinABodyJacOpChanIso
  private :: InitABodyJacOpChanIso
  private :: Set4BodyJacOpChanIso
  private :: SetABodyJacOpChanIso
  private :: WriteABodyJacScalarChanIso
  private :: ReadABodyJacScalarChanIso

  private :: set_4_body_from_3_body_scalar_isospin
  private :: set_4_body_from_3_body_tensor_isospin
  private :: set_A_body_from_A_1_body_scalar_isospin
  private :: set_A_body_from_A_1_body_tensor_isospin

  type, extends(DMat) :: ABodyJacOpChanIso
    type(ABodyJacIsoChan), pointer :: jacobi_ch_bra, jacobi_ch_ket
    type(str) :: OpName
    logical :: is_zero = .true.
    logical :: is_init = .false.
  contains
    procedure :: FinABodyJacOpChanIso
    procedure :: InitABodyJacOpChanIso
    procedure :: Set4BodyJacOpChanIso
    procedure :: SetABodyJacOpChanIso
    procedure :: WriteABodyJacScalarChanIso
    procedure :: ReadABodyJacScalarChanIso
    procedure :: CopyABodyJacOpChanIso
    procedure :: GetChanMEFromIndex
    procedure :: GetChanMEFromQNs

    generic :: release => FinABodyJacOpChanIso
    generic :: init => InitABodyJacOpChanIso
    generic :: set => Set4BodyJacOpChanIso, SetABodyJacOpChanIso
    generic :: assignment(=) => CopyABodyJacOpChanIso
    generic :: GetME => GetChanMEFromIndex, GetChanMEFromQNs
  end type ABodyJacOpChanIso

  type, extends(OperatorDef) :: ABodyJacOpIso
    type(ABodyJacOpChanIso), allocatable :: MatCh(:,:)
    type(ABodyJacIsoSpace), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: FinABodyJacOpIso
    procedure :: InitABodyJacOpIso
    procedure :: InitABodyJacOpIsoFromString
    procedure :: Set4BodyJacOpIso
    procedure :: SetABodyJacOpIso
    procedure :: CopyABodyJacOpIso
    procedure :: SetModelSpace
    procedure :: GetOpMEFromIndex
    procedure :: GetOpMEFromQNs

    generic :: fin => FinABodyJacOpIso
    generic :: init => InitABodyJacOpIso, InitABodyJacOpIsoFromString
    generic :: set => Set4BodyJacOpIso, SetABodyJacOpIso
    generic :: assignment(=) => CopyABodyJacOpIso
    generic :: GetME => GetOpMEFromIndex, GetOpMEFromQNs
  end type ABodyJacOpIso

  type(DMat), private :: cfp1, work, cfp2
  type(sys), private :: sy
contains

  function GetOpMEFromIndex(this, chbra, chket, ibra, iket) result(me)
    class(ABodyJacOpIso), intent(in) :: this
    integer, intent(in) :: chbra, chket, ibra, iket
    real(8) :: me
    me = 0.d0
    if(chbra * chket == 0) return
    me = this%MatCh(chbra,chket)%GetME(ibra,iket)
  end function GetOpMEFromIndex

  function GetOpMEFromQNs(this, Nbra, ibra, jbra, tbra, Nket, iket, jket, tket) result(me)
    class(ABodyJacOpIso), intent(in) :: this
    integer, intent(in) :: Nbra, ibra, jbra, tbra
    integer, intent(in) :: Nket, iket, jket, tket
    integer :: chbra, chket
    real(8) :: me
    me = 0.d0
    chbra = this%ms%GetIndex(jbra, (-1)**Nbra, tbra)
    chket = this%ms%GetIndex(jket, (-1)**Nket, tket)
    if(chbra * chket == 0) return
    me = this%MatCh(chbra,chket)%GetME(Nbra,ibra,Nket,iket)
  end function GetOpMEFromQNs

  function GetChanMEFromIndex(this, bra, ket) result(me)
    class(ABodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: bra, ket
    real(8) :: me
    me = 0.d0
    if(bra*ket == 0) return
    me = this%m(bra,ket)
  end function GetChanMEFromIndex

  function GetChanMEFromQNs(this, Nbra, ibra, Nket, iket) result(me)
    class(ABodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: Nbra, ibra, Nket, iket
    integer :: bra, ket
    real(8) :: me
    bra = this%jacobi_ch_bra%GetASIndex(Nbra,ibra)
    ket = this%jacobi_ch_ket%GetASIndex(Nket,iket)
    me = this%GetME(bra,ket)
  end function GetChanMEFromQNs

  subroutine FinABodyJacOpIso(this)
    class(ABodyJacOpIso), intent(inout) :: this
    integer :: chbra, chket
    if(.not. this%is_init) return
    do chbra = 1, size(this%MatCh,1)
      do chket = 1, size(this%MatCh,2)
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    this%ms => null()
    call this%FinOperatorDef()
    deallocate(this%MatCh)
    this%is_init = .false.
  end subroutine FinABodyJacOpIso

  subroutine InitABodyJacOpIsoFromString(this, ms, oprtr)
    class(ABodyJacOpIso), intent(out) :: this
    type(ABodyJacIsoSpace), intent(in) :: ms
    type(str), intent(in) :: oprtr

    call this%InitOpDef(oprtr, .false.)
    call this%init(ms,this%GetOpJ(),this%GetOpP(),this%GetOpT())
  end subroutine InitABodyJacOpIsoFromString

  subroutine InitABodyJacOpIso(this, ms, jr, pr, tr)
    use MyLibrary, only: triag
    class(ABodyJacOpIso), intent(inout) :: this
    type(ABodyJacIsoSpace), target, intent(in) :: ms
    integer, intent(in) :: jr, pr, tr
    integer :: chbra, jbra, pbra, tbra, nbra
    integer :: chket, jket, pket, tket, nket


    if(allocated(this%MatCh)) call this%fin()
    this%ms => ms
    call this%SetOpJ(jr)
    call this%SetOpP(pr)
    call this%SetOpT(tr)
    allocate(this%MatCh(ms%GetNumberChannels(),ms%GetNumberChannels()))
    do chbra = 1, ms%GetNumberChannels()
      jbra = ms%jpt(chbra)%GetJ()
      pbra = ms%jpt(chbra)%GetParity()
      tbra = ms%jpt(chbra)%GetT()
      nbra = ms%jpt(chbra)%GetNumberAStates()
      do chket = 1, ms%GetNumberChannels()
        jket = ms%jpt(chket)%GetJ()
        pket = ms%jpt(chket)%GetParity()
        tket = ms%jpt(chket)%GetT()
        nket = ms%jpt(chket)%GetNumberAStates()
        if(triag(jbra, jket, 2*jr)) cycle
        if(pbra * pr /= pket) cycle
        if(triag(tbra, tket, 2*tr)) cycle
        call this%MatCh(chbra,chket)%init(ms%jpt(chbra),ms%jpt(chket),this%GetOpName())
      end do
    end do
    this%is_init = .true.
  end subroutine InitABodyJacOpIso

  subroutine CopyABodyJacOpIso(a, b)
    class(ABodyJacOpIso), intent(inout) :: a
    type(ABodyJacOpIso), intent(in) :: b
    integer :: chbra, chket
    if(a%is_init) call a%fin()

    a%ms => b%ms
    call a%CopyOperatorDef(b%OperatorDef)
    allocate(a%MatCh(a%ms%GetNumberChannels(),a%ms%GetNumberChannels()))
    do chbra = 1, b%ms%GetNumberChannels()
      do chket = 1, b%ms%GetNumberChannels()
        if(b%MatCh(chbra,chket)%is_zero) cycle
        a%MatCh(chbra,chket) = b%MatCh(chbra,chket)
      end do
    end do
    a%is_init = .true.
  end subroutine CopyABodyJacOpIso

  subroutine SetModelSpace(this, ms)
    class(ABodyJacOpIso), intent(inout) :: this
    type(ABodyJacIsoSpace), target, intent(in) :: ms
    integer :: chbra, chket
    this%ms => ms
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
      if(this%MatCh(chbra,chket)%is_zero) cycle
      this%MatCh(chbra,chket)%jacobi_ch_bra => ms%jpt(chbra)
      this%MatCh(chbra,chket)%jacobi_ch_ket => ms%jpt(chket)
      end do
    end do
  end subroutine SetModelSpace


  subroutine Set4BodyJacOpIso(this, op)
    use ThreeBodyJacOpsIso
    class(ABodyJacOpIso), intent(inout) :: this
    type(ThreeBodyJacOpIso), intent(in) :: op
    integer :: chbra, chket

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling Set4BodyJacOpIso"
      return
    end if


    if(.not. op%is_init) then
      write(*,*) "Initialize 'op' before calling Set4BodyJacOpIso"
      return
    end if
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%set(op)
      end do
    end do
  end subroutine Set4BodyJacOpIso

  subroutine SetABodyJacOpIso(this, op)
    use ThreeBodyJacOpsIso
    class(ABodyJacOpIso), intent(inout) :: this
    type(ABodyJacOpIso), intent(in) :: op
    integer :: chbra, chket

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling Set4BodyJacOpIso"
      return
    end if


    if(.not. op%is_init) then
      write(*,*) "Initialize 'op' before calling Set4BodyJacOpIso"
      return
    end if
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%set(op)
      end do
    end do
  end subroutine SetABodyJacOpIso

  subroutine FinABodyJacOpChanIso(this)
    class(ABodyJacOpChanIso), intent(inout) :: this
    if(.not. this%is_init) return
    call this%fin()
    this%jacobi_ch_bra => null()
    this%jacobi_ch_ket => null()
    this%is_Zero = .true.
    this%is_init = .false.
  end subroutine FinABodyJacOpChanIso

  subroutine InitABodyJacOpChanIso(this, chbra, chket, oprtr)
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ABodyJacIsoChan), target, intent(in) :: chbra, chket
    type(str), intent(in) :: oprtr
    this%opname = oprtr
    this%jacobi_ch_bra => chbra
    this%jacobi_ch_ket => chket
    if(oprtr%val == "UT") call this%eye(chbra%GetNumberAStates())
    if(oprtr%val /= "UT") call this%zeros(chbra%GetNumberAStates(),chket%GetNumberAStates())
    this%is_zero = .false.
    this%is_init = .true.
  end subroutine InitABodyJacOpChanIso

  subroutine CopyABodyJacOpChanIso(a, b)
    class(ABodyJacOpChanIso), intent(inout) :: a
    type(ABodyJacOpChanIso), intent(in) :: b
    if(.not. a%is_zero) call a%release()
    a%opname = b%opname
    a%jacobi_ch_bra => b%jacobi_ch_bra
    a%jacobi_ch_ket => b%jacobi_ch_ket
    a%DMat = b%DMat
    a%is_zero = b%is_zero
    a%is_init = b%is_init
  end subroutine CopyABodyJacOpChanIso

  subroutine WriteABodyJacScalarChanIso(this, iunit)
    use Profiler, only: timer
    class(ABodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: iunit
    real(8) :: ti
    ti = omp_get_wtime()
    write(iunit) size(this%m,1)
    write(iunit) this%m
    call timer%Add(sy%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteABodyJacScalarChanIso

  subroutine ReadABodyJacScalarChanIso(this, jacl, iusc)
    use Profiler, only: timer
    class(ABodyJacOpChanIso), intent(inout) :: this
    class(ABodyJacIsoChan), intent(in) :: jacl
    integer, intent(in) :: iusc
    class(ABodyJacIsoChan), pointer :: jac
    type(ABodyJacOpChanIso) :: scalar
    integer :: n, n_cut
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling ReadThreeBodyJacScalarChannel"
      return
    end if

    jac => this%jacobi_ch_ket
    n_cut = jac%GetNumberAStates()
    call scalar%init(jacl,jacl,this%opname)
    ti = omp_get_wtime()
    read(iusc) n
    read(iusc) scalar%m
    call timer%Add(sy%str('Read from file'), omp_get_wtime() - ti)

    this%m = scalar%m(:n_cut, :n_cut)
    call scalar%release()
    this%is_Zero = .false.
  end subroutine ReadABodyJacScalarChanIso

  function GetFileNameABodyJacScalarChanIso(this, ms, oprtr, &
        &  genuine3bf, renorm, lambda, hw, cd, ce) result(f)
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ABodyJacIsoChan), intent(in) :: ms
    logical, intent(in) :: genuine3bf
    type(str), intent(in) :: oprtr, renorm
    real(8), intent(in) :: lambda, hw, cd, ce
    type(str) :: f
    type(str) :: pari, dir
    type(sys) :: s

    this%opname = oprtr
    dir = s%str('./A4/')
    f = dir + oprtr + s%str('_')
    if(oprtr%val == 'hamil' .or. oprtr%val == 'Hamil') f = dir
    if(genuine3bf) then
      f = f + s%str('NNN-full')
    else
      f = f + s%str('NNN-ind')
    end if

    if(ms%GetParity() == -1) pari = '-'
    if(ms%GetParity() ==  1) pari = '+'
    f = f + s%str('_j') + s%str(ms%GetJ()) + &
        & s%str('p') + pari + s%str('t') + s%str(ms%GetT()) + &
        & s%str('_Nmax') + s%str(ms%GetNmax())
    if(genuine3bf) f = f + s%str('_cd') + s%str(cd) + s%str('ce') + s%str(ce)
    f = f + s%str('_') + renorm
    if(renorm%val == 'srg') f = f + s%str(lambda)
    f = f + s%str('_hw') + s%str(hw) + s%str('.bin')
  end function GetFileNameABodyJacScalarChanIso

  subroutine Set4bodyJacOpChanIso(this, op)
    use ThreeBodyJacOpsIso
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacOpIso), intent(in) :: op

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling Set4BodyJacOpChanIso"
      return
    end if


    if(.not. op%is_init) then
      write(*,*) "Initialize 'op' before calling Set4BodyJacOpChanIso"
      return
    end if

    if(.not. op%reduced_me()) call set_4_body_from_3_body_scalar_isospin(this, op)
    if(op%reduced_me()) call set_4_body_from_3_body_tensor_isospin(this, op)
  end subroutine Set4bodyJacOpChanIso

  subroutine SetAbodyJacOpChanIso(this, op)
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ABodyJacOpIso), intent(in) :: op

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling SetABodyJacOpChanIso"
      return
    end if


    if(.not. op%is_init) then
      write(*,*) "Initialize 'op' before calling SetABodyJacOpChanIso"
      return
    end if
    if(.not. op%reduced_me()) call set_A_body_from_A_1_body_scalar_isospin(this, op)
    if(op%reduced_me()) call set_A_body_from_A_1_body_tensor_isospin(this, op)
  end subroutine SetAbodyJacOpChanIso

  subroutine set_4_body_from_3_body_scalar_isospin(this, op)
    use Profiler, only: timer
    use ThreeBodyJacOpsIso
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacOpIso), intent(in) :: op
    type(ABodyJacIsoChan), pointer :: jacobi_ch_A
    integer :: north, nphys, ibra, iket
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: ti, elem

    jacobi_ch_A => this%jacobi_ch_ket
    north = jacobi_ch_A%GetNumberAStates()
    nphys = jacobi_ch_A%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphys,nphys)
    cfp1 = jacobi_ch_A%GetCFPMat()

    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, elem) schedule(dynamic)
    do ibra = 1, nphys
      bra => jacobi_ch_A%GetNAS(ibra)
      do iket = 1, ibra
        ket => jacobi_ch_A%GetNAS(iket)

        elem = get_me_3(op, bra, ket)
        work%m(ibra,iket) = elem
        work%m(iket,ibra) = work%m(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%DMat = cfp1%T() * work * cfp1
    call work%fin()
    call cfp1%fin()
    call timer%add(sy%str('SetABodyScalarChannel'), omp_get_wtime() - ti)
  end subroutine set_4_body_from_3_body_scalar_isospin

  subroutine set_4_body_from_3_body_tensor_isospin(this, op)
    use Profiler, only: timer
    use MyLibrary, only: triag, sjs, hat
    use ThreeBodyJacOpsIso
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacOpIso), intent(in) :: op
    type(ABodyJacIsoChan), pointer :: chbra, chket
    integer :: north_bra, nphys_bra, north_ket, nphys_ket, ibra, iket
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: ti, elem, jfac, tfac

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    north_bra = chbra%GetNumberAStates()
    nphys_bra = chbra%GetNumberNAStates()
    north_ket = chket%GetNumberAStates()
    nphys_ket = chket%GetNumberNAStates()
    if(north_bra < 1 .or. nphys_bra < 1) return
    if(north_ket < 1 .or. nphys_ket < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphys_bra,nphys_ket)
    cfp1 = chbra%GetCFPMat()
    cfp2 = chket%GetCFPMat()

    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, jfac, tfac, elem) schedule(dynamic)
    do ibra = 1, nphys_bra
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphys_ket
        ket => chket%GetNAS(iket)
        jfac = (-1.d0) ** ((bra%Getja() + bra%Getjb() + chket%GetJ())/2 + op%GetOpJ()) * &
            & hat(chbra%GetJ()) * hat(chket%GetJ()) * &
            & sjs(bra%Getja(), chbra%GetJ(), bra%Getjb(), chket%GetJ(), ket%Getja(), 2*op%GetOpJ())

        tfac = (-1.d0) ** ((bra%Getta() + 1 + chket%GetT())/2 + op%GetOpT()) * &
            & hat(chbra%GetT()) * hat(chket%GetT()) * &
            & sjs(bra%Getta(), chbra%GetT(),  1, chket%GetT(), ket%Getta(), 2*op%GetOpT())
        elem = get_me_3(op, bra, ket)
        work%m(ibra,iket) = elem * jfac * tfac
      end do
    end do
    !$omp end do
    !$omp end parallel


    this%DMat = cfp1%T() * work * cfp2
    call work%fin()
    call cfp1%fin()
    call timer%add(sy%str('SetABodyTensorChannel'), omp_get_wtime() - ti)
  end subroutine set_4_body_from_3_body_tensor_isospin

  subroutine set_A_body_from_A_1_body_scalar_isospin(this, op)
    use Profiler, only: timer
    use ThreeBodyJacOpsIso
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ABodyJacOpIso), intent(in) :: op
    type(ABodyJacIsoChan), pointer :: jacobi_ch_A
    integer :: north, nphys, ibra, iket
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: ti, elem

    jacobi_ch_A => this%jacobi_ch_ket
    north = jacobi_ch_A%GetNumberAStates()
    nphys = jacobi_ch_A%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphys,nphys)
    cfp1 = jacobi_ch_A%GetCFPMat()

    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, elem) schedule(dynamic)
    do ibra = 1, nphys
      bra => jacobi_ch_A%GetNAS(ibra)
      do iket = 1, ibra
        ket => jacobi_ch_A%GetNAS(iket)

        elem = get_me_A1(op, bra, ket)
        work%m(ibra,iket) = elem
        work%m(iket,ibra) = work%m(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%DMat = cfp1%T() * work * cfp1
    call work%fin()
    call cfp1%fin()
    call timer%add(sy%str('SetABodyScalarChannel'), omp_get_wtime() - ti)
  end subroutine set_A_body_from_A_1_body_scalar_isospin

  subroutine set_A_body_from_A_1_body_tensor_isospin(this, op)
    use Profiler, only: timer
    use MyLibrary, only: triag, sjs, hat
    class(ABodyJacOpChanIso), intent(inout) :: this
    type(ABodyJacOpIso), intent(in) :: op
    type(ABodyJacIsoChan), pointer :: chbra, chket
    integer :: north_bra, nphys_bra, north_ket, nphys_ket, ibra, iket
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: ti, elem, jfac, tfac

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    north_bra = chbra%GetNumberAStates()
    nphys_bra = chbra%GetNumberNAStates()
    north_ket = chket%GetNumberAStates()
    nphys_ket = chket%GetNumberNAStates()
    if(north_bra < 1 .or. nphys_bra < 1) return
    if(north_ket < 1 .or. nphys_ket < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphys_bra,nphys_ket)
    cfp1 = chbra%GetCFPMat()
    cfp2 = chket%GetCFPMat()

    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, jfac, tfac, elem) schedule(dynamic)
    do ibra = 1, nphys_bra
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphys_ket
        ket => chket%GetNAS(iket)

        jfac = 1.d0
        tfac = 1.d0
        if(op%GetOpJ() /= 0) then
          jfac = (-1.d0) ** ((bra%Getja() + bra%Getjb() + chket%GetJ())/2 + op%GetOpJ()) * &
              & hat(chbra%GetJ()) * hat(chket%GetJ()) * &
              & sjs(bra%Getja(), chbra%GetJ(), bra%Getjb(), chket%GetJ(), ket%Getja(), 2*op%GetOpJ())
        end if

        if(op%GetOpT() /= 0) then
          tfac = (-1.d0) ** ((bra%Getta() + 1 + chket%GetT())/2 + op%GetOpT()) * &
              & hat(chbra%GetT()) * hat(chket%GetT()) * &
              & sjs(bra%Getta(), chbra%GetT(),  1, chket%GetT(), ket%Getta(), 2*op%GetOpT())
        end if
        elem = get_me_A1(op, bra, ket)
        work%m(ibra,iket) = elem * jfac * tfac
      end do
    end do
    !$omp end do
    !$omp end parallel

    this%DMat = cfp1%T() * work * cfp2
    call work%fin()
    call cfp1%fin()
    call timer%add(sy%str('SetABodyTensorChannel'), omp_get_wtime() - ti)
  end subroutine set_A_body_from_A_1_body_tensor_isospin

  function get_me_3(op, bra, ket) result(me)
    use MyLibrary, only: triag
    use ThreeBodyJacOpsIso
    type(ThreeBodyJacOpIso), intent(in) :: op
    type(NonAntisymmetrizedQNsA), intent(in) :: bra, ket
    real(8) :: me

    me = 0.d0
    if(bra%Getnb() /= ket%Getnb() ) return
    if(bra%Getlb() /= ket%Getlb() ) return
    if(bra%Getjb() /= ket%Getjb() ) return
    if(triag(bra%Getja(),ket%GetJa(),2*op%GetOpJ())) return
    if((-1)**(bra%GetNa()+ket%GetNa()) /= op%GetOpP()) return
    if(triag(bra%GetTa(),ket%GetTa(),2*op%GetOpT())) return
    me = op%GetME( bra%Getna(), bra%Getia(), bra%Getja(), bra%Getta(), &
        &          ket%Getna(), ket%Getia(), ket%Getja(), ket%Getta())
  end function get_me_3

  function get_me_A1(op, bra, ket) result(me)
    use MyLibrary, only: triag
    type(ABodyJacOpIso), intent(in) :: op
    type(NonAntisymmetrizedQNsA), intent(in) :: bra, ket
    real(8) :: me

    me = 0.d0
    if(bra%Getnb() /= ket%Getnb() ) return
    if(bra%Getlb() /= ket%Getlb() ) return
    if(bra%Getjb() /= ket%Getjb() ) return
    if(triag(bra%Getja(),ket%GetJa(),2*op%GetOpJ())) return
    if((-1)**(bra%GetNa()+ket%GetNa()) /= op%GetOpP()) return
    if(triag(bra%GetTa(),ket%GetTa(),2*op%GetOpT())) return
    me = op%GetME( bra%Getna(), bra%Getia(), bra%Getja(), bra%Getta(), &
        &          ket%Getna(), ket%Getia(), ket%Getja(), ket%Getta())
  end function get_me_A1
end module ABodyJacOpsIso
