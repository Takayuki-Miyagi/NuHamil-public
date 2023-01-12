module ThreeBodyJacOpsIso
  use omp_lib
  use ClassSys
  use OperatorDefinitions
  use ThreeBodyJacobiSpace
  use ThreeBodyJacOpsChanIso
  implicit none

  public :: ThreeBodyJacOpIso
  private
  private :: FinThreeBodyJacOp
  private :: InitThreeBodyJacOpFromString
  private :: InitThreeBodyJacOp
  private :: SetThreeBodyJacOp
  private :: set_three_body_jac_operator
  private :: GetMEFromIndex
  private :: GetMEFromQNs
  private :: GetMemory
  private :: ReducedToNonReduced
  private :: NonReducedToReduced
  private :: NonAntisymmetrize
  private :: PrintNonAntisymmetrizedME
  private :: ReadOperatorFile
  private :: WriteOperatorFile
  private :: GetOperatorFileName

  type, extends(OperatorDef) :: ThreeBodyJacOpIso
    type(ThreeBodyJacOpChanIso), allocatable :: MatCh(:,:)
    type(ThreeBodyJacIsoSpace), pointer :: ms
    logical :: is_init = .false.
    logical :: Antisymmetrized = .true.
  contains
    procedure :: InitThreeBodyJacOp
    procedure :: InitThreeBodyJacOpFromString
    procedure :: GetMEFromIndex
    procedure :: GetMEFromQNs
    procedure :: GetMemory
    procedure :: ReducedToNonReduced
    procedure :: NonReducedToReduced
    procedure :: NonAntisymmetrize
    procedure :: Antisymmetrize
    procedure :: PrintNonAntisymmetrizedME
    procedure :: ReadOperatorFile
    procedure :: WriteOperatorFile
    procedure :: GetOperatorFileName
    procedure :: SpinTensorDecomposition

    procedure :: CopyThreeBodyJacOpIso
    procedure :: SumThreeBodyJacOpIso
    procedure :: SubtractThreeBodyJacOpIso
    procedure :: ScaleThreeBodyJacOpIso
    generic :: assignment(=) => CopyThreeBodyJacOpIso
    generic :: operator(+) => SumThreeBodyJacOpIso
    generic :: operator(-) => SubtractThreeBodyJacOpIso
    generic :: operator(*) => ScaleThreeBodyJacOpIso

    generic :: init => InitThreeBodyJacOp, InitThreeBodyJacOpFromString
    procedure :: fin => FinThreeBodyJacOp
    procedure :: SetThreeBodyJacOp
    procedure :: ChangeModelSpace
    generic :: set => SetThreeBodyJacOp
    generic :: GetME => GetMEFromIndex, GetMEFromQNs
    generic :: writef => WriteOperatorFile
    generic :: readf => ReadOperatorFile
    generic :: GetFileName => GetOperatorFileName
  end type ThreeBodyJacOpIso
contains

  subroutine CopyThreeBodyJacOpIso(a, b)
    class(ThreeBodyJacOpIso), intent(inout) :: a
    class(ThreeBodyJacOpIso), intent(in) :: b
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
        a%MatCh(ichb,ichk)%OpName = b%MatCh(ichb,ichk)%OpName
        a%MatCh(ichb,ichk)%is_Zero = .false.
        a%MatCh(ichb,ichk)%is_init = b%MatCh(ichb,ichk)%is_init
        a%MatCh(ichb,ichk)%DMat = b%MatCh(ichb, ichk)%DMat
        a%MatCh(ichb,ichk)%jacobi_ch_bra => b%MatCh(ichb,ichk)%jacobi_ch_bra
        a%MatCh(ichb,ichk)%jacobi_ch_ket => b%MatCh(ichb,ichk)%jacobi_ch_ket
      end do
    end do
    a%is_init = b%is_init
    a%Antisymmetrized = b%Antisymmetrized
  end subroutine CopyThreeBodyJacOpIso

  function SumThreeBodyJacOpIso(a, b) result(c)
    use LinAlgLib
    type(ThreeBodyJacOpIso) :: c
    class(ThreeBodyJacOpIso), intent(in) :: a, b
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
  end function SumThreeBodyJacOpIso

  function SubtractThreeBodyJacOpIso(a, b) result(c)
    type(ThreeBodyJacOpIso) :: c
    class(ThreeBodyJacOpIso), intent(in) :: a, b
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
  end function SubtractThreeBodyJacOpIso

  function ScaleThreeBodyJacOpIso(a, b) result(c)
    type(ThreeBodyJacOpIso) :: c
    class(ThreeBodyJacOpIso), intent(in) :: a
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
  end function ScaleThreeBodyJacOpIso

  function GetMEFromIndex(this, chbra, chket, ibra, iket) result(me)
    class(ThreeBodyJacOpIso), intent(in) :: this
    integer, intent(in) :: chbra, chket, ibra, iket
    real(8) :: me
    me = 0.d0
    if(chbra * chket == 0) return
    me = this%MatCh(chbra,chket)%GetME(ibra,iket)
  end function GetMEFromIndex

  function GetMEFromQNs(this, Nbra, ibra, Jbra, Tbra, Nket, iket, Jket, Tket) result(me)
    class(ThreeBodyJacOpIso), intent(in) :: this
    integer, intent(in) :: Nbra, ibra, Jbra, Tbra, Nket, iket, Jket, Tket
    integer :: chbra, chket
    real(8) :: me
    me = 0.d0
    chbra = this%ms%GetIndex(Jbra, (-1)**Nbra, Tbra)
    chket = this%ms%GetIndex(Jket, (-1)**Nket, Tket)
    if(chbra*chket == 0) return
    me = this%MatCh(chbra,chket)%GetME(Nbra,ibra,Nket,iket)
  end function GetMEFromQNs

  subroutine FinThreeBodyJacOp(this)
    class(ThreeBodyJacOpIso), intent(inout) :: this
    integer :: chbra, chket
    if(.not. this%is_init) return
    do chbra = 1, size(this%MatCh,1)
      do chket = 1, size(this%MatCh,2)
        call this%MatCh(chbra,chket)%release()
      end do
    end do
    deallocate(this%MatCh)
    this%is_init = .false.
    call this%FinOperatorDef()
    this%ms => null()
    this%Antisymmetrized = .true.
  end subroutine FinThreeBodyJacOp

  subroutine InitThreeBodyJacOpFromString(this, ms, oprtr, Antisymmetrized)
    class(ThreeBodyJacOpIso), intent(out) :: this
    type(ThreeBodyJacIsoSpace), intent(in) :: ms
    type(str), intent(in) :: oprtr
    logical, intent(in), optional :: Antisymmetrized

    call this%InitOpDef(oprtr, .false.)
    if( present( Antisymmetrized )) this%Antisymmetrized = Antisymmetrized
    call this%init(ms,this%GetOpJ(),this%GetOpP(),this%GetOpT(), Antisymmetrized)
  end subroutine InitThreeBodyJacOpFromString

  subroutine InitThreeBodyJacOp(this, ms, jr, pr, tr, Antisymmetrized)
    use MPIFunction, only: myrank
    use MyLibrary, only: triag
    class(ThreeBodyJacOpIso), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), target, intent(in) :: ms
    integer, intent(in) :: jr, pr, tr
    logical, intent(in), optional :: Antisymmetrized
    integer :: chbra, jbra, pbra, tbra, nbra
    integer :: chket, jket, pket, tket, nket
    type(str) :: OpName

    if( present( Antisymmetrized )) this%Antisymmetrized = Antisymmetrized
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
        call this%MatCh(chbra,chket)%init(ms%jpt(chbra),ms%jpt(chket),this%GetOpName(), Antisymmetrized )
      end do
    end do
    this%is_init = .true.
    if(myrank == 0) then
      OpName = this%GetOpName()
      write(*,"(3a,f12.6,a)") "# ", OpName%val, " in Antisymmetrized 3-body basis: ", this%GetMemory(), " GB"
    end if
  end subroutine InitThreeBodyJacOp

  subroutine SetThreeBodyJacOp(this, params, subtract)
    use MPIFunction, only: myrank
    use NuHamilInput, only: InputParameters
    class(ThreeBodyJacOpIso), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: subtract

    if(.not. this%is_init) then
      write(*,*) "Initialize 'this' before calling SetThreeBodyJacOp"
      return
    end if

    if(myrank == 0) then
      write(*,'(4x, "Set three-body operator in jacobi basis ...")')
    end if
    call set_three_body_jac_operator(this, params, subtract)
  end subroutine SetThreeBodyJacOp

  subroutine set_three_body_jac_operator(this, params, subtract)
    use ClassSys, only: sys
    use MPIFunction
    use NuHamilInput, only: InputParameters
    use TwoBodyRelativeSpace
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpIso), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    logical, intent(in), optional :: subtract
    logical :: subtract_mode
    type(ThreeBodyJacIsoSpace), pointer :: ms
    integer, allocatable :: slranks(:)
    type(TwoBodyRelSpaceIsoHOBasis) :: rel_original
    type(TwoBodyRelOpIso) :: op_original
    type(str) :: OpName
    real(8) :: time
    integer :: ibra, iket, cnt, nloops
    integer, allocatable :: i2bra(:), i2ket(:)
#ifdef MPI
    integer :: nbra, nket
#endif
    type(sys) :: s

    if(s%isfile(this%GetFileName(params))) then
      call this%readf(params)
      return
    end if

    subtract_mode = .true.
    if( present(subtract) ) subtract_mode = subtract
    ms => this%ms
    OpName = this%GetOpName()
    select case(OpName%val)
    case("hamil", "Hamil","NNNint")

      call parent_child_procedure(inside_hamil, ms%GetNumberChannels(), slranks, time)
      call timer%Add(s%str("MPI parent-child, set three-body Jacobi op"), time)
#ifdef MPI
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(slranks(1), ms%GetNumberChannels(), mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      do iket = 1, ms%GetNumberChannels()
        nket = ms%jpt(iket)%GetNumberAStates()
        if(nket < 1) cycle
        call mpi_bcast(this%MatCh(iket,iket)%m(1,1), nket**2, mpi_real8, &
            &          slranks(iket), mpi_comm_world, ierr)
      end do
      call mpi_barrier(mpi_comm_world, ierr)
#endif
      deallocate(slranks)

    case default

      cnt = 0
      do ibra = 1, ms%GetNumberChannels()
        do iket = 1, ibra
          if( this%MatCh(ibra,iket)%is_zero ) cycle
          cnt = cnt + 1
        end do
      end do
      nloops = cnt
      allocate(i2bra(nloops))
      allocate(i2ket(nloops))
      cnt = 0
      do ibra = 1, ms%GetNumberChannels()
        do iket = 1, ibra
          if( this%MatCh(ibra,iket)%is_zero ) cycle
          cnt = cnt + 1
          i2bra(cnt) = ibra
          i2ket(cnt) = iket
        end do
      end do

      call rel_original%init(params%hw, params%N2max, params%Jmax2)
      call op_original%init( rel_original, rel_original, this%GetOpName() )
      call op_original%set()

      call parent_child_procedure(inside_ops, nloops, slranks, time)
      call timer%Add(s%str("MPI parent-child, set three-body Jacobi op"), time)
#ifdef MPI
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(slranks(1), nloops, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      do cnt = 1, nloops
        ibra = i2bra(cnt)
        iket = i2ket(cnt)
        nbra = ms%jpt(ibra)%GetNumberAStates()
        nket = ms%jpt(iket)%GetNumberAStates()
        if(nbra * nket < 1) cycle
        if(this%MatCh(ibra,iket)%is_zero) cycle
        call mpi_bcast(this%MatCh(ibra,iket)%m(1,1), nbra*nket, mpi_real8, &
            &          slranks(cnt), mpi_comm_world, ierr)
      end do
      call mpi_barrier(mpi_comm_world, ierr)
#endif
      call op_original%fin()
      call rel_original%fin()
      deallocate(i2bra, i2ket, slranks)
      do ibra = 1, ms%GetNumberChannels()
        do iket = 1, ibra
          if( this%MatCh(ibra,iket)%is_zero ) cycle
          this%MatCh(iket,ibra)%DMat = this%MatCh(ibra,iket)%DMat%T()
        end do
      end do
    end select

    if(params%fname_jacobi_op_for_trans2lab%val /= "none") then
      call this%writef(params)
      write(*,*) "End of calc."
      stop
    end if

  contains
    subroutine inside_hamil(loop)
      integer, intent(in) :: loop
      integer :: ch, Nmax_file
      type(str) :: fn, fn_v
      type(ThreeBodyJacIsoChan) :: jch
      type(ThreeBodyJacOpChanIso) :: ut ! dummy for file name
      integer :: iunit=101
      logical :: ex
      ch = loop
      if(ms%jpt(ch)%GetNumberAStates() < 1) return
      Nmax_file = GetRampNmax(ms%jpt(ch)%GetJ(), params%ramp)
      fn = jch%GetFileName(ms%jpt(ch)%GetJ(), &
          &  ms%jpt(ch)%GetParity(), ms%jpt(ch)%GetT(), Nmax_file, params%path_to_tmp_dir)
      ex = s%isfile(fn, s%str('set_three_body_jac_operator'))
      open(iunit, file=fn%val, form='unformatted', access='stream')
      call jch%readf(ms%GetFrequency(),iunit)
      close(iunit)
      fn_v = ut%GetFileName(jch,params%hw_target,this%GetOpName(),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      if(params%hw_conversion) then
        fn_v = ut%GetFileName(jch,params%hw,this%GetOpName(),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
      end if
      ex = s%isfile(fn_v, s%str('set_three_body_jac_operator'))
      open(iunit, file=fn_v%val, form='unformatted', access='stream')
      call this%MatCh(ch,ch)%readf(jch,&
          & iunit, Nmax_file, params%hw, params%hw_target, &
          & params%hw_conversion)
      close(iunit)
      call jch%fin()
    end subroutine inside_hamil

    subroutine inside_ops(loop)
      integer, intent(in) :: loop
      integer :: chbra, chket
      integer :: Nmax_file_bra, Nmax_file_ket
      type(str) :: fn_jacobi, fn_ut_bra, fn_ut_ket
      type(ThreeBodyJacIsoChan) :: jbra, jket
      type(ThreeBodyJacOpChanIso) :: ut ! dummy for file name
      logical :: ex
      integer :: iunit=101, iunit_bra=102, iunit_ket=103
      type(TwoBodyRelSpaceIsoHOBasis) :: relbra, relket
      type(TwoBodyRelChanIsoHOBasis), pointer :: rchbra, rchket, rchbra_, rchket_
      type(TwoBodyRelOpIso) :: op, opeff
      type(InputParameters) :: input
      !character(:), allocatable :: OpName
      type(HarmonicOscillator), pointer :: ho_bra, ho_ket
      integer :: ichbra, ichket, ichbra_, ichket_, bra, ket, bra_, ket_

      chbra = i2bra(loop)
      chket = i2ket(loop)
      if(ms%jpt(chbra)%GetNumberAStates() * ms%jpt(chket)%GetNumberAStates() < 1) return
      if(this%MatCh(chbra,chket)%is_zero) return
      Nmax_file_bra = GetRampNmax(ms%jpt(chbra)%GetJ(), params%ramp)
      fn_jacobi = jbra%GetFileName(ms%jpt(chbra)%GetJ(), &
          &  ms%jpt(chbra)%GetParity(), ms%jpt(chbra)%GetT(), Nmax_file_bra, params%path_to_tmp_dir)
      ex = s%isfile(fn_jacobi, s%str('In set_three_body_jac_operator'))
      open(iunit, file=fn_jacobi%val, form='unformatted', access='stream')
      call jbra%readf(ms%GetFrequency(),iunit)
      close(iunit)

      Nmax_file_ket = GetRampNmax(ms%jpt(chket)%GetJ(), params%ramp)
      fn_jacobi = jket%GetFileName(ms%jpt(chket)%GetJ(), &
          &  ms%jpt(chket)%GetParity(), ms%jpt(chket)%GetT(), Nmax_file_ket, params%path_to_tmp_dir)
      ex = s%isfile(fn_jacobi, s%str('In set_three_body_jac_operator'))
      open(iunit, file=fn_jacobi%val, form='unformatted', access='stream')
      call jket%readf(ms%GetFrequency(),iunit)
      close(iunit)

      call relbra%init(params%hw, Nmax_file_bra, params%Jmax2)
      call relket%init(params%hw, Nmax_file_ket, params%Jmax2)

      call op%init(relbra,relket,this%GetOpName() )
      call opeff%init(relbra,relket,this%GetOpName() )

      do ichbra = 1, relbra%GetNumberChannels()
        rchbra => relbra%GetChannel(ichbra)
        ichbra_ = rel_original%GetIndex( rchbra%GetJ(), rchbra%GetParity(), rchbra%GetT() )
        rchbra_ => rel_original%GetChannel(ichbra_)
        do ichket = 1, relket%GetNumberChannels()
          rchket => relket%GetChannel(ichket)
          ichket_ = rel_original%GetIndex( rchket%GetJ(), rchket%GetParity(), rchket%GetT() )
          rchket_ => rel_original%GetChannel(ichket_)
          if( op%MatCh(ichbra,ichket)%is_zero ) cycle
          do bra = 1, rchbra%GetNumberStates()
            ho_bra => rchbra%getp(bra)
            bra_ = rchbra_%GetIndex( ho_bra%GetN(), ho_bra%GetL(), ho_bra%GetS() )
            do ket = 1, rchket%GetNumberStates()
              ho_ket => rchket%getp(ket)
              ket_ = rchket_%GetIndex( ho_ket%GetN(), ho_ket%GetL(), ho_ket%GetS() )
              op%MatCh(ichbra,ichket)%m(bra,ket) = &
                  & op_original%MatCh(ichbra_,ichket_)%m(bra_,ket_)
              opeff%MatCh(ichbra,ichket)%m(bra,ket) = &
                  & op_original%MatCh(ichbra_,ichket_)%m(bra_,ket_)
            end do
          end do

        end do
      end do

      input = params
      call opeff%evolve(input, Nmax_file_bra, Nmax_file_ket)

      fn_ut_bra = ut%GetFileName(jbra,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      fn_ut_ket = ut%GetFileName(jket,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)

      if(params%hw_conversion) then
        fn_ut_bra = ut%GetFileName(jbra,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
        fn_ut_ket = ut%GetFileName(jket,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
      end if
      ex = s%isfile(fn_ut_bra, s%str('In set_three_body_jac_operator'))
      ex = s%isfile(fn_ut_ket, s%str('In set_three_body_jac_operator'))
      if(chbra == chket) then
        open(iunit_ket, file=fn_ut_ket%val, form='unformatted', access='stream')
        call this%MatCh(chbra,chket)%set(jbra, jket, iunit_ket, iunit_ket, op, opeff, &
            & Nmax_file_bra,Nmax_file_ket, params%hw, params%hw_target, params%hw_conversion, subtract_mode)
        close(iunit_ket)
        call jbra%fin()
        call jket%fin()
        call opeff%fin()
        call op%fin()
        call relbra%fin()
        call relket%fin()
        return
      end if

      open(iunit_bra, file=fn_ut_bra%val, form='unformatted', access='stream')
      open(iunit_ket, file=fn_ut_ket%val, form='unformatted', access='stream')
      call this%MatCh(chbra,chket)%set(jbra, jket, iunit_bra, iunit_ket, op, opeff, &
          & Nmax_file_bra,Nmax_file_ket, params%hw, params%hw_target, params%hw_conversion, subtract_mode)
      close(iunit_bra)
      close(iunit_ket)
      call jbra%fin()
      call jket%fin()
      call opeff%fin()
      call op%fin()
      call relbra%fin()
      call relket%fin()
    end subroutine inside_ops
  end subroutine set_three_body_jac_operator

  function GetMemory(this) result(r)
    class(ThreeBodyJacOpIso), intent(in) :: this
    real(8) :: r
    integer :: chbra, chket
    r = 0.d0
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        r = r + this%MatCh(chbra,chket)%GetMemory()
      end do
    end do
  end function GetMemory

  subroutine ReducedToNonReduced(this, convJ, convT)
    class(ThreeBodyJacOpIso), intent(inout) :: this
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
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%ReducedToNonReduced(convJ,convT)
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReduced(this, convJ, convT)
    class(ThreeBodyJacOpIso), intent(inout) :: this
    logical, intent(in), optional :: convJ, convT
    logical :: cJ = .true.
    logical :: cT = .true.
    integer :: chbra, chket
    if( present(convJ) ) cJ = convJ
    if( present(convJ) ) cT = convT
    if( cJ .and. this%GetOpJ()/=0 ) then
      write(*,*) "Error, J has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    if( cT .and. this%GetOpT()/=0 ) then
      write(*,*) "Error, T has to be 0 at ", __LINE__, " in ", __FILE__
    end if
    do chbra = 1, this%ms%GetNumberChannels()
      do chket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(chbra,chket)%is_zero) cycle
        call this%MatCh(chbra,chket)%NonReducedToReduced(convJ,convT)
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NonReducedToReduced

  subroutine NonAntisymmetrize( this )
    class(ThreeBodyJacOpIso), intent(inout) :: this
    integer :: ichbra, ichket
    if( .not. this%Antisymmetrized ) then
      write(*,*) "Oh, you already have non-antisymmetrized matrix elements, see at ", &
          & __LINE__, " in ", __FILE__
      return
    end if
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        call this%MatCh(ichbra,ichket)%NonAntisymmetrize()
      end do
    end do
    this%Antisymmetrized = .false.
  end subroutine NonAntisymmetrize

  subroutine Antisymmetrize( this )
    class(ThreeBodyJacOpIso), intent(inout) :: this
    integer :: ichbra, ichket
    if(this%Antisymmetrized ) then
      write(*,*) "Oh, you already have antisymmetrized matrix elements, see at ", &
          & __LINE__, " in ", __FILE__
      return
    end if
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        call this%MatCh(ichbra,ichket)%Antisymmetrize()
      end do
    end do
    this%Antisymmetrized = .true.
  end subroutine Antisymmetrize

  subroutine PrintNonAntisymmetrizedME( this, iunit )
    class(ThreeBodyJacOpIso), intent(inout) :: this
    integer, intent(in), optional :: iunit
    integer :: ichbra, ichket
    integer :: wunit
    if( this%Antisymmetrized ) then
      write(*,*) "Oh, you should have non-antisymmetrized ME first", __LINE__, " in ", __FILE__
      return
    end if
    wunit = 6
    if( present(iunit) ) wunit = iunit
    write(wunit,'(a,15x,a)') &
        & " n12,l12,s12,j12,t12, n3, l3, j3,  J,  T,n12,l12,s12,j12,t12, n3, l3, j3,  J,  T,",&
        & "ME"
    do ichbra = 1, this%ms%GetNumberChannels()
      do ichket = 1, this%ms%GetNumberChannels()
        if(this%MatCh(ichbra,ichket)%is_zero) cycle
        call this%MatCh(ichbra,ichket)%PrintNonAntisymmetrizedME(iunit)
      end do
    end do
  end subroutine PrintNonAntisymmetrizedME

  subroutine WriteOperatorFile(this, params)
    use MPIFunction, only: myrank
    use NuHamilInput, only: InputParameters
    use MyLibrary, only: triag
    class(ThreeBodyJacOpIso), intent(in) :: this
    type(InputParameters), intent(in) :: params
    type(str) :: fname
    type(ThreeBodyJacIsoSpace), pointer :: ms
    integer :: chbra, jbra, pbra, tbra, nbra
    integer :: chket, jket, pket, tket, nket
    integer :: wunit = 99

    if(myrank /= 0) return
    fname = this%GetFileName(params)
    write(*,'(3a)') "writing to ", trim(fname%val), " ..."
    open(wunit, file=fname%val, form="unformatted", access="stream", action="write")
    ms => this%ms
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
        if(triag(jbra, jket, 2*this%GetOpJ())) cycle
        if(pbra * this%GetOpP() /= pket) cycle
        if(triag(tbra, tket, 2*this%GetOpT())) cycle
        write(wunit) jbra, pbra, tbra, nbra, jket, pket, tket, nket
        if(nbra*nket>0) write(wunit) this%MatCh(chbra,chket)%m(:,:)
      end do
    end do
    close(wunit)
  end subroutine WriteOperatorFile

  subroutine ReadOperatorFile(this, params)
    use NuHamilInput, only: InputParameters
    use MyLibrary, only: triag
    class(ThreeBodyJacOpIso), intent(inout) :: this
    type(InputParameters), intent(in) :: params
    type(str) :: fname
    type(ThreeBodyJacIsoSpace), pointer :: ms
    integer :: chbra, jbra, pbra, tbra, nbra
    integer :: chket, jket, pket, tket, nket
    integer :: jbra_read, pbra_read, tbra_read, nbra_read
    integer :: jket_read, pket_read, tket_read, nket_read
    integer :: runit = 100

    fname = this%GetFileName(params)
    write(*,'(3a)') "reading from ", trim(fname%val), " ..."
    open(runit, file=fname%val, form="unformatted", access="stream", action="read")
    ms => this%ms
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
        if(triag(jbra, jket, 2*this%GetOpJ())) cycle
        if(pbra * this%GetOpP() /= pket) cycle
        if(triag(tbra, tket, 2*this%GetOpT())) cycle
        read(runit) jbra_read, pbra_read, tbra_read, nbra_read, jket_read, pket_read, tket_read, nket_read
        if(jbra /= jbra_read .or. pbra /= pbra_read .or. tbra /= tbra_read .or. nbra /= nbra_read .or. &
            & jket /= jket_read .or. pket /= pket_read .or. tket /= tket_read .or. nket /= nket_read ) then
            write(*,*) "Error occurs at", __LINE__, __FILE__
            stop
        end if
        if(nbra*nket>0) read(runit) this%MatCh(chbra,chket)%m(:,:)
      end do
    end do
    close(runit)
  end subroutine ReadOperatorFile

  function GetOperatorFileName(this, params) result(fname)
    use NuHamilInput, only: InputParameters
    use ClassSys, only: sys
    class(ThreeBodyJacOpIso), intent(in) :: this
    type(InputParameters), intent(in) :: params
    type(str) :: fname
    type(sys) :: s
    fname = params%path_to_tmp_dir + s%str("/ops/") + this%GetOpName()
    fname = fname + s%str("_") + params%fname_jacobi_op_for_trans2lab + s%str(".bin")
  end function GetOperatorFileName

  function ChangeModelSpace(this, ms) result(op)
    class(ThreeBodyJacOpIso), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), intent(in) :: ms
    type(ThreeBodyJacOpIso) :: op
    integer :: ichbra, ichket, ibra, iket
    integer :: ichbra_old, ichket_old, ibra_old, iket_old
    integer :: Ebra, Eket, abra, aket
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket, chbra_old, chket_old
    type(DMat) :: tmp

    call op%init(ms, this%GetOpName(), Antisymmetrized=.false.)
    do ichbra = 1, ms%GetNumberChannels()
      chbra => ms%GetChannel(ichbra)
      if(chbra%GetJ() > this%ms%GetJmax()) cycle
      ichbra_old = this%ms%GetIndex(chbra%GetJ(), chbra%GetParity(), chbra%GetT())
      chbra_old => this%ms%GetChannel(ichbra_old)
      do ichket = 1, ms%GetNumberChannels()
        chket => ms%GetChannel(ichket)
        if(chket%GetJ() > this%ms%GetJmax()) cycle
        ichket_old = this%ms%GetIndex(chket%GetJ(), chket%GetParity(), chket%GetT())
        chket_old => this%ms%GetChannel(ichket_old)
        if(op%MatCh(ichbra,ichket)%is_zero) cycle
        if(this%MatCh(ichbra_old,ichket_old)%is_zero) cycle
        if( chbra%GetNumberAStates() < 1 .or. chbra%GetNumberNAStates() < 1) return
        if( chket%GetNumberAStates() < 1 .or. chket%GetNumberNAStates() < 1) return
        if( chbra_old%GetNumberAStates() < 1 .or. chbra_old%GetNumberNAStates() < 1) return
        if( chket_old%GetNumberAStates() < 1 .or. chket_old%GetNumberNAStates() < 1) return
        call this%MatCh(ichbra_old,ichket_old)%NonAntisymmetrize()
        !$omp parallel
        !$omp do private(ibra, Ebra, abra, ibra_old, iket, Eket, aket, iket_old)
        do ibra = 1, chbra%GetNumberNAStates()
          Ebra = chbra%GetENAS(ibra)
          abra = chbra%GetAlphaNAS(ibra)
          if(Ebra > chbra_old%GetNmax()) cycle
          ibra_old = chbra_old%GetNASIndex(Ebra,abra)
          do iket = 1, chket%GetNumberNAStates()
            Eket = chket%GetENAS(iket)
            aket = chket%GetAlphaNAS(iket)
            if(Eket > chket_old%GetNmax()) cycle
            iket_old = chket_old%GetNASIndex(Eket,aket)
            op%MatCh(ichbra,ichket)%m(ibra,iket) = &
                & this%MatCh(ichbra_old,ichket_old)%m(ibra_old,iket_old)
          end do
        end do
        !$omp end do
        !$omp end parallel
        tmp = this%MatCh(ichbra_old,ichket_old)%DMat
        call this%MatCh(ichbra_old,ichket_old)%Antisymmetrize(tmp)
        tmp = op%MatCh(ichbra,ichket)%DMat
        call op%MatCh(ichbra,ichket)%Antisymmetrize(tmp)
      end do
    end do
  end function ChangeModelSpace

  function SpinTensorDecomposition(this, rank) result(op)
    use MyLibrary, only: triag
    use StoreCouplings
    class(ThreeBodyJacOpIso), intent(inout) :: this
    integer, intent(in) :: rank
    type(ThreeBodyJacIsoSpace), pointer :: ms
    type(ThreeBodyJacOpIso) :: op, op_tmp
    type(SixJsStore) :: sixj
    type(NineJsStore) :: ninej
    integer :: ich, Nmax
    integer, allocatable :: slranks(:)
    real(8) :: time
    type(sys) :: s
    type(str) :: opname
#ifdef MPI
    integer :: n_nas
#endif

    if( this%GetOpJ() /= 0 ) then
      write(*,*) "SpinTensorDecompsition is impolemented only for scalar operator"
      write(*,*) __LINE__, __FILE__
      opname = this%GetOpName()
      write(*,*) "OpName=", trim(opname%val), ", J=", this%GetOpJ(), ", Parity=", this%GetOpP(), ", T=", this%GetOpT()
      return
    end if
    if( rank > 3) then
      write(*,*) "Three-body spin-tensor decomposition is valid only up to rank=3"
      return
    end if

    op = this
    call op%NonAntisymmetrize()
    op_tmp = op
    Nmax = -1
    ms => this%ms
    do ich = 1, ms%GetNumberChannels()
      Nmax = max(Nmax, ms%jpt(ich)%GetNmax())
    end do
    call sixj%init(0,2*Nmax, .false., 1, 3, .true., 1, 3, .true.)
    call ninej%init(0,2*Nmax, .false., 0, 2, .false., 0,2*Nmax, .false., 1, 1, .true., j13dmin_in=0, j13dmax_in=2*Nmax)

    call parent_child_procedure(inside, ms%GetNumberChannels(), slranks, time)
    call timer%Add(s%str("3-body spin-tensor decomposition"), time)
#ifdef MPI
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(slranks(1), ms%GetNumberChannels(), mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    do ich = 1, ms%GetNumberChannels()
      n_nas = ms%jpt(ich)%GetNumberNAStates()
      if(n_nas < 1) cycle
      call mpi_bcast(op%MatCh(ich,ich)%m(1,1), n_nas**2, mpi_real8, slranks(ich), mpi_comm_world, ierr)
    end do
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(slranks)
    call sixj%fin()
    call ninej%fin()
    call op%Antisymmetrize()
  contains

    subroutine inside(loop)
      integer, intent(in) :: loop
      type(ThreeBodyJacIsoChan), pointer :: jac1, jac2
      integer :: J, JJ, ich2
      integer :: j12_2, j3_2, j45_2, j6_2
      integer :: ibra, iket, Lbra, Sbra, Lket, Sket, i_bra, i_ket
      type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
      real(8) :: me_rank, me

      jac1 => ms%GetChannel(loop)
      J = jac1%GetJ()
      if(jac1%GetNumberAStates() < 1) return

      !$omp parallel
      !$omp do private(ibra, bra, iket, ket, me_rank, Lbra, Lket, Sbra, Sket, me, JJ, jac2, ich2, &
      !$omp &  j12_2, j3_2, j45_2, j6_2, i_bra, i_ket)
      do ibra = 1, jac1%GetNumberNAStates()
        bra => jac1%GetNAS(ibra)
        do iket = 1, ibra
          ket => jac1%GetNAS(iket)

          me_rank = 0.d0
          do Lbra = abs(bra%GetL12()-bra%GetL3()), bra%GetL12()+bra%GetL3()
            do Lket = max(abs(ket%GetL12()-ket%GetL3()), abs(Lbra-rank)), min(ket%GetL12()+ket%GetL3(), Lbra+rank)
              do Sbra = max(abs(2*bra%GetS12()-1), abs(J-2*Lbra)), min(2*bra%GetS12()+1, J+2*Lbra), 2
                do Sket = max(abs(2*rank-Sbra), abs(2*ket%GetS12()-1), abs(J-2*Lket)), &
                      & min(2*rank+Sbra, 2*ket%GetS12()+1, J+2*Lket), 2

                  me = 0.d0
                  do JJ = max(abs(2*Lbra-Sbra),abs(2*Lket-Sket)), min(2*Lbra+Sbra,2*Lket+Sket), 2
                    ich2 = ms%GetIndex(JJ, jac1%GetParity(), jac1%GetT())
                    if(ich2==0) cycle
                    jac2 => ms%GetChannel(ich2)
                    if(2*(bra%GetN12()+bra%GetN3())+bra%GetL12()+bra%GetL3() > jac2%GetNmax()) cycle
                    if(2*(ket%GetN12()+ket%GetN3())+ket%GetL12()+ket%GetL3() > jac2%GetNmax()) cycle
                    do j12_2 = abs(bra%GetL12()-bra%GetS12()), bra%GetL12()+bra%GetS12()
                      do j3_2 = abs(2*bra%GetL3()-1), 2*bra%GetL3()+1, 2
                        if(triag(2*j12_2, j3_2, JJ)) cycle
                        do j45_2 = abs(ket%GetL12()-ket%GetS12()), ket%GetL12()+ket%GetS12()
                          do j6_2 = abs(2*ket%GetL3()-1), 2*ket%GetL3()+1, 2
                            if(triag(2*j45_2, j6_2, JJ)) cycle
                            i_bra = jac2%GetNASIndex(bra%GetN12(),bra%GetL12(),bra%GetS12(),j12_2,bra%GetT12(),&
                                &bra%GetN3(),bra%GetL3(),j3_2)
                            i_ket = jac2%GetNASIndex(ket%GetN12(),ket%GetL12(),ket%GetS12(),j45_2,ket%GetT12(),&
                                &ket%GetN3(),ket%GetL3(),j6_2)
                            if(i_bra * i_ket ==0) cycle
                            me = me + (-1.d0)**((JJ+J+2*Sbra)/2) * dble(JJ+1) * &
                                & sixj%Get(2*Lbra, Sbra, JJ, Sket, 2*Lket, 2*rank) * &
                                & ls_to_jj(bra%GetL12(), bra%GetS12(), j12_2, bra%GetL3(), j3_2, Lbra, Sbra, JJ) * &
                                & ls_to_jj(ket%GetL12(), ket%GetS12(), j45_2, ket%GetL3(), j6_2, Lket, Sket, JJ) * &
                                & op_tmp%MatCh(ich2,ich2)%m(i_bra,i_ket)
                          end do
                        end do
                      end do
                    end do
                  end do

                  me_rank = me_rank + me * dble(2*rank+1) * &
                      & sixj%Get(2*Lbra, Sbra, J, Sket, 2*Lket, 2*rank) * &
                      & ls_to_jj(bra%GetL12(), bra%GetS12(), bra%GetJ12(), bra%GetL3(), bra%GetJ3(), Lbra, Sbra, J) * &
                      & ls_to_jj(ket%GetL12(), ket%GetS12(), ket%GetJ12(), ket%GetL3(), ket%GetJ3(), Lket, Sket, J)
                end do
              end do

            end do
          end do

          op%MatCh(loop,loop)%m(ibra,iket) = me_rank
          op%MatCh(loop,loop)%m(iket,ibra) = me_rank
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine inside

    function ls_to_jj(l12, s12, j12, l3, j3, L, S, J) result(res)
      integer, intent(in) :: l12, s12, j12, l3, j3, L, S, J
      real(8) :: res
      res = sqrt(dble((2*j12+1)*(j3+1)*(2*L+1)*(S+1))) * &
          & ninej%get(2*l12, 2*s12, 2*j12, 2*l3, 1, j3, 2*L, S, J)
    end function ls_to_jj
  end function SpinTensorDecomposition
end module ThreeBodyJacOpsIso
