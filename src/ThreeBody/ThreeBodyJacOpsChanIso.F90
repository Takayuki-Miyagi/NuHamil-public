module ThreeBodyJacOpsChanIso
  use MPIFunction
  use omp_lib
  use LinAlgLib
  use ClassSys
  use ThreeBodyJacChanIso
  implicit none

  public :: ThreeBodyJacOpChanIso

  private :: InitThreeBodyJacOpChan
  private :: FinThreeBodyJacOpChan
  private :: WriteThreeBodyJacScalarChan
  private :: SetThreeBodyJacOpChanFromTwoBody
  private :: SetThreeBodyJacOpChanEvolved
  private :: FrequencyConversion
  private :: ReadThreeBodyJacScalarChan
  private :: GetMEFromIndex
  private :: GetMEFromEi
  private :: set_three_body_from_two_body_scalar_isospin
  private :: set_three_body_from_two_body_scalar_lim_isospin
  private :: set_three_body_from_two_body_tensor_isospin
  private :: set_three_body_from_two_body_tensor_lim_isospin
  private :: set_three_body_evolved_scalar_isospin
  private :: set_three_body_evolved_tensor_isospin
  private :: frequency_conversion_scalar_isospin
  private :: frequency_conversion_tensor_isospin
  private :: GetFileNameThreeBodyJacScalarChan
  private :: init_overlaps
  private :: fin_overlaps
  private :: GetMemory
  private :: ReducedToNonReduced
  private :: NonReducedToReduced
  private :: PrintNonAntisymmetrizedME
  private :: NonAsym_to_Asym
  private :: Asym_to_NonAsym

  type, extends(DMat) :: ThreeBodyJacOpChanIso
    type(str) :: OpName
    type(ThreeBodyJacIsoChan), pointer :: jacobi_ch_bra, jacobi_ch_ket
    logical :: is_Zero = .true.
    logical :: is_init = .false.
    logical :: Antisymmetrized = .true.
  contains
    procedure :: InitThreeBodyJacOpChan
    procedure :: FinThreeBodyJacOpChan
    procedure :: SetThreeBodyJacOpChanFromOneBody
    procedure :: SetThreeBodyJacOpChanFromTwoBody
    procedure :: SetThreeBodyJacOpChanFromTwoBodyRelCM
    procedure :: SetThreeBodyJacOpChanEvolved
    procedure :: WriteThreeBodyJacScalarChan
    procedure :: ReadThreeBodyJacScalarChan
    procedure :: GetFileNameThreeBodyJacScalarChan
    procedure :: FrequencyConversion
    procedure :: GetMEFromIndex
    procedure :: GetMEFromEi
    procedure :: GetMemory
    procedure :: get_norm_op_chan
    procedure :: ReducedToNonReduced
    procedure :: NonReducedToReduced
    procedure :: PrintNonAntisymmetrizedME
    procedure :: NonAsym_to_Asym
    procedure :: from_NonAsym_to_Asym
    procedure :: Asym_to_NonAsym

    generic :: init => InitThreeBodyJacOpChan
    generic :: release => FinThreeBodyJacOpChan
    generic :: readf => ReadThreeBodyJacScalarChan
    generic :: writef => WriteThreeBodyJacScalarChan
    generic :: GetFileName => GetFileNameThreeBodyJacScalarChan
    generic :: set => SetThreeBodyJacOpChanFromTwoBody, SetThreeBodyJacOpChanFromTwoBodyRelCM, &
        &      SetThreeBodyJacOpChanEvolved, SetThreeBodyJacOpChanFromOneBody
    generic :: FreqConv => FrequencyConversion
    generic :: GetME => GetMEFromIndex, GetMEFromEi
    generic :: Antisymmetrize => NonAsym_to_Asym, from_NonAsym_to_Asym
    generic :: NonAntisymmetrize => Asym_to_NonAsym
  end type ThreeBodyJacOpChanIso

  type(DMat), private :: cfp1, work, cfp2
  type(ThreeBodyJacOpChanIso), private :: ovlap1, ovlap2
  real(8), private, allocatable :: ovlap_rad(:,:,:)
  real(8), private, allocatable :: Radhwin(:,:,:), Radhw(:,:,:)
  real(8), private, allocatable :: meshp(:), weight(:)
  real(8), private, parameter :: pmax = 25.d0
  integer, private, parameter :: nmesh = 500
  type(sys), private :: sy
contains

  function GetMEFromIndex(this, bra, ket) result(me)
    class(ThreeBodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: bra, ket
    real(8) :: me
    me = 0.d0
    if(bra*ket == 0) return
    me = this%m(bra,ket)
  end function GetMEFromIndex

  function GetMEFromEi(this, Nbra, ibra, Nket, iket) result(me)
    class(ThreeBodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: Nbra, ibra, Nket, iket
    real(8) :: me
    me = this%GetME( this%jacobi_ch_bra%GetASIndex(Nbra,ibra), this%jacobi_ch_ket%GetASIndex(Nket,iket) )
  end function GetMEFromEi

  ! isospin formalism methods
  subroutine InitThreeBodyJacOpChan(this, chbra, chket, oprtr, antisymmetrized)
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), target, intent(in) :: chbra, chket
    type(str), intent(in) :: oprtr
    logical, intent(in), optional :: antisymmetrized
    this%jacobi_ch_bra => chbra
    this%jacobi_ch_ket => chket
    this%OpName = oprtr
    if( present(antisymmetrized) ) this%Antisymmetrized=Antisymmetrized
    if(oprtr%val == 'UT') call this%eye(chbra%GetNumberAStates() )
    if(oprtr%val /= 'UT') call this%zeros( chbra%GetNumberAStates(), chket%GetNumberAStates() )
    if( .not. this%Antisymmetrized ) call this%zeros( chbra%GetNumberNAStates(), chket%GetNumberNAStates() )
    this%is_Zero = .false.
    this%is_init = .true.
  end subroutine InitThreeBodyJacOpChan

  subroutine FinThreeBodyJacOpChan(this)
    class(ThreeBodyJacOpChanIso), intent(inout) :: this

    this%is_Zero = .true.
    this%jacobi_ch_bra => null()
    this%jacobi_ch_ket => null()
    this%OpName = 'Unknown'
    call this%fin()
    this%is_init = .false.
    this%antisymmetrized = .true.
  end subroutine FinThreeBodyJacOpChan

  subroutine WriteThreeBodyJacScalarChan(this, iunit)
    use Profiler, only: timer
    class(ThreeBodyJacOpChanIso), intent(in) :: this
    integer, intent(in) :: iunit
    real(8) :: ti
    ti = omp_get_wtime()
    write(iunit) size(this%m,1)
    write(iunit) this%m
    call timer%Add(sy%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyJacScalarChan

  subroutine SetThreeBodyJacOpChanFromTwoBody(this,op)
    use TwoBodyRelOpsIso
    use TwoBodyRelativeSpace
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op

    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(.not. op%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(op%restricted()) then
      if(.not. op%reduced_me() ) call set_three_body_from_two_body_scalar_lim_isospin(this, op)
      if(op%reduced_me() ) call set_three_body_from_two_body_tensor_lim_isospin(this, op)
      return
    else
      !if(this%OpName%val == 'kinetic') then
      !  call set_three_body_kinetic_isospin(this)
      !  return
      !end if
      if(.not. op%reduced_me()) call set_three_body_from_two_body_scalar_isospin(this, op)
      if(op%reduced_me()) call set_three_body_from_two_body_tensor_isospin(this, op)
      return
    end if
  end subroutine SetThreeBodyJacOpChanFromTwoBody

  subroutine SetThreeBodyJacOpChanFromOneBody(this,op)
    use OneBodyLabOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(OneBodyLabOpIso), intent(in) :: op
    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if
    call set_three_body_from_one_body_tensor_isospin(this, op)
  end subroutine SetThreeBodyJacOpChanFromOneBody

  subroutine SetThreeBodyJacOpChanFromTwoBodyRelCM(this, op)
    use TwoBodyRelCMIsoOps
    use TwoBodyRelativeSpace
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelCMIsoOp), intent(in) :: op

    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(.not. op%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    call set_three_body_from_two_body_relcm_tensor_isospin(this, op)
  end subroutine SetThreeBodyJacOpChanFromTwoBodyRelCM

  subroutine SetThreeBodyJacOpChanEvolved(&
        & this, chbral, chketl, iutbra, iutket, &
        & op, opeff, Nmax_bra, Nmax_ket, hw_orig, hw_target, hw_conv, &
        & subtract)
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: chbral, chketl
    integer, intent(in) :: iutbra, iutket  ! unit number for three-body SRG unitary transformation
    type(TwoBodyRelOpIso), intent(in) :: op, opeff
    integer, intent(in) :: Nmax_bra, Nmax_ket
    real(8), intent(in) :: hw_orig, hw_target
    logical, intent(in) :: hw_conv
    logical, intent(in), optional :: subtract

    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(.not. op%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(.not. opeff%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    if(.not. op%reduced_me()) call set_three_body_evolved_scalar_isospin(&
        & this, chketl, iutbra, op, opeff, Nmax_ket, hw_orig, hw_target, hw_conv, subtract)
    if(op%reduced_me()) call set_three_body_evolved_tensor_isospin(&
        & this, chbral, chketl, iutbra, iutket, op, opeff, &
        & Nmax_bra, Nmax_ket, hw_orig, hw_target, hw_conv, subtract)
  end subroutine SetThreeBodyJacOpChanEvolved

  subroutine FrequencyConversion(this, Nmax_bra, Nmax_ket, hw_orig, hw_target)
    use Profiler, only: timer
    use OperatorDefinitions
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    integer, intent(in) :: Nmax_bra, Nmax_ket
    real(8), intent(in) :: hw_orig, hw_target
    type(OperatorDef) :: opdef
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    !write(*,'(a,f7.3,a,f7.3)') "# doing frequency conversion from ", hw_orig, &
    !    & " to ", hw_target
    call opdef%InitOpDef(this%opname, .false.)
    call init_overlaps(hw_orig, hw_target, max(Nmax_bra,Nmax_ket))
    if(opdef%GetOpJ() == 0 .and. opdef%GetOpT() == 0) then
      ti = omp_get_wtime()
      call frequency_conversion_scalar_isospin(this)
      call timer%Add(sy%str("frequency conversion scalar"), omp_get_wtime() - ti)
    end if

    if(opdef%GetOpJ() /= 0 .or. opdef%GetOpT() /= 0) then
      ti = omp_get_wtime()
      call frequency_conversion_tensor_isospin(this)
      call timer%Add(sy%str("frequency conversion tensor"), omp_get_wtime() - ti)
    end if

    call fin_overlaps()
  end subroutine FrequencyConversion

  subroutine ReadThreeBodyJacScalarChan(this, jacl, iusc, Nmax, hw_orig, hw_target, hw_conv)
    use Profiler, only: timer
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    class(ThreeBodyJacIsoChan), intent(in) :: jacl
    integer, intent(in) :: iusc, Nmax
    real(8), intent(in) :: hw_orig, hw_target
    logical, intent(in) :: hw_conv
    class(ThreeBodyJacIsoChan), pointer :: jac
    type(ThreeBodyJacOpChanIso) :: scalar
    integer :: n, n_cut
    real(8) :: ti

    if(.not. this%is_init) then
      write(*,*) "Error: at ", __LINE__, " in ", __FILE__
      return
    end if

    jac => this%jacobi_ch_ket
    n_cut = jac%GetNumberAStates()
    call scalar%init(jacl,jacl,this%opname)
    ti = omp_get_wtime()
    read(iusc) n
    read(iusc) scalar%m
    call timer%Add(sy%str('Read from file'), omp_get_wtime() - ti)

    if(scalar%opname%val == 'hamil' .or. scalar%opname%val == 'Hamil' .or. scalar%OpName%val=="NNNint") then
      if(hw_conv) call scalar%FreqConv(Nmax,Nmax,hw_orig,hw_target)
    end if

    this%m = scalar%m(:n_cut, :n_cut)
    call scalar%release()
    this%is_Zero = .false.
  end subroutine ReadThreeBodyJacScalarChan

  subroutine set_three_body_from_one_body_tensor_isospin(this, op)
    use SingleParticleState
    use OneBodyLabOpsIso
    use MyLibrary, only: triag, sjs
    use Profiler, only: timer
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(OneBodyLabOpIso), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    type(OrbitsIsospin), pointer :: isps
    real(8) :: ti, elem
    integer :: northb, nphysb, northk, nphysk, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    integer :: i1, i2

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    if(northb < 1 .or. nphysb < 1) return
    if(northk < 1 .or. nphysk < 1) return

    isps => op%sps
    ti = omp_get_wtime()
    call work%zeros(nphysb,nphysk)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, i1, i2, elem) schedule(dynamic)
    do ibra = 1, nphysb
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphysk
        ket => chket%GetNAS(iket)
        if( bra%GetN12() /= ket%GetN12() ) cycle
        if( bra%GetL12() /= ket%GetL12() ) cycle
        if( bra%GetS12() /= ket%GetS12() ) cycle
        if( bra%GetJ12() /= ket%GetJ12() ) cycle
        if( bra%GetT12() /= ket%GetT12() ) cycle
        i1 = isps%nlj2idx(bra%GetN3(), bra%GetL3(), bra%GetJ3())
        i2 = isps%nlj2idx(ket%GetN3(), ket%GetL3(), ket%GetJ3())
        elem = op%mat%m(i1,i2) * &
            & (-1.d0)**((chbra%GetJ()+ket%GetJ3())/2+ket%GetJ12()+op%GetOpJ()) * &
            & sqrt(dble((chbra%GetJ()+1)*(chket%GetJ()+1))) * &
            & sjs(bra%GetJ3(), chbra%GetJ(), 2*ket%GetJ12(), chket%GetJ(), ket%GetJ3(), 2*op%GetOpJ()) * &
            & (-1.d0)**((chbra%GetT()+          1)/2+ket%GetT12()+op%GetOpT()) * &
            & sqrt(dble((chbra%GetT()+1)*(chket%GetT()+1))) * &
            & sjs(1, chbra%GetT(), 2*ket%GetT12(), chket%GetT(), 1, 2*op%GetOpT())
        work%m(ibra,iket) = elem
      end do
    end do
    !$omp end do
    !$omp end parallel

    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
      this%DMat = this%DMat * 3.d0
    else
      this%DMat = work
      call work%fin()
    end if
    call timer%add(sy%str('SetThreeBodyJacOpChanOneScalar'), omp_get_wtime() - ti)
  end subroutine set_three_body_from_one_body_tensor_isospin

  subroutine set_three_body_from_two_body_scalar_lim_isospin(this, op)
    use MyLibrary, only: triag
    use Profiler, only: timer
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: jac
    real(8) :: ti, elem
    integer :: north, nphys, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphys,nphys)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, elem) schedule(dynamic)
    do ibra = 1, nphys
      bra => jac%GetNAS(ibra)
      do iket = 1, ibra
        ket => jac%GetNAS(iket)
        if( bra%GetN3() /= ket%GetN3() ) cycle
        if( bra%GetL3() /= ket%GetL3() ) cycle
        if( bra%GetJ3() /= ket%GetJ3() ) cycle
        if( triag(bra%GetJ12(), ket%GetJ12(), op%GetOpJ()) ) cycle
        if( (-1)**(bra%GetL12() + ket%GetL12()) * op%GetOpP() == -1) cycle
        if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) cycle
        elem = get_me12_from_qns(op, bra, ket, .true.)
        work%m(ibra,iket) = elem
        work%m(iket,ibra) = work%m(ibra, iket)
      end do
    end do
    !$omp end do
    !$omp end parallel

    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
      this%DMat = this%DMat * 3.d0
    else
      this%DMat = work
      call work%fin()
    end if
    call timer%add(sy%str('SetThreeBodyJacOpChanScalar'), omp_get_wtime() - ti)
  end subroutine set_three_body_from_two_body_scalar_lim_isospin

  subroutine set_three_body_from_two_body_scalar_isospin(this, op)
    use MyLibrary, only: triag
    use Profiler, only: timer
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: jac
    real(8) :: ti, elem
    integer :: north, nphys, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    ti = omp_get_wtime()
    call work%zeros(nphys,nphys)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, elem) schedule(dynamic)
    do ibra = 1, nphys
      bra => jac%GetNAS(ibra)
      do iket = 1, ibra
        ket => jac%GetNAS(iket)
        if( bra%GetN3() /= ket%GetN3() ) cycle
        if( bra%GetL3() /= ket%GetL3() ) cycle
        if( bra%GetJ3() /= ket%GetJ3() ) cycle
        if( triag(bra%GetJ12(), ket%GetJ12(), op%GetOpJ()) ) cycle
        if( (-1)**(bra%GetL12() + ket%GetL12()) * op%GetOpP() == -1) cycle
        if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) cycle
        elem = get_me12_from_qns(op, bra, ket, .false.)
        work%m(ibra,iket) = elem
        work%m(iket,ibra) = work%m(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel

    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
    else
      this%DMat = work
    end if

    select case(this%OpName%val )
    case('kinetic','Kinetic')
      this%DMat = this%DMat * 2.d0
    case default
      this%DMat = this%DMat * 3.d0
    end select
    call work%fin()
    call timer%add(sy%str('SetThreeBodyJacOpChanScalar'), omp_get_wtime() - ti)
  end subroutine set_three_body_from_two_body_scalar_isospin

  subroutine set_three_body_kinetic_isospin(this)
    use MyLibrary, only: triag
    use Profiler, only: timer
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), pointer :: jac
    real(8) :: ti, elem
    integer :: north, nphys, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

    write(*,*) '3-body kinetic operator in', __LINE__
    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    ti = omp_get_wtime()
    call work%zeros(nphys,nphys)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, elem) schedule(dynamic)
    do ibra = 1, nphys
      bra => jac%GetNAS(ibra)
      do iket = 1, ibra
        ket => jac%GetNAS(iket)
        if( bra%GetS12() /= ket%GetS12() ) cycle
        if( bra%GetJ12() /= ket%GetJ12() ) cycle
        if( bra%GetT12() /= ket%GetT12() ) cycle
        if( bra%GetJ3() /= ket%GetJ3() ) cycle
        elem = 0.d0
        if( bra%GetN3() == ket%GetN3() .and. bra%GetL3() == ket%GetL3() ) then
          elem = elem + kinetic(bra%GetN12(), bra%GetL12(), ket%GetN12(), ket%GetL12(), jac%GetFrequency())
        end if
        if( bra%GetN12() == ket%GetN12() .and. bra%GetL12() == ket%GetL12()) then
          elem = elem + kinetic( bra%GetN3(),  bra%GetL3(),  ket%GetN3(),  ket%GetL3(),  jac%GetFrequency())
        end if
        work%m(ibra,iket) = elem
        work%m(iket,ibra) = work%m(ibra,iket)
      end do
    end do
    !$omp end do
    !$omp end parallel

    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
    else
      this%DMat = work
      call work%fin()
    end if
    call timer%add(sy%str('SetThreeBodyJacOpChanScalar'), omp_get_wtime() - ti)
  contains
    function kinetic(n1, l1, n2, l2, hw) result(r)
      integer, intent(in) :: n1, l1, n2, l2
      real(8) :: hw, r
      r = 0.d0
      if(l1 /= l2) return
      if(abs(n1-n2) > 1) return
      if(n1 == n2) r = dble(2 * n1 + l1) + 1.5d0
      if(n1 == n2-1) r = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
      if(n1 == n2+1) r = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
      r = r * hw * 0.5d0
    end function kinetic
  end subroutine set_three_body_kinetic_isospin

  subroutine set_three_body_from_two_body_tensor_lim_isospin(this, op)
    use MyLibrary, only: triag, sjs, hat, dcg
    use Profiler, only: timer
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    real(8) :: ti, elem, jfac, tfac
    integer :: northb, nphysb, northk, nphysk, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    !integer :: zbra, zket

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    if(northb < 1 .or. nphysb < 1) return
    if(northk < 1 .or. nphysk < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphysb,nphysk)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, jfac, tfac, elem) schedule(dynamic)
    do ibra = 1, nphysb
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphysk
        ket => chket%GetNAS(iket)
        if( bra%GetN3() /= ket%GetN3() ) cycle
        if( bra%GetL3() /= ket%GetL3() ) cycle
        if( bra%GetJ3() /= ket%GetJ3() ) cycle
        if( triag(bra%GetJ12(), ket%GetJ12(), op%GetOpJ()) ) cycle
        if( (-1)**(bra%GetL12() + ket%GetL12()) * op%GetOpP() == -1) cycle
        if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) cycle
        jfac = (-1.d0) ** (bra%GetJ12() + (bra%GetJ3() + chket%GetJ())/2 + op%GetOpJ()) * &
            & hat(chbra%GetJ()) * hat(chket%GetJ()) * &
            & sjs(2*bra%GetJ12(), chbra%GetJ(), ket%GetJ3(), chket%GetJ(), 2*ket%GetJ12(), 2*op%GetOpJ())
        tfac = (-1.d0) ** (bra%GetT12() + (          1 + chket%GetT())/2 + op%GetOpT()) * &
            & hat(chbra%GetT()) * hat(chket%GetT()) * &
            & sjs(2*bra%GetT12(), chbra%GetT(),           1, chket%GetT(), 2*ket%GetT12(), 2*op%GetOpT())
        ! benchmarck with P.Gysbers pp <- nn
        !zbra=-1; zket=3
        !if( chbra%GetT()==3 .and. chket%GetT()==1 ) cycle
        !jfac = 1.d0 / sqrt(dble(2*ket%GetJ12() + 1))
        !tfac = dcg( 2*bra%GetT12(),-2*bra%GetT12(), 1, 1,  chbra%GetT(), zbra ) * &
        !    &  dcg( 2*ket%GetT12(), 2*ket%GetT12(), 1, 1,  chket%GetT(), zket ) / sqrt(5.d0)
        ! benchmark with P.Gysbers
        elem = get_me12_from_qns(op, bra, ket, .true.)
        work%m(ibra,iket) = elem * jfac * tfac
      end do
    end do
    !$omp end do
    !$omp end parallel
    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
      this%DMat = this%DMat * 3.d0
    else
      this%DMat = work
      call work%fin()
    end if
    call timer%add(sy%str('SetThreeBodyJacOpChanTensor'), omp_get_wtime() - ti)
  end subroutine set_three_body_from_two_body_tensor_lim_isospin

  subroutine set_three_body_from_two_body_tensor_isospin(this, op)
    use MyLibrary, only: triag, sjs, hat
    use Profiler, only: timer
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelOpIso), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    real(8) :: ti, elem, jfac, tfac
    integer :: northb, nphysb, northk, nphysk, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    if(northb < 1 .or. nphysb < 1) return
    if(northk < 1 .or. nphysk < 1) return

    ti = omp_get_wtime()
    call work%zeros(nphysb,nphysk)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, jfac, tfac, elem) schedule(dynamic)
    do ibra = 1, nphysb
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphysk
        ket => chket%GetNAS(iket)
        if( bra%GetN3() /= ket%GetN3() ) cycle
        if( bra%GetL3() /= ket%GetL3() ) cycle
        if( bra%GetJ3() /= ket%GetJ3() ) cycle
        if( triag(bra%GetJ12(), ket%GetJ12(), op%GetOpJ()) ) cycle
        if( (-1)**(bra%GetL12() + ket%GetL12()) * op%GetOpP() == -1) cycle
        if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) cycle

        jfac = (-1.d0) ** (bra%GetJ12() + (bra%GetJ3() + chket%GetJ())/2 + op%GetOpJ()) * &
            & hat(chbra%GetJ()) * hat(chket%GetJ()) * &
            & sjs(2*bra%GetJ12(), chbra%GetJ(), ket%GetJ3(), chket%GetJ(), 2*ket%GetJ12(), 2*op%GetOpJ())
        tfac = (-1.d0) ** (bra%GetT12() + ( 1 + chket%GetT())/2 + op%GetOpT()) * &
            & hat(chbra%GetT()) * hat(chket%GetT()) * &
            & sjs(2*bra%GetT12(), chbra%GetT(),           1, chket%GetT(), 2*ket%GetT12(), 2*op%GetOpT())
        elem = get_me12_from_qns(op, bra, ket, .false.)
        work%m(ibra,iket) = elem * jfac * tfac
      end do
    end do
    !$omp end do
    !$omp end parallel

    if(this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
      this%DMat = this%DMat * 3.d0
    else
      this%DMat = work
      call work%fin()
    end if
    call timer%add(sy%str('SetThreeBodyJacOpChanTensor'), omp_get_wtime() - ti)
  end subroutine set_three_body_from_two_body_tensor_isospin

  function get_me12_from_qns(op, bra, ket, mode_int) result(me)
    !
    ! mode_int: interaction mode, return 0 if "me" isn't in "op"
    !
    use MyLibrary, only: triag
    use TwoBodyRelOpsIso
    type(TwoBodyRelOpIso), intent(in) :: op
    type(NonAntisymmetrizedIsoQNs3), intent(in) :: bra, ket
    logical, intent(in) :: mode_int
    real(8) :: me
    type(TwoBodyRelSpaceIsoHOBasis), pointer :: rel_bra, rel_ket
    logical :: use_stored_op
    rel_bra => op%rel_sp_bra
    rel_ket => op%rel_sp_ket
    me = 0.d0
    if( bra%GetN3() /= ket%GetN3() ) return
    if( bra%GetL3() /= ket%GetL3() ) return
    if( bra%GetJ3() /= ket%GetJ3() ) return
    if( triag(bra%GetJ12(), ket%GetJ12(), op%GetOpJ()) ) return
    if( (-1)**(bra%GetL12() + ket%GetL12()) * op%GetOpP() == -1) return
    if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) return
    use_stored_op = .true.
    if( bra%GetJ12() > rel_bra%GetJmax()   ) use_stored_op = .false.
    if( ket%GetJ12() > rel_ket%GetJmax()   ) use_stored_op = .false.
    if( bra%GetL12() > rel_bra%GetJmax()+1 ) use_stored_op = .false.
    if( ket%GetL12() > rel_ket%GetJmax()+1 ) use_stored_op = .false.
    if( use_stored_op ) then
      me = op%Get2ME([bra%GetN12(), bra%GetL12(), bra%GetS12(), bra%GetJ12(), bra%GetT12()],&
          &          [ket%GetN12(), ket%GetL12(), ket%GetS12(), ket%GetJ12(), ket%GetT12()])
      return
    end if
    if( mode_int ) return
    me = CalcMERel( op%OperatorDef, &
        &          [bra%GetN12(), bra%GetL12(), bra%GetS12(), bra%GetJ12(), bra%GetT12()],&
        &          [ket%GetN12(), ket%GetL12(), ket%GetS12(), ket%GetJ12(), ket%GetT12()])
  end function get_me12_from_qns

  subroutine set_three_body_from_two_body_relcm_tensor_isospin(this, op)
    use MyLibrary, only: triag, sjs, hat, dcg
    use Profiler, only: timer
    use TwoBodyRelCMIsoOps
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(TwoBodyRelCMIsoOp), intent(in) :: op
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    real(8) :: ti, elem, tfac
    integer :: northb, nphysb, northk, nphysk, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    integer :: NNmax
    type(SixJsStore) :: sixj1, sixj2, sixj3
    type(TMbracketStore) :: tmbk

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    if(northb < 1 .or. nphysb < 1) return
    if(northk < 1 .or. nphysk < 1) return

    NNmax = max(chbra%GetNmax(),chket%GetNmax())
    call sixj1%init(0, 2*NNmax+2, .false., 0, 2*NNmax, .false., 1, 2*NNmax+1, .true.)
    call sixj2%init(0, 2*NNmax, .false., 0, 2*NNmax, .false., 1, 1, .true.)
    call sixj3%init(0, 2*NNmax+2, .false., chbra%GetJ(), chbra%GetJ(), .true., chket%GetJ(), chket%GetJ(), .true.)
    call tmbk%init(NNmax, 2.d0)

    ti = omp_get_wtime()
    call work%zeros(nphysb,nphysk)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, tfac, elem) schedule(dynamic)
    do ibra = 1, nphysb
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphysk
        ket => chket%GetNAS(iket)
        tfac = (-1.d0) ** (bra%GetT12() + (1 + chket%GetT())/2 + op%GetOpT()) * &
            & hat(chbra%GetT()) * hat(chket%GetT()) * &
            & sjs(2*bra%GetT12(), chbra%GetT(), 1, chket%GetT(), 2*ket%GetT12(), 2*op%GetOpT())
        elem = get_relcm_me12_from_qns(op, chbra%GetJ(), chket%GetJ(), bra, ket)
        work%m(ibra,iket) = elem * tfac
      end do
    end do
    !$omp end do
    !$omp end parallel
    if( this%Antisymmetrized ) then
      call this%Antisymmetrize(work)
      this%DMat = this%DMat * 3.d0
    else
      this%DMat = work
      call work%fin()
    end if

    call sixj1%fin()
    call sixj2%fin()
    call tmbk%fin()
    call timer%add(sy%str('SetThreeBodyJacOpChanCMTensor'), omp_get_wtime() - ti)
  contains
    function get_relcm_me12_from_qns(op, Jbra, Jket, bra, ket) result(me)
      use MyLibrary, only: triag, sjs, gmosh
      use TwoBodyRelCMSpaceIsoHO
      use TwoBodyRelCMIsoOps
      type(TwoBodyRelCMIsoOp), intent(in) :: op
      integer, intent(in) :: Jbra, Jket
      type(NonAntisymmetrizedIsoQNs3), intent(in) :: bra, ket
      real(8) :: me, me_tmp
      type(TwoBodyRelCMSpaceIsoHOBasis), pointer :: ms
      integer :: Npq_bra, Lpq_bra, Jpq_bra, N3_bra, Nmax_bra
      integer :: Npq_ket, Lpq_ket, Jpq_ket, N3_ket, Nmax_ket
      integer :: nr, lr, jr
      integer :: jrmin, jrmax, Lpq_bra_min, Lpq_bra_max, Lpq_ket_min, Lpq_ket_max
      integer :: Jpq_bra_min, Jpq_bra_max, Jpq_ket_min, Jpq_ket_max
      ms => op%ms
      me = 0.d0
      if( triag(bra%GetT12(), ket%GetT12(), op%GetOpT()) ) return
      if( bra%GetJ12() > ms%GetJRelmax()) return
      if( ket%GetJ12() > ms%GetJRelmax()) return
      N3_bra = 2*bra%GetN3() + bra%GetL3()
      N3_ket = 2*ket%GetN3() + ket%GetL3()
      Nmax_bra = 2*bra%GetN12() + bra%GetL12() + 2*bra%GetN3() + bra%GetL3()
      Nmax_ket = 2*ket%GetN12() + ket%GetL12() + 2*ket%GetN3() + ket%GetL3()
      do lr = 0, min(N3_bra, N3_ket)
        do nr = 0, (min(N3_bra, N3_ket)-lr)/2
          do jr = abs(2*lr-1), 2*lr+1, 2

            Lpq_bra_min = max(abs(bra%GetJ3()-jr)/2, abs(bra%GetL3()-lr))
            Lpq_bra_max = min(abs(bra%GetJ3()+jr)/2, abs(bra%GetL3()+lr), N3_bra-2*nr-lr)
            if(Lpq_bra_min > Lpq_bra_max) cycle
            Lpq_ket_min = max(abs(ket%GetJ3()-jr)/2, abs(ket%GetL3()-lr))
            Lpq_ket_max = min(abs(ket%GetJ3()+jr)/2, abs(ket%GetL3()+lr), N3_ket-2*nr-lr)
            if(Lpq_ket_min > Lpq_ket_max) cycle
            do Lpq_bra = Lpq_bra_min, Lpq_bra_max
              Npq_bra = (N3_bra - 2*nr - lr - Lpq_bra)/2
              do Lpq_ket = Lpq_ket_min, Lpq_ket_max
                Npq_ket = (N3_ket - 2*nr - lr - Lpq_ket)/2

                Jpq_bra_min = max(abs(Lpq_bra-bra%GetJ12()), abs(jr-Jbra)/2)
                Jpq_bra_max = min(abs(Lpq_bra+bra%GetJ12()), abs(jr+Jbra)/2, Nmax_bra+1)
                if(Jpq_bra_min > Jpq_bra_max) cycle
                Jpq_ket_min = max(abs(Lpq_ket-ket%GetJ12()), abs(jr-Jket)/2)
                Jpq_ket_max = min(abs(Lpq_ket+ket%GetJ12()), abs(jr+Jket)/2, Nmax_ket+1)
                if(Jpq_ket_min > Jpq_ket_max) cycle

                do Jpq_bra = Jpq_bra_min, Jpq_bra_max
                  do Jpq_ket = Jpq_ket_min, Jpq_ket_max
                    if(triag(Jpq_bra, Jpq_ket, op%GetOpJ())) cycle
                    me = me + (-1.d0)**(Lpq_bra+Lpq_ket-Jpq_ket+1+(bra%GetJ3()+Jbra+ket%GetJ3()+jr)/2+op%GetOpJ()) * &
                        & sqrt(dble( (2*Jpq_bra+1)*(bra%GetJ3()+1)*(2*bra%GetL3()+1)*(Jbra+1)*(Jket+1)* &
                        & (2*Jpq_ket+1)*(ket%GetJ3()+1)*(2*ket%GetL3()+1))) * dble(jr+1) * &
                        & sixj1%get(2*bra%GetJ12(), 2*Lpq_bra, 2*Jpq_bra, jr, Jbra, bra%GetJ3()) * &
                        & sixj2%get(2*Lpq_bra, 2*lr, 2*bra%GetL3(), 1, bra%GetJ3(), jr) * &
                        & sixj1%get(2*ket%GetJ12(), 2*Lpq_ket, 2*Jpq_ket, jr, Jket, ket%GetJ3()) * &
                        & sixj2%get(2*Lpq_ket, 2*lr, 2*ket%GetL3(), 1, ket%GetJ3(), jr) * &
                        & sixj3%get(2*Jpq_bra, Jbra, jr, Jket, 2*Jpq_ket, 2*op%GetOpJ()) * &
                        & tmbk%get(0,0,bra%GetN3(),bra%GetL3(),Npq_bra,Lpq_bra,nr,lr,bra%GetL3()) * &
                        & tmbk%get(0,0,ket%GetN3(),ket%GetL3(),Npq_ket,Lpq_ket,nr,lr,ket%GetL3()) * &
                        & op%Get2ME([Npq_bra,Lpq_bra,bra%GetN12(),bra%GetL12(),bra%GetS12(),bra%GetJ12(),Jpq_bra,bra%GetT12()], &
                        & [Npq_ket,Lpq_ket,ket%GetN12(),ket%GetL12(),ket%GetS12(),ket%GetJ12(),Jpq_ket,ket%GetT12()])
                  end do
                end do

              end do
            end do

          end do
        end do
      end do
    end function get_relcm_me12_from_qns
  end subroutine set_three_body_from_two_body_relcm_tensor_isospin

  subroutine set_three_body_evolved_scalar_isospin(this, jacl, iut, &
        & op, opeff, Nmax, hw_orig, hw_target, hw_conv, subtract)
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: jacl
    type(TwoBodyRelOpIso), intent(in) :: op, opeff
    integer, intent(in) :: iut, Nmax
    real(8), intent(in) :: hw_orig, hw_target
    logical, intent(in) :: hw_conv
    logical, intent(in), optional :: subtract
    type(ThreeBodyJacOpChanIso) :: UT, op3, op3sub
    type(ThreeBodyJacIsoChan), pointer :: jac
    integer :: north, nphys
    logical :: subtract_mode
    type(sys) :: s

    subtract_mode = .true.
    if( present(subtract) ) subtract_mode = subtract
    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    if(north < 1 .or. nphys < 1) return
    call UT%init(jacl, jacl, s%str('UT'))
    call UT%readf(jacl, iut, Nmax, hw_orig, hw_target, hw_conv)

    call op3%init(jacl,jacl,this%opname)
    call op3%set(op)

    call op3sub%init(jacl, jacl, this%opname)
    call op3sub%set(opeff)
    op3%DMat = UT%DMat%T() * op3%DMat * UT%DMat
    if( subtract_mode ) op3%DMat = op3%DMat - op3sub%DMat
    if(hw_conv) call op3%FreqConv(Nmax,Nmax,hw_orig,hw_target)
    this%m(:,:) = op3%m(:jac%GetNumberAStates(),:jac%GetNumberAStates())
    call UT%release()
    call op3%release()
    call op3sub%release()
  end subroutine set_three_body_evolved_scalar_isospin

  subroutine set_three_body_evolved_tensor_isospin(&
        & this, chbral, chketl, iutbra, iutket, &
        & op, opeff, Nmax_bra, Nmax_ket, hw_orig, hw_target, hw_conv, &
        & subtract)
    use TwoBodyRelOpsIso
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: chbral, chketl
    type(TwoBodyRelOpIso), intent(in) :: op, opeff
    integer, intent(in) :: iutbra, iutket, Nmax_bra, Nmax_ket
    real(8), intent(in) :: hw_orig, hw_target
    logical, intent(in) :: hw_conv
    logical, intent(in), optional :: subtract
    logical :: subtract_mode
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    type(ThreeBodyJacOpChanIso) :: UTb, UTk, op3, op3sub
    integer :: northb, nphysb, northk, nphysk
    type(sys) :: s

    subtract_mode = .true.
    if( present(subtract) ) subtract_mode = subtract
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    if(northb < 1 .or. nphysb < 1) return
    if(northk < 1 .or. nphysk < 1) return

    call UTb%init(chbral, chbral, s%str('UT'))
    call UTb%readf(chbral, iutbra, Nmax_bra, hw_orig, hw_target, hw_conv)

    if( iutbra == iutket ) then
      UTk = UTb
    else
      call UTk%init(chketl, chketl, s%str('UT'))
      call UTk%readf(chketl, iutket, Nmax_ket, hw_orig, hw_target, hw_conv)
    end if

    call op3%init(chbral, chketl, this%opname)
    call op3%set(op)

    call op3sub%init(chbral, chketl, this%opname)
    call op3sub%set(opeff)
    op3%DMat = UTb%DMat%T() * op3%DMat * UTk%DMat
    if( subtract_mode ) op3%DMat = op3%DMat - op3sub%DMat
    call UTb%release()
    call UTk%release()
    call op3sub%release()
    if(hw_conv) call op3%FreqConv(Nmax_bra,Nmax_ket,hw_orig,hw_target)
    this%m(:,:) = op3%m(:chbra%GetNumberAStates(), :chket%GetNumberAStates())
    call op3%release()
  end subroutine set_three_body_evolved_tensor_isospin

  subroutine frequency_conversion_scalar_isospin(this)
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    class(ThreeBodyJacIsoChan),  pointer :: jac
    integer :: nphys, north, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    real(8) :: mem
    type(sys) :: s

    jac => this%jacobi_ch_ket
    north = jac%GetNumberAStates()
    nphys = jac%GetNumberNAStates()
    mem = (dble(nphys)+dble(north))**2
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,f8.4,a)') "# Scalar frequency conversion: myrank=", myrank, &
        & ", J=",jac%GetJ(), "/2, P=",jac%GetParity(), &
        & ", T=",jac%GetT(), "/2, ", &
        & 8.d0*mem/1024.d0**3, " GB will be needed at least..."
    call work%zeros(nphys,nphys)

    !$omp parallel
    !$omp do private(ibra, bra, iket, ket)
    do ibra = 1, nphys
      bra => jac%GetNAS(ibra)
      do iket = 1, nphys
        ket => jac%GetNAS(iket)
        if(bra%GetL12() /= ket%GetL12()) cycle
        if(bra%GetS12() /= ket%GetS12()) cycle
        if(bra%GetJ12() /= ket%GetJ12()) cycle
        if(bra%GetT12() /= ket%GetT12()) cycle
        if(bra%GetL3()  /= ket%GetL3() ) cycle
        if(bra%GetJ3()  /= ket%GetJ3() ) cycle
        work%m(ibra,iket) = ovlap_rad(bra%GetN12(),ket%GetN12(),ket%GetL12()) *&
            &               ovlap_rad(bra%GetN3( ),ket%GetN3( ),ket%GetL3( ))
      end do
    end do
    !$omp end do
    !$omp end parallel

    call ovlap1%init(jac,jac,s%str("overlap"))
    call ovlap1%Antisymmetrize(work)
    call work%fin()
    this%DMat = ovlap1%DMat%T() * this%DMat * ovlap1%DMat
    call ovlap1%release()
  end subroutine frequency_conversion_scalar_isospin

  subroutine frequency_conversion_tensor_isospin(this)
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    class(ThreeBodyJacIsoChan), pointer :: chbra,chket
    integer :: nphysb, northb, nphysk, northk, ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    type(sys) :: s
    real(8) :: mem

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    northb = chbra%GetNumberAStates()
    nphysb = chbra%GetNumberNAStates()
    northk = chket%GetNumberAStates()
    nphysk = chket%GetNumberNAStates()
    mem = max( dble(nphysb)**2+2.d0*dble(nphysb)*dble(northb), dble(nphysk)**2+2.d0*dble(nphysk)*dble(northk) ) + &
        & dble(northb)**2 + dble(northk)**2
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f8.4,a)') "# Tensor frequency conversion: myrank=", myrank, &
        & ", Jbra=",chbra%GetJ(), "/2, Pbra=",chket%GetParity(), ", Tbra=",chbra%GetT(), &
        & "/2, Jket=",chket%GetJ(), "/2, Pket=",chket%GetParity(), ", Tket=",chket%GetT(), "/2, ", &
        & 8.d0*mem/1024.d0**3, " GB will be needed at least..."
    call work%zeros(nphysb,nphysb)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket)
    do ibra = 1, nphysb
      bra => chbra%GetNAS(ibra)
      do iket = 1, nphysb
        ket => chbra%GetNAS(iket)
        if(bra%GetL12() /= ket%GetL12()) cycle
        if(bra%GetS12() /= ket%GetS12()) cycle
        if(bra%GetJ12() /= ket%GetJ12()) cycle
        if(bra%GetT12() /= ket%GetT12()) cycle
        if(bra%GetL3()  /= ket%GetL3() ) cycle
        if(bra%GetJ3()  /= ket%GetJ3() ) cycle
        work%m(ibra,iket) = ovlap_rad(bra%GetN12(),ket%GetN12(),ket%GetL12()) *&
            &               ovlap_rad(bra%GetN3( ),ket%GetN3( ),ket%GetL3( ))
      end do
    end do
    !$omp end do
    !$omp end parallel
    call ovlap1%init(chbra,chbra,s%str("overlap"))
    call ovlap1%Antisymmetrize(work)

    call work%zeros(nphysk,nphysk)
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket)
    do ibra = 1, nphysk
      bra => chket%GetNAS(ibra)
      do iket = 1, nphysk
        ket => chket%GetNAS(iket)
        if(bra%GetL12() /= ket%GetL12()) cycle
        if(bra%GetS12() /= ket%GetS12()) cycle
        if(bra%GetJ12() /= ket%GetJ12()) cycle
        if(bra%GetT12() /= ket%GetT12()) cycle
        if(bra%GetL3()  /= ket%GetL3() ) cycle
        if(bra%GetJ3()  /= ket%GetJ3() ) cycle
        work%m(ibra,iket) = ovlap_rad(bra%GetN12(),ket%GetN12(),ket%GetL12()) *&
            &               ovlap_rad(bra%GetN3( ),ket%GetN3( ),ket%GetL3( ))
      end do
    end do
    !$omp end do
    !$omp end parallel
    call ovlap2%init(chket,chket,s%str("overlap"))
    call ovlap2%Antisymmetrize(work)
    call work%fin()

    this%DMat = ovlap1%DMat%T() * this%DMat * ovlap2%DMat
    call ovlap1%release()
    call ovlap2%release()
  end subroutine frequency_conversion_tensor_isospin

  subroutine NonAsym_to_Asym( this, Mat_NonAsym )
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(DMat), intent(inout) :: Mat_NonAsym
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    type(DMat) :: tmp
    integer :: k,l,m,n

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    k = chbra%GetNumberAStates()
    l = chbra%GetNumberNAStates()
    m = chket%GetNumberNAStates()
    n = chket%GetNumberAStates()
    call this%DMat%fin() ! release once to save RAM; it can save ~10 GB in a maximum case.
    if(chbra%GetJ() == chket%GetJ() .and. &
        & chbra%GetParity() == chket%GetParity() .and. &
        & chbra%GetT() == chket%GetT() ) then
      ! scalar
      cfp1 = chket%GetCFPMat()
      call tmp%ini( l, n )
      call dgemm('n','n', l, n, m, 1.d0, Mat_NonAsym%m, l, cfp1%m, m, 0.d0, tmp%m, l)
      call Mat_NonAsym%fin()

      call this%DMat%ini( k,n )
      call dgemm('t','n', k, n, l, 1.d0, cfp1%m, l, tmp%m, l, 0.d0, this%DMat%m, k)
      call tmp%fin()
      call cfp1%fin()
    else
      ! tensor
      cfp1 = chket%GetCFPMat()
      call tmp%ini( l, n )
      call dgemm('n','n', l, n, m, 1.d0, Mat_NonAsym%m, l, cfp1%m, m, 0.d0, tmp%m, l)
      call Mat_NonAsym%fin()

      call cfp1%fin()
      cfp1 = chbra%GetCFPMat()
      call this%DMat%ini( k,n )
      call dgemm('t','n', k, n, l, 1.d0, cfp1%m, l, tmp%m, l, 0.d0, this%DMat%m, k)
      call tmp%fin()
      call cfp1%fin()
    end if
  end subroutine NonAsym_to_Asym

  function GetFileNameThreeBodyJacScalarChan(this, ms, hw, oprtr, &
        &  genuine3bf, Regulator, RegulatorPower, renorm, &
        &  lambda_3nf_local, lambda_3nf_nonlocal, lambda, c1, c3, c4, cd, ce, &
        &  J3max_initial_3nf, path_to_dir) result(f)
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), intent(in) :: ms
    logical, intent(in) :: genuine3bf
    integer, intent(in) :: RegulatorPower, J3max_initial_3nf
    type(str), intent(in) :: oprtr, renorm, Regulator, path_to_dir
    real(8), intent(in) :: lambda, lambda_3nf_local, lambda_3nf_nonlocal
    real(8), intent(in) :: hw, c1, c3, c4, cd, ce
    type(str) :: f, dir, pari
    type(sys) :: s

    this%opname = oprtr
    dir = path_to_dir + s%str('/ops')
    call s%mkdir(dir%val)
    select case(oprtr%val)
    case('hamil',"Hamil","NNNint")
      f = dir + s%str("/")
    case default
      f = dir + s%str("/") + oprtr +  s%str('_')
    end select

    if(genuine3bf .and. ms%GetJ() <= J3max_initial_3nf) then
      f = f + s%str('NNN-full')
    else
      f = f + s%str('NNN-ind')
    end if

    if(ms%GetParity() == -1) pari = '-'
    if(ms%GetParity() ==  1) pari = '+'
    f = f + s%str('_j') + s%str(ms%GetJ()) + &
        & s%str('p') + pari + s%str('t') + s%str(ms%GetT()) + &
        & s%str('_Nmax') + s%str(ms%GetNmax())
    if(genuine3bf .and. ms%GetJ() <= J3max_initial_3nf) then
      f = f + s%str("_") + Regulator + s%str(RegulatorPower)
      select case(Regulator%val)
      case("LNL","lnl",'local-non-local',"Local-Non-Local")
        f = f + s%str("_") + s%str(lambda_3nf_local) + s%str("_") + s%str(lambda_3nf_nonlocal)
      case('Local', 'local')
        f = f + s%str("_") + s%str(lambda_3nf_local)
      case('NonLocal', 'nonlocal', 'Nonlocal')
        f = f + s%str("_") + s%str(lambda_3nf_nonlocal)
      end select
      f = f + s%str("_c1_") + s%str(c1) + s%str("_c3_") + &
          & s%str(c3) + s%str("_c4_") + s%str(c4) + &
          & s%str('_cD_') + s%str(cd) + s%str('_cE_') + s%str(ce)
    end if
    f = f + s%str('_') + renorm
    if(renorm%val == 'srg') f = f + s%str(lambda)
    f = f + s%str('_hw') + s%str(hw) + s%str('.bin')
  end function GetFileNameThreeBodyJacScalarChan

  subroutine init_overlaps(hw_orig, hw_target, Nmax)
    use MyLibrary, only: m_proton, m_neutron, hc, gauss_legendre, ho_radial_wf_norm
    real(8), intent(in) :: hw_orig, hw_target
    integer, intent(in) :: Nmax
    integer :: i, n, l, n1, n2
    real(8) :: r, a, a_target, a_orig, wave_orig, wave_target
    allocate(Radhwin(nmesh, 0:nmax/2, 0:nmax))
    allocate(Radhw(nmesh, 0:nmax/2, 0:nmax))
    allocate(meshp(nmesh), weight(nmesh))
    allocate(ovlap_rad(0:nmax/2, 0:nmax/2, 0:nmax))
    call gauss_legendre(0.d0, pmax, meshp, weight, nmesh)
    !$omp parallel
    !$omp do private(i, n, l, r, a_orig, a_target, wave_orig, wave_target)
    do i = 1, nmesh
      do n = 0, nmax/2
        do l = 0, nmax - 2 * n
          r = meshp(i)
          a_target = 2.d0 * (m_proton * m_neutron) / (m_proton + m_neutron) / hc ** 2 * hw_target
          a_orig = 2.d0 * (m_proton * m_neutron) / (m_proton + m_neutron) / hc ** 2 * hw_orig
          Radhw(i, n, l) = ho_radial_wf_norm(n,l,a_target,r)
          Radhwin(i, n, l) = ho_radial_wf_norm(n,l,a_orig,r)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    do n1 = 0, Nmax/2
      do n2 = 0, Nmax/2
        do l = 0, Nmax
          if(2 * n1 + l > nmax) cycle
          if(2 * n2 + l > nmax) cycle
          a = 0.d0
          do i = 1, nmesh
            a = a + Radhwin(i, n1, l) * Radhw(i, n2, l) * weight(i)
          end do
          ovlap_rad(n1, n2, l) = a
        end do
      end do
    end do
  end subroutine init_overlaps

  subroutine fin_overlaps()
    deallocate(Radhwin, Radhw)
    deallocate(ovlap_rad, meshp, weight)
  end subroutine fin_overlaps

  function get_norm_op_chan(this, mode) result(fmat)
    class(ThreeBodyJacOpChanIso), intent(in) :: this
    type(str), optional, intent(in) :: mode
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    type(NmaxChannelIso3), pointer :: nchbra, nchket
    type(DMat) :: fmat
    real(8), allocatable :: v(:,:)
    integer :: Nbra, Nket, n_states_bra, n_states_ket
    integer :: bra, ket

    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    call fmat%zeros(chbra%GetNmax()+1, chket%GetNmax()+1)
    do Nbra = 0, chbra%GetNmax()
      do Nket = 0, chket%GetNmax()
        nchbra => chbra%GetNmaxChannel(Nbra)
        nchket => chbra%GetNmaxChannel(Nket)
        n_states_bra = nchbra%GetNumberAStates()
        n_states_ket = nchket%GetNumberAStates()
        if(n_states_bra * n_states_ket == 0) cycle
        allocate(v( n_states_bra, n_states_ket ))

        do bra = 1, n_states_bra
          do ket = 1, n_states_ket
            v( bra, ket ) = this%m(chbra%GetASIndex(Nbra,bra), chket%GetASIndex(Nket,ket))
          end do
        end do
        if(present(mode)) then
          if(mode%val == "ave" .or. mode%val == "average") then
            fmat%m(Nbra+1, Nket+1) = get_frobenius_norm(v) / &
                & sqrt(dble( n_states_bra ) * dble( n_states_ket ))
          end if
        end if
        if(.not. present(mode)) fmat%m(Nbra+1, Nket+1) = get_frobenius_norm(v)
        deallocate(v)
      end do
    end do
  end function get_norm_op_chan

  function get_frobenius_norm(v) result(r)
    real(8), intent(in) :: v(:,:)
    real(8) :: r
    integer :: i, j
    r = 0.d0
    !$omp parallel
    !$omp do private(i,j) reduction(+:r)
    do i = 1, size(v,1)
      do j = 1, size(v,2)
        r = r + v(i,j)**2
      end do
    end do
    !$omp end do
    !$omp end parallel
    r = sqrt(r)
  end function get_frobenius_norm

  function GetMemory(this) result(r)
    class(ThreeBodyJacOpChanIso), intent(in) :: this
    real(8) :: r
    r = 8.d0 * this%n_row * this%n_col / 1024.d0**3
  end function GetMemory

  subroutine ReducedToNonReduced(this, convJ, convT)
    class(ThreeBodyJacOpChanIso), intent(inout), target :: this
    logical, intent(in), optional :: convJ, convT
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    real(8) :: fact
    logical :: cJ = .true.
    logical :: cT = .true.
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
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
  end subroutine ReducedToNonReduced

  subroutine NonReducedToReduced(this, convJ, convT)
    class(ThreeBodyJacOpChanIso), intent(inout), target :: this
    logical, intent(in), optional :: convJ, convT
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    real(8) :: fact
    logical :: cJ = .true.
    logical :: cT = .true.
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
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
  end subroutine NonReducedToReduced

  subroutine Asym_to_NonAsym( this )
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    if( .not. this%Antisymmetrized ) then
      write(*,*) "Oh, you already have non-antisymmetrized matrix elements, see at ", &
          & __LINE__, " in ", __FILE__
      return
    end if
    this%antisymmetrized = .false.
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    if( chbra%GetNumberAStates() < 1 .or. chbra%GetNumberNAStates() < 1) return
    if( chket%GetNumberAStates() < 1 .or. chket%GetNumberNAStates() < 1) return
    call work%ini( chbra%GetNumberAStates(), chket%GetNumberAStates() )
    work = this%DMat
    call this%DMat%ini( chbra%GetNumberNAStates(), chket%GetNumberNAStates() )
    cfp1 = chbra%GetCFPMat()
    cfp2 = chket%GetCFPMat()
    this%DMat = cfp1 * work * cfp2%T()
    call cfp1%fin()
    call cfp2%fin()
    call work%fin()
  end subroutine Asym_to_NonAsym

  subroutine from_NonAsym_to_Asym( this )
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    if( this%Antisymmetrized ) then
      write(*,*) "Oh, you already have antisymmetrized matrix elements, see at ", &
          & __LINE__, " in ", __FILE__
      return
    end if

    this%antisymmetrized = .true.
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    if( chbra%GetNumberAStates() < 1 .or. chbra%GetNumberNAStates() < 1) return
    if( chket%GetNumberAStates() < 1 .or. chket%GetNumberNAStates() < 1) return
    call work%ini( chbra%GetNumberNAStates(), chket%GetNumberNAStates() )
    work = this%DMat
    call this%DMat%ini( chbra%GetNumberAStates(), chket%GetNumberAStates() )
    cfp1 = chbra%GetCFPMat()
    cfp2 = chket%GetCFPMat()
    this%DMat = cfp1%T() * work * cfp2
    call cfp1%fin()
    call cfp2%fin()
    call work%fin()
  end subroutine from_NonAsym_to_Asym

  subroutine PrintNonAntisymmetrizedME( this, iunit )
    class(ThreeBodyJacOpChanIso), intent(inout) :: this
    integer, intent(in), optional :: iunit
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    integer :: ibra, iket
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra, ket
    integer :: wunit
    if( this%Antisymmetrized ) then
      write(*,*) "Oh, you should have non-antisymmetrized ME first", __LINE__, " in ", __FILE__
      return
    end if
    wunit = 6
    if( present(iunit)) wunit = iunit
    chbra => this%jacobi_ch_bra
    chket => this%jacobi_ch_ket
    if( chbra%GetNumberNAStates() < 1) return
    if( chket%GetNumberNAStates() < 1) return
    do ibra = 1, chbra%GetNumberNAStates()
      bra => chbra%GetNAS(ibra)
      do iket = 1, chket%GetNumberNAStates()
        ket => chket%GetNAS(iket)
        if( abs(this%m(ibra,iket)) < 1.d-8 ) cycle
        write(wunit,'(20i4,f18.8)') bra%GetN12(), bra%GetL12(), bra%GetS12(), bra%GetJ12(), bra%GetT12(), &
            & bra%GetN3(), bra%GetL3(), bra%GetJ3(), chbra%GetJ(), chbra%GetT(), &
            & ket%GetN12(), ket%GetL12(), ket%GetS12(), ket%GetJ12(), ket%GetT12(), &
            & ket%GetN3(), ket%GetL3(), ket%GetJ3(), chket%GetJ(), chket%GetT(), this%m(ibra,iket)
      end do
    end do
  end subroutine PrintNonAntisymmetrizedME
end module ThreeBodyJacOpsChanIso
