module ThreeBodyJacChanIso
  use omp_lib
  use StoreCouplings
  use LinAlgLib
  use ThreeBodyJacobiQuantumNumbers
  use ClassSys
  implicit none

  public :: ThreeBodyJacIsoChan
  public :: NmaxChannelIso3
  !-- methods for ThreeBodyJacIsoChan --
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetNmax
  private :: GetFrequency
  private :: GetNumberNAStates
  private :: GetNumberAStates
  private :: GetNumberNAStatesNmax
  private :: GetNumberAStatesNmax
  private :: GetNASFromIndex
  private :: GetNASFromEalpha
  private :: GetNASFromI
  private :: GetNASIndexFromEalpha
  private :: GetNASIndexFromJacobi
  private :: GetNASIndexFromQNs
  private :: GetASIndexFromEi
  private :: GetNmaxChannel
  private :: FinThreeBodyJacIsoChan
  private :: InitThreeBodyJacIsoChan
  private :: set_indices
  private :: FinNmaxChannelIso3
  private :: InitNmaxChannelIso3
  private :: set_jacobi_nas
  private :: get_dim_as
  private :: set_cfp
  private :: anti_op_isospin_wrap
  private :: anti_op_isospin
  private :: GetFileNameThreeBodyJacIsoChan
  private :: WriteThreeBodyJacIsoChan
  private :: ReadThreeBodyJacIsoChan
  private :: GetCFPMatThreeBodyJacIsoChan
  private :: FindIndexJacobiChannelFromQNs
  private :: FindIndexJacobiChannel
  private :: SetIndexJacobiChannelFromQNs
  private :: SetIndexJacobiChannel
  private :: FinJs
  private :: InitJs
  private :: FinSpinIsospin
  private :: InitSpinIsospin
  private :: FinIndexJacobiChannel
  private :: InitIndexJacobiChannel
  private :: GetMemoryChannel
  private :: GetMemoryNmaxChannel

  type, private :: Js
    integer, allocatable :: idx(:)
    integer, allocatable :: jj(:,:)
    integer :: nidx
  contains
    procedure, private :: init => initJs
    procedure, private :: fin => FinJs
  end type Js

  type, private :: SpinIsospin
    type(Js), allocatable :: idx(:)
    integer :: nidx
    integer, allocatable :: st(:,:)
  contains
    procedure, private :: init => InitSpinIsospin
    procedure, private :: fin => FinSpinIsospin
  end type SpinIsospin

  type, private :: IndexJacobiChannel
    type(SpinIsospin), allocatable :: idx(:)
    integer :: Nmax
    integer :: nidx
    integer, allocatable :: nl(:,:)
    integer, allocatable :: nls(:,:)
  contains
    procedure :: InitIndexJacobiChannel
    procedure :: FinIndexJacobiChannel
    procedure :: FindIndexJacobiChannel
    procedure :: SetIndexJacobiChannel
    procedure :: FindIndexJacobiChannelFromQNs
    procedure :: SetIndexJacobiChannelFromQNs
    generic :: init => InitIndexJacobiChannel
    generic :: fin => FinIndexJacobiChannel
    generic :: find => FindIndexJacobiChannel, FindIndexJacobiChannelFromQNs
    generic :: set => SetIndexJacobiChannel, SetIndexJacobiChannelFromQNs
  end type IndexJacobiChannel

  type, extends(DMat) :: NmaxChannelIso3
    type(NonAntisymmetrizedIsoQNs3), allocatable :: NAS(:)
    integer, private, allocatable :: NASIndex(:)
    integer, private, allocatable :: ASIndex(:)
    integer, private :: n_states_nas = 0
    integer, private :: n_states_as = 0
  contains
    procedure :: InitNmaxChannelIso3
    procedure :: FinNmaxChannelIso3
    procedure :: GetNumberNAStatesNmax
    procedure :: GetNumberAStatesNmax
    procedure :: GetNASFromI
    procedure :: GetMemoryNmaxChannel
    generic :: init => InitNmaxChannelIso3
    generic :: release => FinNmaxChannelIso3
    generic :: GetNumberNAStates => GetNumberNAStatesNmax
    generic :: GetNumberAStates => GetNumberAStatesNmax
    generic :: GetNAS => GetNASFromI
    generic :: GetMemory => GetMemoryNmaxChannel
  end type NmaxChannelIso3

  type :: ThreeBodyJacIsoChan
    real(8), private :: hw = -1.d0
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: t = -1
    integer, private :: Nmax = -1
    ! CFP
    type(NmaxChannelIso3), allocatable :: cfp(:)
    ! Non-antisymmetrized states
    type(IndexJacobiChannel) :: Jacobi2Index
    integer, private :: n_states_nas = 0
    integer, private, allocatable :: N_nas(:)
    integer, private, allocatable :: alpha_nas(:)
    ! Antisymmetrized states
    integer, private :: n_states_as = 0
    integer, private, allocatable :: N_as(:)
    integer, private, allocatable :: i_as(:)
  contains
    procedure :: InitThreeBodyJacIsoChan
    procedure :: FinThreeBodyJacIsoChan
    procedure :: GetFileNameThreeBodyJacIsoChan
    procedure :: WriteThreeBodyJacIsoChan
    procedure :: ReadThreeBodyJacIsoChan
    procedure :: GetCFPMatThreeBodyJacIsoChan
    procedure :: GetNASFromIndex
    procedure :: GetNASFromEalpha
    procedure :: GetNASIndexFromEalpha
    procedure :: GetNASIndexFromJacobi
    procedure :: GetNASIndexFromQNs
    procedure :: GetASIndexFromEi
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetNmax
    procedure :: GetFrequency
    procedure :: GetNumberNAStates
    procedure :: GetNumberAStates
    procedure :: GetNmaxChannel
    procedure :: GetENAS
    procedure :: GetEAS
    procedure :: GetiAS
    procedure :: GetAlphaNAS
    procedure :: GetMemoryChannel

    generic :: init => InitThreeBodyJacIsoChan
    generic :: fin => FinThreeBodyJacIsoChan
    generic :: GetFileName => GetFileNameThreeBodyJacIsoChan
    generic :: writef => WriteThreeBodyJacIsoChan
    generic :: readf => ReadThreeBodyJacIsoChan
    generic :: GetCFPMat => GetCFPMatThreeBodyJacIsoChan
    generic :: GetNAS => GetNASFromIndex, GetNASFromEalpha
    generic :: GetNASIndex => GetNASIndexFromEalpha,&
        & GetNASIndexFromJacobi, GetNASIndexFromQNs
    generic :: GetASIndex => GetASIndexFromEi
    generic :: GetMemory => GetMemoryChannel
  end type ThreeBodyJacIsoChan

  ! ------- store angular momentum couplings
  type(SixJsStore), private :: sixjs
  type(NineJsStore), private :: ninejs

  type(sys), private :: sy
contains
  function GetJ(this) result(j)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(t)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetNmax(this) result(Nmax)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetENAS(this,i) result(E)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: i
    integer :: E
    E = this%N_nas(i)
  end function GetENAS

  function GetEAS(this,i) result(E)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: i
    integer :: E
    E = this%N_as(i)
  end function GetEAS

  function GetAlphaNAS(this,i) result(alpha)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: i
    integer :: alpha
    alpha = this%alpha_nas(i)
  end function GetAlphaNAS

  function GetiAS(this,i) result(idx)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%i_as(i)
  end function GetiAS

  function GetFrequency(this) result(hw)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetNumberNAStates(this) result(Num)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: Num
    Num = this%n_states_nas
  end function GetNumberNAStates

  function GetNumberAStates(this) result(Num)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer :: Num
    Num = this%n_states_as
  end function GetNumberAStates

  function GetNumberNAStatesNmax(this) result(Num)
    class(NmaxChannelIso3), intent(in) :: this
    integer :: Num
    Num = this%n_states_nas
  end function GetNumberNAStatesNmax

  function GetNumberAStatesNmax(this) result(Num)
    class(NmaxChannelIso3), intent(in) :: this
    integer :: Num
    Num = this%n_states_as
  end function GetNumberAStatesNmax

  function GetNASFromI(this, i) result(r)
    class(NmaxChannelIso3), intent(in), target :: this
    integer, intent(in) :: i
    type(NonAntisymmetrizedIsoQNs3), pointer :: r
    r => null()
    if(i == 0) return
    r => this%NAS(i)
  end function GetNASFromI

  function GetNASFromIndex(this, i) result(r)
    class(ThreeBodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: i
    type(NonAntisymmetrizedIsoQNs3), pointer :: r
    r => this%GetNAS( this%N_nas(i), this%alpha_nas(i) )
  end function GetNASFromIndex

  function GetNASFromEalpha(this, E, alpha) result(r)
    class(ThreeBodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: E, alpha
    type(NonAntisymmetrizedIsoQNs3), pointer :: r
    r => this%cfp(E)%GetNAS(alpha)
  end function GetNASFromEalpha

  function GetNASIndexFromEalpha(this, E, alpha) result(idx)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: E, alpha
    integer :: idx
    idx = this%cfp(E)%NASIndex(alpha)
  end function GetNASIndexFromEalpha

  function GetNASIndexFromJacobi(this, jacobi) result(idx)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    type(NonAntisymmetrizedIsoQNs3), intent(in) :: jacobi
    integer :: idx
    idx = this%Jacobi2Index%find(jacobi)
  end function GetNASIndexFromJacobi

  function GetNASIndexFromQNs(this,n12,l12,s12,j12,t12,n3,l3,j3) result(idx)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: n12,l12,s12,j12,t12,n3,l3,j3
    integer :: idx
    idx = this%Jacobi2Index%find(n12,l12,s12,j12,t12,n3,l3,j3)
  end function GetNASIndexFromQNs

  function GetASIndexFromEi(this, E, i) result(idx)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: E, i
    integer :: idx
    idx = this%cfp(E)%ASIndex(i)
  end function GetASIndexFromEi

  function GetNmaxChannel(this, Nmax) result(r)
    class(ThreeBodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: Nmax
    type(NmaxChannelIso3), pointer :: r
    r => this%cfp(Nmax)
  end function GetNmaxChannel

  function GetMemoryChannel(this) result(r)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    real(8) :: r
    integer :: Nmax
    integer :: idx1, idx2
    r = 0.d0
    ! Index
    r = r + dble( size(this%Jacobi2Index%nl) ) * 4.d0 / 1024.d0**3
    r = r + dble( size(this%Jacobi2Index%nls) ) * 4.d0 / 1024.d0**3
    do idx1 = 1, size(this%Jacobi2Index%idx)
      r = r + dble( size(this%Jacobi2Index%idx(idx1)%st) ) * 4.d0 / 1024.d0**3
      do idx2 = 1, size(this%Jacobi2Index%idx(idx1)%idx)
        r = r + dble( size(this%Jacobi2Index%idx(idx1)%idx(idx2)%jj) ) * 4.d0 / 1024.d0**3
        r = r + dble( size(this%Jacobi2Index%idx(idx1)%idx(idx2)%idx) ) * 4.d0 / 1024.d0**3
      end do
    end do

    do Nmax = 0, this%GetNmax()
      r = r + this%cfp(Nmax)%GetMemory()
    end do
    r = r + ( dble( this%GetNumberNAStates() ) * 8.d0 / 1024.d0**3 )
    r = r + ( dble( this%GetNumberAStates() ) * 8.d0 / 1024.d0**3 )
  end function GetMemoryChannel

  function GetMemoryNmaxChannel(this) result(r)
    class(NmaxChannelIso3), intent(in) :: this
    real(8) :: r
    r = 0.d0
    r = dble(this%n_row) * dble(this%n_col) * 8.d0 / 1024.d0**3
    r = r + dble( this%GetNumberNAStates() ) * 36.d0 / 1024.d0**3
    r = r + dble( this%GetNumberAStates() ) * 4.d0 / 1024.d0**3
  end function GetMemoryNmaxChannel

  subroutine FinThreeBodyJacIsoChan(this)
    class(ThreeBodyJacIsoChan), intent(inout) :: this
    integer :: N
    do N = 0, this%GetNmax()
      call this%cfp(N)%release()
    end do
    deallocate(this%cfp)
    call this%Jacobi2Index%fin()
    if(this%GetNumberNAStates() * this%GetNumberAStates() == 0) return
    deallocate(this%N_nas)
    deallocate(this%alpha_nas)
    deallocate(this%N_as)
    deallocate(this%i_as)
    this%n_states_as = 0
    this%n_states_nas = 0
    this%hw = -1.d0
    this%j = -1
    this%p = -100
    this%t = -1
    this%Nmax = -1
  end subroutine FinThreeBodyJacIsoChan

  subroutine InitThreeBodyJacIsoChan(this, hw, j, p, t, Nmax, path_to_dir)
    use ClassSys, only: sys, str
    use Profiler, only: timer
    use MPIFunction, only: myrank
    class(ThreeBodyJacIsoChan), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: j, p, t, Nmax
    type(str) :: path_to_dir
    integer :: N, i
    real(8) :: ti
    type(sys) :: s
    type(str) :: f
    integer :: runit = 30
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas

    this%hw = hw
    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%t = t

    f = this%GetFileName(j,p,t,Nmax,path_to_dir)
    if(s%isfile(f)) then
      open(runit,file=f%val,status='old',form='unformatted',access='stream')
      call this%readf(hw,runit)
      close(runit)
      return
    end if
    call sixjs%init(1,1,.true., 1,1,.true., 1,1,.true.)
    call ninejs%init(1,1,.true., 1,2*Nmax+1,.true., 1,3,.true., this%j,this%j,.true.) ! avoid overflow

    ti = omp_get_wtime()
    allocate(this%cfp(0:this%Nmax))
    do N = 0, Nmax
      call this%cfp(N)%init(j, p, t, N)
    end do

    call ninejs%fin()
    call sixjs%fin()
    call set_indices(this)
    call this%Jacobi2Index%init(this%GetNmax(), this%GetJ(), this%GetParity(), this%GetT())
    do i = 1, this%GetNumberNAStates()
      nas => this%GetNAS(i)
      call this%Jacobi2Index%set(nas, i)
    end do
    write(*,"(a,i4,a,i3,a,i2,a,i2,a,2i8)") "myrank=", myrank, ", J=", this%j, "/2, P=", this%p, ", T=", this%t, &
        & "/2, Number of physical and orthonormal states: ", this%n_states_nas, this%n_states_as
    call timer%add(sy%str('InitThreeBodyJacIsoChan'), omp_get_wtime()-ti)
  end subroutine InitThreeBodyJacIsoChan

  subroutine set_indices(this)
    class(ThreeBodyJacIsoChan), intent(inout) :: this
    integer :: N, i, n_states_nas, n_states_as
    this%n_states_nas = 0
    this%n_states_as = 0
    do N = 0, this%GetNmax()
      allocate(this%cfp(N)%NASIndex(this%cfp(N)%n_states_nas))
      allocate(this%cfp(N)%ASIndex( this%cfp(N)%n_states_as ))
      this%n_states_nas= this%n_states_nas + this%cfp(N)%n_states_nas
      this%n_states_as = this%n_states_as  + this%cfp(N)%n_states_as
    end do

    if(this%n_states_nas * this%n_states_as == 0) return
    allocate(this%N_nas(    this%n_states_nas))
    allocate(this%alpha_nas(this%n_states_nas))
    allocate(this%N_as(this%n_states_as))
    allocate(this%i_as(this%n_states_as))

    n_states_nas = 0
    do N = 0, this%GetNmax()
      if(this%cfp(N)%n_states_nas == 0) cycle
      do i = 1, this%cfp(N)%n_states_nas
        n_states_nas = n_states_nas + 1
        this%N_nas(    n_states_nas) = N
        this%alpha_nas(n_states_nas) = i
        this%cfp(N)%NASIndex(i) = n_states_nas
      end do
    end do

    n_states_as = 0
    do N = 0, this%GetNmax()
      if(this%cfp(N)%n_states_as == 0) cycle
      do i = 1, this%cfp(N)%n_states_as
        n_states_as = n_states_as + 1
        this%N_as(n_states_as) = N
        this%i_as(n_states_as) = i
        this%cfp(N)%ASIndex(i) = n_states_as
      end do
    end do
  end subroutine set_indices

  subroutine FinNmaxChannelIso3(this)
    class(NmaxChannelIso3), intent(inout) :: this
    call this%fin()
    if(this%n_states_nas * this%n_states_as == 0) return
    deallocate(this%NAS)
    deallocate(this%NASIndex)
    deallocate(this%ASIndex)
    this%n_states_nas = 0
    this%n_states_as = 0
  end subroutine FinNmaxChannelIso3

  subroutine InitNmaxChannelIso3(this, j, p, t, Nmax)
    class(NmaxChannelIso3), intent(inout) :: this
    integer, intent(in) :: j, p, t, Nmax
    call set_jacobi_nas(this, j, p, t, Nmax, this%n_states_nas, .false.)
    if(this%n_states_nas == 0) return
    call set_jacobi_nas(this, j, p, t, Nmax, this%n_states_nas, .true.)

    this%n_states_as = get_dim_as(this, j, t)
    if(this%n_states_as == 0) then
      this%n_states_nas = 0
      deallocate(this%NAS)
      return
    end if
    call this%ini(this%GetNumberAStates(), this%GetNumberNAStates())
    call set_cfp(this, j, t)
  end subroutine InitNmaxChannelIso3

  subroutine set_jacobi_nas(this, j, p, t, N, n_states_nas, set_mode)
    use MyLibrary, only: triag
    class(NmaxChannelIso3), intent(inout) :: this
    integer, intent(in) :: j, p, t, N
    integer, intent(inout) :: n_states_nas
    logical, intent(in) :: set_mode
    integer :: n_states
    integer :: n12, n3, l12, l3, s12, j12, t12, j3
    if(set_mode) allocate(this%NAS( this%n_states_nas ) )
    n_states = 0
    do n12 = 0, N / 2
      do n3 = 0, (N - 2 * n12) / 2
        do l12 = 0, (N - 2 * n12 - 2 * n3)
          l3 = N - 2 * n12 - 2 * n3 - l12
          if((-1) ** (l12 + l3) /= p) cycle
          do s12 = 0, 1
            do j12 = iabs(l12 - s12), l12 + s12
              do t12 = 0, 1
                if(triag(2 * t12,  1, t)) cycle
                if((-1) ** (l12 + s12 + t12) /= -1) cycle
                do j3 = iabs(2 * l3 - 1), 2 * l3 + 1, 2
                  if(triag(2 * j12, j3, j)) cycle
                  n_states = n_states + 1
                  if(set_mode) call this%NAS(n_states)%set(n12,l12,s12,j12,t12,n3,l3,j3)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    if(.not. set_mode) n_states_nas = n_states
  end subroutine set_jacobi_nas

  function get_dim_as(this, j, t) result(n_states)
    class(NmaxChannelIso3), intent(inout), target :: this
    integer, intent(in) :: j, t
    integer :: n_states
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas
    real(8), allocatable :: diag(:)
    integer :: n_nas
    allocate(diag( this%n_states_nas ))
    diag(:) = 0.d0
    do n_nas = 1, this%n_states_nas
      nas => this%NAS(n_nas)
      diag(n_nas) = anti_op_isospin_wrap(nas, nas, j, t)
    end do
    n_states = int(sum(diag) + 1.d-2)
    if(abs(sum(diag) - dble(n_states)) > 1.d-4) then
      write(*,*) "error in trace of A operator: tr(A) = ", sum(diag)
      stop
    end if
    deallocate(diag)
  end function get_dim_as

  subroutine set_cfp(this, j, t)
    use MyLibrary, only: triag
    use LinAlgLib
    class(NmaxChannelIso3), intent(inout), target :: this
    integer, intent(in) :: j, t
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas1, nas2
    integer :: bra, ket, i, k, ii
    real(8) :: anti
    type(DMat) :: A
    type(EigenSolSymD) :: sol

    call A%ini(this%n_states_nas, this%n_states_nas)
    !$omp parallel
    !$omp do private(bra, nas1, ket, nas2, anti) schedule(dynamic)
    do bra = 1, this%n_states_nas
      nas1 => this%NAS(bra)
      do ket = 1, bra
        nas2 => this%NAS(ket)
        anti = anti_op_isospin_wrap(nas1, nas2, j, t)
        A%m(bra,ket) = anti
        A%m(ket,bra) = anti
      end do
    end do
    !$omp end do
    !$omp end parallel
    call sol%init(A)
    call sol%DiagSym(A)
    k = 0
    i = 0
    do ii = 1, this%n_states_nas
      if(abs(1.d0 - sol%eig%v(ii)) < 1.d-2) then
        k = k + 1
        this%m(k,:) = sol%vec%m(:,ii)
      else
        if(abs(sol%eig%v(ii)) > 1.d-2) then
          i = i + 1
        end if
      end if
    end do
    if( i /= 0 ) then
      write(*,'(a)') 'Warning: SetThreeBodyIsospinJacNmax'
      call sol%eig%prt()
    end if
    call sol%fin()
    call A%fin()
  end subroutine set_cfp

  real(8) function anti_op_isospin_wrap(bra, ket, j, t) result(anti)
    type(NonAntisymmetrizedIsoQNs3), intent(in) :: bra, ket
    integer, intent(in) :: j, t
    integer :: n12, l12, s12, j12, t12, n3, l3, j3
    integer :: n45, l45, s45, j45, t45, n6, l6, j6
    n12 = bra%GetN12()
    l12 = bra%GetL12()
    s12 = bra%GetS12()
    j12 = bra%GetJ12()
    t12 = bra%GetT12()
    n3  = bra%GetN3()
    l3  = bra%GetL3()
    j3  = bra%GetJ3()
    n45 = ket%GetN12()
    l45 = ket%GetL12()
    s45 = ket%GetS12()
    j45 = ket%GetJ12()
    t45 = ket%GetT12()
    n6  = ket%GetN3()
    l6  = ket%GetL3()
    j6  = ket%GetJ3()
    anti = anti_op_isospin(n12, l12, s12, j12, t12, n3, l3, j3, &
        &                  n45, l45, s45, j45, t45, n6, l6, j6, &
        &                  j, t)
  end function anti_op_isospin_wrap

  ! <(n12,l12,s12,j12,t12; n3, l3, j3): JT| A |(n45,l45,s45,j45,t45; n6, l6, j6): JT>
  ! <(n12,l12,s12,j12,t12; n3, l3, j3): JT| 1 - 2 * T |(n45,l45,s45,j45,t45; n6, l6, j6): JT> / 3
  ! <(n12,l12,s12,j12,t12; n3, l3, j3): JT| T |(n45,l45,s45,j45,t45; n6, l6, j6): JT>
  !          = delta_{N123, N456} sqrt(2 * t12 + 1) * sqrt(2 * t45 + 1)
  !            x {1/2  1/2  t12}  sum_{LS} (2 * L + 1) * (2 * S + 1)
  !              {1/2   T   t45}
  !            x sqrt(2 * j12 + 1) * sqrt(2 * j45 + 1) * sqrt(2 * j3 + 1) * sqrt(2 * j6 + 1)
  !            x sqrt(2 * s12 + 1) * sqrt(2 * s45 + 1) {1/2  1/2  s12}
  !                                                    {1/2  S    s45}
  !                       {l12  s12  j12}  {l45  s45  j45}
  !            x (-1)^{L} {l3   1/2   j3}  {l6   1/2   j6} <<n12 l12 n3 l3; L| n6 l6 n45 l45; L>>_{3}
  !                       {L    S      J}  {L    S      J}
  ! A: Antisymmetrizer
  ! N123: 2 * n12 + l12 + 2 * n3 + l3
  ! << | >>_{x} : Transformation braket for mass ratio x
  ! reference: Progress in Particle and Nuclear Physics 69 (2013) 131–181
  !            http://www.sciencedirect.com/science/article/pii/S0146641012001184
  real(8) function anti_op_isospin(n12, l12, s12, j12, t12, n3, l3, j3, &
        &                      n45, l45, s45, j45, t45, n6, l6, j6, &
        &                      jtot, ttot) result(anti)
    use MyLibrary, only: gmosh
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3
    integer, intent(in) :: n45, l45, s45, j45, t45, n6, l6, j6
    integer, intent(in) :: jtot, ttot
    integer :: lambda, lambda_min, lambda_max, stot, is_min, is_max
    real(8) :: ex

    anti = 0.d0; ex = 0.d0
    if(2 * n12 + l12 + 2 * n3 + l3 /= 2 * n45 + l45 + 2 * n6 + l6) return

    is_min = max(iabs(2 * s12 - 1), iabs(2 * s45 - 1))
    is_max = min(2 * s12 + 1, 2 * s45 + 1)

    do stot = is_min, is_max, 2
      lambda_min = max(iabs(l12 - l3), iabs(l45 - l6), iabs(jtot - stot)/2)
      lambda_max = min(l12 + l3, l45 + l6, (jtot + stot)/2)
      do lambda = lambda_min, lambda_max
        ex = ex + &
            &  dsqrt(dble(2*j12+1) * dble(j3+1) * dble(2*lambda+1) * dble(stot+1)) * &
            &  dsqrt(dble(2*j45+1) * dble(j6+1) * dble(2*lambda+1) * dble(stot+1)) * &
            !&  ninejs%get(2*l12,2*s12,2*j12,2*l3,1,j3,2*lambda,stot,jtot) * &
            !&  ninejs%get(2*l45,2*s45,2*j45,2*l6,1,j6,2*lambda,stot,jtot) * &
            &  ninejs%get(1,j3,2*l3,stot,jtot,2*lambda,2*s12,2*j12,2*l12) * &
            &  ninejs%get(1,j6,2*l6,stot,jtot,2*lambda,2*s45,2*j45,2*l45) * &
            &  dsqrt(dble(2*s12+1) * dble(2*s45+1)) * &
            &  sixjs%get(1,1,2*s12,1,stot,2*s45) * &
            &  gmosh(n12, l12, n3, l3, n45, l45, n6, l6, lambda, 1.d0/3.d0)
        !ex = ex + &
        !    &  init4j(l12, s12, l3, 1)%inter4j(j12, j3, lambda, stot)%v(jtot) * &
        !    &  init4j(l45, s45, l6, 1)%inter4j(j45, j6, lambda, stot)%v(jtot) * &
        !    &  init3j(1, 1, 1)%inter2j(s12, s45)%v(stot) * &
        !    &  gmosh(n12, l12, n3, l3, n45, l45, n6, l6, lambda, 1.d0/3.d0)
      end do
    end do
    !ex = ex * init3j(1, 1, 1)%inter2j(t12, t45)%v(ttot) * (-1.d0) ** (s12 + t12 + s45 + t45)
    ex = ex * (-1.d0) ** (s12 + t12 + s45 + t45) * &
            &  dsqrt(dble(2*t12+1) * dble(2*t45+1)) * &
            &  sixjs%get(1,1,2*t12,1,ttot,2*t45)
    anti = - 2.d0 * ex / 3.d0
    if(n12 == n45 .and. l12 == l45 .and. s12 == s45 .and. j12 == j45 .and. &
        & t12 == t45 .and. n3 == n6 .and. l3 == l6 .and. j3 == j6) then
      anti = (1.d0 - 2.d0 * ex ) / 3.d0
    end if
  end function anti_op_isospin

  function GetFileNameThreeBodyJacIsoChan(this, j, p, t, Nmax, path_to_dir) result(f)
    use ClassSys, only: sys, str
    class(ThreeBodyJacIsoChan), intent(inout) :: this
    integer, intent(in) :: j, p, t, Nmax
    type(str), intent(in) :: path_to_dir
    type(str) :: f, dir, pari
    type(sys) :: s
    this%j = j
    this%p = p
    this%t = t
    this%Nmax = Nmax
    dir = path_to_dir + s%str('/cfp')
    call s%mkdir(dir%val)
    if(p == -1) pari = '-'
    if(p ==  1) pari = '+'
    f = dir + s%str('/cfp_A3_j') + s%str(j) + s%str('p') + pari + &
        & s%str('t') + s%str(t) + s%str('_Nmax') + s%str(Nmax) + s%str('.bin')
  end function GetFileNameThreeBodyJacIsoChan

  subroutine WriteThreeBodyJacIsoChan(this, iunit)
    use Profiler, only: timer
    class(ThreeBodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: iunit
    type(NmaxChannelIso3), pointer :: ch_nmax
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas
    integer :: N, i
    real(8) :: ti

    ti = omp_get_wtime()
    write(iunit) this%GetJ(), this%GetParity(), this%GetT(), this%GetNmax()
    do N = 0, this%GetNmax()
      ch_nmax => this%GetNmaxChannel(N)
      write(iunit) ch_nmax%GetNumberNAStates(), ch_nmax%GetNumberAStates(), N
      !if(ch_nmax%GetNumberNAStates() * ch_nmax%GetNumberAStates() == 0) cycle
      if( ch_nmax%GetNumberNAStates() < 1 ) cycle
      do i = 1, ch_nmax%n_states_nas
        nas => ch_nmax%GetNAS(i)
        write(iunit) nas%GetN12()
        write(iunit) nas%GetL12()
        write(iunit) nas%GetS12()
        write(iunit) nas%GetJ12()
        write(iunit) nas%GetT12()
        write(iunit) nas%GetN3()
        write(iunit) nas%GetL3()
        write(iunit) nas%GetJ3()
      end do
      if( ch_nmax%GetNumberAStates() < 1 ) cycle
      do i = 1, ch_nmax%n_states_as
        write(iunit) this%cfp(N)%m(i,:)
      end do
    end do
    call timer%Add(sy%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteThreeBodyJacIsoChan

  subroutine ReadThreeBodyJacIsoChan(this, hw, iunit, Nmax)
    use Profiler, only: timer
    class(ThreeBodyJacIsoChan), intent(inout) :: this
    integer, intent(in) :: iunit
    real(8), intent(in) :: hw
    integer, intent(in), optional :: Nmax
    integer :: N, i, Nmax_read, N_read
    integer :: nphys, north
    integer :: n12, l12, s12, j12, t12, n3, l3, j3
    type(NmaxChannelIso3), pointer :: ch_nmax
    type(NonAntisymmetrizedIsoQNs3), pointer :: nas
    real(8) :: ti

    if(allocated(this%cfp)) call this%fin()
    ti = omp_get_wtime()
    this%hw = hw

    read(iunit) this%j, this%p, this%t, Nmax_read

    if(present(Nmax)) then; this%Nmax = min(Nmax,Nmax_read)
    else; this%Nmax = Nmax_read; end if
    allocate(this%cfp( 0:this%GetNmax() ))

    do N = 0, min(this%GetNmax(),Nmax_read)
      read(iunit) nphys, north, N_read
      ch_nmax => this%GetNmaxChannel(N)
      if(N /= N_read) cycle
      ch_nmax%n_states_nas = nphys
      ch_nmax%n_states_as  = north

      if( ch_nmax%GetNumberNAStates() < 1 ) cycle
      allocate( ch_nmax%NAS(nphys) )
      do i = 1, ch_nmax%GetNumberNAStates()
        read(iunit) n12, l12, s12, j12, t12, n3, l3, j3
        call ch_nmax%NAS(i)%set(n12,l12,s12,j12,t12,n3,l3,j3)
      end do

      if( ch_nmax%GetNumberAStates() < 1 ) cycle
      call this%cfp(N)%ini( north, nphys )
      do i = 1, ch_nmax%GetNumberAStates()
        read(iunit) this%cfp(N)%m(i,:)
      end do
    end do

    call set_indices(this)
    call this%Jacobi2Index%init(this%GetNmax(), this%GetJ(), this%GetParity(), this%GetT())
    do i = 1, this%GetNumberNAStates()
      nas => this%GetNAS(i)
      call this%Jacobi2Index%set(nas, i)
    end do
    call timer%Add(sy%str('Read from file'), omp_get_wtime() - ti)
  end subroutine ReadThreeBodyJacIsoChan

  function GetCFPMatThreeBodyJacIsoChan(this) result(cfp)
    class(ThreeBodyJacIsoChan), intent(in) :: this
    type(DMat) :: cfp
    integer :: n_nas, n_as
    integer :: bra, ket, Nmax, l, i
    n_nas = this%GetNumberNAStates()
    n_as  = this%GetNumberAStates()
    call cfp%zeros(n_nas, n_as)
    !$omp parallel
    !$omp do private(bra, Nmax, i, ket, l)
    do bra = 1, n_as
      Nmax = this%N_as(bra)
      i    = this%i_as(bra)
      do ket = 1, n_nas
        if(Nmax /= this%N_nas(ket)) cycle
        l = this%alpha_nas(ket)
        cfp%m(ket,bra) = this%cfp(Nmax)%m(i,l)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end function GetCFPMatThreeBodyJacIsoChan

  function FindIndexJacobiChannelFromQNs(this, jacobi) result(idx)
    class(IndexJacobiChannel), intent(in) :: this
    type(NonAntisymmetrizedIsoQNs3), intent(in) :: jacobi
    integer :: n12, l12, s12, j12, t12, n3, l3, j3
    integer :: idx
    n12 = jacobi%GetN12()
    l12 = jacobi%GetL12()
    s12 = jacobi%GetS12()
    j12 = jacobi%GetJ12()
    t12 = jacobi%GetT12()
    n3  = jacobi%GetN3()
    l3  = jacobi%GetL3()
    j3  = jacobi%GetJ3()
    idx = this%find(n12,l12,s12,j12,t12,n3,l3,j3)
  end function FindIndexJacobiChannelFromQNs

  function FindIndexJacobiChannel(this, n12, l12, s12, j12, t12, n3, l3, j3) result(idx)
    class(IndexJacobiChannel), intent(in) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3
    integer :: idx12, idx3, idx_nls, idx_st, idx_jj
    integer :: idx
    idx = 0
    idx12 = this%nl(n12,l12)
    idx3  = this%nl(n3 ,l3 )
    idx_nls = this%nls(idx12, idx3)
    if(idx_nls==0) return
    idx_st = this%idx(idx_nls)%st(s12,t12)
    if(idx_st==0) return
    idx_jj = this%idx(idx_nls)%idx(idx_st)%jj(j12,j3)
    if(idx_jj==0) return
    idx = this%idx(idx_nls)%idx(idx_st)%idx(idx_jj)
  end function FindIndexJacobiChannel

  subroutine SetIndexJacobiChannelFromQNs(this, jacobi, idx)
    class(IndexJacobiChannel), intent(inout) :: this
    type(NonAntisymmetrizedIsoQNs3), intent(in) :: jacobi
    integer, intent(in) :: idx
    integer :: n12, l12, s12, j12, t12, n3, l3, j3
    n12 = jacobi%GetN12()
    l12 = jacobi%GetL12()
    s12 = jacobi%GetS12()
    j12 = jacobi%GetJ12()
    t12 = jacobi%GetT12()
    n3  = jacobi%GetN3()
    l3  = jacobi%GetL3()
    j3  = jacobi%GetJ3()
    call this%set(n12, l12, s12, j12, t12, n3, l3, j3, idx)
  end subroutine SetIndexJacobiChannelFromQNs

  subroutine SetIndexJacobiChannel(this, n12, l12, s12, j12, t12, n3, l3, j3, idx)
    class(IndexJacobiChannel), intent(inout) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3, idx
    integer :: idx12, idx3, idx_nls, idx_st, idx_jj
    idx12 = this%nl(n12,l12)
    idx3  = this%nl(n3 ,l3 )
    if(idx12 * idx3 == 0) then
      write(*,"(a,4i4)") "Error in SetIndexJacobiChannelIsospin, n12, l12, n3, l3:", n12, l12, n3, l3
      return
    end if
    idx_nls = this%nls(idx12, idx3)
    if(idx_nls == 0) then
      write(*,"(a,4i4)") "Error in SetIndexJacobiChannelIsospin, n12, l12, n3, l3:", n12, l12, n3, l3
      return
    end if
    idx_st = this%idx(idx_nls)%st(s12,t12)
    if(idx_st == 0) then
      write(*,"(a,2i4)") "Error in SetIndexJacobiChannelIsospin, s12, t12:", s12, t12
      return
    end if
    idx_jj = this%idx(idx_nls)%idx(idx_st)%jj(j12,j3)
    if(idx_jj == 0) then
      write(*,"(a,2i4)") "Error in SetIndexJacobiChannelIsospin: j12, 2*j3:", j12, j3
      return
    end if
    this%idx(idx_nls)%idx(idx_st)%idx(idx_jj) = idx
  end subroutine SetIndexJacobiChannel

  subroutine FinJs(this)
    class(Js), intent(inout) :: this
    if(allocated(this%idx)) deallocate(this%idx)
    if(allocated(this%jj)) deallocate(this%jj)
  end subroutine FinJs

  subroutine InitJs(this, l12, l3, s12, Jtot)
    use MyLibrary, only: triag
    class(Js), intent(inout) :: this
    integer, intent(in) :: l12, l3, s12, Jtot
    integer :: j12min, j12max, j3min, j3max
    integer :: idx, j12, j3
    j12min = abs(l12-s12)
    j12max = l12+s12
    j3min = abs(2*l3-1)
    j3max = 2*l3+1
    allocate(this%jj(j12min:j12max,j3min:j3max))
    this%jj(:,:) = 0
    idx = 0
    do j12 = j12min, j12max
      do j3 = j3min, j3max, 2
        if(triag(2*j12, j3, Jtot)) cycle
        idx = idx + 1
      end do
    end do
    this%nidx= idx
    allocate(this%idx(this%nidx))
    this%idx(:) = 0
    idx = 0
    do j12 = j12min, j12max
      do j3 = j3min, j3max, 2
        if(triag(2*j12, j3, Jtot)) cycle
        idx = idx + 1
        this%jj(j12,j3) = idx
      end do
    end do
  end subroutine InitJs

  subroutine FinSpinIsospin(this)
    class(SpinIsospin), intent(inout) :: this
    integer :: idx
    do idx = 1, this%nidx
      call this%idx(idx)%fin()
    end do
    deallocate(this%idx)
    deallocate(this%st)
  end subroutine FinSpinIsospin

  subroutine InitSpinIsospin(this, l12, l3, Jtot, Ttot)
    use MyLibrary, only: triag
    use MyLibrary, only: triag
    class(SpinIsospin), intent(inout) :: this
    integer, intent(in) :: l12, l3, Jtot, Ttot
    integer :: s12, t12, idx

    allocate(this%st(0:1,0:1))
    this%st(:,:) = 0
    idx = 0
    do s12 = 0, 1
      do t12 = 0, 1
        if(mod(l12+s12+t12, 2) == 0) cycle
        if(triag(2*t12, 1, Ttot)) cycle
        idx = idx + 1
      end do
    end do
    this%nidx = idx
    allocate(this%idx(this%nidx))
    idx = 0
    do s12 = 0, 1
      do t12 = 0, 1
        if(mod(l12+s12+t12, 2) == 0) cycle
        if(triag(2*t12, 1, Ttot)) cycle
        idx = idx + 1
        this%st(s12,t12) = idx
        call this%idx(idx)%init(l12, l3, s12, Jtot)
      end do
    end do
  end subroutine InitSpinIsospin

  subroutine FinIndexJacobiChannel(this)
    class(IndexJacobiChannel), intent(inout) :: this
    integer :: idx
    do idx = 1, this%nidx
      call this%idx(idx)%fin()
    end do
    deallocate(this%nl)
    deallocate(this%nls)
    deallocate(this%idx)
  end subroutine FinIndexJacobiChannel

  subroutine InitIndexJacobiChannel(this, Nmax, Jtot, Ptot, Ttot)
    class(IndexJacobiChannel), intent(inout) :: this
    integer, intent(in) :: Nmax, Jtot, Ptot, Ttot
    integer :: n, l, idx, nidx_nl
    integer :: n12, l12, n3, l3
    allocate(this%nl(0:Nmax/2,0:Nmax))
    this%nl(:,:) = 0
    idx = 0
    do n = 0, Nmax/2
      do l = 0, Nmax-2*n
        idx = idx + 1
        this%nl(n,l) = idx
      end do
    end do
    nidx_nl = idx
    allocate(this%nls(nidx_nl, nidx_nl))
    this%nls(:,:) = 0
    idx = 0
    do n12 = 0, Nmax/2
      do l12 = 0, Nmax-2*n12
        do n3 = 0, Nmax/2
          do l3 = 0, Nmax-2*n3
            if((-1)**(l12+l3) /= Ptot) cycle
            idx = idx + 1
          end do
        end do
      end do
    end do
    this%nidx = idx
    allocate(this%idx(this%nidx))
    idx = 0
    do n12 = 0, Nmax/2
      do l12 = 0, Nmax-2*n12
        do n3 = 0, Nmax/2
          do l3 = 0, Nmax-2*n3
            if((-1)**(l12+l3) /= Ptot) cycle
            idx = idx + 1
            this%nls(this%nl(n12,l12), this%nl(n3,l3)) = idx
            call this%idx(idx)%init(l12, l3, Jtot, Ttot)
          end do
        end do
      end do
    end do
  end subroutine InitIndexJacobiChannel
end module ThreeBodyJacChanIso
