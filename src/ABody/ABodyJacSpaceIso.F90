module ABodyJacSpaceIso
  use omp_lib
  use ClassSys
  use Profiler, only: timer
  use LinAlgLib
  use ThreeBodyJacobiSpace
  use ABodyJacQuantumNumbers
  implicit none

  public :: ABodyJacIsoSpace
  public :: ABodyJacIsoChan
  public :: NmaxChannelA

  private :: FinNmaxChannelA
  private :: InitNmaxChannel4
  private :: InitNmaxChannelA
  private :: GetNumberNAStatesNmax
  private :: GetNumberAStatesNmax
  private :: GetNASFromI
  private :: GetASIndex
  private :: GetNASIndex
  private :: set_jacobi_nas4
  private :: get_dim_as4
  private :: set_cfp4
  private :: set_jacobi_nasA
  private :: get_dim_asA
  private :: set_cfpA
  private :: anti_op4_wrap
  private :: Aop4
  private :: anti_op_wrap
  private :: Aop

  private :: FinABodyJacIsoChan
  private :: Init4BodyJacIsoChanFrom3
  private :: InitABodyJacIsoChan
  private :: ReadABodyJacIsoChan
  private :: WriteABodyJacIsoChan
  private :: GetCFPMatABodyJacIsoChan
  private :: GetNmaxChannel
  private :: GetFrequencyChan
  private :: GetJ
  private :: GetParity
  private :: GetT
  private :: GetNmaxChan
  private :: GetNumberParticlesChan
  private :: GetNumberNAStates
  private :: GetNumberAStates
  private :: GetNASIndexFromEi
  private :: GetASIndexFromEi
  private :: GetNASFromEalpha
  private :: GetNASFromIndex

  private :: FinABodyJacIsoSpace
  private :: Init4BodyJacIsoSpace
  private :: InitABodyJacIsoSpace
  private :: GetNmaxSpace
  private :: GetFrequencySpace
  private :: GetJmaxSpace
  private :: GetJminSpace
  private :: GetTmaxSpace
  private :: GetTminSpace
  private :: GetNumberParticlesSpace
  private :: GetNumberChannels
  private :: GetJPTIndex
  private :: GetChannelFromIndex
  private :: GetChannelFromJPT

  type, extends(DMat) :: NmaxChannelA
    type(NonAntisymmetrizedQNsA), allocatable :: NAS(:)
    integer, private, allocatable :: NASIndex(:)
    integer, private, allocatable :: ASIndex(:)
    integer, private :: n_states_nas = 0
    integer, private :: n_states_as  = 0
  contains
    procedure :: FinNmaxChannelA
    procedure :: InitNmaxChannel4
    procedure :: InitNmaxChannelA
    procedure :: GetNumberNAStatesNmax
    procedure :: GetNumberAStatesNmax
    procedure :: GetNASFromI
    procedure :: GetASIndex
    procedure :: GetNASIndex
    generic :: release => FinNmaxChannelA
    generic :: init => InitNmaxChannel4, InitNmaxChannelA
    generic :: GetNumberNAStates => GetNumberNAStatesNmax
    generic :: GetNumberAStates => GetNumberAStatesNmax
    generic :: GetNAS => GetNASFromI
  end type NmaxChannelA

  type :: AbodyJacIsoChan
    real(8), private :: hw = -1.d0
    integer, private :: j = -1
    integer, private :: p = -100
    integer, private :: t = -1
    integer, private :: Nmax = -1
    integer, private :: n_particle = 0
    ! CFP
    type(NmaxChannelA), allocatable :: cfp(:)
    ! Non-Antisymmetrized states
    integer, private :: n_states_nas = 0
    integer, private, allocatable :: N_nas(:)
    integer, private, allocatable :: alpha_nas(:)
    ! Antisymmetrized states
    integer, private :: n_states_as = 0
    integer, private, allocatable :: N_as(:)
    integer, private, allocatable :: i_as(:)
  contains
    procedure :: FinABodyJacIsoChan
    procedure :: Init4BodyJacIsoChanFrom3
    procedure :: InitABodyJacIsoChan
    procedure :: ReadABodyJacIsoChan
    procedure :: WriteABodyJacIsoChan
    procedure :: GetCFPMatABodyJacIsoChan
    procedure :: GetNmaxChannel
    procedure :: GetFrequencyChan
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetT
    procedure :: GetNmaxChan
    procedure :: GetNumberParticlesChan
    procedure :: GetNumberNAStates
    procedure :: GetNumberAStates
    procedure :: GetNASIndexFromEi
    procedure :: GetASIndexFromEi
    procedure :: GetNASFromEalpha
    procedure :: GetNASFromIndex

    generic :: fin => FinABodyJacIsoChan
    generic :: init => Init4BodyJacIsoChanFrom3, InitABodyJacIsoChan
    generic :: readf => ReadABodyJacIsoChan
    generic :: writef => WriteABodyJacIsoChan
    generic :: GetCFPMat => GetCFPMatABodyJacIsoChan
    generic :: GetFrequency => GetFrequencyChan
    generic :: GetNmax => GetNmaxChan
    generic :: GetNumberParticles => GetNumberParticlesChan
    generic :: GetNASIndex => GetNASIndexFromEi
    generic :: GetASIndex => GetASIndexFromEi
    generic :: GetNAS => GetNASFromEalpha, GetNASFromIndex
  end type AbodyJacIsoChan

  type :: AbodyJacIsoSpace
    real(8), private :: hw = -1.d0
    integer, private :: Jmax = -1
    integer, private :: Jmin = -1
    integer, private :: Tmax = -1
    integer, private :: Tmin = -1
    integer, private :: Nmax = -1
    integer, private :: n_particle = 0
    integer, private :: NChan = 0
    type(AbodyJacIsoChan), allocatable :: jpt(:)
    integer, private, allocatable :: jpt2idx(:,:,:)
  contains
    procedure :: FinABodyJacIsoSpace
    procedure :: Init4BodyJacIsoSpace
    procedure :: InitABodyJacIsoSpace
    procedure :: GetNmaxSpace
    procedure :: GetFrequencySpace
    procedure :: GetJmaxSpace
    procedure :: GetJminSpace
    procedure :: GetTmaxSpace
    procedure :: GetTminSpace
    procedure :: GetNumberParticlesSpace
    procedure :: GetNumberChannels
    procedure :: GetJPTIndex
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromJPT

    generic :: fin => FinABodyJacIsoSpace
    generic :: init => Init4BodyJacIsoSpace, InitABodyJacIsoSpace
    generic :: GetNmax => GetNmaxSpace
    generic :: GetFrequency => GetFrequencySpace
    generic :: GetJmax => GetJmaxSpace
    generic :: GetJmin => GetJminSpace
    generic :: GetTmax => GetTmaxSpace
    generic :: GetTmin => GetTminSpace
    generic :: GetNumberParticles => GetNumberParticlesSpace
    generic :: GetIndex => GetJPTIndex
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromJPT
  end type AbodyJacIsoSpace

  type(SixJsStore), private :: sixj_3l, sixj_2l1j, sixj_isospin
  type(NineJsStore), private :: ninej
  type(TMbracketStore), private :: tmbk
  type(sys), private :: sy
contains

  function GetNmaxSpace(this) result(Nmax)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmaxSpace

  function GetFrequencySpace(this) result(hw)
    class(ABodyJacIsoSpace), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequencySpace

  function GetJmaxSpace(this) result(Jmax)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmaxSpace

  function GetJminSpace(this) result(Jmin)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: Jmin
    Jmin = this%Jmin
  end function GetJminSpace

  function GetTmaxSpace(this) result(Tmax)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: Tmax
    Tmax = this%Tmax
  end function GetTmaxSpace

  function GetTminSpace(this) result(Tmin)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: Tmin
    Tmin = this%Tmin
  end function GetTminSpace

  function GetNumberParticlesSpace(this) result(n_particle)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: n_particle
    n_particle = this%n_particle
  end function GetNumberParticlesSpace

  function GetNumberChannels(this) result(NChan)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetJPTIndex(this, j, p, t) result(idx)
    class(ABodyJacIsoSpace), intent(in) :: this
    integer, intent(in) :: j, p, t
    integer :: idx
    idx = this%jpt2idx(j,p,t)
  end function GetJPTIndex

  function GetChannelFromIndex(this, idx) result(r)
    class(ABodyJacIsoSpace), intent(in), target :: this
    integer :: idx
    type(ABodyJacIsoChan), pointer :: r
    r => null()
    if(idx == 0) return
    r => this%jpt(idx)
  end function GetChannelFromIndex

  function GetChannelFromJPT(this, j, p, t) result(r)
    class(ABodyJacIsoSpace), intent(in), target :: this
    integer, intent(in) :: j, p, t
    type(ABodyJacIsoChan), pointer :: r
    r => this%GetChannel( this%GetIndex(j,p,t) )
  end function GetChannelFromJPT

  function GetNmaxChannel(this,Nmax) result(r)
    class(ABodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: Nmax
    type(NmaxChannelA), pointer :: r
    r => this%cfp(Nmax)
  end function GetNmaxChannel

  function GetFrequencyChan(this) result(hw)
    class(ABodyJacIsoChan), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequencyChan

  function GetJ(this) result(j)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetT(this) result(t)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: t
    t = this%t
  end function GetT

  function GetNmaxChan(this) result(Nmax)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmaxChan

  function GetNumberParticlesChan(this) result(np)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: np
    np = this%n_particle
  end function GetNumberParticlesChan

  function GetNumberNAStates(this) result(n_states)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states_nas
  end function GetNumberNAStates

  function GetNumberAStates(this) result(n_states)
    class(ABodyJacIsoChan), intent(in) :: this
    integer :: n_states
    n_states = this%n_states_as
  end function GetNumberAStates

  function GetNASIndexFromEi(this, E, i) result(idx)
    class(ABodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: E, i
    integer :: idx
    idx = this%cfp(E)%GetNASIndex(i)
  end function GetNASIndexFromEi

  function GetASIndexFromEi(this, E, i) result(idx)
    class(ABodyJacIsoChan), intent(in) :: this
    integer, intent(in) :: E, i
    integer :: idx
    idx = this%cfp(E)%GetASIndex(i)
  end function GetASIndexFromEi

  function GetNASFromEalpha(this,E,alpha) result(r)
    class(ABodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: E, alpha
    type(NonAntisymmetrizedQNsA), pointer :: r
    r => this%cfp(E)%NAS(alpha)
  end function GetNASFromEalpha

  function GetNASFromIndex(this,i) result(r)
    class(ABodyJacIsoChan), intent(in), target :: this
    integer, intent(in) :: i
    type(NonAntisymmetrizedQNsA), pointer :: r
    r => this%GetNAS( this%N_nas(i), this%alpha_nas(i) )
  end function GetNASFromIndex

  function GetNumberNAStatesNmax(this) result(n_states)
    class(NmaxChannelA), intent(in) :: this
    integer :: n_states
    n_states = this%n_states_nas
  end function GetNumberNAStatesNmax

  function GetNumberAStatesNmax(this) result(n_states)
    class(NmaxChannelA), intent(in) :: this
    integer :: n_states
    n_states = this%n_states_as
  end function GetNumberAStatesNmax

  function GetNASFromI(this,i) result(r)
    class(NmaxChannelA), intent(in), target :: this
    integer, intent(in) :: i
    type(NonAntisymmetrizedQNsA), pointer :: r
    r => null()
    if(i == 0) return
    r => this%NAS(i)
  end function GetNASFromI

  function GetNASIndex(this, i) result(idx)
    class(NmaxChannelA), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%NASIndex(i)
  end function GetNASIndex

  function GetASIndex(this, i) result(idx)
    class(NmaxChannelA), intent(in) :: this
    integer, intent(in) :: i
    integer :: idx
    idx = this%ASIndex(i)
  end function GetASIndex

  subroutine FinAbodyJacIsoSpace(this)
    class(AbodyJacIsoSpace), intent(inout) :: this
    integer :: ich
    do ich = 1, this%GetNumberChannels()
      call this%jpt(ich)%fin()
    end do
    deallocate(this%jpt)
    deallocate(this%jpt2idx)
    this%hw = -1.d0
    this%Jmax = -1
    this%Jmin = -1
    this%Tmax = -1
    this%Tmin = -1
    this%Nmax = -1
    this%n_particle = 0
    this%NChan = 0
  end subroutine FinAbodyJacIsoSpace

  subroutine Init4bodyJacIsoSpace(this, Nmax, jacobi3, path_to_dir)
    class(ABodyJacIsoSpace), intent(inout) :: this
    integer, intent(in) :: Nmax
    type(str), intent(in) :: path_to_dir
    type(ThreeBodyJacIsoSpace), intent(in) :: jacobi3
    integer :: j, p, t, ch

    if(allocated(this%jpt)) call this%fin()

    this%n_particle = 4
    this%Nmax = Nmax
    this%Jmax = 2*Nmax + 4
    this%Tmax = 4
    this%Jmin = 0
    this%Tmin = 0
    this%hw = jacobi3%GetFrequency()

    this%NChan = (this%Jmax/2 + 1) * (this%Tmax/2 + 1) * 2
    allocate(this%jpt(this%NChan))
    allocate(this%jpt2idx(this%Jmin:this%Jmax,-1:1,this%Tmin:this%Tmax))
    this%jpt2idx(:,:,:) = 0
    ch = 0
    do j = this%Jmin, this%Jmax, 2
      do p = 1, -1, -2
        do t = this%Tmin, this%Tmax, 2
          ch = ch + 1
          this%jpt2idx(j,p,t) = ch
          call this%jpt(ch)%init(j, p, t, Nmax, jacobi3, path_to_dir)
        end do
      end do
    end do
  end subroutine Init4bodyJacIsoSpace

  subroutine InitAbodyJacIsoSpace(this, Nmax, n_p, jacobiA_1, path_to_dir)
    class(ABodyJacIsoSpace), intent(inout) :: this
    integer, intent(in) :: Nmax, n_p
    type(str), intent(in) :: path_to_dir
    type(ABodyJacIsoSpace), intent(in) :: jacobiA_1
    integer :: j, p, t, ch

    if(allocated(this%jpt)) call this%fin()

    this%Nmax = Nmax
    this%n_particle = n_p
    this%Jmax = 2*Nmax + n_p
    this%Tmax = n_p
    this%hw = jacobiA_1%hw

    if(mod(n_p, 2) == 0) then
      this%Jmin = 0
      this%Tmin = 0
      this%NChan = (this%Jmax/2 + 1) * (this%Tmax/2 + 1) * 2
    end if

    if(mod(n_p, 2) == 1) then
      this%Jmin = 1
      this%Tmin = 1
      this%NChan = (this%Jmax + 1) * (this%Tmax + 1) / 2
    end if
    allocate(this%jpt(this%NChan))
    allocate(this%jpt2idx(this%Jmin:this%Jmax,-1:1,this%Tmin:this%Tmax))
    this%jpt2idx(:,:,:) = 0
    ch = 0
    do j = this%Jmin, this%Jmax, 2
      do p = 1, -1, -2
        do t = this%Tmin, this%Tmax, 2
          ch = ch + 1
          this%jpt2idx(j,p,t) = ch
          call this%jpt(ch)%init(j, p, t, Nmax, n_p, jacobiA_1, path_to_dir)
        end do
      end do
    end do
  end subroutine InitAbodyJacIsoSpace

  function get_file_name(n_particle, j_in, p, t_in, Nmax, path_to_dir) result(f)
    integer, intent(in) :: n_particle, j_in, p, t_in, Nmax
    type(str), intent(in) :: path_to_dir
    integer :: j, t
    type(str) :: f
    type(str) :: dir, pari
    type(sys) :: s

    j = j_in
    t = t_in
    if(mod(n_particle, 2) == 0) then
      j = j_in / 2
      t = t_in / 2
    end if

    dir = path_to_dir
    call s%mkdir(dir%val)
    if(p == -1) pari = '-'
    if(p ==  1) pari = '+'
    f = dir + s%str('/cfp/cfp_A') + s%str(n_particle) + s%str('_j') + s%str(j) + s%str('p') + &
        & pari + s%str('t') + s%str(t) + s%str('_Nmax') + s%str(Nmax) + s%str('.bin')
  end function get_file_name

  subroutine FinAbodyJacIsoChan(this)
    class(AbodyJacIsoChan), intent(inout) :: this
    integer :: NAmax
    do NAmax = 0, this%Nmax
      call this%cfp(NAmax)%release()
    end do
    deallocate(this%cfp)
    if( this%GetNumberNAStates()*this%GetNumberAStates() == 0) return
    deallocate(this%N_nas)
    deallocate(this%alpha_nas)
    deallocate(this%N_as)
    deallocate(this%i_as)
    this%n_states_nas = 0
    this%n_states_as = 0
    this%hw = -1.d0
    this%j = -1
    this%p = -100
    this%t = -1
    this%Nmax = -1
    this%n_particle = 0
  end subroutine FinAbodyJacIsoChan

  subroutine Init4BodyJacIsoChanFrom3(this, j, p, t, Nmax, jac, path_to_dir)
    class(ABodyJacIsoChan), intent(inout) :: this
    integer, intent(in) :: j, p, t, Nmax
    type(str), intent(in) :: path_to_dir
    type(ThreeBodyJacIsoSpace), intent(in) :: jac
    integer :: N
    real(8) :: ti
    type(str) :: f
    type(sys) :: s
    integer :: runit=21, wunit = 20

    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%t = t
    this%n_particle = 4
    this%hw = jac%GetFrequency()
    f = get_file_name(this%n_particle, j, p, t, Nmax, path_to_dir)
    if(s%isfile(f)) then
      write(*,"(a,i3,a,i2,a,i2,a,i2,a)") "# Reading cfp...  A= 4, Nmax= ", Nmax, ", J=", j, "/2, P=", p, ", T=", t,  "/2"
      open(runit, file=f%val, form="unformatted", access="stream")
      call this%readf(this%GetFrequency(), runit, Nmax)
      close(runit)
      return
    end if

    write(*,"(a,i3,a,i2,a,i2,a,i2,a)") "# Calculating cfp...  A= 4, Nmax= ", Nmax, ", J=", j, "/2, P=", p, ", T=", t,  "/2"
    call sixj_3l%init(0,2*Nmax, .false., 0,2*Nmax, .false., 0,2*Nmax, .false.)
    call sixj_2l1j%init(0,2*Nmax,.false., 0,2*Nmax,.false., 1,2*Nmax+1, .true.)
    call sixj_isospin%init(1,1,.true., 0,2,.false., 1,1, .true.)
    !call ninej%init(0,2*Nmax+2,.false., 1,2*Nmax+1, .true., 1,2*Nmax+1,.true., 0,2*Nmax,.false., jdmin_in=j, jdmax_in=j)
    call ninej%init(0,2*Nmax+2,.false., 1,2*Nmax+3, .true., 1,2*Nmax+3,.true., j,j,.false.)
    call tmbk%init(Nmax, 1.d0 / 8.d0)
    !write(*,"(a,f12.6,a)") "SixJs is using: ", sixj_3l%GetMemory()+ sixj_2l1j%GetMemory(), " GB"
    !write(*,"(a,f12.6,a)") "NineJ is using:  ", ninej%GetMemory(), " GB"
    !write(*,"(a,f12.6,a)") "TMbk is using:  ", tmbk%GetMemory(), " GB"

    ti = omp_get_wtime()
    allocate(this%cfp(0:this%Nmax))
    do N = 0, this%Nmax
      call this%cfp(N)%init(j, p, t, N, jac)
    end do
    call set_indices(this)
    call timer%add(sy%str("Init4BodyJacIsoChanFrom3"), omp_get_wtime() - ti)
    write(*,"(a,i3,a,i3,a,i2,a,i2,a,2i8)") "A=", this%n_particle, ", J=", &
        & this%j, "/2, P=", this%p, ", T=", this%t, &
        & "/2, Number of physical and orthonormal states: ", &
        & this%GetNumberNAStates(), this%GetNumberAStates()

    call tmbk%fin()
    call sixj_3l%fin()
    call sixj_2l1j%fin()
    call sixj_isospin%fin()
    call ninej%fin()
    open(wunit, file=f%val, form="unformatted",access="stream")
    call this%writef(wunit)
    close(wunit)
  end subroutine Init4BodyJacIsoChanFrom3

  subroutine InitABodyJacIsoChan(this, j, p, t, Nmax, n_p, jac, path_to_dir)
    ! j and t are doubled
    class(ABodyJacIsoChan), intent(inout) :: this
    integer, intent(in) :: j, p, t, Nmax, n_p
    type(ABodyJacIsoSpace), intent(in) :: jac
    type(str), intent(in) :: path_to_dir
    integer :: N
    real(8) :: ti
    type(str) :: f
    type(sys) :: s
    integer :: runit=21, wunit = 20

    this%Nmax = Nmax
    this%j = j
    this%p = p
    this%t = t
    this%n_particle = n_p
    this%hw = jac%GetFrequency()
    f = get_file_name(this%n_particle, j, p, t, Nmax, path_to_dir)
    if(s%isfile(f)) then
      write(*,"(a,i2,a,i3,a,i2,a,i2,a,i2,a)") &
          &"# Reading cfp...  A=", n_p,", Nmax= ", Nmax, ", J=", j, "/2, P=", p, ", T=", t,  "/2"
      open(runit, file=f%val, form="unformatted",access="stream")
      call this%readf(this%GetFrequency(), runit, Nmax)
      close(runit)
      return
    end if

    write(*,"(a,i2,a,i3,a,i2,a,i2,a,i2,a)") &
        &"# Calculating cfp...  A=", n_p,", Nmax= ", Nmax, ", J=", j, "/2, P=", p, ", T=", t,  "/2"
    call sixj_3l%init(0,2*Nmax, .false., 0,2*Nmax, .false., 0,2*Nmax, .false.)
    call sixj_2l1j%init(0,2*Nmax,.false., 0,2*Nmax,.false., 1,2*Nmax+1, .true.)
    if(mod(n_p, 2) == 0) then
      call sixj_isospin%init(1,1,.true., 0,n_p-2,.false., 1,1, .true.)
      call ninej%init(0,2*Nmax+n_p-2,.false., 1,2*Nmax+1, .true., 1,2*Nmax+1,.true., 0,2*Nmax,.false., jdmin_in=j, jdmax_in=j)
      !call ninej%init(0,2*Nmax+n_p-2,.false., 1,2*Nmax+n_p-1, .true., 1,2*Nmax+n_p-1,.true., j,j,.false.)
    end if

    if(mod(n_p, 2) == 1) then
      call sixj_isospin%init(1,1,.true., 1,n_p-2,.true., 1,1, .true.)
      call ninej%init(1,2*Nmax+n_p-2,.true., 1,2*Nmax+1, .true., 1,2*Nmax+1,.true., 0,2*Nmax,.false., jdmin_in=j, jdmax_in=j)
      !call ninej%init(1,2*Nmax+n_p-2,.true., 0,2*Nmax+n_p-1, .false., 0,2*Nmax+n_p-1,.false., j,j,.true.)
    end if
    call tmbk%init(Nmax, 1.d0 / dble(n_p * (n_p - 2)))

    ti = omp_get_wtime()
    allocate(this%cfp(0:this%Nmax))
    do N = 0, this%Nmax
      call this%cfp(N)%init(j, p, t, N, n_p, jac)
    end do
    call set_indices(this)
    call timer%add(sy%str("InitABodyJacIsoChan"), omp_get_wtime() - ti)
    write(*,"(a,i3,a,i3,a,i2,a,i2,a,2i8)") "A=", this%n_particle, ", J=", &
        & this%j, "/2, P=", this%p, ", T=", this%t, &
        & "/2, Number of physical and orthonormal states: ", &
        & this%GetNumberNAStates(), this%GetNumberAStates()

    call tmbk%fin()
    call sixj_3l%fin()
    call sixj_2l1j%fin()
    call sixj_isospin%fin()
    call ninej%fin()
    open(wunit, file=f%val, form="unformatted",access="stream")
    call this%writef(wunit)
    close(wunit)
  end subroutine InitABodyJacIsoChan

  subroutine WriteAbodyJacIsoChan(this, iunit)
    use Profiler, only: timer
    class(AbodyJacIsoChan), intent(inout) :: this
    integer, intent(in) :: iunit
    type(NmaxChannelA), pointer :: nchjac
    type(NonAntisymmetrizedQNsA), pointer :: nas
    integer :: nAmax, i
    real(8) :: ti

    ti = omp_get_wtime()
    write(iunit) this%GetJ(), this%GetParity(), this%GetT(), this%GetNmax()
    do nAmax = 0, this%GetNmax()
      nchjac => this%GetNmaxChannel(nAmax)
      write(iunit) nchjac%GetNumberNAStates(), nchjac%GetNumberAStates(), nAmax
      if( nchjac%GetNumberNAStates() < 1) cycle
      do i = 1, nchjac%GetNumberNAStates()
        nas => nchjac%GetNAS(i)
        write(iunit) nas%Getna()
        write(iunit) nas%Getia()
        write(iunit) nas%Getja()
        write(iunit) nas%Getta()
        write(iunit) nas%Getnb()
        write(iunit) nas%Getlb()
        write(iunit) nas%Getjb()
      end do
      if( nchjac%GetNumberAStates() < 1) cycle
      do i = 1, nchjac%GetNumberAStates()
        write(iunit) nchjac%m(i,:)
      end do
    end do
    call timer%Add(sy%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteAbodyJacIsoChan

  subroutine ReadAbodyJacIsoChan(this, hw, iunit, Nmax)
    use Profiler, only: timer
    class(AbodyJacIsoChan), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: iunit
    integer, intent(in), optional :: Nmax
    type(NmaxChannelA), pointer :: nchjac
    integer :: nAmax, i, Nmax_read, nAmax_read
    integer :: nphys, north
    integer :: na, ia, ja, ta, nb, lb, jb
    real(8) :: ti

    if(allocated(this%cfp)) call this%fin()
    ti = omp_get_wtime()
    this%hw = hw
    read(iunit) this%j, this%p, this%t, Nmax_read

    if(present(Nmax)) then; this%Nmax = Nmax
    else; this%Nmax = Nmax_read; end if

    allocate(this%cfp(0:this%GetNmax() ))

    do nAmax = 0, this%GetNmax()
      nchjac => this%GetNmaxChannel(nAmax)
      read(iunit) nphys, north, nAmax_read
      if(nAmax /= nAmax_read) cycle
      nchjac%n_states_nas = nphys
      nchjac%n_states_as  = north
      if( nphys < 1) cycle
      allocate( nchjac%NAS(nphys) )
      do i = 1, nchjac%GetNumberNAStates()
        read(iunit) na,ia,ja,ta,nb,lb,jb
        call nchjac%NAS(i)%set(na,ia,ja,ta,nb,lb,jb)
      end do

      if( north < 1) cycle
      call nchjac%ini( north, nphys )
      do i = 1, nchjac%GetNumberAStates()
        read(iunit) nchjac%m(i,:)
      end do
    end do
    call timer%Add(sy%str('Read from file'), omp_get_wtime() - ti)
    call set_indices(this)
  end subroutine ReadAbodyJacIsoChan

  function GetCFPMatABodyJacIsoChan(this) result(cfp)
    class(ABodyJacIsoChan), intent(in) :: this
    type(DMat) :: cfp
    integer :: north, nphys, bra, ket, Nmax, l, i
    north = this%GetNumberAStates()
    nphys = this%GetNumberNAStates()
    call cfp%zeros(nphys,north)
    !$omp parallel
    !$omp do private(bra, Nmax, i, ket, l)
    do bra = 1, north
      Nmax = this%N_as(bra)
      i    = this%i_as(bra)
      do ket = 1, nphys
        if(Nmax /= this%N_nas(ket)) cycle
        l = this%alpha_nas(ket)
        cfp%m(ket,bra) = this%cfp(Nmax)%m(i,l)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end function GetCFPMatABodyJacIsoChan

  subroutine set_indices(this)
    class(ABodyJacIsoChan), intent(inout) :: this
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

  subroutine FinNmaxChannelA(this)
    class(NmaxChannelA), intent(inout) :: this
    if(this%GetNumberAStates() * this%GetNumberAStates() == 0) return
    call this%fin()
    deallocate(this%NAS)
    deallocate(this%NASIndex)
    deallocate(this%ASIndex)
    this%n_states_nas = 0
    this%n_states_as  = 0
  end subroutine FinNmaxChannelA

  subroutine InitNmaxChannel4(this, j, p, t, N, jspace)
    class(NmaxChannelA), intent(inout) :: this
    integer, intent(in) :: j, p, t, N
    type(ThreeBodyJacIsoSpace), intent(in) :: jspace
    call set_jacobi_nas4(this,j,p,t,N,jspace,this%n_states_nas,.false.)
    if(this%GetNumberNAStates() == 0) return
    call set_jacobi_nas4(this,j,p,t,N,jspace,this%n_states_nas,.true.)

    this%n_states_as = get_dim_as4(this, j, t, jspace)
    if(this%GetNumberAStates() == 0) return
    call this%ini( this%GetNumberAStates(), this%GetNumberNAStates() )
    call set_cfp4(this, j, t, jspace)
  end subroutine InitNmaxChannel4

  subroutine InitNmaxChannelA(this, j, p, t, N, n_p, jspace)
    class(NmaxChannelA), intent(inout) :: this
    integer, intent(in) :: j, p, t, N, n_p
    type(ABodyJacIsoSpace), intent(in) :: jspace
    call set_jacobi_nasA(this,j,p,t,N,this%n_states_nas,jspace,.false.)
    if(this%GetNumberNAStates() == 0) return
    call set_jacobi_nasA(this,j,p,t,N,this%n_states_nas,jspace,.true.)

    this%n_states_as = get_dim_asA(this, j, t, jspace, n_p)
    if(this%GetNumberAStates() == 0) return
    call this%ini( this%GetNumberAStates(), this%GetNumberNAStates() )
    call set_cfpA(this, j, t, jspace, n_p)
  end subroutine InitNmaxChannelA

  subroutine set_jacobi_nas4(this, j, p, t, N, jspace, n_states_nas, set_mode)
    use MyLibrary, only: triag
    class(NmaxChannelA), intent(inout) :: this
    integer, intent(in) :: j, p, t, N
    integer, intent(inout) :: n_states_nas
    type(ThreeBodyJacIsoSpace), intent(in) :: jspace
    logical, intent(in) :: set_mode
    type(ThreeBodyJacIsoChan), pointer :: chjac
    integer :: n_states, i
    integer :: na, pa, ia, ja, ta, nb, lb, jb

    if(set_mode) allocate( this%NAS( this%GetNumberNAStates() ))
    n_states = 0
    do ja = 1, 2*N+3, 2
      do pa = 1, -1, -2
        do ta = 1, 3, 2
          if(triag(ta, 1, t)) cycle
          chjac => jspace%GetChannel(ja,pa,ta)
          if(.not. associated(chjac) ) cycle

          do i = 1, chjac%GetNumberAStates()
            na = chjac%GetEAS(i)
            ia = chjac%GetiAS(i)

            do nb = 0, (N-na)/2
              lb = N - na - 2*nb
              if(lb < 0) cycle
              if(pa * (-1)**lb /= p) cycle
              do jb = abs(2*lb-1), 2*lb+1, 2
                if(triag(ja,jb,j)) cycle
                n_states = n_states + 1
                if(set_mode) call this%NAS(n_states)%set(na,ia,ja,ta,nb,lb,jb)
              end do
            end do
          end do

        end do
      end do
    end do
    if(.not. set_mode) n_states_nas = n_states
  end subroutine set_jacobi_nas4

  function get_dim_as4(this, j, t, jspace) result(n_states)
    class(NmaxChannelA), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), intent(in) :: jspace
    integer, intent(in) :: j, t
    integer :: n_states
    type(NonAntisymmetrizedQNsA), pointer :: nas
    real(8), allocatable :: diag(:)
    integer :: n_nas
    allocate( diag(this%GetNumberNAStates() ))
    diag(:) = 0.d0
    do n_nas = 1, this%GetNumberNAStates()
      nas => this%GetNAS(n_nas)
      diag(n_nas) = anti_op4_wrap(nas,nas,jspace,j,t)
    end do
    n_states = int(sum(diag) + 1.d-2)
    if(abs(sum(diag) - dble(n_states)) > 1.d-4) then
      write(*,*) "error in trace of A operator: tr(A) = ", sum(diag)
      stop
    end if
    deallocate(diag)
  end function get_dim_as4

  subroutine set_cfp4(this, j, t, jspace)
    class(NmaxChannelA), intent(inout) :: this
    type(ThreeBodyJacIsoSpace), intent(in) :: jspace
    integer, intent(in) :: j, t
    integer :: ibra, iket, i, k, ii
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: anti
    type(DMat) :: A
    type(EigenSolSymD) :: sol

    call A%ini( this%GetNumberNAStates(), this%GetNumberNAStates() )
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, anti) schedule(dynamic)
    do ibra = 1, this%GetNumberNAStates()
      bra => this%GetNAS(ibra)
      do iket = 1, ibra
        ket => this%GetNAS(iket)
        anti = anti_op4_wrap(bra, ket, jspace, j, t)
        A%m(ibra,iket) = anti
        A%m(iket,ibra) = anti
      end do
    end do
    !$omp end do
    !$omp end parallel
    call sol%init(A)
    call sol%DiagSym(A)
    k = 0
    i = 0
    do ii = 1, this%GetNumberNAStates()
      if(abs(1.d0 - sol%eig%v(ii)) < 1.d-4) then
        k = k + 1
        this%m(k,:) = sol%vec%m(:,ii)
      else
        if(abs(sol%eig%v(ii)) > 1.d-4) then
          i = i + 1
        end if
      end if
    end do
    if( i /= 0 ) then
      write(*,'(a)') 'Warning: cfp', __LINE__,__FILE__
      call sol%eig%prt()
    end if
    call sol%fin()
    call A%fin()
  end subroutine set_cfp4

  subroutine set_jacobi_nasA(this, j, p, t, N, n_states_nas, jac, set_mode)
    use MyLibrary, only: triag
    class(NmaxChannelA), intent(inout) :: this
    integer, intent(in) :: j, p, t, N
    integer, intent(inout) :: n_states_nas
    type(ABodyJacIsoSpace), intent(in) :: jac
    logical, intent(in) :: set_mode
    type(ABodyJacIsoChan), pointer :: chjac
    integer :: n_states, i
    integer :: na, pa, ia, ja, ta, nb, lb, jb

    if(set_mode) allocate( this%NAS( this%GetNumberNAStates() ))
    n_states = 0
    do ja = jac%jmin, jac%jmax, 2
      do pa = 1, -1, -2
        do ta = jac%tmin, jac%tmax, 2
          if(triag(ta, 1, t)) cycle
          chjac => jac%GetChannel(ja,pa,ta)
          if(.not. associated(chjac) ) cycle

          do i = 1, chjac%GetNumberAStates()
            na = chjac%N_as(i)
            ia = chjac%i_as(i)

            do nb = 0, (N-na)/2
              lb = N - na - 2*nb
              if(lb < 0) cycle
              if(pa * (-1)**lb /= p) cycle
              do jb = abs(2*lb-1), 2*lb+1, 2
                if(triag(ja,jb,j)) cycle
                n_states = n_states + 1
                if(set_mode) call this%NAS(n_states)%set(na,ia,ja,ta,nb,lb,jb)
              end do
            end do
          end do

        end do
      end do
    end do
    if(.not. set_mode) n_states_nas = n_states
  end subroutine set_jacobi_nasA

  function get_dim_asA(this, j, t, jspace, np) result(n_states)
    class(NmaxChannelA), intent(inout) :: this
    type(ABodyJacIsoSpace), intent(in) :: jspace
    integer, intent(in) :: j, t, np
    integer :: n_states
    type(NonAntisymmetrizedQNsA), pointer :: nas
    real(8), allocatable :: diag(:)
    integer :: n_nas
    allocate( diag(this%GetNumberNAStates() ))
    diag(:) = 0.d0
    do n_nas = 1, this%GetNumberNAStates()
      nas => this%GetNAS(n_nas)
      diag(n_nas) = anti_op_wrap(nas,nas,jspace,j,t, np)
    end do
    n_states = int(sum(diag) + 1.d-2)
    if(abs(sum(diag) - dble(n_states)) > 1.d-4) then
      write(*,*) "error in trace of A operator: tr(A) = ", sum(diag)
      stop
    end if
    deallocate(diag)
  end function get_dim_asA

  subroutine set_cfpA(this, j, t, jspace, np)
    class(NmaxChannelA), intent(inout) :: this
    type(ABodyJacIsoSpace), intent(in) :: jspace
    integer, intent(in) :: j, t, np
    integer :: ibra, iket, i, k, ii
    type(NonAntisymmetrizedQNsA), pointer :: bra, ket
    real(8) :: anti
    type(DMat) :: A
    type(EigenSolSymD) :: sol

    call A%ini( this%GetNumberNAStates(), this%GetNumberNAStates() )
    !$omp parallel
    !$omp do private(ibra, bra, iket, ket, anti) schedule(dynamic)
    do ibra = 1, this%GetNumberNAStates()
      bra => this%GetNAS(ibra)
      do iket = 1, ibra
        ket => this%GetNAS(iket)
        anti = anti_op_wrap(bra, ket, jspace, j, t, np)
        A%m(ibra,iket) = anti
        A%m(iket,ibra) = anti
      end do
    end do
    !$omp end do
    !$omp end parallel
    call sol%init(A)
    call sol%DiagSym(A)
    k = 0
    i = 0
    do ii = 1, this%GetNumberNAStates()
      if(abs(1.d0 - sol%eig%v(ii)) < 1.d-4) then
        k = k + 1
        this%m(k,:) = sol%vec%m(:,ii)
      else
        if(abs(sol%eig%v(ii)) > 1.d-4) then
          i = i + 1
        end if
      end if
    end do
    if( i /= 0 ) then
      write(*,'(a)') 'Warning: cfp', __LINE__,__FILE__
      call sol%eig%prt()
    end if
    call sol%fin()
    call A%fin()
  end subroutine set_cfpA

  function anti_op4_wrap(bra, ket, jspace, j, t) result(anti)
    type(NonAntisymmetrizedQNsA), intent(in) :: bra, ket
    type(ThreeBodyJacIsoSpace), intent(in) :: jspace
    type(ThreeBodyJacIsoChan), pointer :: chbra, chket
    integer, intent(in) :: j, t
    real(8) :: anti
    integer :: Nabra, iabra, jabra, tabra, nbbra, lbbra, jbbra
    integer :: Naket, iaket, jaket, taket, nbket, lbket, jbket
    Nabra = bra%Getna()
    iabra = bra%Getia()
    jabra = bra%Getja()
    tabra = bra%Getta()
    nbbra = bra%Getnb()
    lbbra = bra%Getlb()
    jbbra = bra%Getjb()

    Naket = ket%Getna()
    iaket = ket%Getia()
    jaket = ket%Getja()
    taket = ket%Getta()
    nbket = ket%Getnb()
    lbket = ket%Getlb()
    jbket = ket%Getjb()
    chbra => jspace%GetChannel(jabra, (-1)**nabra, tabra)
    chket => jspace%GetChannel(jaket, (-1)**naket, taket)
    anti = Aop4( chbra%GetNmaxChannel(nabra), [Nabra, iabra, jabra, tabra, nbbra, lbbra, jbbra], &
        &        chket%GetNmaxChannel(naket), [Naket, iaket, jaket, taket, nbket, lbket, jbket], j, t)
  end function anti_op4_wrap

  real(8) function Aop4(jacb, bra, jack, ket, jtot, ttot) result(anti)
    use MyLibrary, only: hat
    integer, intent(in) :: bra(7), ket(7)
    type(NmaxChannelIso3), intent(in) :: jacb, jack
    type(NonAntisymmetrizedIsoQNs3), pointer :: bra3, ket3
    integer, intent(in) :: jtot, ttot
    integer :: j3L, t3L, nmax3L, iadqnL, n4L, l4L, j4L
    integer :: j3R, t3R, nmax3R, iadqnR, n4R, l4R, j4R
    real(8) :: ex, ex1, ex2, gms
    integer :: n, l, s, j, t, jindL, jindR
    integer :: nitL, litL, jit3L, nitR, litR, jit3R, K, Lp
    integer :: kmin, kmax, Lpmin, Lpmax

    Nmax3L = bra(1); Nmax3R = ket(1)
    iadqnL = bra(2); iadqnR = ket(2)
    j3l    = bra(3); j3r    = ket(3)
    t3l    = bra(4); t3r    = ket(4)
    n4l    = bra(5); n4r    = ket(5)
    l4l    = bra(6); l4r    = ket(6)
    j4l    = bra(7); j4r    = ket(7)
    ex = 0.d0
    do jindL = 1, jacb%GetNumberNAStates()
      bra3 => jacb%GetNAS(jindL)
      n = bra3%GetN12()
      l = bra3%GetL12()
      s = bra3%GetS12()
      j = bra3%GetJ12()
      t = bra3%GetT12()
      nitL = bra3%GetN3()
      litL = bra3%GetL3()
      jit3L = bra3%GetJ3()
      do jindR = 1, jack%GetNumberNAStates()
        ket3 => jack%GetNAS(jindR)
        if(n /= ket3%GetN12() ) cycle
        if(l /= ket3%GetL12() ) cycle
        if(s /= ket3%GetS12() ) cycle
        if(j /= ket3%GetJ12() ) cycle
        if(t /= ket3%GetT12() ) cycle
        nitR = ket3%GetN3()
        litR = ket3%GetL3()
        jit3R = ket3%GetJ3()

        kmin = max(abs(jit3L - j4R) / 2, abs(jit3R - j4L) / 2, abs(litR - l4L), abs(litL - l4R))
        kmax = min((jit3L + j4R) / 2, (jit3R + j4L) / 2, litR + l4L, litL + l4R)
        Lpmin = max(abs(litL - l4L), abs(LitR - l4R))
        Lpmax = min(litL + l4L, LitR + l4R)
        if(2*n4l + l4l + 2*NitL + LitL /= 2*NitR + LitR + 2*n4R + l4R) cycle

        ex1 = 0.d0
        do Lp = Lpmin, Lpmax
          ex2 = 0.d0
          do K = kmin, kmax
            ex2 = ex2 + &
                & sixj_3l%get(2*l4l, 2*litR, 2*K, 2*l4R, 2*LitL, 2*Lp) * &
                !& ninej%get(2*j, jit3L, j3L, jit3R, 2*K, j4L, j3R, j4R, jtot) * &
                & ninej%get(2*j, j3L, jit3L, j3R, Jtot, j4R, jit3R, j4L, 2*K) * &
                & dble(2*K+1) * &
                & sixj_2l1j%get(2*litR, 2*l4L, 2*K, j4L, jit3R, 1) * &
                & sixj_2l1j%get(2*litL, 2*l4R, 2*K, j4R, jit3L, 1)
          end do
          gms = tmbk%get(NitL,LitL,n4L,l4L, NitR,LitR,n4R,l4R, Lp)
          ex1 = ex1 +  dble(2*Lp+1) * gms * ex2
        end do
        ex = ex + (-1.d0) ** ((t3L + t3R + jit3L + jit3R)/2 + LitL + LitR) * &
          & hat(jit3L) * hat(jit3R) * hat(j4L) * hat(j4R) * hat(j3L) * hat(j3R) * &
          & sixj_isospin%get(1, 2*t, t3R, 1, ttot, t3L) * &
          & hat(t3L) * hat(t3R) * &
          & jacb%m(iadqnL, jindL) * &
          & jack%m(iadqnR, jindR) * ex1
      end do
    end do
    anti = - 3.d0 * ex / 4.d0
    if(j3L == j3R .and. t3L == t3R .and. nmax3L == nmax3R &
      &           .and. iadqnL == iadqnR .and. n4L == n4R .and. l4L == l4R &
      &           .and. j4L == j4R) then
      anti = (1.d0 - 3.d0 * ex) / 4.d0
    end if
  end function Aop4

  function anti_op_wrap(bra, ket, jspace, j, t, np) result(anti)
    type(NonAntisymmetrizedQNsA), intent(in) :: bra, ket
    type(ABodyJacIsoSpace), intent(in) :: jspace
    type(ABodyJacIsoChan), pointer :: chbra, chket
    integer, intent(in) :: j, t, np
    real(8) :: anti
    integer :: Nabra, iabra, jabra, tabra, nbbra, lbbra, jbbra
    integer :: Naket, iaket, jaket, taket, nbket, lbket, jbket
    Nabra = bra%Getna()
    iabra = bra%Getia()
    jabra = bra%Getja()
    tabra = bra%Getta()
    nbbra = bra%Getnb()
    lbbra = bra%Getlb()
    jbbra = bra%Getjb()

    Naket = ket%Getna()
    iaket = ket%Getia()
    jaket = ket%Getja()
    taket = ket%Getta()
    nbket = ket%Getnb()
    lbket = ket%Getlb()
    jbket = ket%Getjb()
    chbra => jspace%GetChannel(jabra, (-1)**nabra, tabra)
    chket => jspace%GetChannel(jaket, (-1)**naket, taket)
    anti = Aop( chbra%GetNmaxChannel(nabra), [Nabra, iabra, jabra, tabra, nbbra, lbbra, jbbra], &
        &       chket%GetNmaxChannel(naket), [Naket, iaket, jaket, taket, nbket, lbket, jbket], np, j, t)
  end function anti_op_wrap

  function Aop(chbra, bra, chket, ket, n_p, jtot, ttot)
    ! Note: all j ant t are double of the actual values
    integer, intent(in) :: bra(7), ket(7)
    type(NmaxChannelA), intent(in) :: chbra, chket
    type(NonAntisymmetrizedQNsA), pointer :: bra1, ket1
    integer, intent(in) :: n_p, jtot, ttot
    real(8) :: Aop

    integer :: nl, il, jl, tl, nzl, lzl, jzl
    integer :: nr, ir, jr, tr, nzr, lzr, jzr
    integer :: jacl, jacr
    integer :: nl1, il1, jl1, tl1, nl2, ll2, jl2
    integer :: nr1, ir1, jr1, tr1, nr2, lr2, jr2
    integer :: lmin, lmax, l, kmin, kmax, k
    real(8) :: cfpb, cfpk, st, njtot, sjl1, sjl2, sll, gms, ex1, ex2, ex

    nl = bra(1); nr = ket(1)
    il = bra(2); ir = ket(2)
    jl = bra(3); jr = ket(3)
    tl = bra(4); tr = ket(4)
    nzl= bra(5); nzr= ket(5)
    lzl= bra(6); lzr= ket(6)
    jzl= bra(7); jzr= ket(7)
    ex = 0.d0

    ex = 0.d0
    do jacl = 1, chbra%GetNumberNAStates()
      bra1 => chbra%GetNAS(jacl)
      nl1 = bra1%Getna() ! (A-2)-body
      il1 = bra1%Getia() ! (A-2)-body
      jl1 = bra1%Getja() ! (A-2)-body
      tl1 = bra1%Getta() ! (A-2)-body
      nl2 = bra1%Getnb()
      ll2 = bra1%Getlb()
      jl2 = bra1%Getjb()

      do jacr = 1, chket%GetNumberNAStates()
        ket1 => chket%GetNAS(jacr)
        nr1 = ket1%Getna()
        ir1 = ket1%Getia()
        jr1 = ket1%Getja()
        tr1 = ket1%Getta()
        nr2 = ket1%Getnb()
        lr2 = ket1%Getlb()
        jr2 = ket1%Getjb()
        if(nl1 /= nr1) cycle
        if(il1 /= ir1) cycle
        if(jl1 /= jr1) cycle
        if(tl1 /= tr1) cycle
        kmin = max(abs(jl2 - jzr) / 2, abs(jr2 - jzl) / 2, abs(ll2 - lzr), abs(lr2 - lzl))
        kmax = min((jl2 + jzr) / 2, (jr2 + jzl) / 2, (ll2 + lzr), (lr2 + lzl))
        lmin = max(abs(ll2 - lzl), abs(lr2 - lzr))
        lmax = min((ll2 + lzl), (lr2 + lzr))
        if(2*nzl + lzl + 2*nl2 + ll2 /= 2*nzr + lzr + 2*nr2 + lr2) cycle

        ex1 = 0.d0
        do l = lmin, lmax
          ex2 = 0.d0
          do k = kmin, kmax
            njtot = ninej%get(jl1, jl2, jl, jr2, 2*K, jzl, jr, jzr, jtot) * &
                & sqrt(dble(jl+1) * dble(jzl+1) * dble(jr+1) * dble(jzr+1))
            !njtot = ninej%get(jl1, jl, jl2, jr, jtot, jzr, jr2, jzl, 2*K) * &
            !    & sqrt(dble(jl+1) * dble(jzl+1) * dble(jr+1) * dble(jzr+1))
            sjl1 = sixj_2l1j%get(2*ll2, 2*lzr, 2*K, jzr, jl2, 1)
            sjl2 = sixj_2l1j%get(2*lr2, 2*lzl, 2*K, jzl, jr2, 1)
            sll  = sixj_3l%get(2*lzl, 2*lr2, 2*K, 2*lzr, 2*ll2, 2*l)
            ex2 = ex2 + njtot * sjl1 * sjl2 * sll * dble(2*K+1)
          end do
          gms = tmbk%get(nl2, ll2, nzl, lzl, nr2, lr2, nzr, lzr, l)
          ex1 = ex1 + gms * dble(2*l+1) * ex2 * (-1.d0) ** (ll2 + lr2)
        end do
        cfpb = chbra%m(il, jacl)
        cfpk = chket%m(ir, jacr)
        st = sixj_isospin%get(1, tl1, tr, 1, ttot, tl) * sqrt(dble(tr+1) * dble(tl+1))
        ex = ex + (-1.d0) ** ((jl2 + jr2 + tl + tr)/2) * &
            & sqrt(dble(jl2+1) * dble(jr2+1)) * st * cfpb * cfpk * ex1
      end do
    end do
    Aop = - dble(n_p-1) * ex / dble(n_p)
    if(nl == nr .and. il == ir .and. jl == jr .and. tl == tr &
      & .and. nzl == nzr .and. lzl == lzr .and. jzl == jzr) then
      Aop = (1.d0 - dble(n_p-1) * ex) / dble(n_p)
    end if
  end function Aop
end module ABodyJacSpaceIso
