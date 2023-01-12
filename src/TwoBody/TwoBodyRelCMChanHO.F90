module TwoBodyRelCMChanHO
  use TwoBodyRelCMRadQNumbers
  implicit none
  public :: TwoBodyRelCMChanHOBasis
  private

  type :: TwoBodyRelCMChanHOBasis
    integer, private :: j, p, z, j_rel, l_cm, n_states, Nmax, Nmax_rel, Nmax_cm
    integer, private :: l_rel_min, l_rel_max, s_min, s_max
    real(8), private :: hw
    integer, private, allocatable :: nls2idx(:,:,:,:) ! n_rel, n_cm, l_rel, spin
    type(RelCMHarmonicOscillator), allocatable :: HOn(:)
  contains
    procedure :: InitTwoBodyRelCMChanHOBasis
    procedure :: FinTwoBodyRelCMChanHOBasis
    procedure :: GetRelCMHO
    procedure :: GetRelCMHOFromQNumbers
    procedure :: GetNmax
    procedure :: GetJ
    procedure :: GetParity
    procedure :: GetZ
    procedure :: GetJRel
    procedure :: GetLCM
    procedure :: GetNumberStates
    procedure :: GetIndexFromQNs
    procedure :: GetIndexFromHO
    procedure :: GetNmaxRel
    procedure :: GetNmaxCM
    procedure :: GetFrequency
    ! overriding
    generic :: init => InitTwoBodyRelCMChanHOBasis  ! constructor
    generic :: fin => FinTwoBodyRelCMChanHOBasis    ! destructor
    generic :: GetHOs => GetRelCMHO, GetRelCMHOFromQNumbers
    generic :: GetIndex => GetIndexFromQNs, GetIndexFromHO
  end type TwoBodyRelCMChanHOBasis

contains

  function GetJ(this) result(j)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: j
    j = this%j
  end function GetJ

  function GetParity(this) result(p)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: p
    p = this%p
  end function GetParity

  function GetZ(this) result(z)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: z
    z = this%z
  end function GetZ

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJRel(this) result(j_rel)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: j_rel
    j_rel = this%j_rel
  end function GetJRel

  function GetLCM(this) result(l_cm)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: l_cm
    l_cm = this%l_cm
  end function GetLCM

  function GetIndexFromQNs(this,n_rel,l_rel,spin,n_cm) result(idx)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer, intent(in) :: n_rel, l_rel, spin, n_cm
    integer :: idx
    idx = this%nls2idx(n_rel,n_cm,l_rel,spin)
  end function GetIndexFromQNs

  function GetIndexFromHO(this,ho) result(idx)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    type(RelCMHarmonicOscillator), intent(in) :: ho
    integer :: idx
    idx = this%GetIndex(ho%GetNRel(),ho%GetLRel(),ho%GetSpin(),ho%GetNCM())
  end function GetIndexFromHO

  function GetNmaxRel(this) result(Nmax)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax_rel
  end function GetNmaxRel

  function GetNmaxCM(this) result(Nmax)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax_cm
  end function GetNmaxCM

  function GetNumberStates(this) result(n_states)
    class(TwoBodyRelCMChanHOBasis), intent(in) :: this
    integer :: n_states
    n_states = this%n_states
  end function GetNumberStates

  function GetFrequency(this) result(hw)
    class(TwoBodyRelCMChanHOBasis), target, intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetRelCMHO(this, idx) result(r)
    class(TwoBodyRelCMChanHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(RelCMHarmonicOscillator), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%HOn(idx)
  end function GetRelCMHO

  function GetRelCMHOFromQNumbers(this, n_rel, l_rel, spin, n_cm) result(r)
    class(TwoBodyRelCMChanHOBasis), target, intent(in) :: this
    integer, intent(in) :: n_rel, l_rel, spin, n_cm
    type(RelCMHarmonicOscillator), pointer :: r
    r => this%GetHOs( this%GetIndex(n_rel, l_rel, spin, n_cm) )
  end function GetRelCMHOFromQNumbers

  subroutine FinTwoBodyRelCMChanHOBasis(this)
    class(TwoBodyRelCMChanHOBasis), intent(inout) :: this
    deallocate(this%nls2idx)
    deallocate(this%HOn)
  end subroutine FinTwoBodyRelCMChanHOBasis

  subroutine InitTwoBodyRelCMChanHOBasis(this, hw, Nmax, Nmax_rel, j, p, z, j_rel, l_cm, Nmax_cm)
    class(TwoBodyRelCMChanHOBasis), intent(inout) :: this
    real(8), intent(in) :: hw
    integer, intent(in) :: Nmax, Nmax_rel, j, p, z, j_rel, l_cm
    integer, intent(in), optional :: Nmax_cm
    integer :: l_rel, spin
    integer :: s_min, s_max, l_rel_min, l_rel_max
    integer :: cnt, n_rel, n_cm

    this%hw = hw
    this%Nmax = Nmax
    this%Nmax_rel = Nmax_rel
    this%Nmax_cm = Nmax
    if(present(Nmax_cm)) this%Nmax_cm = Nmax_cm
    this%j = j
    this%p = p
    this%z = z
    this%j_rel = j_rel
    this%l_cm = l_cm
    s_min = 9999
    s_max =-9999
    l_rel_min = 9999
    l_rel_max =-9999
    cnt = 0
    do spin = 0, 1
      do l_rel = abs(j_rel-spin), j_rel+spin
        if((-1)**(l_rel+l_cm) /= p) cycle
        if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
        s_min = min(s_min, spin)
        s_max = max(s_max, spin)
        l_rel_min = min(l_rel, l_rel_min)
        l_rel_max = max(l_rel, l_rel_max)
        do n_rel = 0, (Nmax_rel-l_rel)/2
          do n_cm = 0, (this%Nmax_cm-l_cm)/2
            if(2*(n_rel+n_cm)+l_rel+l_cm > this%Nmax) cycle
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%n_states = cnt
    this%l_rel_min = l_rel_min
    this%l_rel_max = l_rel_max
    this%s_min = s_min
    this%s_max = s_max
    allocate(this%HOn(cnt))
    allocate(this%nls2idx(0:Nmax_rel/2, 0:(this%Nmax_cm-l_cm)/2, l_rel_min:l_rel_max, s_min:s_max))
    this%nls2idx(:,:,:,:) = 0

    cnt = 0
    do spin = 0, 1
      do l_rel = abs(j_rel-spin), j_rel+spin
        if((-1)**(l_rel+l_cm) /= p) cycle
        if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
        do n_rel = 0, (Nmax_rel-l_rel)/2
          do n_cm = 0, (this%Nmax_cm-l_cm)/2
            if(2*(n_rel+n_cm)+l_rel+l_cm > this%Nmax) cycle
            cnt = cnt + 1
            call this%HOn(cnt)%set(n_rel,l_rel,spin,n_cm)
            this%nls2idx(n_rel,n_cm,l_rel,spin) = cnt
            !write(*,'(a,8i6)') 'num, |n_rel,l_rel,spin,j_rel,n_cm,l_cm,tz>: ', &
            !    & cnt,n_rel,l_rel,spin,j_rel,n_cm,l_cm,z
          end do
        end do
      end do
    end do
  end subroutine InitTwoBodyRelCMChanHOBasis

end module TwoBodyRelCMChanHO
