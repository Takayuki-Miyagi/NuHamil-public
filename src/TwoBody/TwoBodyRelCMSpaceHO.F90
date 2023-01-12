module TwoBodyRelCMSpaceHO
  use TwoBodyRelCMChanHO
  implicit none
  public :: TwoBodyRelCMSpaceHOBasis
  private

  type :: TwoBodyRelCMSpaceHOBasis
    type(TwoBodyRelCMChanHOBasis), allocatable :: jpz_jrel_lcm(:)
    real(8), private :: hw
    integer, private :: Nmax, Nmax_rel, NcmMax, Jmax, NChan, j_rel_max, LcmMax
    integer, private, allocatable :: jpz_jrel_lcm2idx(:,:,:,:,:)
    logical, private :: initialized = .false.
  contains
    procedure :: InitTwoBodyRelCMSpaceHOBasis
    procedure :: FinTwoBodyRelCMSpaceHOBasis
    procedure :: GetNumberChannels
    procedure :: GetIndex
    procedure :: GetJmax
    procedure :: GetJRelmax
    procedure :: GetNmax
    procedure :: GetNmaxRel
    procedure :: GetNmaxCM
    procedure :: GetLcmMax
    procedure :: GetChannelFromIndex
    procedure :: GetChannelFromQNumbers
    procedure :: GetFrequency
    ! overriding
    generic :: init => InitTwoBodyRelCMSpaceHOBasis
    generic :: fin => FinTwoBodyRelCMSpaceHOBasis
    generic :: GetChannel => GetChannelFromIndex, GetChannelFromQNumbers
  end type TwoBodyRelCMSpaceHOBasis

contains

  function GetFrequency(this) result(hw)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    real(8) :: hw
    hw = this%hw
  end function GetFrequency

  function GetNmax(this) result(Nmax)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: Nmax
    Nmax = this%Nmax
  end function GetNmax

  function GetJmax(this) result(Jmax)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: Jmax
    Jmax = this%Jmax
  end function GetJmax

  function GetNmaxRel(this) result(Nmax_rel)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: Nmax_rel
    Nmax_rel = this%Nmax_rel
  end function GetNmaxRel

  function GetNmaxCM(this) result(Nmax_cm)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: Nmax_cm
    Nmax_cm = this%NcmMax
  end function GetNmaxCM

  function GetLcmMax(this) result(LcmMax)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: LcmMax
    LcmMax = this%LcmMax
  end function GetLcmMax

  function GetNumberChannels(this) result(NChan)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: NChan
    NChan = this%NChan
  end function GetNumberChannels

  function GetJRelMax(this) result(j_rel_max)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer :: j_rel_max
    j_rel_max = this%j_rel_max
  end function GetJrelMax

  function GetIndex(this, j, p, z, jrel, lcm) result(idx)
    class(TwoBodyRelCMSpaceHOBasis), intent(in) :: this
    integer, intent(in) :: j, p, z, jrel, lcm
    integer :: idx
    idx = this%jpz_jrel_lcm2idx(j,p,z,jrel,lcm)
  end function GetIndex

  function GetChannelFromIndex(this, idx) result(r)
    class(TwoBodyRelCMSpaceHOBasis), target, intent(in) :: this
    integer, intent(in) :: idx
    type(TwoBodyRelCMChanHOBasis), pointer :: r
    r => null()
    if(idx < 1) return
    r => this%jpz_jrel_lcm(idx)
  end function GetChannelFromIndex

  function GetChannelFromQNumbers(this, j, p, z, jrel, lcm) result(r)
    class(TwoBodyRelCMSpaceHOBasis), target, intent(in) :: this
    integer, intent(in) :: j, p, z, jrel, lcm
    type(TwoBodyRelCMChanHOBasis), pointer :: r
    r => this%GetChannel(this%GetIndex(j,p,z,jrel,lcm))
  end function GetChannelFromQNumbers

  subroutine FinTwoBodyRelCMSpaceHOBasis(this)
    class(TwoBodyRelCMSpaceHOBasis), intent(inout) :: this
    integer :: ich
    do ich = 1, this%NChan
      call this%jpz_jrel_lcm(ich)%fin()
    end do
    deallocate(this%jpz_jrel_lcm)
    deallocate(this%jpz_jrel_lcm2idx)
    this%initialized = .false.
  end subroutine FinTwoBodyRelCMSpaceHOBasis

  subroutine InitTwoBodyRelCMSpaceHOBasis(this, hw, Nmax, Jmax, Nmax_rel, j_rel_max, NcmMax, LcmMax)
    class(TwoBodyRelCMSpaceHOBasis), intent(inout) :: this
    integer, intent(in) :: Nmax, Jmax, Nmax_rel, j_rel_max
    integer, intent(in), optional :: NcmMax, LcmMax
    real(8) :: hw
    integer :: j, p, z, j_rel, l_cm, spin, l_rel, cnt, nch, loop

    this%initialized = .true.
    this%Nmax = Nmax
    this%Jmax = Jmax
    this%Nmax_rel = Nmax_rel
    this%J_rel_max = j_rel_max
    this%NcmMax = Nmax
    this%LcmMax = Jmax + j_rel_max
    if(present(NcmMax)) this%NcmMax = NcmMax
    if(present(LcmMax)) this%LcmMax = LcmMax
    this%hw = hw

    do loop = 1, 2
      nch = 0
      do z = -1, 1
        do j = 0, this%GetJmax()
          do p = 1, -1, -2

            do j_rel = 0, this%GetJRelMax()
              do l_cm = abs(j-j_rel), min(j+j_rel, this%LcmMax)
                if((this%NcmMax-l_cm)/2 < 0) cycle

                cnt = 0
                do spin = 0, 1
                  do l_rel = abs(j_rel-spin), j_rel+spin
                    if((-1)**(l_rel+l_cm) /= p) cycle
                    if(abs(z)==1 .and. (-1)**(l_rel+spin)==-1) cycle
                    if(l_cm + l_rel > this%Nmax) cycle
                    cnt = cnt + 1
                  end do
                end do
                if(cnt > 0) then
                  nch = nch + 1
                  if(loop==2) then
                    this%jpz_jrel_lcm2idx(j,p,z,j_rel,l_cm) = nch
                    call this%jpz_jrel_lcm(nch)%init(hw, Nmax, Nmax_rel, j, p, z, j_rel, l_cm, this%NcmMax)
                  end if
                end if

              end do
            end do

          end do
        end do
      end do
      if(loop==1) then
        this%NChan = nch
        allocate(this%jpz_jrel_lcm(nch))
        allocate(this%jpz_jrel_lcm2idx(0:Jmax,-1:1,-1:1,0:j_rel_max,0:this%LcmMax))
        this%jpz_jrel_lcm2idx(:,:,:,:,:) = 0
      end if
    end do

  end subroutine InitTwoBodyRelCMSpaceHOBasis
end module TwoBodyRelCMSpaceHO

!program test
!  use RelativeSpaceHarmonicOscillator
!  implicit none
!  type(TwoBodyRelSpaceHOBasis) :: ho
!  type(TwoBodyRelativeChannelHOBasis), pointer :: ch
!  type(HarmonicOscillator), pointer :: m
!  integer :: ich, n
!
!  call ho%init(20,8)
!  do ich = 1, ho%NChan
!    ch => ho%getp(ich)
!    write(*,'(3i3)') ch%j, ch%p, ch%z
!  end do
!  call ho%fin()
!
!  call ho%init(20,8,1,1,0)
!  do ich = 1, ho%NChan
!    ch => ho%getp(ich)
!    write(*,'(3i3)') ch%j, ch%p, ch%z
!    do n = 1, ch%n_state
!      m => ch%getp(n)
!      write(*,'(3i3)') m%n, m%l, m%s
!    end do
!  end do
!  call ho%fin()
!end program test
