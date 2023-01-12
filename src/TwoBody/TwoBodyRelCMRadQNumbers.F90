module TwoBodyRelCMRadQNumbers
  implicit none
  public :: RelCMHarmonicOscillator
  public :: RelCMMesh
  private

  type :: RelCMHarmonicOscillator
    integer, private :: n_rel = -1 ! radial quantum number (rel)
    integer, private :: l_rel = -1 ! orbital angular momentum (rel)
    integer, private :: spin = -1  ! total spin
    integer, private :: n_cm = -1  ! radial quantum number (cm)
  contains
    procedure :: SetRelCMHarmonicOscillator
    procedure :: GetRelCMHarmonicOscillatorNRel
    procedure :: GetRelCMHarmonicOscillatorLRel
    procedure :: GetRelCMHarmonicOscillatorSpin
    procedure :: GetRelCMHarmonicOscillatorNCM
    generic :: set => SetRelCMHarmonicOscillator
    generic :: GetNRel => GetRelCMHarmonicOscillatorNRel
    generic :: GetLRel => GetRelCMHarmonicOscillatorLRel
    generic :: GetSpin => GetRelCMHarmonicOscillatorSpin
    generic :: GetNCM  => GetRelCMHarmonicOscillatorNCM
  end type RelCMHarmonicOscillator

  type :: RelCMMesh
    integer, private :: i_rel = -1    ! index for mesh
    integer, private :: l_rel = -1    ! orbital angular momentum
    integer, private :: spin = -1     ! total spin
    integer, private :: i_cm = -1     ! index for mesh
    real(8), private :: p_rel = -1.d0 ! mesh point
    real(8), private :: w_rel = -1.d0 ! weight for integral
    real(8), private :: p_cm = -1.d0  ! mesh point
    real(8), private :: w_cm = -1.d0  ! weight for integral
  contains
    procedure :: SetRelCMMesh
    procedure :: GetRelCMMeshIRel
    procedure :: GetRelCMMeshLRel
    procedure :: GetRelCMMeshSpin
    procedure :: GetRelCMMeshICM
    procedure :: GetRelCMMeshPRel
    procedure :: GetRelCMMeshPCM
    procedure :: GetRelCMMeshWRel
    procedure :: GetRelCMMeshWCM
    generic :: set => SetRelCMMesh
    generic :: GetIRel => GetRelCMMeshIRel
    generic :: GetLRel => GetRelCMMeshLRel
    generic :: GetSpin => GetRelCMMeshSpin
    generic :: GetICM  => GetRelCMMeshICM
    generic :: GetPRel => GetRelCMMeshPRel
    generic :: GetPCM  => GetRelCMMeshPCM
    generic :: GetWRel => GetRelCMMeshWRel
    generic :: GetWCM  => GetRelCMMeshWCM
  end type RelCMMesh

contains

  subroutine SetRelCMHarmonicOscillator(this, n_rel, l_rel, spin, n_cm)
    class(RelCMHarmonicOscillator), intent(inout) :: this
    integer, intent(in) :: n_rel, l_rel, spin, n_cm
    this%n_rel = n_rel
    this%l_rel = l_rel
    this%spin = spin
    this%n_cm = n_cm
  end subroutine SetRelCMHarmonicOscillator

  function GetRelCMHarmonicOscillatorNRel(this) result(n_rel)
    class(RelCMHarmonicOscillator), intent(in) :: this
    integer :: n_rel
    n_rel = this%n_rel
  end function GetRelCMHarmonicOscillatorNRel

  function GetRelCMHarmonicOscillatorLRel(this) result(l_rel)
    class(RelCMHarmonicOscillator), intent(in) :: this
    integer :: l_rel
    l_rel = this%l_rel
  end function GetRelCMHarmonicOscillatorLRel

  function GetRelCMHarmonicOscillatorSpin(this) result(spin)
    class(RelCMHarmonicOscillator), intent(in) :: this
    integer :: spin
    spin = this%spin
  end function GetRelCMHarmonicOscillatorSpin

  function GetRelCMHarmonicOscillatorNCM(this) result(n_cm)
    class(RelCMHarmonicOscillator), intent(in) :: this
    integer :: n_cm
    n_cm = this%n_cm
  end function GetRelCMHarmonicOscillatorNCM

  subroutine SetRelCMMesh(this, i_rel, l_rel, spin, i_cm, p_rel, w_rel, p_cm, w_cm)
    class(RelCMMesh), intent(inout) :: this
    integer, intent(in) :: i_rel, l_rel, spin, i_cm
    real(8), intent(in) :: p_rel, w_rel, p_cm, w_cm
    this%i_rel = i_rel
    this%l_rel = l_rel
    this%spin = spin
    this%i_cm = i_cm
    this%p_rel = p_rel
    this%w_rel = w_rel
    this%p_cm = p_cm
    this%w_cm = w_cm
  end subroutine SetRelCMMesh

  function GetRelCMMeshIRel(this) result(i_rel)
    class(RelCMMesh), intent(in) :: this
    integer :: i_rel
    i_rel = this%i_rel
  end function GetRelCMMeshIRel

  function GetRelCMMeshLRel(this) result(l_rel)
    class(RelCMMesh), intent(in) :: this
    integer :: l_rel
    l_rel = this%l_rel
  end function GetRelCMMeshLRel

  function GetRelCMMeshSpin(this) result(spin)
    class(RelCMMesh), intent(in) :: this
    integer :: spin
    spin = this%spin
  end function GetRelCMMeshSpin

  function GetRelCMMeshICM(this) result(i_cm)
    class(RelCMMesh), intent(in) :: this
    integer :: i_cm
    i_cm = this%i_cm
  end function GetRelCMMeshICM

  function GetRelCMMeshPRel(this) result(p_rel)
    class(RelCMMesh), intent(in) :: this
    real(8) :: p_rel
    p_rel = this%p_rel
  end function GetRelCMMeshPRel

  function GetRelCMMeshWRel(this) result(w_rel)
    class(RelCMMesh), intent(in) :: this
    real(8) :: w_rel
    w_rel = this%w_rel
  end function GetRelCMMeshWRel

  function GetRelCMMeshPCM(this) result(p_cm)
    class(RelCMMesh), intent(in) :: this
    real(8) :: p_cm
    p_cm = this%p_cm
  end function GetRelCMMeshPCM

  function GetRelCMMeshWCM(this) result(w_cm)
    class(RelCMMesh), intent(in) :: this
    real(8) :: w_cm
    w_cm = this%w_cm
  end function GetRelCMMeshWCM

end module TwoBodyRelCMRadQNumbers
