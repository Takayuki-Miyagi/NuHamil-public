module ThreeBodyJacobiQuantumNumbers
  public :: NonAntisymmetrizedIsoQNs3
  public :: NonAntisymmetrizedQNs3
  public :: AntisymmetrizedQNs3

  private :: SetNonAntisymmetrizedQNs3
  private :: GetN12
  private :: GetL12
  private :: GetS12
  private :: GetJ12
  private :: GetZ12
  private :: GetN3
  private :: GetL3
  private :: GetJ3
  private :: GetZ3

  private :: SetNonAntisymmetrizedIsoQNs3
  private :: GetIsoN12
  private :: GetIsoL12
  private :: GetIsoS12
  private :: GetIsoJ12
  private :: GetIsoT12
  private :: GetIsoN3
  private :: GetIsoL3
  private :: GetIsoJ3

  type :: NonAntisymmetrizedQNs3
    integer, private :: n12 = -100
    integer, private :: l12 = -100
    integer, private :: s12 = -100
    integer, private :: j12 = -100
    integer, private :: z12 = -100
    integer, private :: n3  = -100
    integer, private :: l3  = -100
    integer, private :: j3  = -100
    integer, private :: z3  = -100
  contains
    procedure :: SetNonAntisymmetrizedQNs3
    procedure :: GetN12
    procedure :: GetL12
    procedure :: GetS12
    procedure :: GetJ12
    procedure :: GetZ12
    procedure :: GetN3
    procedure :: GetL3
    procedure :: GetJ3
    procedure :: GetZ3
    generic :: Set => SetNonAntisymmetrizedQNs3
  end type NonAntisymmetrizedQNs3

  type :: NonAntisymmetrizedIsoQNs3
    integer, private :: n12 = -100
    integer, private :: l12 = -100
    integer, private :: s12 = -100
    integer, private :: j12 = -100
    integer, private :: t12 = -100
    integer, private :: n3  = -100
    integer, private :: l3  = -100
    integer, private :: j3  = -100
  contains
    procedure :: SetNonAntisymmetrizedIsoQNs3
    procedure :: GetIsoN12
    procedure :: GetIsoL12
    procedure :: GetIsoS12
    procedure :: GetIsoJ12
    procedure :: GetIsoT12
    procedure :: GetIsoN3
    procedure :: GetIsoL3
    procedure :: GetIsoJ3
    generic :: Set => SetNonAntisymmetrizedIsoQNs3
    generic :: GetN12 => GetIsoN12
    generic :: GetL12 => GetIsoL12
    generic :: GetS12 => GetIsoS12
    generic :: GetJ12 => GetIsoJ12
    generic :: GetT12 => GetIsoT12
    generic :: GetN3 => GetIsoN3
    generic :: GetL3 => GetIsoL3
    generic :: GetJ3 => GetIsoJ3
  end type NonAntisymmetrizedIsoQNs3

contains

  subroutine SetNonAntisymmetrizedIsoQNs3(this,n12,l12,s12,j12,t12,n3,l3,j3)
    class(NonAntisymmetrizedIsoQNs3), intent(inout) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3
    this%n12 = n12
    this%l12 = l12
    this%s12 = s12
    this%j12 = j12
    this%t12 = t12
    this%n3  = n3
    this%l3  = l3
    this%j3  = j3
  end subroutine SetNonAntisymmetrizedIsoQNs3

  subroutine SetNonAntisymmetrizedQNs3(this,n12,l12,s12,j12,z12,n3,l3,j3,z3)
    class(NonAntisymmetrizedQNs3), intent(inout) :: this
    integer, intent(in) :: n12, l12, s12, j12, z12, n3, l3, j3, z3
    this%n12 = n12
    this%l12 = l12
    this%s12 = s12
    this%j12 = j12
    this%z12 = z12
    this%n3  = n3
    this%l3  = l3
    this%j3  = j3
    this%z3  = z3
  end subroutine SetNonAntisymmetrizedQNs3

  function GetN12(this) result(n12)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: n12
    n12 = this%n12
  end function GetN12

  function GetL12(this) result(l12)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: l12
    l12 = this%l12
  end function GetL12

  function GetS12(this) result(s12)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: s12
    s12 = this%s12
  end function GetS12

  function GetJ12(this) result(j12)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: j12
    j12 = this%j12
  end function GetJ12

  function GetZ12(this) result(z12)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: z12
    z12 = this%z12
  end function GetZ12

  function GetN3(this) result(n3)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: n3
    n3 = this%n3
  end function GetN3

  function GetL3(this) result(l3)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: l3
    l3 = this%l3
  end function GetL3

  function GetJ3(this) result(j3)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: j3
    j3 = this%j3
  end function GetJ3

  function GetZ3(this) result(z3)
    class(NonAntisymmetrizedQNs3), intent(in) :: this
    integer :: z3
    z3 = this%z3
  end function GetZ3

  function GetIsoN12(this) result(n12)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: n12
    n12 = this%n12
  end function GetIsoN12

  function GetIsoL12(this) result(l12)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: l12
    l12 = this%l12
  end function GetIsoL12

  function GetIsoS12(this) result(s12)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: s12
    s12 = this%s12
  end function GetIsoS12

  function GetIsoJ12(this) result(j12)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: j12
    j12 = this%j12
  end function GetIsoJ12

  function GetIsoT12(this) result(t12)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: t12
    t12 = this%t12
  end function GetIsoT12

  function GetIsoN3(this) result(n3)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: n3
    n3 = this%n3
  end function GetIsoN3

  function GetIsoL3(this) result(l3)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: l3
    l3 = this%l3
  end function GetIsoL3

  function GetIsoJ3(this) result(j3)
    class(NonAntisymmetrizedIsoQNs3), intent(in) :: this
    integer :: j3
    j3 = this%j3
  end function GetIsoJ3

end module ThreeBodyJacobiQuantumNumbers
