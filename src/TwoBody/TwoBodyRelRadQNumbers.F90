module TwoBodyRelRadQNumbers
  implicit none
  public :: HarmonicOscillator
  public :: Mesh
  private :: SetHarmonicOscillator
  private :: SetMesh
  private :: GetHarmonicOscillatorN
  private :: GetHarmonicOscillatorL
  private :: GetHarmonicOscillatorS
  private :: GetMeshI
  private :: GetMeshL
  private :: GetMeshS
  private :: GetMeshP
  private :: GetMeshW

  type :: HarmonicOscillator
    integer, private :: n = -1 ! radial quantum number
    integer, private :: l = -1 ! orbital angular momentum
    integer, private :: s = -1 ! total spin
  contains
    procedure :: SetHarmonicOscillator
    procedure :: GetHarmonicOscillatorN
    procedure :: GetHarmonicOscillatorL
    procedure :: GetHarmonicOscillatorS
    generic :: set => SetHarmonicOscillator
    generic :: GetN => GetHarmonicOscillatorN
    generic :: GetL => GetHarmonicOscillatorL
    generic :: GetS => GetHarmonicOscillatorS
  end type HarmonicOscillator

  type :: Mesh
    integer, private :: i = -1    ! index for mesh
    integer, private :: l = -1    ! orbital angular momentum
    integer, private :: s = -1    ! total spin
    real(8), private :: p = -1.d0 ! mesh point
    real(8), private :: w = -1.d0 ! weight for integral
  contains
    procedure :: SetMesh
    procedure :: GetMeshI
    procedure :: GetMeshL
    procedure :: GetMeshS
    procedure :: GetMeshP
    procedure :: GetMeshW
    generic :: set => SetMesh
    generic :: GetI => GetMeshI
    generic :: GetL => GetMeshL
    generic :: GetS => GetMeshS
    generic :: GetP => GetMeshP
    generic :: GetW => GetMeshW
  end type Mesh

contains

  function GetHarmonicOscillatorN(this) result(N)
    class(HarmonicOscillator), intent(in) :: this
    integer :: N
    N = this%n
  end function GetHarmonicOscillatorN

  function GetHarmonicOscillatorL(this) result(L)
    class(HarmonicOscillator), intent(in) :: this
    integer :: L
    L = this%l
  end function GetHarmonicOscillatorL

  function GetHarmonicOscillatorS(this) result(S)
    class(HarmonicOscillator), intent(in) :: this
    integer :: S
    S = this%s
  end function GetHarmonicOscillatorS

  subroutine SetHarmonicOscillator(this, n, l, s)
    class(HarmonicOscillator), intent(inout) :: this
    integer, intent(in) :: n, l, s
    this%n = n
    this%l = l
    this%s = s
  end subroutine SetHarmonicOscillator

  subroutine SetMesh(this, i, l, s, p, w)
    class(Mesh), intent(inout) :: this
    integer, intent(in) :: i, l, s
    real(8), intent(in) :: p, w
    this%i = i
    this%l = l
    this%s = s
    this%p = p
    this%w = w
  end subroutine SetMesh

  function GetMeshI(this) result(I)
    class(Mesh), intent(in) :: this
    integer :: I
    I = this%i
  end function GetMeshI

  function GetMeshL(this) result(L)
    class(Mesh), intent(in) :: this
    integer :: L
    L = this%l
  end function GetMeshL

  function GetMeshS(this) result(S)
    class(Mesh), intent(in) :: this
    integer :: S
    S = this%s
  end function GetMeshS

  function GetMeshP(this) result(P)
    class(Mesh), intent(in) :: this
    real(8) :: P
    P = this%p
  end function GetMeshP

  function GetMeshW(this) result(W)
    class(Mesh), intent(in) :: this
    real(8) :: W
    W = this%w
  end function GetMeshW
end module TwoBodyRelRadQNumbers
