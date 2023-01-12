module ABodyJacQuantumNumbers
  !
  ! Jacobi basis for A-body system
  !     _______
  !    /       \ na, ia, ja, ta
  !   |        |
  !   |   A    |---------o B
  !   |        |         nb, lb, jb
  !    \______/
  !
  !
  type :: NonAntisymmetrizedQNsA
    integer, private :: na = -100
    integer, private :: ia = -100
    integer, private :: ja = -100
    integer, private :: ta = -100
    integer, private :: nb = -100
    integer, private :: lb = -100
    integer, private :: jb = -100
  contains
    procedure :: SetNonAntisymmetrizedQNsA
    procedure :: Getna
    procedure :: GetIa
    procedure :: Getja
    procedure :: Getta
    procedure :: Getnb
    procedure :: Getlb
    procedure :: Getjb
    generic :: set => SetNonAntisymmetrizedQNsA
  end type NonAntisymmetrizedQNsA
contains

  subroutine SetNonAntisymmetrizedQNsA(this, na, ia, ja, ta, nb, lb, jb)
    class(NonAntisymmetrizedQNsA), intent(inout) :: this
    integer, intent(in) :: na, ia, ja, ta, nb, lb, jb
    this%na = na
    this%ia = ia
    this%ja = ja
    this%ta = ta
    this%nb = nb
    this%lb = lb
    this%jb = jb
  end subroutine SetNonAntisymmetrizedQNsA

  function Getna(this) result(na)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: na
    na = this%na
  end function Getna

  function Getia(this) result(ia)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: ia
    ia = this%ia
  end function Getia

  function Getja(this) result(ja)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: ja
    ja = this%ja
  end function Getja

  function Getta(this) result(ta)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: ta
    ta = this%ta
  end function Getta

  function Getnb(this) result(nb)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: nb
    nb = this%nb
  end function Getnb

  function Getlb(this) result(lb)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: lb
    lb = this%lb
  end function Getlb

  function Getjb(this) result(jb)
    class(NonAntisymmetrizedQNsA), intent(in) :: this
    integer :: jb
    jb = this%jb
  end function Getjb

end module ABodyJacQuantumNumbers
