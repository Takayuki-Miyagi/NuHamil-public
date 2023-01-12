module OperatorDefinitions
  use NuHamilInput, only: InputParameters
  use ClassSys
  implicit none
  public :: OperatorDef
  public :: CalcMERel
  public :: CalcMERelCM
  public :: CalcMEOneBody
  public :: SetParameters
  public :: SetFrequency
  public :: PrintParameters
  public :: parse_operator_string
  public :: parse_opname
  public :: set_pv_couplings
  public :: set_pvtv_couplings
  private

  type :: OperatorDef
    type(str), private :: OpName
    integer, private :: jr=-1
    integer, private :: pr=-100
    integer, private :: tr=-1
    integer, private :: zr=-100
    real(8), private :: Q=0.d0
    logical, private :: pn_formalism = .false.
    logical, private :: reduced_matrix_element = .true.
    logical, private :: subtract_onebody=.false.
    logical, private :: channel_restriction=.true.
    real(8), private :: op_args(6)
  contains
    procedure :: InitOperatorDefFromJPT
    procedure :: InitOperatorDefFromName
    procedure :: FinOperatorDef
    procedure :: CopyOperatorDef
    procedure :: GetOpJ
    procedure :: GetOpP
    procedure :: GetOpT
    procedure :: GetOpZ
    procedure :: GetQ
    procedure :: GetSubtract
    procedure :: SetOpJ
    procedure :: SetOpP
    procedure :: SetOpT
    procedure :: SetOpZ
    procedure :: SetQ
    procedure :: SetOpName
    procedure :: SetSubtract
    procedure :: pn
    procedure :: reduced_me
    procedure :: SetReduced
    procedure :: GetOpName
    procedure :: restricted
    generic :: assignment(=) => CopyOperatorDef
    generic :: InitOpDef => InitOperatorDefFromJPT, InitOperatorDefFromName
    generic :: FinOpDef => FinOperatorDef
  end type OperatorDef
  type(InputParameters), private :: parameters
contains

  subroutine SetParameters(params)
    type(InputParameters), intent(in) :: params
    parameters = params
  end subroutine SetParameters

  subroutine SetFrequency(hw)
    real(8), intent(in) :: hw
    parameters%hw = hw
  end subroutine SetFrequency

  subroutine PrintParameters()
    call parameters%PrintParametersALL()
  end subroutine PrintParameters

  subroutine InitOperatorDefFromJPT(this, pn_form, jr, pr, t_z_r)
    class(OperatorDef), intent(inout) :: this
    logical, intent(in) :: pn_form
    integer, intent(in) :: jr, pr, t_z_r

    this%pn_formalism = pn_form
    this%jr = jr
    this%pr = pr
    if(      this%pn_formalism) this%zr = t_z_r
    if(.not. this%pn_formalism) this%tr = t_z_r
    if(      this%pn_formalism .and. this%jr==0 .and. this%pr==1 .and. this%zr==0) this%reduced_matrix_element=.false.
    if(.not. this%pn_formalism .and. this%jr==0 .and. this%pr==1 .and. this%tr==0) this%reduced_matrix_element=.false.
    if(.not. allocated(this%OpName%val)) this%OpName = "Unknown"
  end subroutine InitOperatorDefFromJPT

  subroutine InitOperatorDefFromName(this, opname, pn_form)
    class(OperatorDef), intent(inout) :: this
    type(str), intent(in) :: opname
    logical, intent(in) :: pn_form
    this%opname = opname
    this%pn_formalism = pn_form
    if(this%pn_formalism) call get_operator_rank_pn(this)
    if(.not. this%pn_formalism) call get_operator_rank_iso(this)
  end subroutine InitOperatorDefFromName

  subroutine FinOperatorDef(this)
    class(OperatorDef), intent(inout) :: this
    deallocate(this%opname%val)
    this%pn_formalism = .false.
    this%reduced_matrix_element = .true.
    this%jr = -1
    this%pr = -100
    this%zr = -100
    this%tr = -1
  end subroutine FinOperatorDef

  subroutine CopyOperatorDef(a, b)
    class(OperatorDef), intent(inout) :: a
    type(OperatorDef), intent(in) :: b
    a%opname = b%opname
    a%pn_formalism= b%pn_formalism
    a%reduced_matrix_element = b%reduced_matrix_element
    a%jr= b%jr
    a%pr= b%pr
    a%zr= b%zr
    a%tr= b%tr
  end subroutine CopyOperatorDef

  function GetOpName(this) result(opname)
    class(OperatorDef), intent(in) :: this
    type(str) :: opname
    opname = this%opname
  end function GetOpName

  subroutine SetOpName(this, opname)
    class(OperatorDef), intent(inout) :: this
    type(str), intent(in) :: opname
    this%OpName = OpName
  end subroutine SetOpName

  function GetOpJ(this) result(J)
    class(OperatorDef), intent(in) :: this
    integer :: J
    J = this%jr
  end function GetOpJ

  function GetOpP(this) result(P)
    class(OperatorDef), intent(in) :: this
    integer :: P
    P = this%Pr
  end function GetOpP

  function GetOpT(this) result(T)
    class(OperatorDef), intent(in) :: this
    integer :: T
    T = this%Tr
  end function GetOpT

  function GetOpZ(this) result(Z)
    class(OperatorDef), intent(in) :: this
    integer :: Z
    Z = this%Zr
  end function GetOpZ

  function GetQ(this) result(Q)
    class(OperatorDef), intent(in) :: this
    real(8) :: Q
    Q = this%Q
  end function GetQ

  function GetSubtract(this) result(r)
    class(OperatorDef), intent(in) :: this
    logical :: r
    r = this%subtract_onebody
  end function GetSubtract

  subroutine SetOpJ(this,J)
    class(OperatorDef), intent(inout) :: this
    integer, intent(in) :: J
    this%jr = J
  end subroutine SetOpJ

  subroutine SetOpP(this,P)
    class(OperatorDef), intent(inout) :: this
    integer, intent(in) :: P
    this%Pr = P
  end subroutine SetOpP

  subroutine SetOpT(this,T)
    class(OperatorDef), intent(inout) :: this
    integer, intent(in) :: T
    this%Tr = T
  end subroutine SetOpT

  subroutine SetOpZ(this,Z)
    class(OperatorDef), intent(inout) :: this
    integer, intent(in) :: Z
    this%Zr = Z
  end subroutine SetOpZ

  subroutine SetQ(this, Q)
    class(OperatorDef), intent(inout) :: this
    real(8), intent(in) :: Q
    this%Q = Q
  end subroutine SetQ

  subroutine SetSubtract(this,subtract)
    class(OperatorDef), intent(inout) :: this
    logical, intent(in) :: subtract
    this%subtract_onebody = subtract
  end subroutine SetSubtract

  function pn(this)
    class(OperatorDef), intent(in) :: this
    logical :: pn
    pn = this%pn_formalism
  end function pn

  function reduced_me(this)
    class(OperatorDef), intent(in) :: this
    logical :: reduced_me
    reduced_me = this%reduced_matrix_element
  end function reduced_me

  subroutine SetReduced(this, reduced)
    class(OperatorDef), intent(inout) :: this
    logical, intent(in) :: reduced
    this%reduced_matrix_element = reduced
  end subroutine SetReduced

  function restricted(this) result(r)
    class(OperatorDef), intent(in) :: this
    logical :: r
    r = this%channel_restriction
  end function restricted

  subroutine get_operator_rank_pn(this)
    use ClassSys, only: sys
    class(OperatorDef), intent(inout) :: this
    type(str) :: OpName
    type(str), allocatable :: strings(:)
    type(sys) :: s
    integer :: J_read
    real(8) :: Q_read

    select case(this%opname%val)

    case('Hamil','hamil', "NNint", "NNNint", &
          &'UnitaryTransformation','UT',"NNNfrom2N","NNNfrom2N_central",&
          &"NNNfrom2N_spinorbit","NNNfrom2N_tensor","NNNfromTkin", "NNNinduced", &
          &"NNN_c1","NNN_c3","NNN_c4","NNN_TPE","NNN_cD","NNN_OPE","NNN_cE",&
          &"NNN_Contact","NNN_Genuine", "NNNinduced_N3LO_EM500_OPE", &
          &"NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
      this%jr = 0
      this%pr = 1
      this%zr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction=.true.
    case('kinetic','Kinetic','hopot','HOpot', "Hcm", "HOHamil", &
          &"T_l0", "T_radial", "T_Gj", "T_Gl", 'R2','r2')
      this%jr = 0
      this%pr = 1
      this%zr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction=.false.
    case('Sigma','Sigma_Tauz','Spin')
      this%jr = 1
      this%pr = 1
      this%zr = 0
      this%subtract_onebody=.true.
      this%channel_restriction=.false.
    case('M1_L_IS','M1_L_IV','M1_S_IS','M1_S_IV')
      this%jr = 1
      this%pr = 1
      this%zr = 0
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case('M1_2BC', 'M1_2BC_intr', 'M1_2BC_Sachs')
      this%jr = 1
      this%pr = 1
      this%zr = 0
      this%channel_restriction=.true.
    case('E2_IS','E2_IV')
      this%jr = 2
      this%pr = 1
      this%zr = 0
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case('E1',"TDM_mag", "TDM_conv","E1cm")
      this%jr = 1
      this%pr = -1
      this%zr = 0
      this%channel_restriction=.false.
    case("GamowTeller")
      this%jr = 1
      this%pr = 1
      this%zr = 1
      !this%subtract_onebody=.true.
      this%channel_restriction=.false.
    case("Fermi")
      this%jr = 0
      this%pr = 1
      this%zr = 1
      this%subtract_onebody=.true.
      this%channel_restriction=.false.
    case("DFermi","DGamowTeller0")
      this%jr = 0
      this%pr = 1
      this%zr = 2
      this%channel_restriction=.false.
    case("DGamowTeller2")
      this%jr = 2
      this%pr = 1
      this%zr = 2
      this%channel_restriction=.false.
    case default
      if( this%OpName%val(:4) == "Sp_E") then
        OpName = this%OpName
        call s%split(OpName, s%str("Sp_E"), strings)
        read(strings(2)%val,*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read)
        this%zr = 0
        return
      end if

      if( s%find(this%OpName,s%str("Sp_M"))) then
        OpName = this%OpName
        call s%split(OpName, s%str("Sp_M"), strings)
        read(strings(2)%val,*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        this%zr = 0
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('M5_'))) then
        ! format example:
        !   M5_1B_J0_Tz0_Q10
        !   M5_2B_J0_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('M_'))) then
        ! format example:
        !   M_1B_J0_Tz0_Q10
        !   M_2B_J0_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if


      if( s%find(this%OpName,s%str('L5_'))) then
        ! format example:
        !   L5_1B_J1_Tz0_Q10
        !   L5_2B_J1_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('L_'))) then
        ! format example:
        !   L_1B_J1_Tz0_Q10
        !   L_2B_J1_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tel_'))) then
        ! format example:
        !   Tel_1B_J1_Tz0_Q10
        !   Tel_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tmag_'))) then
        ! format example:
        !   Tmag_1B_J1_Tz0_Q10
        !   Tmag_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tel5_'))) then
        ! format example:
        !   Tel5_1B_J1_Tz0_Q10
        !   Tel5_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tmag5_'))) then
        ! format example:
        !   Tmag5_1B_J1_Tz0_Q10
        !   Tmag5_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read)
        read(strings(4)%val(3:),*) J_read
        this%zr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str("0vbb"))) then
        !
        ! Format "0vbb(op type)_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "0vbbFermi_Ec_7.72-N2LO-NonLocal2-500", "0vbbContact_Ec_0-LO-NonLocal2-500"
        !
        this%jr = 0
        this%pr = 1
        this%zr = 2
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("AxialV_Tz0"))) then
        !
        ! Format "AxialV_Tz0_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "AxialV_Tz0-N2LO", "AxialV_Z0-N2LO-NonLocal2-500"
        !
        this%jr = 1
        this%pr = 1
        this%zr = 0
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("AxialV_Tz1"))) then
        !
        ! Format "AxialV_Tz1_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "AxialV_Tz1-N2LO", "AxialV_Tz1-N2LO-NonLocal2-500"
        !
        this%jr = 1
        this%pr = 1
        this%zr = 1
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("Vector_Tz0"))) then
        !
        ! Format "Vector_Tz0_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "Vector_Tz0-NLO", "Vector_Z0-NLO-NonLocal2-500"
        !
        this%jr = 1
        this%pr =-1
        this%zr = 0
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("Vector_Tz1"))) then
        !
        ! Format "Vector_Tz1_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "Vector_Tz1-NLO", "Vector_Tz1-NLO-NonLocal2-500"
        !
        this%jr = 1
        this%pr =-1
        this%zr = 1
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("PVTCint"))) then
        !
        ! Format "PVTCint_(op parameter name1)_(op parameter value1)-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "PVTCint_ChEFT_hpi1-N2LO-NonLocal2-500", "PVTCint_DDH_best-LO-NonLocal2-500"
        !
        this%jr = 0
        this%pr = -1
        this%zr = 0
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("PVTVint"))) then
        !
        ! Format "PVTVint_Term_(Term)-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "PVTVint_Term_hpi1-N2LO-NonLocal2-500"
        !
        this%jr = 0
        this%pr = -1
        this%zr = 0
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("HOHamil"))) then
        ! Format "HOHamil(frequency)"
        this%jr = 0
        this%pr = 1
        this%zr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction=.false.
        return
      end if

      if( s%find(this%OpName,s%str("T_Gl"))) then
        ! Format "T_Gl(width)": Tkin * exp( -(l/width)**2 )
        this%jr = 0
        this%pr = 1
        this%zr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction=.false.
        call s%split( this%OpName, s%str("_Gl"), strings)
        this%op_args(1) = 1.d0
        if(size(strings)>1) read(strings(2)%val, *) this%op_args(1)
        return
      end if

      if( s%find(this%OpName,s%str("T_Gj"))) then
        ! Format "T_Gl(width)": Tkin * exp( -(j/width)**2 )
        this%jr = 0
        this%pr = 1
        this%zr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction=.false.
        call s%split( this%OpName, s%str("_Gj"), strings)
        this%op_args(1) = 1.d0
        if(size(strings)>1) read(strings(2)%val, *) this%op_args(1)
        return
      end if

      if( s%find(this%OpName,s%str("Gaussian"))) then
        this%jr = 0
        this%pr = 1
        this%zr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction=.false.
        call s%split( this%OpName, s%str("_"), strings)
        this%op_args(1) = 1.d0
        if(size(strings)>1) read(strings(2)%val, *) this%op_args(1)
        return
      end if

      write(*,'(2a)') 'In get_operator_rank_pn, selected operator has not been implemented: oprtr=', trim(this%opname%val)
      write(*,"(a)") "Assuming scalar J=0, P=1, Z=0"
      this%jr = 0
      this%pr = 1
      this%zr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction=.true.
      stop
    end select
  end subroutine get_operator_rank_pn

  subroutine get_operator_rank_iso(this)
    use ClassSys, only: sys, str
    class(OperatorDef), intent(inout) :: this
    type(sys) :: s
    type(str), allocatable :: strings(:)
    type(str) :: OpName
    integer :: J_read
    real(8) :: Q_read

    select case(this%opname%val)

    case('Hamil','hamil', "NNint", "NNNint", &
          &'UnitaryTransformation','UT',"NNNfrom2N","NNNfrom2N_central",&
          &"NNNfrom2N_spinorbit","NNNfrom2N_tensor","NNNfromTkin","NNNinduced",&
          &"NNN_c1","NNN_c3","NNN_c4","NNN_TPE","NNN_cD","NNN_OPE","NNN_cE",&
          &"NNN_Contact","NNN_Genuine", "NNNinduced_N3LO_EM500_OPE", &
          &"NNNinduced_N3LO_EM500_TPE", "NNNinduced_N3LO_EM500_Contacts")
      this%jr = 0
      this%pr = 1
      this%tr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction = .true.
    case('kinetic','Kinetic','hopot','HOpot','R2','r2',&
          &"HOHamil", "T_l0", "T_radial","T_Gj", "T_Gl")
      this%jr = 0
      this%pr = 1
      this%tr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction = .false.
    case('R2_IV')
      this%jr = 0
      this%pr = 1
      this%tr = 1
      this%channel_restriction = .false.
    case('Sigma','Spin')
      this%jr = 1
      this%pr = 1
      this%tr = 0
      this%channel_restriction = .false.
    case('Sigma_Tauz')
      this%jr = 1
      this%pr = 1
      this%tr = 1
      this%channel_restriction = .false.
    case('E1', "TDM_mag_st")
      this%jr = 1
      this%pr = -1
      this%tr = 1
      this%channel_restriction = .false.
    case('M1_2BC', 'M1_2BC_intr', 'M1_2BC_Sachs')
      this%jr = 1
      this%pr = 1
      this%tr = 1
      this%channel_restriction = .true.
    case('M1_L_IS','M1_S_IS')
      this%jr = 1
      this%pr = 1
      this%tr = 0
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case('M1_L_IV','M1_S_IV')
      this%jr = 1
      this%pr = 1
      this%tr = 1
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case('E2_IS')
      this%jr = 2
      this%pr = 1
      this%tr = 0
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case('E2_IV')
      this%jr = 2
      this%pr = 1
      this%tr = 1
      this%subtract_onebody=.false.
      this%channel_restriction=.false.
    case("TDM_mag_s")
      this%jr = 1
      this%pr = -1
      this%tr = 0
      this%channel_restriction = .false.
    case("Fermi")
      this%jr = 0
      this%pr = 1
      this%tr = 1
      this%channel_restriction = .false.
    case("GamowTeller")
      this%jr = 1
      this%pr = 1
      this%tr = 1
      this%channel_restriction = .false.
    case("DFermi","DGamowTeller0")
      this%jr = 0
      this%pr = 1
      this%tr = 2
      this%channel_restriction = .false.
    case("DGamowTeller2")
      this%jr = 2
      this%pr = 1
      this%tr = 2
      this%channel_restriction = .false.
    case default
      if( s%find(this%OpName,s%str("0vbb"))) then
        this%jr = 0
        this%pr = 1
        this%tr = 2
        this%reduced_matrix_element = .true.
        this%channel_restriction = .true.
        return
      end if

      if( s%find(this%OpName,s%str("AxialV_T1"))) then
        !
        ! Format "AxialV_T1_(intrinsic pararameter name1)_(intrinsic parameter value1)_(name2)_(value2)...-(order)-(regulator)(regulator_power)-(cutoff)"
        ! Ex. "AxialV_T1-N2LO", "AxialV_T1-N2LO-NonLocal2-500"
        !
        this%jr = 1
        this%pr = 1
        this%tr = 1
        this%reduced_matrix_element = .true.
        this%channel_restriction=.true.
        return
      end if

      if( s%find(this%OpName,s%str("PVTCint"))) then
        this%jr = 0
        this%pr = -1
        if(s%find(this%OpName, s%str("hpi1"))) this%tr = 1
        if(s%find(this%OpName, s%str("Ct1")))  this%tr = 0
        if(s%find(this%OpName, s%str("Ct2")))  this%tr = 0
        if(s%find(this%OpName, s%str("Ct3")))  this%tr = 1
        if(s%find(this%OpName, s%str("Ct4")))  this%tr = 1
        if(s%find(this%OpName, s%str("Ct5")))  this%tr = 2
        if(s%find(this%OpName, s%str("hV0")))  this%tr = 0
        if(s%find(this%OpName, s%str("hV1")))  this%tr = 1
        if(s%find(this%OpName, s%str("hV2")))  this%tr = 2
        if(s%find(this%OpName, s%str("hA1")))  this%tr = 1
        if(s%find(this%OpName, s%str("hA2")))  this%tr = 2
        if(s%find(this%OpName, s%str("rho0")))  this%tr = 0
        if(s%find(this%OpName, s%str("rho1")))  this%tr = 1
        if(s%find(this%OpName, s%str("rho1p")))  this%tr = 1
        if(s%find(this%OpName, s%str("rho2")))  this%tr = 2
        if(s%find(this%OpName, s%str("omega0")))  this%tr = 0
        if(s%find(this%OpName, s%str("omega1")))  this%tr = 1
        if(s%find(this%OpName, s%str("beta"))) then
          write(*,*) "Isospin cannot be assigned."
          stop
        end if
        this%reduced_matrix_element = .true.
        this%channel_restriction = .true.
        return
      end if

      if( s%find(this%OpName,s%str("PVTVint"))) then
        this%jr = 0
        this%pr = -1
        if(s%find(this%OpName, s%str("gpi0"))) this%tr = 0
        if(s%find(this%OpName, s%str("gpi1"))) this%tr = 1
        if(s%find(this%OpName, s%str("gpi2"))) this%tr = 2
        if(s%find(this%OpName, s%str("delta"))) this%tr = 1
        if(s%find(this%OpName, s%str("Ct1")))  this%tr = 0
        if(s%find(this%OpName, s%str("Ct2")))  this%tr = 0
        if(s%find(this%OpName, s%str("Ct3")))  this%tr = 1
        if(s%find(this%OpName, s%str("Ct4")))  this%tr = 1
        if(s%find(this%OpName, s%str("Ct5")))  this%tr = 2
        if(s%find(this%OpName, s%str("grho0"))) this%tr = 0
        if(s%find(this%OpName, s%str("grho1"))) this%tr = 1
        if(s%find(this%OpName, s%str("grho2"))) this%tr = 2
        if(s%find(this%OpName, s%str("geta0"))) this%tr = 0
        if(s%find(this%OpName, s%str("geta1"))) this%tr = 1
        if(s%find(this%OpName, s%str("gomega0"))) this%tr = 0
        if(s%find(this%OpName, s%str("gomega1"))) this%tr = 1
        this%reduced_matrix_element = .true.
        this%channel_restriction = .true.
        return
      end if

      if( s%find(this%OpName,s%str("HOHamil"))) then
        ! Format "HOHamil(frequency)"
        this%jr = 0
        this%pr = 1
        this%tr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction = .true.
        return
      end if

      if( s%find(this%OpName,s%str("T_Gl"))) then
        ! Format "T_Gl(width)": Tkin * exp( -(l/width)**2 )
        this%jr = 0
        this%pr = 1
        this%tr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction = .true.
        call s%split( this%OpName, s%str("_Gl"), strings)
        this%op_args(1) = 1.d0
        if(size(strings)>1) read(strings(2)%val, *) this%op_args(1)
        return
      end if

      if( s%find(this%OpName,s%str("T_Gj"))) then
        ! Format "T_Gl(width)": Tkin * exp( -(j/width)**2 )
        this%jr = 0
        this%pr = 1
        this%tr = 0
        this%reduced_matrix_element = .false.
        this%channel_restriction = .true.
        call s%split( this%OpName, s%str("_Gj"), strings)
        this%op_args(1) = 1.d0
        if(size(strings)>1) read(strings(2)%val, *) this%op_args(1)
        return
      end if

      if( s%find(this%OpName,s%str('M5_'))) then
        ! format example:
        !   M5_1B_J0_Tz0_Q10
        !   M5_2B_J0_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('M_'))) then
        ! format example:
        !   M_1B_J0_Tz0_Q10
        !   M_2B_J0_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('L5_'))) then
        ! format example:
        !   L5_1B_J1_Tz0_Q10
        !   L5_2B_J1_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('L_'))) then
        ! format example:
        !   L_1B_J1_Tz0_Q10
        !   L_2B_J1_Tz0_Q10

        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tel_'))) then
        ! format example:
        !   Tel_1B_J1_Tz0_Q10
        !   Tel_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**J_read
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tmag_'))) then
        ! format example:
        !   Tmag_1B_J1_Tz0_Q10
        !   Tmag_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tel5_'))) then
        ! format example:
        !   Tel5_1B_J1_Tz0_Q10
        !   Tel5_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read+1)
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      if( s%find(this%OpName,s%str('Tmag5_'))) then
        ! format example:
        !   Tmag5_1B_J1_Tz0_Q10
        !   Tmag5_2B_J1_Tz0_Q10
        OpName = this%OpName
        call s%split(OpName, s%str("_"), strings)
        read(strings(3)%val(2:),*) J_read
        this%jr = J_read
        this%pr = (-1)**(J_read)
        read(strings(4)%val(3:),*) J_read
        this%tr = J_read
        read(strings(5)%val(2:),*) Q_read
        this%Q = Q_read
        this%reduced_matrix_element = .true.
        return
      end if

      write(*,'(2a)') 'In get_operator_rank_iso, selected operator has not been implemented: oprtr=', trim(this%opname%val)
      write(*,"(a)") "Assuming scalar J=0, P=1, T=0"
      this%jr = 0
      this%pr = 1
      this%tr = 0
      this%reduced_matrix_element = .false.
      this%channel_restriction = .true.
      stop
    end select
  end subroutine get_operator_rank_iso

  function CalcMEOneBody(this, bra, ket, isospin, e_charge) result(r)
    use MyLibrary
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(4), ket(4)
    integer, intent(in), optional :: isospin
    real(8), intent(in), optional :: e_charge
    integer :: nbra, lbra, jbra, zbra
    integer :: nket, lket, jket, zket
    real(8) :: r, b2, e
    type(str) :: OpName
    type(sys) :: s

    OpName = this%OpName
    nbra = bra(1)
    lbra = bra(2)
    jbra = bra(3)
    zbra = bra(4)

    nket = ket(1)
    lket = ket(2)
    jket = ket(3)
    zket = ket(4)
    r = 0.d0
    b2 = hc**2/( (m_proton+m_neutron)*0.5*parameters%hw)
    e = 1.d0
    if(present(e_charge)) e = e_charge
    if(present(isospin)) then
      if(isospin==0) e = 0.5d0 * e
      if(isospin==1 .and. zket == -1) e = 0.5d0 * e
      if(isospin==1 .and. zket ==  1) e =-0.5d0 * e
    end if

    select case(OpName%val)
    case( "HOHamil" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_kinetic(nbra,lbra,nket,lket,parameters%hw) + &
          & two_body_hopot(nbra,lbra,nket,lket,parameters%hw)
      return
    case( "kinetic", "Kinetic" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_kinetic(nbra,lbra,nket,lket,parameters%hw)
      return
    case( "T_radial" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_kinetic(nbra,lbra,nket,lket,parameters%hw) - &
          & kin_rotation_part(nbra,lbra,nket,lket,parameters%hw)
      return
    case( "T_Gl" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_kinetic(nbra,lbra,nket,lket,parameters%hw) * exp( - (dble(lket)/this%op_args(1))**2 )
      return
    case( "T_Gj" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_kinetic(nbra,lbra,nket,lket,parameters%hw) * exp( - (dble(jket)/this%op_args(1))**2 )
      return
    case( "hopot", "HOpot" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_hopot(nbra,lbra,nket,lket,parameters%hw)
      return
    case( "R2", "r2" )
      if(zbra /= zket) return
      if(jbra /= jket) return
      r = two_body_hopot(nbra,lbra,nket,lket,parameters%hw)
      r = r * 0.5d0 * b2 / parameters%hw
      return
    case( "Sigma", 'Spin' )
      if(zbra /= zket) return
      if(nbra /= nket) return
      if(lbra /= lket) return
      r = j_to_ls_1( lbra, jbra, lket, jket, 0, 1, 1) * &
          & sqrt(dble(2*lket+1)) * sqrt(6.d0)
      return
    case( "Sigma_Tauz" )
      if(zbra /= zket) return
      if(nbra /= nket) return
      if(lbra /= lket) return
      r = j_to_ls_1( lbra, jbra, lket, jket, 0, 1, 1) * &
          & sqrt(dble(2*lket+1)) * sqrt(6.d0) * dble(zket)
      return
    case('E1')
      if(zbra /= zket) return
      if(lbra == lket) return
      if( iabs(nbra-nket) > 1 ) return
      if( iabs(lbra-lket) > 1 ) return
      if(.not. present(isospin) .and. zket == 1) return

      r = j_to_ls_1( lbra, jbra, lket, jket, 1, 0, 1) * &
          & red_r_l(nbra,lbra,nket,lket) * sqrt(2.d0) * &
          & sqrt(b2) * sqrt(3.d0 /(4.d0 * pi)) * e
      return
    case('E1cm')
      if(zbra /= zket) return
      if(lbra == lket) return
      if( iabs(nbra-nket) > 1 ) return
      if( iabs(lbra-lket) > 1 ) return
      if(.not. present(isospin) .and. zket == 1) return
      r = j_to_ls_1( lbra, jbra, lket, jket, 1, 0, 1) * &
          & red_r_l(nbra,lbra,nket,lket) * sqrt(2.d0) * &
          & sqrt(b2) * sqrt(3.d0 /(4.d0 * pi)) * e
      return
    case('GamowTeller')
      if(zbra == zket) return
      if(lbra /= lket) return
      if(nbra /= nket) return
      if( iabs(jbra-jket) > 2 ) return
      r = j_to_ls_1( lbra, jbra, lket, jket, 0, 1, 1 ) * &
          & sqrt(dble(2*lket+1)) * sqrt(6.d0) * sqrt(2.d0) * dble(zket)
    case('Fermi')
      if(zbra == zket) return
      if(lbra /= lket) return
      if(nbra /= nket) return
      if(jbra /= jket) return
      r = j_to_ls_1( lbra, jbra, lket, jket, 0, 0, 0 ) * &
          & sqrt(dble(2*lket+1)) * sqrt(2.d0) * sqrt(2.d0) * dble(zket)
    case('DFermi','DGamowTeller0','DGamowTeller2')
      return
    case default
      if( s%find(this%OpName,s%str("Sp_E"))) then
        if(triag(jbra, jket, 2*this%jr)) return
        if(zbra /= zket) return
        if(.not. present(isospin) .and. zket == 1) return
        r = single_particle_electric_multipole(nbra, lbra, jbra, nket, lket, jket, this%jr, sqrt(b2)) * e
        return
      end if

      if( s%find(this%OpName,s%str("Sp_M"))) then
        if(triag(jbra, jket, 2*this%jr)) return
        if(zbra /= zket) return
        r = single_particle_magnetic_multipole(nbra, lbra, jbra, nket, lket, jket, zket, this%jr, sqrt(b2), &
          & isospin, e_charge)
        return
      end if

      write(*,'(2a)') 'In CalcMEOneBody, selected operator has not been implemented. ', trim(OpName%val)
    end select
  end function CalcMEOneBody

  function CalcMERel(this, bra, ket) result(r)
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(5), ket(5)
    real(8) :: r

    if(this%pn()) then
      r = get_two_body_me_pn(this, bra, ket)
      return
    end if

    if(.not. this%pn()) then
      r = get_two_body_me_iso(this, bra, ket)
      return
    end if
  end function CalcMERel

  function CalcMERelCM(this, bra, ket) result(r)
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(8), ket(8)
    real(8) :: r

    if(this%pn()) then
      r = get_two_body_me_pn_cm(this, bra, ket)
      return
    end if

    if(.not. this%pn()) then
      r = get_two_body_me_iso_cm(this, bra, ket)
      return
    end if
  end function CalcMERelCM

  recursive function get_two_body_me_pn(this, bra, ket) result(r)
    !
    ! Here, users can define their own operators
    ! In this function, the HO matrix of operator. For non-scalar operator, use the reduced matrix element.
    !
    ! bra and ket are arrays
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,z1]
    ! ket = [n2,l2,s2,j2,z2]
    ! Note:
    ! | Tz= 1 > = | nlSJ:nn > (J have to be even)
    ! | Tz= 0 > = ( | nlSJ:pn > + (-1)^(l+S) | nlSJ:np > ) / sqrt(2)
    ! | Tz=-1 > = | nlSJ:pp > (J have to be even)
    use MyLibrary
    use M1_2b_current
    use ClassSys, only: sys
    real(8) :: r
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(5), ket(5)
    integer :: j1, tz1, n1, l1, s1, lam
    integer :: j2, tz2, n2, l2, s2
    integer :: NMesh = 500
    real(8), allocatable :: p(:), w(:)
    real(8) :: anu, integral
    integer :: ibra
    type(str), allocatable :: strings(:)
    type(str) :: OpName
    type(sys) :: csys
    real(8) :: hw_pot, b2

    n1 = bra(1); n2 = ket(1)
    l1 = bra(2); l2 = ket(2)
    s1 = bra(3); s2 = ket(3)
    j1 = bra(4); j2 = ket(4)
    tz1= bra(5); tz2= ket(5)
    r = 0.d0
    b2 = hc**2 / (m_red_pn * parameters%hw)
    OpName = this%GetOpName()

    select case(OpName%val)
    case('kinetic', 'Kinetic')
      ! kinetic energy term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw)
      return

    case('T_radial')
      ! kinetic energy term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) - &
          & kin_rotation_part(n1,l1,n2,l2,parameters%hw)
      return

    case('T_Gl')
      ! kinetic energy term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) * exp(- (dble(l1)/this%op_args(1))**2)
      return

    case('T_Gj')
      ! kinetic energy term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) * exp(- (dble(j1)/this%op_args(1))**2)
      return

    case('T_l0')
      ! kinetic energy term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      if(l1 /= l2) return
      if(n1 == n2) r = dble(2 * n1) + 1.5d0
      if(n1 == n2-1) r = dsqrt(dble(n1 + 1) * (dble(n1) + 1.5d0))
      if(n1 == n2+1) r = dsqrt(dble(n2 + 1) * (dble(n2) + 1.5d0))
      r = r * parameters%hw * 0.5d0
      return

    case('hopot', 'HOpot')
      ! HO potential term
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_hopot(n1,l1,n2,l2,parameters%hw)
      return

    case('R2', 'r2')
      ! r^{2} operator
      if(tz1 /= tz2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_r2(n1,l1,n2,l2,parameters%hw)
      return

    case('Sigma', 'Spin')
      ! (sigma_{1,z} + sigma_{2,z})
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & tau1_plus_tau2_iso( s1, s2 ) * sqrt(dble(2*l1+1))
      return

    case('Sigma_Tauz')
      ! (sigma_{1,z} tau_{1,z} + sigma_{2,z} tau_{2,z})
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * &
          & (tau1_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_m,1,0) + &
          &  tau2_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau2_m,1,0))
      return

    case('E1')
      ! E1 = - (tau1z - tau2z) * (z1 - z2) / 2
      if(l1 == l2) return
      if(iabs(l1 - l2) > 1) return
      if(iabs(n1 - n2) > 1) return
      if(s1 /= s2) return
      if(tz1 /= tz2) return
      r = - j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1 ) * &
          & red_r_l(n1, l1, n2, l2) * sqrt(dble(2*s1+1)) * &
          & asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_minus_tau2,1,0) * &
          & sqrt(b2) * 0.5d0 * sqrt(3.d0 / (4.d0*pi))
      return
    case('E1cm')
      return
    case("M1_2BC_intr")
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current_intr(n1, l1, s1, j1, tz1, n2, l2, s2, j2, tz2, this%pn())
      return

    case('M1_L_IS')
      ! < n'l'S'J'Tz || L || nlSJTz >
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(s1 /= s2) return
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1) * &
          & sqrt(dble(2*s1+1)) * sqrt( dble( l1 * (l1+1) * (2*l1+1) )) * sqrt(3.d0 / (4.d0 * pi))
    case('M1_L_IV')
      ! < n'l'S'J'Tz || (tau1,z + tau2,z) L || nlSJTz >
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(s1 /= s2) return
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1) * &
          & sqrt(dble(2*s1+1)) * sqrt( dble( l1 * (l1+1) * (2*l1+1) )) * &
          & asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_plus_tau2,1,0) * &
          & sqrt(3.d0 / (4.d0*pi))
    case('M1_S_IS')
      ! < n'l'S'J'Tz || s1 + s2 || nlSJTz > x 0.880
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(tz1 /= tz2) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & tau1_plus_tau2_iso( s1, s2 ) * sqrt(dble(2*l1+1)) * &
          & gs * 0.5d0 * sqrt(3.d0 / (4.d0*pi))
      return
    case('M1_S_IV')
      ! < n'l'S'J'Tz || s1 tau1,z + s2 tau2,z || nlSJTz > x -4.706
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(tz1 /= tz2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * &
          & (tau1_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_m,1,0) + &
          &  tau2_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau2_m,1,0)) * &
          & gv * (-1.d0) * 0.5d0 * sqrt(3.d0 / (4.d0*pi))
      return
    case('E2_IS')
      ! < n1 l1 s1 j1 tz1 || r^2 Y^2 || n2 l2 s2 j2 tz2 >
      if(tz1 /= tz2) return
      if(s1 /= s2) return
      if((-1)**(l1+l2) /= 1) return
      if(triag(j1,j2,2)) return
      if(triag(l1,l2,2)) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 2, 0, 2) * &
          & sqrt(dble(2*s1+1)) * red_l_Y(l1, 2, l2) * &
          & radius_power(2, n1, l1, n2, l2) * b2
      return
    case('E2_IV')
      ! < n1 l1 s1 j1 tz1 || r^2 Y^2 (tau1,z + tau2,z) || n2 l2 s2 j2 tz2 >
      if(tz1 /= tz2) return
      if(s1 /= s2) return
      if((-1)**(l1+l2) /= 1) return
      if(triag(j1,j2,2)) return
      if(triag(l1,l2,2)) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 2, 0, 2) * &
          & sqrt(dble(2*s1+1)) * red_l_Y(l1, 2, l2) * &
          & radius_power(2, n1, l1, n2, l2) * b2 * &
          & asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_plus_tau2,1,0)
      return
    case('GamowTeller')
      ! 1/g_A A -> - sigma tau / 2 (Q->0)
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(abs(j1 - j2) > 1) return
      if(tz2 == tz1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * &
          & (tau1_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_m,1,tz1-tz2) + &
          &  tau2_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau2_m,1,tz1-tz2)) * (-0.5d0)
      return

    case('Fermi')
      ! rho -> tau/2 (Q->0)
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(s1 /= s2) return
      if(j1 /= j2) return
      if(tz2 == tz1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 0, 0 ) * &
          & sqrt(dble(2*l1+1)) * sqrt(dble(2*s1+1)) * &
          & asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_plus_tau2,1,tz1-tz2) * 0.5d0
      return
    case('DFermi')
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(s1 /= s2) return
      if(j1 /= j2) return
      if(abs(tz1-tz2) /= 2) return
      r = sqrt(dble(2*j1+1)) * 2.d0
      return
    case('DGamowTeller0')
      lam = 0
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(triag(s1,s2,lam)) return
      if(triag(j1,j2,lam)) return
      if(abs(tz1-tz2) /= 2) return
      lam = 0
      r = 12.d0 * j_to_ls( l1, s1, j1, l2, s2, j2, 0, lam, lam ) * &
          & sqrt( dble(2*lam+1) * dble(2*l1+1) * dble(2*s1+1) * dble(2*s2+1) ) * &
          & snj( 1, 1, 2, 1, 1, 2, 2*s1, 2*s2, 2*lam )
      return
    case('DGamowTeller2')
      lam = 2
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(triag(s1,s2,lam)) return
      if(triag(j1,j2,lam)) return
      if(abs(tz1-tz2) /= 2) return
      r = 12.d0 * j_to_ls( l1, s1, j1, l2, s2, j2, 0, lam, lam ) * &
          & sqrt( dble(2*lam+1) * dble(2*l1+1) * dble(2*s1+1) * dble(2*s2+1) ) * &
          & snj( 1, 1, 2, 1, 1, 2, 2*s1, 2*s2, 2*lam )
      return

    case("TDM_mag")
      ! NOTE: derived from -(1/4) int dx x^2 j(x)
      ! [ g1 r x (s1-s2) - g2 r x (s1t1 - s2t2) ]/8m
      ! -sqrt(2.0) comes from r x s = -i sqrt(2.0) [r s]^1
      ! i is not taken into account here
      if( l1==l2 ) return
      if(iabs(l1-l2) > 1) return
      if(iabs(n1-n2) > 1) return
      if(tz1 /= tz2) return
      r = - sqrt(2.d0) * hc / (4.d0 * (m_proton+m_neutron)) * &
          & j_to_ls( l1, s1, j1, l2, s2, j2, 1, 1, 1) * &
          & red_r_l(n1,l1,n2,l2) * &
          & sqrt(b2) * &
          & (tau1_minus_tau2_iso(s1,s2) * gs - &
          & (tau1_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_m,1,0) - &
          &  tau2_iso(s1,s2) * asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau2_m,1,0)) * gv)

    case("TDM_conv")
      ! - e[ r + r2 nabla ]/8m (tau1z - tau2z)
      if( l1==l2 ) return
      if( s1/=s2 ) return
      if(iabs(l1-l2) > 1) return
      if(iabs(n1-n2) > 2) return
      if(tz1 /= tz2) return
      r = - 1.d0 * hc / (4.d0 * (m_proton+m_neutron)) * sqrt(dble(2*s1+1)) * &
          & j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1 ) * &
          & ( red_r_l(n1,l1,n2,l2) + red_r2_nab_l(n1,l1,n2,l2) ) * &
          & sqrt(b2) * &
          & asym_isospin_func_pn(l1,s1,tz1,l2,s2,tz2,tau1_minus_tau2,1,0)

    case default

      if( csys%find(OpName,csys%str("Gaussian"))) then
        if( j1 /= j2 ) return
        if( s1 /= s2 ) return
        if(tz1 /=tz2 ) return
        if( l1 /= l2 ) return
        call gauss_legendre(0.d0, 25.d0, p, w, NMesh)
        anu = hc**2 / (m_red_pn * parameters%hw)
        r = 0.d0
        integral = 0.d0

        ! < n'l | f(p) | nl >
        do ibra = 1, NMesh
            r = r + w(ibra)*exp( -p(ibra)**2 / (this%op_args(1)/hc)**2 ) * &
                & ho_radial_wf_norm(n1,l1,anu,p(ibra)) * ho_radial_wf_norm(n2,l2,anu,p(ibra))
        end do
        r = r * (-1)**(n1+n2)
        return
      end if

      if( csys%find(this%OpName, csys%str("HOHamil") )) then
        call csys%split(this%OpName, csys%str("Hamil"), strings )
        if(tz1 /= tz2) return
        if(l1 /= l2) return
        if(j1 /= j2) return
        hw_pot = parameters%hw
        if(size(strings)>1) read(strings(2)%val,*) hw_pot
        r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) + &
            & two_body_hopot(n1,l1,n2,l2,hw_pot)
        return
      end if


      write(*,'(2a)') 'In get_twobody_me_pn, selected operator has not been implemented. ', trim(OpName%val)
      stop
    end select
  end function get_two_body_me_pn

  function get_two_body_me_pn_cm(this, bra, ket) result(r)
    !
    ! Here, users can define their own operators
    ! In this function, the HO matrix of operator. For non-scalar operator, use the reduced matrix element.
    !
    ! bra and ket are arrays
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,jr1,Ncm1,Lcm1,j1,z1]
    ! ket = [n2,l2,s2,jr2,Ncm2,Lcm2,j2,z2]
    ! Note:
    ! | Tz= 1 > = | nlSJ:nn > (J have to be even)
    ! | Tz= 0 > = ( | nlSJ:pn > + (-1)^(l+S) | nlSJ:np > ) / sqrt(2)
    ! | Tz=-1 > = | nlSJ:pp > (J have to be even)
    use MyLibrary
    use M1_2b_current
    real(8) :: r
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(8), ket(8)
    integer :: n1, l1, s1, jr1, Ncm1, Lcm1, j1, z1
    integer :: n2, l2, s2, jr2, Ncm2, Lcm2, j2, z2
    type(str) :: OpName
    n1   = bra(1); n2   = ket(1)
    l1   = bra(2); l2   = ket(2)
    s1   = bra(3); s2   = ket(3)
    jr1  = bra(4); jr2  = ket(4)
    Ncm1 = bra(5); Ncm2 = ket(5)
    Lcm1 = bra(6); Lcm2 = ket(6)
    j1   = bra(7); j2   = ket(7)
    z1   = bra(8); z2   = ket(8)
    OpName = this%GetOpName()
    select case(OpName%val)
    case("M1_2BC")
      if(z1 /= z2) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current(n1, l1, s1, jr1, ncm1, lcm1, j1, z1, n2, l2, s2, jr2, ncm2, lcm2, j2, z2, this%pn())
      return
    case("M1_2BC_Sachs")
      if(z1 /= z2) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current_sachs(n1, l1, s1, jr1, ncm1, lcm1, j1, z1, n2, l2, s2, jr2, ncm2, lcm2, j2, z2, this%pn())
      return
    end select
  end function get_two_body_me_pn_cm

  recursive function get_two_body_me_iso(this, bra, ket) result(r)
    !
    ! Here, users can define their own operators
    ! In this function, the HO matrix of operator. For non-scalar operator, use the reduced matrix element.
    !
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,j1,z1]
    ! ket = [n2,l2,s2,j2,z2]
    !
    use MyLibrary
    use M1_2b_current
    real(8) :: r
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(5), ket(5)
    integer :: j1, t1, n1, l1, s1, lam
    integer :: j2, t2, n2, l2, s2
    type(str), allocatable :: strings(:)
    real(8) :: s, hw_pot, b2
    type(sys) :: csys
    type(str) :: OpName

    n1 = bra(1); n2 = ket(1)
    l1 = bra(2); l2 = ket(2)
    s1 = bra(3); s2 = ket(3)
    j1 = bra(4); j2 = ket(4)
    t1 = bra(5); t2 = ket(5)
    r = 0.d0
    s = 1.d0
    b2 = hc**2 / (m_red_pn * parameters%hw)

    if(t1 == -1 .or. t2 == -1) return
    if((-1) ** (l1 + s1 + t1) == 1) return
    if((-1) ** (l2 + s2 + t2) == 1) return
    OpName = this%GetOpName()
    select case(OpName%val)

    case('kinetic', 'Kinetic')
      ! kinetic energy term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw)
      return

    case('T_radial')
      ! kinetic energy term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) - &
          & kin_rotation_part(n1,l1,n2,l2,parameters%hw)
      return

    case('T_Gl')
      ! kinetic energy term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) * exp(- (dble(l1)/this%op_args(1))**2)
      return

    case('T_Gj')
      ! kinetic energy term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) * exp(- (dble(j1)/this%op_args(1))**2)
      return

    case('T_l0')
      ! kinetic energy term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      if(l1 /= l2) return
      if(n1 == n2) r = dble(2 * n1) + 1.5d0
      if(n1 == n2-1) r = dsqrt(dble(n1 + 1) * (dble(n1) + 1.5d0))
      if(n1 == n2+1) r = dsqrt(dble(n2 + 1) * (dble(n2) + 1.5d0))
      r = r * parameters%hw * 0.5d0
      return

    case('hopot', 'HOpot')
      ! HO potential term
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_hopot(n1,l1,n2,l2,parameters%hw)
      return

    case( "HOHamil" )
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) + &
          & two_body_hopot(n1,l1,n2,l2,parameters%hw)
      return

    case('R2', 'r2')
      ! r^{2} operator
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_r2(n1,l1,n2,l2,parameters%hw)
      return

    case('R2_IV')
      ! r^{2} operator
      if(t1 /= t2) return
      if(j1 /= j2) return
      if(s1 /= s2) return
      r = two_body_r2(n1,l1,n2,l2,parameters%hw) * sqrt(dble(2*j1+1)) * tau1_plus_tau2_iso(t1,t2)
      return

    case('Spin','Sigma')
      ! (sigma_{1,z} + sigma_{2,z})
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(t1 /= t2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble((2*l1+1)*(2*t1+1))) * (tau1_plus_tau2_iso(s1,s2))
      return

    case('Sigma_Tauz')
      ! (sigma_{1,z} tau_{1,z} + sigma_{2,z} tau_{2,z})
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * (tau1_iso(s1, s2)*tau1_iso(t1, t2) + tau2_iso(s1, s2)*tau2_iso(t1, t2))
      return

    case('E1')
      ! E1 = - (tau1z - tau2z) * (z1 - z2) / 2
      if(l1 == l2) return
      if(iabs(l1 - l2) > 1) return
      if(iabs(n1 - n2) > 1) return
      if(s1 /= s2) return
      if(t1 == t2) return
      r = - j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1 ) * &
          & red_r_l(n1, l1, n2, l2) * sqrt(dble(2*s1+1)) * &
          & tau1_minus_tau2_iso( t1,t2 ) * &
          & sqrt(b2) * 0.5d0 * sqrt(3.d0 / (4.d0 * pi))
      return

    case("M1_2BC_intr")
      if(iabs(t1 - t2) > 1) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current_intr(n1, l1, s1, j1, t1, n2, l2, s2, j2, t2, this%pn())
      return

    case('M1_L_IS')
      ! < n'l'S'J'Tz || L || nlSJTz >
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(s1 /= s2) return
      if(t1 /= t2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1) * &
          & sqrt(dble(2*s1+1)) * sqrt( dble( l1 * (l1+1) * (2*l1+1) * (2*t1+1) )) * sqrt(3.d0 / (4.d0*pi))

    case('M1_L_IV')
      ! < n'l'S'J'Tz || (tau1,z + tau2,z) L || nlSJTz >
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(s1 /= s2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 1, 0, 1) * &
          & sqrt(dble(2*s1+1)) * sqrt( dble( l1 * (l1+1) * (2*l1+1) )) * tau1_plus_tau2_iso(t1,t2) * sqrt(3.d0/(4.d0*pi))
    case('M1_S_IS')
      ! < n'l'S'J'Tz || s1 + s2 || nlSJTz > x 0.880
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(t1 /= t2) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & tau1_plus_tau2_iso(s1,s2) * sqrt(dble((2*l1+1)*(2*t1+1))) * &
          & gs * 0.5d0 * sqrt(3.d0/(4.d0*pi))
      return
    case('M1_S_IV')
      ! < n'l'S'J'Tz || s1 tau1,z + s2 tau2,z || nlSJTz > x -4.706
      if(n1 /= n2) return
      if(l1 /= l2) return
      if(iabs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * &
          & (tau1_iso(s1, s2)*tau1_iso(t1, t2) + tau2_iso(s1,s2)*tau2_iso(t1,t2)) * &
          & gv * (-1.d0) * 0.5d0 * sqrt(3.d0/(4.d0*pi))
      return

    case('E2_IS')
      ! < n1 l1 s1 j1 t1 || r^2 Y^2 || n2 l2 s2 j2 t2 >
      if(t1 /= t2) return
      if(s1 /= s2) return
      if((-1)**(l1+l2) /= 1) return
      if(triag(j1,j2,2)) return
      if(triag(l1,l2,2)) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 2, 0, 2) * &
          & sqrt(dble(2*s1+1)*dble(2*t1+1)) * red_l_Y(l1, 2, l2) * &
          & radius_power(2, n1, l1, n2, l2) * b2
      return
    case('E2_IV')
      ! < n1 l1 s1 j1 t1 || r^2 Y^2 (tau1,z + tau2,z) || n2 l2 s2 j2 t2 >
      if(s1 /= s2) return
      if((-1)**(l1+l2) /= 1) return
      if(iabs(t1-t2)>1) return
      if(triag(j1,j2,2)) return
      if(triag(l1,l2,2)) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 2, 0, 2) * &
          & sqrt(dble(2*s1+1)) * red_l_Y(l1, 2, l2) * &
          & radius_power(2, n1, l1, n2, l2) * b2 * &
          & tau1_plus_tau2_iso(t1,t2)
      return

    case("GamowTeller")
      ! 1/g_A A -> - sigma tau / 2 (Q->0)
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(abs(j1 - j2) > 1) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 1, 1 ) * &
          & sqrt(dble(2*l1+1)) * &
          & (tau1_iso(s1,s2) * tau1_iso(t1,t2) + &
          &  tau2_iso(s1,s2) * tau2_iso(t1,t2)) * (-0.5d0)
      return
    case("Fermi")
      ! rho -> tau / 2 (Q->0)
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(s1 /= s2) return
      if(j1 /= j2) return
      r = j_to_ls( l1, s1, j1, l2, s2, j2, 0, 0, 0 ) * &
          & sqrt(dble(2*l1+1)) * sqrt(dble(2*s1+1)) * &
          & tau1_plus_tau2_iso(t1,t2) * 0.5d0
      return
    case('DFermi')
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(s1 /= s2) return
      if(j1 /= j2) return
      if(triag(t1,t2,2)) return
      r = sqrt(dble(2*j1+1)) * 2.d0 * sqrt(5.d0)
      return
    case('DGamowTeller0')
      lam = 0
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(triag(s1,s2,lam)) return
      if(triag(j1,j2,lam)) return
      if(triag(t1,t2,2)) return
      r = 12.d0 * j_to_ls( l1, s1, j1, l2, s2, j2, 0, lam, lam ) * &
          & sqrt( dble(2*lam+1) * dble(2*l1+1) * dble(2*s1+1) * dble(2*s2+1) ) * &
          & snj( 1, 1, 2, 1, 1, 2, 2*s1, 2*s2, 2*lam ) * sqrt(5.d0)
      return
    case('DGamowTeller2')
      lam = 2
      if(l1 /= l2) return
      if(n1 /= n2) return
      if(triag(s1,s2,lam)) return
      if(triag(j1,j2,lam)) return
      if(triag(t1,t2,2)) return
      r = 12.d0 * j_to_ls( l1, s1, j1, l2, s2, j2, 0, lam, lam ) * &
          & sqrt( dble(2*lam+1) * dble(2*l1+1) * dble(2*s1+1) * dble(2*s2+1) ) * &
          & snj( 1, 1, 2, 1, 1, 2, 2*s1, 2*s2, 2*lam ) * sqrt(5.d0)
      return

    case("TDM_mag_s")
      ! gs r x (sigma_1 - sigma_2) / 8m / i
      ! tauz = -1 (proton), tauz= 1 (neutron)
      if(l1==l2) return
      if(iabs(l1-l2) > 1) return
      if(iabs(n1-n2) > 1) return
      if(t1 /= t2) return
      r = red_r_l(n1,l1,n2,l2) * sqrt(b2) * &
          & 6.d0 * sqrt( dble( (2*j1+1) * (2*j2+1) * (2*s1+1) * (2*s2+1) * (2*t1+1) ) ) * &
          & snj(2*l1, 2*s1, 2*j1, 2*l2, 2*s2, 2*j2, 2, 2, 2) * &
          & ( (-1.d0)**s2 - (-1.d0)**s1 ) * &
          & sjs(1, 1, 2, 2*s2, 2*s1, 1) * hc / (4.d0 * (m_proton+m_neutron)) * gs

    case("TDM_mag_st")
      ! gv r x (sigma_1 tauz_1 - sigma_2 tauz_2) / 8m / i
      ! tauz = -1 (proton), tauz= 1 (neutron)
      if(l1==l2) return
      if(iabs(l1-l2) > 1) return
      if(iabs(n1-n2) > 1) return
      r = red_r_l(n1,l1,n2,l2) * sqrt(b2) * &
          & 6.d0 * sqrt(6.d0 * dble( (2*j1+1) * (2*j2+1) * (2*s1+1) * (2*s2+1) * &
          & (2*t1+1) * (2*t2+1) ) ) * &
          & snj(2*l1, 2*s1, 2*j1, 2*l2, 2*s2, 2*j2, 2, 2, 2) * &
          & sjs(1, 1, 2, 2*s2, 2*s1, 1) * sjs(1, 1, 2, 2*t2, 2*t1, 1) * &
          & ( (-1.d0)**(s2+t2) - (-1.d0)**(s1+t1) ) * hc / (4.d0 * (m_proton+m_neutron) ) * gv

    case default
      if( csys%find(OpName, csys%str("HOHamil") )) then
        call csys%split(OpName, csys%str("Hamil"), strings )
        if(t1 /= t2) return
        if(l1 /= l2) return
        if(j1 /= j2) return
        hw_pot = parameters%hw
        if(size(strings)>1) read(strings(2)%val,*) hw_pot
        r = two_body_kinetic(n1,l1,n2,l2,parameters%hw) + &
            & two_body_hopot(n1,l1,n2,l2,hw_pot)
        return
      end if

      write(*,'(2a)') 'In get_twobody_me_iso, selected operator has not been implemented. ', trim(OpName%val)
      stop
    end select
  end function get_two_body_me_iso

  function get_two_body_me_iso_cm(this, bra, ket) result(r)
    !
    ! Here, users can define their own operators
    ! In this function, the HO matrix of operator. For non-scalar operator, use the reduced matrix element.
    !
    ! bra and ket are arrays
    ! < bra | O | ket >
    ! bra = [n1,l1,s1,jr1,Ncm1,Lcm1,j1,t1]
    ! ket = [n2,l2,s2,jr2,Ncm2,Lcm2,j2,t2]
    use MyLibrary
    use M1_2b_current
    real(8) :: r
    type(OperatorDef), intent(in) :: this
    integer, intent(in) :: bra(8), ket(8)
    integer :: n1, l1, s1, jr1, Ncm1, Lcm1, j1, t1
    integer :: n2, l2, s2, jr2, Ncm2, Lcm2, j2, t2
    type(str) :: OpName
    n1   = bra(1); n2   = ket(1)
    l1   = bra(2); l2   = ket(2)
    s1   = bra(3); s2   = ket(3)
    jr1  = bra(4); jr2  = ket(4)
    Ncm1 = bra(5); Ncm2 = ket(5)
    Lcm1 = bra(6); Lcm2 = ket(6)
    j1   = bra(7); j2   = ket(7)
    t1   = bra(8); t2   = ket(8)
    OpName = this%GetOpName()
    select case(OpName%val)
    case("M1_2BC")
      if(iabs(t1 - t2) > 1) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current(n1, l1, s1, jr1, ncm1, lcm1, j1, t1, n2, l2, s2, jr2, ncm2, lcm2, j2, t2, this%pn())
      return
    case("M1_2BC_Sachs")
      if(iabs(t1 - t2) > 1) return
      if(iabs(j1 - j2) > 1) return
      r = me_M1_2b_current_sachs(n1, l1, s1, jr1, ncm1, lcm1, j1, t1, n2, l2, s2, jr2, ncm2, lcm2, j2, t2, this%pn())
      return
    end select
  end function get_two_body_me_iso_cm

  function two_body_kinetic(n1, l1, n2, l2, hw) result(r)
    !
    ! p^2 / 2m
    !
    integer, intent(in) :: n1, l1, n2, l2
    real(8) :: hw, r
    r = 0.d0
    if(l1 /= l2) return
    if(abs(n1-n2) > 1) return
    if(n1 == n2) r = dble(2 * n1 + l1) + 1.5d0
    if(n1 == n2-1) r = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
    if(n1 == n2+1) r = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
    r = r * hw * 0.5d0
  end function two_body_kinetic

  function two_body_hopot(n1, l1, n2, l2, hw) result(r)
    !
    ! m w^2 r^2 / 2
    !
    integer, intent(in) :: n1, l1, n2, l2
    real(8) :: hw, r
    r = 0.d0
    if(l1 /= l2) return
    if(abs(n1-n2) > 1) return
    if(n1 == n2) r = dble(2 * n1 + l1) + 1.5d0
    if(n1 == n2-1) r = -dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
    if(n1 == n2+1) r = -dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
    r = r * hw * 0.5d0
  end function two_body_hopot

  function two_body_r2(n1, l1, n2, l2, hw) result(r)
    use MyLibrary, only: hc, m_red_pn
    !
    ! (r_i - r_j)^2
    !
    integer, intent(in) :: n1, l1, n2, l2
    real(8) :: hw, r, b2
    r = two_body_hopot(n1, l1, n2, l2, hw)
    b2 = hc**2  / (m_red_pn * hw**2)
    r = r * (2.d0 * b2)
  end function two_body_r2

  !
  ! isospin functions pn
  !
  !function tau_dot_tau_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{1} dot tau_{2} | lket, sket >
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra /= zket) return
  !  r = 1.d0
  !  if(abs(zket) == 1) return
  !  r = (-1.d0)**(lbra+sbra) + (-1.d0)**(lket+sket) - 0.5d0 * (1.d0 + (-1.d0)**(lbra+lket+sbra+sket))
  !end function tau_dot_tau_pn

  !function tau1_tau2_tensor_0_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! 3 tau_{1,z}tau_{2,z} - tau_{1} dot tau_{2}
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 3.d0 * tau1z_tau2z_pn(lbra,sbra,zbra,lket,sket,zket) - &
  !      & tau_dot_tau_pn(lbra,sbra,zbra,lket,sket,zket)
  !end function tau1_tau2_tensor_0_pn

  !function tau1_tau2_vector_0_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! tau_{1,x}tau_{2,y} - tau_{1,y}tau_{2,x}
  !  ! i is not taken into account here
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra /= zket) return
  !  if(abs(zket) == 1) return
  !  r = (-1.d0)**(lbra+sbra) - (-1.d0)**(lket+sket)
  !end function tau1_tau2_vector_0_pn

  !function tau1z(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{1,z} | lket, sket >
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if( zbra /= zket ) return
  !  if( abs(zket) == 1 ) r = dble(zket)
  !  if( zket == 0 ) r = -0.5d0 + 0.5d0 * (-1.d0)**(lbra+sbra+lket+sket)
  !end function tau1z

  !function tau2z(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{2,z} | lket, sket >
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if( zbra /= zket ) return
  !  if( abs(zket) == 1 ) r = dble(zket)
  !  if( zket == 0 ) r = 0.5d0 - 0.5d0 * (-1.d0)**(lbra+sbra+lket+sket)
  !end function tau2z

  !function tau1z_tau2z_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{1,z} tau_{2,z} | lket, sket >
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra /= zket) return
  !  r = 1.d0
  !  if(abs(zket) == 1) return
  !  r = - 0.5d0 * (1.d0 + (1.d0)**(lbra+lket+sbra+sket))
  !end function tau1z_tau2z_pn

  !function tau1z_minus_tau2z_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{1}_z - tau_{2}_z | lket, sket >
  !  integer, intent(in) :: lbra, sbra, zbra, lket, sket, zket
  !  real(8) :: r
  !  r = tau1z(lbra,sbra,zbra,lket,sket,zket) - &
  !      & tau2z(lbra,sbra,zbra,lket,sket,zket)
  !end function tau1z_minus_tau2z_pn

  !function tau1z_plus_tau2z_pn(zbra,zket) result(r)
  !  ! < lbra, sbra| tau_{1}_z + tau_{2}_z | lket, sket >
  !  integer, intent(in) :: zbra, zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra /= zket) return
  !  if(zket == 0) return
  !  r = dble(2*zket)
  !end function tau1z_plus_tau2z_pn

  !function tau1_plus_minus(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{1}_+ | lket, sket >
  !  ! < lbra, sbra| tau_{1}_- | lket, sket >
  !  ! here, < n | tau_{+} | p > = 1 and < p | tau_{-} | n > = 1 are used.
  !  integer, intent(in) :: lbra,sbra,zbra,lket,sket,zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra == zket) return
  !  if(zbra == 0 .and. zket ==-1) r = (-1.d0)**(lbra+sbra) / sqrt(2.d0)
  !  if(zbra == 1 .and. zket == 0) r = 1.d0 / sqrt(2.d0)
  !  if(zbra == 0 .and. zket == 1) r = 1.d0 / sqrt(2.d0)
  !  if(zbra ==-1 .and. zket == 0) r = (-1.d0)**(lket+sket) / sqrt(2.d0)
  !end function tau1_plus_minus

  !function tau2_plus_minus(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  ! < lbra, sbra| tau_{2}_+ | lket, sket >
  !  ! < lbra, sbra| tau_{2}_- | lket, sket >
  !  ! here, < n | tau_{+} | p > = 1 and < p | tau_{-} | n > = 1 are used.
  !  integer, intent(in) :: lbra,sbra,zbra,lket,sket,zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra == zket) return
  !  if(zbra == 0 .and. zket ==-1) r = 1.d0 / sqrt(2.d0)
  !  if(zbra == 1 .and. zket == 0) r = (-1.d0)**(lket+sket) / sqrt(2.d0)
  !  if(zbra == 0 .and. zket == 1) r = (-1.d0)**(lbra+sbra) / sqrt(2.d0)
  !  if(zbra ==-1 .and. zket == 0) r = 1.d0 / sqrt(2.d0)
  !end function tau2_plus_minus

  !!
  !! isospin functions iso
  !!
  !function tau_dot_tau_iso(tbra, tket, reduced) result(r)
  !  ! < lbra, sbra| tau_{1} dot tau_{2} | lket, sket >
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: tbra, tket
  !  logical, intent(in) :: reduced
  !  real(8) :: r
  !  r = 0.d0
  !  if(tbra /= tket) return
  !  r = 6.d0 * sqrt(dble(2*tket+1)) * (-1.d0)**(1+tket) * sjs(1, 1, 2*tket, 1, 1, 2)
  !  if(reduced) return
  !  if(tket == 0) r = -3.d0
  !  if(tket == 1) r =  1.d0
  !end function tau_dot_tau_iso

  !function tau1_tau2_tensor_iso(tbra, tket, Trank) result(r)
  !  !  < Tbra || [tau1 x tau2]^{Trank} || Tket >
  !  use MyLibrary, only: snj
  !  integer, intent(in) :: tbra, tket, Trank
  !  real(8) :: r
  !  r = 6.d0 * sqrt(dble( (2*Tbra+1) * (2*Tket+1) * (2*Trank+1) ) ) * &
  !      & snj(1, 1, 2*Tbra, 1, 1, 2*Tket, 2, 2, 2*Trank)
  !end function tau1_tau2_tensor_iso

  !function tau1z_minus_tau2z_iso(Tbra,Tket) result(r)
  !  ! < Tbra || tau_{1}_z - tau_{2}_z || Tket >
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: Tbra, Tket
  !  real(8) :: r
  !  r = sigma1(Tbra,Tket) - sigma2(Tbra,Tket)
  !end function tau1z_minus_tau2z_iso

  !function tau1z_plus_tau2z_iso(tbra,tket) result(r)
  !  ! < Tbra || tau_{1}_z + tau_{2}_z || Tket >
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: Tbra, Tket
  !  real(8) :: r
  !  r = sigma1(Tbra,Tket) + sigma2(Tbra,Tket)
  !end function tau1z_plus_tau2z_iso

  !!
  !! spin functions
  !!
  !function sigma_dot_sigma(sbra, sket, reduced) result(r)
  !  integer, intent(in) :: sbra, sket
  !  logical, intent(in) :: reduced
  !  real(8) :: r
  !  r = tau_dot_tau_iso(sbra, sket, reduced)
  !end function sigma_dot_sigma

  !function sigma1_sigma2_tensor(sbra, sket, Srank) result(r)
  !  integer, intent(in) :: sbra, sket, Srank
  !  real(8) :: r
  !  r = tau1_tau2_tensor_iso(sbra,sket,Srank)
  !end function sigma1_sigma2_tensor

  !function sigma1(Sbra,Sket) result(r)
  !  ! < Sbra || sigma_{1} || Sket >
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: Sbra, Sket
  !  real(8) :: r
  !  r = sqrt(6.d0) * (-1.d0)**Sket * &
  !      & sqrt(dble( (2*Sbra+1) * (2*Sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*Sket, 2*Sbra, 1)
  !end function sigma1

  !function sigma2(Sbra,Sket) result(r)
  !  ! < Sbra || sigma_{2} || Sket >
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: Sbra, Sket
  !  real(8) :: r
  !  r = sqrt(6.d0) * (-1.d0)**Sbra * &
  !      & sqrt(dble( (2*Sbra+1) * (2*Sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*Sket, 2*Sbra, 1)
  !end function sigma2


  !function sigma1z_minus_sigma2z(Sbra,Sket) result(r)
  !  integer, intent(in) :: Sbra, Sket
  !  real(8) :: r
  !  r = Sigma1( Sbra, Sket ) - Sigma2( Sbra, Sket )
  !end function sigma1z_minus_sigma2z

  !function sigma1z_plus_sigma2z(Sbra,Sket) result(r)
  !  integer, intent(in) :: Sbra, Sket
  !  real(8) :: r
  !  r = Sigma1( Sbra, Sket ) + Sigma2( Sbra, Sket )
  !end function sigma1z_plus_sigma2z

  !function sigtau1z_plus_sigtau2z_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: lbra,sbra,zbra,lket,sket,zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra/=zket) return
  !  r =   Sigma1( Sbra, Sket ) * Tau1z( lbra, sbra, zbra, lket, sket, zket ) + &
  !      & Sigma2( Sbra, Sket ) * Tau2z( lbra, sbra, zbra, lket, sket, zket )
  !end function sigtau1z_plus_sigtau2z_pn

  !function sigtau1z_minus_sigtau2z_pn(lbra,sbra,zbra,lket,sket,zket) result(r)
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: lbra,sbra,zbra,lket,sket,zket
  !  real(8) :: r
  !  r = 0.d0
  !  if(zbra/=zket) return
  !  r =   Sigma1( Sbra, Sket ) * Tau1z( lbra, sbra, zbra, lket, sket, zket ) - &
  !      & Sigma2( Sbra, Sket ) * Tau2z( lbra, sbra, zbra, lket, sket, zket )
  !end function sigtau1z_minus_sigtau2z_pn

  !function sigtau1z_plus_sigtau2z_iso(sbra,tbra,sket,tket) result(r)
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: sbra,tbra,sket,tket
  !  real(8) :: r
  !  r = (-1.d0)**sket * sqrt( dble( 6*(2*sbra+1)*(2*sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*sket, 2*sbra, 1) * &
  !      & (-1.d0)**tket * sqrt( dble( 6*(2*tbra+1)*(2*tket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*tket, 2*tbra, 1) + &
  !      & (-1.d0)**sbra * sqrt( dble( 6*(2*sbra+1)*(2*sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*sket, 2*sbra, 1) * &
  !      & (-1.d0)**tbra * sqrt( dble( 6*(2*tbra+1)*(2*tket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*tket, 2*tbra, 1)
  !end function sigtau1z_plus_sigtau2z_iso

  !function sigtau1z_minus_sigtau2z_iso(sbra,tbra,sket,tket) result(r)
  !  use MyLibrary, only: sjs
  !  integer, intent(in) :: sbra,tbra,sket,tket
  !  real(8) :: r
  !  r = (-1.d0)**sket * sqrt( dble( 6*(2*sbra+1)*(2*sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*sket, 2*sbra, 1) * &
  !      & (-1.d0)**tket * sqrt( dble( 6*(2*tbra+1)*(2*tket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*tket, 2*tbra, 1) - &
  !      & (-1.d0)**sbra * sqrt( dble( 6*(2*sbra+1)*(2*sket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*sket, 2*sbra, 1) * &
  !      & (-1.d0)**tbra * sqrt( dble( 6*(2*tbra+1)*(2*tket+1) ) ) * &
  !      & sjs(1, 1, 2, 2*tket, 2*tbra, 1)
  !end function sigtau1z_minus_sigtau2z_iso

  !
  !
  !
  function j_to_ls( lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank) result(r)
    use MyLibrary, only: snj
    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket, lrank, srank, jrank
    real(8) :: r
    r = sqrt(dble( (2*jbra+1)*(2*jket+1)*(2*jrank+1) ) ) * &
        & snj(2*lbra,2*sbra,2*jbra,2*lket,2*sket,2*jket,2*lrank,2*srank,2*jrank)
  end function j_to_ls

  function j_to_ls_1( lbra, jbra, lket, jket, lrank, srank, jrank) result(r)
    use MyLibrary, only: snj
    integer, intent(in) :: lbra, jbra, lket, jket, lrank, srank, jrank
    real(8) :: r
    r = sqrt(dble( (jbra+1)*(jket+1)*(2*jrank+1) ) ) * &
        & snj(2*lbra,1,jbra,2*lket,1,jket,2*lrank,2*srank,2*jrank)
  end function j_to_ls_1

  function red_r2_nab_l(n1,l1,n2,l2) result(r)
    integer, intent(in) :: n1, l1, n2, l2
    real(8) :: r
    r = 0.d0
    if(l1 == l2-1) then
      if(n1 == n2-1) r = sqrt( dble(l2*n2) * (dble(n2+l2)-0.5d0) * (dble(n2+l2)+0.5d0) )
      if(n1 == n2  ) r =-sqrt( dble(l2) * (dble(n2+l2)-0.5d0)**2 * (dble(n2+l2)+0.5d0) )
      if(n1 == n2+1) r =-sqrt( dble(l2) * dble(n2+1) * dble(n2+2)**2 )
      if(n1 == n2+2) r = sqrt( dble(l2) * dble(n2+1) * dble(n2+2)* (dble(n2+l2)+1.5d0) )
      return
    end if

    if(l1 == l2+1) then
      if(n1 == n2-2) r = sqrt( dble(l2+1) * dble(n2-1) * dble(n2) * (dble(n2+l2)+0.5d0) )
      if(n1 == n2-1) r =-sqrt( dble(l2+1) * dble(n2) * dble(n2-1)**2 )
      if(n1 == n2  ) r =-sqrt( dble(l2+1) * (dble(n2+l2)+1.5d0) * (dble(n2+l2)+2.5d0)**2 )
      if(n1 == n2+1) r = sqrt( dble(l2+1) * dble(n2+1) * (dble(n2+l2)+1.5d0)* (dble(n2+l2)+2.5d0) )
      return
    end if
  end function red_r2_nab_l

  function kin_rotation_part(n1,l1,n2,l2,hw) result(r)
    use MyLibrary, only: ln_gamma
    integer, intent(in) :: n1, l1, n2, l2
    real(8), intent(in) :: hw
    real(8) :: r
    real(8) :: s
    r = 0.d0
    if(l1 /= l2) return
    if( n1 < n2 ) then
      s = 0.5d0*(ln_gamma(dble(n2+1))+ln_gamma(dble(n1+l1)+1.5d0)-ln_gamma(dble(n1+1))-ln_gamma(dble(n2+l1)+1.5d0))
      r = exp(s) * 0.5d0 * hw * dble( l1*(l1+1) )
      return
    end if

    if( n1 >= n2 ) then
      s = 0.5d0*(ln_gamma(dble(n1+1))+ln_gamma(dble(n2+l1)+1.5d0)-ln_gamma(dble(n2+1))-ln_gamma(dble(n1+l1)+1.5d0))
      r = exp(s) * 0.5d0 * hw * dble( l1*(l1+1) )
      return
    end if
  end function kin_rotation_part

  function single_particle_electric_multipole(nbra, lbra, jbra, nket, lket, jket, lam, b) result(r)
    use MyLibrary, only: pi, tjs, triag
    integer, intent(in) :: nbra, lbra, jbra, nket, lket, jket, lam
    real(8), intent(in) :: b
    real(8) :: r
    r = 0.d0
    if((-1)**(lbra + lket + lam) == -1) return
    r = 1/sqrt(4.d0 * pi) * (-1.d0)**((jket-1)/2 + lam) * &
        & sqrt( dble( (jbra+1)*(2*lam+1)*(jket+1) )) * tjs(jbra,jket,2*lam,1,-1,0) * &
        & radius_power(lam, nbra, lbra, nket, lket) * b**lam
  end function single_particle_electric_multipole

  function single_particle_magnetic_multipole(nbra, lbra, jbra, nket, lket, jket, z, lam, b, &
      & isospin, e_charge) result(r)
    use MyLibrary, only: pi, tjs, triag
    integer, intent(in) :: nbra, lbra, jbra, nket, lket, jket, z, lam
    real(8), intent(in) :: b
    integer, intent(in), optional :: isospin
    real(8), intent(in), optional :: e_charge
    real(8) :: r, e
    real(8) :: kappa, g_spin, g_orb
    r = 0.d0
    if((-1)**(lbra + lket + lam + 1) == -1) return
    e = 1.d0
    if(present(e_charge)) e = e_charge
    if(present(isospin)) then
      if(isospin==0) then
        g_spin = 0.880d0
        g_orb = 0.5d0 * e
      else if(isospin==1) then
        g_spin = 4.706d0 * dble(-z)
        g_orb = 0.5d0 * e * dble(-z)
      else
        write(*,*) "Error", __LINE__, __FILE__
        stop
      end if
    else
      if(z == -1) then
        g_spin = 5.586d0
        g_orb = e

      else if(z== 1) then
        g_spin =-3.826d0
        g_orb = 0.d0
      else
        write(*,*) "Error", __LINE__, __FILE__
        stop
      end if
    end if
    kappa = ((-1.d0)**(lbra + (jbra+1)/2) * dble(jbra+1) + (-1.d0)**(lket + (jket+1)/2) * dble(jket+1)) * 0.5d0
    r = 1/sqrt(4.d0 * pi) * (-1.d0)**((jket-1)/2 + lam) * &
        & sqrt( dble( (jbra+1)*(2*lam+1)*(jket+1) )) * tjs(jbra,jket,2*lam,1,-1,0) * &
        & (dble(lam)-kappa) * (g_orb * ( 1.d0 + kappa / dble(lam+1)) - 0.5d0 * g_spin) * &
        & radius_power(lam-1, nbra, lbra, nket, lket) * b**(lam-1)
  end function single_particle_magnetic_multipole

  function red_l_Y(lbra, lam, lket) result(r)
    ! < lbra || Y^lam || lket >
    use MyLibrary, only: pi, tjs
    integer, intent(in) :: lbra, lket, lam
    real(8) :: r
    r = (-1.d0)**lbra * sqrt(dble((2*lbra+1)*(2*lam+1)*(2*lket+1))/(4.d0*pi)) * tjs(2*lbra, 2*lam, 2*lket, 0, 0, 0)
  end function red_l_Y

  function radius_power(k, n1, l1, n2, l2) result(s)
    !
    !  radial integral for the harmonic oscillator wave function
    !
    !  radius_power(k, n1, l1, n2, l2) =  <n1, l1|r^k/b^k|n2, l2>
    !
    !  n1, n2: the number of nodes (n1, n2 = 0, 1, 2, ...)
    !  integral of r only
    !
    integer, intent(in) :: k, n1, l1, n2, l2
    real(8) :: s
    integer :: ll1, ll2, ll, i, imin, imax

    ll = l1 + l2 + k
    ll1 = l2 - l1 + k
    ll2 = l1 - l2 + k
    s = 0.0d0
    if (mod(ll, 2) == 1 .or. ll1 < 0 .or. ll2 < 0) then
      ! direct integral should be done instead
      s = numerical_integral()
      return
    end if


    ll1 = ll1/2
    ll2 = ll2/2

    imin = max(0, n1-ll1, n2-ll2)
    imax = min(n1, n2)
    do i = imin, imax
      s = s + dbinomial(ll1, n1-i) * dbinomial(ll2, n2-i) &
          & * (double_factorial(ll+2*i+1)/double_factorial(2*i))
    end do

    s = s * sqrt(double_factorial(2*n1)/double_factorial(2*n1+2*l1+1)) &
        & * sqrt(double_factorial(2*n2)/double_factorial(2*n2+2*l2+1)) &
        & * sqrt(1.0d0/2**k) * (-1)**(n1-n2)
  contains
    function numerical_integral() result(r)
      use MyLibrary
      integer :: NMesh, i
      real(8) :: xmax
      real(8), allocatable :: x(:), w(:)
      real(8) :: r
      NMesh = 200
      xmax = 25.d0
      call gauss_legendre(0.d0, xmax, x, w, NMesh)
      r = 0.d0
      do i = 1, NMesh
        r = r + w(i) * ho_radial_wf_norm(n1,l1,1.d0,x(i)) * ho_radial_wf_norm(n2,l2,1.d0,x(i)) * x(i)**k
      end do
      deallocate(x,w)
    end function numerical_integral
  end function radius_power

  function double_factorial(n) result(s)
    !
    !  double factorial n!!
    !
    !  output: double precision (not an integer)
    !
    integer, intent(in) :: n
    real(8) :: s
    integer :: i

    if (n > 300) then
       write(*,'(1a,1i5,1x,1a)') &
            & 'error [double_factorial]: n =', n, 'is too large'
       stop
    end if

    s = 1.0d0
    do i = n, 1, -2
       s = s * i
    end do

  end function double_factorial

  function dbinomial(n, m) result(s)
    !
    !  binomial coefficient: n_C_m
    !  s: double precision
    !
    integer, intent(in) :: n, m
    real(8) :: s, s1, s2
    integer :: i, m1

    s = 1.0d0
    m1 = min(m, n-m)
    if (m1 == 0) return
    if (n > 1000) then
       write(*,'(1a, 1i6, 1a)') '[dbinomial]: n =', n, ' is too large'
       stop
    end if

    if (n < 250) then
       s1 = 1.0d0
       s2 = 1.0d0
       do i = 1, m1
          s1 = s1 * (n-i+1)
          s2 = s2 * (m1-i+1)
       end do
       s = s1 / s2
    else
       do i = 1, m1
          s = (s * (n-i+1)) / (m1-i+1)
       end do
    endif

  end function dbinomial

  subroutine parse_operator_string(op_string, OpName, order, regulator, regulator_power, regulator_cutoff)
    !
    ! (OpName)-(Order)-(Regulator)(Regulator Power)-(Regulator Cutoff)
    !
    type(str), intent(in) :: op_string
    type(str), intent(out) :: OpName, regulator
    integer, intent(out) :: order, regulator_power
    real(8), intent(out) :: regulator_cutoff
    type(str), allocatable :: strings(:)
    type(sys) :: s
    integer :: i

    call s%split(op_string, s%str("-"), strings)
    OpName = strings(1)
    order = -1
    regulator = "None"
    regulator_power = 2
    regulator_cutoff = 500.d0
    if(size(strings)>1) then
      do i = 2, size(strings)
        if(strings(i)%val=="LO")   order=0
        if(strings(i)%val=="NLO")  order=1
        if(strings(i)%val=="N2LO") order=2
        if(strings(i)%val=="N3LO") order=3
        if(strings(i)%val=="N4LO") order=4
        if(s%find(strings(i), s%str("NonLocal"))) then
          regulator = "NonLocal"
          read(strings(i)%val(9:),*) regulator_power
          regulator_cutoff = s%dbl(strings(i+1))
        else if(s%find(strings(i), s%str("Local"))) then
          regulator = "Local"
          read(strings(i)%val(6:),*) regulator_power
          regulator_cutoff = s%dbl(strings(i+1))
        end if
      end do
    end if
    deallocate(strings)
  end subroutine parse_operator_string

  subroutine parse_opname(OpName, OpNameBase, par_names, par_vals)
    !
    ! (OpNameBase)_(prameter name1)_(parameter1)_(parameter name2)_(parameter2)_...
    !
    type(str), intent(in) :: OpName
    type(str), intent(out) :: OpNameBase
    type(str), intent(out), allocatable :: par_names(:), par_vals(:)
    type(str), allocatable :: strings(:)
    type(sys) :: s
    integer :: i, n_par

    call s%split(OpName, s%str("_"), strings)
    OpNameBase = strings(1)
    n_par = (size(strings)-1)/2
    if(n_par==0) return
    allocate(par_names(n_par), par_vals(n_par))
    do i = 1, n_par
      par_names= strings(2*i)
      par_vals = strings(2*i+1)
    end do
    deallocate(strings)
  end subroutine parse_opname

  subroutine set_pv_couplings(pv_couplings, pot_type, LECname)
    real(8), intent(inout) :: pv_couplings(:)
    character(*), intent(in) :: pot_type, LECname
    if(pot_type=="ChEFT") then
      if(LECname=="hpi1") pv_couplings(1 ) = 1.d0
      if(LECname=="Ct1")  pv_couplings(2 ) = 1.d0
      if(LECname=="Ct2")  pv_couplings(3 ) = 1.d0
      if(LECname=="Ct3")  pv_couplings(4 ) = 1.d0
      if(LECname=="Ct4")  pv_couplings(5 ) = 1.d0
      if(LECname=="Ct5")  pv_couplings(6 ) = 1.d0
      if(LECname=="hV0")  pv_couplings(7 ) = 1.d0
      if(LECname=="hV1")  pv_couplings(8 ) = 1.d0
      if(LECname=="hV2")  pv_couplings(9 ) = 1.d0
      if(LECname=="hA1")  pv_couplings(10) = 1.d0
      if(LECname=="hA2")  pv_couplings(11) = 1.d0
      return
    end if

    if(pot_type=="DDH") then
      if(LECname=="hpi1")   pv_couplings(1) = 1.d0
      if(LECname=="hrho0")  pv_couplings(2) = 1.d0
      if(LECname=="hrho1")  pv_couplings(3) = 1.d0
      if(LECname=="hrho1p") pv_couplings(4) = 1.d0
      if(LECname=="hrho2")  pv_couplings(5) = 1.d0
      if(LECname=="homega0")pv_couplings(6) = 1.d0
      if(LECname=="homega1")pv_couplings(7) = 1.d0
      if(LECname=="best") then
        ! "best value" set is from Ann. Phys. (N. Y). 124, 449 (1980).
        pv_couplings(1:7) = [12.d0, -30.d0, -0.5d0, 0.d0, -25.d0, -5.d0, -3.d0]
        pv_couplings(:) = pv_couplings(:) * 3.8d0 * 1.d-8
      end if
      return
    end if
  end subroutine set_pv_couplings

  subroutine set_pvtv_couplings(pvtv_couplings, pot_type, LECname)
    real(8), intent(inout) :: pvtv_couplings(:)
    character(*), intent(in) :: pot_type, LECname
    if(pot_type=="ChEFT") then
      if(LECname=="gpi0")  pvtv_couplings(1) = 1.d0
      if(LECname=="gpi1")  pvtv_couplings(2) = 1.d0
      if(LECname=="gpi2")  pvtv_couplings(3) = 1.d0
      if(LECname=="delta") pvtv_couplings(4) = 1.d0
      if(LECname=="Ct1")   pvtv_couplings(5) = 1.d0
      if(LECname=="Ct2")   pvtv_couplings(6) = 1.d0
      if(LECname=="Ct3")   pvtv_couplings(7) = 1.d0
      if(LECname=="Ct4")   pvtv_couplings(8) = 1.d0
      if(LECname=="Ct5")   pvtv_couplings(9) = 1.d0
      return
    end if

    if(pot_type=="OBE") then
      if(LECname=="gpi0")    pvtv_couplings( 1) = 1.d0
      if(LECname=="gpi1")    pvtv_couplings( 2) = 1.d0
      if(LECname=="gpi2")    pvtv_couplings( 3) = 1.d0
      if(LECname=="grho0")   pvtv_couplings( 4) = 1.d0
      if(LECname=="grho1")   pvtv_couplings( 5) = 1.d0
      if(LECname=="grho2")   pvtv_couplings( 6) = 1.d0
      if(LECname=="geta0")   pvtv_couplings( 7) = 1.d0
      if(LECname=="geta1")   pvtv_couplings( 8) = 1.d0
      if(LECname=="gomega0") pvtv_couplings( 9) = 1.d0
      if(LECname=="gomega1") pvtv_couplings(10) = 1.d0
      return
    end if
  end subroutine set_pvtv_couplings


end module OperatorDefinitions
