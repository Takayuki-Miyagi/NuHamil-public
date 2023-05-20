module TwoBodyManager
  use Profiler, only: timer
  use ClassSys
  use NuHamilInput, only: InputParameters
  use MPIFunction, only: myrank
  use omp_lib
  implicit none

contains
  subroutine TwoBodyManagers()
    use NuHamilInput, only: params
    use NNForce, only: phase_shift_analysis
    type(InputParameters) :: input

    input = params
    if(myrank == 0 .and. params%particle_rank == 2) then
      write(*,'(a)') '############################################################'
      write(*,'(a)') '           Calculations for NN sector'
      write(*,'(a)') '############################################################'
      write(*,*)
    end if

    if(.not. input%pn_form) then
      write(*,*) 'NN calculations with isospin formalism is not supported.'
      write(*,*) 'pn_form is turned to .true.'
      input%pn_form = .true.
    end if

    if(input%file_convert) then
      call file_convert_nn(input)
      return
    end if

    if(allocated(input%files_combined) .and. allocated(input%weights_combined) &
        &.and. size(input%files_combined)==size(input%weights_combined)) then
      call file_combine_nn(input)
      return
    end if

    if(input%file_name_phase_shift%val /= 'none') then
      call phase_shift_analysis(input)
      return
    end if

    if( input%file_nn_cartesian%val /= "none" ) then
      call Vrel_Cartesian( input )
      return
    end if

    if( input%NNrel_op_output%val /= "none" ) then
      call Write_NN_op_file(input)
      return
    end if

    if( input%NNrel_op_input%val /= "none" ) then
      call TMtransformation_from_file(input)
      return
    end if

    if( input%output_nn_file%val /= "none" ) then
      call nn_srg_analysis(input)
      return
    end if


    if(input%trans2lab) then

      call manage_trans_to_lab(input)

    else

      call calc_two_body_system(input)

    end if
  end subroutine TwoBodyManagers

  subroutine manage_trans_to_lab(params)
    use SingleParticleState
    use TwoBodyLabSpacePN
    use TwoBodyTransCoef, only: TransRel2LabSpace
    type(InputParameters), intent(in) :: params
    type(Orbits) :: sps
    type(TwoBodyLabPNSpace) :: lab
    type(TransRel2LabSpace) :: rel2lab
    type(sys) :: s
    integer :: i

    if(params%count_memory) then
      call count_memory_2bme(params)
      return
    end if

    do i = 1, size(params%Operators)
      if(s%find(params%Operators(i), s%str("L5_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tel5_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tmag5_2B")) .or. &
          & s%find(params%Operators(i), s%str("L_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tel_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tmag_2B")) .or. &
          & s%find(params%Operators(i), s%str("M1_2BC_Sachs"))) then
          call trans_to_lab_relcm_op(params%Operators(i))
          cycle
      end if

      select case(params%Operators(i)%val)
      case("Hcm","hcm")
        if(params%renorm%val/="bare") then
          write(*,*) "Only bare", params%Operators(i)%val
          stop
        end if
        call cm_operator(params%Operators(i))
      case default
        call trans_to_lab_op(params%Operators(i))
      end select
    end do

    call rel2lab%fin()
    call lab%fin()
    call sps%fin()

  contains
    subroutine trans_to_lab_op(oprtr)
      use TwoBodyRelativeSpace
      use NNForce, only: NNForceHO
      use TwoBodyRelOps
      use TwoBodyLabOps
      type(str), intent(in) :: oprtr
      type(TwoBodyRelSpaceHOBasis) :: rel, rel_for_trans
      type(TwoBodyRelSpaceSpinHOBasis) :: relspin
      type(NNForceHO) :: vnnspin, Uspin
      type(TwoBodyRelOp) :: vnn, U, oprel, optmp
      type(TwoBodyLabOp) :: oplab
      type(InputParameters) :: input
      type(str) :: filename
      type(sys) :: s
      logical :: ex


      if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
      if(.not. lab%is_Constructed) call lab%init(params%hw, sps, params%e2max, params%snt_mass)
      call oplab%init(lab,oprtr)

      filename = oplab%GetFile(params%file_name_nn, &
          & params%NNInt, params%renorm, params%lambda, params%hw, &
          & params%emax, params%e2max, params%coul, &
          & params%snt_mass)

      inquire(file=filename%val, exist=ex)
      if(ex) then
        call oplab%fin()
        call lab%fin()
        call sps%fin()
        write(*,'(2a)') trim(filename%val), ' already exists.'
        return
      end if

      if(.not. rel2lab%is_Constructed) then
        call rel2lab%init(lab, params%e2max, params%Jmax2)
        write(*,"(a,f12.6,a)") "Estimated memory for two-body Tcoefs: ", rel2lab%GetMemory(), " GB"
      end if

      call relspin%init(params%hw, params%N2max, params%J2max_NNint)
      call vnnspin%init(relspin)
      call Uspin%init(relspin)

      input = params; input%N_ls = input%e2max
      call vnnspin%setNNForceHO(Uspin, input)
      call rel%init(params%hw, params%N2max, params%Jmax2)
      call vnn%init(rel,rel,s%str('hamil'),params%pn_same_mass)
      call U%init(rel,rel,s%str('UT'),params%pn_same_mass)
      call vnn%SetTwoBodyScalarSpinBreaking(vnnspin)
      call U%SetTwoBodyScalarSpinBreaking(Uspin)

      call vnnspin%fin()
      call Uspin%fin()
      call relspin%fin()

      if(oprtr%val == 'hamil' .or. oprtr%val == 'Hamil') then

        oprel = vnn

      else

        call oprel%init(rel, rel, oprtr, params%pn_same_mass)
        if(s%find(oprtr,s%str("PVTC")) .or. &
            & s%find(oprtr,s%str("PVTV")) .or. &
            & s%find(oprtr,s%str("AxialV")) .or. &
            & s%find(oprtr,s%str("0vbb"))) then
          call oprel%set(params%NNint, pmax=params%pmax2, NMesh=params%NMesh2)
          !call oprel%prt()
        else
          call oprel%set()
        end if
        call oprel%UT(U,U)

      end if


      call vnn%fin()
      call U%fin()
      call rel_for_trans%init(params%hw, params%e2max, params%Jmax2)
      oprel = oprel%truncate(rel_for_trans, rel_for_trans)
      call rel%fin()
      if(params%spin_tensor_decomposition /= -1) then
        select case(params%spin_tensor_decomposition)
        case(0, 1, 2)
          optmp = oprel%SpinTensorDecomposition(params%spin_tensor_decomposition)
          oprel = optmp
        case(3) ! C + LS
          optmp = oprel%SpinTensorDecomposition(0) + oprel%SpinTensorDecomposition(1)
          oprel = optmp
        case(4) ! C + T
          optmp = oprel%SpinTensorDecomposition(0) + oprel%SpinTensorDecomposition(2)
          oprel = optmp
        case(-100) ! Total = C + LS + T
          optmp = oprel%SpinTensorDecomposition(0) + oprel%SpinTensorDecomposition(1) + oprel%SpinTensorDecomposition(2)
          oprel = optmp
        case default
          write(*,*) 'Unknown value of spin_tensor_decomposition', __LINE__, __FILE__
        end select
      end if

      call oplab%TMtrans(rel2lab,oprel)
      if(params%svd_rank_op_lab/=-1) call oplab%SVD(params%svd_rank_op_lab)
      if(params%averaged_file_for_test%val == "none" ) then
        call oplab%writef(filename)
      else
        call oplab%writef(filename, params%averaged_file_for_test)
      end if

      call oprel%fin()
      call oplab%fin()
      call rel_for_trans%fin()
    end subroutine trans_to_lab_op

    subroutine trans_to_lab_relcm_op(oprtr)
      use MyLibrary, only: gauss_legendre, gv, m_nucleon
      use TwoBodyRelativeSpace
      use TwoBodyRelOps
      use TwoBodyRelCMOps
      use TwoBodyLabOps
      use TwoBodyMultiPoleOperator
      type(str), intent(in) :: oprtr
      type(TwoBodyRelCMSpaceHOBasis) :: relcm_ho
      type(TwoBodyRelCMSpaceMBasis) :: relcm_mom
      type(TwoBOdyRelCMOp) :: oprelcm, tmp
      type(TwoBodyLabOp) :: oplab
      real(8), allocatable :: x(:), w(:)
      type(TwoBodyMultiPoleOp) :: mltop_2b
      real(8) :: c6
      integer :: J2maxLab, Lcm2Max

      type(str) :: filename
      logical :: ex
      type(sys) :: s


      if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
      if(.not. lab%is_Constructed) call lab%init(params%hw, sps, params%e2max, params%snt_mass)
      call oplab%init(lab,oprtr)

      filename = oplab%GetFile(params%file_name_nn, &
          & params%NNInt, params%renorm, params%lambda, params%hw, &
          & params%emax, params%e2max, params%coul, &
          & params%snt_mass)

      J2maxLab = 2*params%e2max+1
      Lcm2max = 2*params%e2max
      if(params%J2maxLab /= -1) J2maxLab = params%J2maxLab
      if(params%Lcm2max /= -1) Lcm2Max = params%Lcm2max

      call relcm_ho%init(params%hw, params%e2max, J2maxLab, params%e2max, params%jmax2, LcmMax=Lcm2Max)
      call oprelcm%init(relcm_ho, oprtr)
      write(*,'(a,f10.6,a)') 'Memory estimation for rel-cm operator: ', oprelcm%GetMemory()*2.d0, ' GB'

      if(s%find(params%Operators(i), s%str("L5_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tel5_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tmag5_2B"))) then

        call relcm_mom%init(0.d0, params%pmax2, params%NMesh2, 0.d0, params%pmax2, params%NMesh2, &
            & J2maxLab, params%jmax2, Lcm_Max=Lcm2Max)
        call gauss_legendre(0.d0, params%pmax2, x, w, params%NMesh2)
        call relcm_mom%SetMeshWeight(x, w, x, w)

        c6 = gv / m_nucleon * 1.d3
        c6 = 0.d0
        call mltop_2b%init(oprelcm, relcm_mom, c1=params%c1, c3=params%c3, c4=params%c4, c6=c6, cD=params%cD)
        call tmp%init(relcm_ho, oprtr)
        call mltop_2b%SetHOMatrix(tmp, s%str('c1'))
        oprelcm = tmp
        call tmp%fin()

        call tmp%init(relcm_ho, oprtr)
        call mltop_2b%SetHOMatrix(tmp, s%str('c3'))
        oprelcm = oprelcm + tmp
        call tmp%fin()

        call tmp%init(relcm_ho, oprtr)
        call mltop_2b%SetHOMatrix(tmp, s%str('c4'))
        oprelcm = oprelcm + tmp
        call tmp%fin()

        call tmp%init(relcm_ho, oprtr)
        call mltop_2b%SetHOMatrix(tmp, s%str('cD'))
        oprelcm = oprelcm + tmp
        call tmp%fin()

        call tmp%init(relcm_ho, oprtr)
        call mltop_2b%SetHOMatrix(tmp, s%str('c6'))
        oprelcm = oprelcm + tmp
        call tmp%fin()
        call mltop_2b%fin()
        call relcm_mom%fin()
      else if(s%find(params%Operators(i), s%str("L_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tel_2B")) .or. &
          & s%find(params%Operators(i), s%str("Tmag_2B"))) then
        call relcm_mom%init(0.d0, params%pmax2, params%NMesh2, 0.d0, params%pmax2, params%NMesh2, &
            & J2maxLab, params%jmax2, Lcm_Max=Lcm2Max)
        call gauss_legendre(0.d0, params%pmax2, x, w, params%NMesh2)
        call relcm_mom%SetMeshWeight(x, w, x, w)

        call mltop_2b%init(oprelcm, relcm_mom, c1=params%c1, c3=params%c3, c4=params%c4, c6=c6, cD=params%cD)
        call mltop_2b%SetHOMatrix(oprelcm, s%str(''))

        call mltop_2b%fin()
        call relcm_mom%fin()
      else

        call oprelcm%set()

      end if

      call oprelcm%IsZero()
      if(.not. rel2lab%is_Constructed) then
        call rel2lab%init(lab, params%e2max, params%Jmax2)
        write(*,"(a,f12.6,a)") "Estimated memory for two-body Tcoefs: ", rel2lab%GetMemory(), " GB"
      end if

      call oplab%TMtrans(rel2lab,oprelcm)
      call oplab%writef(filename)

      call oprelcm%fin()
      call relcm_ho%fin()
      call oplab%fin()
    end subroutine trans_to_lab_relcm_op

    subroutine cm_operator(opname)
      use TwoBodyLabOps
      type(str), intent(in) :: opname
      type(TwoBodyLabOp) :: op
      type(str) :: filename
      if(.not. sps%is_Constructed) call sps%init(params%emax, params%lmax)
      if(.not. lab%is_Constructed) call lab%init(params%hw, sps, params%e2max, params%snt_mass)
      call op%init(lab,opname)

      filename = op%GetFile(params%file_name_nn, &
          & params%NNInt, params%renorm, params%lambda, params%hw, &
          & params%emax, params%e2max, params%coul, &
          & params%snt_mass)
      call set_two_body_hcm(op)
      call op%writef(filename)
      call op%fin()
    end subroutine cm_operator
  end subroutine manage_trans_to_lab

  subroutine TMtransformation_from_file(params)
    use SingleParticleState
    use TwoBodyTransCoef, only: TransRel2LabSpace
    use TwoBodyRelOps
    use TwoBodyLabOps
    type(InputParameters), intent(in) :: params
    type(Orbits) :: sps
    type(TwoBodyLabPNSpace) :: lab
    type(TransRel2LabSpace) :: rel2lab
    type(str) :: OpName
    type(TwoBodyRelSpaceHOBasis) :: rel, rel_for_trans
    type(TwoBodyRelOp) :: oprel
    type(TwoBodyLabOp) :: oplab
    type(str) :: filename
    integer :: ichbra, ichket

    OpName = params%Operators(1)
    call sps%init(params%emax, params%lmax)
    call lab%init(params%hw, sps, params%e2max, params%snt_mass)
    call oplab%init(lab,OpName)
    write(*,*) "Doing for ", OpName%val
    call rel2lab%init(lab, params%e2max, params%Jmax2)

    filename = oplab%GetFile(params%file_name_nn, &
        & params%NNInt, params%renorm, params%lambda, params%hw, &
        & params%emax, params%e2max, params%coul, &
        & params%snt_mass)
    call rel%init(params%hw, params%N2max, params%Jmax2)
    call oprel%init(rel, rel, OpName, params%pn_same_mass)
    call oprel%ReadFile( params%NNrel_op_input )

    call rel_for_trans%init(params%hw, params%e2max, params%Jmax2)
    oprel = oprel%truncate(rel_for_trans, rel_for_trans)
    call rel%fin()

    call oplab%TMtrans(rel2lab,oprel)
    call oplab%writef(filename)

    call oprel%fin()
    call oplab%fin()
    call rel_for_trans%fin()
    call rel2lab%fin()
    call lab%fin()
    call sps%fin()
  end subroutine TMtransformation_from_file

  subroutine Write_NN_op_file(params)
    use TwoBodyRelOps
    type(InputParameters), intent(in) :: params
    type(str) :: OpName
    type(TwoBodyRelSpaceHOBasis) :: rel
    type(TwoBodyRelOp) :: oprel
    integer :: i
    type(sys) :: s

    write(*,*)
    write(*,*)
    write(*,*) "Write NNrel files: "
    write(*,*)
    call rel%init(params%hw, params%N2max, params%Jmax2)
    do i = 1, size(params%Operators)
      OpName = params%Operators(i)
      write(*,*) "Doing for ", OpName%val
      call oprel%init(rel, rel, OpName, params%pn_same_mass)
      if(s%find(OpName,s%str("PVTC")) .or. s%find(OpName,s%str("PVTV")) .or. s%find(OpName,s%str("0vbb"))) then
        call oprel%set(params%NNint, pmax=params%pmax2, NMesh=params%NMesh2)
      else
        call oprel%set()
      end if
      call oprel%WriteFile(params%NNrel_op_output)
      call oprel%fin()
    end do
    call rel%fin()
  end subroutine Write_NN_op_file

  subroutine count_memory_2bme(params)
    use SingleParticleState
    use TwoBodyLabSpacePN
    use OperatorDefinitions
    use MyLibrary, only: triag
    type(InputParameters), intent(in) :: params
    type(TwoBodyLabPNSpace) :: lab
    type(Orbits) :: sps
    integer :: i, ichbra, ichket
    integer :: nb, jb, pb, tzb, nk, jk, pk, tzk
    type(OperatorDef) :: opdef
    integer(8) :: nelm
    real(8) :: mem

    call sps%init(params%emax)
    call lab%init(params%hw, sps, params%e2max, params%snt_mass)

    do i = 1, size(params%Operators)
      call opdef%InitOpDef(params%Operators(i), .true.)
      nelm = 0
      mem = 0.d0
      do ichbra = 1, lab%GetNumberChannels()
        nb = lab%jpz(ichbra)%GetNumberStates()
        jb = lab%jpz(ichbra)%GetJ()
        pb = lab%jpz(ichbra)%GetParity()
        tzb= lab%jpz(ichbra)%GetZ()
        do ichket = 1, ichbra
          nk = lab%jpz(ichket)%GetNumberStates()
          jk = lab%jpz(ichket)%GetJ()
          pk = lab%jpz(ichket)%GetParity()
          tzk= lab%jpz(ichket)%GetZ()
          if(triag(jb, jk, opdef%GetOpJ())) cycle
          if(pb * opdef%GetOpP() /= pk) cycle
          if(abs(tzk-tzb) /= opdef%GetOpZ()) cycle
          mem = mem + &
              & dble(lab%jpz(ichbra)%GetNumberStates()) * &
              & dble(lab%jpz(ichket)%GetNumberStates()) * 8.d0 / 1024.d0**3
          nelm = nelm + lab%jpz(ichbra)%GetNumberStates() * lab%jpz(ichket)%GetNumberStates()
        end do
      end do
      write(*,"(3a,i12)") "Matrix-element number of ", params%Operators(i)%val, ": ", nelm
      write(*,"(f18.6,a)") mem, " GB"
    end do

    call lab%fin()
    call sps%fin()

  end subroutine count_memory_2bme

  subroutine calc_two_body_system(params)
    type(InputParameters), intent(in) :: params
    integer :: i

    if(abs(params%tzbra) > params%tbra) stop 'Warning: tz_bra have to be smaller than t_bra'
    if(abs(params%tzket) > params%tket) stop 'Warning: Tz_ket have to be smaller than t_ket'

    do i = 1, size(params%Operators)

      call calc_two_body(params, params%Operators(i))

    end do

  end subroutine calc_two_body_system

  subroutine calc_two_body(params, oprtr)
    use LinAlgLib
    use MyLibrary, only: triag, pi, geometry_part
    use OperatorDefinitions
    use TwoBodyRelOps
    type(InputParameters), intent(in) :: params
    type(str), intent(in) :: oprtr
    type(OperatorDef) :: opdef
    type(DMat) :: bra_states, ket_states, UTbra, UTket
    type(DVec) :: bra_energies, ket_energies
    integer :: jbra, pbra, tzbra
    integer :: jket, pket, tzket
    type(str) :: f
    logical :: calc_bra = .true., ex
    real(8) :: vals(2), vals0(2), vals1(2), vals2(2), vals3(2)
    integer :: A, Z, N
    type(sys) :: s

    jbra = params%jbra
    pbra = params%pbra
    tzbra= params%tzbra

    jket = params%jket
    pket = params%pket
    tzket= params%tzket

    select case(oprtr%val)
    case("EDM","TDM","AlphaD")
      if( params%renorm%val /= "bare" ) return
    case("Mmoment_IA","Qmoment_IA","Radius")
      if(jbra/=jket) return
      if(pbra/=pket) return
      if(tzbra/=tzket) return
    case default
      call opdef%InitOpDef(oprtr, .true.)
      if(triag(jbra,jket,opdef%GetOpJ())) return
      if(pbra * pket * opdef%GetOpP() /= 1) return
      if(abs(tzket-tzbra) /= opdef%GetOpZ() ) return
    end select
    A=2
    Z=(2-2*tzket)/2
    N=(2+2*tzket)/2

    if(abs(tzbra) == 1 .and. pbra == 1 .and. mod(jbra,2) == 1) then
      write(*,"(a,i2,a,i2,a,i2)") "We do not have jbra=", jbra, ", pbra=",pbra,", tzbra=",tzbra
      stop
    end if

    if(abs(tzket) == 1 .and. pket == 1 .and. mod(jket,2) == 1) then
      write(*,"(a,i2,a,i2,a,i2)") "We do not have jket=", jket, ", pket=",pket,", tzket=",tzket
      stop
    end if

    if(jbra == jket .and. pbra == pket .and. tzbra == tzket) calc_bra = .false.
    call solve_two_body_hamil(params,jket,pket,tzket,ket_energies,ket_states, UTbra)
    if(calc_bra) then
      call solve_two_body_hamil(params,jbra,pbra,tzbra,bra_energies,bra_states, UTket)
    else
      UTket = UTbra
      bra_energies = ket_energies
      bra_states = ket_states
    end if

    select case(oprtr%val)
    case("hamil","Hamil")
      write(*,'(f18.8)') ket_energies%v(1)

    case("EDM")
      call calculation_electric_dipole_moment(params, &
          & ket_energies%v(1), ket_states%m(:,1))

    case("TDM")
      call calculation_toroidal_dipole_moment(params, &
          & ket_energies%v(1), ket_states%m(:,1))

    case("AlphaD")
      call calculation_dipole_polarizability(params, &
          & ket_energies%v(1), ket_states%m(:,1))
    case("Mmoment_IA")
      vals0 = calc_exp_vals(s%str("M1_S_IS"), params, bra_states, ket_states, UTbra, UTket)
      vals1 = calc_exp_vals(s%str("M1_S_IV"), params, bra_states, ket_states, UTbra, UTket)
      vals2 = calc_exp_vals(s%str("M1_L_IS"), params, bra_states, ket_states, UTbra, UTket)
      vals3 = calc_exp_vals(s%str("M1_L_IV"), params, bra_states, ket_states, UTbra, UTket)
      vals0 = vals0 * geometry_part(2*jbra,2,2*jket,2*jbra,0,2*jket) / dble(A-1)
      vals1 = vals1 * geometry_part(2*jbra,2,2*jket,2*jbra,0,2*jket) / dble(A-1)
      vals2 = vals2 * geometry_part(2*jbra,2,2*jket,2*jbra,0,2*jket) * dble(2*N) / dble(A)**2
      vals3 = vals3 * geometry_part(2*jbra,2,2*jket,2*jbra,0,2*jket) / dble(A)
      vals = vals0 + vals1 + vals2 - vals3
      vals = vals * sqrt(4.d0 * pi/3.d0)
      write(*,'(a12, 4f18.8)') trim(oprtr%val), bra_energies%v(1), ket_energies%v(1), vals(:)

    case("Qmoment_IA")
      vals0 = calc_exp_vals(s%str("E2_IS"), params, bra_states, ket_states, UTbra, UTket)
      vals1 = calc_exp_vals(s%str("E2_IV"), params, bra_states, ket_states, UTbra, UTket)
      vals0 = vals0 * geometry_part(2*jbra,4,2*jket,2*jbra,0,2*jket) * dble(N) / dble(A)**2
      vals1 = vals1 * geometry_part(2*jbra,4,2*jket,2*jbra,0,2*jket) / dble(2*A)
      vals = vals0 - vals1
      vals = vals * sqrt(16.d0 * pi/5.d0)
      write(*,'(a12,4f18.8)') trim(oprtr%val),bra_energies%v(1), ket_energies%v(1), vals(:)

    case("Radius")
      vals = calc_exp_vals(s%str("R2"), params, bra_states, ket_states, UTbra, UTket)
      vals = sqrt(vals / dble(A)**2)
      write(*,'(a12,4f18.8)') trim(oprtr%val),bra_energies%v(1), ket_energies%v(1), vals(:)

    case default
      vals = calc_exp_vals(oprtr,params,bra_states,ket_states,UTbra,UTket)
      write(*,'(a12, 4f18.8)') trim(oprtr%val), bra_energies%v(1), ket_energies%v(1), vals(:)
    end select

    call bra_states%fin()
    call ket_states%fin()
    call bra_energies%fin()
    call ket_energies%fin()

  end subroutine calc_two_body

  function calc_exp_vals(OpName,params,bra_states,ket_states,UTbra,UTket) result(r)
    use LinAlgLib
    use TwoBodyRelOps
    type(str), intent(in) :: OpName
    type(InputParameters), intent(in) :: params
    type(DMat), intent(in) :: bra_states, ket_states, UTbra, UTket
    real(8) :: r(2)
    integer :: jbra, pbra, tzbra
    integer :: jket, pket, tzket
    real(8) :: expctnv, expctnv_renorm
    type(TwoBodyRelChanHOBasis) :: bra_rel, ket_rel
    type(TwoBodyRelOpChan) :: op, op_renorm
    type(sys) :: s

    jbra = params%jbra
    pbra = params%pbra
    tzbra= params%tzbra

    jket = params%jket
    pket = params%pket
    tzket= params%tzket
    call bra_rel%init(params%N2max, jbra, pbra, tzbra, params%hw)
    call ket_rel%init(params%N2max, jket, pket, tzket, params%hw)

    call op%init(OpName,bra_rel,ket_rel,params%pn_same_mass)
    if( s%find(OpName,s%str("0vbb") )) then
      call op%set([jbra,pbra,tzbra],[jket,pket,tzket], params%NNint, params%pmax2, params%NMesh2)
    else
      call op%set([jbra,pbra,tzbra],[jket,pket,tzket])
    end if
    op_renorm = op
    op_renorm%DMat = UTbra%T() * op%DMat * UTket
    expctnv = dot_product(bra_states%m(:,1), matmul(op%m, ket_states%m(:,1)))
    expctnv_renorm = dot_product(bra_states%m(:,1), matmul(op_renorm%m, ket_states%m(:,1)))
    r(1) = expctnv
    r(2) = expctnv_renorm

    call op%fin()
    call op_renorm%fin()
    call bra_rel%fin()
    call ket_rel%fin()
  end function calc_exp_vals

  subroutine solve_two_body_hamil(params,j,p,tz,energies,states,UT)
    use LinAlgLib
    use NNForce, only: NNForceHO
    use TwoBodyRelOps
    use OperatorDefinitions
    type(InputParameters), intent(in) :: params
    integer, intent(in) :: j, p, tz
    type(DVec), intent(out) :: energies
    type(DMat), intent(out) :: states, UT
    type(OperatorDef) :: def_trel
    type(TwoBodyRelSpaceHOBasis) :: rel
    type(TwoBodyRelSpaceSpinHOBasis) :: relspin
    type(NNForceHO) :: vnnspin, Uspin
    type(TwoBodyRelOp) :: vnn, Trel, UnitaryTrans
    type(InputParameters) :: input
    type(DMat) :: h
    type(EigenSolSymD) :: sol
    real(8) :: ti
    type(sys) :: s
    integer :: ich

    input = params
    call def_trel%InitOpDef(s%str('kinetic'), .true.)

    call relspin%init(params%hw, params%N2max, params%J2max_NNint)
    call vnnspin%init(relspin)
    call Uspin%init(relspin)
    call vnnspin%setNNForceHO(Uspin, input)

    call rel%init(params%hw, params%N2max, params%Jmax2)
    call vnn%init(rel, rel, s%str("NNint"), params%pn_same_mass)
    call UnitaryTrans%init(rel, rel, s%str("UT"), params%pn_same_mass)
    call vnn%SetTwoBodyScalarSpinBreaking(vnnspin)
    call UnitaryTrans%SetTwoBodyScalarSpinBreaking(Uspin)

    if(params%spin_tensor_decomposition /= -1) then
      select case(params%spin_tensor_decomposition)
      case(0, 1, 2)
        vnn = vnn%SpinTensorDecomposition(params%spin_tensor_decomposition)
      case(3) ! C + LS
        vnn = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(1)
      case(4) ! C + T
        vnn = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(2)
      case(-100) ! Total = C + LS + T
        vnn = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(1) + vnn%SpinTensorDecomposition(2)
      case default
        write(*,*) 'Unknown value of spin_tensor_decomposition', __LINE__, __FILE__
      end select
    end if

    call vnnspin%fin()
    call Uspin%fin()
    call relspin%fin()
    call Trel%init(rel, rel, s%str("kinetic"), params%pn_same_mass)
    ich = rel%GetIndex(j, p, tz)
    call Trel%set()
    h = Trel%MatCh(ich,ich)%DMat + vnn%MatCh(ich,ich)%DMat

    ti = omp_get_wtime()
    call sol%init(h)
    call sol%DiagSym(h)
    call timer%add(s%str('solve_two_body_hamil'), omp_get_wtime() - ti)

    call states%ini(size(h%m,1),size(h%m,2))
    call energies%ini(size(h%m,1))
    states = sol%vec
    energies = sol%eig
    UT = UnitaryTrans%MatCh(ich,ich)%DMat

    call sol%fin()
    call h%fin()
    call Trel%fin()
    call vnn%fin()
    call rel%fin()
    call def_trel%FinOperatorDef()
  end subroutine solve_two_body_hamil

  subroutine calculation_electric_dipole_moment(params, egs, gs)
    use LinAlgLib
    use TwoBodyRelOps
    use MyLibrary, only: geometry_part, pi
    type(InputParameters), intent(in) :: params
    integer :: j, p, tz
    real(8), intent(in) :: egs, gs(:)
    type(DVec) :: inter_energies, groundstate
    type(DMat) :: inter_states, propagator, UT
    type(TwoBodyRelChanHOBasis), target :: rel, rel_flip
    type(TwoBodyRelOpChan) :: pvtv, e1
    integer :: i
    real(8) :: me
    type(sys) :: s
    real(8) :: edm_g1_lo, edm_c3, edm_c4, edm_delta_nlo, edm_delta_n2lo

    edm_g1_lo = edm_lec(s%str("PVTVint_ChEFT_gpi1-LO-Local2-500"))
    edm_delta_nlo = edm_lec(s%str("PVTVint_ChEFT_delta-NLO-Local2-500"))
    edm_c3 = edm_lec(s%str("PVTVint_ChEFT_Ct3-N2LO-Local2-500"))
    edm_c4 = edm_lec(s%str("PVTVint_ChEFT_Ct4-N2LO-Local2-500"))
    edm_delta_n2lo = edm_lec(s%str("PVTVint_ChEFT_delta-N2LO-Local2-500"))

    write(*, '(a)') '# Egs, g1 LO, D NLO, C3 NLO, C4 NLO, D N2LO'
    write(*, '(7f18.8)') egs, edm_g1_lo, edm_delta_nlo, edm_c3, edm_c4, edm_delta_n2lo
  contains
    function edm_lec(opname) result(edm)
      type(str), intent(in) :: opname
      real(8) :: edm
      call groundstate%ini(size(gs))
      groundstate%v = gs
      j = params%jket
      p = params%pket
      tz= params%tzket
      call rel%init(params%N2max,j,p,tz,params%hw)
      call rel_flip%init(params%N2max,j,-p,tz,params%hw)

      call solve_two_body_hamil(params,j,-p,tz,inter_energies,inter_states, UT)
      call propagator%zeros(rel_flip%GetNumberStates(), rel_flip%GetNumberStates() )
      do i = 1, rel_flip%GetNumberStates()
        propagator%m(i,i) = 1.d0 / (egs - inter_energies%v(i))
      end do
      propagator = inter_states * propagator * inter_states%T()

      call pvtv%init(opname, rel, rel_flip, params%pn_same_mass)
      call e1%init(s%str("E1"),rel_flip, rel, params%pn_same_mass)

      call pvtv%set( [j,p,tz], [j,-p,tz], params%NNint, params%pmax2, params%NMesh2)
      call e1%set([j,-p,tz],[j,p,tz])
      call inter_energies%fin()
      call inter_states%fin()
      me = groundstate * pvtv%DMat * propagator * e1%DMat * groundstate
      me = me * geometry_part(2*params%jket, 2, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from D operator
      me = me * geometry_part(2*params%jket, 0, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from PVTV operator
      me = me * 2.d0 ! c.c. term
      me = me / dble(params%particle_rank)
      edm = me * sqrt(4.d0 * pi / 3.d0)

      call e1%fin()
      call pvtv%fin()
      call rel_flip%fin()
      call rel%fin()
      call propagator%fin()
      call groundstate%fin()
    end function edm_lec

  end subroutine calculation_electric_dipole_moment

  subroutine calculation_electric_dipole_moment_isospin(params)
    use LinAlgLib
    use TwoBodyRelativeSpace
    use TwoBodyRelOpsIso
    use NNForceIsospin
    use MyLibrary, only: geometry_part, pi
    type(InputParameters), intent(in) :: params
    type(InputParameters) :: input
    integer :: j, p, t, z, ch, ch_flip, tt, i
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: vnn_rel, kin_rel
    type(TwoBodyRelOpIso) :: E1, pvtv_g1, pvtv_delta, pvtv_c3, pvtv_c4
    type(DVec) :: gs
    type(DMat) :: h, propagator
    type(EigenSolSymD) :: sol
    real(8) :: me, egs, edm, edm_g1, edm_delta, edm_c3, edm_c4
    type(sys) :: s

    j = params%jket
    p = params%pket
    t = params%tket
    z = params%tzket

    call relspin%init(params%hw, params%N2max, params%J2max_NNint)
    call rel%init(params%hw, params%N2max, params%Jmax2)

    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%renorm = "bare"
    call vnn_relspin%SetNNForceHOIsospin(U_relspin, input, 1, 0)
    call vnn_rel%init(rel,rel,s%str('NNint'))
    call vnn_rel%SetTwoBodyScalarSpinBreaking(vnn_relspin)

    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    call kin_rel%init(rel,rel,s%str('kinetic'))
    call kin_rel%set()

    call E1%init(rel,rel,s%str("E1"))
    call pvtv_g1%init(   rel,rel,s%str("PVTVint_ChEFT_g1-N2LO"))
    call pvtv_delta%init(rel,rel,s%str("PVTVint_ChEFT_delta-N2LO"))
    call pvtv_c3%init(   rel,rel,s%str("PVTVint_ChEFT_Ct3-N2LO"))
    call pvtv_c4%init(   rel,rel,s%str("PVTVint_ChEFT_Ct4-N2LO"))

    call E1%set()
    call pvtv_g1%set()
    call pvtv_delta%set()
    call pvtv_c3%set()
    call pvtv_c4%set()

    ch = rel%GetIndex(j,p,t)
    h = kin_rel%MatCh(ch,ch)%DMat + vnn_rel%MatCh(ch,ch)%DMat
    call sol%init(h)
    call sol%DiagSym(h)
    gs = sol%eig
    gs%v(:) = sol%vec%m(:,1)
    egs = sol%eig%v(1)
    edm = 0.d0
    tt = 1
    ch_flip = rel%GetIndex(j,-p,tt)
    h = kin_rel%MatCh(ch_flip,ch_flip)%DMat + vnn_rel%MatCh(ch_flip,ch_flip)%DMat
    call sol%init(h)
    call sol%DiagSym(h)
    call propagator%zeros( h%n_col, h%n_col )
    do i = 1, h%n_col
      propagator%m(i,i) = 1.d0 / (egs - sol%eig%v(i))
    end do
    propagator = sol%vec * propagator * sol%vec%T()

    me = gs * e1%MatCh(ch,ch_flip)%DMat * propagator * pvtv_g1%MatCh(ch_flip,ch)%DMat * gs
    me = me*geometry_part(2*tt, 2, 2*t,-2*z, 0,-2*z) ! from PVTV interaction
    me = me*geometry_part(2*t, 2, 2*tt, 2*z, 0, 2*z) ! from D operator
    edm = edm + me
    edm = edm*geometry_part(2*j, 2, 2*j, 2*j, 0, 2*j) ! from D operator
    edm = edm*geometry_part(2*j, 0, 2*j, 2*j, 0, 2*j) ! from PVTV operator
    edm = edm * 2.d0 / dble(params%particle_rank)
    edm_g1 = edm

    me = gs * e1%MatCh(ch,ch_flip)%DMat * propagator * pvtv_delta%MatCh(ch_flip,ch)%DMat * gs
    me = me*geometry_part(2*tt, 2, 2*t,-2*z, 0,-2*z) ! from PVTV interaction
    me = me*geometry_part(2*t, 2, 2*tt, 2*z, 0, 2*z) ! from D operator
    edm = edm + me
    edm = edm*geometry_part(2*j, 2, 2*j, 2*j, 0, 2*j) ! from D operator
    edm = edm*geometry_part(2*j, 0, 2*j, 2*j, 0, 2*j) ! from PVTV operator
    edm = edm * 2.d0 / dble(params%particle_rank)
    edm_delta = edm * sqrt(4.d0 * pi / 3.d0 )

    me = gs * e1%MatCh(ch,ch_flip)%DMat * propagator * pvtv_c3%MatCh(ch_flip,ch)%DMat * gs
    me = me*geometry_part(2*tt, 2, 2*t,-2*z, 0,-2*z) ! from PVTV interaction
    me = me*geometry_part(2*t, 2, 2*tt, 2*z, 0, 2*z) ! from D operator
    edm = edm + me
    edm = edm*geometry_part(2*j, 2, 2*j, 2*j, 0, 2*j) ! from D operator
    edm = edm*geometry_part(2*j, 0, 2*j, 2*j, 0, 2*j) ! from PVTV operator
    edm = edm * 2.d0 / dble(params%particle_rank)
    edm_c3 = edm * sqrt(4.d0 * pi / 3.d0 )

    me = gs * e1%MatCh(ch,ch_flip)%DMat * propagator * pvtv_c4%MatCh(ch_flip,ch)%DMat * gs
    me = me*geometry_part(2*tt, 2, 2*t,-2*z, 0,-2*z) ! from PVTV interaction
    me = me*geometry_part(2*t, 2, 2*tt, 2*z, 0, 2*z) ! from D operator
    edm = edm + me
    edm = edm*geometry_part(2*j, 2, 2*j, 2*j, 0, 2*j) ! from D operator
    edm = edm*geometry_part(2*j, 0, 2*j, 2*j, 0, 2*j) ! from PVTV operator
    edm = edm * 2.d0 / dble(params%particle_rank)
    edm_c4 = edm * sqrt(4.d0 * pi / 3.d0)

    write(*, '(5f18.8)') egs, edm_g1, edm_delta, edm_c3, edm_c4
    call vnn_rel%fin()
    call E1%fin()
    call PVTV_g1%fin()
    call PVTV_delta%fin()
    call PVTV_c3%fin()
    call PVTV_c4%fin()
    call rel%fin()
  end subroutine calculation_electric_dipole_moment_isospin

  subroutine calculation_toroidal_dipole_moment(params, egs, gs)
    use LinAlgLib
    use TwoBodyRelOps
    use MyLibrary, only: geometry_part
    use OperatorDefinitions
    type(InputParameters), intent(in) :: params
    integer :: j, p, tz
    real(8), intent(in) :: egs, gs(:)
    type(DVec) :: inter_energies, groundstate
    type(DMat) :: inter_states, propagator, UT
    type(TwoBodyRelChanHOBasis) :: rel, rel_flip
    type(TwoBodyRelOpChan) :: pvtc, tdm_mag, tdm_conv
    integer :: i
    real(8) :: me, me_conv
    type(sys) :: s

    call groundstate%ini(size(gs))
    groundstate%v = gs

    j = params%jket
    p = params%pket
    tz= params%tzket
    call rel%init(params%N2max,j,p,tz,params%hw)
    call rel_flip%init(params%N2max,j,-p,tz,params%hw)

    call solve_two_body_hamil(params,j,-p,tz,inter_energies,inter_states, UT)
    call propagator%zeros(rel_flip%GetNumberStates(), rel_flip%GetNumberStates())
    do i = 1, rel_flip%GetNumberStates()
      propagator%m(i,i) = 1.d0 / (egs - inter_energies%v(i))
    end do
    propagator = inter_states * propagator * inter_states%T()
    call inter_energies%fin()
    call inter_states%fin()

    call pvtc%init(s%str("PVTCint_DDH_hpi1"),rel, rel_flip, params%pn_same_mass)
    call tdm_mag%init(s%str("TDM_mag"),rel_flip, rel, params%pn_same_mass)
    call tdm_conv%init(s%str("TDM_conv"),rel_flip, rel, params%pn_same_mass)
    call pvtc%set( [j,p,tz], [j,-p,tz], params%NNint, params%pmax2, params%NMesh2)
    call tdm_mag%set([j,-p,tz],[j,p,tz])
    call tdm_conv%set([j,-p,tz],[j,p,tz])
    call tdm_mag%prt()
    tdm_mag%DMat = tdm_mag%DMat * (1.d0/ dble(params%particle_rank))
    tdm_conv%DMat = tdm_conv%DMat * ((params%particle_rank-1)**2/ dble(params%particle_rank)**3)

    me = groundstate * pvtc%DMat * propagator * tdm_mag%DMat * groundstate
    me = me * geometry_part(2*params%jket, 2, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from TDM operator
    me = me * geometry_part(2*params%jket, 0, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from PVTC operator
    me = me * 2.d0 ! c.c. term

    me_conv = groundstate * pvtc%DMat * propagator * tdm_conv%DMat * groundstate
    me_conv = me_conv * geometry_part(2*params%jket, 2, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from TDM operator
    me_conv = me_conv * geometry_part(2*params%jket, 0, 2*params%jket, 2*params%jket, 0, 2*params%jket) ! from PVTC operator
    me_conv = me_conv * 2.d0 ! c.c. term

    write(*,'(2a)') '# Nmax, hw, Energy (MeV), tdm mag, tdm conv'
    write(*,'(f18.8,5es14.6)') egs, me, me_conv

    call tdm_mag%fin()
    call rel_flip%fin()
    call rel%fin()
    call propagator%fin()
    call groundstate%fin()
  end subroutine calculation_toroidal_dipole_moment

  subroutine calculation_dipole_polarizability(params, egs, gs)
    use LinAlgLib
    use TwoBodyRelOps
    use MyLibrary, only: geometry_part, alpha, hc, pi
    use OperatorDefinitions
    type(InputParameters), intent(in) :: params
    integer :: j, p, tz, j_partner
    real(8), intent(in) :: egs, gs(:)
    type(DVec) :: inter_energies, groundstate
    type(DMat) :: inter_states, propagator, UT
    type(TwoBodyRelChanHOBasis) :: rel, rel_flip
    type(TwoBodyRelOpChan) :: e1
    integer :: i
    real(8) :: me, alphaD
    type(sys) :: s

    call groundstate%ini(size(gs))
    groundstate%v = gs

    j = params%jket
    p = params%pket
    tz= params%tzket
    call rel%init(params%N2max,j,p,tz,params%hw)
    alphaD = 0.d0
    do j_partner = abs( j-1 ), j+1
      call rel_flip%init(params%N2max,j_partner,-p,tz,params%hw)
      call solve_two_body_hamil(params,j_partner,-p,tz,inter_energies,inter_states, UT)
      call propagator%zeros(rel_flip%GetNumberStates(), rel_flip%GetNumberStates())
      do i = 1, rel_flip%GetNumberStates()
        propagator%m(i,i) = 1.d0 / (inter_energies%v(i)-egs)
      end do
      propagator = inter_states * propagator * inter_states%T()
      call inter_energies%fin()
      call inter_states%fin()
      call e1%init(s%str("E1"),rel, rel_flip, params%pn_same_mass)
      call e1%set( [j,p,tz], [j_partner,-p,tz] )
      me = groundstate * e1%DMat * propagator * e1%DMat%T() * groundstate
      me = me * (-1.d0)**( j_partner - j )
      me = me * geometry_part(2*params%jket, 2, 2*j_partner, 2*params%jket, 0, 2*params%jket)
      me = me * geometry_part(2*j_partner, 2, 2*params%jket, 2*params%jket, 0, 2*params%jket)
      alphaD = alphaD + me
      call e1%fin()
      call rel_flip%fin()
      call propagator%fin()
    end do
    alphaD = alphaD * 2.d0 * hc / (alpha * dble(params%particle_rank)**2) * 4.d0 * pi / 3.d0

    write(*, '(2a)') '# Energy (MeV), alphaD (fm3)'
    write(*, '(f18.8,5es14.6)') egs, alphaD

    call rel%fin()
    call groundstate%fin()
  end subroutine calculation_dipole_polarizability

  subroutine file_convert_nn(params)
    use SingleParticleState
    use TwoBodyLabOps
    type(InputParameters), intent(in) :: params
    type(Orbits) :: sps, sps_convert
    type(TwoBodyLabPNSpace) :: lab, lab_convert
    type(TwoBodyLabOp) :: oplab, oplab_convert
    type(str) :: filename

    if(params%file_name_nn_original%val == "none") then
      write(*,*) "Error, set the original file name."
      return
    end if

    write(*,*) ""
    write(*,"(a)") "#################### NN file convert mode ####################"

    call sps%init(params%emax, params%lmax)
    call lab%init(params%hw, sps, params%e2max, params%snt_mass)
    call oplab%init( lab, params%Operators(1) )

    filename = trim(params%file_name_nn_converted%val)

    write(*,*) "Reading..."
    call oplab%readf(params%file_name_nn_original)
    write(*,*) "Done."

    call sps_convert%init(params%emax_convert, params%lmax_convert)
    call lab_convert%init(params%hw, sps_convert, params%e2max_convert, params%snt_mass)
    oplab_convert = oplab%Truncate( lab_convert )
    call oplab%fin()
    call lab%fin()
    call sps%fin()

    write(*,*) "Writing..."
    call oplab_convert%writef(params%file_name_nn_converted)
    call oplab_convert%fin()
    call lab_convert%fin()
    call sps_convert%fin()

  end subroutine file_convert_nn

  subroutine file_combine_nn(params)
    use SingleParticleState
    use TwoBodyLabOps
    type(InputParameters), intent(in) :: params
    type(Orbits) :: sps
    type(TwoBodyLabPNSpace) :: lab
    type(TwoBodyLabOp) :: oplab, oplab_tmp
    integer :: i

    write(*,*) ""
    write(*,"(a)") "#################### NN file combine mode ####################"

    call sps%init(params%emax, params%lmax)
    call lab%init(params%hw, sps, params%e2max, params%snt_mass)

    write(*,*) "Reading ", trim(params%files_combined(1)%val)
    call oplab_tmp%init( lab, params%Operators(1) )
    call oplab_tmp%readf(params%files_combined(1))
    oplab = oplab_tmp * params%weights_combined(1)
    call oplab_tmp%fin()
    do i = 2, size(params%files_combined)
      write(*,*) "Reading ", trim(params%files_combined(i)%val)
      call oplab_tmp%init( lab, params%Operators(1) )
      call oplab_tmp%readf(params%files_combined(i))
      oplab = oplab + oplab_tmp * params%weights_combined(i)
      call oplab_tmp%fin()
    end do


    write(*,*) "Writing..."
    call oplab%writef(params%file_name_nn)
    call oplab%fin()
    call lab%fin()
    call sps%fin()
  end subroutine file_combine_nn

  subroutine Vrel_Cartesian( params )
    use MyLibrary, only: sys, pi
    use NNForceCartesian
    type(InputParameters), intent(in) :: params
    type(NNCartesian) :: vrel
    type(sys) :: s
    type(str), allocatable :: strings(:)
    type(str) :: tmp
    integer :: Nmax
    real(8) :: Lbox

    call s%split( params%file_nn_cartesian, s%str("Lbox"), strings)
    if( size(strings) /= 2 ) then
      write(*,*) "Error file_nn_cartesian has to be like: xxx_Lbox@@_Nmax@@_xxx"
      stop
    end if
    tmp = strings(2)
    call s%split( tmp, s%str("_"), strings)
    tmp = strings(1)
    read(tmp%val,*) Lbox

    call s%split( params%file_nn_cartesian, s%str("Nmax"), strings)
    if( size(strings) /= 2 ) then
      write(*,*) "Error file_nn_cartesian has to be like: xxx_Lbox@@_Nmax@@_xxx"
      stop
    end if
    tmp = strings(2)
    call s%split( tmp, s%str("_"), strings)
    tmp = strings(1)
    read(tmp%val,*) Nmax
    write(*,'(a,f7.3,a,i3)') "Generating the NN interaction file on Cartesian momentum grids: Lbox=", Lbox, " fm, Nmax=", Nmax
    call vrel%init( Lbox=Lbox, Nmax=Nmax )
    call vrel%SetNNCartesian( params )
    call vrel%WriteFile(filename=params%file_nn_cartesian%val)
    call vrel%fin()
  end subroutine Vrel_Cartesian

  function get_operator_in_lab_frame(params, lab) result(op)
    use SingleParticleState
    use TwoBodyLabSpacePN
    use TwoBodyTransCoef, only: TransRel2LabSpace
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceHO
    use TwoBodyRelOps
    use TwoBodyLabOps
    type(TwoBodyLabOp) :: op
    type(InputParameters), intent(in) :: params
    type(TwoBodyLabPNSpace), intent(in), target :: lab
    type(str) :: OpName
    type(Orbits), pointer :: sps
    type(TransRel2LabSpace) :: rel2lab
    type(TwoBodyRelSpaceHOBasis) :: rel, rel_for_trans
    type(TwoBodyRelSpaceSpinHOBasis) :: relspin
    type(NNForceHO) :: vnnspin, Uspin
    type(TwoBodyRelOp) :: vnn, U, oprel
    type(InputParameters) :: input
    type(sys) :: s

    OpName = params%Operators(1)
    sps => lab%sps
    call rel2lab%init(lab, params%e2max, params%Jmax2)
    call op%init(lab,OpName)

    ! --- calc vNN and UNN
    call relspin%init(params%hw, params%N2max, params%J2max_NNint)
    call vnnspin%init(relspin)
    call Uspin%init(relspin)
    input = params; input%N_ls = input%e2max
    call vnnspin%setNNForceHO(Uspin, input)
    call rel%init(params%hw, params%N2max, params%Jmax2)
    call vnn%init(rel,rel,s%str('hamil'), params%pn_same_mass)
    call U%init(rel,rel,s%str('UT'), params%pn_same_mass)
    call vnn%SetTwoBodyScalarSpinBreaking(vnnspin)
    call U%SetTwoBodyScalarSpinBreaking(Uspin)
    call vnnspin%fin()
    call Uspin%fin()
    call relspin%fin()
    ! ---

    if(OpName%val == 'hamil' .or. OpName%val == 'Hamil') then
      oprel = vnn
    else
      call oprel%init(rel, rel, OpName, params%pn_same_mass)
      if(s%find(OpName,s%str("PVTC")) .or. &
          & s%find(OpName,s%str("PVTV")) .or. &
          & s%find(OpName,s%str("AxialV")) .or. &
          & s%find(OpName,s%str("0vbb"))) then
        call oprel%set(params%NNint, pmax=params%pmax2, NMesh=params%NMesh2)
      else
        call oprel%set()
      end if
      call oprel%UT(U,U)
    end if
    call vnn%fin()
    call U%fin()

    call rel_for_trans%init(params%hw, params%e2max, params%Jmax2)
    oprel = oprel%truncate(rel_for_trans, rel_for_trans)
    call rel%fin()
    call op%TMtrans(rel2lab,oprel)

    call oprel%fin()
    call rel_for_trans%fin()
    call rel2lab%fin()
  end function get_operator_in_lab_frame

  subroutine nn_srg_analysis(params)
    use TwoBodyRelativeSpace
    use NNForce, only: NNForceHO
    type(InputParameters), intent(in) :: params
    type(TwoBodyRelSpaceSpinHOBasis) :: rel
    type(NNForceHO) :: vnn, U
    call rel%init(params%hw, params%N2max, params%J2max_NNint)
    call vnn%init(rel)
    call U%init(rel)

    call vnn%setNNForceHO(U, params, fn=params%output_nn_file)

    call vnn%fin()
    call U%fin()
    call rel%fin()
  end subroutine nn_srg_analysis

end module TwoBodyManager
