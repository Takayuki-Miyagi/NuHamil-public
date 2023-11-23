module ThreeBodyManager
  use LinAlgLib
  use Profiler, only: timer
  use NuHamilInput, only: InputParameters
  use ClassSys
  use MPIFunction, only: myrank
  use omp_lib
  implicit none
contains
  subroutine ThreeBodyManagers()
    use NuHamilInput, only: params
    use TMTransformManager
    use ThreeBodyLabFile
    type(sys) :: s
    type(InputParameters) :: input

    input = params
    if(myrank == 0 .and. params%particle_rank == 3) then
      write(*,'(a)') '############################################################'
      write(*,'(a)') '           Calculations for NNN sector'
      write(*,'(a)') '############################################################'
      write(*,*)
    end if

    if(params%pn_form) then
      if(myrank == 0) then
        write(*,"(a)") "You are trying to calcule three-body system."
        write(*,"(a)") "Then, 'pn_form' is turned to false."
      end if
      input%pn_form = .false.
    end if

    if( s%find( input%averaged_file_for_test, s%str("Jacobi")) ) then
      call check_norms(input)
      return
    end if

    if(input%file_convert) then
      if( input%only_no2b_elements ) call file_convert_no2b(input)
      if(.not. input%only_no2b_elements .and. .not. input%only_hf_monopole) call file_convert_3n(input)
      return
    end if

    if(allocated(input%files_combined) .and. allocated(input%weights_combined) &
        &.and. size(input%files_combined)==size(input%weights_combined)) then
      if(input%lab_3bme_precision%val=="half") then
        call file_combineHalf(params, input%files_combined, input%file_name_3n, input%weights_combined)
      else if(input%lab_3bme_precision%val=="single") then
        call file_combineSingle(params, input%files_combined, input%file_name_3n, input%weights_combined)
      else if(input%lab_3bme_precision%val=="double") then
        call file_combineDouble(params, input%files_combined, input%file_name_3n, input%weights_combined)
      else
        write(*,*) "Unknown precision: ", trim(input%lab_3bme_precision%val)
      end if
      return
    end if

    if(input%trans2lab) then
      if(input%coul) then
        if(myrank == 0) then
          write(*,"(a)") "You are trying to calcule three-body matrix element in Laboratory frame."
          write(*,"(a)") "Then, 'coul' is turned to false."
        end if
        input%coul = .false.
      end if

      if(input%NN_only) then
        if(myrank == 0) then
          write(*,"(a)") "You are trying to calcule three-body matrix element in Laboratory frame."
          write(*,"(a)") "Then, 'NN_only' is turned to false."
        end if
        input%NN_only = .false.
      end if

      call manage_trans_to_lab(input)

    else

      if(input%hw /= input%hw_target) then
        if(myrank == 0) then
          write(*,'(a)') "You are trying to calculate three-body eigen value problem."
          write(*,'(a,f6.2,a)') "Then, hw is ", input%hw, " MeV."
        end if
        input%hw_target = input%hw
      end if
      call calc_three_body_system(input)

    end if
  end subroutine ThreeBodyManagers

  subroutine calc_NNN_jacobi_isospin_channel( params )
    use ClassSys, only: sys
    use ThreeBodyJacobiSpace
    use NNNForceHOIsospin
    type(InputParameters), intent(in) :: params
    type(ThreeBodyJacIsoChan) :: jac
    integer :: j, p, t
    type(sys) :: s
    type(ThreeBodyJacOpChanIso) :: U
    type(NNNForceIsospin) :: v3
    type(str) :: f, fv, fut
    integer :: wunit = 21, runit = 22, Nmax

    j = params%jket
    p = params%pket
    t = params%tket
    Nmax = GetRampNmax(j, params%ramp)
    f = jac%GetFileName(j,p,t,Nmax,params%path_to_tmp_dir)
    write(*,*) "Calculating NNN force and transformation matrix..."
    write(*,'(a, i4, a, i2, a, i2, a, i2, a, i2)') &
        &  'myrank=', myrank, ',  :Calculating for J = ', &
        &   j, '/2,   P = ', p, ',   T = ', t, &
        &   '/2, Nmax = ', Nmax

    if(.not. s%isfile(f)) then
      call jac%init(params%hw,j,p,t,Nmax,params%path_to_tmp_dir)
      ! -- write to file
      open(wunit,file=f%val,status='replace',form='unformatted',access='stream')
      call jac%writef(wunit)
      close(wunit)
      call jac%fin()
    end if

    ! -- read from file
    open(runit,file=f%val,status='old',form='unformatted',access='stream')
    call jac%readf(params%hw,runit,Nmax)
    close(runit)

    if(jac%GetNumberNAStates() < 1 .or. jac%GetNumberAStates() < 1) then
      call jac%fin()
      return
    end if

    if(params%NN_only) then
      call jac%fin()
      write(*,'(a)') 'To calculate the 3BME, turn NN_only to .false.'
      return
    end if

    if(s%isfile(fv) .and. s%isfile(fut)) then
      call jac%fin()
      return
    end if
    call v3%InitNNNForce(jac)
    call U%init(jac,jac,s%str("UT"))
    call v3%SetNNNForce(U, params)

    call v3%fin()
    call U%fin()
    call jac%fin()
  end subroutine calc_NNN_jacobi_isospin_channel

  subroutine calc_three_body_system(params)
    type(InputParameters), intent(in) :: params
    if(mod(params%jbra, 2) == 0) stop 'Warning: Jbra have to be odd number'
    if(mod(params%jket, 2) == 0) stop 'Warning: Jket have to be odd number'
    if(mod(params%tbra, 2) == 0) stop 'Warning: tbra have to be odd number'
    if(mod(params%tket, 2) == 0) stop 'Warning: tket have to be odd number'
    if(mod(params%tzbra, 2) == 0) stop 'Warning: tzbra have to be odd number'
    if(mod(params%tzket, 2) == 0) stop 'Warning: tzket have to be odd number'
    if(abs(params%tzbra) > params%tbra) stop 'Warning: tz_bra have to be smaller than t_bra'
    if(abs(params%tzket) > params%tket) stop 'Warning: Tz_ket have to be smaller than t_ket'

    call calc_three_body(params)
  end subroutine calc_three_body_system

  subroutine calc_three_body(params)
    type(InputParameters), intent(in) :: params
    call calc_three_body_isospin(params)
  end subroutine calc_three_body

  subroutine calc_three_body_isospin(params)
    use MyLibrary, only: triag, geometry_part, pi
    use ClassSys, only: sys
    use OperatorDefinitions
    type(InputParameters), intent(in) :: params
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    type(str) :: f, oprtr
    logical :: calc_bra = .true.
    type(DMat) :: bra_states, ket_states
    type(DVec) :: bra_energies, ket_energies
    integer :: iloop

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket

    if(jbra == jket .and. pbra == pket .and. tbra == tket .and. zbra==zket) calc_bra = .false.
    call solve_three_body_hamil_isospin(params,jket,pket,tket,zket,ket_energies,ket_states)
    if(calc_bra) then
      call solve_three_body_hamil_isospin(params,jbra,pbra,tbra,zbra,bra_energies,bra_states)
    else
      bra_energies = ket_energies
      bra_states = ket_states
    end if

    do iloop = 1, size(params%Operators)
      oprtr = params%Operators(iloop)
      call inside()
    end do

    call bra_states%fin()
    call ket_states%fin()
    call bra_energies%fin()
    call ket_energies%fin()
  contains
    subroutine inside()
      use MyLibrary
      type(sys) :: s
      integer :: i
      integer :: A, Z, N
      type(OperatorDef) :: opdef
      type(str) :: tmp
      real(8) :: val, vals(3), vals0(3), vals1(3), vals2(3), vals3(3), vals4(3)

      A=3
      Z=(3-zket)/2
      N=(3+zket)/2
      select case(oprtr%val)
      case('EDM',"TDM","Rm2","Rp2","Rn2","AlphaD","Mmoment_IA","Mmoment_2B","Mmoment_full","Qmoment_IA")
        if(jbra/=jket) return
        if(pbra/=pket) return
        if(tbra/=tket) return
        if(zbra/=zket) return
      case("BetaDecay")
        if(zbra==zket) return
      case default
        call opdef%InitOpDef(oprtr, .false.)
        if(triag(jbra,jket,opdef%GetOpJ())) return
        if(pbra * pket * opdef%GetOpP() /= 1) return
        if(triag(tbra,tket,opdef%GetOpT())) return
      end select

      select case(oprtr%val)
      case('hamil', 'Hamil')
        write(*,'(a18)') "Energy"
        write(*,'(f18.8)') ket_energies%v(1)

      case("EDM")
        if(params%renorm%val /= "bare") then
          write(*,*) "Not implemented yet, use renorm='bare'."
          return
        end if
        call calc_edm_isospin(params, ket_energies%v(1), ket_states%m(:,1))

      case("AlphaD")
        if(params%renorm%val /= "bare") then
          write(*,*) "Not implemented yet, use renorm='bare'."
          return
        end if
        vals0= calc_dipole_polarizability_isospin(params, ket_energies%v(1), ket_states%m(:,1))
        write(*,'(2a18)') "Energy", "alphaD"
        write(*,'(2f18.8)') ket_energies%v(1), vals0

      case("TDM")
        if(params%renorm%val /= "bare") then
          write(*,*) "Not implemented yet, use renorm='bare'."
          return
        end if
        call calc_tdm_isospin(params, ket_energies%v(1), ket_states%m(:,1))

      case("Rm2")
        vals = calc_expectation_val_isospin(s%str("R2"), params, bra_states, ket_states)
        vals = vals / dble(A**2)
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,*)

      case("Rp2")
        vals0 = calc_expectation_val_isospin(s%str("R2"), params, bra_states, ket_states)
        vals1 = calc_expectation_val_isospin(s%str("R2_IV"), params, bra_states, ket_states)
        vals0 = vals0 * dble(N) / dble(Z*A**2)
        vals1 = vals1 * geometry_part(tbra,2,tket,zbra,0,zket)*geometry_part(jbra,0,jket,jbra,0,jket) / dble(2*A*Z)
        vals = vals0 - vals1
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,*)

      case("Rn2")
        vals0 = calc_expectation_val_isospin(s%str("R2"), params, bra_states, ket_states)
        vals1 = calc_expectation_val_isospin(s%str("R2_IV"), params, bra_states, ket_states)
        vals0 = vals0 * dble(Z) / dble(N*A**2)
        vals1 = vals1 * geometry_part(tbra,2,tket,zbra,0,zket)*geometry_part(jbra,0,jket,jbra,0,jket) / dble(2*A*N)
        vals = vals0 + vals1
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,*)

      case("Mmoment_IA")
        vals0 = calc_expectation_val_isospin(s%str("M1_S_IS"), params, bra_states, ket_states)
        vals1 = calc_expectation_val_isospin(s%str("M1_S_IV"), params, bra_states, ket_states)
        vals2 = calc_expectation_val_isospin(s%str("M1_L_IS"), params, bra_states, ket_states)
        vals3 = calc_expectation_val_isospin(s%str("M1_L_IV"), params, bra_states, ket_states)

        vals0 = vals0 * geometry_part(tbra,0,tket,zbra,0,zket)*geometry_part(jbra,2,jket,jbra,0,jket) / dble(A-1)
        vals1 = vals1 * geometry_part(tbra,2,tket,zbra,0,zket)*geometry_part(jbra,2,jket,jbra,0,jket) / dble(A-1)
        vals2 = vals2 * geometry_part(tbra,0,tket,zbra,0,zket)*geometry_part(jbra,2,jket,jbra,0,jket) * dble(2*N) / dble(A)**2
        vals3 = vals3 * geometry_part(tbra,2,tket,zbra,0,zket)*geometry_part(jbra,2,jket,jbra,0,jket) / dble(A)
        vals0 = vals0 * sqrt(4.d0*pi/3.d0)
        vals1 = vals1 * sqrt(4.d0*pi/3.d0)
        vals2 = vals2 * sqrt(4.d0*pi/3.d0)
        vals3 = vals3 * sqrt(4.d0*pi/3.d0) * (-1.d0)
        vals = vals0 + vals1 + vals2 + vals3
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,*)

      case("Qmoment_IA")
        vals0 = calc_expectation_val_isospin(s%str("E2_IS"), params, bra_states, ket_states)
        vals1 = calc_expectation_val_isospin(s%str("E2_IV"), params, bra_states, ket_states)
        vals0 = vals0 * geometry_part(tbra,0,tket,zbra,0,zket)*geometry_part(jbra,4,jket,jbra,0,jket) * dble(N) / dble(A)**2
        vals1 = vals1 * geometry_part(tbra,2,tket,zbra,0,zket)*geometry_part(jbra,4,jket,jbra,0,jket) / dble(2*A)
        vals = vals0 - vals1
        vals = vals * sqrt(16.d0 * pi/5.d0)
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,*)

      case("Mmoment_2B")
        vals = calc_expectation_val_isospin(s%str("M1_2BC_intr"), params, bra_states, ket_states)
        val = calc_expectation_val_cm_isospin(s%str("M1_2BC_Sachs"), params, bra_states, ket_states)
        vals = vals * geometry_part(tbra,2,tket,zbra,0,zket) * geometry_part(jbra,2,jket,jbra,0,jket)
        vals = vals * sqrt(4.d0 * pi/3.d0)
        val = val * geometry_part(tbra,2,tket,zbra,0,zket) * geometry_part(jbra,2,jket,jbra,0,jket)
        val = val * sqrt(4.d0 * pi/3.d0)
        write(*,'(a16,6a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>", &
            & "<Op (2b ev)>", "<Op (3b ev)>"
        write(*,'(a16,5f18.8)') "M1_2BC_intr",bra_energies%v(1), ket_energies%v(1), vals(:)
        write(*,'(a16,3a18)') "Operator", "Energy bra", "Energy ket", "<Op (bare)>"
        write(*,'(a16,3f18.8)') "M1_2BC_Sachs",bra_energies%v(1), ket_energies%v(1), val
        write(*,*)

      case("BetaDecay")
        if(s%find(params%Regulator, s%str("NonLocal"))) then
          tmp = "AxialV_T1-N2LO-" + params%Regulator + s%str(params%RegulatorPower) + "-" + s%str(params%lambda_3nf_nonlocal)
        else if(s%find(params%Regulator, s%str("Local"))) then
          tmp = "AxialV_T1-N2LO-" + params%Regulator + s%str(params%RegulatorPower) + "-" + s%str(params%lambda_3nf_local)
        else
          tmp = "AxialV_T1-N2LO"
        end if
        vals0 = calc_expectation_val_isospin(s%str("Fermi"), params, bra_states, ket_states)
        vals0 = vals0 * geometry_part(tbra,2,tket,zbra,zbra-zket,zket) / dble(A-1)
        vals1 = calc_expectation_val_isospin(s%str("GamowTeller"), params, bra_states, ket_states)
        vals1 = vals1 * geometry_part(tbra,2,tket,zbra,zbra-zket,zket) / dble(A-1)
        vals2 = calc_expectation_val_isospin(tmp, params, bra_states, ket_states)
        vals2 = vals2 * geometry_part(tbra,2,tket,zbra,zbra-zket,zket)  / g_A * (-1.d0) * sqrt(2.d0)
        write(*,'(a16,5f18.8)') "<|| Fermi ||>: ",bra_energies%v(1), ket_energies%v(1), vals0(:)
        write(*,'(a16,5f18.8)') "<|| GT ||>: ",bra_energies%v(1), ket_energies%v(1), vals1(:)
        write(*,'(a16,5f18.8)') "<|| A(2BC) ||>: ",bra_energies%v(1), ket_energies%v(1), vals2(:)
        write(*,*)


      case default
        if(s%find(oprtr, s%str("Tmag_1B"))) then
          val = calc_expectation_val_onebody_isospin(oprtr, params, bra_states, ket_states, 0) * &
            & geometry_part(tbra,0,tket,zbra,0,zket)
          val = val + calc_expectation_val_onebody_isospin(oprtr, params, bra_states, ket_states, 1) * &
            & geometry_part(tbra,2,tket,zbra,0,zket)
          val = val / opdef%GetQ() * 2.d0 * m_nucleon * (-1.d0) * sqrt(pi)
          val = val * dble(opdef%GetOpP()) ! this comes from the embedding phase (-1)^l3'+l3
          write(*,'(a16,3f18.8,es18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), opdef%GetQ(), val

        else if(s%find(oprtr, s%str("Tmag_2B"))) then
          val = calc_expectation_val_cm_isospin(oprtr, params, bra_states, ket_states) * &
            & geometry_part(tbra,2,tket,zbra,0,zket)
          val = val / opdef%GetQ() * 2.d0 * m_nucleon * (-1.d0) * sqrt(pi)
          val = val * dble(opdef%GetOpP()) ! this comes from the embedding phase (-1)^l3'+l3
          write(*,'(a16,3f18.8,es18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), opdef%GetQ(), val

        else if(s%find(oprtr, s%str("M_1B"))) then
          val = calc_expectation_val_onebody_isospin(oprtr, params, bra_states, ket_states, 0) * &
            & geometry_part(tbra,0,tket,zbra,0,zket)
          val = val + calc_expectation_val_onebody_isospin(oprtr, params, bra_states, ket_states, 1) * &
            & geometry_part(tbra,2,tket,zbra,0,zket)
          val = val * dble(opdef%GetOpP()) ! this comes from the embedding phase (-1)^l3'+l3
          val = val * sqrt(2.d0 * pi) / dble(Z)
          write(*,'(a16,3f18.8,es18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), opdef%GetQ(), val

        else
          vals = calc_expectation_val_isospin(oprtr, params, bra_states, ket_states)
          write(*,'(a16,5f18.8)') oprtr%val,bra_energies%v(1), ket_energies%v(1), vals(:)
        end if
      end select
    end subroutine inside

  end subroutine calc_three_body_isospin

  function calc_expectation_val_isospin(opname, params, bra_vecs, ket_vecs) result(vals)
    use ClassSys, only: sys
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    type(str), intent(in) :: opname
    type(InputParameters), intent(in) :: params
    type(DMat), intent(in) :: bra_vecs, ket_vecs
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: U, op, opeff
    type(ThreeBodyJacIsoChan) :: chbra, chket
    type(ThreeBodyJacOpChanIso) :: opjac_bare, opjac_2b, opjac_3b, ut
    type(InputParameters) :: input
    type(sys) :: s
    type(str) :: futbra, futket, fn
    integer :: iutbra = 101, iutket = 102, wunit = 103
    real(8) :: exp_bare, exp_2b, exp_3b
    real(8) :: vals(3)
    logical :: ex, calc_bra = .true.
    integer :: A, Z, N

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket
    if(jbra == jket .and. pbra == pket .and. tbra == tket .and. zbra==zket) calc_bra = .false.

    call chbra%init(params%hw,jbra,pbra,tbra,params%N3max,params%path_to_tmp_dir)
    call chket%init(params%hw,jket,pket,tket,params%N3max,params%path_to_tmp_dir)

    call relspin%init(params%hw,params%N3max, params%J2max_NNint)
    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = params%N3max
    call vnn_relspin%setNNForceHOIsospin(U_relspin, input, params%Tket, params%Tzket)

    call rel%init(params%hw, params%N3max, params%Jmax2)
    call U%init(rel,rel,s%str("UT"))
    call U%SetTwoBodyScalarSpinBreaking(U_relspin)
    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    call op%init(rel,rel,opname)
    A = 3; Z = (3-zket)/2; N = (3+zket)/2
    call op%set()

    call opeff%init(rel,rel,opname)
    call opeff%set()
    call opeff%UT(U,U)
    call U%fin()

    ! --- bare ---
    call opjac_bare%init(chbra,chket,opname)
    call opjac_bare%set(op)

    ! --- 2b evolved ---
    call opjac_2b%init(chbra,chket,opname)
    call opjac_2b%set(opeff)

    ! --- 3b evolved ---
    call op%set()
    call opjac_3b%init(chbra,chket,opname)
    if(.not. params%NN_only .and. params%renorm%val /= "bare") then

      futbra = ut%GetFileName(chbra,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      futket = ut%GetFileName(chket,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      if(params%hw_conversion) then
        futbra = ut%GetFileName(chbra,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
        futket = ut%GetFileName(chket,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
            & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
            & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
            & params%path_to_tmp_dir)
      end if
      ex = s%isfile(futbra, s%str('In set_three_body_jac_operator_isospin'))
      ex = s%isfile(futket, s%str('In set_three_body_jac_operator_isospin'))
      if(.not. calc_bra) then
        open(iutket, file=futket%val, form='unformatted', access='stream')
        call opjac_3b%set(chbra, chket, &
            & iutket, iutket, op, opeff, &
            & params%N3max, params%N3max, params%hw, params%hw_target, params%hw_conversion)
        close(iutket)
      else

        open(iutbra, file=futbra%val, form='unformatted', access='stream')
        open(iutket, file=futket%val, form='unformatted', access='stream')
        call opjac_3b%set(chbra, chket, iutbra, iutket, &
            & op, opeff, &
            & params%N3max, params%N3max, params%hw, params%hw_target, params%hw_conversion)
        close(iutbra)
        close(iutket)
      end if
    end if

    exp_bare = dot_product(bra_vecs%m(:,1), matmul(opjac_bare%m, ket_vecs%m(:,1)))
    exp_2b = dot_product(bra_vecs%m(:,1), matmul(opjac_2b%m, ket_vecs%m(:,1)))
    exp_3b = dot_product(bra_vecs%m(:,1), matmul(opjac_3b%m, ket_vecs%m(:,1)))
    if(params%fn_fewbody_wave_function%val /= 'none') then
      fn = opname + "_" + params%fn_fewbody_wave_function + "_" + &
          & s%str(jbra) + s%str(pbra) + s%str(tbra) + s%str(zbra) + "_" + &
          & s%str(jket) + s%str(pket) + s%str(tket) + s%str(zket) + ".op"
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) opjac_bare%m(:,:)
      close(wunit)
    end if
    call opjac_3b%fin()
    call opjac_2b%fin()
    call opjac_bare%fin()
    call opeff%fin()
    call op%fin()
    call rel%fin()
    vals(:) = [exp_bare, exp_2b, exp_2b+exp_3b]
    call chket%fin()
    call chbra%fin()

  end function calc_expectation_val_isospin

  function calc_expectation_val_cm_isospin(opname, params, bra_vecs, ket_vecs) result(val)
    use MyLibrary, only: gauss_legendre
    use ClassSys, only: sys
    use TwoBodyRelativeSpace
    use TwoBodyRelCMIsoOps
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    use TwoBodyMultiPoleOperator
    use TwoBodyRelCMOps
    type(str), intent(in) :: opname
    type(InputParameters), intent(in) :: params
    type(DMat), intent(in) :: bra_vecs, ket_vecs
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    type(TwoBodyRelCMSpaceIsoHOBasis) :: relcm
    type(TwoBodyRelCMIsoOp) :: op
    type(ThreeBodyJacIsoChan) :: chbra, chket
    type(ThreeBodyJacOpChanIso) :: opjac_bare
    type(InputParameters) :: input
    type(sys) :: s
    type(str) :: futbra, futket
    integer :: iutbra = 101, iutket = 102
    real(8) :: exp_bare, exp_2b, exp_3b
    real(8) :: val
    logical :: ex, calc_bra = .true.
    integer :: A, Z, N
    type(TwoBodyRelCMSpaceHOBasis) :: relcm_pn
    type(TwoBodyRelCMSpaceMBasis) :: relcm_mom
    type(TwoBOdyRelCMOp) :: op_pn
    type(TwoBodyMultiPoleOp) :: mltop_2b
    real(8), allocatable :: x(:), w(:)

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket
    if(jbra == jket .and. pbra == pket .and. tbra == tket .and. zbra==zket) calc_bra = .false.

    call chbra%init(params%hw,jbra,pbra,tbra,params%N3max,params%path_to_tmp_dir)
    call chket%init(params%hw,jket,pket,tket,params%N3max,params%path_to_tmp_dir)

    call relcm%init(params%hw, params%N3max, params%Jmax2, params%N3max, params%Jmax2)
    call op%init(relcm,opname)


    if(s%find(opname,s%str("Tmag_2B"))) then
      call op%fin()
      call op%init(relcm, 1, 1, 1)
      call relcm_pn%init(params%hw, params%N3max, params%J2maxLab, params%N3max, params%jmax2, LcmMax=params%Lcm2Max)
      call op_pn%init(relcm_pn, opname)
      call relcm_mom%init(0.d0, params%pmax2, params%NMesh2, 0.d0, params%pmax2, params%NMesh2, &
        & params%J2maxLab, params%jmax2, Lcm_Max=params%Lcm2Max)
      call gauss_legendre(0.d0, params%pmax2, x, w, params%NMesh2)
      call relcm_mom%SetMeshWeight(x, w, x, w)
      call mltop_2b%init(op_pn, relcm_mom, c1=0.d0, c3=0.d0, c4=0.d0, c6=0.d0, cD=0.d0)
      call mltop_2b%SetHOMatrix(op_pn, s%str(''))
      call op%MakeIsospinSymmetric(op_pn)

      call mltop_2b%fin()
      call relcm_mom%fin()
      call op_pn%fin()
      call relcm_pn%fin()
    else

      call op%set()
    end if

    call opjac_bare%init(chbra,chket,opname)
    call opjac_bare%set(op)
    exp_bare = dot_product(bra_vecs%m(:,1), matmul(opjac_bare%m, ket_vecs%m(:,1)))
    call opjac_bare%fin()
    call op%fin()
    call relcm%fin()
    val = exp_bare
    call chket%fin()
    call chbra%fin()

  end function calc_expectation_val_cm_isospin

  function calc_expectation_val_onebody_isospin(opname, params, bra_vecs, ket_vecs, isospin) result(val)
    use ClassSys, only: sys
    use SingleParticleState
    use OneBodyLabOps
    use OneBodyLabOpsIso
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    type(str), intent(in) :: opname
    type(InputParameters), intent(in) :: params
    type(DMat), intent(in) :: bra_vecs, ket_vecs
    integer :: isospin
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    type(Orbits) :: sps_pn
    type(OrbitsIsospin) :: sps
    type(OneBodyLabOp) :: op_pn
    type(OneBodyLabOpIso) :: op
    type(ThreeBodyJacIsoChan) :: chbra, chket
    type(ThreeBodyJacOpChanIso) :: opjac
    type(InputParameters) :: input
    type(sys) :: s
    real(8) :: val
    logical :: ex, calc_bra = .true.
    integer :: A, Z, N

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket
    if(jbra == jket .and. pbra == pket .and. tbra == tket .and. zbra==zket) calc_bra = .false.

    call chbra%init(params%hw,jbra,pbra,tbra,params%N3max,params%path_to_tmp_dir)
    call chket%init(params%hw,jket,pket,tket,params%N3max,params%path_to_tmp_dir)

    call sps%init(params%N3max)
    call sps_pn%init(params%N3max)
    call op_pn%init(sps_pn, params%hw, opname)
    if(opname%val == "Sp_M1") then
      call op_pn%set(isospin, e_charge=2.d0/3.d0)
    else if(s%find(opname,s%str("M_1B")) .or. s%find(opname,s%str("Tmag_1B"))) then
      op_pn%hw = params%hw * 3.d0 / 2.d0
      call op_pn%set(isospin)
      op_pn%hw = params%hw
    else
      call op_pn%set(isospin)
    end if

    call op%init(sps, op_pn%GetOpJ(), op_pn%GetOpP(), isospin, params%hw)
    call op%MakeIsospinSymmetric(op_pn)

    call op_pn%fin()
    call sps_pn%fin()

    call opjac%init(chbra,chket,s%str(""))
    call opjac%set(op)

    val = dot_product(bra_vecs%m(:,1), matmul(opjac%m, ket_vecs%m(:,1)))
    call opjac%fin()
    call op%fin()
    call sps%fin()
    call chket%fin()
    call chbra%fin()

  end function calc_expectation_val_onebody_isospin

  subroutine solve_three_body_hamil_isospin(params,j,p,t,tz,energies,states)
    use LinAlgLib
    use ClassSys, only: sys
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    use NNNForceHOIsospin
    type(InputParameters), intent(in) :: params
    integer, intent(in) :: j, p, t, tz
    type(DVec), intent(out) :: energies
    type(DMat), intent(out) :: states
    type(DMat) :: h
    type(sys) :: s
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelOpIso) :: vnn, kin_rel, vnn_tmp
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(ThreeBodyJacIsoChan) :: jac
    type(ThreeBodyJacOpChanIso) :: kin_jac, vnn_jac, U_jac
    type(NNNForceIsospin) :: v3n_jac
    type(EigenSolSymD) :: sol
    type(InputParameters) :: input
    type(DVec) :: gs_wf
    type(str) :: fjac
    type(str) :: fn
    integer :: wjac = 100, rjac = 101, wunit=20
    real(8) :: ti

    call relspin%init(params%hw,params%N2max, params%J2max_NNint)

    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = params%N3max
    call vnn_relspin%setNNForceHOIsospin(U_relspin, input, t, tz)

    call rel%init(params%hw, params%N2max, params%Jmax2)
    call vnn%init(rel,rel,s%str('NNint'))
    call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)
    if(params%spin_tensor_decomposition /= -1) then
      select case(params%spin_tensor_decomposition)
      case(0, 1, 2)
        vnn_tmp = vnn%SpinTensorDecomposition(params%spin_tensor_decomposition)
        vnn = vnn_tmp
      case(3) ! C + LS
        vnn_tmp = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(1)
        vnn = vnn_tmp
      case(4) ! C + T
        vnn_tmp = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(2)
        vnn = vnn_tmp
      case(-100) ! C + LS + T
        vnn_tmp = vnn%SpinTensorDecomposition(0) + vnn%SpinTensorDecomposition(1) + vnn%SpinTensorDecomposition(2)
        vnn = vnn_tmp
      case default
        write(*,*) 'Unknown value of spin_tensor_decomposition', __LINE__, __FILE__
      end select
    end if
    call kin_rel%init(rel,rel,s%str('kinetic'))
    call kin_rel%set()
    fjac = jac%GetFileName(j,p,t,params%N3max,params%path_to_tmp_dir)
    if(.not. s%isfile(fjac)) then
      call jac%init(params%hw,j,p,t,params%N3max, params%path_to_tmp_dir)
      open(wjac,file=fjac%val,status='replace',form='unformatted',access='stream')
      call jac%writef(wjac)
      close(wjac)
      call jac%fin()
    end if
    open(rjac,file=fjac%val,status='old',form='unformatted',access='stream')
    call jac%readf(params%hw,rjac)
    close(rjac)

    call kin_jac%init(jac,jac,s%str('kinetic'))
    call vnn_jac%init(jac,jac,s%str('NNint'))
    call kin_jac%set(kin_rel)
    call vnn_jac%set(vnn)

    call U_jac%init(jac,jac,s%str("UT"))
    call v3n_jac%init(jac)
    if(.not. params%NN_only) then
      call v3n_jac%SetNNNForce(U_jac, params)
    end if

    ti = omp_get_wtime()

    h = kin_jac%DMat + vnn_jac%DMat + v3n_jac%DMat
    call sol%init(h)
    call sol%DiagSym(h)
    states = sol%vec
    energies = sol%eig

    if(params%fn_fewbody_wave_function%val /= 'none') then
      fn = "wf_" + params%fn_fewbody_wave_function + "_" + s%str(j) + s%str(p) + s%str(t) + s%str(tz) + ".wav"
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) states%m(:,1)
      close(wunit)
      fn = "Hamil_" + params%fn_fewbody_wave_function + "_" + s%str(j) + s%str(p) + s%str(t) + s%str(tz) + ".op"
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) h%m(:,:)
      close(wunit)
    end if

    call sol%fin()
    call h%fin()
    call v3n_jac%fin()
    call vnn_jac%fin()
    call kin_jac%fin()
    call kin_rel%fin()
    call vnn%fin()
    call rel%fin()
    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()
    call timer%add(s%str('solve_three_body_hamil_isospin'), omp_get_wtime() - ti)
  end subroutine solve_three_body_hamil_isospin

  subroutine calc_edm_isospin(params, egs, gs)
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    real(8) :: g0, g1, g2, delta, ct1, ct2, ct3, ct4, ct5
    type(sys) :: s
    g0    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_OBE_gpi0") )
    g1    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_OBE_gpi1") )
    g2    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_OBE_gpi2") )
    !g0    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_gpi0-N2LO-Local2-500") )
    !g1    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_gpi1-N2LO-Local2-500") )
    !g2    = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_gpi2-N2LO-Local2-500") )
    delta = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_delta-N2LO-Local2-500") )
    Ct1   = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_Ct1-N2LO-Local2-500") )
    Ct2   = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_Ct2-N2LO-Local2-500") )
    Ct3   = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_Ct3-N2LO-Local2-500") )
    Ct4   = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_Ct4-N2LO-Local2-500") )
    Ct5   = calc_edm_isospin_lec( params, egs, gs, s%str("PVTVint_ChEFT_Ct5-N2LO-Local2-500") )
    write(*,*) "egs, g0, g1, g2, delta, ct1, ct2, ct3, ct4, ct5"
    write(*,'(10f18.8)') egs, g0, g1, g2, delta, ct1, ct2, ct3, ct4, ct5
  end subroutine calc_edm_isospin

  function calc_edm_isospin_lec(params, egs, gs, opname) result(edm)
    use LinAlgLib
    use ClassSys, only: sys
    use MyLibrary, only: geometry_part, pi
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use TwoBodyRelOps
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    use NNNForceHOIsospin
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    type(str), intent(in) :: opname
    type(DVec) :: groundstate, inter_energies
    type(DMat) :: propagator, inter_states
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelSpaceHOBasis) :: rel_pn
    type(ThreeBodyJacIsoChan) :: jacobi_ch, jacobi_ch_flip
    type(TwoBodyRelOpIso) :: E1_2, PVTV_2
    type(ThreeBodyJacOpChanIso) :: E1_3, PVTV_3
    integer :: ttot_flip, i
    real(8) :: me, edm
    type(sys) :: s
    ! for details
    !real(8) :: me_e1, me_pv
    !type(DVec) :: tmp

    call groundstate%ini(size(gs))
    groundstate%v = gs
    call jacobi_ch%init(params%hw, params%jket, params%pket, params%tket, params%N3max, params%path_to_tmp_dir)
    call rel%init(params%hw, params%N3max, params%Jmax2)
    call rel_pn%init(params%hw, params%N3max, params%Jmax2)

    call E1_2%init(rel,rel,s%str("E1"))
    call PVTV_2%init(rel,rel,opname)
    ! set two-body operators
    call E1_2%set()
    call PVTV_2%set()

    edm = 0.d0
    do ttot_flip = 1, 3, 2
      call jacobi_ch_flip%init(params%hw,params%jket,-params%pket,ttot_flip,params%N3max,params%path_to_tmp_dir)
      call solve_three_body_hamil_isospin(params,params%jket,-params%pket,ttot_flip,params%tzket,inter_energies,inter_states)
      call propagator%zeros( jacobi_ch_flip%GetNumberAStates(), jacobi_ch_flip%GetNumberAStates() )
      do i = 1, jacobi_ch_flip%GetNumberAStates()
        propagator%m(i,i) = 1.d0 / (egs - inter_energies%v(i))
      end do
      propagator = inter_states * propagator * inter_states%t()

      call E1_3%init(jacobi_ch, jacobi_ch_flip, s%str("E1"))
      call E1_3%set(E1_2)
      call PVTV_3%init(jacobi_ch_flip, jacobi_ch, s%str("PVTVint"))
      call PVTV_3%set(PVTV_2)

      ! detail
      !call tmp%ini( jacobi_ch_flip%GetNumberAStates() )
      !do i = 1, jacobi_ch_flip%GetNumberAStates()
      !  tmp%v(:) = inter_states%m(:,i)
      !  me_e1 = groundstate * E1_3%DMat * tmp * &
      !      &   geometry_part(params%tket, 2*E1_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket)
      !  me_pv = tmp * PVTV_3%DMat * groundstate * &
      !      &   geometry_part(ttot_flip, 2*PVTV_2%GetOpT(), params%tket, params%tzket, 0, params%tzket)
      !  !call tmp%prt("inter")
      !  !call PVTV_3%prt("pv int")
      !  !call groundstate%prt("gs")
      !  write(*,'(a,3f12.6)') "energy denominator, E1, PVTV int.:", egs-inter_energies%v(i), me_e1, me_pv
      !end do
      !call tmp%fin()

      me = groundstate * E1_3%DMat * propagator * PVTV_3%DMat * groundstate
      me = me*geometry_part(ttot_flip, 2*PVTV_2%GetOpT(), params%tket, params%tzket, 0, params%tzket) ! from PVTV interaction
      me = me*geometry_part(params%tket, 2*E1_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket) ! from D operator
      edm = edm + me

      call E1_3%fin(); call PVTV_3%fin()
      call propagator%fin()
      call inter_energies%fin()
      call inter_states%fin()
      call jacobi_ch_flip%fin()
    end do
    edm = edm * 2.d0 / dble(params%particle_rank)
    write(*,"(a,f18.8,a)") " Reduced EDM = ", edm, " efm"
    edm = edm*geometry_part(params%jket, 2, params%jket, params%jket, 0, params%jket) ! from D operator
    edm = edm*geometry_part(params%jket, 0, params%jket, params%jket, 0, params%jket) ! from PVTV operator
    edm = edm* sqrt(4.d0 * pi / 3.d0)
    write(*,"(a,f18.8,a)") " EDM = ", edm, " efm"
    call E1_2%fin()
    call PVTV_2%fin()
    call jacobi_ch%fin()
    call rel%fin()
  end function calc_edm_isospin_lec

  subroutine calc_tdm_isospin(params, egs, gs)
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    real(8) :: hpi1_lo, hpi1_nlo, hpi1_n2lo, Ct1, Ct2, Ct3, Ct4, Ct5, hV0, hV1, hV2, hA1, hA2
    type(sys) :: s
    hpi1_lo   = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hpi1-LO-Local2-500") )
    hpi1_nlo  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hpi1-NLO-Local2-500") )
    hpi1_n2lo = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hpi1-N2LO-Local2-500") )
    Ct1  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_Ct1-NLO-Local2-500") )
    Ct2  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_Ct2-NLO-Local2-500") )
    Ct3  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_Ct3-NLO-Local2-500") )
    Ct4  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_Ct4-NLO-Local2-500") )
    Ct5  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_Ct5-NLO-Local2-500") )
    hV0  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hV0-N2LO-Local2-500") )
    hV1  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hV1-N2LO-Local2-500") )
    hV2  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hV2-N2LO-Local2-500") )
    hA1  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hA1-N2LO-Local2-500") )
    hA2  = calc_tdm_isospin_lec( params, egs, gs, s%str("PVTCint_ChEFT_hA2-N2LO-Local2-500") )
    write(*,*) "egs, hpi1_lo, hpi1_nlo, hpi1_n2lo, Ct1, Ct2, Ct3, Ct4, Ct5, hV0, hV1, hV2, hA1, hA2"
    write(*,'(14f16.6)') egs, hpi1_lo, hpi1_nlo, hpi1_n2lo, Ct1, Ct2, Ct3, Ct4, Ct5, &
        & hV0, hV1, hV2, hA1, hA2
  end subroutine calc_tdm_isospin

  function calc_tdm_isospin_lec(params, egs, gs, opname) result(tdm)
    use LinAlgLib
    use ClassSys, only: sys
    use MyLibrary, only: geometry_part, pi
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use TwoBodyRelOps
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    use NNNForceHOIsospin
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    type(str), intent(in) :: opname
    type(DVec) :: groundstate, inter_energies
    type(DMat) :: propagator, inter_states
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelSpaceHOBasis) :: rel_pn
    type(ThreeBodyJacIsoChan) :: jacobi_ch, jacobi_ch_flip
    type(TwoBodyRelOpIso) :: TDM_is_2, TDM_iv_2, PVTC_2
    type(ThreeBodyJacOpChanIso) :: TDM_is_3, TDM_iv_3, PVTC_3
    integer :: ttot_flip, i
    real(8) :: me_is, me_iv, tdm
    type(sys) :: s
    ! for details
    !real(8) :: me_tmp1, me_tmp2
    !type(DVec) :: tmp

    call groundstate%ini(size(gs))
    groundstate%v = gs
    call jacobi_ch%init(params%hw, params%jket, params%pket, params%tket, params%N3max, params%path_to_tmp_dir)
    call rel%init(params%hw, params%N3max, params%Jmax2)
    call rel_pn%init(params%hw, params%N3max, params%Jmax2)

    call TDM_is_2%init(rel,rel,s%str("TDM_mag_s"))
    call TDM_iv_2%init(rel,rel,s%str("TDM_mag_st"))
    call PVTC_2%init(rel,rel,opname)
    ! set two-body operators
    call TDM_is_2%set()
    call TDM_iv_2%set()
    call PVTC_2%set()

    tdm = 0.d0
    do ttot_flip = 1, 3, 2
      call jacobi_ch_flip%init(params%hw,params%jket,-params%pket,ttot_flip,params%N3max,params%path_to_tmp_dir)
      call solve_three_body_hamil_isospin(params,params%jket,-params%pket,ttot_flip,params%tzket,inter_energies,inter_states)
      call propagator%zeros( jacobi_ch_flip%GetNumberAStates(), jacobi_ch_flip%GetNumberAStates() )
      do i = 1, jacobi_ch_flip%GetNumberAStates()
        propagator%m(i,i) = 1.d0 / (egs - inter_energies%v(i))
      end do
      propagator = inter_states * propagator * inter_states%t()

      call TDM_is_3%init(jacobi_ch, jacobi_ch_flip, s%str("TDM_mag_s"))
      call TDM_iv_3%init(jacobi_ch, jacobi_ch_flip, s%str("TDM_mag_st"))
      call TDM_is_3%set(TDM_is_2)
      call TDM_iv_3%set(TDM_iv_2)
      call PVTC_3%init(jacobi_ch_flip, jacobi_ch, s%str("PVTCint"))
      call PVTC_3%set(PVTC_2)

      ! detail
      !call tmp%ini( jacobi_ch_flip%GetNumberAStates() )
      !do i = 1, jacobi_ch_flip%GetNumberAStates()
      !  tmp%v(:) = inter_states%m(:,i)
      !  me_tmp1 = groundstate * TDM_iv_3%DMat * tmp * &
      !      &   geometry_part(params%tket, 2*TDM_iv_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket)
      !  me_tmp2 = tmp * PVTC_3%DMat * groundstate * &
      !      &   geometry_part(ttot_flip, 2*PVTC_2%GetOpT(), params%tket, params%tzket, 0, params%tzket)
      !  !call tmp%prt("inter")
      !  !call PVTC_3%prt("pv int")
      !  !call groundstate%prt("gs")
      !  write(*,'(a,3f12.6)') "energy denominator, IV, PVTC int.:", egs-inter_energies%v(i), me_tmp1, me_tmp2
      !end do
      !call tmp%fin()

      me_is = groundstate * TDM_is_3%DMat * propagator * PVTC_3%DMat * groundstate
      me_iv = groundstate * TDM_iv_3%DMat * propagator * PVTC_3%DMat * groundstate
      me_is = me_is*geometry_part(ttot_flip, 2*PVTC_2%GetOpT(), params%tket, params%tzket, 0, params%tzket) ! from PVTC interaction
      me_is = me_is*geometry_part(params%tket, 2*TDM_is_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket) ! from operator
      me_iv = me_iv*geometry_part(ttot_flip, 2*PVTC_2%GetOpT(), params%tket, params%tzket, 0, params%tzket) ! from PVTC interaction
      me_iv = me_iv*geometry_part(params%tket, 2*TDM_iv_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket) ! from operator
      tdm = tdm + me_is + me_iv

      call TDM_is_3%fin(); call TDM_iv_3%fin(); call PVTC_3%fin()
      call propagator%fin()
      call inter_energies%fin()
      call inter_states%fin()
      call jacobi_ch_flip%fin()
    end do
    tdm = tdm * 2.d0 / dble(params%particle_rank)
    write(*,"(a,f18.8,a)") " Reduced spin TDM = ", tdm, " efm2"
    tdm = tdm*geometry_part(params%jket, 2, params%jket, params%jket, 0, params%jket) ! from operator
    tdm = tdm*geometry_part(params%jket, 0, params%jket, params%jket, 0, params%jket) ! from PVTC operator
    tdm = tdm
    write(*,"(a,f18.8,a)") " TDM = ", tdm, " efm2"
    call TDM_is_2%fin()
    call TDM_iv_2%fin()
    call PVTC_2%fin()
    call jacobi_ch%fin()
    call rel%fin()
  end function calc_tdm_isospin_lec

  function calc_dipole_polarizability_isospin(params, egs, gs) result(alphaD)
    use LinAlgLib
    use ClassSys, only: sys
    use MyLibrary, only: geometry_part, hc, alpha, pi
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use TwoBodyRelOps
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsChanIso
    use NNNForceHOIsospin
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    type(DVec) :: groundstate, inter_energies
    type(DMat) :: propagator, inter_states
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelSpaceHOBasis) :: rel_pn
    type(ThreeBodyJacIsoChan) :: jacobi_ch, jacobi_ch_flip
    type(TwoBodyRelOpIso) :: E1_2
    type(ThreeBodyJacOpChanIso) :: E1_3
    integer :: ttot_flip, j_partner, i
    real(8) :: me, alphaD
    type(sys) :: s

    call groundstate%ini(size(gs))
    groundstate%v = gs
    call jacobi_ch%init(params%hw, params%jket, params%pket, params%tket, params%N3max,params%path_to_tmp_dir)
    call rel%init(params%hw, params%N3max, params%Jmax2)
    call rel_pn%init(params%hw, params%N3max, params%Jmax2)

    call E1_2%init(rel,rel,s%str("E1"))
    call E1_2%set()

    alphaD = 0.d0
    do j_partner = abs(params%jket - 2), abs(params%jket + 2), 2
      do ttot_flip = 1, 3, 2
        call jacobi_ch_flip%init(params%hw,j_partner,-params%pket,ttot_flip,params%N3max,params%path_to_tmp_dir)
        call solve_three_body_hamil_isospin(params,j_partner,-params%pket,ttot_flip,params%tzket,inter_energies,inter_states)
        call propagator%zeros( jacobi_ch_flip%GetNumberAStates(), jacobi_ch_flip%GetNumberAStates() )
        do i = 1, jacobi_ch_flip%GetNumberAStates()
          propagator%m(i,i) = 1.d0 / ( inter_energies%v(i) - egs )
        end do

        propagator = inter_states * propagator * inter_states%t()

        call E1_3%init(jacobi_ch, jacobi_ch_flip, s%str("E1"))
        call E1_3%set(E1_2)
        me = groundstate * E1_3%DMat * propagator * E1_3%DMat%T() * groundstate
        me = me * (-1.d0)**( (j_partner-params%jket+ttot_flip-params%tket)/2 )
        me = me*geometry_part(params%tket, 2*E1_2%GetOpT(), ttot_flip, params%tzket, 0, params%tzket) ! from D operator
        me = me*geometry_part(ttot_flip, 2*E1_2%GetOpT(), params%tket, params%tzket, 0, params%tzket) ! from D operator
        me = me*geometry_part(params%jket, 2*E1_2%GetOpJ(), j_partner, params%jket, 0, params%jket) ! from D operator
        me = me*geometry_part(j_partner, 2*E1_2%GetOpJ(), params%jket, params%jket, 0, params%jket) ! from D operator
        alphaD = alphaD + me

        call E1_3%fin()
        call propagator%fin()
        call inter_energies%fin()
        call inter_states%fin()
        call jacobi_ch_flip%fin()
      end do
    end do
    alphaD = alphaD * 2.d0 * hc / (alpha * dble(params%particle_rank)**2) * 4.d0 * pi / 3.d0
    write(*,"(a,f18.8,a)") " alphaD = ", alphaD, " fm3"
    call E1_2%fin()
    call jacobi_ch%fin()
    call rel%fin()
  end function calc_dipole_polarizability_isospin

  subroutine check_induced_nnn_part(params)
    use MyLibrary, only: triag
    use LinAlgLib
    use ClassSys, only: sys
    use OperatorDefinitions
    use NNNForceHOIsospin
    use TMTransformManager
    type(InputParameters), intent(in) :: params
    type(sys) :: s
    integer :: jbra, pbra, tbra, tzbra
    integer :: jket, pket, tket, tzket
    type(ThreeBodyJacOpChanIso) :: opjac_vs0, opjac_vs1, opjac_vs2, U_jac, opjac_T
    type(NNNForceIsospin) :: v3n_jac
    integer :: wut = 20
    logical :: calc_bra=.true.
    type(ThreeBodyJacIsoChan) :: chbra, chket
    type(str) :: f
    integer :: bra, ket
    type(DVec) :: energies
    type(DMat) :: states

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    tzbra= params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    tzket= params%tzket

    f = "test.dat"
    if(s%isfile(f)) then
      write(*,'(2a)') trim(f%val), ' already exists.'
      return
    end if

    if(jbra == jket .and. pbra == pket .and. tbra == tket) calc_bra = .false.
    if(calc_bra) then
      write(*,*) "return: check_induced_nnn_part"
      return
    end if
    call solve_three_body_hamil_isospin(params,jket,pket,tket,tzket,energies,states)
    call chbra%init(params%hw,jbra,pbra,tbra,params%N3max,params%path_to_tmp_dir)
    call chket%init(params%hw,jket,pket,tket,params%N3max,params%path_to_tmp_dir)

    ! --- 3b evolved ---
    call U_jac%init(chket,chket,s%str("UT"))
    call v3n_jac%init(chket)
    call v3n_jac%SetNNNForce(U_jac,params)
    call opjac_T%init(chbra,chket,s%str("NNNfromTkin"))
    call opjac_vs0%init(chbra,chket,s%str("NNNfrom2N_central"))
    call opjac_vs1%init(chbra,chket,s%str("NNNfrom2N_spinorbit"))
    call opjac_vs2%init(chbra,chket,s%str("NNNfrom2N_tensor"))

    call set_induced_nnn_ch(opjac_T, params, s%str("NNNfromTkin"))
    call set_induced_nnn_ch(opjac_vs0, params, s%str("NNNfrom2N_central"))
    call set_induced_nnn_ch(opjac_vs1, params, s%str("NNNfrom2N_spinorbit"))
    call set_induced_nnn_ch(opjac_vs2, params, s%str("NNNfrom2N_tensor"))

    open(wut, file=f%val, status="replace")
    do bra = 1, chbra%GetNumberAStates()
      do ket = 1, chket%GetNumberAStates()
        write(wut,"(2i6,5es18.6)") bra, ket, &
            & opjac_T%m(bra,ket), opjac_vs0%m(bra,ket), opjac_vs1%m(bra,ket), opjac_vs2%m(bra,ket), &
            & v3n_jac%m(bra,ket)
      end do
    end do
    close(wut)

    call opjac_vs0%fin()
    call opjac_vs1%fin()
    call opjac_vs2%fin()
    call opjac_T%fin()
  end subroutine check_induced_nnn_part

  subroutine check_norms(params)
    use ClassSys, only: sys
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use NNNForceHOIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace, only: ThreeBodyJacIsoChan, GetRampNmax
    use ThreeBodyJacOpsChanIso
    !
    type(InputParameters), intent(in) :: params
    type(InputParameters) :: input
    type(ThreeBodyJacIsoChan) :: jac
    type(ThreeBodyJacOpChanIso) :: U, T_jac, vnn_jac, T_sub_jac, vnn_sub_jac, Generator
    type(ThreeBodyJacOpChanIso) :: T_tr, v2_tr, v3_tr
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: vnn, vnn_sub, t2, t2_sub
    type(DMat) :: mat_t, mat_tind, mat_v2, mat_v2ind, mat_v3, mat_v3tr, h, eta0, rhs0, h0
    type(NNNForceIsospin) :: v3n_jac
    type(str) :: fn
    type(sys) :: s
    integer :: wunit=21
    integer :: Nmax
    integer :: Nbra, Nket
    if(params%jbra /= params%jket) return
    if(params%pbra /= params%pket) return
    if(params%tbra /= params%tket) return
    Nmax = params%N3max
    call jac%init(params%hw, params%jket, params%pket, params%tket, Nmax, params%path_to_tmp_dir)

    call relspin%init(params%hw, Nmax, params%J2max_NNint)
    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = Nmax; input%renorm="bare"
    call vnn_relspin%setNNForceHOIsospin(U_relspin,input)

    call rel%init(params%hw, Nmax, params%Jmax2)
    call vnn%init(rel,rel,s%str("NNint"))
    call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)

    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    call t2%init(rel,rel,s%str("kinetic"))
    call t2%set()
    t2_sub = t2
    vnn_sub = vnn
    input%renorm=params%renorm
    call vnn_sub%evolve(input, jac%GetNmax(), jac%GetNmax())
    call t2_sub%evolve(input, jac%GetNmax(), jac%GetNmax())

    call v3n_jac%Init(jac)
    call U%init(jac,jac,s%str("UT"))
    call T_jac%init(jac,jac,s%str("kinetic"))
    call T_jac%set(t2)
    call T_sub_jac%init(jac,jac,s%str("kinetic"))
    call T_sub_jac%set(t2_sub)
    call vnn_jac%init(jac,jac,s%str("NNint"))
    call vnn_jac%set(vnn)
    call vnn_sub_jac%init(jac,jac,s%str("NNint"))
    call vnn_sub_jac%set(vnn_sub)
    if(params%genuine_3bf) then
      call set_genuine_three_body_interaction(v3n_jac,params,verbose=.true.)
    end if
    mat_t = T_jac%get_norm_op_chan(  mode=s%str("average"))
    mat_v2= vnn_jac%get_norm_op_chan(mode=s%str("average"))
    mat_v3= v3n_jac%get_norm_op_chan(mode=s%str("average"))
    h = T_jac%DMat + vnn_jac%DMat + v3n_jac%DMat
    h0 = h
    eta0 = (T_jac%DMat * h) - (h * T_jac%DMat)
    rhs0 = (eta0*h) - (h*eta0)
    if(params%renorm%val == 'srg') then
      select case( params%srg_generator%val )
      case( "kinetic" )
        Generator = T_jac
      case default
        Generator = set_srg_generator( jac, params%srg_generator )
      end select
      call three_body_srg_evolution(h, Generator, U%DMat, params%lambda, params%Nmax_srg_edge)
    end if

    call T_tr%init(jac,jac,s%str("kinetic"))
    call v2_tr%init(jac,jac,s%str("NNint"))
    call v3_tr%init(jac,jac,s%str("NNNint"))
    t_tr%DMat = ( U%DMat%T() * T_jac%DMat * U%DMat ) - (1.5d0 * T_sub_jac%DMat) + (0.5d0 * T_jac%DMat)
    v2_tr%DMat = ( U%DMat%T() * vnn_jac%DMat * U%DMat ) - vnn_sub_jac%DMat
    v3_tr%DMat = U%DMat%T() * v3n_jac%DMat * U%DMat

    mat_tind  = t_tr%get_norm_op_chan( mode=s%str("average"))
    mat_v2ind = v2_tr%get_norm_op_chan(mode=s%str("average"))
    mat_v3tr  = v3_tr%get_norm_op_chan(mode=s%str("average"))
    fn = "Averaged_" + params%averaged_file_for_test + ".txt"
    open(wunit, file=fn%val)
    write(wunit,"(a)") "# E',   E,   T bare, NN bare, 3N bare, Tind, NNind, 3N"
    do Nbra = 0, Nmax
      if((-1)**Nbra /= params%pket) cycle
      do Nket = 0, Nmax
        if((-1)**Nket /= params%pket) cycle
        write(wunit,"(2i4,6es18.6)") Nbra, Nket, &
            & mat_t%m(Nbra+1, Nket+1), mat_v2%m(Nbra+1, Nket+1), &
            & mat_v3%m(Nbra+1, Nket+1), mat_tind%m(Nbra+1, Nket+1), &
            & mat_v2ind%m(Nbra+1, Nket+1), mat_v3tr%m(Nbra+1, Nket+1)
      end do
    end do
    close(wunit)

    fn = "eta0_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) eta0%m(:,:)
    close(wunit)

    fn = "rhs0_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) rhs0%m(:,:)
    close(wunit)

    fn = "H0_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) H0%m(:,:)
    close(wunit)

    fn = "Hevol_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) h%m(:,:)
    close(wunit)

    fn = "T_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) T_jac%DMat%m(:,:)
    close(wunit)
    fn = "Tind_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) t_tr%m(:,:)
    close(wunit)
    fn = "NN_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) vnn_jac%DMat%m(:,:)
    close(wunit)
    fn = "NNind_" + params%averaged_file_for_test + ".bin"
    open(wunit, file=fn%val, form="unformatted", access="stream")
    write(wunit) v2_tr%m(:,:)
    close(wunit)
    if(params%genuine_3bf) then
      fn = "NNN_" + params%averaged_file_for_test + ".bin"
      open(wunit, file=fn%val, form="unformatted", access="stream")
      write(wunit) v3n_jac%DMat%m(:,:)
      close(wunit)
      fn = "NNNEvol_" + params%averaged_file_for_test + ".bin"
      open(wunit, file=fn%val, form="unformatted", access="stream")
      write(wunit) v3_tr%m(:,:)
      close(wunit)
    end if
    call T_jac%fin()
    call T_sub_jac%fin()
    call vnn_jac%fin()
    call vnn_sub_jac%fin()
    call v3n_jac%fin()
    call u%fin()
    call jac%fin()
  end subroutine check_norms

  subroutine check_norms_old(params)
    use ClassSys, only: sys
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace, only: ThreeBodyJacIsoChan, GetRampNmax
    use ThreeBodyJacOpsChanIso
    !
    type(InputParameters), intent(in) :: params
    type(sys) :: s
    type(InputParameters) :: input
    type(ThreeBodyJacIsoChan) :: jac
    type(ThreeBodyJacOpChanIso) :: ut, vtot, v3, tmp, v2, t, t_sub
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: vnn, vnn_sub, t2, t2_sub
    type(DMat) :: fmat_tind, fmat_v2, fmat_v3, fmat_tot
    type(str) :: f, fv, fut
    integer :: runit=20, wunit=21
    integer :: Nmax
    integer :: Nbra, Nket
    if(params%jbra /= params%jket) return
    if(params%pbra /= params%pket) return
    if(params%tbra /= params%tket) return
    Nmax = GetRampNmax(params%jket, params%ramp)
    f = jac%GetFileName(params%jket, params%pket, params%tket, Nmax, params%path_to_tmp_dir)
    fv = ut%GetFileName(jac,params%hw,s%str("NNNint"),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)
    fut = ut%GetFileName(jac,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)

    if(.not. s%isfile(f)) then
      write(*,*) "cfp file is not found! Do actual calcualtions first!"
      return
    end if
    if(.not. s%isfile(fv)) then
      write(*,*) "3N int file is not found! Do actual calcualtions first!"
      return
    end if
    if(.not. s%isfile(fut)) then
      write(*,*) "3N unitrary transformation file is not found! Do actual calcualtions first!"
      return
    end if
    call jac%init(params%hw, params%jket, params%pket, params%tket, Nmax, params%path_to_tmp_dir)
    open(runit,file=f%val,status='old',form='unformatted',access='stream')
    call jac%readf(params%hw,runit,Nmax)
    close(runit)

    call vtot%init(jac,jac,s%str("hamil"))
    call ut%init(jac,jac,s%str("UT"))
    open(runit,file=fv%val,status='old',form='unformatted',access='stream')
    call vtot%readf(jac, runit, jac%GetNmax(), params%hw, params%hw_target, params%hw_conversion)
    close(runit)
    open(runit,file=fut%val,status='old',form='unformatted',access='stream')
    call ut%readf(jac, runit, jac%GetNmax(), params%hw, params%hw_target, params%hw_conversion)
    close(runit)

    call relspin%init(params%hw, Nmax, params%J2max_NNint)
    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = Nmax; input%renorm="bare"
    call vnn_relspin%setNNForceHOIsospin(U_relspin,input)

    call rel%init(params%hw, Nmax, params%Jmax2)
    call vnn%init(rel,rel,s%str("NNint"))
    call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)

    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    call t2%init(rel,rel,s%str("kinetic"))
    call t2%set()
    t2_sub = t2
    vnn_sub = vnn
    input%renorm=params%renorm
    call vnn_sub%evolve(input, jac%GetNmax(), jac%GetNmax())
    call t2_sub%evolve(input, jac%GetNmax(), jac%GetNmax())

    call t%init(jac, jac, s%str("kinetic"))
    call tmp%init(jac, jac, s%str("kinetic"))
    call tmp%set(t2)
    t%DMat = UT%DMat%T() * tmp%DMat * UT%DMat
    call t_sub%init(jac, jac, s%str("kinetic"))
    call t_sub%set(t2_sub)
    t%DMat = t%DMat - (t_sub%DMat * 1.5d0)
    t%DMat = t%DMat + 0.5d0 * tmp%DMat
    if(params%hw_conversion) call t%FreqConv(jac%GetNmax(), jac%GetNmax(), params%hw, params%hw_target)

    open(runit,file=fut%val,status='old',form='unformatted',access='stream')
    call v2%init(jac,jac,s%str("NNint"))
    call v2%set(jac, jac, runit, runit, vnn, vnn_sub, &
        & Nmax, Nmax, params%hw, params%hw_target, params%hw_conversion)
    close(runit)
    call v3%init(jac,jac,s%str("NNNint"))
    v3%DMat = vtot%DMat - v2%DMat - t%DMat

    fmat_tind= t%get_norm_op_chan( mode=s%str("average"))
    fmat_v2  = v2%get_norm_op_chan(mode=s%str("average"))
    fmat_v3  = v3%get_norm_op_chan(mode=s%str("average"))
    fmat_tot = vtot%get_norm_op_chan(mode=s%str("average"))
    f = params%file_name_3n
    open(wunit, file=f%val)
    write(wunit,"(a)") "# E',   E,   Tind, NNind, 3N,  Total"
    do Nbra = 0, Nmax
      if((-1)**Nbra /= params%pket) cycle
      do Nket = 0, Nmax
        if((-1)**Nket /= params%pket) cycle
        write(wunit,"(2i4,4es18.6)") Nbra, Nket, &
            & fmat_tind%m(Nbra+1, Nket+1), fmat_v2%m(Nbra+1, Nket+1), &
            & fmat_v3%m(Nbra+1, Nket+1), fmat_tot%m(Nbra+1, Nket+1)
      end do
    end do
    close(wunit)
    call t%fin()
    call v2%fin()
    call v3%fin()
    call vtot%fin()
    call ut%fin()
    call jac%fin()
  end subroutine check_norms_old

  subroutine file_convert_3n(params)
    use SingleParticleState
    use ThreeBodyLabSpace
    ! temporary
    use ThreeBodyLabOpsIsoSingle, only: ThreeBodyLabOpIso => ThreeBodyLabOpIsoSingle
    type(InputParameters), intent(in) :: params
    type(ThreeBodyLabOpIso) :: oplab, oplab_convert
    type(OrbitsIsospin) :: sps, sps_convert
    type(ThreeBodyLabIsoSpace) :: lab, lab_convert
    type(str) :: filename

    if(params%file_name_3n_original%val == "none") then
      write(*,*) "Error, set the original file name."
      return
    end if

    write(*,*) ""
    write(*,"(a)") "#################### 3N file convert mode ####################"
    write(*,*) ""

    call sps%init(params%emax, params%lmax)
    call lab%init(sps, params%e2max, params%e3max, params%hw)
    call oplab%init(lab, params%Operators(1))

    filename = params%file_name_3n_converted
    write(*,*) "Reading..."
    call oplab%readf(params%file_name_3n_original)
    write(*,*) "Done."

    call sps_convert%init(params%emax_convert, params%lmax_convert)
    call lab_convert%init(sps_convert, params%e2max_convert, params%e3max_convert, params%hw)
    oplab_convert = oplab%Truncate( lab_convert )
    call oplab%fin()
    call lab%fin()
    call sps%fin()
    write(*,*) "Writing..."
    if(params%averaged_file_for_test%val == "none") call oplab_convert%writef(filename)
    if(params%averaged_file_for_test%val /= "none") call oplab_convert%writef(filename, params%averaged_file_for_test)
    call oplab_convert%fin()
    call lab_convert%fin()
    call sps_convert%fin()

  end subroutine file_convert_3n

  subroutine file_convert_no2b(params)
    use SingleParticleState
    use ThreeBodyLabSpaceIsoNO2B
    !temp
    use ThreeBodyNO2BIsospinSingle, only: &
        & ThreeBodyNO2BChIso => ThreeBodyNO2BChIsoSingle, &
        & ThreeBodyNO2BIso => ThreeBodyNO2BIsoSingle, &
        & ThreeBodyLabOpIsoChan => ThreeBodyLabOpIsoChanSingle
    !temp
    type(InputParameters), intent(in) :: params
    type(ThreeBodyNO2BIso) :: oplab
    type(OrbitsIsospin) :: sps
    type(NO2BThreeBodyIsoSpace) :: lab
    type(str) :: filename

    if(params%file_name_3n_original%val == "none") then
      write(*,*) "Error, set the original file name."
      return
    end if

    write(*,*) ""
    write(*,"(a)") "#################### 3N file convert mode ####################"
    write(*,"(a)") "# Don't use for tensor operators"

    call sps%init(params%emax, params%lmax)
    call lab%init(sps, params%e2max, params%e3max)
    call oplab%init(lab)

    filename = params%file_name_3n_converted

    call oplab%readf(params%file_name_3n_original)
    call oplab%writef(filename)

    call oplab%fin()
    call lab%fin()
    call sps%fin()
  end subroutine file_convert_no2b

  subroutine set_induced_nnn_ch(v, params, opname)
    use TwoBodyRelativeSpace
    use TwoBodyRelOpsIso
    use NNForceIsospin
    use ThreeBodyJacOpsChanIso
    use ThreeBodyJacOpsIso
    use ThreeBodyJacobiSpace
    type(ThreeBodyJacOpChanIso), intent(inout) :: v
    type(InputParameters), intent(in) :: params
    type(str), intent(in) :: opname
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: vnn, op2, op2_sub, t_rel
    type(ThreeBodyJacIsoChan) :: chjac
    type(InputParameters) :: input
    type(ThreeBodyJacOpChanIso) :: tmp, ut, t_orig, t_eff, vnn_jac
    integer :: iut = 101, Nmax_file
    logical :: ex
    type(sys) :: s
    type(str) :: fut, fn_jacobi

    if(params%renorm%val == "bare") return
    if(v%jacobi_ch_bra%GetNumberAStates() < 1) return
    if(v%jacobi_ch_ket%GetNumberAStates() < 1) return
    Nmax_file = GetRampNmax(v%jacobi_ch_ket%GetJ(), params%ramp)
    fn_jacobi = chjac%GetFileName(v%jacobi_ch_ket%GetJ(), &
        &  v%jacobi_ch_ket%GetParity(), v%jacobi_ch_ket%GetT(), Nmax_file, params%path_to_tmp_dir)
    ex = s%isfile(fn_jacobi, s%str('In set_three_body_jac_operator_isospin'))
    open(iut, file=fn_jacobi%val, form='unformatted', access='stream')
    call chjac%readf(params%hw,iut)
    close(iut)

    call relspin%init(params%hw, params%N3max, params%J2max_NNint)
    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = params%N3max; input%renorm="bare"
    call vnn_relspin%setNNForceHOIsospin(U_relspin,input)

    call rel%init(params%hw, params%N3max, params%Jmax2)
    call vnn%init(rel,rel,s%str("NNint"))
    call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)

    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    select case(opname%val)
    case("NNNfrom2N","NNNinduced")
      call op2%init(rel,rel,opname)
      op2 = vnn
    case("NNNfrom2N_central")
      call op2%init(rel,rel,opname)
      op2 = vnn%SpinTensorDecomposition(0)
    case("NNNfrom2N_spinorbit")
      call op2%init(rel,rel,opname)
      op2 = vnn%SpinTensorDecomposition(1)
    case("NNNfrom2N_tensor")
      call op2%init(rel,rel,opname)
      op2 = vnn%SpinTensorDecomposition(2)
    case("NNNfromTkin")
      call op2%init(rel,rel,s%str('kinetic'))
      call op2%set()
    case default
      return
    end select
    op2_sub = op2
    input%renorm=params%renorm
    call op2_sub%evolve(input, chjac%GetNmax(), chjac%GetNmax())
    fut = ut%GetFileName(chjac,params%hw_target,s%str("UT"),params%genuine_3bf, params%regulator, &
        & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
        & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
        & params%path_to_tmp_dir)
    if(params%hw_conversion) then
      fut = ut%GetFileName(chjac,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
    end if
    ex = s%isfile(fut, s%str('In set_three_body_jac_operator_isospin'))
    open(iut, file=fut%val, form='unformatted', access='stream')
    select case(opname%val)
    case("NNNfromTkin")
      call UT%init(chjac, chjac, s%str('UT'))
      call UT%readf(chjac, iut, params%N3max, params%hw, params%hw_target, params%hw_conversion)

      call t_orig%init(chjac, chjac, s%str("kinetic"))
      call t_orig%set(op2)
      t_orig%DMat = UT%DMat%T() * t_orig%DMat * UT%DMat
      call UT%release()

      call t_eff%init(chjac, chjac, s%str("kinetic"))
      call t_eff%set(op2_sub)
      t_orig%DMat = t_orig%DMat - (t_eff%DMat * 1.5d0)
      call t_eff%release()

      call tmp%init(chjac, chjac, s%str("kinetic"))
      call tmp%set(op2)
      t_orig%DMat = t_orig%DMat + 0.5d0 * tmp%DMat
      call tmp%release()
      if(params%hw_conversion) call t_orig%FreqConv(chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target)
      v%m(:,:) = t_orig%m(:v%jacobi_ch_bra%GetNumberAStates(), :v%jacobi_ch_ket%GetNumberAStates())
      call t_orig%release()

    case("NNNinduced")
      call UT%init(chjac, chjac, s%str('UT'))
      call UT%readf(chjac, iut, params%N3max, params%hw, params%hw_target, params%hw_conversion)
      call t_rel%init(rel,rel,s%str('kinetic'))
      call t_rel%set()
      ! -- 3N evolution
      call t_orig%init(chjac, chjac, s%str("kinetic"))
      call t_orig%set(t_rel)
      call vnn_jac%init(chjac, chjac, s%str("NNint"))
      call vnn_jac%set(op2)
      call tmp%init(chjac, chjac, s%str("hamil"))
      tmp%DMat = UT%DMat%T() * (t_orig%DMat + vnn_jac%DMat) * UT%DMat
      call vnn_jac%release()
      call vnn_jac%init(chjac, chjac, s%str("NNint"))
      !

      ! -- subtract
      op2_sub = op2 + t_rel
      call op2_sub%evolve(input, chjac%GetNmax(), chjac%GetNmax())
      op2_sub = op2_sub - t_rel
      call vnn_jac%set(op2_sub)
      tmp%DMat = tmp%DMat - vnn_jac%DMat - t_orig%DMat
      !
      call t_orig%release()
      call vnn_jac%release()
      if(params%hw_conversion) call tmp%FreqConv(chjac%GetNmax(), chjac%GetNmax(), params%hw, params%hw_target)
      v%m(:,:) = tmp%m(:v%jacobi_ch_bra%GetNumberAStates(), :v%jacobi_ch_ket%GetNumberAStates())
      call tmp%release()

    case default
      call v%set(chjac, chjac, iut, iut, op2, op2_sub, &
          & params%N3max, params%N3max, params%hw, params%hw_target, params%hw_conversion)
    end select
    close(iut)
    call vnn%fin()
    call op2%fin()
    call op2_sub%fin()
    call chjac%fin()
  end subroutine set_induced_nnn_ch
end module ThreeBodyManager
