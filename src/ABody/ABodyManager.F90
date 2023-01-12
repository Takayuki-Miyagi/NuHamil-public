module ABodyManager
  use omp_lib
  use ClassSys, only: sys, str
  use LinAlgLib
  use Profiler, only: timer
  use NuHamilInput, only: InputParameters, params
  use ThreeBodyJacobiSpace
  use ABodyJacSpaceIso
  use ABodyJacOpsIso
  implicit none
contains
  subroutine ABodyManagers()
    integer :: i
    type(InputParameters) :: input

    input = params
    if(params%pn_form) then
      write(*,"(a)") "For A-body calculations, turn off pn formalism."
      input%pn_form = .false.
    end if

    if(params%trans2lab) then
      write(*,"(2a)") "For more than 2-body and 3-body calculations", &
          & ", the transformation from Jacobi to Lab systems is not suppoted."
      stop
    end if

    if(mod(params%particle_rank,2) == 1) then
      if(mod(params%jbra, 2) == 0) stop 'Warning: Jbra have to be odd number'
      if(mod(params%jket, 2) == 0) stop 'Warning: Jket have to be odd number'
      if(mod(params%tbra, 2) == 0) stop 'Warning: tbra have to be odd number'
      if(mod(params%tket, 2) == 0) stop 'Warning: tket have to be odd number'
      if(mod(params%tzbra, 2) == 0) stop 'Warning: tzbra have to be odd number'
      if(mod(params%tzket, 2) == 0) stop 'Warning: tzket have to be odd number'
    end if
    if(abs(params%tzbra) > params%tbra) stop 'Warning: tz_bra have to be smaller than t_bra'
    if(abs(params%tzket) > params%tket) stop 'Warning: Tz_ket have to be smaller than t_ket'

    do i = 1, size(params%Operators)
      call calc_A_body(params, params%Operators(i))
    end do
  end subroutine ABodyManagers

  subroutine calc_A_body(params, oprtr)
    use MyLibrary, only: triag, geometry_part
    use OperatorDefinitions
    type(InputParameters), intent(in) :: params
    type(str), intent(in) :: oprtr
    type(sys) :: s
    type(OperatorDef) :: opdef
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    integer :: A, Z, N
    logical :: calc_bra
    type(DMat) :: bra_states, ket_states
    type(DVec) :: bra_energies, ket_energies
    real(8) :: vals(3), vals0(3), vals1(3)

    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket

    if(mod(params%particle_rank, 2) == 0) then
      jbra = 2*params%jbra
      tbra = 2*params%tbra
      zbra = 2*params%tzbra

      jket = 2*params%jket
      tket = 2*params%tket
      zket = 2*params%tzket

    end if
    A = params%particle_rank
    Z = (A-zket)/2
    N = (A+zket)/2

    select case(oprtr%val)
    case("Rp2", "Rn2","EDM","TDM","AlphaD")
      if(jbra/=jket) return
      if(pbra/=pket) return
      if(tbra/=tket) return
      if(zbra/=zket) return
    case default
      call opdef%InitOpDef(oprtr, .false.)
      if(triag(jbra,jket,2*opdef%GetOpJ())) return
      if(pbra * pket * opdef%GetOpP() /= 1) return
      if(triag(tbra,tket,2*opdef%GetOpT())) return
    end select

    calc_bra = .true.
    if(jbra==jket .and. pbra==pket .and. tbra==tket) calc_bra=.false.
    call solve_A_body_hamil(params,jket,pket,tket,zket,ket_energies,ket_states)
    bra_energies = ket_energies
    bra_states = ket_states
    if(calc_bra) then
      call solve_A_body_hamil(params,jbra,pbra,tbra,zbra,bra_energies,bra_states)
    end if

    select case(oprtr%val)
    case('hamil', 'Hamil')
      write(*,'(a,f18.8)') "Energy: ", ket_energies%v(1)
    case('EDM')
      if(jbra /= jket) return
      if(pbra /= pket) return
      if(tbra /= tket) return
      write(*,"(2a)") "Not implemented yet: operator =", trim(oprtr%val)
    case('TDM')
      if(jbra /= jket) return
      if(pbra /= pket) return
      if(tbra /= tket) return
      write(*,"(2a)") "Not implemented yet: operator =", trim(oprtr%val)
    case('AlphaD')
      if(jbra /= jket) return
      if(pbra /= pket) return
      if(tbra /= tket) return
      vals0(1) = calc_dipole_polarizability( params, ket_energies%v(1), ket_states%m(:,1))
      write(*,'(a,2f18.8)') "AlphaD: ", ket_energies%v(1), vals0(1)
    case("Rp2")
      vals0 = calc_expectation_val(s%str("R2"), params, bra_states, ket_states)
      vals1 = calc_expectation_val(s%str("R2_IV"), params, bra_states, ket_states)
      vals1 = vals1 * geometry_part(tbra, 2, tket, zbra, 0, zket)
      vals = dble(N) / dble(A**2 * Z) * vals0 - 1.d0 / dble(2 * A * Z) *  vals1
      write(*,'(a,5f18.8)') "Rp2: ", bra_energies%v(1), ket_energies%v(1), vals(:)
    case("Rn2")
      vals0 = calc_expectation_val(s%str("R2"), params, bra_states, ket_states)
      vals1 = calc_expectation_val(s%str("R2_IV"), params, bra_states, ket_states)
      vals1 = vals1 * geometry_part(tbra, 2, tket, zbra, 0, zket)
      vals = dble(Z) / dble(A**2 * N) * vals0 + 1.d0 / dble(2 * A * N) *  vals1
      write(*,'(a,5f18.8)') "Rn2: ", bra_energies%v(1), ket_energies%v(1), vals(:)
    case default
      vals = calc_expectation_val(oprtr, params, bra_states, ket_states)
      write(*,'(a,5f18.8)') trim(oprtr%val), bra_energies%v(1), ket_energies%v(1), vals(:)
    end select
    call bra_states%fin()
    call ket_states%fin()
    call bra_energies%fin()
    call ket_energies%fin()
  end subroutine calc_A_body

  function calc_expectation_val(opname, params, bra_vecs, ket_vecs) result(vals)
    use ClassSys, only: sys
    use MyLibrary, only: geometry_part
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsIso
    type(str), intent(in) :: opname
    type(InputParameters), intent(in) :: params
    type(DMat), intent(in) :: bra_vecs, ket_vecs
    integer :: jbra, pbra, tbra, zbra
    integer :: jket, pket, tket, zket
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(TwoBodyRelOpIso) :: U, op, opeff
    type(ThreeBodyJacOpIso) :: bare_3, eff2_3, eff3_3
    type(ABodyJacOpChanIso) :: bare_A, eff2_A, eff3_A
    type(InputParameters) :: input
    type(ThreeBodyJacIsoSpace) :: jacobi3
    type(ABodyJacIsoChan) :: chbra, chket
    real(8) :: exp_bare, exp_2b, exp_3b, red_to_non_red
    real(8) :: vals(3)
    logical :: calc_bra = .true.
    type(sys) :: s

    if(params%particle_rank > 4) then
      write(*,*) "Not implemented yet for the system A>4"
      return
    end if
    jbra = params%jbra
    pbra = params%pbra
    tbra = params%tbra
    zbra = params%tzbra

    jket = params%jket
    pket = params%pket
    tket = params%tket
    zket = params%tzket

    if(mod(params%particle_rank, 2) == 0) then
      jbra = 2*params%jbra
      tbra = 2*params%tbra
      zbra = 2*params%tzbra

      jket = 2*params%jket
      tket = 2*params%tket
      zket = 2*params%tzket
    end if

    vals = 0.d0
    if(jbra == jket .and. pbra == pket .and. tbra == tket .and. zbra==zket) calc_bra = .false.
    call relspin%init(params%hw,params%NAmax, params%J2max_NNint)
    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = params%N3max
    call vnn_relspin%setNNForceHOIsospin(U_relspin, input, params%Tket, params%Tzket)

    call rel%init(params%hw, params%NAmax, params%Jmax2)
    call U%init(rel,rel,s%str("UT"))
    call U%SetTwoBodyScalarSpinBreaking(U_relspin)
    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()

    call op%init(rel,rel,opname)
    !op%A = params%particle_rank; op%Z = (op%A-zket)/2; op%N = (op%A+zket)/2
    !call op%set(params%NNint, params%pmax2, params%NMesh2)
    call op%set()

    call opeff%init(rel,rel,opname)
    !opeff%A = params%particle_rank; opeff%Z = (opeff%A-zket)/2; opeff%N = (opeff%A+zket)/2
    !call opeff%set(params%NNint, params%pmax2, params%NMesh2)
    call opeff%set()
    call opeff%UT(U,U)
    call U%fin()

    call jacobi3%init(params%hw, 2*params%NAmax+3, params%NAmax, params%ramp, params%N3max, params%path_to_tmp_dir)
    if(params%particle_rank == 4) then
      call chbra%init(jbra,pbra,tbra,params%NAmax,jacobi3,params%path_to_tmp_dir)
      call chket%init(jket,pket,tket,params%NAmax,jacobi3,params%path_to_tmp_dir)
    end if

    if(params%particle_rank > 4) then
    end if

    ! --- bare ---
    bare_3 = get_3_body_jacobi_op_isospin(jacobi3, op, opname)
    eff2_3 = get_3_body_jacobi_op_isospin(jacobi3, opeff, opname)
    call eff3_3%init(jacobi3, opname)
    if(params%renorm%val /= "bare") call eff3_3%set(params)

    call bare_A%init(chbra, chket, opname)
    call eff2_A%init(chbra, chket, opname)
    call eff3_A%init(chbra, chket, opname)
    if(params%particle_rank == 4) then
      call bare_A%set(bare_3)
      call eff2_A%set(eff2_3)
      call eff3_A%set(eff3_3)
      bare_A%m = (dble( params%particle_rank * (params%particle_rank-1)) / 6.d0) * bare_A%m
      eff2_A%m = (dble( params%particle_rank * (params%particle_rank-1)) / 6.d0) * eff2_A%m
      eff3_A%m = (dble( params%particle_rank * (params%particle_rank-1) * (params%particle_rank-2)) / 6.d0) * eff3_A%m
    end if

    if(params%particle_rank > 4) then
    end if

    exp_bare= dot_product(bra_vecs%m(:,1), matmul(bare_A%m, ket_vecs%m(:,1)))
    exp_2b  = dot_product(bra_vecs%m(:,1), matmul(eff2_A%m, ket_vecs%m(:,1)))
    exp_3b  = dot_product(bra_vecs%m(:,1), matmul(eff3_A%m, ket_vecs%m(:,1)))
    red_to_non_red = 1.d0
    !red_to_non_red = geometry_part(jbra, 0, jket, jbra, 0, jket) * &
    !    &            geometry_part(tbra, 4, tket, zbra, (zbra-zket), zket)
    !red_to_non_red = geometry_part(tbra, 4, tket, zbra, (zbra-zket), zket)
    exp_bare= exp_bare * red_to_non_red
    exp_2b  = exp_2b   * red_to_non_red
    exp_3b  = exp_3b   * red_to_non_red
    vals(:) = [exp_bare, exp_2b, exp_2b+exp_3b]

    call bare_A%fin()
    call eff2_A%fin()
    call eff3_A%fin()

    call bare_3%fin()
    call eff2_3%fin()
    call eff3_3%fin()
    call jacobi3%fin()

    call opeff%fin()
    call op%fin()
    call rel%fin()
    call chket%fin()
    call chbra%fin()


  end function calc_expectation_val

  function calc_dipole_polarizability(params, egs, gs) result(alphaD)
    use ClassSys, only: sys
    use MyLibrary, only: geometry_part, alpha, hc, pi
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacobiSpace
    use ThreeBodyJacOpsIso
    type(InputParameters), intent(in) :: params
    real(8), intent(in) :: egs, gs(:)
    type(DVec) :: groundstate, inter_energies
    type(DMat) :: propagator, inter_states
    integer :: j, p, t, z, i, j_partner, t_partner
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelOpIso) :: op
    type(ThreeBodyJacOpIso) :: bare_3
    type(ABodyJacOpChanIso) :: bare_A
    type(ThreeBodyJacIsoSpace) :: jacobi3
    type(ABodyJacIsoChan) :: chjac, chjac_partner
    type(sys) :: s
    real(8) :: exp_bare
    real(8) :: alphaD

    AlphaD = 0.d0
    if(params%particle_rank > 4) then
      write(*,*) "Not implemented yet for the system A>4"
      return
    end if
    j = 2*params%jket
    p =   params%pket
    t = 2*params%tket
    z = 2*params%tzket

    call groundstate%ini(size(gs))
    groundstate%v = gs
    call rel%init(params%hw, params%NAmax, params%Jmax2)
    call op%init(rel,rel,s%str("E1"))
    call op%set()

    call chjac%init(j,p,t,params%NAmax,jacobi3,params%path_to_tmp_dir)
    call jacobi3%init(params%hw, 2*params%NAmax+3, params%NAmax, params%ramp, params%N3max, params%path_to_tmp_dir)
    bare_3 = get_3_body_jacobi_op_isospin(jacobi3, op, s%str("E1"))

    alphaD = 0.d0
    do j_partner = abs(j-2*op%GetOpJ()), j+2*op%GetOpJ(), 2 ! assuming A=even
      do t_partner = abs(t-2*op%GetOpT()), t+2*op%GetOpT(), 2 ! assuming A=even
        if( j_partner < 0 ) cycle
        if( t_partner < 0 ) cycle
        call chjac_partner%init(j_partner,-p,t_partner,params%NAmax,jacobi3,params%path_to_tmp_dir)

        call solve_A_body_hamil(params,j_partner,-p,t_partner,z,inter_energies,inter_states)
        call propagator%zeros( chjac_partner%GetNumberAStates(), chjac_partner%GetNumberAStates() )
        do i = 1, chjac_partner%GetNumberAStates()
          propagator%m(i,i) = 1.d0 / ( inter_energies%v(i) - egs )
        end do
        propagator = inter_states * propagator * inter_states%t()

        call bare_A%init(chjac, chjac_partner, s%str("E1"))
        call bare_A%set(bare_3)
        bare_A%m = (dble( params%particle_rank * (params%particle_rank-1)) / 6.d0) * bare_A%m
        exp_bare = groundstate * bare_A%DMat * propagator * bare_A%DMat%T() * groundstate
        exp_bare = exp_bare * (-1.d0)**( (j_partner-j+t_partner-t) ) ! assumed A=even
        exp_bare = exp_bare*geometry_part(t, 2, t_partner, z, 0, z) ! from D operator
        exp_bare = exp_bare*geometry_part(t_partner, 2, t, z, 0, z) ! from D operator
        exp_bare = exp_bare*geometry_part(j, 2, j_partner, j, 0, j) ! from D operator
        exp_bare = exp_bare*geometry_part(j_partner, 2, j, j, 0, j) ! from D operator
        alphaD = alphaD + exp_bare

        call propagator%fin()
        call inter_energies%fin()
        call bare_A%fin()
        call chjac_partner%fin()
      end do
    end do
    alphaD = alphaD * 2.d0 * hc / ( alpha * dble(params%particle_rank)**2 ) * 4.d0 * pi / 3.d0
    write(*,"(a,f18.8,a)") " alphaD = ", alphaD, " fm3"
    call bare_3%fin()
    call jacobi3%fin()
    call op%fin()
    call rel%fin()
    call chjac%fin()
  end function calc_dipole_polarizability

  subroutine solve_A_body_hamil(params, j, p, t, z, energies, states)
    use TwoBodyRelativeSpace
    use NNForceIsospin
    use TwoBodyRelOpsIso
    use ThreeBodyJacOpsIso
    type(InputParameters), intent(in) :: params
    integer, intent(in) :: j, p, t, z
    type(DVec), intent(out) :: energies
    type(DMat), intent(out) :: states
    type(InputParameters) :: input
    type(DMat) :: h
    type(EigenSolSymD) :: sol
    type(TwoBodyRelSpaceSpinIsoHOBasis) :: relspin
    type(TwoBodyRelSpaceIsoHOBasis) :: rel
    type(TwoBodyRelOpIso) :: vnn, kin_rel
    type(NNForceHOIsospin) :: vnn_relspin, U_relspin
    type(ThreeBodyJacIsoSpace) :: jacobi3, jacobi3_3n
    type(ABodyJacIsoSpace), target :: jacobiA, jacobiA_1
    type(ABodyJacIsoChan) :: jacobiA_ch
    type(ThreeBodyJacOpIso) :: T_3, vnn_3, v3n_3_tmp, v3n_3
    type(ABodyJacOpChanIso) :: T_A_ch, vnn_A_ch, v3n_A_ch
    type(ABodyJacOpIso) :: T_A1, vnn_A1, v3n_A1, T_A, vnn_A, v3n_A
    type(str) :: ramp, fn
    type(sys) :: s
    integer :: n_p, wunit=20
    logical :: NN_only=.false.

    if(params%NN_only) NN_only=.true.
    if(.not. params%genuine_3bf .and. params%renorm%val=="bare") NN_only=.true.

    call relspin%init(params%hw,params%N2max, params%J2max_NNint)

    call vnn_relspin%init(relspin)
    call U_relspin%init(relspin)
    input = params; input%N_ls = params%NAmax
    call vnn_relspin%setNNForceHOIsospin(U_relspin, input, t, z)

    call rel%init(params%hw, params%N2max, params%Jmax2)
    call vnn%init(rel,rel,s%str('NNint'))
    call vnn%SetTwoBodyScalarSpinBreaking(vnn_relspin)
    call kin_rel%init(rel,rel,s%str('kinetic'))
    call kin_rel%set()


    ! *** cfp ***
    write(*,*)
    write(*,*) "CFP part:"
    ramp = "flat" + s%str(params%N3max)
    write(*,"(2a)") "Space for CFP calc.: ", ramp%val
    call jacobi3%init(params%hw, 2*params%NAmax+3, params%NAmax, ramp, params%N3max, params%path_to_tmp_dir)
    if(params%particle_rank == 4) then
      call jacobiA_ch%init(j, p, t, params%NAmax, jacobi3, params%path_to_tmp_dir)
      call jacobi3%fin()
    end if

    if(params%particle_rank > 4) then
      do n_p = 4, params%particle_rank-1
        if(n_p == 4) then
          call jacobiA%init(params%NAmax, jacobi3, params%path_to_tmp_dir)
          call jacobi3%fin()
          if(params%particle_rank - 1 > 4) call jacobiA%fin()
        end if

        if(n_p > 4) then
          call jacobiA_1%init(params%NAmax, n_p-1, jacobiA, params%path_to_tmp_dir)
          call jacobiA%init(params%NAmax, n_p, jacobiA_1, params%path_to_tmp_dir)
          call jacobiA_1%fin()
          if(n_p /= params%particle_rank-1) call jacobiA%fin()
        end if
      end do
      call jacobiA_ch%init(j, p, t, params%NAmax, params%particle_rank, jacobiA, params%path_to_tmp_dir)
    end if
    ! *** cfp ***

    ! *** Hamiltonian ***
    write(*,*)
    write(*,*) "Hamiltonian part:"
    call jacobi3%init(params%hw, 2*params%NAmax+3, params%NAmax, ramp, params%N3max, params%path_to_tmp_dir)
    T_3 = get_3_body_jacobi_op_isospin(jacobi3, kin_rel, s%str("kinetic"))
    vnn_3 = get_3_body_jacobi_op_isospin(jacobi3, vnn, s%str("NNint"))
    if(.not. NN_only) then
      write(*,"(2a)") "Space for NNN: ",params%ramp%val
      call jacobi3_3n%init(params%hw, params%jmax3, params%NAmax, params%ramp, params%N3max, params%path_to_tmp_dir)
      call calc_3body_files(params)
      call v3n_3_tmp%init(jacobi3_3n, s%str("NNNint"))
      call v3n_3_tmp%set(params)

      v3n_3 = v3n_3_tmp%ChangeModelSpace(jacobi3)
      call v3n_3_tmp%fin()
      call jacobi3_3n%fin()
    end if

    call T_A_ch%init(jacobiA_ch, jacobiA_ch, s%str("kinetic"))
    call vnn_A_ch%init(jacobiA_ch, jacobiA_ch, s%str("NNint"))
    if(.not. NN_only) call v3n_A_ch%init(jacobiA_ch, jacobiA_ch, s%str("NNNint"))

    if(params%particle_rank == 4) then
      call T_A_ch%set(T_3)
      call vnn_A_ch%set(vnn_3)
      if(.not. NN_only) call v3n_A_ch%set(v3n_3)
      n_p = params%particle_rank
      h = dble(n_p-1) * 0.5d0 * T_A_ch%DMat
      h = h + (dble(n_p * (n_p-1)) / 6.d0) * vnn_A_ch%DMat
      if(.not. NN_only) h = h + (dble(n_p * (n_p-1) * (n_p-2))/6.d0) * v3n_A_ch%DMat
      call T_3%fin()
      call vnn_3%fin()
      call v3n_3%fin()
      call jacobi3%fin()
    end if

    if(params%particle_rank > 4) then
      do n_p = 4, params%particle_rank - 1
        if(n_p == 4) then
          call jacobiA%init(params%NAmax, jacobi3, params%path_to_tmp_dir)
          T_A = get_4_body_jacobi_op_isospin(jacobiA, T_3)
          vnn_A = get_4_body_jacobi_op_isospin(jacobiA, vnn_3)
          if(.not. NN_only) v3n_A = get_4_body_jacobi_op_isospin(jacobiA, v3n_3)
          call jacobi3%fin()
          call T_3%fin()
          call vnn_3%fin()
          if(.not. NN_only) call v3n_3%fin()
          call jacobi3%fin()
          cycle
        end if

        if(n_p > 4) then
          call jacobiA_1%init(params%NAmax, n_p-1, jacobiA, params%path_to_tmp_dir)
          T_A1 = T_A
          vnn_A1 = vnn_A
          call T_A1%SetModelSpace(jacobiA_1)
          call vnn_A1%SetModelSpace(jacobiA_1)
          if(.not. NN_only) then
            v3n_A1 = v3n_A
            call v3n_A1%SetModelSpace(jacobiA_1)
          end if
          call jacobiA%init(params%NAmax, n_p, jacobiA_1, params%path_to_tmp_dir)
          T_A = get_a_body_jacobi_op_isospin(jacobiA, T_A1)
          vnn_A = get_a_body_jacobi_op_isospin(jacobiA, vnn_A1)
          if(.not. NN_only) v3n_A = get_a_body_jacobi_op_isospin(jacobiA, v3n_A1)
          cycle
        end if
      end do
      call T_A_ch%set(T_A)
      call vnn_A_ch%set(vnn_A)
      if(.not. NN_only) call v3n_A_ch%set(v3n_A)
      n_p = params%particle_rank
      h = dble(n_p-1) * 0.5d0 * T_A_ch%DMat + &
          & (dble(n_p * (n_p-1)) / 6.d0) * vnn_A_ch%DMat
      if(.not. NN_only) then
        h = h + (dble(n_p * (n_p-1) * (n_p-2))/6.d0) * v3n_A_ch%DMat
      end if
    end if


    call T_A_ch%fin()
    call vnn_A_ch%fin()
    if(.not. NN_only) call v3n_A_ch%fin()
    call jacobiA_ch%fin()
    ! *** Hamiltonian ***

    call sol%init(h)
    call sol%DiagSym(h)
    states = sol%vec
    energies = sol%eig
    if(params%fn_fewbody_wave_function%val /= 'none') then
      fn = "wf_" + params%fn_fewbody_wave_function + ".wav"
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) states%m(:,1)
      close(wunit)
      fn = 'Hamil_' + params%fn_fewbody_wave_function + ".op"
      open(wunit, file=fn%val, form='unformatted', access='stream')
      write(wunit) h%m(:,:)
      close(wunit)
    end if

    call sol%fin()
    call h%fin()
    call vnn%fin()
    call rel%fin()
    call vnn_relspin%fin()
    call U_relspin%fin()
    call relspin%fin()
  end subroutine solve_A_body_hamil

  function get_3_body_jacobi_op_isospin(jacobi3, op2, opname) result(op)
    use TwoBodyRelOpsIso
    use ThreeBodyJacOpsIso
    type(ThreeBodyJacOpIso) :: op
    type(TwoBodyRelOpIso), intent(in) :: op2
    type(ThreeBodyJacIsoSpace), intent(in) :: jacobi3
    type(str), intent(in) :: opname
    integer :: ch, chbra, chket


    if(opname%val == "kinetic") then
      call op%init(jacobi3, opname)
      do ch = 1, jacobi3%GetNumberChannels()
        if(op%MatCh(ch,ch)%is_zero) cycle
        call op%MatCh(ch,ch)%set(op2)
      end do
      return
    end if

    if(opname%val == "NNint") then
      call op%init(jacobi3, opname)
      do ch = 1, jacobi3%GetNumberChannels()
        if(op%MatCh(ch,ch)%is_zero) cycle
        call op%MatCh(ch,ch)%set(op2)
      end do
      return
    end if

    call op%init(jacobi3, opname)
    do chbra = 1, jacobi3%GetNumberChannels()
      do chket = 1, jacobi3%GetNumberChannels()
        if(op%MatCh(chbra,chket)%is_zero) cycle
        call op%MatCh(chbra,chket)%set(op2)
      end do
    end do
  end function get_3_body_jacobi_op_isospin

  function get_4_body_jacobi_op_isospin(jacobi4, op3) result(op4)
    use ThreeBodyJacOpsIso
    type(ABodyJacOpIso) :: op4
    type(ThreeBodyJacOpIso), intent(in) :: op3
    type(ABodyJacIsoSpace), intent(in) :: jacobi4
    integer :: chbra, chket

    call op4%init(jacobi4, op3%GetOpName())
    do chbra = 1, jacobi4%GetNumberChannels()
      do chket = 1, jacobi4%GetNumberChannels()
        if(op4%MatCh(chbra,chket)%is_zero) cycle
        call op4%MatCh(chbra,chket)%set(op3)
      end do
    end do
  end function get_4_body_jacobi_op_isospin

  function get_a_body_jacobi_op_isospin(jacobiA, opA1) result(opA)
    type(ABodyJacOpIso) :: opA
    type(ABodyJacOpIso), intent(in) :: opA1
    type(ABodyJacIsoSpace), intent(in) :: jacobiA
    integer :: chbra, chket

    call opA%init(jacobiA, opA1%GetOpName())
    do chbra = 1, jacobiA%GetNumberChannels()
      do chket = 1, jacobiA%GetNumberChannels()
        if(opA%MatCh(chbra,chket)%is_zero) cycle
        call opA%MatCh(chbra,chket)%set(opA1)
      end do
    end do
  end function get_a_body_jacobi_op_isospin

  subroutine calc_3body_files(params)
    type(InputParameters), intent(in) :: params
    integer, allocatable :: jpt(:,:)
    integer :: n, j, p, t, nch, ch
    nch = ((params%jmax3 + 1) / 2) * 4
    allocate(jpt(3,nch))
    n = 0
    do t = 1, 3, 2
      do j = 1, params%jmax3, 2
        do p = 1, -1, -2
          n = n + 1
          jpt(:,n) = [j,p,t]
        end do
      end do
    end do

    do ch = 1, nch
      call calc_each_channel(ch)
    end do
    deallocate(jpt)
  contains
    subroutine calc_each_channel(ich)
      use MPIfunction
      use ClassSys, only: sys
      use ThreeBodyJacobiSpace
      use NNNForceHOIsospin
      integer, intent(in) :: ich
      type(ThreeBodyJacIsoChan) :: jac
      type(sys) :: s
      type(ThreeBodyJacOpChanIso) :: U
      type(NNNForceIsospin) :: v3
      type(str) :: f, fv, fut
      integer :: wunit = 21, runit = 22, Nmax

      Nmax = GetRampNmax(jpt(1,ich), params%ramp)
      f = jac%GetFileName(jpt(1,ich),jpt(2,ich),jpt(3,ich),Nmax,params%path_to_tmp_dir)

      write(*,'(a, i4, a, i2, a, i2, a, i2, a, i2)') &
          &  'myrank=', myrank, ',  :Calculating for J = ', &
          &   jpt(1,ich), '/2,   P = ', jpt(2,ich), ',   Tz = ', jpt(3,ich), &
          &   '/2, Nmax = ', Nmax

      if(.not. s%isfile(f)) then
        call jac%init(params%hw,jpt(1,ich),jpt(2,ich),jpt(3,ich),Nmax,params%path_to_tmp_dir)

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

      fv = u%GetFileName(jac,params%hw,s%str("NNNint"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)
      fut= u%GetFileName(jac,params%hw,s%str("UT"),params%genuine_3bf, params%regulator, &
          & params%regulatorpower, params%renorm, params%lambda_3nf_local, params%lambda_3nf_nonlocal, &
          & params%lambda, params%c1, params%c3, params%c4, params%cd, params%ce, params%J3max_initial_3nf, &
          & params%path_to_tmp_dir)

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
    end subroutine calc_each_channel
  end subroutine calc_3body_files

end module ABodyManager
